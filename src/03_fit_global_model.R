#--------------------------------------------
#-------------- Load packages ---------------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")

packages <- c( "dplyr", "stringr", "here", "qs","CoordinateCleaner","terra", "raster", "rnaturalearth", "rnaturalearthdata", "ggplot2","tidyterra", "dismo", "sdm", "caret", "viridisLite", "kableExtra","future", "future.apply","randomForest","earth", "progressr", "sf", "gbm", "PresenceAbsence","geosphere","arm", "RStoolbox", "ecospat", "viridis", "patchwork", "grid")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}

sdm::installAll()


#--------------------------------------------
#- Load helper functions and configurations -
#--------------------------------------------
source(here::here("src", "helper_functions.R"))
source(here::here("src", "00_configurations.R"))


#--------------------------------------------
#--- Load global occurrences and taxa info---
#--------------------------------------------
global<-qs::qread( paste0("./data/projects/",project,"/",project,"_occurrences.qs"))
taxa_info<-read.csv2(paste0("./data/projects/",project,"/",project,"_taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  dplyr::pull(acceptedTaxonKey)%>%
  unique()


#--------------------------------------------
#-------- Filter global occurrences----------
#--------------------------------------------
#remove unverified records
identificationVerificationStatus_to_discard <- c( "unverified",
                                                  "unvalidated",
                                                  "not validated",
                                                  "under validation",
                                                  "not able to validate",
                                                  "control could not be conclusive due to insufficient knowledge",
                                                  "1",
                                                  "uncertain",
                                                  "unconfirmed",
                                                  "douteux",
                                                  "invalide",
                                                  "non r\u00E9alisable",
                                                  "verification needed" ,
                                                  "probable",
                                                  "unconfirmed - not reviewed",
                                                  "validation requested",
                                                  "unconfirmed - plausible")

#enter value for max coordinate uncertainty in meters, default = 5000
#Reasoning: there are several invasive species datasets with default uncertainty of 5km
global.occ<-global %>%
  dplyr::filter(acceptedTaxonKey%in%accepted_taxonkeys) %>%   
  dplyr::filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters <= 5000) %>% 
  dplyr::filter(!str_to_lower(identificationVerificationStatus) %in% identificationVerificationStatus_to_discard)

#Remove coordinates that for both lon and lat values, have less than 4 decimal places
global.occ$lon_dplaces<-sapply(global.occ$decimalLongitude, function(x) decimalplaces(x))
global.occ$lat_dplaces<-sapply(global.occ$decimalLatitude, function(x) decimalplaces(x))
global.occ[global.occ$lon_dplaces < 4 & global.occ$lat_dplaces < 4 , ]<-NA
global.occ<-global.occ[ which(!is.na(global.occ$lon_dplaces)),]
global.occ<-within(global.occ,rm("lon_dplaces","lat_dplaces")) 


#--------------------------------------------
#------------ Define species group-----------
#--------------------------------------------
fish_orders <- c(
  "Acipenseriformes", "Albuliformes", "Amiiformes", "Anguilliformes", "Atheriniformes",
  "Batrachoidiformes", "Beloniformes", "Beryciformes", "Blenniiformes", "Carangiformes",
  "Characiformes", "Chimaeriformes", "Clupeiformes", "Cypriniformes", "Cyprinodontiformes",
  "Elopiformes", "Esociformes", "Gadiformes", "Gasterosteiformes", "Gonorynchiformes",
  "Gymnotiformes", "Holocentriformes", "Istiophoriformes", "Lampriformes", "Lophiiformes",
  "Mugiliformes", "Myliobatiformes", "Myxiniformes", "Notacanthiformes", "Ophidiiformes",
  "Osmeriformes", "Osteoglossiformes", "Perciformes", "Percopsiformes", "Petromyzontiformes",
  "Pleuronectiformes", "Polypteriformes", "Pristiformes", "Rajiformes", "Salmoniformes",
  "Scorpaeniformes", "Siluriformes", "Squaliformes", "Stomiiformes", "Syngnathiformes",
  "Tetraodontiformes", "Torpediniformes", "Trachiniformes", "Zeiformes"
)

global.occ <- global.occ %>%
  dplyr::mutate(Group = case_when(
    kingdom == "Plantae" ~ "Plants",
    class == "Aves" ~ "Birds",
    phylum == "Mollusca" ~ "Molluscs",
    class == "Amphibia" ~ "Amphibians",
    class == "Mammalia" ~ "Mammals",
    class %in% c("Crocodylia", "Testudines", "Sphenodontia", "Squamata") ~ "Reptiles",
    class == "Malacostraca" ~ "Malacostraca",
    class == "Insecta" ~ "Insects",
    order %in% fish_orders ~ "Fish",
    TRUE ~ NA_character_
  ))


#--------------------------------------------
#-------Prepare occurrence dataset-----------
#--------------------------------------------
global.occ.LL<-global.occ%>%
  rename(species= acceptedScientificName)%>%
  dplyr::select(decimalLongitude, decimalLatitude, species, acceptedTaxonKey, Group, coordinateUncertaintyInMeters) 
rm(global.occ, global)


#--------------------------------------------
#-----------Do coordinate cleaning-----------
#--------------------------------------------
# Clean coordinates based on their proximity to country centroids, capitals, biodiversity institutions, GBIF headquarters, and the 0/0 point
cleaned<-global.occ.LL%>%
  CoordinateCleaner::cc_cen(buffer=100) %>% # remove points within a buffer of 100m around country centroids, default 1km
  CoordinateCleaner::cc_cap(buffer=100) %>% # remove capitals centroids (buffer 100m), default 10km
  CoordinateCleaner::cc_inst(buffer=100) %>% # remove zoo and herbaria records buffer of 100 m around biodiversity institutes, default 100m
  CoordinateCleaner::cc_gbif(buffer=100)%>% #remove around GBIF headquarters in Copenhagen (buffer 100m), default 100m
  CoordinateCleaner::cc_zero() #Remove around the 0/0 point (buffer 0.5 degrees)


#---------------------------------------------
#--- Prepare data at 1km coord uncertainty ---
#---------------------------------------------
#These data will be used for the European model
cleaned_1km<-cleaned%>%
  filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters <= 1000)
 

#--------------------------------------------
#------- Load global climate rasters --------
#--------------------------------------------
# Only include files that start with "scaled_layer_" and end with .tif: 
scaled_files <- list.files(
  here::here("data", "external","climate","trias_CHELSA"),
  pattern = "^scaled_layer.*\\.tif$",
  full.names = TRUE)

# Load and stack
globalclimpreds_terra <- terra::rast(scaled_files)

#Decrease resolution to match coordinate uncertainty of global occurrences: use around 5km at equator by averaging
globalclimpreds_terra<- aggregate(globalclimpreds_terra, fact=5, fun=mean, na.rm=TRUE)


#--------------------------------------------
#--------Load European climate rasters-------
#--------------------------------------------
euboundary<-sf::st_read(here::here("data", "external", "GIS", "Europe", "EUROPE.shp")) 
# Convert sf boundary to SpatVector
euboundary_vect <- terra::vect(euboundary)


#---------------------------------------------
#- Remove NA pixels and mask to European ext -
#---------------------------------------------
# This ensures all rasters have the same NA structure
na_mask_globalclimpreds_terra <- anyNA(globalclimpreds_terra)
globalclimpreds_terra  <- terra::mask(
  globalclimpreds_terra,
  na_mask_globalclimpreds_terra,
  maskvalue = 1,
  filename = "data/external/climate/masked_globalclimpreds.tif",
  overwrite = TRUE,
  wopt = list(gdal = c("COMPRESS=LZW"))
)

# Crop and mask scaled_stack to European extent
eu_climpreds.10 <- terra::crop(globalclimpreds_terra, euboundary_vect)
eu_climpreds.10 <- terra::mask(eu_climpreds.10, euboundary_vect)


#--------------------------------------------
#------ Load future climate rasters ----------
#--------------------------------------------

for (period in c("2041-2070","2071-2100")){
  for(scenario in c("ssp126", "ssp370", "ssp585")){
    
    # List future raster files
    future_files<- list.files( here::here("data", "external", "climate", "chelsa_future", period,scenario), pattern = "\\.tif$", full.names = TRUE)
    
    #Stack them together
    future_stack <- rast(future_files)
    
    #Aggregate at a resolution of 5km
    future_stack<- aggregate(future_stack, fact=5, fun=mean, na.rm=TRUE)
    
    #Reproject future stack to match CRS and resolution of eu_climpreds.10
    #TODO: Regenerate these rasters so this is not necessary
    future_aligned <- terra::project(
      future_stack,
      eu_climpreds.10,
      method = "bilinear")
    
    # Get mask from future stack NA structure
    na_mask_future <- anyNA(future_aligned)
    
    # Mask future stack with its own NA structure
    future_aligned_masked <- terra::mask(
      future_aligned,
      na_mask_future,
      maskvalue = 1)  
    
    #Assign correct name to raster stack
    assign(paste0(scenario,"_",period), future_aligned_masked)
    
    #Clean up
    rm(future_aligned, future_aligned_masked, future_stack)
    
  }
}


#--------------------------------------------
#--------- Load shape of the world ----------
#--------------------------------------------
world<-rnaturalearth::ne_countries(scale=50)


#--------------------------------------------
#--------------Load ecoregions --------------
#--------------------------------------------
wwf_eco_biome<-sf::st_read(here::here("./data/external/GIS/official/newRealms.shp")) 


# Dissolve by BIOME
#wwf_eco_biome <- wwf_eco %>%
#  group_by(BIOME) %>%
#  summarise(geometry = sf::st_union(geometry), .groups = "drop") %>%
#  sf::st_as_sf()

# Optionally, make geometry valid
#wwf_eco_biome <- sf::st_make_valid(wwf_eco)

#--------------------------------------------
#-------Load file paths to bias grids -------
#--------------------------------------------
bias_grid_folder<-here::here("data","external", "bias_grids")
bias_grid_paths <- list(
  Plants = here::here(bias_grid_folder, "plants_10km_bias_layer_log.tif"), #0-13.24
  Amphibians = here::here(bias_grid_folder, "amphibians_10km_bias_layer_log.tif"),#0-12.06
  Birds = here::here(bias_grid_folder, "birds_1deg_min5.tif"),#5-1703018
  Mammals = here::here(bias_grid_folder, "mammals_10km_bias_layer_log.tif"), #0-13.36
  Molluscs = here::here(bias_grid_folder, "mollusca_10km_bias_layer_log.tif"),#0-12.48
  Reptiles = here::here(bias_grid_folder, "reptiles_10km_bias_layer_log.tif"), #0-11.34
  Fish = here::here(bias_grid_folder, "fish_10km_bias_layer_log.tif"),#0-14.67
  Malacostraca = here::here(bias_grid_folder, "malacostraca_10km_bias_layer_log.tif"),#0-13.12
  Insects = here::here(bias_grid_folder, "insects_10km_bias_layer_log.tif")) #0-15.78


#--------------------------------------------
#------- Split dataframe by taxonkey --------
#--------------------------------------------
sort(unique(cleaned$species))
split_df<-split(cleaned,cleaned$species)
split_df_all_occs <- split_df                                 #keep a copy with all occurrences
#split_df <- thinOccurrences50(split_df_all_occs)  #split_df as filename for the rarefied occurrences


#--------------------------------------------
#---------------- Clean up ------------------
#--------------------------------------------
rm(global.occ.LL,cleaned)
gc()

names(split_df)


#--------------------------------------------
#-------Start loop for SDM modelling --------
#--------------------------------------------
with_progress({
  p <- progressr::progressor(along = 1:length(split_df)) 
  for(i in 1:length(split_df)) {
    
    #--------------------------------------------
    #------------- Track progress ---------------
    #--------------------------------------------
    p()
    
    #--------------------------------------------
    #----------Load species data ----------------
    #--------------------------------------------
    species <- names(split_df)[i]
    print(species)
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)  # Extract first two words of species name
    global.occ.LL.cleaned<-split_df[[i]]
    taxonkey<-unique(global.occ.LL.cleaned$acceptedTaxonKey)
    speciesgroup<-unique(global.occ.LL.cleaned$Group)
    global.occ.LL.cleaned<-global.occ.LL.cleaned %>%
      dplyr::select(c(decimalLongitude,decimalLatitude))
    global.occ_1KM<-cleaned_1km %>%
      filter(acceptedTaxonKey == taxonkey)
    
    #Generate file for informing PA selection containing all occurrences (no thinning, in case we thinned split_df)
    for_PA_selection <- split_df_all_occs[[i]] %>%
      dplyr::select(c(decimalLongitude, decimalLatitude))%>%
      sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                   crs = 4326)
    
    
    #--------------------------------------------
    #-------------Create folders-----------------
    #--------------------------------------------
    # Define the folder paths
    folder_paths<-list(list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "Rasters", "Interim"),
                            "name"= "Rasters/Interim"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "Rasters", "Global", "2041-2070","ssp126"),
                            "name"= "Rasters/Global/2041-2070/ssp126"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "Rasters", "Global", "2041-2070","ssp370"),
                            "name"= "Rasters/Global/2041-2070/ssp370"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "Rasters", "Global", "2041-2070","ssp585"),
                            "name"= "Rasters/Global/2041-2070/ssp585"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "Rasters", "Global", "2071-2100","ssp126"),
                            "name"= "Rasters/Global/2071-2100/ssp126"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "Rasters", "Global", "2071-2100","ssp370"),
                            "name"= "Rasters/Global/2071-2100/ssp370"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "Rasters", "Global", "2071-2100","ssp585"),
                            "name"= "Rasters/Global/2071-2100/ssp585"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "PDFs", "Global"),
                            "name"= "PDFs"),
                       list("path"=here::here("data", "projects", project, paste0(first_two_words, "_", taxonkey), "PNGs", "Global"),
                            "name"= "PNGs"))
    
    # Check and create each folder if necessary
    lapply(folder_paths, function(folder){
      create_folder(folder$path, folder$name)
    })
    
    
    #--------------------------------------------
    #------ Process occurrence data -----
    #--------------------------------------------
    #Keep only one occurrence per grid cell
    global.occ.LL.cleaned <- remove_duplicates(occurrences = global.occ.LL.cleaned, rast_template = globalclimpreds_terra)
    
    #Remove occurrences within grid cells with NA values
    global.occ.sf <- remove_nodata_occurrences(occurrences = global.occ.LL.cleaned, rast_template=globalclimpreds_terra, crs=4326)
    
    #add column indicating species presence (1) for modeling
    global.occ.sf$species <- rep(1, nrow(global.occ.sf)) 
   
    
    #-----------------------------------------------
    #------ Limit to 10,000 occupied grid cells ----
    #-----------------------------------------------
    if(nrow(global.occ.sf) > 10000){
      if(occurrence_thinning_method == "random"){
      set.seed(101) 
      global.occ.sf <- global.occ.sf[sample(nrow(global.occ.sf), 10000, replace=FALSE), ]
      }else if (occurrence_thinning_method == "kmeans_clustering"){
        #Extract environmental data in each occurrence grid cell
        env_data <- terra::extract(globalclimpreds_terra, global.occ.sf, ID = FALSE)
        
        # K-means clustering
        set.seed(101)
        clust <- kmeans(env_data, centers = n_clusters,iter.max = 10, nstart = 1)$cluster
        occ_env <- cbind(global.occ.sf, env_data, clust)
        
        # Sample within clusters
        max_per_cluster <- 10000/n_clusters
        row_sample <- sapply(1:10000, function(x) {
          rowids <- which(clust == x)
          sample(rowids, min(max_per_cluster, length(rowids)), replace = FALSE)
        })
        
        #Get final dataframe
        global.occ.sf <- occ_env[row_sample,] %>%
          dplyr::select(decimalLongitude, decimalLatitude, geometry, species)

      }
    }
    
    
    #-------------------------------------------------------
    #-Don't fit model if less than 20 global presences -----
    #-------------------------------------------------------   
    if(nrow(global.occ.sf)<20){
      warning(paste0("Skipping species ", species, " because the number of occurrences is less than 20 (n =",nrow(global.occ.sf),")"))
      next  # Skip the rest of the loop and move to the next iteration
    }
    
    
    #--------------------------------------------
    #------ Plot distribution of occurrences ----
    #--------------------------------------------
    # ggplot()+
    # geom_sf(data = world,  colour = "black", fill = NA)+
    # geom_point(data=global.occ.sf, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
    # labs(x="Longitude", y="Latitude")+
    # theme_bw()
    
    
    #--------------------------------------------
    #- Select ecoregions containing occurrences -
    #--------------------------------------------
    
    # Ensure valid geometries
    #wwf_eco_biome <- sf::st_make_valid(wwf_eco_biome)
    
    # Disable S2 geometry engine to avoid topological issues
    sf::sf_use_s2(FALSE)
    
    # Keep only biome polygons that intersect at least one occurrence point
    has_occurrence <- lengths(sf::st_intersects(wwf_eco_biome, global.occ.sf)) > 0
    wwf_ecoSub1 <- wwf_eco_biome[has_occurrence, ]
    sf::sf_use_s2(TRUE)
    
    #--------------------------------------------
    #------------- Plot ecoregions --------------
    --------------------------------------------
    # ggplot()+
    # geom_sf(data = world,  colour = "black", fill = NA)+
    # geom_sf(data=wwf_ecoSub1, fill="#f7786f")+
    # geom_point(data=global.occ.sf, aes(decimalLongitude, decimalLatitude), color="blue")+
    # labs(x="Longitude", y="Latitude")+
    # theme_bw()

    
    #--------------------------------------------
    #------ Import right bias grid --------------
    #--------------------------------------------
    if (speciesgroup %in% names(bias_grid_paths)) {
      biasgrid <- terra::rast(bias_grid_paths[[speciesgroup]])
      if(speciesgroup %in% c("Amphibians", "Molluscs", "Mammals", "Reptiles","Birds","Plants","Fish","Malacostraca","Insects")){
        # Resample biasgrid to match the resolution of globalclimpreds_terra
        biasgrid <- terra::resample(biasgrid, globalclimpreds_terra[[1]], method="bilinear")
      }
    } else {
      next("No bias grid available for this species. Species has to be one of the following: Amphibians, Molluscs, Mammals, Reptiles, Birds, Plants, Fish, Malacostraca, or Insects.")
    }

    
    #--------------------------------------------
    #------------ Process  bias grid ------------
    #--------------------------------------------
    #Mask biasgrid with climate layers (no PA can be selected in NA climate pixels)
    biasgrid_log <- terra::mask(biasgrid, globalclimpreds_terra[[1]])
    
    # Rescale raster values to range from 1 to 20
    min_val <- global(biasgrid_log, fun = "min", na.rm = TRUE)[[1]]
    max_val <- global(biasgrid_log, fun = "max", na.rm = TRUE)[[1]]
    biasgrid <- ((biasgrid_log - min_val) / (max_val - min_val)) * 19 + 1
    
    #Mask biasgrid with biomes with occurrences
    wwf_ecoSub1_ext<-terra::ext(wwf_ecoSub1) 
    wwf_ecoSub1_vector <- terra::vect(wwf_ecoSub1) 
    biasgrid_crop <- terra::crop(biasgrid, wwf_ecoSub1_ext) 
    biasgrid_sub <- terra::mask(biasgrid_crop, wwf_ecoSub1_vector)
    
    
    #--------------------------------------------
    #---------------Visualize biasgrid-----------
    #--------------------------------------------
    # ggplot()+
    # geom_sf(data = world,  colour = "black", fill = NA)+
    # geom_spatraster(data=biasgrid_sub)+
    # scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange")+
    # labs(x="Longitude", y="Latitude")+
    # theme_bw()


    #--------------------------------------------
    #---------- Generate pseudoabsences ---------
    #--------------------------------------------
    
    #LIMIT TO AREA OF 100KM AROUND OCCURRENCES
    #global.occ.sf_buffer_100km <- st_buffer(global.occ.sf, dist = 100000)
    #global.occ.sf_buffer_100km <- vect(st_union(global.occ.sf_buffer_100km))
    # Convert raster from 'raster' to 'terra' format if needed
    #if (inherits(biasgrid_sub_raster, "Raster")) {
    #  biasgrid_sub_raster <- rast(biasgrid_sub_raster)
    #}
    #biasgrid_sub_raster<-terra::mask(biasgrid_sub_raster,global.occ.sf_buffer_100km)
    #biasgrid_sub_raster_raster <- raster::raster(biasgrid_sub_raster)
    
    
    #---------Mask cells that contain occurrences---------
    # Convert sf to terra-compatible vector
    for_PA_vect <- terra::vect(for_PA_selection)
    # Identify raster cells that correspond to these points
    cells_to_exclude <- terra::cellFromXY(biasgrid_sub, terra::crds(for_PA_vect))
    # Set those cells to NA
    biasgrid_sub[cells_to_exclude] <- NA
    
    #--------------Generate pseudoabsences-----------------
    #TODO check if alternative is needed when not enough points can be selected
    #Before we used ecoregions grid, without biasgrid mask as alternative
    set.seed(728)
    global_points <- terra::spatSample(
      biasgrid_sub,
      size = 10000,
      method = "weights",     # weighted random sampling
      as.points = TRUE,       # return SpatVector of points
      na.rm = TRUE            # ignore NA pixels
    )
    
    
    #--------------------------------------------
    #--- Create presence-pseudoabsence dataset---
    #--------------------------------------------
    # Extract coordinates from SpatVector
    coords <- terra::crds(global_points)
    
    # Add coordinates and convert
    global_pseudoAbs <- coords %>%
      as.data.frame() %>%
      mutate(species = 0) %>%
      rename(decimalLongitude = x,
             decimalLatitude = y) %>%
      sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326, remove = F) 
    
    #Combine with presence data
    global_presabs <- rbind(global.occ.sf, global_pseudoAbs)
    
    #Clean up
    rm(global_points, global_pseudoAbs)
    
    
    #--------------------------------------------
    #--Visualize presence-pseudoabsence dataset--
    #--------------------------------------------
    # ggplot()+
    #   geom_sf(data = world,  colour = "black", fill = NA)+
    #   geom_spatraster(data=biasgrid_sub)+
    #   geom_point(data=global_presabs, aes(x=decimalLongitude, y=decimalLatitude, color=factor(species)))+
    #   scale_color_manual(values=c("red", "green"), name="Occurrence type")+
    #   scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange", name="Weight")+
    #   labs(x="Longitude", y="Latitude")+
    #   theme_bw()

    
    #--------------------------------------------
    #---- Extract climate data for modelling-----
    #--------------------------------------------
    global.data.df <- sdm::sdmData(species~.,train=vect(global_presabs),predictors=globalclimpreds_terra)%>%
      as.data.frame()
    
    
    #--------------------------------------------
    #--- Remove highly correlated predictors---
    #--------------------------------------------
    # Calculate correlation matrix (excluding rID and species)
    correlationMatrix <- cor(global.data.df[, -c(1, 2)])
    
    # Identify highly correlated variables (cutoff = 0.7)
    highlyCorrelated <- caret::findCorrelation(correlationMatrix, cutoff = 0.7, exact = TRUE, names = TRUE)
    
    # Remove highly correlated predictors and rID, and prepare species factor
    global.data.df.uncor <- global.data.df %>%
      dplyr::select(-all_of(highlyCorrelated), -rID) %>%
      dplyr::mutate(
        species = as.factor(species),
        species = recode_factor(species, '0' = "absent", '1' = "present"),
        species = relevel(species, ref = "present")
      )
    
    # Remove them from climate stack
    globalclimpreds_terra_selection <- globalclimpreds_terra %>%
      subset(!names(globalclimpreds_terra) %in% highlyCorrelated)
    
    eu_climpreds.10_selection <- eu_climpreds.10 %>%
      subset(!names(eu_climpreds.10) %in% highlyCorrelated)
      
    
    #--------------------------------------------
    #--- Run multiple machine learning models ---
    #--------------------------------------------
    
    #Define prevalence ratio
    n1 <- nrow(global.occ.sf)  # presences
    n0 <- 10000                # pseudoabsences (adjust to your setup if different)
    prev_ratio <- n1 / n0
    
    #define methods and data
    sdm_data <- sdm::sdmData(species~.,train=vect(global_presabs),predictors=globalclimpreds_terra_selection) 
    methods <- c("glm", "gam", "bioclim", "brt", "rf", "glmpoly", "mars", "maxent", "fda", "cart") 
    
    #run model
    model <- sdm(
      species ~ ., data = sdm_data,
      methods = methods  # 10 models
    )
    print(model)
   
    
    #--------------------------------------------
    #---  Make predictions using each model  ---
    #-------------------------------------------- 
    # Get model info
    info <- sdm::getModelInfo(model)
    
    #Get presence data and their associated climatological values
    pres_data <- subset(global.data.df.uncor, species == "present")
    pres_features <- pres_data[, -which(names(pres_data) == "species")]
    
    #Create empty list to store models in
    modeloutput<-list()
    
    system.time({
      for(modelmethod in methods){
        
        print(modelmethod)
        
        #Create raster with predictions for Europe
        pred_raster <- predict(model,
                               newdata = eu_climpreds.10_selection,
                               method = modelmethod)
        
        # Get model IDs
        model_ids <- info$modelID[info$method == modelmethod]
        
        # Subset using those IDs
        method_model <- model[[model_ids]]  
        
        #Apply the transformation to the raster
        fav_raster <- favourability_from_prob(pred_raster, prev_ratio)
        
        #Get threshold
        fitted_model <- predict(method_model, newdata = pres_features, type = "response")
        fitted_favourability <- favourability_from_prob(fitted_model[[1]], prev_ratio)
        
        # 1% minimum training presence threshold
        threshold_1pct <- quantile(fitted_favourability, probs = 0.01, na.rm = TRUE)
        
        #5% minimum training presence threshold
        threshold_5pct <- quantile(fitted_favourability, probs = 0.05, na.rm = TRUE)
        
        # Binarize rasters using the thresholds
        binary_1pct <- fav_raster >= threshold_1pct 
        binary_5pct <- fav_raster >= threshold_5pct
        
        # Plot
       plot(fav_raster, main = paste0("Favourability (", modelmethod,")"))
       plot(binary_1pct, main = paste0("Binary map ",modelmethod," (1% threshold)"))
       plot(binary_5pct, main = paste0("Binary map ",modelmethod," (5% threshold)"))
        
        #Store
        modeloutput[[modelmethod]]<-list(fav_raster=fav_raster,
                                         binary1pct=binary_1pct,
                                         binary5pct=binary_5pct,
                                         model=method_model)
        
        rm(fav_raster, binary_1pct, binary_5pct, method_model)
      }
    })
   
    
    #---------------------------------------------
    #-- Create Ensemble model using PCAm method --
    #---------------------------------------------

    # List favourability rasters
    fav_rasters_list <- lapply(modeloutput, function(x) x$fav_raster)
    
    # Combine into a SpatRaster stack
    fav_stack <- terra::rast(fav_rasters_list)
    
    # Assign layer names based on model methods
    names(fav_stack) <- names(modeloutput)
    
    #make PCA
    pca_result <- rasterPCA(fav_stack, nSamples = NULL, spca = FALSE, maskCheck = TRUE)
    
    #-----------------GET TOP 5 variance models----------------
    # Step 1: Recover original raster stack used in rasterPCA
    fav_stack <- eval(pca_result$call$img)
    
    # Step 2: Extract PC1 loadings from princomp object
    loadings <- pca_result$model$loadings[, 1]  # Comp.1 = PC1
    names(loadings) <- rownames(pca_result$model$loadings)
    
    # Step 3: Convert raster stack to matrix (rows = pixels, cols = models)
    fav_matrix <- as.matrix(fav_stack)
    
    # Step 4: Calculate variance along PC1 for each model
    model_variances <- setNames(numeric(nlyr(fav_stack)), names(fav_stack))
    
    for (lyr in 1:nlyr(fav_stack)) {
      model_vals <- fav_matrix[, lyr]
      centered <- model_vals - mean(model_vals, na.rm = TRUE)
      projection <- centered * loadings[lyr]
      model_variances[lyr] <- var(projection, na.rm = TRUE)
    }
    
    # Step 5: Select top 5 models with highest variance on PC1
    top5_models <- names(sort(model_variances, decreasing = TRUE))[1:5]
    cat("Top 5 models by variance along PC1:\n")
    
    # Get model IDs
    top_ids <- info$modelID[info$method %in% top5_models]
    
    # Subset using those IDs
    top5models <- model[[top_ids]]  
    
    # Step 6: Subset fav_stack to top 5 layers
    top5_stack <- subset(fav_stack, top5_models)
    
    # Step 7: Compute pixel-wise median = consensus model
    consensus_median <- app(top5_stack, median)
    
    # Step 8: Plot result
    plot(consensus_median, main = "Consensus (Median of Top 5 Models by Variance on PC1)")
    
    
    #---------------------------------------------
    #-- Get response curves of 5 selected models -
    #---------------------------------------------
    response_list<-list()
    varimp_list<-list()
    
    for(topmethod in top5_models){
      # Get model id
      id <- info$modelID[info$method == topmethod]
      
      #Get response curve
      response_curves<-sdm::getResponseCurve(model,id)@response
      
      #Get variable importance
      varimp<-sdm::getVarImp(model,id)@varImportance
      
      #Store
      response_list[[topmethod]]<-response_curves
      varimp_list[[topmethod]]<-varimp
    }
    
    # Convert list to a dataframe
    response_df <- imap_dfr(response_list, function(model_list, model_name) {
      imap_dfr(model_list, function(df, var_name) {
        response_df <- df %>%
          setNames(c("Predictor_value", "Response"))%>%
          mutate( Algorithm = model_name,
                  Predictor = var_name)})
    }) %>%
      dplyr::select(Algorithm,Predictor, Predictor_value, Response)
    
    
    varimp_df <- imap_dfr(varimp_list, function(df, model_name) {
      df %>%
        setNames(c("Predictor", "corTest" , "AUCtest"))%>%
        dplyr::mutate(Algorithm = model_name)
    })%>%
      dplyr::select(Algorithm,Predictor, corTest, AUCtest)
    
    
    # Plot response curves
    response_plot <- ggplot(response_df, aes(x = Predictor_value,
                                             y = Response, 
                                             color = Algorithm)) +
      geom_line(size=0.8) +
      facet_wrap(~ Predictor, scales = "free_x")+
      labs(title= "Climatological response curves" ,x= "Predictor value")+
      theme_bw()
    
    # Plot variable importance 
    varimp_plot <- ggplot(varimp_df, aes(x = Predictor, y = corTest)) +
      geom_col(fill = "steelblue") +
      coord_flip() +  #horizontal bars
      facet_wrap(~ Algorithm) +  
      geom_hline(yintercept = 0, color = "black") + 
      labs(
        x = "Variable",
        y = "Importance",
        title = "Variable importance per model"
      ) +
      theme_bw()
    
    #Write ggplots to pngs
    PNG_folder <- here::here("data","projects", project, paste0(first_two_words, "_", taxonkey), "PNGs")
    
    for(plot in c("varimp_plot", "response_plot")){
      
      #Define filename
      type <- switch(plot,
                     "varimp_plot" = "variable_importance",
                     "response_plot" =  "response_curves")
      PNG_file <- paste(first_two_words, "_", taxonkey,"_",type, "_global.png", sep = "")
      
      #Get plot
      gg_plot<-get(plot)
      
      #Save plot
      ggplot2::ggsave(filename = PNG_file, plot = gg_plot, 
                      device = "png", width =8.27 , height = 5.845, path= PNG_folder )
      #Print
      print(paste(PNG_file," has been created.", sep="")) 
    }
    
    
    #------------------------------------------
    #-- Create map with ensemble suitability --
    #------------------------------------------
    #Add EU occurrence points and EU Boyce-index value
    
    #Extract all raster values (excluding NAs)
    all_suit_vals <- values(consensus_median)
    all_suit_vals <- all_suit_vals[!is.na(all_suit_vals)]
    
    #Extract suitability values at occurrence locations
    occ_suit_vals <- terra::extract(consensus_median, vect(global.occ.sf))[,2]
    occ_suit_vals <- occ_suit_vals[!is.na(occ_suit_vals)]
    
    #Compute Boyce only if there are enough occurrences
    if (length(occ_suit_vals) > 0) {
      boyce_result <- ecospat.boyce(
        fit = all_suit_vals,
        obs = occ_suit_vals,
        nclass = 0
      )
      boyce_val <- round(boyce_result$cor, 3)
    } else {
      warning(paste("No EU occurrences available to calculate Boyce index for", species))
      boyce_val <- "NA (no EU data)"
      }
    
    #Prepare title and species extension
    species_title <- gsub("_", " ", first_two_words)
    rest_of_name <- if (grepl("^\\S+\\s+\\S+$", species)) "" else sub("^\\S+\\s+\\S+\\s+", "", species)
    
    #Define PNG and PDF folders
    PDF_folder <- file.path("./data/projects", project, paste0(first_two_words, "_", taxonkey), "PDFs")
    PNG_folder<-file.path("./data/projects", project, paste0(first_two_words, "_", taxonkey), "PNGs")
    
    for (occs in list(NULL, global.occ.sf)){
      # Export PDFs with and without occurrences plotted
      exportPDF(predictions = consensus_median,
                taxonName = first_two_words,
                nameExtension= rest_of_name,
                dataType = "Suit",
                taxonNameTitle = species_title,
                taxonKey = taxonkey,
                scenario = "hist",
                regionName = "Global",
                returnPredictions = FALSE,
                returnPNG = FALSE,
                occ_data=occs,
                exportPNG=TRUE,
                LabelValue= boyce_val,
                LabelName="Boyce index")
    }

    
    #------------------------------------------
    #------------ Create binary map -----------
    #------------------------------------------
    # Get predictor values at occurrence points
    predictors_only <- global.data.df.uncor%>%
      dplyr::filter(species=="present")%>%
      dplyr::select(-species)

    
    # Predict for top 5 models
    pred_vals <- list()
    for (method in top5_models) {
      pred_vals[[method]] <- predict(model, newdata = predictors_only, method = tolower(method))
    }
    
    # Favourability transformation
    fav_vals <- lapply(pred_vals, function(p) favourability_from_prob(p[[1]], prev_ratio))
    
    binary_maps<-list()
    for (probs in c(0.05, 0.01)){
    
    #Define mtp_pct
    mtp_pct <- switch(as.character(probs),
                     "0.05" = "5%",
                     "0.01" = "1%")  
      
    # Thresholds
    thr <- sapply(fav_vals, function(fv) quantile(fv, probs = probs, na.rm = TRUE))
    
    # Averages
    mean_pct <- mean(thr, na.rm = TRUE)
    cat(paste0("Mean ",mtp_pct," minimum training presence threshold: ", round(mean_pct, 4), "\n"))
    
    # Binary raster using MTP threshold
    binary_map_pct <- consensus_median >= mean_pct  
    binary_map_pct <- as.factor( binary_map_pct*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
    levels( binary_map_pct) <- data.frame(ID = c(0, 1),
                                          class = c("Absent", "Present"))

    # Calculate sensitivity in Europe
    occ_values <- terra::extract(binary_map_pct, vect(global.occ.sf))[,2]  # Extract raster values (2nd column) at occurrence locations
    cat("Occurrences in suitable cells:", sum(occ_values == "Present", na.rm = TRUE), "\n")
    cat("Occurrences in unsuitable cells:", sum(occ_values == "Absent", na.rm = TRUE), "\n")
    cat("Occurrences in NA cells:", sum(is.na(occ_values)), "\n")
    global_EU_sensitivity <- sum(occ_values == "Present", na.rm = TRUE) / sum(occ_values %in% c("Present", "Absent"), na.rm = TRUE)
    
    
    # export as PDF and PNG with and without occurrences plotted and return as PNG
    for (occs in list(NULL, global.occ.sf)){
      exportPDF(predictions = binary_map_pct,
                taxonName = first_two_words,
                nameExtension= rest_of_name,
                dataType = "Binary",
                taxonNameTitle = species_title,
                taxonKey = taxonkey,
                scenario = "hist",
                regionName = "Global",
                returnPredictions = FALSE,
                returnPNG = TRUE,
                occ_data=occs,
                exportPNG=TRUE,
                LabelValue= mtp_pct,
                LabelName="MTP threshold",
                Label2Value=round(global_EU_sensitivity,3),
                Label2Name="Sensitivity")
    }
    
    binary_maps[[mtp_pct]]<-list(binary_raster=binary_map_pct,
                                 EU_sensitivity=global_EU_sensitivity,
                                 mean_MTP=mean_pct)
    rm(binary_map_pct)
    }
    

    #--------------------------------------------
    #-- Prepare global_presabs for export--------
    #--------------------------------------------  
    #Decimallon and decimalLat are converted to x and y, geometry is dropped
    global_presabs<-global_presabs%>%
      dplyr::select(decimalLongitude, decimalLatitude)%>%
      dplyr::rename("x"= decimalLongitude,
             "y"= decimalLatitude)%>%
      sf::st_drop_geometry()
    
    
    #--------------------------------------------
    #-- Export results as .qs list
    #--------------------------------------------
    globalmodels <-list(species = species, #Species name
                        taxonkey = taxonkey, #Species taxonkey
                        global_data_df_uncor=global.data.df.uncor, #Data used to fit the model (climate data for each presence/pseudoabsence)
                        global_presabs = global_presabs,#xy coordinates of presences and pseudoabsences used to fit the models
                        occurrences=for_PA_selection, #All raw occurrences
                        occurrences5km = global.occ.sf, #Processed occurrences used to fit the models
                        occurrences1km = global.occ_1KM,
                        mtp_5_threshold = binary_maps$`5%`$mean_MTP,
                        mtp_1_threshold = binary_maps$`1%`$mean_MTP,
                        sdm_model = model,
                        pca_result = pca_result,
                        top5_models = top5_models,
                        final_boyce = boyce_val,
                        global_EU_sensitivity_5pct_threshold =  binary_maps$`5%`$EU_sensitivity,
                        global_EU_sensitivity_1pct_threshold =  binary_maps$`1%`$EU_sensitivity,
                        response_df = response_df,
                        varimp_df = varimp_df,
                        top5models = top5models)
                        
    qs::qsave(globalmodels, paste0("./data/projects/",project,"/",first_two_words,"_",taxonkey,"/Global_model_",first_two_words,"_",taxonkey,".qs"))
    
    
    #--------------------------------------------
    #--Export raster layers in folder "rasters"--
    #--------------------------------------------
    #We don't store them in .qs file as some important metadate would be stored in a temp folder, which would be removed after a while 
    biasgrid_file<- file.path("./data/projects",project,paste0(first_two_words,"_",taxonkey),"Rasters","Interim",paste0("Biasgrid_",first_two_words,"_",taxonkey,".tif"))
    ensemble_median_file <- file.path( "./data/projects", project, paste0(first_two_words, "_", taxonkey),"Rasters", "Global", paste0("Ensemble_median_", first_two_words, "_", taxonkey, ".tif"))
    
    terra::writeRaster(biasgrid_sub, filename = biasgrid_file, overwrite = TRUE)
    terra::writeRaster(consensus_median, filename = ensemble_median_file, overwrite = TRUE)
  
    
    #--------------------------------------------
    #------------------ Clean up-----------------
    #--------------------------------------------
    print(paste("Global model has been created for", species))
    
    rm(list = setdiff(ls(), c("p","wwf_eco_biome","eu_climpreds.10", "split_df",  "decimalplaces", "globalclimpreds_terra","bias_grid_paths", "i", "world", "project", "create_folder", "split_df_all_occs", "exportPDF", "remove_duplicates", "remove_nodata_occurrences", "favourability_from_prob", "cleaned_1km")))
   
  }
})


#--------------------------------------------
#---------- Clean R environment--------------
#--------------------------------------------
rm(list = ls())

