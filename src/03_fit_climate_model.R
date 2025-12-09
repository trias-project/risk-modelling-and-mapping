#--------------------------------------------
#-------------- Load packages ---------------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")#Prevent terra from writing aux.xml files

packages <- c( "dplyr", "stringr", "here", "qs","CoordinateCleaner","terra", "raster", "rnaturalearth", "rnaturalearthdata", "ggplot2","tidyterra", "dismo", "sdm", "caret", "viridisLite", "kableExtra","future", "future.apply","randomForest","earth", "progressr", "sf", "gbm", "PresenceAbsence","geosphere","arm", "RStoolbox", "ecospat", "viridis", "patchwork", "grid", "purrr", "magick")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}

sdm::installAll()


#--------------------------------------------
#- Load helper functions and configurations -
#--------------------------------------------
source(file.path("src", "helper_functions.R"))
source(file.path("src", "00_configurations.R"))


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
  dplyr::rename(species= acceptedScientificName)%>%
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
  dplyr::filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters <= 1000)


#--------------------------------------------
#------- Load global climate rasters --------
#--------------------------------------------
# Only include files that start with "scaled_layer_" and end with .tif: 
scaled_files <- list.files(
  file.path("data", "external","climate","chelsa_current"),
  pattern = "^scaled_layer.*\\.tif$",
  full.names = TRUE)

# Load and stack
globalclimpreds_terra <- terra::rast(scaled_files)

#Decrease resolution to match coordinate uncertainty of global occurrences: use around 5km at equator by averaging
globalclimpreds_terra<- terra::aggregate(globalclimpreds_terra, fact=5, fun=mean, na.rm=TRUE)


#--------------------------------------------
#--------Load European boundary layer -------
#--------------------------------------------
euboundary <- terra::rast(file.path("data", "external", "habitat", "Agriculture.tif"))%>%
  terra::project(globalclimpreds_terra[[1]])


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
eu_climpreds.10 <- terra::crop(globalclimpreds_terra, euboundary)
eu_climpreds.10 <- terra::mask(eu_climpreds.10, euboundary)%>%
  terra::crop(terra::ext(-38, 50,  24.29152732065, 72.66652712715))


#--------------------------------------------
#------ Load future climate rasters ----------
#--------------------------------------------

for (period in c("2041-2070","2071-2100")){
  for(scenario in c("ssp126", "ssp370", "ssp585")){
    
    # List future raster files
    future_files<- list.files( file.path("data", "external", "climate", "chelsa_future", period,scenario), pattern = "^scaled_layer.*\\.tif$", full.names = TRUE)
    
    #Stack them together
    future_stack <- terra::rast(future_files)
    
    #Aggregate at a resolution of 5km
    future_stack <- terra::aggregate(future_stack, fact=5, fun = mean, na.rm=TRUE)
    
    #Reproject future stack to match CRS and resolution of eu_climpreds.10
    future_aligned <- terra::project(
      future_stack,
      eu_climpreds.10,
      method = "bilinear")
    
    future_aligned <- future_aligned %>%
      terra::mask(eu_climpreds.10) %>%
      terra::crop(terra::ext(-38, 50,  24.29152732065, 72.66652712715))
    
    
    
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
wwf_eco_biome<-sf::st_read(file.path("./data/external/GIS/official/newRealms.shp")) 

# Optionally, make geometry valid
#wwf_eco_biome <- sf::st_make_valid(wwf_eco)


#--------------------------------------------
#-------Load file paths to bias grids -------
#--------------------------------------------
bias_grid_folder<-file.path("data","external", "bias_grids")
bias_grid_paths <- list(
  Plants = file.path(bias_grid_folder, "log_plants_1degree_layer.tif"), #0-13.24
  Amphibians = file.path(bias_grid_folder, "log_amphibians_1degree_layer.tif"),#0-12.06
  Birds = file.path(bias_grid_folder, "birds_1deg_min5.tif"),#5-1703018
  Mammals = file.path(bias_grid_folder, "log_mammals_1degree_layer.tif"), #0-13.36
  Molluscs = file.path(bias_grid_folder, "log_mollusca_1degree_layer.tif"),#0-12.48
  Reptiles = file.path(bias_grid_folder, "log_reptiles_1degree_layer.tif"), #0-11.34
  Fish = file.path(bias_grid_folder, "log_fish_1degree_layer.tif"),#0-14.67
  Malacostraca = file.path(bias_grid_folder, "log_malacostraca_1degree_layer.tif"),#0-13.12
  Insects = file.path(bias_grid_folder, "log_insects_1degree_layer.tif")) #0-15.78


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
    #----------- Load species details -----------
    #--------------------------------------------
    species <- names(split_df)[i]
    taxonkey<- unique(split_df[[i]]$acceptedTaxonKey)
    speciesName <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)  # Extract first two words of species name
    speciesgroup<-unique(split_df[[i]]$Group)
    
    
    #--------------------------------------------
    #----------- Load occurrence data -----------
    #--------------------------------------------
    global.occ.LL.cleaned<-split_df[[i]]%>%
      dplyr::select(c(decimalLongitude,decimalLatitude))
    global.occ_1KM<-cleaned_1km %>%
      dplyr::filter(acceptedTaxonKey == taxonkey)
    
    #Generate file for informing PA selection containing all occurrences (no thinning, in case we thinned split_df)
    for_PA_selection <- split_df_all_occs[[i]] %>%
      dplyr::select(c(decimalLongitude, decimalLatitude))%>%
      sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),crs = 4326)
    
    
    #---------------------------------------------
    #-- Prepare filenames and titles for export --
    #---------------------------------------------
    #Prepare PDF title 
    nameExtension<- if (grepl("^\\S+\\s+\\S+$", species)) "" else sub("^\\S+\\s+\\S+\\s+", "", species)
    PDF_title<-bquote(italic(.(gsub("_", " ", speciesName))) ~ .(nameExtension) ~ "(" * .(taxonkey) * ")")
    
    #Prepare current and future basefile
    basefile<-  paste0(speciesName,"_Climate_")
    
    
    #--------------------------------------------
    #-------------Create folders-----------------
    #--------------------------------------------
    # Define base project folder
    base_dir <- file.path("data", "projects", project, paste0(speciesName, "_", taxonkey))
    
    # Define outputs, periods, and scenarios
    periods   <- c("Current","2041-2070", "2071-2100")
    scenarios <- c("ssp126", "ssp370", "ssp585")
    outputs   <- c("Rasters", "PDFs", "PNGs")
    
    #Create folders for each combination
    scenario_folders <- list()
    
    for(period in periods){
      for(output in outputs){
        if(period=="Current"){
          loop_list <- list(list(path = file.path(base_dir, "Climate", period,"Predictions",output),
                                 name = paste("Climate", period, "Predictions",output,  sep = "/")),
                            list(path = file.path(base_dir, "Climate", period,"Diagnostics", "Variable_importance"),
                                 name = paste("Climate", period, "Diagnostics", "Variable_importance",  sep = "/")),
                            list(path = file.path(base_dir, "Climate", period,"Diagnostics", "Response_curves"),
                                 name = paste("Climate", period,"Diagnostics", "Response_curves", sep = "/")),
                            list(path = file.path(base_dir, "Climate", period,"Diagnostics", "Confidence_maps",output),
                                 name = paste("Climate", period, "Diagnostics", "Confidence_maps", output,  sep = "/")))
          scenario_folders <- c(scenario_folders, loop_list)  
          
        }else{
          for(scenario in scenarios){
            loop_list <- list(list(path = file.path(base_dir, "Climate", period, scenario, "Predictions", output),
                                   name = paste("Climate", period, scenario, "Predictions", output, sep = "/")),
                              list(path = file.path(base_dir, "Climate", period, scenario, "Diagnostics", "Confidence_maps",output),
                                   name = paste("Climate", period, scenario, "Diagnostics", "Confidence_maps", output,  sep = "/")))
            scenario_folders <- c(scenario_folders, loop_list)
          }
        }
      }
    }
    
    # Add Rasters/Interim folder
    fixed_folders <- list(
      list(path = file.path(base_dir, "Climate", "Current", "Interim"), 
           name = "Interim"))
    
    # Combine 
    folder_paths <- c(fixed_folders, scenario_folders)
    
    # Check and create each folder if necessary
    lapply(folder_paths, function(folder){
      create_folder(folder$path, folder$name)
    })
    
    
    #--------------------------------------------
    #------ Process occurrence data -----
    #--------------------------------------------
    #Keep only one occurrence per grid cell
    global.occ.LL.cleaned <- remove_duplicates(occurrences = global.occ.LL.cleaned, rast_template = globalclimpreds_terra[[1]])
    
    #Remove occurrences within grid cells with NA values
    global.occ.sf <- remove_nodata_occurrences(occurrences = global.occ.LL.cleaned, rast_template=globalclimpreds_terra[[1]], crs=4326)
    
    #add column indicating species presence (1) for modeling
    global.occ.sf$species <- rep(1, nrow(global.occ.sf)) 
    
    
    #-----------------------------------------------
    #------ Limit to 10,000 occupied grid cells ----
    #-----------------------------------------------
    if(nrow(global.occ.sf) > 10000){
      if(occurrence_thinning_method == "random"){
        print("Thinning occurrences randomly")
        set.seed(101) 
        global.occ.sf <- global.occ.sf[sample(nrow(global.occ.sf), 10000, replace=FALSE), ]
      }else if (occurrence_thinning_method == "kmeans_clustering"){
        
        print("Thinning occurrences based on k-means clustering")
        #Extract environmental data in each occurrence grid cell
        env_data <- terra::extract(globalclimpreds_terra, global.occ.sf, ID = FALSE)
        
        #Check how many unique rows there are and set centers to lowest of either 10000 or #unique rows
        unique_centers<-nrow(unique(env_data))
        center_number<-min(unique_centers, 10000)
        
        # K-means clustering
        set.seed(101)
        clust <- kmeans(env_data, centers = center_number,iter.max = 10, nstart = 1)$cluster
        occ_env<- cbind(global.occ.sf, env_data, clust)%>%
          dplyr::mutate(rID =row_number())
        
        # Keep 1 occurrence per cluster
        sampled <- occ_env %>%
          dplyr::group_by(clust) %>%
          dplyr::slice_sample(n = 1) %>%
          dplyr::ungroup()
        
        # How many presences do we still need
        remaining <- 10000 - nrow(sampled)
        
        # sample extra occurrences if fewer than 10000
        if (remaining > 0) {
          # Randomly sample additional presences excluding already chosen ones
          extra_occ <- occ_env %>%
            dplyr::filter(!rID %in% sampled$rID)%>%
            dplyr::slice_sample(n = remaining) 
          
          global.occ.sf <- bind_rows(sampled, extra_occ)
          rm(extra_occ)
          
        } else {
          global.occ.sf <- sampled
        }
        
        # Keep only occurrence columns
        global.occ.sf <- global.occ.sf %>%
          dplyr::select(decimalLongitude, decimalLatitude, geometry, species)
        
        rm(env_data, occ_env, sampled, remaining, unique_centers, center_number, clust)
        
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
    #-------------- Plot biomes -----------------
    #--------------------------------------------
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
    
    #---------Mask cells that contain occurrences---------
    # Convert sf to terra-compatible vector
    for_PA_vect <- terra::vect(for_PA_selection)
    # Identify raster cells that correspond to these points
    cells_to_exclude <- terra::cellFromXY(biasgrid_sub, terra::crds(for_PA_vect))
    # Set those cells to NA
    biasgrid_sub[cells_to_exclude] <- NA
    
    #--------------Generate pseudoabsences-----------------
    set.seed(728)
    global_points <- terra::spatSample(
      biasgrid_sub,
      size = 30000, #three times the number we need
      method = "weights",     # weighted random sampling
      as.points = TRUE,       # return SpatVector of points
      na.rm = TRUE            # ignore NA pixels
    )
    
    #Extract environmental data in each occurrence and presence absence grid cell
    occ_climate_data <- terra::extract(globalclimpreds_terra, global.occ.sf, ID = FALSE)
    pa_climate_data <- terra::extract(globalclimpreds_terra, global_points, ID = FALSE)
    
    # Find which combinations of climate values (rows) in pa_climate_df are already in occ_climate_df
    rows_to_keep <- !do.call(paste, pa_climate_data) %in% do.call(paste, occ_climate_data)
    
    #Filter points that have the exact same combo of environmental variables as some presence records
    global_points <- global_points[rows_to_keep, ]
    
    #Convert to sf dataframe
    global_points <- sf::st_as_sf(global_points) %>%
      dplyr::mutate(decimalLongitude = sf::st_coordinates(.)[,1],
                    decimalLatitude  = sf::st_coordinates(.)[,2]) %>%
      dplyr::select(-sum)  
    
    #If after filtering we have less than 10000 PA left, repeat with higher initial sampling
    if(nrow(global_points)< 10000) {
      set.seed(728)
      global_points <- terra::spatSample(
        biasgrid_sub,
        size = 50000, #Five times the number we need
        method = "weights",     # weighted random sampling
        as.points = TRUE,       # return SpatVector of points
        na.rm = TRUE            # ignore NA pixels
      )
      
      #Extract environmental data in each occurrence and presence absence grid cell
      pa_climate_data <- terra::extract(globalclimpreds_terra, global_points, ID = FALSE)
      
      # Find which combo of values in pa_climate_data are already in occ_climate_data
      rows_to_keep <- !do.call(paste, pa_climate_data) %in% do.call(paste, occ_climate_data)
      
      #Filter points that have the exact same combo of environmental variables as some presence records
      global_points <- global_points[rows_to_keep, ]
      
      #Convert to sf dataframe
      global_points <- sf::st_as_sf(global_points) %>%
        dplyr::mutate(decimalLongitude = sf::st_coordinates(.)[,1],
                      decimalLatitude  = sf::st_coordinates(.)[,2]) %>%
        dplyr::select(-sum)  
      
      #If still less than 10000 points, skip species 
      if(nrow(global_points) < 10000) {   
        warning(paste0(
          "Skipping species ", species, 
          " because too many pseudoabsences with a combination of environmental values ",
          "occurring in the presence dataset were selected."
        ))
        next  # Skip the rest of the loop and move to the next iteration
      }
    }
    
    #Select 10000 pseudoabsences
    if(nrow(global_points) > 10000){
      if(pseudoabsence_thinning_method == "random"){
        print("Thinning pseudoabsences randomly")
        set.seed(101) 
        global_points <- global_points[sample(nrow(global_points), 10000, replace=FALSE), ]
      }else if (pseudoabsence_thinning_method == "kmeans_clustering"){
        print("Thinning pseudoabsences based on k-means clustering")
        
        #Extract environmental data from filtered pseudoabsences
        pa_climate_data <- terra::extract(globalclimpreds_terra, global_points, ID = FALSE)
        
        #Check how many unique rows there are and set centers to lowest of either 10000 or #unique rows
        unique_centers<-nrow(unique(pa_climate_data))
        center_number<-min(unique_centers, 10000)
        
        # K-means clustering
        set.seed(101)
        clust <- kmeans(pa_climate_data, centers = center_number,iter.max = 10, nstart = 1)$cluster
        pa_climate <- cbind(global_points, pa_climate_data, clust)%>%
          dplyr::mutate(rID =row_number())
        
        # Keep 1 pseudoabsence per cluster
        sampled <- pa_climate %>%
          dplyr::group_by(clust) %>%
          dplyr::slice_sample(n = 1) %>%
          dplyr::ungroup()
        
        # How many pseudoabsences do we still need
        remaining <- 10000 - nrow(sampled)
        
        # sample extra pseudoabsences if fewer than 10000
        if (remaining > 0) {
          # Randomly sample additional pseudoabsences excluding already chosen ones
          extra_pa <- pa_climate %>%
            dplyr::filter(!rID %in% sampled$rID)%>%
            dplyr::slice_sample(n = remaining) 
          
          global_points <- bind_rows(sampled, extra_pa)
          rm(extra_pa)
          
        } else {
          global_points <- sampled
        }
        
        # Keep only occurrence columns
        global_points <- global_points %>%
          dplyr::select(decimalLongitude, decimalLatitude, geometry)
        
        rm(pa_climate_data, pa_climate, sampled, remaining, unique_centers, center_number, clust)
        
      }
    }
    
    
    #--------------------------------------------
    #--- Create presence-pseudoabsence dataset---
    #--------------------------------------------
    
    # Add coordinates and convert
    global_pseudoAbs <- global_points %>%
      dplyr::mutate(species = 0)
    
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
    global.data.df <- sdm::sdmData(species~.,train=vect(global_presabs),
                                   predictors = globalclimpreds_terra)%>%
      as.data.frame()
    
    
    #--------------------------------------------
    #--- Remove highly correlated predictors---
    #--------------------------------------------
    # Calculate correlation matrix (excluding rID and species)
    correlationMatrix <- cor(global.data.df[, -c(1, 2)])
    
    # Identify highly correlated variables (cutoff = 0.7)
    highlyCorrelated <- caret::findCorrelation(correlationMatrix, 
                                               cutoff = 0.7, 
                                               exact = TRUE, 
                                               names = TRUE)
    
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
        
        #Store
        modeloutput[[modelmethod]]<-list(fav_raster=fav_raster,
                                         model=method_model)
        
        rm(fav_raster, binary_1pct, binary_5pct, method_model)
      }
    })
    
    # List favourability rasters
    fav_rasters_list <- lapply(modeloutput, function(x) x$fav_raster)
    
    # Combine into a SpatRaster stack
    fav_stack <- terra::rast(fav_rasters_list)
    
    # Assign layer names based on model methods
    names(fav_stack) <- names(modeloutput)
    
    
    #---------------------------------------------
    #-- Create Ensemble model using PCAm method --
    #---------------------------------------------
    #Step 0: make PCA
    pca_result <- rasterPCA(fav_stack, nSamples = NULL, spca = FALSE, maskCheck = TRUE)
    
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
    cat("Top 5 models by variance along PC1:\n", top5_models)
    
    # Get model IDs
    top_ids <- info$modelID[info$method %in% top5_models]
    
    # Subset using those IDs
    top5models <- model[[top_ids]]  
    
    # Step 6: Subset fav_stack to top 5 layers
    top5_stack <- subset(fav_stack, top5_models)
    
    # Step 7: Compute pixel-wise median = consensus model
    consensus_median <- app(top5_stack, median)
    
    # Step 8: Compute pixel-wise mean for SD calculation
    consensus_mean <- mean(top5_stack, na.rm=TRUE)
    
    # Step 9: Compute pixel-wise population SD
    consensus_sd <- stdev(top5_stack, pop=TRUE)
    
    
    #------------------------------------------
    #-- Create map with ensemble suitability --
    #------------------------------------------
    #Define name of files
    base_file <- paste0(basefile, "current_ensemble")
    
    #Export PDFs with and without occurrences plotted
    for (occs in list(NULL, global.occ.sf)){
      filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))

      exportPDF(predictions = consensus_median,
                dataType = "Suit",
                scenario = "Current",
                returnPredictions = FALSE,
                returnPNG = FALSE,
                occ_data=occs,
                exportPNG=TRUE,
                PDF_title = PDF_title,
                PNG_folder=file.path(base_dir, "Climate", "Current", "Predictions", "PNGs"),
                PDF_folder=file.path(base_dir, "Climate", "Current", "Predictions","PDFs"),
                filename = filename)
    }
    
    
    #------------------------------------------
    #-- Create map with ensemble SD --
    #------------------------------------------
    #Define name of files
    filename <- paste0(basefile, "current_ensemble_SD")
    
    exportPDF(predictions = consensus_sd,
              dataType = "Stdev",
              scenario = "Current",
              returnPredictions = FALSE,
              returnPNG = FALSE,
              occ_data=NULL,
              exportPNG=TRUE,
              PDF_title = PDF_title,
              PNG_folder=file.path(base_dir, "Climate", "Current", "Diagnostics", "Confidence_maps", "PNGs"),
              PDF_folder=file.path(base_dir, "Climate", "Current", "Diagnostics", "Confidence_maps","PDFs"),
              filename = filename)
    
    
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
    
    #Create one df with the median favorability value for each occurrence
    fav_vals <- fav_vals %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(median = apply(., 1, median, na.rm = TRUE)) %>% #1 = apply to rows
      dplyr::select(median)
    
    # Create binary maps
    binary_maps<-list()
    for (probs in mtp_probabilities){
      
      #Define mtp_pct and mtp_value
      mtp_value <- probs*100
      mtp_pct <- paste0(mtp_value, "%")
      
      # Obtain threshold
      to_omit <- floor(probs * nrow(fav_vals)) #Define how many of lowest ranked occs to omit based on mtp threshold
      thr <- sort(fav_vals$median)[to_omit + 1]
      cat(paste0("Mean ",mtp_pct," minimum training presence threshold: ", round(thr, 4), "\n"))
      
      # Create binary raster using MTP threshold
      binary_map_pct <- consensus_median >= thr 
      binary_map_pct <- as.factor( binary_map_pct*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
      levels( binary_map_pct) <- data.frame(ID = c(0, 1),
                                            class = c("Absent", "Present"))
      
      # Calculate sensitivity in Europe
      occ_values <- terra::extract(binary_map_pct, vect(global.occ.sf))[,2]  
      global_EU_sensitivity <- sum(occ_values == "Present", na.rm = TRUE) / sum(occ_values %in% c("Present", "Absent"), na.rm = TRUE)
      
      #Store raster
      raster_folder <- file.path(base_dir, "Climate","Current", "Predictions", "Rasters")
      binary_file <- file.path (raster_folder, paste0(basefile,"current_binary",mtp_value,"pct.tif"))
      terra::writeRaster(binary_map_pct, filename = binary_file, overwrite = TRUE)
      
      # export as PDF and PNG with and without occurrences plotted 
      base_file <- paste0(basefile, "current_binary",mtp_value,"pct")
      for (occs in list(NULL, global.occ.sf)){
        filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
        exportPDF(predictions = binary_map_pct,
                  dataType = "Binary",
                  scenario = "Current",
                  returnPredictions = FALSE,
                  returnPNG = TRUE,
                  occ_data=occs,
                  exportPNG=TRUE,
                  LabelValue= round(thr,3),
                  LabelName=paste0(mtp_pct, " MTP threshold"),
                  Label2Value=round(global_EU_sensitivity,3),
                  Label2Name="Sensitivity",
                  PDF_title = PDF_title,
                  PNG_folder=file.path(base_dir, "Climate","Current", "Predictions", "PNGs"),
                  PDF_folder=file.path(base_dir,"Climate" ,"Current", "Predictions", "PDFs"),
                  filename = filename)
      }
      
      assign(paste0(mtp_value,"pct"), thr)
      binary_maps[[mtp_pct]] <- list(binary_raster=binary_map_pct,
                                     EU_sensitivity=global_EU_sensitivity,
                                     mean_MTP= thr)
      rm(binary_map_pct, binary_file)
    }
    
    
    #--------------------------------------------
    #--- Create maps with future projections ----
    #--------------------------------------------
    for (period in c("2041-2070","2071-2100")){
      for(scenario in c("ssp126", "ssp370", "ssp585")){
        
        print(paste("[FUTURE] Projecting:", period,scenario))
        
        #Get climate data for specific period and scenario
        future_rast <- get(paste0(scenario, "_", period))
        
        # Keep relevant predictors in the raster stack
        eu_future_selection <- future_rast %>%
          subset(names(eu_climpreds.10_selection))
        
        # Project each of the top 5 models
        future_modeloutput <- list()
        for(modelmethod in top5_models){
          pred_raster_future <- predict(model,
                                        newdata = eu_future_selection,
                                        method = modelmethod)
          fav_raster_future <- favourability_from_prob(pred_raster_future, prev_ratio)
          future_modeloutput[[modelmethod]] <- fav_raster_future
          rm(fav_raster_future, pred_raster_future)
        }
        
        # Create Ensemble predictions for future
        future_fav_stack <- terra::rast(future_modeloutput)
        future_consensus_median <- app(future_fav_stack, median)
        future_consensus_mean <- mean(future_fav_stack, na.rm=TRUE)
        future_consensus_sd <- stdev(future_fav_stack, pop=TRUE)
        
        # Export future ensemble raster (favorability) 
        future_folder <- file.path(base_dir, "Climate", period, scenario, "Predictions", "Rasters")
        ensemble_file <- file.path(future_folder, paste0(basefile, period,"_",scenario,"_ensemble.tif"))
        terra::writeRaster(future_consensus_median, filename = ensemble_file, overwrite = TRUE)
        
        # Export future sd raster 
        future_sd_folder <- file.path(base_dir, "Climate", period, scenario, "Diagnostics", "Confidence_maps", "Rasters")
        ensemble_sd_file <- file.path(future_sd_folder, paste0(basefile, period,"_",scenario,"_ensemble_SD.tif"))
        terra::writeRaster(future_consensus_sd, filename = ensemble_sd_file, overwrite = TRUE)
        
        # Export future mean raster 
        future_mean_folder <- file.path(base_dir, "Climate", "Current", "Interim")
        ensemble_mean_file <- file.path(future_mean_folder, paste0(basefile, period,"_",scenario,"_ensemble_mean.tif"))
        terra::writeRaster(future_consensus_mean, filename = ensemble_mean_file, overwrite = TRUE)
        
        # # Export future single-model rasters
        # for (mod in top5_models) {
        #   singlemodfile <- file.path(future_folder,
        #                               paste0(basefile, period, "_",scenario,"_", mod, ".tif"))
        #   terra::writeRaster(future_fav_stack[[mod]], filename = singlemodfile, overwrite = TRUE)
        # }
        
        # Export ensemble predictions as PDF and PNG with and without occurrences
        base_file <- paste0(basefile, scenario,"_", period,"_ensemble")
        
        for (occs in list(NULL, global.occ.sf)){
          filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))

          exportPDF(predictions = future_consensus_median,
                    dataType = "Suit",
                    period = period,
                    scenario = scenario,
                    returnPredictions = FALSE,
                    returnPNG = TRUE,
                    occ_data=occs,
                    exportPNG=TRUE,
                    PDF_title=PDF_title,
                    PNG_folder=file.path(base_dir, "Climate", period, scenario, "Predictions", "PNGs"),
                    PDF_folder=file.path(base_dir, "Climate", period, scenario, "Predictions", "PDFs"),
                    filename = filename)
        }
        
        # Export ensemble SD predictions as PDF and PNG 
        filename<- paste0(basefile, scenario,"_", period,"_ensemble_SD")
        
        exportPDF(predictions = future_consensus_sd,
                  dataType = "Stdev",
                  period = period,
                  scenario = scenario,
                  returnPredictions = FALSE,
                  returnPNG = TRUE,
                  occ_data=NULL,
                  exportPNG=TRUE,
                  PDF_title=PDF_title,
                  PNG_folder=file.path(base_dir, "Climate", period, scenario, "Diagnostics", "Confidence_maps", "PNGs"),
                  PDF_folder=file.path(base_dir, "Climate", period, scenario, "Diagnostics", "Confidence_maps", "PDFs"),
                  filename = filename)

        
        # Create binarized ensemble predictions for future
        for(probs in mtp_probabilities){
          #Define mtp_pct and mtp_value
          mtp_label <- paste0(probs*100, "%")
          mtp_text <- paste0(probs*100,"pct")
          
          #Get threshold value and apply to consensus predictions
          threshold<-get(mtp_text)
          binary_map_future <- future_consensus_median  >= threshold
          binary_map_future <- as.factor( binary_map_future*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
          levels( binary_map_future) <- data.frame(ID = c(0, 1),
                                                   class = c("Absent", "Present"))
          
          #Store raster
          binary_file <- file.path(future_folder, paste0(basefile, period,"_",scenario,"_binary",mtp_text,".tif"))
          terra::writeRaster(binary_map_future, filename = binary_file, overwrite = TRUE)
          
          # Export binarized ensemble predictions as PDF and PNG with and without occurrences 
          base_file <- paste0(basefile, period,"_", scenario, "_binary",mtp_text)
          
          for (occs in list(NULL, global.occ.sf)){

            filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
            exportPDF(predictions = binary_map_future,
                      dataType = "Binary",
                      period = period,
                      scenario = scenario,
                      occ_data=occs,
                      returnPredictions = FALSE,
                      returnPNG = TRUE,
                      exportPNG=TRUE,
                      LabelValue= round(threshold,3),
                      LabelName=paste0(mtp_label, " MTP threshold"),
                      PDF_title=PDF_title,
                      PNG_folder=file.path(base_dir,"Climate", period, scenario, "Predictions", "PNGs"),
                      PDF_folder=file.path(base_dir, "Climate",period, scenario, "Predictions", "PDFs"),
                      filename=filename)
          }
        }
      }
    }
    
    
    #-----------------------------------------------
    #- Get response curves and variable importance -
    #-----------------------------------------------
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
    response_df <- purrr::imap_dfr(response_list, function(model_list, model_name) {
      purrr::imap_dfr(model_list, function(df, var_name) {
        response_df <- df %>%
          setNames(c("Predictor_value", "Response"))%>%
          mutate( Algorithm = model_name,
                  Predictor = var_name)})
    }) %>%
      dplyr::select(Algorithm,Predictor, Predictor_value, Response)
    
    
    varimp_df <- purrr::imap_dfr(varimp_list, function(df, model_name) {
      df %>%
        setNames(c("Predictor", "corTest" , "AUCtest"))%>%
        dplyr::mutate(Algorithm = model_name)
    })%>%
      dplyr::select(Algorithm,Predictor, corTest, AUCtest)
    
    
    # Plot response curves
    response_plot <- ggplot(response_df, aes(x = Predictor_value,
                                             y = Response, 
                                             color = Algorithm)) +
      geom_line(linewidth=0.8) +
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
    
    #Save plot
    PNG_folder <- file.path(base_dir, "Climate", "Current", "Diagnostics")
    
    ggplot2::ggsave(filename = paste(basefile, "variable_importance.png"), plot = varimp_plot ,  device = "png", width =8.27 , height = 5.845, path= file.path(PNG_folder, "Variable_importance") )
    ggplot2::ggsave(filename = paste(basefile, "response_curves.png"), plot = response_plot,  device = "png", width =8.27 , height = 5.845, path=  file.path(PNG_folder, "Response_curves") )
    
    
    
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
    climatemodel <-list(species = species, #Species name
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
                        global_EU_sensitivity_5pct_threshold =  binary_maps$`5%`$EU_sensitivity,
                        global_EU_sensitivity_1pct_threshold =  binary_maps$`1%`$EU_sensitivity,
                        response_df = response_df,
                        varimp_df = varimp_df,
                        top5models = top5models,
                        selected_predictors = names(eu_future_selection),
                        future_consensus_median = future_consensus_median)
    
    qs::qsave(climatemodel, file.path(base_dir, "Climate", paste0("Climate_model_",speciesName,"_",taxonkey,".qs")))
    
    
    #--------------------------------------------
    #--Export raster layers in folder "Rasters"--
    #--------------------------------------------
    #We don't store them in .qs file as some important metadata would be stored in a temp folder, which would be removed after a while 
    biasgrid_file<- file.path(base_dir,"Climate", "Current", "Interim", 
                              paste0("Biasgrid_",speciesName,"_",taxonkey,".tif"))
    ensemble_median_file <- file.path( base_dir,"Climate", "Current", "Predictions", "Rasters",
                                       paste0(basefile, "current_ensemble.tif"))
    ensemble_mean_file <- file.path( base_dir,"Climate", "Current", "Interim",
                                     paste0(basefile, "current_ensemble_mean.tif"))
    ensemble_sd_file <- file.path( base_dir,"Climate", "Current","Diagnostics", "Confidence_maps", "Rasters",
                                   paste0(basefile, "current_ensemble_SD.tif"))
    
    terra::writeRaster(biasgrid_sub, filename = biasgrid_file, overwrite = TRUE)
    terra::writeRaster(consensus_median, filename = ensemble_median_file, overwrite = TRUE)
    terra::writeRaster(consensus_mean, filename = ensemble_mean_file, overwrite = TRUE)
    terra::writeRaster(consensus_sd, filename = ensemble_sd_file, overwrite = TRUE)
    
    
    #--------------------------------------------
    #------------------ Clean up-----------------
    #--------------------------------------------
    rm(list = setdiff(ls(), c("p","wwf_eco_biome","eu_climpreds.10", "split_df",  "decimalplaces", "globalclimpreds_terra","bias_grid_paths", "i", "world", "project", "create_folder", "split_df_all_occs", "exportPDF", "remove_duplicates", "remove_nodata_occurrences", "favourability_from_prob", "cleaned_1km", "occurrence_thinning_method", "n_clusters", "ssp126_2041-2070","ssp370_2041-2070","ssp585_2041-2070","ssp126_2071-2100","ssp370_2071-2100","ssp585_2071-2100", "mtp_probabilities", "pseudoabsence_thinning_method")))
    
  }
})


#--------------------------------------------
#---------- Clean R environment--------------
#--------------------------------------------
rm(list = ls())

