#--------------------------------------------
#-----------To do: specify project ----------
#--------------------------------------------
#specify project name
projectname <- "PA prob & Alternative Treshold & Ensemble Boyce"


#--------------------------------------------
#-------------- Load packages ---------------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")

packages <- c( "dplyr", "stringr", "here", "qs","CoordinateCleaner","terra", "raster", "rnaturalearth", "rnaturalearthdata", 
               "ggplot2","tidyterra","mapview", "dismo", "sdm", "caret", "viridisLite", "kableExtra","future", "future.apply",
               "randomForest","earth", "progressr", "sf", "gbm", "PresenceAbsence","geosphere","arm")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#--- Load right version of caretEnsemble-----
#--------------------------------------------
# Define the required version
desired_version <- "2.0.3"

if ("caretEnsemble" %in% rownames(installed.packages())) {
  current_version <- packageVersion("caretEnsemble")
  if (as.character(current_version) != desired_version) {
    remove.packages("caretEnsemble")
    devtools::install_github("zachmayer/caretEnsemble@2.0.3")
    library(caretEnsemble)
    rm(current_version, desired_version)
  } else {
    library(caretEnsemble)
    rm(current_version, desired_version)
  }
  
} else {
  devtools::install_github("zachmayer/caretEnsemble@2.0.3")
  library(caretEnsemble)
  rm(desired_version)
}


#--------------------------------------------
#-------- Source helper functions------------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#--- Load global occurrences and taxa info---
#--------------------------------------------
global<-qs::qread( paste0("./data/projects/",projectname,"/",projectname,"_occurrences.qs"))
taxa_info<-read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  dplyr::pull(speciesKey)%>%
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
#Reasoning: several invasive species that have datasets from designs with default uncertainty of 5km
global.occ<-global %>%
  dplyr::filter(speciesKey%in%accepted_taxonkeys) %>%   
  dplyr::filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters<= 5000) %>% 
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
  dplyr::select(decimalLongitude, decimalLatitude, species, speciesKey, Group) 
rm(global.occ, global)


#--------------------------------------------
#-----------Do coordinate cleaning-----------
#--------------------------------------------
# OPTIONAL: Coordinates are tested for several things: whether they are in capitals, whether ... . For each coordinate a column per test is added indicating wether the result is potentially problematic (FALSE) or a clean coordinate (TRUE)
#flags_report<-clean_coordinates(x = global.occ.LL, lon= "decimalLongitude", lat= "decimalLatitude",
#  tests = c("capitals", 
# "centroids","gbif", "institutions", 
# "seas", "zeros"))

# Clean coordinates based on their proximity to country centroids, capitals, biodiversity institutions, GBIF headquarters, and the 0/0 point
cleaned<-global.occ.LL%>%
  CoordinateCleaner::cc_cen(buffer=100) %>% # remove points within a buffer of 100m around country centroids, default 1km
  CoordinateCleaner::cc_cap(buffer=100) %>% # remove capitals centroids (buffer 100m), default 10km
  CoordinateCleaner::cc_inst(buffer=100) %>% # remove zoo and herbaria records buffer of 100 m around biodiversity institutes, default 100m
  CoordinateCleaner::cc_gbif(buffer=100)%>% #remove around GBIF headquarters in Copenhagen (buffer 100m), default 100m
  CoordinateCleaner::cc_zero() #Remove around the 0/0 point (buffer 0.5 degrees)


#--------------------------------------------
#--------Load global climate rasters --------
#--------------------------------------------
# Only include files that start with "scaled_layer_" and end with .tif: TODO: NEED TO CREATE THEM FIRST!
scaled_files <- list.files(
  here::here("data/external/climate/scaled_layers"),
  pattern = "^scaled_layer_\\d+\\.tif$",
  full.names = TRUE)

# Load and stack
globalclimpreds_terra <- terra::rast(scaled_files)

# Optional: check
print(globalclimpreds_terra)


#--------------------------------------------
#--------Load European climate rasters-------
#--------------------------------------------
euboundary<-sf::st_read(here("./data/external/GIS/Europe/EUROPE.shp")) 
# Convert sf boundary to SpatVector
euboundary_vect <- terra::vect(euboundary)

# Crop and mask scaled_stack
scaled_stack_europe <- terra::crop(globalclimpreds_terra, euboundary_vect)
eu_climpreds.10 <- terra::mask(scaled_stack_europe, euboundary_vect)


#---------------------------------------------
#-- Remove NA pixels from climate rasters ----
#---------------------------------------------
# This ensures all rasters have the same NA structure
na_mask_globalclimpreds_terra <- terra::sum(is.na(globalclimpreds_terra)) > 0
globalclimpreds_terra <- terra::mask(globalclimpreds_terra, na_mask_globalclimpreds_terra, maskvalue = 1)
na_mask_eu_climpreds.10 <- terra::sum(is.na(eu_climpreds.10)) > 0
eu_climpreds.10 <- terra::mask(eu_climpreds.10, na_mask_eu_climpreds.10, maskvalue = 1)


#--------------------------------------------
#--------- Load shape of the world ----------
#--------------------------------------------
world<-rnaturalearth::ne_countries(scale=50)


#--------------------------------------------
#--------------Load ecoregions --------------
#--------------------------------------------
wwf_eco_biome<-sf::st_read(here::here("./data/external/GIS/official/newRealms.shp")) #TODO load file
#wwf_eco_biome<-sf::st_transform(wwf_eco, 4326) %>%
#  sf::st_make_valid()
plot(wwf_eco_biome)


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
bias_grid_paths <- list(
  Plants = here::here("./data/external/bias_grids/final/trias/plants_10km_bias_layer_log.tif"), #0-13.24
  Amphibians = here::here("./data/external/bias_grids/final/trias/amphibians_10km_bias_layer_log.tif"),#0-12.06
  Birds = here::here("./data/external/bias_grids/final/trias/birds_1deg_min5.tif"),#5-1703018
  Mammals = here::here("./data/external/bias_grids/final/trias/mammals_10km_bias_layer_log.tif"), #0-13.36
  Molluscs = here::here("./data/external/bias_grids/final/trias/mollusca_10km_bias_layer_log.tif"),#0-12.48
  Reptiles = here::here("./data/external/bias_grids/final/trias/reptiles_10km_bias_layer_log.tif"), #0-11.34
  Fish = here::here("./data/external/bias_grids/final/trias/fish_10km_bias_layer_log.tif"),#0-14.67
  Malacostraca = here::here("./data/external/bias_grids/final/trias/malacostraca_10km_bias_layer_log.tif"),#0-13.12
  Insects = here::here("./data/external/bias_grids/final/trias/insects_10km_bias_layer_log.tif")) #0-15.78


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
    species<-names(split_df)[i]
    print(species)
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)  # Extract first two words of species name
    global.occ.LL.cleaned<-split_df[[i]]
    taxonkey<-unique(global.occ.LL.cleaned$speciesKey)
    speciesgroup<-unique(global.occ.LL.cleaned$Group)
    global.occ.LL.cleaned<-global.occ.LL.cleaned %>%
      dplyr::select(c(decimalLongitude,decimalLatitude))
    
    #replicate code above to generate file for informing PA selection - no PA here
    sp <- names(split_df_all_occs)[i]
    sp_short <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", sp)  # Extract first two words
    sp_data <- split_df_all_occs[[i]]
    sp_key <- unique(sp_data$speciesKey)
    sp_group <- unique(sp_data$Group)
    sp_data <- sp_data %>%
      dplyr::select(c(decimalLongitude, decimalLatitude))
    
    for_PA_selection <- sf::st_as_sf(global.occ.LL.cleaned,
                                     coords = c("decimalLongitude", "decimalLatitude"),
                                     crs = 4326)
    
    
    #--------------------------------------------
    #-------------Create folders-----------------
    #--------------------------------------------
    # Define the folder paths
    folder_paths<-list(list("path"=file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Rasters", "Interim"),
                            "name"= "Rasters/Interim"),
                       list("path"=file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Rasters", "Global"),
                            "name"= "Rasters/Global"),
                       list("path"=file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PDFs"),
                            "name"= "PDFs"),
                       list("path"=file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PNGs"),
                            "name"= "PNGs"))
    
    # Check and create each folder if necessary
    lapply(folder_paths, function(folder){
      create_folder(folder$path, folder$name)
    })
    
    
    #--------------------------------------------
    #------ Remove duplicates per grid cell -----
    #--------------------------------------------
    global.occ.LL.cleaned$cell<-terra::cellFromXY( globalclimpreds_terra, global.occ.LL.cleaned) #Indicate for each occurrence point in which cell of the raster it falls
    global.occ.LL.cleaned<-global.occ.LL.cleaned[!is.na(global.occ.LL.cleaned$cell),]
    unique_occurrences <- !duplicated(global.occ.LL.cleaned$cell)# Identify unique occurrences
    global.occ.LL.cleaned <- global.occ.LL.cleaned[unique_occurrences, 1:2] # Subset the occurrence points to keep only one occurrence per raster cell 
    
    global.occ.LL.cleaned<- terra::extract(globalclimpreds_terra, global.occ.LL.cleaned, xy = T, ID=F)%>%
      dplyr::filter(rowSums(is.na(.[, 1:(ncol(.) - 2)])) == 0)%>% #Keep rows that do not have any NA values in column 1- 3rd last 
      dplyr::select(c(x,y))%>%
      dplyr::rename(decimalLongitude=x,
                    decimalLatitude=y) #Extract climatic values of occurrence points from each raster layer and remove occurrence points that fall in cells with NA values in at least one rasterlayer
    
    #Convert to sf dataframe
    global.occ.LL.cleaned$species<- rep(1,length(global.occ.LL.cleaned$decimalLongitude)) #adds columns indicating species presence (1) needed for modeling
    global.occ.sf<-sf::st_as_sf(global.occ.LL.cleaned, coords=c("decimalLongitude", "decimalLatitude"), crs=4326, remove= FALSE)
    
    
    #--------------------------------------------
    #------ Define number of pseudoabsences -----
    #--------------------------------------------
    numb.global.pseudoabs <-length(global.occ.sf$decimalLongitude) #sets the number of pseudoabsences equal to number of unique presences
    rm(global.occ.LL.cleaned)
    
    
    #-------------------------------------------------------
    #-Don't fit model if less than 20 global presences -----
    #-------------------------------------------------------   
    if(numb.global.pseudoabs<20){
      warning(paste0("Skipping species ", species, " because the number of occurrences is less than 20 (n =",numb.global.pseudoabs,")"))
      next  # Skip the rest of the loop and move to the next iteration
    }
    
    
    #--------------------------------------------
    #------ Plot distribution of occurrences ----
    #--------------------------------------------
    #ggplot()+ 
    #geom_sf(data = world,  colour = "black", fill = NA)+
    #geom_point(data=global.occ.sf, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
    #labs(x="Longitude", y="Latitude")+
    #theme_bw()
    
    
    #--------------------------------------------
    #- Select ecoregions containing occurrences -
    #--------------------------------------------
    
    # Ensure valid geometries
    #wwf_eco_biome <- sf::st_make_valid(wwf_eco_biome)
    
    # Disable S2 geometry engine to avoid topological issues
    #sf::sf_use_s2(FALSE)
    
    # Keep only biome polygons that intersect at least one occurrence point
    has_occurrence <- lengths(sf::st_intersects(wwf_eco_biome, global.occ.sf)) > 0
    wwf_ecoSub1 <- wwf_eco_biome[has_occurrence, ]
    
    # Plot result
    #plot(wwf_ecoSub1$geometry)
    #points(global.occ.sf, col = "red", pch = 16, cex = 0.5)
    
    
    #--------------------------------------------
    #------------- Plot ecoregions --------------
    #--------------------------------------------
    #ggplot()+ 
    #geom_sf(data = world,  colour = "black", fill = NA)+
    #geom_sf(data=wwf_ecoSub1, fill="#f7786f")+
    #labs(x="Longitude", y="Latitude")+
    #theme_bw()
    
    
    #--------------------------------------------
    #------ Import right bias grid --------------
    #--------------------------------------------
    if (speciesgroup %in% names(bias_grid_paths)) {
      biasgrid <- terra::rast(bias_grid_paths[[speciesgroup]])
      if(speciesgroup %in% c("Amphibians", "Molluscs", "Mammals", "Reptiles","Birds","Plants","Fish","Malacostraca","Insects")){
        # Resample biasgrid to match the resolution of globalclimpreds_terra
        biasgrid<- terra::resample(biasgrid, globalclimpreds_terra[[1]], method="bilinear")
      }
    } else {
      stop("No bias grid available for this species. Species has to be one of the following: Amphibians, Molluscs, Mammals, Reptiles, Birds, Plants, Fish, Malacostraca, or Insects.")
    }
    
    #SET NAs
    #savegr <- biasgrid
    biasgrid_masked <- terra::mask(biasgrid, globalclimpreds_terra[[1]])
    biasgrid_log <- log(biasgrid_masked+1)
    plot(biasgrid_log)
    
    
    #rescale between 1 and 20
    min_val <- global(biasgrid_log, fun = "min", na.rm = TRUE)[[1]]
    max_val <- global(biasgrid_log, fun = "max", na.rm = TRUE)[[1]]
    print(min_val)
    print(max_val)
    # Rescale raster values to range from 1 to 20
    biasgrid <- ((biasgrid_log - min_val) / (max_val - min_val)) * 19 + 1
    plot(biasgrid_20)
    summary(biasgrid_20)
    
    
    # Define the dynamic output path
    biasgrid_file <- file.path("./data/projects",projectname, paste0(first_two_words, "_", taxonkey),"Rasters", "Interim", paste0("Biasgrid_", first_two_words, "_", taxonkey, ".tif"))
    
    # Ensure the output folder exists
    dir.create(dirname(biasgrid_file), recursive = TRUE, showWarnings = FALSE)
    
    # Save the biasgrid raster
    terra::writeRaster(biasgrid, filename = biasgrid_file, overwrite = TRUE)
    
    # Optional: confirm writing
    cat("✔ Biasgrid saved for", species, "→", biasgrid_file, "\n")
    
    
    #--------------------------------------------
    #Mask biasgrid by ecoregions with occurrences 
    #--------------------------------------------
    wwf_ecoSub1_ext<-terra::ext(wwf_ecoSub1) 
    wwf_ecoSub1_vector <- vect(wwf_ecoSub1) #Convert wwf_ecoSub1 to a SpatVector that can be used for masking
    biasgrid_crop <- terra::crop(biasgrid, wwf_ecoSub1_ext) #Crop biasgrid to extent wwf_ecoSub1
    biasgrid_sub <- terra::mask(biasgrid_crop, wwf_ecoSub1_vector)#Mask cropped biasgrid with SpatVector
    
    #Mask biasgrid with one of the climatic layers, to make sure it doesn't extend beyond them
    climategrid_rast<-terra::crop(globalclimpreds_terra[[1]], wwf_ecoSub1_ext)
    biasgrid_sub<-terra::mask(biasgrid_sub, climategrid_rast) 
    
    #biasgrid_sub_raster <- raster::raster(biasgrid_sub) #Convert SpatRaster back to normal raster object, needed in function generate_pseudoabs(which is based on dismo::randomPoints)
    
    
    #--------------------------------------------
    #---------------Visualize biasgrid-----------
    #--------------------------------------------
    #ggplot()+ 
    #geom_sf(data = world,  colour = "black", fill = NA)+
    #geom_spatraster(data=biasgrid_sub)+
    #scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange")+
    #labs(x="Longitude", y="Latitude")+
    #theme_bw()
    
    
    #--------------------------------------------
    #---------- Generate pseudoabsences ---------
    #--------------------------------------------
    # generates 10.000 pseudo absences, from cells that are not NA in biasgrid_sub and not from cells that have occurrence points as indicated in global.occ.sf. 
    #extf 1.1 increases the size of extent with 5% at each side of the extent (default value) but when ext is NULL it won't do anything.
    #excludep = TRUE indicates that presence points are excluded from the background, prob: if TRUE the values in mask are interpreted as probability weights (only works for rasters with a modest size that can be loaded into RAM)
    
    #Create alternative raster consisting of only ecoregions without biasgrid mask, used only when not enough pseudoabsence points can be generated using biasgrid_sub as mask layer
    ecoregions_crop<-terra::crop(globalclimpreds_terra[[1]], wwf_ecoSub1_ext) #Crop one of the climate rasters to extent ecoregions
    ecoregions_raster<-terra::mask(ecoregions_crop,wwf_ecoSub1_vector) #Mask with ecoregions vector
    
    #LIMIT TO AREA OF 100KM AROUND OCCURRENCES
    #global.occ.sf_buffer_100km <- st_buffer(global.occ.sf, dist = 100000)
    #global.occ.sf_buffer_100km <- vect(st_union(global.occ.sf_buffer_100km))
    # Convert raster from 'raster' to 'terra' format if needed
    #if (inherits(biasgrid_sub_raster, "Raster")) {
    #  biasgrid_sub_raster <- rast(biasgrid_sub_raster)
    #}
    #biasgrid_sub_raster<-terra::mask(biasgrid_sub_raster,global.occ.sf_buffer_100km)
    #biasgrid_sub_raster_raster <- raster::raster(biasgrid_sub_raster)
    
    
    #MASK OUR CELLS THAT CONTAIN OCCURRENCES
    # Convert sf to terra-compatible vector
    for_PA_vect <- vect(for_PA_selection)
    # Identify raster cells that correspond to these points
    cells_to_exclude <- cellFromXY(biasgrid_sub, crds(for_PA_vect))
    # Set those cells to NA
    biasgrid_sub[cells_to_exclude] <- NA
    
    #Generate pseudoabsences
    set.seed(728)
    global_points <- spatSample(
      biasgrid_sub,
      size = 10000,
      method = "weights",     # weighted random sampling
      as.points = TRUE,       # return SpatVector of points
      na.rm = TRUE # ignore NA pixels
    )
    
    
    #--------------------------------------------
    #--- Create presence-pseudoabsence dataset---
    #--------------------------------------------
    # Extract coordinates from SpatVector
    coords <- terra::crds(global_points)
    
    # Add coordinates and convert
    global_pseudoAbs <- global_points %>%
      as.data.frame() %>%
      mutate(x = coords[, 1],
             y = coords[, 2],
             species = 0) %>%
      sf::st_as_sf(coords = c("x", "y"), crs = 4326, remove = FALSE) %>%
      rename(decimalLongitude = x,
             decimalLatitude = y)
    
    # Now combine
    # Optional: remove extra column if present in pseudoabsence data
    global_pseudoAbs_clean <- global_pseudoAbs %>%
      dplyr::select(decimalLongitude, decimalLatitude, species, geometry)
    
    
    global_presabs <- rbind(global.occ.sf, global_pseudoAbs_clean)
    #rm(global_points)
    
    
    #--------------------------------------------
    #--Visualize presence-pseudoabsence dataset--
    #--------------------------------------------
    #mapview(biasgrid_sub_raster, 
    #col.regions = colorRampPalette(c("blue", "orange")),
    #alpha=1, 
    #  na.color = "transparent", 
    #layer.name = "Bias Grid") +
    #mapview(global_presabs, zcol = "species", 
    #col.regions = c("red", "yellow"),
    #layer.name = "Species distribution")
    
    
    #--------------------------------------------
    #---- Extract climate data for modelling-----
    #--------------------------------------------
    global.data <- sdm::sdmData(species~.,train=vect(global_presabs),predictors=globalclimpreds_terra) 
    global.data.df<-as.data.frame(global.data)
    
    
    #--------------------------------------------
    #--- Remove highly correlated predictors---
    #--------------------------------------------
    # Identify highly correlated predictors
    correlationMatrix<-cor(global.data.df[,-c(1,2)]) #Calculate pearson correlation among environmental values
    highlyCorrelated <- caret::findCorrelation(correlationMatrix, cutoff=0.7,exact=TRUE,names=TRUE)#Returns names of environmental variables to be removed because they are correlated more than 0.7 with other variables.  If two variables have a high correlation, the function removes the variable with the largest mean absolute correlation.
    preds<-as.data.frame(highlyCorrelated)
    
    # Remove highly correlated predictors from dataframe 
    global.data.df.subset<- global.data.df %>%
      dplyr::select (-all_of(highlyCorrelated), -rID) %>% 
      dplyr::mutate(species = as.factor(species)) %>%
      dplyr::mutate(species = recode_factor(species, 
                                            '0' = "absent",
                                            '1' = "present")) %>%  # Later steps require non numeric dependent variable
      dplyr::mutate(species = relevel(species, ref = "present")) 
    
    #Remove them from climate stack
    eu_climpreds.10_selection<-eu_climpreds.10%>%
      subset(!names(eu_climpreds.10) %in% highlyCorrelated)
    
    
    #--------------------------------------------
    #--Correct climate data from integer format--
    #--------------------------------------------
    #Divide all climate data by 10
    global.data.df.uncor<-cbind("species"=  global.data.df.subset$species,divide10(global.data.df.subset[,-c(1)]))
    
    
    #--------------------------------------------
    #--- Run multiple machine learning models ---
    #--------------------------------------------
    #preProc: preprocessing of predictors (environmental data). method = "center" subtracts the mean of the predictor's data (again from the data in x) from the predictor values while method = "scale" divides by the standard deviation.
    control <- caret::trainControl(method="cv",
                                   number=10,
                                   savePredictions="final", 
                                   preProc=c("center","scale"),
                                   classProbs=TRUE)
    classList1 <- c("glm","gbm","rf","earth")
    set.seed(457)
    global_train <- caretEnsemble::caretList(species~., 
                                             data= global.data.df.uncor,
                                             trControl=control,
                                             methodList=classList1,
                                             metric="Accuracy")
    
    
    #--------------------------------------------
    #--Return accurracy, kappa and correlation --
    #--------------------------------------------
    #Return the results for each model
    GlobalModelResults<-caret::resamples(global_train)  
    
    # Display accuracy of each model
    Global.Mod.Accuracy<-summary(GlobalModelResults)
    
    #Display correlation among models.
    #Weakly correlated algorithms are persuasive for stacking them in ensemble.
    Global.Mod.Cor<-caret::modelCor(resamples(global_train))
    
    
    #--------------------------------------------
    #---------- Create ensemble model -----------
    #--------------------------------------------
    #combine individual models into one
    set.seed(478)
    global_stack <- caretEnsemble(
      global_train, 
      metric="Accuracy",
      trControl=trainControl(method="cv",
                             number=10,
                             savePredictions= "final",
                             classProbs=TRUE))
    print(global_stack)
    
    
    #--------------------------------------------
    #Identify best threshold and get accurracy
    #--------------------------------------------
    #Identify threshold that maximizes spec=sens
    global.ens.thresh<-findThresh(global_stack$ens_model$pred)
    
    #Return accurracy
    ensemble_accurracy<-accuracyStats(global_stack$ens_model$pred,global.ens.thresh$predicted)
    
    
    #--------------------------------------------
    #-- Get variable importance of global model--
    #--------------------------------------------
    variableImportance_global<-caret::varImp(global_stack)
    
    
    #--------------------------------------------
    #-------- Make predictions for Europe--------
    #--------------------------------------------
    global_model <- terra::predict(eu_climpreds.10_selection,global_stack,type="prob", na.rm = TRUE) 
    
    
    #--------------------------------------------
    #-------------- Plot predictions-------------
    #--------------------------------------------
    #brks <- seq(0, 1, by=0.1) 
    #nb <- length(brks)-1 
    # Generate Viridis palette
    #viridis_palette <- viridis(nb)
    
    #ggplot() + 
    #geom_sf(data = world,  colour = "grey", fill = NA)+
    #geom_spatraster(data = global_model) +
    #scale_fill_gradientn(colors = viridis_palette, breaks = brks, labels = brks, na.value = NA) +
    #geom_sf(data = global.occ.sf, color = "black", fill = "red", size =1.5, shape = 21) +
    #coord_sf(xlim = c(-10, 40), ylim = c(35, 72)) + 
    #labs(fill = "Suitability")+
    #theme_bw()
    
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
                        global_ensemble_model = global_stack, #Global ensemble model 
                        global_data_df_uncor=global.data.df.uncor, #Data used to fit the global ensemble model (climate data for each presence/pseudoabsence)
                        global_presabs=global_presabs,#xy coordinates of presences and pseudoabsences used to fit the models
                        model_accuracy = ensemble_accurracy, #Accuracy of ensemble model
                        variable_importance = variableImportance_global, #Variable importance in each separate model and overall
                        occurrences=global.occ.sf, #Sf dataframe of occurrence data used to fit the models
                        model_correlation = Global.Mod.Cor #Correlation between the separate models
    )
    
    qs::qsave(globalmodels, paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey,"/Global_model_",first_two_words,"_",taxonkey,".qs"))
    
    
    #--------------------------------------------
    #--Export raster layers in folder "rasters"--
    #--------------------------------------------
    #We don't store them in .qs file as some important metadate would be stored in a temp folder, which would be removed after a while 
    biasgrid_file<- file.path("./data/projects",projectname,paste0(first_two_words,"_",taxonkey),"Rasters","Interim",paste0("Biasgrid_",first_two_words,"_",taxonkey,".tif"))
    global_model_file<- file.path("./data/projects",projectname,paste0(first_two_words,"_",taxonkey),"Rasters","Global",paste0("Global_model_",first_two_words,"_",taxonkey,".tif"))
    euclimpreds_file<- file.path("./data/projects",projectname,paste0(first_two_words,"_",taxonkey),"Rasters","Interim",paste0("EU_climpreds10_",first_two_words,"_",taxonkey,".tif"))
    
    terra::writeRaster(biasgrid_sub, filename = biasgrid_file, overwrite = TRUE)
    terra::writeRaster(global_model, filename = global_model_file, overwrite = TRUE)
    terra::writeRaster(eu_climpreds.10_selection, filename = euclimpreds_file, overwrite = TRUE)
    
    
    #--------------------------------------------
    #------------------ Clean up-----------------
    #--------------------------------------------
    print(paste("Global model has been created for", species))
    rm(list = setdiff(ls(), c("p","wwf_eco","eu_climpreds.10", "split_df", "accuracyStats", "decimalplaces", "divide10", "findThresh", "predict_large_raster", "globalclimpreds_terra","bias_grid_paths", "i", "world", "projectname", "generate_pseudoabs", "create_folder")))
    
  }
  
})

#--------------------------------------------
#---------- Clean R environment--------------
#--------------------------------------------
rm(list = ls())

