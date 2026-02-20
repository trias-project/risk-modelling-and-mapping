#---------------------------------------------------------------
#---------------------------------------------------------------
# This script processes the CORINE land use dataset
# (version U2018_CLC2018_V2020_20u1). We start with the
# 100m-resolution GeoTIFF file and aggregate its 100m x 100m
# pixels into 1km x 1km cells. Each cell represents the
# percentage coverage of a specific land cover cateogry.
# These maps were created for eight land cover categories:
# Agriculture, Deciduous forest, Coniferous forest, mixed forest,
# Shrubs and herbaceous, Inland wetland, Coastal wetland, and
# artificial structures. 1km grid cells based on more than 5% NA's
# or more than 50% Sea and ocean cells were set to NA.
#---------------------------------------------------------------
#---------------------------------------------------------------
 

#--------------------------------------------
#-------------- Load packages ---------------
#--------------------------------------------
packages <- c("terra", "dplyr")

for (package in packages) {
  print(package)
  if (!package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#-------- Source helper functions -----------
#--------------------------------------------
source(here::here("src", "helper_functions.R"))


#--------------------------------------------
#------- Load original CORINE raster --------
#--------------------------------------------
corine_layers <- terra::rast(here::here("data", "external", "habitat", "u2018_clc2018_v2020_20u1_raster100m", "DATA", "U2018_CLC2018_V2020_20u1.tif" ))
terra::plot(corine_layers)
levels(corine_layers)


#--------------------------------------------
#------- Convert NODATA to NA --------
#--------------------------------------------
corine_layers[corine_layers == "NODATA"] <- NA


#--------------------------------------------
#------- Define land use categories ---------
#--------------------------------------------
#AGRICULTURE:(https://www.eea.europa.eu/publications/COR0-part2/land_coverPart2.2.pdf)
Agriculture <- c("Non-irrigated arable land", "Permanently irrigated land", "Rice fields", "Vineyards", "Fruit trees and berry plantations", "Olive groves", "Pastures", "Annual crops associated with permanent crops", "Complex cultivation patterns", "Land principally occupied by agriculture, with significant areas of natural vegetation", "Agro-forestry areas")

#CONIFEROUS
Coniferous_forest <- c("Coniferous forest")

#DECIDUOUS
Deciduous_forest <- c("Broad-leaved forest")

#MIXED FOREST
Mixed_forest <- c("Mixed forest")

#Shrub and/or herbaceous vegetation associations
Shrub_and_herbaceous <- c("Natural grasslands", "Moors and heathland", "Sclerophyllous vegetation", "Transitional woodland-shrub")

#INLAND WETLAND
Inland_wetland <- c("Inland marshes", "Peat bogs")

#COASTAL WETLANDS
Coastal_wetland <- c("Salt marshes", "Salines", "Intertidal flats")

#ARTIFICIAL SURFACES
Artificial <- c("Continuous urban fabric", "Discontinuous urban fabric", "Industrial or commercial units", "Road and rail networks and associated land", "Port areas", "Airports","Mineral extraction sites", "Dump sites", "Construction sites", "Green urban areas", "Sport and leisure facilities")

categories <- c("Agriculture", "Coniferous_forest", "Deciduous_forest", "Mixed_forest", "Shrub_and_herbaceous", "Inland_wetland", "Coastal_wetland", "Artificial")


#-------------------------------------------------
#----------- Create marine and NA mask -----------
#-------------------------------------------------
# Marine cells get a 1, others 0
ocean_mask <- terra::ifel(corine_layers == "Sea and ocean", 1, 0)

# Aggregate from 100m -> 1km
ocean_mask <- terra::aggregate(ocean_mask, fact = 10, fun = mean, na.rm = TRUE) * 100

#Convert masking cells (with more than 50% ocean pixels) to 1, rest to 0
ocean_mask <- terra::ifel(ocean_mask > 50, 1, 0)

#Define function to calculate proportion of NAs
na_fraction <- function(x) {
  sum(is.na(x)) / length(x)
}

# Aggregate from 100m -> 1km
NA_mask <- terra::aggregate(corine_layers, fact = 10, fun = na_fraction) * 100

#Convert masking cells (with more than 5% NA pixels) to 1, rest to 0
NA_mask <- terra::ifel(NA_mask > 5, 1, 0)

#-------------------------------------------------
#- Create rasters with % cover for each category -
#-------------------------------------------------
habitat_list <- list()

for (landuse_category in categories){
  
  subcategories <- get(landuse_category)
  
  # Landuse category cells get a 1, others 0
  bin_raster <- terra::ifel(corine_layers %in% subcategories, 1, 0)
  
  # Aggregate from 100m -> 1km
  bin_raster_pct <- aggregate(bin_raster, fact = 10, fun = mean, na.rm = TRUE) * 100
  
  #Mask 1km cells that consist for more than 50% of 100m marine grid cells
  bin_raster_mask <- terra::mask(bin_raster_pct, ocean_mask, maskvalues=1)
  
  #Mask 1km cells that consist for more than 5% of 100m NA grid cells
  bin_raster_mask <- terra::mask(bin_raster_mask, NA_mask, maskvalues=1)
  
  habitat_list[[landuse_category]] <- bin_raster_mask
  
  print(paste0("Habitat layer created for ", landuse_category))
}

# Combine into multilayer SpatRaster
habitat_rasters <- rast(habitat_list)
names(habitat_rasters) <- categories


#-------------------------------------------------
#------------ Export raster layers ---------------
#-------------------------------------------------
#Quality check in QGIS: everything looked OK + alignment with other layers is also ok
habitat_folder <- here::here("data","external","habitat", "new_zenodo_layers")
create_folder(habitat_folder, "new_zenodo_layers")

for (i in 1:nlyr(habitat_rasters)) {
  terra::writeRaster(habitat_rasters[[i]], filename = here::here(habitat_folder, paste0(names(habitat_rasters[[i]]), ".tif")), overwrite=TRUE)
}
