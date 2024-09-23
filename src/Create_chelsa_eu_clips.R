library(rnaturalearth)
library(dplyr)
library(sf)
library(terra)

#Prepare vector of Europe for cropping
europe<-rnaturalearth::ne_countries(continent="europe", scale=10)%>%
        dplyr::filter(sovereignt!="Russia")
mapview(europe[1])
europe <- sf::st_transform(europe, st_crs = "EPSG:4326" )%>%
          sf::st_crop(xmin=-50, ymin=0, xmax=60, ymax=90)%>%
          terra::vect()

#Check the extent of the vector, and make it 10 larger on each side for nicer visualisation
terra::ext(europe)
crop_ext <- terra::ext(-31.2849014959999-10, 40.1595430910001 +10, 27.6422386740001-5, 80.7700869810001+10)

#Load CHELSA climate layers and stack them
globalclimrasters <- list.files((here::here("./data/external/climate/trias_CHELSA")),pattern='tif',full.names = T) #import CHELSA data
globalclimpreds <- terra::rast(globalclimrasters)

#Crop raster files to extent of Europe and then mask them 
ras_crop <- terra::crop(globalclimpreds, crop_ext)
ras_mask <- terra::mask(ras_crop, europe) 
terra::plot(ras_mask[[1]])

# Save each raster layer as an individual GeoTIFF
output_dir <- "./data/external/climate/chelsa_eu_clips/" 

for (i in 1:nlyr(ras_mask)) {
  # Generate a filename for each layer
  layer_name <- names(ras_mask)[i]  
  output_filepath <- paste0(output_dir, layer_name, "_eu.tif")
  
  # Export the layer as a GeoTIFF
  terra::writeRaster(ras_mask[[i]], filename = output_filepath,  overwrite = TRUE)
  
  #Follow progress
  print(paste0("Exported raster layer ",i," of ",nlyr(ras_mask),": ",layer_name,"_eu.tif"))
  
  #Remove the file
  unlink(ras_mask[[i]])
}
