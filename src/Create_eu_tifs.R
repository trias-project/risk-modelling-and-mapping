#Load packages
library(readr)
library(dplyr)
library(raster)
library(pbapply)

#Import coordinate files and store them as vector format
lat_1km <- readr::read_csv(here::here("./data/external/climate/EU_data/latlon_data/lat_1km.csv"),  col_select=c(2:6072), col_types = "n")
lon_1km <- readr::read_csv(here::here("./data/external/climate/EU_data/latlon_data/lon_1km.csv"),  col_select=c(2:6072), col_types = "n") #Not possible to use row.names 1 due to NA fields

lat_vec <- as.vector(as.matrix(lat_1km))
lon_vec <- as.vector(as.matrix(lon_1km))
rm(lat_1km, lon_1km)

#Read in csv's and convert to raster file, export as tif
eudata <- list.files(("./data/external/climate/EU_data"),pattern='csv',full.names = T) #import CHELSA data

#initiate progress bar
i<-1

for(file in eudata){
  
 # Read file
 df<-suppressWarnings(suppressMessages(read_csv(file, col_select=c(2:6072), col_types = "n")))
 
 #Convert to vector
 df<- as.vector(as.matrix(df))
 
 # Create a data frame with lon, lat, and variable data
 lonlatvar <- data.frame(lon = lon_vec, lat = lat_vec, var1 = df)
 
 # Create a raster object and set crs to WGS84
 r <- raster::rasterFromXYZ(lonlatvar, crs="epsg:4326")
 
 #Reproject raster to epsg:3035
 r<-terra::project(rast(r), terra::crs(habitatexample))
 r<-crop(r, ext(habitatexample))
 #Create filepath
 filename <- sub("\\.csv$", "", basename(file))
 filepath<- paste0("./data/external/climate/EU_data/tifs/",filename,".tif")
 
 # Explicitly remove existing file
 if (file.exists(filepath)) {
   file.remove(filepath)
 }
 
 #Export raster
 terra::writeRaster(r, filepath , filetype = "GTiff", overwrite = TRUE)
 
 #Set progressbar
print(paste0("Exported layer ",i," of ",length(eudata),": ",filename))
 
 i<- i+1
 
 #Remove redundant objects
 rm(df, lonlatvar, r)
 
}

rm(lat_vec, lon_vec)





habitat<-list.files((here("./data/external/habitat")),pattern='tif',full.names = T)
habitat_stack<-stack(habitat[1:5])
fullstack<-stack(raster(r),habitat_stack) #combine uncorrelated climate variable selected earlier with habitat
