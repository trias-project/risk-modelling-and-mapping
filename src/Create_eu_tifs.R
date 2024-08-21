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
 
 # Create a raster object and set crs to a Lambert Azimuthal Equal-Area projection 
 r <- raster::rasterFromXYZ(lonlatvar, crs="+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs")
 
 #Create filepath
 filename <- sub("\\.csv$", "", basename(file))
 filepath<- paste0("./data/external/climate/EU_data/tifs/",filename,".tif")
 
 #Export raster
 raster::writeRaster(r, filepath , format = "GTiff", overwrite = TRUE)
 
 #Set progressbar
print(paste0("Exported layer ",i," of ",length(eudata),": ",filename))
 
 i<- i+1
 
 #Remove redundant objects
 rm(df, lonlatvar, r)
 
}


