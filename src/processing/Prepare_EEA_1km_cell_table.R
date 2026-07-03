#--------------------------------------------
#---------------- Load packages -------------
#--------------------------------------------
library(sf)
library(terra)
library(arrow)


#--------------------------------------------
#------Download EEA grid geopackage ---------
#--------------------------------------------
fistools::download_gdrive_if_missing(gfileID = "1XX9eYv9N_aBCquUm7SkN2CvVizgVvm29",
                                     destfile = "data/external/GIS/Europe/EEA-1km-GRID-COUNTRIES-2013.gpkg",
                                     update_always = TRUE,
                                     email = Sys.getenv("email"))


#--------------------------------------------
#--------- List  layers for processing ------
#--------------------------------------------
EEA_grid_Europe <- file.path("data","external", "GIS", "Europe", "EEA-1km-GRID-COUNTRIES-2013.gpkg")
layers_info <- sf::st_layers(EEA_grid_Europe)%>%
  filter(name!=rs_1km_polygon) #This throws an error, add manually

layers_info<-layers_info[38:42,]


#--------------------------------------------
#---------Create template raster  ------
#--------------------------------------------
#Create template raster with extent and CRS of EEA shapefile
eea_template <- rast(xmin = -200000,
                     xmax = 8400000,
                     ymin = 700000,
                     ymax = 7500000,
                     resolution = 1000,
                     crs = "EPSG:3035")


#--------------------------------------------
#---------Process each country ------
#--------------------------------------------
#Create df 
eea_cells<-data.frame()

#Store EEA cells of every country
for(name in layers_info$name){
  #Load country shapefile
  country <- sf::st_read(EEA_grid_Europe, layer = name)
  
  #Force geometry to multipolygon
  sf::st_geometry(country) <- sf::st_cast(
    sf::st_cast(sf::st_geometry(country), "GEOMETRY"),
    "MULTIPOLYGON")
  
  #Convert to spatvector
  country <- terra::vect(country)
  
  #Extract centroids of each country EEA polygon
  pts <- terra::centroids(country)
  
  #Link centroid of country EEA cell to a cell in the EEA template raster
  cells <- terra::cellFromXY(eea_template, terra::crds(pts))
  
  #Link EEA template raster cell to real EEA cellcode
  lookup <- data.frame(
    CELLID = cells,
    CELLCODE = country$CELLCODE)
  
  #Store
  eea_cells<-bind_rows(eea_cells,lookup)
  
}


#--------------------------------------------
#---Process Portugal and Serbia separately----
#--------------------------------------------
#Serbia throws an error in the .gpk and Pt is not included
for(name in c("rs", "pt")){
  
  #Load country shapefile
  country <- sf::st_read(file.path("data", "external", "GIS", "Europe",
                                   paste0(name,"_1km.shp")))
  
  sf::st_geometry(country) <- sf::st_cast(
    sf::st_cast(sf::st_geometry(country), "GEOMETRY"),
    "MULTIPOLYGON")
  
  #Convert to spatvector
  country <- terra::vect(country)
  
  #Extract centroids of each country EEA polygon
  pts <- terra::centroids(country)
  
  #Link centroid of country EEA cell to a cell in the EEA template raster
  cells <- terra::cellFromXY(eea_template, terra::crds(pts))
  
  #Link EEA template raster cell to real EEA cellcode
  lookup <- data.frame(
    CELLID = cells,
    CELLCODE = country$CELLCODE)
  
  #Store
  eea_cells<-bind_rows(eea_cells,lookup)
  
}


#--------------------------------------------------------
#--Create and export final list and template raster------
#--------------------------------------------------------
#Keep only unique cellcodes
eea_cells<-eea_cells%>%
  dplyr::distinct(CELLCODE, .keep_all = TRUE)

#Export lookup table
arrow::write_parquet(eea_cells,"data/external/GIS/EEA-1km-GRID/EEA_1km_link.parquet")

#Export template raster
values(eea_template) <- NA
writeRaster(eea_template, "data/external/GIS/EEA-1km-GRID/EEA-1km_template.tif", overwrite=T)
