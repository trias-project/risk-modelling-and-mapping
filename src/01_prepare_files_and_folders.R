#--------------------------------------------
#-----------  Load packages  ----------------
#--------------------------------------------
packages <- c("curl", "zen4R", "here", "sf", "terra", "rnaturalearth", "devtools")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#---  Install and load rnaturalearthhires  --
#--------------------------------------------
if ("rnaturalearthhires" %in% rownames(installed.packages())) {
  library(rnaturalearthhires)
} else {
  devtools::install_github("ropensci/rnaturalearthhires")
  library(rnaturalearthhires)
}


#--------------------------------------------
#-Source helper functions and configurations-
#--------------------------------------------
source(here::here("src", "helper_functions.R"))
source(here::here("src", "00_configurations.R"))


#-------------------------------------------------
#--------------- Create folders ---------------
#-------------------------------------------------
# Define the folder paths
habitat_folder <- here::here("data","external","habitat")
country_folder <- here::here("data","external","GIS","Country")
chelsa_current_folder <- here::here("data","external","climate","chelsa_current")
chelsa_mask_folder <- here::here("data","external","climate","chelsa_mask")
biasgrids_folder <- here::here("data","external","bias_grids")

#Store in a vector
folders<-c(habitat_folder, country_folder, chelsa_current_folder, chelsa_mask_folder, biasgrids_folder)

# Check and create each folder if necessary
for(folder in folders){
  if(!dir.exists(folder)) dir.create(folder, recursive=TRUE)
}


#-------------------------------------------------
#--------- Store the CHELSA layers  --------------
#-------------------------------------------------
options(timeout = 600) #set time-out to 10 min 

for(i in c("1", "4", "5", "6","7", "12","13","14","15")){
  
  # Define CHELSA layer name
  layer_name <- switch(i,
                       "1" = "meantemp",
                       "4" = "temp_seasonality",
                       "5" = "maxTmpWarmestMon",
                       "6"= "minTmpColdestMon",
                       "7"="temp_annRange",
                       "12"="annPrecip",
                       "13"="precipWettestMon",
                       "14"="precipDriestMon",
                       "15"="precipSeasonality")
  
  destfile <- here::here(chelsa_current_folder,paste0("CHELSA_",layer_name,"_",i,".tif"))
  
  if( update_files_logic(dest_file = destfile,
                         update_files)){
    if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
      download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",i,"_1981-2010_V.2.1.tif"),
                    mode = "wb",
                    destfile = destfile)
    }else{
      download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",i,"_1981-2010_V.2.1.tif"),
                    destfile = destfile)
    }
  }
}


#--------------------------------------------------------------------
#----- Store CHELSA v1 layer as mask template for marine pixels  ----
#--------------------------------------------------------------------
#Download a V1 chelsa layer (check for marine pixels: none seem to be present)
destfile <- here::here(chelsa_mask_folder,paste0("CHELSA_meantemp1.tif"))
if(update_files_logic(dest_file = destfile,
                      update_files)){
  if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav1/climatologies/bio/CHELSA_bio10_01.tif"),
                  mode = "wb",
                  destfile = destfile)
  }else{
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav1/climatologies/bio/CHELSA_bio10_",i,".tif"),
                  destfile = destfile)
  }
  
  chelsa_mask <- terra::rast(destfile)
}

#-------------------------------------------------
#----- Scale and mask current CHELSA layer  ------
#-------------------------------------------------
#List files
chelsa_current <- list.files(here::here(chelsa_current_folder), pattern = "^CHELSA_.*\\.tif$", full.names = TRUE)

for (file in chelsa_current){
  
  #Create layer name 
  layer_name<-sub("\\.tif$", "", basename(file))
  print(paste0("Processing ", layer_name))
  
  #Mask raster
  masked_r<- terra::mask(terra::rast(file), chelsa_mask)
  
  #Obtain mean and sd of raster and use for scaling
  m<-global(masked_r,"mean",na.rm=TRUE)$mean
  s<-global(masked_r,"sd",na.rm=TRUE)$sd
  scaled_r<-(masked_r - m) / s
  
  #Round the scaled layer
  scaled_r<-terra::round(scaled_r, 2)
  
  #Assign name to raster
  names(scaled_r)<-layer_name
  
  # #Convert units of temp seasonality layer to °C: not necessary when you scale afterwards
  # if(names(rast_file) == "CHELSA_temp_seasonality_4"){
  #   rast_file <- rast_file/100
  #   print("Converted the unit of layer bio 4 (temperature seasonality) to °C") 
  # }
  
  # #Write raster to disk
  out_name <- here::here(chelsa_current_folder, paste0("scaled_layer_", layer_name,".tif"))
  terra::writeRaster(scaled_r, filename = out_name, overwrite = TRUE)
  
  #Print write statement
  print(paste0("Created rasterlayer ", basename(out_name)," in folder ", basename(chelsa_current_folder)))
  
  #Clean up
  rm(masked_r, m, s, scaled_r, layer_name, out_name)
  gc()
}


#---------------------------------------------------------
#---------- Store future CHELSA layers  ---------
#---------------------------------------------------------
#Note that there are values missing in all future layers of Precipitation driest month (bio14), this is not a problem as they fall outside of the EU!
for (period in c("2041-2070","2071-2100")) {
  for (scenario in c("ssp126","ssp370","ssp585")) {
    
    future_folder <- here::here("data","external", "climate", "chelsa_future", period, scenario)
    if(!dir.exists(future_folder)) dir.create(future_folder, recursive=TRUE)
    
    zen4R::download_zenodo(
      doi = "10.5281/zenodo.17724735",
      path = future_folder,
      files = c(paste0("scaled_layer_CHELSA_meantemp_1_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_temp_seasonality_4_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_maxTmpWarmestMon_5_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_minTmpColdestMon_6_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_temp_annRange_7_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_annPrecip_12_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_precipWettestMon_13_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_precipDriestMon_14_",period,"_",scenario,".tif"), 
                paste0("scaled_layer_CHELSA_precipSeasonality_15_",period,"_",scenario,".tif")
      ),
      timeout=600,
      quiet = FALSE
    )
  }
}


# #-------------------------------------------------
# #-- Store habitat layers for the European model --
# #-------------------------------------------------

zen4R::download_zenodo(doi = "10.5281/zenodo.17724735", 
                       path = habitat_folder, 
                       files = list("Agriculture.tif",
                                    "Artificial.tif",
                                    "Coastal_wetland.tif",
                                    "Coniferous_forest.tif",
                                    "Deciduous_forest.tif",
                                    "Inland_wetland.tif",
                                    "Mixed_forest.tif",
                                    "Shrub_and_herbaceous.tif",
                                    "log_distance_to_water.tif",
                                    "log_total_water_length.tif",
                                    "proportion_total_water_polygon_cover.tif"),
                       quiet=FALSE)


# #-------------------------------------------------
# #---------------- Store biasgrids  ---------------
# #-------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.17724735", 
                       path=biasgrids_folder,
                       files=list("log_amphibians_1degree_layer.tif",
                                  "log_birds_1degree_layer.tif",
                                  "log_fish_1degree_layer.tif",
                                  "log_hydrozoa_1degree_layer.tif",
                                  "log_insects_1degree_layer.tif",
                                  "log_malacostraca_1degree_layer.tif",
                                  "log_mammals_1degree_layer.tif",
                                  "log_mollusca_1degree_layer.tif",
                                  "log_plants_1degree_layer.tif",
                                  "log_reptiles_1degree_layer.tif"),
                       quiet=FALSE)


#-------------------------------------------------
#----- Store the country boundary shapefile  -----
#-------------------------------------------------
#This may take some time!
country <- rnaturalearth::ne_countries(country=country_of_interest, scale=10)[1]
country_vector <- terra::vect(country) #Convert to a SpatVector, used for masking
country_ext <- terra::ext(country_vector) 
sf::write_sf(country, here::here(country_folder,"country.shp"))


#--------------------------------------------
#---------- Clean R environment--------------
#--------------------------------------------
rm(list = ls())
