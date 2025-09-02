#--------------------------------------------
#-----  TO DO: define country of interest----
#--------------------------------------------
country_of_interest <- "Belgium"


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
#--------- Source helper functions ----------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#------------ Define folder paths -----------
#--------------------------------------------
habitat_folder <- here::here("data","external","habitat")
country_folder <- here::here("data","external","GIS","Country")
Europe_folder <- here::here("data","external","GIS","Europe")
chelsa_eu_folder <- here::here("data","external","climate","chelsa_eu_clips")
rcp26_country_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","country_rcps","rcp26")
rcp70_country_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","country_rcps","rcp70")
rcp85_country_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","country_rcps","rcp85")
rcp26_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","rcp26")
rcp70_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","rcp70")
rcp85_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","rcp85")
global_climate_folder <- here::here("data","external","climate","trias_CHELSA")
chelsa_mask_folder <- here::here("data","external","climate","chelsa_mask")
scaled_layers_folder <- here::here("data","external","climate","scaled_layers")
biasgrids_folder <- here::here("data","external","bias_grids","final","trias")


#-------------------------------------------------
#--------------- Create EU folders ---------------
#-------------------------------------------------
# Define the folder paths
folder_paths<-list(list("path"= habitat_folder,
                        "name"= "habitat"),
                   list("path"= country_folder,
                        "name"= "Country"),
                   list("path"= Europe_folder,
                        "name"= "Europe"),
                   list("path"= rcp26_country_eumodel_folder,
                        "name"= "rcp26"),
                   list("path"= rcp45_country_eumodel_folder,
                        "name"= "rcp45"),
                   list("path"= rcp85_country_eumodel_folder,
                        "name"= "rcp85"),
                   list("path"= chelsa_eu_folder,
                        "name"= "chelsa_eu_clips"),
                   list("path"= rcp26_country_globalmodel_folder,
                        "name"= "rcp26"),
                   list("path"= rcp70_country_globalmodel_folder,
                        "name"= "rcp70"),
                   list("path"= rcp85_country_globalmodel_folder,
                        "name"= "rcp85"),
                   list("path"= rcp26_globalmodel_folder,
                        "name"= "rcp26"),
                   list("path"= rcp70_globalmodel_folder,
                        "name"= "rcp70"),
                   list("path"= rcp85_globalmodel_folder,
                        "name"= "rcp85"),
                   list("path"= global_climate_folder,
                        "name"= "trias_CHELSA"),
                   list("path"= chelsa_mask_folder,
                        "name"= "CHELSA mask"),
                   list("path"= scaled_layers_folder,
                        "name"= "scaled_layers"),
                   list("path"= biasgrids_folder,
                        "name"= "trias")
                   )

# Check and create each folder if necessary
lapply(folder_paths, function(folder){
  create_folder(folder$path, folder$name)
})


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
  
  if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",i,"_1981-2010_V.2.1.tif"),
                  mode = "wb",
                  destfile = here::here(global_climate_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }else{
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",i,"_1981-2010_V.2.1.tif"),
                  destfile = here::here(global_climate_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }
}


#---------------------------------------------------------
#---------- Store future CHELSA layers  ---------
#---------------------------------------------------------
#Note that there are values missing in all future layers of Precipitation driest month bio14, this should not be a problem though as they fall outside of the EU!
#Store future CHELSA layers at global level
options(timeout = 600) #set time-out to 10 min 

for(i in c("1", "4", "5", "6","7", "12","13","14","15")){
  
  # Remove leading zeros during download
  i_download <- as.integer(i)
  
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
  
  if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp126/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp126_V.2.1.tif"),
                  mode = "wb",
                  destfile = here::here(rcp26_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }else{
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp126/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp126_V.2.1.tif"),
                  destfile = here::here(rcp26_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }
  
  if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp370/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp370_V.2.1.tif"),
                  mode = "wb",
                  destfile = here::here(rcp70_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }else{
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp370/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp370_V.2.1.tif"),
                  destfile = here::here(rcp70_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }
  
  if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp585/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp585_V.2.1.tif"),
                  mode = "wb",
                  destfile = here::here(rcp85_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }else{
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp585/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp585_V.2.1.tif"),
                  destfile = here::here(rcp85_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }
}


#--------------------------------------------------------------------
#----- Store CHELSA v1 layer as mask template for marine pixels  ----
#--------------------------------------------------------------------
#Download a V1 chelsa layer (check for marine pixels: none seem to be present)
if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
  download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav1/climatologies/bio/CHELSA_bio10_01.tif"),
                mode = "wb",
                destfile = here::here(chelsa_mask_folder,paste0("CHELSA_meantemp1.tif")))
}else{
  download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav1/climatologies/bio/CHELSA_bio10_",i,".tif"),
                destfile = here::here(chelsa_mask_folder,paste0("CHELSA_meantemp1.tif")))
}

chelsa_mask <- terra::rast(here::here(chelsa_mask_folder,paste0("CHELSA_meantemp1.tif")))


#-------------------------------------------------
#----- Store the European boundary shapefile  ----
#-------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.15102496", 
                       path=Europe_folder, 
                       files=list("EUROPE.shp", 
                                  "EUROPE.dbf",
                                  "EUROPE.shx",
                                  "EUROPE.prj",
                                  "EUROPE.sbn",
                                  "EUROPE.sbx",
                                  "EUROPE.shp.xml"), 
                       quiet=FALSE)

euboundary_vect <- sf::st_read(here::here(Europe_folder,"EUROPE.shp")) %>%
  terra::vect()


#-------------------------------------------------
#-------- Scale and mask CHELSA layers  ----------
#-------------------------------------------------
#List files
chelsa_current_climate <- list.files(here::here(global_climate_folder), pattern = 'tif', full.names = TRUE)

chelsa_climate_rcp26 <- list.files(here::here(rcp26_globalmodel_folder), pattern = 'tif', full.names = TRUE)

chelsa_climate_rcp70 <- list.files(here::here(rcp70_globalmodel_folder), pattern = 'tif', full.names = TRUE)
  
chelsa_climate_rcp85 <- list.files(here::here(rcp85_globalmodel_folder), pattern = 'tif', full.names = TRUE)

list_names <- c("chelsa_current_climate", "chelsa_climate_rcp26", "chelsa_climate_rcp70", "chelsa_climate_rcp85")
 
# Iterate over the list names
for (list_name in list_names) {
  # get the list 
  current_list <- get(list_name)
  
  #Define folder to store the rasters in
  rcp_folder <- switch(list_name,
                       "chelsa_current_climate" = global_climate_folder,
                       "chelsa_climate_rcp26" = rcp26_globalmodel_folder,
                       "chelsa_climate_rcp70" = rcp70_globalmodel_folder,
                       "chelsa_climate_rcp85" = rcp85_globalmodel_folder)
  
  for (file in current_list){
    
  #Open climate layer as spatRaster
  rast_file <- terra::rast(file)
  
  #Mask marine pixels
  rast_file <- terra::mask(rast_file, chelsa_mask)
  
  # #Convert units of temp seasonality layer to °C: not necessary when you scale afterwards
  # if(names(rast_file) == "CHELSA_temp_seasonality_4"){
  #   rast_file <- rast_file/100
  #   print("Converted the unit of layer bio 4 (temperature seasonality) to °C") 
  # }
  
  #Scale layer
  rast_file <- terra::scale(rast_file, center=TRUE, scale=TRUE)
  
  # Create output filename
  out_name <- here::here(rcp_folder, paste0("scaled_layer_", basename(file)))
  
  #Write raster to scaled_layers folder
  terra::writeRaster(rast_file, filename = out_name, overwrite = TRUE)
  
  #Print write statement
  print(paste0("Created rasterlayer ", basename(out_name)," in folder ", basename(rcp_folder)))
  
  #Store current EU layers
  if(list_name == "chelsa_current_climate"){
    
    #Mask current rasters with europe shape
    rast_file_eu <- terra::crop(rast_file, ext(-18.69139, 36.5828, 29.80069, 76.13302))
    rast_file_eu <- terra::mask(rast_file_eu, euboundary_vect)
    
    # Create output filename
    out_name <- here::here(chelsa_eu_folder, paste0("scaled_eu_layer_", basename(file)))
    
    #Write raster to scaled_layers folder
    terra::writeRaster(rast_file_eu, filename = out_name, overwrite = TRUE)
    print(paste0("Created rasterlayer ", basename(out_name)," in folder ", basename(chelsa_eu_folder)))
  }
  
  rm(rast_file, rast_file_eu)
}
}


#------------------------------------------------------------------------
#--Create future climate layers for a specific country (global model)  --
#------------------------------------------------------------------------
#This may take some time!
country <- rnaturalearth::ne_countries(country=country_of_interest, scale=10)[1]
country_vector <- terra::vect(country) #Convert to a SpatVector, used for masking
country_ext <- terra::ext(country_vector) 

#List files
country26 <- list.files(rcp26_globalmodel_folder , pattern = 'tif', full.names = TRUE)
country70 <- list.files(rcp70_globalmodel_folder , pattern = 'tif', full.names = TRUE)
country85 <- list.files(rcp85_globalmodel_folder , pattern = 'tif', full.names = TRUE)

# List scenarios to iterate through
set.seed(123)
list_names <- c("country26", "country70", "country85")

# Iterate over the list names
for (list_name in list_names) {
  # get the list ("country26", "country70", "country85")
  current_list <- get(list_name)
  
  #Define folder to store the rasters in
  rcp_folder <- switch(list_name,
                       "country26" = rcp26_country_globalmodel_folder,
                       "country70" = rcp70_country_globalmodel_folder,
                       "country85" = rcp85_country_globalmodel_folder)
  
  fullstack <- lapply(current_list, function(f) {
    r <- terra::rast(f)
    r <- terra::crop(r, country_ext)
    r <- terra::mask(r, country_vector)
    terra::writeRaster(r, filename = here::here(rcp_folder, basename(f)), overwrite=TRUE)

  })
}
  

#-------------------------------------------------
#----- Store the country boundary shapefile  -----
#-------------------------------------------------
sf::write_sf(country, here::here(country_folder,"country.shp"))


#-------------------------------------------------
#-- Store habitat layers for the European model --
#-------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.7841324", 
                       path=habitat_folder, 
                       quiet=FALSE)


#-------------------------------------------------
#---------------- Store biasgrids  ---------------
#-------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.7556851", 
                       path=biasgrids_folder, 
                       files=list("amphib_1deg_min5.tif",
                                  "birds_1deg_min5.tif",
                                  "mammals_1deg_min5.tif",
                                  "molluscs_1deg_min5.tif",
                                  "plants_1deg_min5.tif",
                                  "reptiles_1deg_min5.tif"), 
                       quiet=FALSE)

