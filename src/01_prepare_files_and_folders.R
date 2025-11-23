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
Europe_folder <- here::here("data","external","GIS","Europe")
chelsa_current_folder <- here::here("data","external","climate","chelsa_current")
chelsa_mask_folder <- here::here("data","external","climate","chelsa_mask")
biasgrids_folder <- here::here("data","external","bias_grids")

#Store in a vector
folders<-c(habitat_folder, country_folder, Europe_folder, chelsa_current_folder, chelsa_mask_folder, biasgrids_folder)

# Check and create each folder if necessary
for(folder in folders){
  if(!dir.exists(folder)) dir.create(folder, recursive=TRUE)
}


#-------------------------------------------------
#----------- Create future folders ---------------
#-------------------------------------------------
for (period in c("2041-2070","2071-2100")){
  for(scenario in c("ssp126", "ssp370", "ssp585")){
    
    folder <- here::here("data","external", "climate", "chelsa_future", period, scenario)
    if(!dir.exists(folder)) dir.create(folder, recursive=TRUE)
  }
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
  
  if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",i,"_1981-2010_V.2.1.tif"),
                  mode = "wb",
                  destfile = here::here(chelsa_current_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }else{
    download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/1981-2010/bio/CHELSA_bio",i,"_1981-2010_V.2.1.tif"),
                  destfile = here::here(chelsa_current_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
  }
}


#---------------------------------------------------------
#---------- Store future CHELSA layers  ---------
#---------------------------------------------------------
#Note that there are values missing in all future layers of Precipitation driest month bio14, this should not be a problem though as they fall outside of the EU!
#Store future CHELSA layers at global level
# options(timeout = 600) #set time-out to 10 min 
# 
# for(i in c("1", "4", "5", "6","7", "12","13","14","15")){
#   
#   # Remove leading zeros during download
#   i_download <- as.integer(i)
#   
#   # Define CHELSA layer name
#   layer_name <- switch(i,
#                        "1" = "meantemp",
#                        "4" = "temp_seasonality",
#                        "5" = "maxTmpWarmestMon",
#                        "6"= "minTmpColdestMon",
#                        "7"="temp_annRange",
#                        "12"="annPrecip",
#                        "13"="precipWettestMon",
#                        "14"="precipDriestMon",
#                        "15"="precipSeasonality")
#   
#   if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
#     download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp126/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp126_V.2.1.tif"),
#                   mode = "wb",
#                   destfile = here::here(rcp26_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
#   }else{
#     download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp126/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp126_V.2.1.tif"),
#                   destfile = here::here(rcp26_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
#   }
#   
#   if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
#     download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp370/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp370_V.2.1.tif"),
#                   mode = "wb",
#                   destfile = here::here(rcp70_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
#   }else{
#     download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp370/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp370_V.2.1.tif"),
#                   destfile = here::here(rcp70_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
#   }
#   
#   if(grepl("windows", Sys.getenv("OS"), ignore.case = TRUE)) {
#     download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp585/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp585_V.2.1.tif"),
#                   mode = "wb",
#                   destfile = here::here(rcp85_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
#   }else{
#     download.file(url = paste0("https://os.zhdk.cloud.switch.ch/chelsav2/GLOBAL/climatologies/2041-2070/GFDL-ESM4/ssp585/bio/CHELSA_bio",i_download,"_2041-2070_gfdl-esm4_ssp585_V.2.1.tif"),
#                   destfile = here::here(rcp85_globalmodel_folder,paste0("CHELSA_",layer_name,"_",i,".tif")))
#   }
# }


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


#-------------------------------------------------
#----- Store the country boundary shapefile  -----
#-------------------------------------------------
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


# #-------------------------------------------------
# #-- Store habitat layers for the European model --
# #-------------------------------------------------
# zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.7841324", 
#                        path=habitat_folder, 
#                        quiet=FALSE)


# #-------------------------------------------------
# #---------------- Store biasgrids  ---------------
# #-------------------------------------------------
# zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.7556851", 
#                        path=biasgrids_folder, 
#                        files=list("amphib_1deg_min5.tif",
#                                   "birds_1deg_min5.tif",
#                                   "mammals_1deg_min5.tif",
#                                   "molluscs_1deg_min5.tif",
#                                   "plants_1deg_min5.tif",
#                                   "reptiles_1deg_min5.tif"), 
#                        quiet=FALSE)
# 
