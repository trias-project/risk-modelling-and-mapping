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
Belgium_folder <- here::here("data","external","GIS","Belgium")
Europe_folder <- here::here("data","external","GIS","Europe")
Ecoregions_folder <- here::here("data","external","GIS")
rcp26_belgium_eumodel_folder <- here::here("data","external","climate","byEEA_finalRCP","belgium_rcps","rcp26")
rcp45_belgium_eumodel_folder <- here::here("data","external","climate","byEEA_finalRCP","belgium_rcps","rcp45")
rcp85_belgium_eumodel_folder <- here::here("data","external","climate","byEEA_finalRCP","belgium_rcps","rcp85")
chelsa_eu_folder <- here::here("data","external","climate","chelsa_eu_clips")
rcp26_country_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","country_rcps","rcp26")
rcp70_country_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","country_rcps","rcp70")
rcp85_country_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","country_rcps","rcp85")
rcp26_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","rcp26")
rcp70_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","rcp70")
rcp85_globalmodel_folder <- here::here("data","external","climate","Global_finalRCP","rcp85")
eu_climate_folder <- here::here("data","external","climate","rmi_corrected")
global_climate_folder <- here::here("data","external","climate","trias_CHELSA")
biasgrids_folder <- here::here("data","external","bias_grids","final","trias")


#-------------------------------------------------
#--------------- Create EU folders ---------------
#-------------------------------------------------
# Define the folder paths
folder_paths<-list(list("path"= habitat_folder,
                        "name"= "habitat"),
                   list("path"= country_folder,
                        "name"= "Country"),
                   list("path"= Belgium_folder,
                        "name"= "Belgium"),
                   list("path"= Europe_folder,
                        "name"= "Europe"),
                   list("path"= Ecoregions_folder,
                        "name"= "GIS"),
                   list("path"= rcp26_belgium_eumodel_folder,
                        "name"= "rcp26"),
                   list("path"= rcp45_belgium_eumodel_folder,
                        "name"= "rcp45"),
                   list("path"= rcp85_belgium_eumodel_folder,
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
                   list("path"=eu_climate_folder,
                        "name"= "rmi_corrected"),
                   list("path"= global_climate_folder,
                        "name"= "trias_CHELSA"),
                   list("path"= biasgrids_folder,
                        "name"= "trias"))

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


#-------------------------------------------------
#-- Store habitat layers for the European model --
#-------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.7841324", 
                       path=habitat_folder, 
                       quiet=FALSE)


#---------------------------------------------------------
#---------- Store future climate layers (CHELSA) ---------
#---------------------------------------------------------
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


#------------------------------------------------------------------------
#------Create future climate layers for Belgium (global model)  ---------
#------------------------------------------------------------------------
#Crop and mask layers to Belgium, this may take some time!
#If you'd like to predict for another country, change the country name
belgium<-rnaturalearth::ne_countries(country="Belgium", scale=10)[1]
belgium_ext<-terra::ext(belgium) 
belgium_vector <- terra::vect(belgium) #Convert to a SpatVector, used for masking

#List files
be26 <- list.files(here("data/external/climate/Global_finalRCP/rcp26"), pattern = 'tif', full.names = TRUE)
be45 <- list.files(here("data/external/climate/Global_finalRCP/rcp45"), pattern = 'tif', full.names = TRUE)
be85 <- list.files(here("data/external/climate/Global_finalRCP/rcp85"), pattern = 'tif', full.names = TRUE)

# List scenarios to iterate through
set.seed(123)
list_names <- c("be26", "be45", "be85")

# Iterate over the list names
for (list_name in list_names) {
  # get the list (be26, be45, and be85)
  current_list <- get(list_name)
  
  fullstack <- lapply(current_list, function(f) {
    r <- terra::rast(f)
    r[r==-32768]<-NA #Set marine pixels to NA
    r <- terra::crop(r, belgium_ext)
    r <- terra::mask(r, belgium_vector)
    r<- r/10 #Divide by 10 to comply with predictors in global model framework
    return(r)
  })
  
  # Combine 
  fullstack <- do.call(c, fullstack)  
  
  # Assign the processed fullstack to a variable dynamically (e.g., fullstack26, fullstack45, fullstack85)
  assign(paste0("fullstack", substr(list_name, 3, 4)), fullstack)
}

#--------------- Export rasters -------------

#RCP 2.6 
for (i in 1:nlyr(fullstack26)) {
  terra::writeRaster(fullstack26[[i]], filename = here::here(rcp26_belgium_globalmodel_folder, paste0(names(fullstack26[[i]]), ".tif")), overwrite=TRUE)
}

#RCP 4.5
for (i in 1:nlyr(fullstack45)) {
  terra::writeRaster(fullstack45[[i]], filename = here::here(rcp45_belgium_globalmodel_folder, paste0(names(fullstack45[[i]]), ".tif")), overwrite=TRUE)
}

#RCP 8.5
for (i in 1:nlyr(fullstack85)) {
  terra::writeRaster(fullstack85[[i]], filename = here::here(rcp85_belgium_globalmodel_folder, paste0(names(fullstack85[[i]]), ".tif")), overwrite=TRUE)
}


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

#-------------------------------------------------
#----- Store the Belgium boundary shapefile  -----
#-------------------------------------------------
zen4R::download_zenodo(doi="https://doi.org/10.5281/zenodo.15102496", 
                       path=Belgium_folder, 
                       files=list("belgium_boundary.shp", 
                                  "belgium_boundary.dbf",
                                  "belgium_boundary.shx",
                                  "belgium_boundary.prj",
                                  "belgium_boundary.sbn",
                                  "belgium_boundary.sbx",
                                  "belgium_boundary.shp.xml"), 
                       quiet=FALSE)


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


