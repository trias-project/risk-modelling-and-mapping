#--------------------------------------------
#-----------  Load packages  ----------------
#--------------------------------------------
packages <- c("curl", "zen4R", "here", "sf", "terra")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}



#--------------------------------------------
#--------- Source helper functions ----------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#------------ Define folder paths -----------
#--------------------------------------------
habitat_folder <- file.path("./data/external/habitat")
Belgium_folder <- file.path("./data/external/GIS/Belgium")
Europe_folder <- file.path("./data/external/GIS/Europe")
Ecoregions_folder <- file.path("./data/external/GIS")
rcp26_belgium_eumodel_folder <- file.path("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26")
rcp45_belgium_eumodel_folder <- file.path("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45")
rcp85_belgium_eumodel_folder <- file.path("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85")
chelsa_eu_folder <- file.path("./data/external/climate/chelsa_eu_clips")
rcp26_belgium_globalmodel_folder <- file.path("./data/external/climate/Global_finalRCP/belgium_rcps/rcp26")
rcp45_belgium_globalmodel_folder <- file.path("./data/external/climate/Global_finalRCP/belgium_rcps/rcp45")
rcp85_belgium_globalmodel_folder <- file.path("./data/external/climate/Global_finalRCP/belgium_rcps/rcp85")
rcp26_globalmodel_folder <- file.path("./data/external/climate/Global_finalRCP/rcp26")
rcp45_globalmodel_folder <- file.path("./data/external/climate/Global_finalRCP/rcp45")
rcp85_globalmodel_folder <- file.path("./data/external/climate/Global_finalRCP/rcp85")
eu_climate_folder <- file.path("./data/external/climate/rmi_corrected")
global_climate_folder <- file.path("./data/external/climate/trias_CHELSA")
biasgrids_folder <- file.path("./data/external/bias_grids/final/trias")


#-------------------------------------------------
#--------------- Create EU folders ---------------
#-------------------------------------------------
# Define the folder paths
folder_paths<-list(list("path"= habitat_folder,
                        "name"= "habitat"),
                   list("path"= Belgium_folder,
                        "name"= "Belgium"),
                   list("path"= Europe_folder,
                        "name"= "Europe"),
                   list("path"= Ecoregions_folder,
                        "name"= "official"),
                   list("path"= rcp26_belgium_eumodel_folder,
                        "name"= "rcp26"),
                   list("path"= rcp45_belgium_eumodel_folder,
                        "name"= "rcp45"),
                   list("path"= rcp85_belgium_eumodel_folder,
                        "name"= "rcp85"),
                   list("path"= chelsa_eu_folder,
                        "name"= "chelsa_eu_clips"),
                   list("path"= rcp26_belgium_globalmodel_folder,
                        "name"= "rcp26"),
                   list("path"= rcp45_belgium_globalmodel_folder,
                        "name"= "rcp45"),
                   list("path"= rcp85_belgium_globalmodel_folder,
                        "name"= "rcp85"),
                   list("path"= rcp26_globalmodel_folder,
                        "name"= "rcp26"),
                   list("path"= rcp45_globalmodel_folder,
                        "name"= "rcp45"),
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
#--------- Download and store data ---------------
#-------------------------------------------------
#The WWF ecoregions file
curl::curl_download("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip?1349272619", destfile = paste0(Ecoregions_folder,"/official.zip"))
unzip(paste0(Ecoregions_folder,"/official.zip"), exdir = paste0(Ecoregions_folder,"/official"))
unlink(paste0(Ecoregions_folder,"/official.zip"))

#The CHELSA layers
