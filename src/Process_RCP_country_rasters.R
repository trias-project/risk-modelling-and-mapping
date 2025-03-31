#--------------------------------------------
#----------  Load packages  --------
#--------------------------------------------
library(sf)
library(terra)
library(dplyr)
library(here)

#--------------------------------------------
#----------  To do: specify country  --------
#--------------------------------------------
#If you'd like to predict for another country, change the country name
belgium<-rnaturalearth::ne_countries(country="Belgium", scale=10)[1]
belgium_ext<-terra::ext(belgium) 
belgium_vector <- terra::vect(belgium) #Convert to a SpatVector, used for masking

#--------------------------------------------
#--Create RCP scenario files for Belgium  ---
#--------------------------------------------
be26 <- list.files(here("./data/external/climate/Global_finalRCP/rcp26"), pattern = 'tif', full.names = TRUE)
be45 <- list.files(here("./data/external/climate/Global_finalRCP/rcp45"), pattern = 'tif', full.names = TRUE)
be85 <- list.files(here("./data/external/climate/Global_finalRCP/rcp85"), pattern = 'tif', full.names = TRUE)


#--------------------------------------------
#--Create stacks country RCP rasters  ---
#--------------------------------------------
# IMPORTANT: note that due to something I don't know, every time you run this, the layers are a little different (tiny changes, but still a bit weird)!
# List scenarios to iterate through
set.seed(123)
list_names <- c("be26", "be45", "be85")

# Iterate over the list names
for (list_name in list_names) {
  # get the list (be26, be45, and be85)
  current_list <- get(list_name)
  
  fullstack <- lapply(current_list, function(f) {
    r <- rast(f)
    r[r==-32768]<-NA #Set marine pixels to NA
    r <- crop(r, belgium_ext)
    r <- mask(r, belgium_vector)
    r<- r/10 #Divide by 10 to comply with predictors in global model framework
    return(r)
  })
  
  fullstack <- do.call(c, fullstack)  # Combine 
  
  #Rename layers: they need to have the same names as the layers used to fit the global model
  layer_names<-c("CHELSA_meantemp_01",
                 "CHELSA_annPrecip_12", 
                 "CHELSA_precipWettestMon_13", 
                 "CHELSA_precipDriestMon_14",
                 "CHELSA_precipSeasonality_15",
                 "CHELSA_temp_seasonality_04",
                 "CHELSA_maxTmpWarmestMon_05",
                 "CHELSA_minTmpColdestMon_06",
                 "CHELSA_temp_annRange_07")
  names(fullstack) <- layer_names
  
  # Assign the processed fullstack to a variable dynamically (e.g., fullstack26, fullstack45, fullstack85)
  assign(paste0("fullstack", substr(list_name, 3, 4)), fullstack)
}


#--------------------------------------------
#--------------- Export rasters -------------
#--------------------------------------------
#RCP 2.6 
for (i in 1:nlyr(fullstack26)) {
  writeRaster(fullstack26[[i]], filename = here(("./data/external/climate/Global_finalRCP/belgium_rcps/rcp26"), paste0(names(fullstack26[[i]]), ".tif")), overwrite=TRUE)
}

#RCP 4.5
for (i in 1:nlyr(fullstack45)) {
  writeRaster(fullstack45[[i]], filename = here(("./data/external/climate/Global_finalRCP/belgium_rcps/rcp45"), paste0(names(fullstack45[[i]]), ".tif")), overwrite=TRUE)
}

#RCP 8.5
for (i in 1:nlyr(fullstack85)) {
  writeRaster(fullstack85[[i]], filename = here(("./data/external/climate/Global_finalRCP/belgium_rcps/rcp85"), paste0(names(fullstack85[[i]]), ".tif")), overwrite=TRUE)
}
