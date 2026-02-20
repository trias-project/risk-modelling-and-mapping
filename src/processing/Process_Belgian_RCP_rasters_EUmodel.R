#--------------------------------------------
#----------  Reason  --------
#--------------------------------------------
#In the original Belgian RCP raster layers, there is a band of pixels along the coast that all have a value of 0 instead of NA
#We will convert these pixels to NA by masking them a max temp raster where these values were converted to NA
#The latter raster was chosen as no other pixels in Belgium were 0, so we could be sure to only convert these erroneous values to NA

#--------------------------------------------
#----------  Load packages  --------
#--------------------------------------------
library(sf)
library(terra)
library(dplyr)
library(here)

#--------------------------------------------
#--Create RCP scenario files for Belgium  ---
#--------------------------------------------
be26 <- list.files(here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26"), pattern = 'tif', full.names = TRUE)
be45 <- list.files(here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45"), pattern = 'tif', full.names = TRUE)
be85 <- list.files(here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85"), pattern = 'tif', full.names = TRUE)


#----------------------------------------------------------------
#--Get raster where this these pixels are 0 and convert to NA ---
#----------------------------------------------------------------
mask_raster<-rast("C:/Users/soria_delva/Documents/GitHub/risk-modelling-and-mapping/./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26/maxtemp.tif")
mask_raster[mask_raster==0]<-NA


#--------------------------------------------
#--Create stacks country RCP rasters  ---
#--------------------------------------------
set.seed(123)

list_names <- c("be26", "be45", "be85")

# Iterate over the list names
for (list_name in list_names) {
  # get the list (be26, be45, and be85)
  current_list <- get(list_name)
  
  fullstack <- lapply(current_list, function(f) {
    r <- rast(f)
    r <- mask(r, mask_raster)
    return(r)
  })
  
  fullstack <- do.call(c, fullstack)  # Combine 
  
  
  # Assign the processed fullstack to a variable dynamically (e.g., fullstack26, fullstack45, fullstack85)
  assign(paste0("fullstack", substr(list_name, 3, 4)), fullstack)
}


#--------------------------------------------
#--------------- Export rasters -------------
#--------------------------------------------
#RCP 2.6 
for (i in 1:nlyr(fullstack26)) {
  writeRaster(fullstack26[[i]], filename = here(("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26"), paste0(names(fullstack26[[i]]), ".tif")), overwrite=TRUE)
}

#RCP 4.5
for (i in 1:nlyr(fullstack45)) {
  writeRaster(fullstack45[[i]], filename = here(("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45"), paste0(names(fullstack45[[i]]), ".tif")), overwrite=TRUE)
}

#RCP 8.5
for (i in 1:nlyr(fullstack85)) {
  writeRaster(fullstack85[[i]], filename = here(("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85"), paste0(names(fullstack85[[i]]), ".tif")), overwrite=TRUE)
}
