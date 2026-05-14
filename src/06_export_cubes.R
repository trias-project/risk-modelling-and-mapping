#--------------------------------------------
#--------------- Load packages --------------
#--------------------------------------------
packages <- c( "dplyr", "sf", "terra", "arrow", "ggplot2", "plyr")
for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#- Load helper functions and configurations -
#--------------------------------------------
source(file.path("src", "helper_functions.R"))
source(file.path("src", "00_configurations.R"))


#--------------------------------------------
#-Only run code if cube needs to be exported-
#--------------------------------------------
if (export_cubes) {
  
  #--------------------------------------------
  #------------- Load species data ------------
  #--------------------------------------------
  #Get taxa info
  taxa_info <- read.csv2(file.path("data", "projects", project, paste0(project, "_taxa_info.csv")))%>%
    dplyr::distinct(acceptedTaxonKey, .keep_all = T)
  
  
  #--------------------------------------------
  #-- EEA 1km template raster and link file---
  #--------------------------------------------
  eea_template<-terra::rast("data/external/GIS/EEA-1km-GRID/EEA-1km_template.tif")
  eea_link <- as.data.frame(arrow::read_parquet("data/external/GIS/EEA-1km-GRID/EEA_1km_link.parquet"))
  cellcode_link <- eea_link$CELLCODE
  names(cellcode_link) <- eea_link$CELLID
  
  
  #--------------------------------------------
  #--Provide scenario names---
  #--------------------------------------------
  scenario_map <- c(ssp126 = "SSP1-2.6",
                    ssp370 = "SSP3-7.0",
                    ssp585 = "SSP5-8.5")
  
  
  #--------------------------------------------
  #----- Create empty list for final cube -----
  #--------------------------------------------
  #Create folders for each combination
  EEA_climate_cube_all <- list()
  EEA_habitat_cube_all <- list()
  EEA_combined_cube_all <- list()
  
  
  for (i in seq_along(taxa_info$acceptedTaxonKey)) {
    
    #--------------------------------------------
    #----------- Load species details -----------
    #--------------------------------------------
    species <- taxa_info$acceptedScientificName[i]
    taxonkey<- taxa_info$acceptedTaxonKey[i]
    speciesName <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)  # Extract first two words of species name
    
    message(
      "\n", strrep("=", 72),
      "\nSPECIES: ", speciesName,
      "  [taxonkey: ", taxonkey, "]",
      "\n", strrep("=", 72)
    )
    
    #--------------------------------------------
    #-----------  Specify base_dir  -------------
    #--------------------------------------------
    base_dir <- file.path("data", "projects", project, paste0(speciesName, "_", taxonkey))
    
    
    #----------------------------------
    #-Only run when climate model ran -
    #----------------------------------
    #Check if habitat model was created, otherwise skip to next species
    climate_qs_file <- file.path(base_dir, "Climate",
                                 paste0("Climate_model_", speciesName, "_", taxonkey, ".qs"))
    
    if (!file.exists(climate_qs_file)) {
      warning("No climate predictions available for ",species,
              ".\nSkipping species.")
      next
    }
    
    
    #--------------------------------------------
    #---------- Define raster paths ------------
    #--------------------------------------------
    # Define outputs, periods, and scenarios
    periods   <- c("Current","2041-2070", "2071-2100")
    scenarios <- c("ssp126", "ssp370", "ssp585")
    
    #Create folders for each combination
    EEA_climate_cubes <- list()
    EEA_habitat_cubes <- list()
    EEA_combined_cubes <- list()
    
    
    #--------------------------------------------
    #-----PART 1: Create climate cubes ----------
    #--------------------------------------------
    for(period in periods){
      if(period=="Current"){
        
        
        #-------------------------------------------------------
        #--- Load climate predictions for region of interest ---
        #-------------------------------------------------------
        climatepath <- file.path(base_dir, "Climate", period, "Predictions", "Rasters", 
                                 paste0(speciesName, "_Climate_current_ensemble.tif"))
        
        if(!file.exists(climatepath)){
          message("No climate predictions available for period: ", period ,". Skipping this period.")
          next
        }else{
          message("Processing climate predictions for period: ", period)
        }
        
        #---------------------------
        #- project to eea template -
        #---------------------------
        climateraster<-terra::rast(climatepath)%>%
          terra::project("EPSG:3035",
                         method = "bilinear")
        
        climateraster<- terra::resample(climateraster,
                                        eea_template,
                                        method = "bilinear")
        
        #Get values of climate raster and CELLIDs of cells that are not NA
        vals <- values(climateraster)
        idx <- which(is.finite(vals))
      
        
        #Generate climate cube
        climate_cube <- data.frame(TaxonKey = taxonkey,
                                   Species = species,
                                   CELLID = idx,
                                   Climate_favourability = vals[idx],
                                   Time = period,
                                   Scenario = "Baseline")%>%
          dplyr::left_join(eea_link,
                           by = "CELLID" )%>%
          dplyr::select(Species, TaxonKey, CELLCODE,Time,Scenario,Climate_favourability)
        
        
        #Store in list
        EEA_climate_cubes[[period]]<-climate_cube
        
        #Clean up
        rm(climate_cube,climateraster, climatepath) 
        
      }else{
        for(scenario in scenarios){
          
          
          #---------------------------
          #--- Define scenario name ---
          #---------------------------
          scenarioTitle<- switch(scenario,
                                 "ssp126" = "SSP1-2.6",
                                 "ssp370" = "SSP3-7.0",
                                 "ssp585" = "SSP5-8.5")
          
          
          #---------------------------
          #--- Load climate raster ---
          #---------------------------
          climatepath <- file.path(base_dir, "Climate", period, scenario,"Predictions", "Rasters", 
                                   paste0(speciesName, "_Climate_",period,"_",scenario,"_ensemble.tif") )
          
          if(!file.exists(climatepath)){
            message("No climate predictions available for period: ", period ," scenario ",scenario,". Skipping this.")
            next
          }else{
            message("Processing climate predictions for period: ", period," and scenario: ", scenario)
          }
          
          
          #---------------------------
          #- project to eea template -
          #---------------------------
          climateraster<-terra::rast(climatepath)%>%
            terra::project("EPSG:3035",
                           method = "bilinear")
          
          climateraster<- terra::resample(climateraster,
                                          eea_template,
                                          method = "bilinear")
          
          #Get values of climate raster and CELLIDs of cells that are not NA
          vals <- values(climateraster)
          idx <- which(is.finite(vals))
          
          #Generate climate cube
          climate_cube <- data.frame(TaxonKey = taxonkey,
                                     Species = species,
                                     CELLID = idx,
                                     Climate_favourability = vals[idx],
                                     Time = period,
                                     Scenario = scenarioTitle)%>%
            dplyr::left_join(eea_link,
                             by = "CELLID" )%>%
            dplyr::select(Species, TaxonKey, CELLCODE,Time,Scenario,Climate_favourability)
          
          
          #Store in list
          EEA_climate_cubes[[paste0(period,"_",scenario)]]<-climate_cube
          
          #Clean up
          rm(climate_cube,climateraster, climatepath) 
        }
      }
    }
    
    #Combine into one cube
    climate_cube<-dplyr::bind_rows(EEA_climate_cubes)
    
    #Plot cube 
    plot_cube(cube= climate_cube,
              type= "Climate",
              path = file.path(base_dir, "Climate"))
    
    #Store and export cube
    EEA_climate_cube_all[[speciesName]]<-climate_cube
    write.csv(climate_cube, file.path(base_dir, "Climate", "Climate_cube.csv"))
    
    
