#--------------------------------------------
#--------------- Load packages --------------
#--------------------------------------------
packages <- c( "dplyr", "sf", "terra", "arrow", "ggplot2", "purrr")
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
                    ssp585 = "SSP5-8.5",
                    Baseline = "Baseline")
  
  
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
    base_dir <- file.path("data", "projects", project, paste0(speciesName, "_", taxonkey))
    
    message(
      "\n", strrep("=", 72),
      "\nSPECIES: ", speciesName,
      "  [taxonkey: ", taxonkey, "]",
      "\n", strrep("=", 72)
    )
    

    #--------------------------------------------
    #-----PART 1: Create climate cubes ----------
    #--------------------------------------------
    
    #Check if climate model was created, otherwise skip to next species
    climate_qs_file <- file.path(base_dir, "Climate",
                                 paste0("Climate_model_", speciesName, "_", taxonkey, ".qs"))
    
    if (!file.exists(climate_qs_file)) {
      warning("No climate predictions available for ",species,
              ".\nSkipping species.")
      next
    }
    
    #Define periods and scenarios to be included
    runs <- data.frame(Period = "Current", 
                       Scenario = "Baseline")%>%
      rbind(expand.grid(Period = c("2041-2070", "2071-2100"),
                        Scenario = c("ssp126", "ssp370", "ssp585"),
                        stringsAsFactors = FALSE))
    
    #Create climate cube
    climate_cube <- purrr::pmap_dfr(runs, function(Period, Scenario) {
      build_cube(base_dir, 
                 speciesName,
                 taxonkey,
                 Period,
                 Scenario,
                 eea_template,
                 eea_link,
                 type = "Climate",
                 scenario_map )})
    
    #Plot cube 
    plot_cube(cube= climate_cube,
              type= "Climate",
              path = file.path(base_dir, "Climate"))
    
    #Store and export cube
    EEA_climate_cube_all[[speciesName]]<-climate_cube
    write.csv(climate_cube, file.path(base_dir, "Climate", "Climate_cube.csv"))
    
    
    #--------------------------------------------
    #-----PART 2: Create habitat cubes ----------
    #--------------------------------------------
    
    #Check if habitat model was created, otherwise skip to next species
    habitat_qs_file <- file.path(base_dir, "Habitat",
                                 paste0("Habitat_model_", speciesName, "_", taxonkey, ".qs"))
    
    if (!file.exists(habitat_qs_file)) {
      warning("No habitat predictions available for ",species,
              ".\nSkipping species.")
      next
    }
    
    #Create habitat cube
    habitat_cube <- purrr::pmap_dfr(runs, function(Period, Scenario) {
      build_cube(base_dir, 
                 speciesName,
                 taxonkey,
                 Period,
                 Scenario,
                 eea_template,
                 eea_link,
                 type = "Habitat",
                 scenario_map )})
    
    #Plot cube 
    plot_cube(cube= habitat_cube,
              type= "Habitat",
              path = file.path(base_dir, "Habitat"))
    
    
    #Store and export cube
    EEA_habitat_cube_all[[speciesName]]<-habitat_cube
    write.csv(habitat_cube, file.path(base_dir, "Habitat", "Habitat_cube.csv"))
    
    
    #--------------------------------------------
    #-----PART 3: Create ensemble cubes ---------
    #--------------------------------------------
    
    #Create combined cube
    combined_cube <- purrr::pmap_dfr(runs, function(Period, Scenario) {
      build_cube(base_dir, 
                 speciesName,
                 taxonkey,
                 Period,
                 Scenario,
                 eea_template,
                 eea_link,
                 type = "Combined",
                 scenario_map )})
    
    #Plot cube 
    plot_cube(cube= combined_cube,
              type= "Combined",
              path = file.path(base_dir, "Combined"))
    
    #Store and export cube
    EEA_combined_cube_all[[speciesName]]<-combined_cube
    write.csv(combined_cube, file.path(base_dir, "Combined", "Combined_cube.csv"))
    
    gc()
  }
  
  #--------------------------------------------
  #----------- Export final cubes--------------
  #--------------------------------------------
  #Combine into one cube
  climate_cube_all<-dplyr::bind_rows(EEA_climate_cube_all)
  habitat_cube_all<-dplyr::bind_rows(EEA_habitat_cube_all)
  combined_cube_all<-dplyr::bind_rows(EEA_combined_cube_all)
  
  #Store and export cube
  arrow::write_parquet(climate_cube_all, file.path("data", "projects", project, "Climate_cube.parquet"))
  arrow::write_parquet(habitat_cube_all, file.path("data", "projects", project, "Habitat_cube.parquet"))
  arrow::write_parquet(combined_cube_all, file.path("data", "projects", project, "Combined_cube.parquet"))
  
}

