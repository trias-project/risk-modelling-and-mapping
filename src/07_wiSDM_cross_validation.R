#-------------------------------------------------------------------------------
#--------------------------    Load packages      ------------------------------
#-------------------------------------------------------------------------------
packages <- c( "dplyr", "qs", "terra", "tidyterra", "sf", "here", "matrixStats",
               "ggplot2", "dismo",  "sdm", "purrr", "ecospat", "blockCV")

installed <- rownames(installed.packages())
for (package in packages) {
  print(package)
  if (!package %in% installed) install.packages(package)
  library(package, character.only = TRUE)
}

suppressWarnings(try(sdm::installAll(), silent = TRUE))


#-------------------------------------------------------------------------------
#------------------------- Set Terra options -----------------------------------
#-------------------------------------------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")#Prevent terra from writing aux.xml files
terraOptions(memfrac = 0.4,
             tempdir = file.path(tempdir()),
             todisk = TRUE)


#-------------------------------------------------------------------------------
#---------------- Load helper functions and configurations ---------------------
#-------------------------------------------------------------------------------
source(here::here("src", "helper_functions.R"))
source(here::here("src", "00_configurations.R"))


#-------------------------------------------------------------------------------
#-------------------------- Define file paths ----------------------------------
#-------------------------------------------------------------------------------
#climate stack
climate_path <- file.path("data", "external", "climate", "chelsa_current","processed", "globalclimpreds.tif")

#EU climate stack
eu_climpreds_path <- file.path("data", "external", "climate", "chelsa_current","processed","euclimpreds.tif")

#habitat stack
habitat_path <- file.path("data", "external", "habitat", "processed", "habitat_stack.tif")

#Biome file path
biome_path<-file.path("data", "external", "GIS", "official", "newRealms.shp")


#--------------------------------------------
#------------- Load species data ------------
#--------------------------------------------
#Get taxa info
taxa_info <- read.csv2(file.path("data", "projects", project, paste0(project, "_taxa_info.csv")))

#Select unique taxonkeys
accepted_taxonkeys <- unique(taxa_info$acceptedTaxonKey)


#--------------------------------------------
#------------ Load euboundary  --------------
#--------------------------------------------
euboundary_path<-file.path("data", "external", "GIS", "Europe", "EUboundary.shp")
euboundary<-sf::st_read(euboundary_path)
euboundary_wgs84<-euboundary%>%
  sf::st_transform(crs = 4326)


#--------------------------------------------
#-----------------Load rasters---------------
#--------------------------------------------
#Load rasters 
climate_stack <- terra::rast(climate_path)
habitat_stack <- terra::rast(habitat_path)
eu_climate_stack<- terra::rast(eu_climpreds_path)%>%
  terra::project(habitat_stack[[1]])


#-----------------------------------------------------------
#- Sample background data once
#-----------------------------------------------------------
#Extract a subsample of European pixels for Boyce calculation 
set.seed(728)
eu_subsample <- terra::spatSample(
  habitat_stack[[1]],
  size = boyce_background_size, 
  method = "random", 
  na.rm = TRUE, #Ignore NA pixels
  as.points = TRUE)

# Extract climate data at eu subsample points
eu_climate_sub <- terra::extract(eu_climate_stack, eu_subsample, ID = FALSE, xy = FALSE)

# Extract habitat data at eu subsample points
eu_habitat_sub <- terra::extract(habitat_stack, eu_subsample, ID = FALSE, xy = FALSE)


#----------------------------------------------
#----------- Define results directory ---------
#----------------------------------------------
validation_dir <- file.path("data", "projects", project, "Model_validation")
dir.create(validation_dir, recursive = TRUE, showWarnings = FALSE)
validation_summary<-list()


#--------------------------------------------
#----------- Start species loop -------------
#--------------------------------------------
for (i in seq_along(accepted_taxonkeys)) {
  
  species_validation_summary<-data.frame()
  
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
  
  
  #==============================================
  #=                                            =
  #=           PART 1: Climate model            =
  #=                                            =
  #==============================================
  
  message("\n--- Part 1: Global climate model ---")
  
  
  #--------------------------------------------
  #------- Define qs file paths --------
  #--------------------------------------------
  climate_qs_file <- file.path(base_dir, "Climate",
                               paste0("Climate_model_", speciesName, "_", taxonkey, ".qs"))
  
  habitat_qs_file <- file.path(base_dir, "Habitat",
                               paste0("Habitat_model_", speciesName, "_", taxonkey, ".qs"))
  
  
  #----------------------------------------------
  #- Only do validation if climate model exists -
  #----------------------------------------------
  if (!file.exists(climate_qs_file)) {
    warning("No climate model was found for ",species,
            "\nRun 03_fit_climate_model.R first.\n Skipping species.")
    next
  }
  
  
  #--------------------------------------------
  #--- Load  data stored in climate model qs --
  #--------------------------------------------
  climatemodel   <- qs::qread(climate_qs_file)
  top5_methods  <- climatemodel$top5_models
  global_presabs <- climatemodel$global_presabs
  climate_predictors <- climatemodel$selected_predictors
  rm(climatemodel)
  
  
  #--------------------------------------------------
  #------------------Define validation types---------
  #--------------------------------------------------
  eu_occ <- global_presabs%>%
    dplyr::filter(species==1)%>%
    sf::st_filter(euboundary_wgs84) 
  
  #Only validate climate model in Europe if 40 or more occs 
  eu_climate_validation <- nrow(eu_occ) >= 40
  
  #Only validate ensemble model if habitat model could be fitted
  ensemble_validation <- file.exists(habitat_qs_file)
  
  
  #---------------------------------------------------------
  #- Select climate rasters used in 03_fit_climate_model.R -
  #---------------------------------------------------------
  climate_selection <- terra::subset(climate_stack,
                                     climate_predictors[climate_predictors %in%
                                                          names(climate_stack)])
  
  
  #-----------------------------------------------------------
  #- Obtain global climate subsample values for selected predictors
  #-----------------------------------------------------------
  #Load biomes
  wwf_eco_biome<-sf::st_read(biome_path) 
  
  # Keep only biome polygons that intersect at least one occurrence point
  global_presences<-dplyr::filter(global_presabs, species==1)
  sf::sf_use_s2(FALSE)
  has_occurrence <- lengths(sf::st_intersects(wwf_eco_biome, global_presences)) > 0
  wwf_ecoSub1 <- wwf_eco_biome[has_occurrence, ]
  rm(wwf_eco_biome)
  sf::sf_use_s2(TRUE)
  
  
  #Mask Chelsa layer with biomes with occurrences
  wwf_ecoSub1_ext<-terra::ext(wwf_ecoSub1) 
  wwf_ecoSub1_vector <- terra::vect(wwf_ecoSub1) 
  climate_sub <- terra::crop(climate_selection[[1]], wwf_ecoSub1_ext) 
  climate_sub <- terra::mask(climate_sub, wwf_ecoSub1_vector)
  
  #Extract a subsample of global pixels for Boyce calculation
  set.seed(728)
  global_subsample<- terra::spatSample(
    climate_sub,
    size = boyce_background_size, 
    method = "random", 
    na.rm = TRUE, #Ignore NA pixels
    as.points = TRUE) 
  
  # Extract climate data at global subsample points
  global_points<- terra::extract(climate_selection, global_subsample, ID = FALSE, xy = FALSE)%>%
    dplyr::mutate(ID = dplyr::row_number())
  
  #Clean up
  rm( wwf_ecoSub1, wwf_ecoSub1_ext, wwf_ecoSub1_vector,climate_sub, global_subsample)
  
  
  #-----------------------------------------------------------
  #- Obtain European climate subsample values for selected predictors
  #-----------------------------------------------------------
  if(eu_climate_validation || ensemble_validation){
    eu_points <- eu_climate_sub %>%
      dplyr::select(any_of(climate_predictors))%>%
      dplyr::mutate(ID = dplyr::row_number())
  }
  
  
  #-----------------------------------------------------------------
  #- Define if cross validation can be done and for how many folds -
  #-----------------------------------------------------------------
  
  #Default is 0 folds and no CV
  cv_folds <- 0L
  use_cv <- FALSE
  use_cv_climate_only <- FALSE
  
  #If ensemble validation is done, let EU data drive number of folds
  if(ensemble_validation){
    
    #Load presabs data of habitat model
    habitatmodel   <- qs::qread(habitat_qs_file)
    eu_presabs <- habitatmodel$eu_presabs
    n_pres_ensemble <- sum(eu_presabs$species == 1)
    rm(habitatmodel)
    
    if (n_pres_ensemble>= 40L) {
      use_cv <- TRUE
      cv_folds <- min(5L, floor(n_pres_ensemble / 20L))
    }else if (nrow(global_presences) >= 40L) {
      use_cv_climate_only <- TRUE
      cv_folds <- min(5L, floor(nrow(global_presences) / 20L))
      eu_presabs <- eu_presabs%>%
        dplyr::mutate(ID = dplyr::row_number())
      
    }else{
      eu_presabs <- eu_presabs%>%
        dplyr::mutate(ID = dplyr::row_number())
    }
    
  }else{
    
    if (nrow(global_presences) >= 40L) {
      use_cv <- TRUE
      cv_folds <- min(5L, floor(nrow(global_presences) / 20L))
    }else{
      global_presabs <- global_presabs%>%
        dplyr::mutate(ID = dplyr::row_number())
    }
  }
  
  
  #-----------------------------------------------------------------
  #            OPTION 1: SPATIAL CROSS VALIDATION
  #-----------------------------------------------------------------
  if (use_cv || use_cv_climate_only) {
    
    #---------------------------------
    #----- Create spatial folds-------
    #---------------------------------
    
    if(ensemble_validation && !use_cv_climate_only){
      
      #--------------------------------------------------------------
      #---Prepare combined dataset of global_presabs and eu_presabs
      #--------------------------------------------------------------
      # Combine data
      all_presabs<-eu_presabs%>%
        sf::st_transform(crs=sf::st_crs(global_presabs))%>%
        dplyr::mutate(decimalLatitude = sf::st_coordinates(.)[, "Y"],
                      decimalLongitude = sf::st_coordinates(.)[, "X"])%>%
        dplyr::bind_rows(global_presabs)
      
      #Remove duplicates
      all_presabs$cell <- terra::cellFromXY( climate_stack[[1]], 
                                             all_presabs%>%
                                               st_coordinates%>%
                                               as.data.frame()) 
      
      all_presabs <- all_presabs %>%
        dplyr::filter(!is.na(cell))%>%
        group_by(cell) %>%
        dplyr::distinct(cell, .keep_all=TRUE)%>%
        dplyr::ungroup()
      
      
      #-----------------------------
      #---Generate spatial folds
      #-----------------------------
      sf::sf_use_s2(FALSE)
      set.seed(123)
      sb <- blockCV::cv_spatial(
        x         = vect(all_presabs),
        column    = "species",
        k         = cv_folds,
        hexagon   = TRUE,
        selection = "random",
        iteration = 200,
        size      = 100000) #100 km
      
      sf::sf_use_s2(TRUE)
      fold_structure<-sb$blocks["folds"]
      
      
      #------------------------------------------
      #- Assign occs of ensemble model to folds -
      #------------------------------------------
      eu_presabs_perfold <- sf::st_join(sf::st_transform(eu_presabs, crs=sf::st_crs(global_presabs)),
                                        fold_structure,  
                                        join = sf::st_within,        
                                        left = TRUE)%>%
        dplyr::filter(!is.na(folds))%>%
        dplyr::mutate(ID = dplyr::row_number())
      
      if(nrow(eu_presabs_perfold)!=nrow(eu_presabs)){
        warning(nrow(eu_presabs)- nrow(eu_presabs_perfold)," Ensemble model point(s) not assigned to a fold and removed from dataset.")
      }
      
    }else{
      
      sf::sf_use_s2(FALSE)
      # Hex, class-balanced spatial folds
      set.seed(123)
      sb <- blockCV::cv_spatial(
        x         = vect(global_presabs),
        column    = "species",
        k         = cv_folds,
        hexagon   = TRUE,
        selection = "random",
        iteration = 200,
        size      = 100000) #100 km
      
      sf::sf_use_s2(TRUE)
      fold_structure<-sb$blocks["folds"]
    }
    
    
    #-----------------------------------------
    #- Assign occs of climate model to folds -
    #-----------------------------------------
    global_presabs_perfold <- sf::st_join(global_presabs,
                                          fold_structure,  
                                          join = sf::st_within,        
                                          left = TRUE)%>%
      dplyr::filter(!is.na(folds))%>%
      dplyr::mutate(ID = dplyr::row_number())
    
    if(nrow(global_presabs_perfold)!=nrow(global_presabs)){
      warning(nrow(global_presabs)- nrow(global_presabs_perfold)," global point(s) not assigned to a fold and removed from dataset.")
    }
    
    
    #-----------------------------------------
    #- Create lists for storing results -
    #-----------------------------------------
    global_validation_climate <- list()
    if(eu_climate_validation) eu_validation_climate <- list()
    median_favourability_climate_perfold <- vector("list", cv_folds)
    
    
    #----------------------------------------------------------
    #-- Fit models on each training set and predict test set --
    #----------------------------------------------------------
    #Start loop per fold
    for (fold in seq_len(cv_folds)) {
      
      message(sprintf("Creating climate validation metrics for fold %d/%d: use folds %s for training", 
                      fold, cv_folds, paste(seq_len(cv_folds)[-fold], collapse = ", ")))
      
      
      #--------------------------------------
      #-          Define train data         -
      #--------------------------------------
      #Create training dataset
      train_data  <- global_presabs_perfold%>%
        dplyr::filter(folds!=fold)
      
      # Prevalence ratio from training data
      pres_train <- sum(train_data$species == 1)
      abs_train  <- sum(train_data$species == 0)
      prev_ratio <- pres_train/abs_train
      
      
      #--------------------------------------
      #-      Fit models on train data      -
      #--------------------------------------
      #Prepare model framework
      sdm_data <- sdm::sdmData(
        species ~ .,
        train      = vect(train_data),
        predictors = climate_selection
      )
      
      #Fit models
      model <- sdm::sdm(species ~ ., data = sdm_data, methods = top5_methods)
      
      
      #-----------------------------------------------------------
      #--    Prepare datasets with climate data for predictions --
      #-----------------------------------------------------------
      #Extract data for global validation
      test_data  <- global_presabs_perfold%>%
        dplyr::filter(folds == fold)
      
      global_env <- extract_env(test_data, climate_selection)
      datasets <- list(global_points = global_points,
                       occ_env       = global_env$presences,
                       abs_env       = global_env$absences)
      
      #Extract data for validation in Europe
      if(eu_climate_validation){
        eu_test_data  <- test_data%>%
          sf::st_filter(euboundary_wgs84)
        
        eu_env<-extract_env(eu_test_data, climate_selection)
        datasets$eu_occ_env <- eu_env$presences
        datasets$eu_abs_env <- eu_env$absences
      }
      
      #Extract data for validation of ensemble model
      if(ensemble_validation & !use_cv_climate_only){
        ensemble_test_data<-eu_presabs_perfold%>%
          dplyr::filter(folds == fold)
        
        ensemble_env<-extract_env(ensemble_test_data, climate_selection)
        datasets$ens_occ_env <- ensemble_env$presences
        datasets$ens_abs_env <- ensemble_env$absences
      }
      
      #Add EU background data for validation of ensemble and Europe
      if (eu_climate_validation || ensemble_validation) {datasets$eu_points  <- eu_points}
      
      
      #---------------------------------------------------------------------
      #---- Make predictions per model algorithm and dataset and get median 
      #----------------------------------------------------------------------
      median_favourability_climate_perfold[[fold]]<- compute_median_favourability(model,
                                                                                  datasets,
                                                                                  top5_methods,
                                                                                  prev_ratio)
      
      #-----------------------------------------
      #------- Compute Boyce, AUC, and TSS -----
      #-----------------------------------------
      climate_fav<-median_favourability_climate_perfold[[fold]]
      
      #Global
      global_validation_climate[[fold]] <- compute_validation_metrics(
        species= speciesName,
        type = "Global_climate",
        fold = fold,
        all_suit_vals = climate_fav$global_points$median_favourability,
        occ_suit_vals = climate_fav$occ_env$median_favourability,
        abs_suit_vals = climate_fav$abs_env$median_favourability)
      
      #EU
      if(eu_climate_validation){
        eu_validation_climate[[fold]] <- compute_validation_metrics(
          species= speciesName,
          type = "Europe_climate",
          fold = fold,
          all_suit_vals = climate_fav$eu_points$median_favourability,
          occ_suit_vals = climate_fav$eu_occ_env$median_favourability,
          abs_suit_vals = climate_fav$eu_abs_env$median_favourability)
      }
      
      #Clean
      terra::tmpFiles(remove = TRUE)
    }
    
    
    #-----------------------------------------
    #---- Store validation metrics in dfs ----
    #-----------------------------------------
    #Set AUC and tss to NA in global validation as these are calculated for different regions per species and, hence, are not comparable
    global_validation_climate <- dplyr::bind_rows(global_validation_climate)
    if(eu_climate_validation){
      eu_validation_climate <- dplyr::bind_rows(eu_validation_climate)
    } 
    
    
  } else if (!use_cv) {
    
    #--------------------------------------------------
    #-          OPTION 2: NO CROSS VALIDATION
    #--------------------------------------------------
    
    #--------------------------------------
    #-      Fit models on full data      -
    #--------------------------------------
    #Prepare model framework
    sdm_data <- sdm::sdmData(
      species ~ .,
      train      = vect(global_presabs),
      predictors = climate_selection
    )
    
    #Fit models
    model <- sdm::sdm(species ~ ., data = sdm_data, methods = top5_methods)
    
    
    #--------------------------------------
    #------- Define prevalence ratio ------
    #--------------------------------------
    pres_total <- sum(global_presabs$species == 1)
    abs_total  <- sum(global_presabs$species == 0)
    prev_ratio <- pres_total / abs_total
    
    
    #-----------------------------------------------------------
    #---- Prepare datasets with climate data for predictions --
    #-----------------------------------------------------------
    datasets<-list()
    
    if(!use_cv_climate_only){
      #Extract data for global validation
      global_env<-extract_env(global_presabs, climate_selection)
      datasets <- list(global_points = global_points,
                       occ_env       = global_env$presences,
                       abs_env       = global_env$absences)
      
      #Extract data for validation in Europe
      if(eu_climate_validation){
        euboundary_presabs  <- global_presabs%>%
          sf::st_filter(euboundary_wgs84)
        
        eu_env<-extract_env(euboundary_presabs, climate_selection)
        datasets$eu_occ_env <- eu_env$presences
        datasets$eu_abs_env <- eu_env$absences
      }
    }
    
    #Extract data for validation of ensemble model
    if(ensemble_validation){
      ensemble_presabs<-eu_presabs%>%
        st_transform(crs=sf::st_crs(global_presabs))
      
      ensemble_env<-extract_env(ensemble_presabs, climate_selection)
      datasets$ens_occ_env <- ensemble_env$presences
      datasets$ens_abs_env <- ensemble_env$absences
    }
    
    #Add EU background data for validation of ensemble and Europe
    if (eu_climate_validation || ensemble_validation) {datasets$eu_points  <- eu_points}
    
    
    #-----------------------------------------------------------
    #---- Make predictions per model algorithm and dataset----
    #-----------------------------------------------------------
    median_fav_climate<- compute_median_favourability(model,
                                                      datasets,
                                                      top5_methods,
                                                      prev_ratio)
    
    
    #-----------------------------------------
    #------- Compute Boyce, AUC, and TSS -----
    #-----------------------------------------
    
    if(!use_cv_climate_only){
      message("Calculating validation metrics (no cross-validation)")
      
      #Global
      global_validation_climate <- compute_validation_metrics(
        species = speciesName,
        type =  "Global_climate",
        fold = "No cross-validation",
        all_suit_vals = median_fav_climate$global_points$median_favourability,
        occ_suit_vals = median_fav_climate$occ_env$median_favourability,
        abs_suit_vals = median_fav_climate$abs_env$median_favourability)
      
      
      #EU
      if(eu_climate_validation){
        eu_validation_climate <- compute_validation_metrics(
          species = speciesName,
          type =  "Europe_climate",
          fold = "No cross-validation",
          all_suit_vals = median_fav_climate$eu_points$median_favourability,
          occ_suit_vals = median_fav_climate$eu_occ_env$median_favourability,
          abs_suit_vals = median_fav_climate$eu_abs_env$median_favourability)
        
      }
    }
    
  }
  
  #--------------------------------------------
  #-------------- Export results --------------
  #--------------------------------------------
  # Define directories
  climate_validation_dir<-file.path(base_dir, "Climate", "Current", "Diagnostics", "Model_validation")
  if(!dir.exists(climate_validation_dir)) dir.create(climate_validation_dir, recursive = TRUE, showWarnings = FALSE)
  
  
  # Export validation summary (mean across folds) when relevant
  if(use_cv || use_cv_climate_only){
    
    #Export per fold validation
    readr::write_csv(global_validation_climate,
                     file.path(climate_validation_dir, paste0(speciesName, "_global_climate_validation_per_fold.csv"))) 
    
    #Export summary
    global_validation_clim_mean <- summarise_validation(df = global_validation_climate, 
                                                        validation = "Cross-validation")
    readr::write_csv(global_validation_clim_mean,
                     file.path(climate_validation_dir, paste0(speciesName, "_global_climate_validation_summary.csv"))) 
    
    #Bind results to validation summary
    species_validation_summary <- species_validation_summary %>%
      dplyr::bind_rows(global_validation_clim_mean)
    
    if (eu_climate_validation) {
      
      #Export per fold validation
      readr::write_csv(eu_validation_climate,
                       file.path(climate_validation_dir, paste0(speciesName, "_eu_climate_validation_per_fold.csv")))
      
      #Export summary
      eu_validation_clim_mean <- summarise_validation(eu_validation_climate, 
                                                      validation ="Cross-validation")
      readr::write_csv(eu_validation_clim_mean,
                       file.path(climate_validation_dir, paste0(speciesName, "_eu_climate_validation_summary.csv")))
      
      #Add summary to species validation df
      species_validation_summary<-species_validation_summary%>%
        dplyr::bind_rows(eu_validation_clim_mean)
      
    }
    
  }else{
    #Export non-cross-validated results
    readr::write_csv(global_validation_climate,
                     file.path(climate_validation_dir, paste0(speciesName, "_global_climate_validation_summary.csv")))
    
    #Add summary to species validation df
    global_validation_clim_mean<-summarise_validation(df = global_validation_climate,
                                                      validation = "No cross-validation")
    
    species_validation_summary<-species_validation_summary%>%
      dplyr::bind_rows(global_validation_clim_mean)
    
    if (eu_climate_validation) {
      #Export non crossvalidated results
      readr::write_csv(eu_validation_climate,
                       file.path(climate_validation_dir, paste0(speciesName, "_eu_climate_validation_summary.csv")))
      
      #Store summary in species validation df
      eu_validation_clim_mean<-summarise_validation(eu_validation_climate,
                                                    validation = "No cross-validation")
      
      species_validation_summary<-species_validation_summary%>%
        dplyr::bind_rows(eu_validation_clim_mean)
      
    }
  }
  
  
  
  #==============================================
  #=                                            =
  #=     PART 2: European landcover model       =
  #=                                            =
  #==============================================
  
  #--------------------------------------------
  #---- Should habitat validation be done? ----
  #--------------------------------------------
  if (!ensemble_validation) {
    warning("No habitat model was fitted for species ", species,
            "\n Skipping habitat model validation.")
    validation_summary[[speciesName]]<-  species_validation_summary
    next
    
  }
  
  
  #--------------------------------------------
  #--- Load  data stored in climate model qs --
  #--------------------------------------------
  habitatmodel   <- qs::qread(habitat_qs_file)
  eu_presabs <- habitatmodel$eu_presabs
  #TOACTIVATE
  #top5_methods  <- habitatmodel$top5_models
  #habitat_predictors <- habitatmodel$selected_predictors
  #TOREMOVE
  top5_habitat_methods  <- c("gam", "rf", "glmpoly", "mars", "maxent")
  habitat_predictors<-unique(habitatmodel[["varimp_df"]][["Predictor"]])
  rm(habitatmodel)
  
  
  #---------------------------------------------------------
  #- Select landcover rasters used in 04_fit_climate_model.R -
  #---------------------------------------------------------
  habitat_selection <- terra::subset(habitat_stack,
                                     habitat_predictors[habitat_predictors %in%
                                                          names(habitat_stack)])
  
  
  #-----------------------------------------------------------
  #- Obtain habitat subsample values for selected predictors
  #-----------------------------------------------------------
  eu_habitat_points<-eu_habitat_sub %>%
    dplyr::select(any_of(habitat_predictors))%>%
    dplyr::mutate(ID = dplyr::row_number())
  
  
  #-----------------------------------------------------------------
  #            OPTION 1: SPATIAL CROSS VALIDATION
  #-----------------------------------------------------------------
  if (use_cv) {
    
    #--------------------------------------------
    #--Put fold assignment data in right CRS ----
    #--------------------------------------------
    eu_presabs_perfold<- eu_presabs_perfold%>%
      sf::st_transform(crs=sf::st_crs(eu_presabs))
    
    
    #----------------------------------------------------------
    #-- Fit models on each training set and predict test set --
    #----------------------------------------------------------
    #Define lists to store validation metrics
    validation_habitat <- list()
    median_favourability_habitat_perfold <- vector("list", cv_folds)
    
    #Start loop per fold
    for (fold in seq_len(cv_folds)) {
      
      message(sprintf("Creating habitat validation metrics for fold %d/%d: use folds %s for training", 
                      fold, cv_folds, paste(seq_len(cv_folds)[-fold], collapse = ", ")))
      
      
      #--------------------------------------
      #-          Define train data         -
      #--------------------------------------
      #Create training dataset
      train_data  <- eu_presabs_perfold%>%
        dplyr::filter(folds!=fold)
      
      # Prevalence ratio from training data
      pres_train <- sum(train_data$species == 1)
      abs_train  <- sum(train_data$species == 0)
      prev_ratio <- pres_train/abs_train
      
      
      #--------------------------------------
      #-      Fit models on train data      -
      #--------------------------------------
      #Prepare model framework
      sdm_data <- sdm::sdmData(
        species ~ .,
        train      = vect(train_data),
        predictors = habitat_selection
      )
      
      #Fit models
      habitat_model <- sdm::sdm(species ~ ., data = sdm_data, methods = top5_habitat_methods)
      
      
      #-----------------------------------------------------------
      #---- Prepare datasets with habitat data for predictions --
      #-----------------------------------------------------------
      #Extract data for validation in Europe
      test_data  <- eu_presabs_perfold%>%
        dplyr::filter(folds == fold)
      
      europe_hab<-extract_env(test_data, habitat_selection)
      habitat_datasets <- list(eu_habitat_points = eu_habitat_points,
                               occ_hab       = europe_hab$presences,
                               abs_hab       = europe_hab$absences)
      
      
      #---------------------------------------------------------------------
      #---- Make predictions per model algorithm and dataset and get median 
      #----------------------------------------------------------------------
      median_favourability_habitat_perfold[[fold]]<- compute_median_favourability(habitat_model,
                                                                                  habitat_datasets,
                                                                                  top5_habitat_methods,
                                                                                  prev_ratio)
      
      
      #-----------------------------------------
      #------- Compute Boyce, AUC, and TSS -----
      #-----------------------------------------
      habitat_fav<-median_favourability_habitat_perfold[[fold]]
      
      #EU
      validation_habitat[[fold]] <- compute_validation_metrics(
        species= speciesName,
        type = "Europe_habitat",
        fold = fold,
        all_suit_vals = habitat_fav$eu_habitat_points$median_favourability,
        occ_suit_vals = habitat_fav$occ_hab$median_favourability,
        abs_suit_vals = habitat_fav$abs_hab$median_favourability)
      
      #Clean
      terra::tmpFiles(remove = TRUE)
    }
    
    
    #-----------------------------------------
    #---- Store validation metrics in dfs ----
    #-----------------------------------------
    eu_validation_habitat <- dplyr::bind_rows(validation_habitat)
    
    
  } else {
    
    #--------------------------------------------------
    #-          OPTION 2: NO CROSS VALIDATION
    #--------------------------------------------------
    
    #--------------------------------------
    #-      Fit models on full data      -
    #--------------------------------------
    #Prepare habitat_model framework
    sdm_data <- sdm::sdmData(
      species ~ .,
      train      = vect(eu_presabs),
      predictors = habitat_selection
    )
    
    #Fit models
    habitat_model <- sdm::sdm(species ~ ., data = sdm_data, methods = top5_habitat_methods)
    
    
    #--------------------------------------
    #------- Define prevalence ratio ------
    #--------------------------------------
    pres_total <- sum(eu_presabs$species == 1)
    abs_total  <- sum(eu_presabs$species == 0)
    prev_ratio <- pres_total / abs_total
    
    
    #-----------------------------------------------------------
    #---- Prepare datasets with habitat data for predictions --
    #-----------------------------------------------------------
    #Extract data for validation in Europe
    europe_hab<-extract_env(eu_presabs, habitat_selection)
    habitat_datasets <- list(eu_habitat_points = eu_habitat_points,
                             occ_hab       = europe_hab$presences,
                             abs_hab       = europe_hab$absences)
    
    
    #---------------------------------------------------------------------
    #---- Make predictions per model algorithm and dataset and get median 
    #----------------------------------------------------------------------
    median_fav_habitat<- compute_median_favourability(habitat_model,
                                                      habitat_datasets,
                                                      top5_habitat_methods,
                                                      prev_ratio)
    
    
    #-----------------------------------------
    #------- Compute Boyce, AUC, and TSS -----
    #-----------------------------------------
    eu_validation_habitat <- compute_validation_metrics(
      species= speciesName,
      type = "Europe_habitat",
      fold = "No cross-validation",
      all_suit_vals = median_fav_habitat$eu_habitat_points$median_favourability,
      occ_suit_vals = median_fav_habitat$occ_hab$median_favourability,
      abs_suit_vals = median_fav_habitat$abs_hab$median_favourability)
    
    
  }
  
  #--------------------------------------------
  #-------------- Export results --------------
  #--------------------------------------------
  # Export validation overview
  habitat_validation_dir<-file.path(base_dir, "Habitat", "Current", "Diagnostics", "Model_validation")
  if(!dir.exists(habitat_validation_dir)) dir.create(habitat_validation_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Export validation summary (mean across folds) when relevant
  if(use_cv){
    
    #Export per fold validation metrics
    readr::write_csv(eu_validation_habitat,
                     file.path(habitat_validation_dir, paste0(speciesName, "_habitat_validation_per_fold.csv")))
    
    #Export summary validation metrics
    validation_hab_mean<-summarise_validation(eu_validation_habitat,
                                              validation = "Cross-validation")
    readr::write_csv(validation_hab_mean,
                     file.path(habitat_validation_dir, paste0(speciesName, "_habitat_validation_summary.csv")))
    
  }else{
    validation_hab_mean<-summarise_validation(eu_validation_habitat,
                                              validation = "No cross-validation")
    
    readr::write_csv(eu_validation_habitat,
                     file.path(habitat_validation_dir, paste0(speciesName, "_habitat_validation_summary.csv")))
  }
  
  species_validation_summary<-species_validation_summary%>%
    dplyr::bind_rows(validation_hab_mean)
  
  
  
  #==============================================
  #=                                            =
  #=         PART 3: Ensemble validation        =
  #=                                            =
  #==============================================
  
  
  #-----------------------------------------------------------------
  #            OPTION 1: SPATIAL CROSS VALIDATION
  #-----------------------------------------------------------------
  if (use_cv) {
    
    #----------------------------------------------------------
    #-- Combine predictions of habitat and climate model --
    #----------------------------------------------------------
    #Define lists to store validation metrics
    validation_ensemble<- list()
    
    #Start loop per fold
    for (fold in seq_len(cv_folds)) {
      
      message(sprintf("Calculating ensemble validation metrics for test fold %d/%d", fold,cv_folds))
      
      # Extract median favourability for the current fold
      hab_fav <- median_favourability_habitat_perfold[[fold]]
      clim_fav <- median_favourability_climate_perfold[[fold]]
      
      #Generate ensemble favourability for background points, occs, and abs.
      ensemble_background_fav <- ensemble_geom_mean(hab_fav$eu_habitat_points, clim_fav$eu_points, type="background")
      ensemble_occ_fav <- ensemble_geom_mean(hab_fav$occ_hab,clim_fav$ens_occ_env, type="occurrence")
      ensemble_abs_fav <- ensemble_geom_mean(hab_fav$abs_hab,clim_fav$ens_abs_env, type="absence")
      
      
      #-----------------------------------------
      #------- Compute Boyce, AUC, and TSS -----
      #-----------------------------------------
      validation_ensemble[[fold]] <- compute_validation_metrics(
        species= speciesName,
        type = "Europe_ensemble",
        fold = fold,
        all_suit_vals = ensemble_background_fav,
        occ_suit_vals =  ensemble_occ_fav ,
        abs_suit_vals = ensemble_abs_fav)
    }
    
    #-----------------------------------------
    #---- Store validation metrics in df ----
    #-----------------------------------------
    eu_validation_ensemble <- dplyr::bind_rows(validation_ensemble)
    
    
  } else{
    
    #--------------------------------------------------
    #-          OPTION 2: NO CROSS VALIDATION
    #--------------------------------------------------
    
    # Extract median favourability for habitat and climate
    hab_fav <- median_favourability_habitat
    clim_fav <- median_favourability_climate
    
    
    #Generate ensemble favourability for background points, occs, and abs.
    ensemble_background_fav <- ensemble_geom_mean(hab_fav$eu_habitat_points, clim_fav$eu_points, type="background")
    ensemble_occ_fav <- ensemble_geom_mean(hab_fav$occ_hab,clim_fav$ens_occ_env, type="occurrence")
    ensemble_abs_fav <- ensemble_geom_mean(hab_fav$abs_hab,clim_fav$ens_abs_env, type="absence")
    
    
    #-----------------------------------------
    #------- Compute Boyce, AUC, and TSS -----
    #-----------------------------------------
    message("Calculating ensemble validation metrics (no cross-validation)")
    
    eu_validation_ensemble <- compute_validation_metrics(
      species= speciesName,
      type = "Europe_ensemble",
      fold = "No cross-validation",
      all_suit_vals = ensemble_background_fav,
      occ_suit_vals =  ensemble_occ_fav ,
      abs_suit_vals = ensemble_abs_fav)
    
  }
  
  
  #--------------------------------------------
  #-------------- Export results --------------
  #--------------------------------------------
  # Define directory
  ensemble_validation_dir<-file.path(base_dir, "Combined", "Current", "Diagnostics", "Model_validation")
  if(!dir.exists(ensemble_validation_dir)) dir.create(ensemble_validation_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Export validation summary (mean across folds) when relevant
  if(use_cv){
    
    #Export per fold validation metrics
    readr::write_csv(eu_validation_ensemble,
                     file.path(ensemble_validation_dir, paste0(speciesName, "_combined_validation_per_fold.csv")))
    
    #Export summary validation metrics
    validation_ens_mean<- summarise_validation(eu_validation_ensemble,
                                               validation = "Cross-validation")
    readr::write_csv(validation_ens_mean,
                     file.path(ensemble_validation_dir, paste0(speciesName, "_combined_validation_summary.csv")))
    
  }else{
    
    #Export summary validation metrics
    validation_ens_mean<-summarise_validation(eu_validation_ensemble,
                                              validation="No cross-validation")
    readr::write_csv(eu_validation_ensemble,
                     file.path(ensemble_validation_dir, paste0(speciesName, "_combined_validation_summary.csv")))
    
  }
  
  species_validation_summary<-species_validation_summary%>%
    dplyr::bind_rows(validation_ens_mean)
  
  #--------------------------------------------
  #-----Store results in validation summary ---
  #--------------------------------------------
  validation_summary[[speciesName]]<-  species_validation_summary
  
}

#--------------------------------------------
#----- Export combined validation results--------------
#--------------------------------------------
final_validation <- bind_rows(validation_summary)
readr::write_csv(final_validation,
                 file.path(validation_dir, "Validation_summary.csv"))
