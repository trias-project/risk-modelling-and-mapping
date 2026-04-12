#-------------------------------------------------------------------------------
#--------------------------    Load packages      ------------------------------
#-------------------------------------------------------------------------------
packages <- c( "dplyr", "qs", "terra", "tidyterra", "sf", "here",
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
#Extract a subsample of global pixels for Boyce calculation: around 40min-1h!
set.seed(728)
global_subsample<- terra::spatSample(
  climate_stack[[1]],
  size = boyce_background_size, 
  method = "random", 
  na.rm = TRUE, #Ignore NA pixels
  as.points = TRUE) 

#Extract a subsample of European pixels for Boyce calculation 
set.seed(728)
eu_subsample <- terra::spatSample(
  habitat_stack[[1]],
  size = boyce_background_size, 
  method = "random", 
  na.rm = TRUE, #Ignore NA pixels
  as.points = TRUE)

# Extract climate data at global subsample points
global_climate_sub <- terra::extract(climate_stack, global_subsample, ID = FALSE, xy = FALSE)

# Extract climate data at eu subsample points
eu_climate_sub <- terra::extract(eu_climate_stack, eu_subsample, ID = FALSE, xy = FALSE)

# Extract habitat data at eu subsample points
eu_habitat_sub <- terra::extract(habitat_stack, eu_subsample, ID = FALSE, xy = FALSE)


#----------------------------------------------
#----------- Define results directory ---------
#----------------------------------------------
results_dir <- file.path("data", "projects", project, "Model_validation_results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


#--------------------------------------------
#----------- Start species loop -------------
#--------------------------------------------
for (i in seq_along(accepted_taxonkeys)) {
  
  
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
  #==============================================
  #=           PART 1: Climate model            =
  #==============================================
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
    warning("No climate model could be fitted for ",species,
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
    st_transform(crs = st_crs(euboundary)) %>%
  
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
  eu_climate_selection  <- terra::subset(eu_climate_stack,
                                         climate_predictors[climate_predictors %in%
                                                              names(eu_climate_stack)])
  
  
  #-----------------------------------------------------------
  #- Obtain climate subsample values for selected predictors
  #-----------------------------------------------------------
  global_points<-global_climate_sub %>%
    dplyr::select(any_of(climate_predictors))
  
  if(predict.eu){
    eu_points <- eu_climate_sub %>%
      dplyr::select(any_of(climate_predictors))
  }
  
  
  #-----------------------------------------------------------------
  #- Define if cross validation can be done and for how many folds -
  #-----------------------------------------------------------------
  n_pres <- sum(global_presabs$species == 1)
  k <- 0L
  use_cv <- FALSE
  
  if (n_pres >= 40L) {
    k <- min(5L, floor(n_pres / 20L))
    use_cv <- k >= 2L
  }
  
  
  #-----------------------------------------------------------------
  #            OPTION 1: SPATIAL CROSS VALIDATION
  #-----------------------------------------------------------------
  if (use_cv) {
    
    #---------------------------------
    #----- Create spatial folds-------
    #---------------------------------
    sf::sf_use_s2(FALSE)
    # Hex, class-balanced spatial folds
    sb <- blockCV::cv_spatial(
      x         = vect(global_presabs),
      column    = "species",
      r         = climate_selection,
      k         = k,
      hexagon   = TRUE,
      selection = "random",
      iteration = 200,
      size      = 100000 #100 km
    )
    sf::sf_use_s2(TRUE)
    
    fold_ids <- sb$folds_ids
    stopifnot(length(fold_ids) == nrow(global_presabs))
    
    # Initiate per-fold list
    fold_fav  <- vector("list", k)  
    
    
    #----------------------------------------------------------
    #-- Fit models on each training set and predict test set --
    #----------------------------------------------------------
    #Define lists to store validation metrics
    global_validation_climate <- list()
    if(predict.eu) eu_validation_climate <- list()
    
    #Start loop per fold
    for (fold in seq_len(k)) {
      message(sprintf("Creating validation metrics for fold %d/%d: use folds %s for training", fold, k, paste(seq_len(k)[-fold], collapse = ", ")))
      
      #--------------------------------------
      #-          Define train data         -
      #--------------------------------------
      #Create training dataset
      train_idx <- which(fold_ids != fold)
      train_data  <- global_presabs[train_idx, ]
      
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
      
      
      #---------------------------------------------------
      #- Get climate values at pres/abs in EU and global -
      #---------------------------------------------------
      #Define test presences and absences
      test_idx <- which(fold_ids == fold)
      test_data  <- global_presabs[test_idx, ]
      test_presences <- test_data[test_data$species==1,]
      test_absences <- test_data[test_data$species==0,]
      
      #Extract global climate data at presences and absences
      occ_env <- terra::extract(climate_selection, terra::vect(test_presences), ID = FALSE, xy = FALSE)
      abs_env <- terra::extract(climate_selection, terra::vect(test_absences), ID = FALSE, xy = FALSE)
      occ_env <- occ_env[complete.cases(occ_env), ]
      abs_env <- abs_env[complete.cases(abs_env), ]
      
      #Extract EU climate data at presences and absences
      if(predict.eu){
        eu_occ_env <- terra::extract(eu_climate_selection, terra::vect(test_presences), ID = FALSE, xy = FALSE)
        eu_abs_env <- terra::extract(eu_climate_selection, terra::vect(test_absences), ID = FALSE, xy = FALSE)
        eu_occ_env <- eu_occ_env[complete.cases(eu_occ_env), ]
        eu_abs_env <- eu_abs_env[complete.cases(eu_abs_env), ]
      }
      
      
      #-----------------------------------------------------------
      #---------------- List datasets for predictions-------------
      #-----------------------------------------------------------
      if(predict.eu){
        datasets <- list(global_points = global_points,
                         occ_env       = occ_env,
                         abs_env       = abs_env,
                         eu_points     = eu_points,
                         eu_occ_env    = eu_occ_env,
                         eu_abs_env    = eu_abs_env)
      }else{
        datasets <- list(global_points = global_points,
                         occ_env       = occ_env,
                         abs_env       = abs_env)
      }
      
      
      #-----------------------------------------------------------
      #---- Make predictions per model algorithm and dataset----
      #-----------------------------------------------------------
      favourability_pred <- list()
      for(modelmethod in top5_methods){
        
        message("Predicting for method: ", modelmethod,".")
        
        for(dataset_name in names(datasets)) {
          
          #Load datasets
          dataset <- datasets[[dataset_name]]
          
          #Predict for dataset
          dataset_suit <- predict(model,
                                  newdata = dataset,
                                  method = modelmethod)
          
          #Convert suitability to favourability
          dataset_fav<- favourability_from_prob(dataset_suit[[1]], prev_ratio)
          
          #Store in list
          favourability_pred[[modelmethod]][[dataset_name]] <- dataset_fav
          
          #Clean up
          rm(dataset_suit, dataset_fav)
          gc()
        }
      }  
      
      
      #-----------------------------------------
      #---- Calculate median favourability  ----
      #-----------------------------------------
      #Global
      global_bg_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "global_points")),
                                      1, median, na.rm = TRUE) 
      global_pres_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "occ_env")),
                                        1, median, na.rm = TRUE) 
      global_abs_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "abs_env")),
                                       1, median, na.rm = TRUE) 
      
      #Europe
      if(predict.eu){
        eu_bg_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_points")),
                                    1, median, na.rm = TRUE) 
        eu_pres_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_occ_env")),
                                      1, median, na.rm = TRUE) 
        eu_abs_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_abs_env")),
                                     1, median, na.rm = TRUE) 
      }
      
      
      #-----------------------------------------
      #------- Compute Boyce, AUC, and TSS -----
      #-----------------------------------------
      message(sprintf("Calculating validation metrics for test fold %d/%d", fold,k))
      
      #Global
      global_validation_climate[[fold]] <- compute_validation_metrics(
        fold = fold,
        all_suit_vals = global_bg_favourability,
        occ_suit_vals = global_pres_favourability,
        abs_suit_vals = global_abs_favourability)
      
      #EU
      if(predict.eu){
        eu_validation_climate[[fold]] <- compute_validation_metrics(
          fold = fold,
          all_suit_vals = eu_bg_favourability,
          occ_suit_vals = eu_pres_favourability,
          abs_suit_vals = eu_abs_favourability)
      }
      
      #Clean
      terra::tmpFiles(remove = TRUE)
    }
    
    
    #-----------------------------------------
    #---- Store validation metrics in dfs ----
    #-----------------------------------------
    global_validation_climate <- as.data.frame(do.call(rbind, global_validation_climate))
    if(predict.eu) eu_validation_climate <- as.data.frame(do.call(rbind, eu_validation_climate))
    
    
  } else {
    
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
    
    
    #--------------------------------------
    #-        Extract climate data        -
    #--------------------------------------
    #Define presences and absences
    presences <- global_presabs[global_presabs$species==1,]
    absences <- global_presabs[global_presabs$species==0,]
    
    #Extract global climate data at presences and absences
    occ_env <- terra::extract(climate_selection, terra::vect(presences), ID = FALSE, xy = FALSE)
    abs_env <- terra::extract(climate_selection, terra::vect(absences), ID = FALSE, xy = FALSE)
    occ_env <- occ_env[complete.cases(occ_env), ]
    abs_env <- abs_env[complete.cases(abs_env), ]
    
    #Extract EU climate data at presences and absences
    if(predict.eu){
      eu_occ_env <- terra::extract(eu_climate_selection, terra::vect(presences), ID = FALSE, xy = FALSE)
      eu_abs_env <- terra::extract(eu_climate_selection, terra::vect(absences), ID = FALSE, xy = FALSE)
      eu_occ_env <- eu_occ_env[complete.cases(eu_occ_env), ]
      eu_abs_env <- eu_abs_env[complete.cases(eu_abs_env), ]
    }
    
    
    #-----------------------------------------------------------
    #---------------- List datasets for predictions-------------
    #-----------------------------------------------------------
    if(predict.eu){
      datasets <- list(global_points = global_points,
                       occ_env       = occ_env,
                       abs_env       = abs_env,
                       eu_points     = eu_points,
                       eu_occ_env    = eu_occ_env,
                       eu_abs_env    = eu_abs_env)
    }else{
      datasets <- list(global_points = global_points,
                       occ_env       = occ_env,
                       abs_env       = abs_env)
    }
    
    
    #-----------------------------------------------------------
    #---- Make predictions per model algorithm and dataset----
    #-----------------------------------------------------------
    favourability_pred <- list()
    for(modelmethod in top5_methods){
      
      message("Predicting for method: ", modelmethod,".")
      
      for(dataset_name in names(datasets)) {
        
        #Load datasets
        dataset <- datasets[[dataset_name]]
        
        #Predict for dataset
        dataset_suit <- predict(model,
                                newdata = dataset,
                                method = modelmethod)
        
        #Convert suitability to favourability
        dataset_fav<- favourability_from_prob(dataset_suit[[1]], prev_ratio)
        
        #Store in list
        favourability_pred[[modelmethod]][[dataset_name]] <- dataset_fav
        
        #Clean up
        rm(dataset_suit, dataset_fav)
        gc()
        
      }
    }  
    
    
    #-----------------------------------------
    #---- Calculate median favourability  ----
    #-----------------------------------------
    #Global
    global_bg_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "global_points")),
                                    1, median, na.rm = TRUE) 
    global_pres_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "occ_env")),
                                      1, median, na.rm = TRUE) 
    global_abs_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "abs_env")),
                                     1, median, na.rm = TRUE) 
    
    #Europe
    if(predict.eu){
      eu_bg_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_points")),
                                  1, median, na.rm = TRUE) 
      eu_pres_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_occ_env")),
                                    1, median, na.rm = TRUE) 
      eu_abs_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_abs_env")),
                                   1, median, na.rm = TRUE) 
    }
    
    
    #-----------------------------------------
    #------- Compute Boyce, AUC, and TSS -----
    #-----------------------------------------
    message("Calculating validation metrics (no cross-validation)")
    
    #Global
    global_validation_climate <- compute_validation_metrics(
      fold = "No cross-validation",
      all_suit_vals = global_bg_favourability,
      occ_suit_vals = global_pres_favourability,
      abs_suit_vals = global_abs_favourability)
    
    #EU
    if(predict.eu){
      eu_validation_climate <- compute_validation_metrics(
        fold = "No cross-validation",
        all_suit_vals = eu_bg_favourability,
        occ_suit_vals = eu_pres_favourability,
        abs_suit_vals = eu_abs_favourability)
    }
    
  }
  
  
  #==============================================
  #=                                            =
  #=     PART 2: European landcover model       =
  #=                                            =
  #==============================================
  
  #--------------------------------------------
  #--Check if a habitat model has been fitted -
  #--------------------------------------------
  habitat_qs_file <- file.path(base_dir, "Habitat",
                               paste0("Habitat_model_", speciesName, "_", taxonkey, ".qs"))
  
  if (!file.exists(habitat_qs_file)) {
    warning("No habitat model was fitted for species ", species,
            "\n Skipping habitat model validation.")
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
  top5_methods  <- c("gam", "rf", "glmpoly", "mars", "maxent")
  habitat_predictors<-unique(habitatmodel[["varimp_df"]][["Predictor"]])
  rm(habitatmodel)
  
  
  #--------------------------------------------
  #--------------Define target CRS ------------
  #--------------------------------------------
  target_crs    <- sf::st_crs(terra::crs(habitat_stack))
  
  
  #---------------------------------------------------------
  #- Select landcover rasters used in 04_fit_climate_model.R -
  #---------------------------------------------------------
  habitat_selection <- terra::subset(habitat_stack,
                                     habitat_predictors[habitat_predictors %in%
                                                          names(habitat_stack)])
  
  
  #-----------------------------------------------------------
  #- Obtain habitat subsample values for selected predictors
  #-----------------------------------------------------------
  eu_points<-eu_habitat_sub %>%
    dplyr::select(any_of(habitat_predictors))
  
  
  #-----------------------------------------------------------------
  #- Define if cross validation can be done and for how many folds -
  #-----------------------------------------------------------------
  n_pres <- sum(eu_presabs$species == 1)
  k <- 0L
  use_cv <- FALSE
  
  if (n_pres >= 40L) {
    k <- min(5L, floor(n_pres / 20L))
    use_cv <- k >= 2L
  }
  
  
  #-----------------------------------------------------------------
  #            OPTION 1: SPATIAL CROSS VALIDATION
  #-----------------------------------------------------------------
  if (use_cv) {
    
    #---------------------------------
    #----- Create spatial folds-------
    #---------------------------------
    sf::sf_use_s2(FALSE)
    # Hex, class-balanced spatial folds
    sb <- blockCV::cv_spatial(
      x         = vect(eu_presabs),
      column    = "species",
      r         = habitat_selection,
      k         = k,
      hexagon   = TRUE,
      selection = "random",
      iteration = 200,
      size      = 100000 #100 km
    )
    sf::sf_use_s2(TRUE)
    
    fold_ids <- sb$folds_ids
    stopifnot(length(fold_ids) == nrow(eu_presabs))
    
    # Initiate per-fold list
    fold_fav  <- vector("list", k)  
    
    
    #----------------------------------------------------------
    #-- Fit models on each training set and predict test set --
    #----------------------------------------------------------
    #Define lists to store validation metrics
    validation_habitat <- list()
    
    #Start loop per fold
    for (fold in seq_len(k)) {
      message(sprintf("Creating validation metrics (habitat) for fold %d/%d: use folds %s for training", fold, k, paste(seq_len(k)[-fold], collapse = ", ")))
      
      #--------------------------------------
      #-          Define train data         -
      #--------------------------------------
      #Create training dataset
      train_idx <- which(fold_ids != fold)
      train_data  <- eu_presabs[train_idx, ]
      
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
      habitat_model <- sdm::sdm(species ~ ., data = sdm_data, methods = top5_methods)
      
      
      #---------------------------------------------------
      #------ Get habitat values at pres/abs in EU  ------
      #---------------------------------------------------
      #Define test presences and absences
      test_idx <- which(fold_ids == fold)
      test_data  <- eu_presabs[test_idx, ]
      test_presences <- test_data[test_data$species==1,]
      test_absences <- test_data[test_data$species==0,]
      
      #Extract EU habitat data at presences and absences
      occ_habitat <- terra::extract(habitat_selection, terra::vect(test_presences), ID = FALSE, xy = FALSE)
      abs_habitat <- terra::extract(habitat_selection, terra::vect(test_absences), ID = FALSE, xy = FALSE)
      occ_habitat <- occ_habitat[complete.cases(occ_habitat), ]
      abs_habitat <- abs_habitat[complete.cases(abs_habitat), ]
      
      
      #-----------------------------------------------------------
      #---------------- List datasets for predictions-------------
      #-----------------------------------------------------------
      datasets <- list(eu_points   = eu_points,
                       occ_habitat = occ_habitat,
                       abs_habitat = abs_habitat)
      
      
      
      #-----------------------------------------------------------
      #---- Make predictions per model algorithm and dataset----
      #-----------------------------------------------------------
      favourability_pred <- list()
      for(modelmethod in top5_methods){
        
        message("Predicting for method: ", modelmethod,".")
        
        for(dataset_name in names(datasets)) {
          
          #Load datasets
          dataset <- datasets[[dataset_name]]
          
          #Predict for dataset
          dataset_suit <- predict(habitat_model,
                                  newdata = dataset,
                                  method = modelmethod)
          
          #Convert suitability to favourability
          dataset_fav<- favourability_from_prob(dataset_suit[[1]], prev_ratio)
          
          #Store in list
          favourability_pred[[modelmethod]][[dataset_name]] <- dataset_fav
          
          #Clean up
          rm(dataset_suit, dataset_fav)
          gc()
        }
      }  
      
      
      #-----------------------------------------
      #---- Calculate median favourability  ----
      #-----------------------------------------
      habitat_bg_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_points")),
                                       1, median, na.rm = TRUE) 
      habitat_pres_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "occ_habitat")),
                                         1, median, na.rm = TRUE) 
      habitat_abs_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "abs_habitat")),
                                        1, median, na.rm = TRUE) 
      
      
      
      #-----------------------------------------
      #------- Compute Boyce, AUC, and TSS -----
      #-----------------------------------------
      message(sprintf("Calculating habitat validation metrics for test fold %d/%d", fold,k))
      
      #EU
      validation_habitat[[fold]] <- compute_validation_metrics(
        fold = fold,
        all_suit_vals = habitat_bg_favourability,
        occ_suit_vals = habitat_pres_favourability,
        abs_suit_vals = habitat_abs_favourability)
      
      #Clean
      terra::tmpFiles(remove = TRUE)
    }
    
    
    #-----------------------------------------
    #---- Store validation metrics in dfs ----
    #-----------------------------------------
    eu_validation_habitat <- as.data.frame(do.call(rbind, validation_habitat))
    
    
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
    habitat_model <- sdm::sdm(species ~ ., data = sdm_data, methods = top5_methods)
    
    
    #--------------------------------------
    #------- Define prevalence ratio ------
    #--------------------------------------
    pres_total <- sum(eu_presabs$species == 1)
    abs_total  <- sum(eu_presabs$species == 0)
    prev_ratio <- pres_total / abs_total
    
    
    #--------------------------------------
    #-        Extract climate data        -
    #--------------------------------------
    #Define presences and absences
    presences <- eu_presabs[eu_presabs$species==1,]
    absences <- eu_presabs[eu_presabs$species==0,]
    
    #Extract global climate data at presences and absences
    occ_habitat <- terra::extract(habitat_selection, terra::vect(presences), ID = FALSE, xy = FALSE)
    abs_habitat <- terra::extract(habitat_selection, terra::vect(absences), ID = FALSE, xy = FALSE)
    occ_habitat <- occ_habitat[complete.cases(occ_habitat), ]
    abs_habitat <- abs_habitat[complete.cases(abs_habitat), ]
    
    
    #-----------------------------------------------------------
    #---------------- List datasets for predictions-------------
    #-----------------------------------------------------------
    datasets <- list(eu_points = eu_points,
                     occ_habitat       = occ_habitat,
                     abs_habitat       = abs_habitat)
    
    
    #-----------------------------------------------------------
    #---- Make predictions per model algorithm and dataset----
    #-----------------------------------------------------------
    favourability_pred <- list()
    for(modelmethod in top5_methods){
      
      message("Predicting for method: ", modelmethod,".")
      
      for(dataset_name in names(datasets)) {
        
        #Load datasets
        dataset <- datasets[[dataset_name]]
        
        #Predict for dataset
        dataset_suit <- predict(habitat_model,
                                newdata = dataset,
                                method = modelmethod)
        
        #Convert suitability to favourability
        dataset_fav<- favourability_from_prob(dataset_suit[[1]], prev_ratio)
        
        #Store in list
        favourability_pred[[modelmethod]][[dataset_name]] <- dataset_fav
        
        #Clean up
        rm(dataset_suit, dataset_fav)
        gc()
        
      }
    }  
    
    
    #-----------------------------------------
    #---- Calculate median favourability  ----
    #-----------------------------------------
    habitat_bg_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "eu_points")),
                                     1, median, na.rm = TRUE) 
    habitat_pres_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "occ_habitat")),
                                       1, median, na.rm = TRUE) 
    habitat_abs_favourability = apply(do.call(cbind, lapply(favourability_pred, `[[`, "abs_habitat")),
                                      1, median, na.rm = TRUE) 
    
    
    #-----------------------------------------
    #------- Compute Boyce, AUC, and TSS -----
    #-----------------------------------------
    eu_validation_habitat <- compute_validation_metrics(
      fold = "No cross-validation",
      all_suit_vals = habitat_bg_favourability,
      occ_suit_vals = habitat_pres_favourability,
      abs_suit_vals = habitat_abs_favourability)
    
  }
  
    #--------------------------------------------------
    #-          OPTION 2: NO CROSS VALIDATION
    #--------------------------------------------------
    
    # Extract median favourability for habitat and climate
    hab_fav <- median_favourability_habitat
    clim_fav <- median_favourability_climate
    
    #Generate ensemble favourability for background points, occs, and abs.
    ensemble_background_fav <- sqrt(hab_fav$eu_habitat_points*clim_fav$eu_points)
    ensemble_occ_fav <-sqrt(hab_fav$occ_hab*clim_fav$ens_occ_env)
    ensemble_abs_fav <- sqrt(hab_fav$abs_hab*clim_fav$ens_abs_env)
    
    
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
    

}