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
climate_path <- file.path("data", "external", "climate", "chelsa_current",
                          "processed", "globalclimpreds.tif")

#EU climate stack
eu_climpreds_path <- file.path("data", "external", "climate", "chelsa_current",
                               "processed","euclimpreds.tif")

#habitat stack
habitat_path <- file.path("data", "external", "habitat", "processed", 
                          "habitat_stack.tif")

#Biome file path
biome_path<-file.path("data", "external", "GIS", "official", "newRealms.shp")

#--------------------------------------------
#------------- Load species data ------------
#--------------------------------------------
#Get taxa info
taxa_info <- read.csv2(file.path("data", "projects", project, 
                                 paste0(project, "_taxa_info.csv")))

#Select unique taxonkeys
accepted_taxonkeys <- unique(taxa_info$acceptedTaxonKey)


#--------------------------------------------
#-----------------Load rasters---------------
#--------------------------------------------
#Load rasters
habitat_stack <- terra::rast(habitat_path)


#--------------------------------------------
#------ Build default validation context ----
#--------------------------------------------
default_euboundary <- load_eu_boundary(
  custom_path = custom_eu_boundary_path,
  reference = habitat_stack[[1]]
)
default_euboundary_vect <- terra::vect(default_euboundary)

if (is.null(custom_eu_boundary_path)) {
  default_eu_sampling_mask <- habitat_stack[[1]]
} else {
  default_eu_sampling_mask <- habitat_stack[[1]] %>%
    terra::crop(default_euboundary_vect) %>%
    terra::mask(default_euboundary_vect)
}

default_eu_sampling_cells <- terra::global(!is.na(default_eu_sampling_mask), 
                                           "sum", na.rm = TRUE)[[1]]

if (is.na(default_eu_sampling_cells) || default_eu_sampling_cells == 0) {
  stop("No non-NA European habitat cells are available 
       for default validation sampling.", call. = FALSE)
}
default_climate_stack <- NULL
default_eu_climate_stack <- NULL
default_eu_subsample <- NULL
default_eu_climate_sub <- NULL
default_eu_habitat_sub <- NULL


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
  climate_input_mode_saved <- climatemodel$climate_input_mode
  
  if (is.null(climate_input_mode_saved) || !nzchar(climate_input_mode_saved)) {
    climate_input_mode_saved <- "chelsa"
  }
  
  climate_manifest_path_saved <- climatemodel$climate_manifest_path
  predictor_crs_saved <- climatemodel$predictor_crs
  
  if (!climate_input_mode_saved %in% c("chelsa", "user_specific")) {
    stop(
      "Unsupported climate input mode stored in ", basename(climate_qs_file),
      ": ", climate_input_mode_saved, call. = FALSE
    )
  }
  
  ensemble_validation <- file.exists(habitat_qs_file)
  habitatmodel <- NULL
  habitat_stack_active <- habitat_stack
  landcover_input_mode_saved <- "default"
  landcover_manifest_path_saved <- NULL
  landcover_predictor_crs_saved <- NULL
  
  if (ensemble_validation) {
    habitatmodel <- qs::qread(habitat_qs_file)
    landcover_input_mode_saved <- habitatmodel$landcover_input_mode
    
    if (is.null(landcover_input_mode_saved) || !nzchar(landcover_input_mode_saved)) {
      landcover_input_mode_saved <- "default"
    }
    
    landcover_manifest_path_saved <- habitatmodel$landcover_manifest_path
    landcover_predictor_crs_saved <- habitatmodel$landcover_predictor_crs
    
    if (!landcover_input_mode_saved %in% c("default", "user_specific")) {
      stop(
        "Unsupported land-cover input mode stored in ",
        basename(habitat_qs_file),
        ": ",
        landcover_input_mode_saved,
        call. = FALSE
      )
    }
    
    if (landcover_input_mode_saved == "user_specific") {
      if (is.null(landcover_manifest_path_saved) || !nzchar(landcover_manifest_path_saved)) {
        stop(
          "The saved habitat model metadata is missing 'landcover_manifest_path'. Rerun 04_fit_habitat_model.R for ",
          speciesName,
          ".",
          call. = FALSE
        )
      }
      
      if (is.null(landcover_predictor_crs_saved) || !nzchar(landcover_predictor_crs_saved)) {
        stop(
          "The saved habitat model metadata is missing 'landcover_predictor_crs'. Rerun 04_fit_habitat_model.R for ",
          speciesName,
          ".",
          call. = FALSE
        )
      }
      
      landcover_manifest_saved <- load_user_specific_landcover_manifest(landcover_manifest_path_saved)
      habitat_stack_active <- materialize_user_specific_landcover_stack(
        stack_rows = landcover_manifest_saved$current_rows,
        period = "current",
        scenario = "current"
      )
      
      saved_landcover_crs <- sf::st_crs(landcover_predictor_crs_saved)
      active_landcover_crs <- sf::st_crs(terra::crs(habitat_stack_active))
      if (!isTRUE(saved_landcover_crs == active_landcover_crs)) {
        stop(
          "The current habitat stack loaded from the saved manifest does not match the CRS stored in the habitat model metadata.",
          call. = FALSE
        )
      }
    }
  }
  
  if (climate_input_mode_saved == "user_specific") {
    if (is.null(climate_manifest_path_saved) || !nzchar(climate_manifest_path_saved)) {
      stop(
        "The saved climate model metadata is missing 'climate_manifest_path'. 
        Rerun 03_fit_climate_model.R for ", speciesName,".", call. = FALSE
      )
    }
    
    if (is.null(predictor_crs_saved) || !nzchar(predictor_crs_saved)) {
      stop(
        "The saved climate model metadata is missing 'predictor_crs'. 
        Rerun 03_fit_climate_model.R for ", speciesName, ".", call. = FALSE
      )
    }
    
    climate_manifest_saved <- load_user_specific_climate_manifest(climate_manifest_path_saved)
    climate_stack_active <- materialize_user_specific_current_stack(climate_manifest_saved)
    
    saved_predictor_crs <- sf::st_crs(predictor_crs_saved)
    active_predictor_crs <- sf::st_crs(terra::crs(climate_stack_active))
    
    if (!isTRUE(saved_predictor_crs == active_predictor_crs)) {
      stop(
        "The current climate stack loaded from the saved manifest does not match 
        the CRS stored in the climate model metadata.",
        call. = FALSE
      )
    }
    
    if (landcover_input_mode_saved == "user_specific") {
      if (!isTRUE(active_predictor_crs == sf::st_crs(terra::crs(habitat_stack_active)))) {
        stop(
          "The saved user-specific climate and land-cover predictors do not share the same CRS for ",
          speciesName,
          ".",
          call. = FALSE
        )
      }
    }
    
    euboundary_active <- load_eu_boundary(
      custom_path = custom_eu_boundary_path,
      reference = climate_stack_active[[1]]
    )
    euboundary_active_vect <- terra::vect(euboundary_active)
    
    eu_climate_stack_active <- climate_stack_active %>%
      terra::crop(euboundary_active_vect) %>%
      terra::mask(euboundary_active_vect)
    
    habitat_stack_on_climate <- align_continuous_raster(habitat_stack_active, 
                                                        eu_climate_stack_active[[1]])
    
    eu_sampling_mask_active <- terra::mask(habitat_stack_on_climate[[1]], 
                                           eu_climate_stack_active[[1]])
    
    eu_sampling_cells_active <- terra::global(!is.na(eu_sampling_mask_active), 
                                              "sum", na.rm = TRUE)[[1]]
    
    if (is.na(eu_sampling_cells_active) || eu_sampling_cells_active == 0) {
      
      stop("No overlapping non-NA European cells remain after aligning habitat 
           coverage to the user-specific climate template.",
        call. = FALSE
      )
    }
    
    set.seed(728)
    eu_subsample_active <- terra::spatSample(
      eu_sampling_mask_active,
      size = boyce_background_size,
      method = "random",
      na.rm = TRUE,
      as.points = TRUE
    )
    
    eu_climate_sub_active <- terra::extract(eu_climate_stack_active, 
                                            eu_subsample_active, 
                                            ID = FALSE, xy = FALSE)
    eu_habitat_sub_active <- terra::extract(habitat_stack_on_climate, 
                                            eu_subsample_active, 
                                            ID = FALSE, xy = FALSE)
    
  } else {
    
    climate_manifest_saved <- NULL
    
    if (is.null(default_climate_stack)) {
      default_climate_stack <- terra::rast(climate_path)
    }
    
    climate_stack_active <- default_climate_stack
    
    if (landcover_input_mode_saved == "user_specific") {
      euboundary_active <- load_eu_boundary(
        custom_path = custom_eu_boundary_path,
        reference = habitat_stack_active[[1]]
      )
      euboundary_active_vect <- terra::vect(euboundary_active)
      
      if (is.null(custom_eu_boundary_path)) {
        eu_sampling_mask_active <- habitat_stack_active[[1]]
        eu_climate_stack_active <- terra::rast(eu_climpreds_path) %>%
          terra::project(habitat_stack_active[[1]])
      } else {
        eu_sampling_mask_active <- habitat_stack_active[[1]] %>%
          terra::crop(euboundary_active_vect) %>%
          terra::mask(euboundary_active_vect)
        eu_climate_stack_active <- terra::rast(eu_climpreds_path) %>%
          terra::project(habitat_stack_active[[1]]) %>%
          terra::crop(euboundary_active_vect) %>%
          terra::mask(euboundary_active_vect)
      }
      
      eu_sampling_cells_active <- terra::global(!is.na(eu_sampling_mask_active),
                                                "sum", na.rm = TRUE)[[1]]
      if (is.na(eu_sampling_cells_active) || eu_sampling_cells_active == 0) {
        stop("No non-NA European habitat cells are available for validation sampling.", call. = FALSE)
      }
      
      set.seed(728)
      eu_subsample_active <- terra::spatSample(
        eu_sampling_mask_active,
        size = boyce_background_size,
        method = "random",
        na.rm = TRUE,
        as.points = TRUE
      )
      
      eu_climate_sub_active <- terra::extract(eu_climate_stack_active,
                                              eu_subsample_active,
                                              ID = FALSE, xy = FALSE)
      eu_habitat_sub_active <- terra::extract(habitat_stack_active,
                                              eu_subsample_active,
                                              ID = FALSE, xy = FALSE)
    } else {
      if (is.null(default_eu_climate_stack)) {
        if (is.null(custom_eu_boundary_path)) {
          default_eu_climate_stack <- terra::rast(eu_climpreds_path) %>%
            terra::project(habitat_stack[[1]])
        } else {
          default_eu_climate_stack <- terra::rast(eu_climpreds_path) %>%
            terra::project(habitat_stack[[1]]) %>%
            terra::crop(default_euboundary_vect) %>%
            terra::mask(default_euboundary_vect)
        }
        
        set.seed(728)
        default_eu_subsample <- terra::spatSample(
          default_eu_sampling_mask,
          size = boyce_background_size,
          method = "random",
          na.rm = TRUE,
          as.points = TRUE
        )
        
        default_eu_climate_sub <- terra::extract(default_eu_climate_stack, 
                                                 default_eu_subsample, 
                                                 ID = FALSE, xy = FALSE)
        default_eu_habitat_sub <- terra::extract(habitat_stack, 
                                                 default_eu_subsample, 
                                                 ID = FALSE, xy = FALSE)
      }
      
      euboundary_active <- default_euboundary
      eu_climate_stack_active <- default_eu_climate_stack
      eu_climate_sub_active <- default_eu_climate_sub
      eu_habitat_sub_active <- default_eu_habitat_sub
    }
  }
  
  missing_climate_predictors <- setdiff(climate_predictors, 
                                        names(climate_stack_active))
  
  if (length(missing_climate_predictors) > 0) {
    stop("The current climate stack is missing predictor(s) stored in ",
      basename(climate_qs_file), ": ",paste(missing_climate_predictors, 
                                            collapse = ", "),
      call. = FALSE
    )
  }
  
  euboundary_occurrence <- euboundary_active %>%
    sf::st_transform(crs = sf::st_crs(global_presabs))
  rm(climatemodel)
  
  
  #--------------------------------------------------
  #------------------Define validation types---------
  #--------------------------------------------------
  eu_occ <- global_presabs%>%
    dplyr::filter(species==1)%>%
    sf::st_filter(euboundary_occurrence) 
  
  #Only validate climate model in Europe if 40 or more occs 
  eu_climate_validation <- nrow(eu_occ) >= 40
  
  #---------------------------------------------------------
  #- Select climate rasters used in 03_fit_climate_model.R -
  #---------------------------------------------------------
  climate_selection <- terra::subset(climate_stack_active, climate_predictors)
  
  
  #-----------------------------------------------------------
  #- Obtain global climate subsample values for selected predictors
  #-----------------------------------------------------------
  #Load biomes
  wwf_eco_biome <- sf::st_read(biome_path, quiet = TRUE) %>%
    sf::st_transform(crs = get_reference_crs(climate_selection))
  
  # Keep only biome polygons that intersect at least one occurrence point
  global_presences <- dplyr::filter(global_presabs, species==1) %>%
    sf::st_transform(crs = sf::st_crs(wwf_eco_biome))
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
    eu_points <- eu_climate_sub_active %>%
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
    if (is.null(habitatmodel)) {
      stop("Habitat model object was not loaded for ensemble validation.", call. = FALSE)
    }
    eu_presabs <- habitatmodel$eu_presabs
    n_pres_ensemble <- sum(eu_presabs$species == 1)
    
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
      all_presabs$cell <- terra::cellFromXY( climate_stack_active[[1]], 
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
      warning(nrow(global_presabs)- nrow(global_presabs_perfold),
              " global point(s) not assigned to a fold and removed from dataset.")
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
          sf::st_filter(euboundary_occurrence)
        
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
      if (eu_climate_validation || ensemble_validation) {
        datasets$eu_points  <- eu_points
        }
      
      
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
        type = "Climate",
        region = "Global",
        fold = fold,
        all_suit_vals = climate_fav$global_points$median_favourability,
        occ_suit_vals = climate_fav$occ_env$median_favourability,
        abs_suit_vals = climate_fav$abs_env$median_favourability)
      
      #EU
      if(eu_climate_validation){
        eu_validation_climate[[fold]] <- compute_validation_metrics(
          species= speciesName,
          type = "Climate",
          region = "Europe",
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
          sf::st_filter(euboundary_occurrence)
        
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
        type =  "Climate",
        region = "Global",
        fold = "No cross-validation",
        all_suit_vals = median_fav_climate$global_points$median_favourability,
        occ_suit_vals = median_fav_climate$occ_env$median_favourability,
        abs_suit_vals = median_fav_climate$abs_env$median_favourability)
      
      
      #EU
      if(eu_climate_validation){
        eu_validation_climate <- compute_validation_metrics(
          species = speciesName,
          type =  "Climate",
          region = "Europe",
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
  climate_validation_dir<-file.path(base_dir, "Climate", "Current", 
                                    "Diagnostics", "Model_validation")
  if(!dir.exists(climate_validation_dir)) dir.create(climate_validation_dir, 
                                                     recursive = TRUE, 
                                                     showWarnings = FALSE)
  
  
  # Export validation summary (mean across folds) when relevant
  if(use_cv || use_cv_climate_only){
    
    #Export per fold validation
    readr::write_csv(global_validation_climate,
                     file.path(climate_validation_dir, 
                               paste0(speciesName, "_global_climate_validation_per_fold.csv"))) 
    
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
  #--- Load data stored in habitat model qs ---
  #--------------------------------------------
  eu_presabs <- habitatmodel$eu_presabs
  top5_habitat_methods  <- habitatmodel$top5_models
  habitat_predictors <- habitatmodel$selected_predictors
  #habitat_predictors<-unique(habitatmodel[["varimp_df"]][["Predictor"]])
  rm(habitatmodel)
  
  
  #---------------------------------------------------------
  #- Select landcover rasters used in 04_fit_climate_model.R -
  #---------------------------------------------------------
  habitat_selection <- terra::subset(habitat_stack_active,
                                     habitat_predictors[habitat_predictors %in%
                                                          names(habitat_stack_active)])
  minimum_habitat_methods <- 3L
  
  
  #-----------------------------------------------------------
  #- Obtain habitat subsample values for selected predictors
  #-----------------------------------------------------------
  eu_habitat_points<-eu_habitat_sub_active %>%
    dplyr::select(any_of(habitat_predictors))%>%
    dplyr::mutate(ID = dplyr::row_number())
  
  
  #-----------------------------------------------------------------
  #            OPTION 1: SPATIAL CROSS VALIDATION
  #-----------------------------------------------------------------
  valid_habitat_folds <- integer(0)
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
      habitat_prediction <- compute_median_favourability_safe(
        habitat_model,
        habitat_datasets,
        top5_habitat_methods,
        prev_ratio,
        min_successful_methods = minimum_habitat_methods
      )
      
      if (length(habitat_prediction$failed_methods) > 0) {
        warning(
          "Habitat validation for species '",
          speciesName,
          "', fold ",
          fold,
          " dropped method(s): ",
          paste(habitat_prediction$failed_methods, collapse = ", "),
          ". Continuing with: ",
          paste(habitat_prediction$successful_methods, collapse = ", "),
          ".",
          call. = FALSE
        )
      }
      
      if (!isTRUE(habitat_prediction$valid)) {
        warning(
          "Skipping habitat validation for species '",
          speciesName,
          "', fold ",
          fold,
          ": only ",
          habitat_prediction$success_count,
          " of ",
          length(top5_habitat_methods),
          " habitat methods predicted successfully.",
          call. = FALSE
        )
        terra::tmpFiles(remove = TRUE)
        next
      }
      
      median_favourability_habitat_perfold[[fold]] <- habitat_prediction$median_favourability
      valid_habitat_folds <- c(valid_habitat_folds, fold)
      
      
      #-----------------------------------------
      #------- Compute Boyce, AUC, and TSS -----
      #-----------------------------------------
      habitat_fav<-median_favourability_habitat_perfold[[fold]]
      
      #EU
      validation_habitat[[fold]] <- compute_validation_metrics(
        species= speciesName,
        type = "Habitat",
        region = "Europe",
        fold = fold,
        all_suit_vals = habitat_fav$eu_habitat_points$median_favourability,
        occ_suit_vals = habitat_fav$occ_hab$median_favourability,
        abs_suit_vals = habitat_fav$abs_hab$median_favourability)
      
      #Clean
      terra::tmpFiles(remove = TRUE)
    }
    
    if (length(valid_habitat_folds) == 0L) {
      warning(
        "Skipping habitat and ensemble validation for species '",
        speciesName,
        "': no habitat cross-validation folds retained at least ",
        minimum_habitat_methods,
        " successful methods.",
        call. = FALSE
      )
      validation_summary[[speciesName]] <- species_validation_summary
      next
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
    habitat_prediction <- compute_median_favourability_safe(
      habitat_model,
      habitat_datasets,
      top5_habitat_methods,
      prev_ratio,
      min_successful_methods = minimum_habitat_methods
    )
    
    if (length(habitat_prediction$failed_methods) > 0) {
      warning(
        "Habitat validation for species '",
        speciesName,
        "' dropped method(s): ",
        paste(habitat_prediction$failed_methods, collapse = ", "),
        ". Continuing with: ",
        paste(habitat_prediction$successful_methods, collapse = ", "),
        ".",
        call. = FALSE
      )
    }
    
    if (!isTRUE(habitat_prediction$valid)) {
      warning(
        "Skipping habitat and ensemble validation for species '",
        speciesName,
        "': only ",
        habitat_prediction$success_count,
        " of ",
        length(top5_habitat_methods),
        " habitat methods predicted successfully.",
        call. = FALSE
      )
      validation_summary[[speciesName]] <- species_validation_summary
      next
    }
    
    median_fav_habitat <- habitat_prediction$median_favourability
    
    
    #-----------------------------------------
    #------- Compute Boyce, AUC, and TSS -----
    #-----------------------------------------
    eu_validation_habitat <- compute_validation_metrics(
      species= speciesName,
      type = "Habitat",
      region = "Europe",
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
    for (fold in valid_habitat_folds) {
      
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
        type = "Ensemble",
        region = "Europe",
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
    hab_fav <- median_fav_habitat
    clim_fav <- median_fav_climate
    
    
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
      type = "Ensemble",
      region= "Europe",
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
