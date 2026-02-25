#--------------------------------------------
#- Load helper functions and configurations -
#--------------------------------------------
source(here::here("src", "helper_functions.R"))
source(here::here("src", "00_configurations.R"))


#--------------------------------------------
#-----------    Load packages      ----------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")

packages <- c(
  "viridis","dplyr","here","qs","terra","tidyterra","sf","ggplot2","RColorBrewer",
  "magick","patchwork","grid","randomForest","progressr","raster","dismo","caret",
  "caretEnsemble","kableExtra","gbm","PresenceAbsence","RStoolbox","sdm",
  "future","future.apply","sp","ecospat","blockCV"
)
for (package in packages) {
  if (!package %in% rownames(installed.packages())) install.packages(package)
  library(package, character.only = TRUE)
}


#--------------------------------------------
#----- Stable project root (works in workers)
#--------------------------------------------
PROJECT_ROOT <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
pr <- function(...) normalizePath(file.path(PROJECT_ROOT, ...), winslash = "/", mustWork = FALSE)


#--------------------------------------------
#---------   Load euboundary  ---------
#--------------------------------------------
euboundary <- terra::rast(file.path("data", "external", "habitat", "Agriculture.tif"))
euboundary<-(euboundary*0+1)
euboundary <- terra::as.polygons(euboundary, dissolve = TRUE)  # merge adjacent cells
euboundary <- sf::st_as_sf(euboundary)  # convert to sf
euboundary_path<-file.path("data", "external", "GIS", "Europe", "EUboundary.shp")
sf::write_sf(euboundary, euboundary_path)
stopifnot(file.exists(euboundary_path))


#---------------------------------------------
#------------ Define habitat path  -----------
#---------------------------------------------
processed_folder<-file.path("data", "external", "habitat", "processed")
habitat_path <- file.path(processed_folder, "habitat_stack.tif")


#---------------------------------------------
#--------- Load WWF ecoregions file ----------
#---------------------------------------------
wwf_path <- pr("data/external/GIS/official/wwf_terr_ecos.shp")
stopifnot(file.exists(wwf_path))


#--------------------------------------------
#------------- Load species data ------------
#--------------------------------------------
taxa_info_path <- pr("data/projects", project, paste0(project, "_taxa_info.csv"))
stopifnot(file.exists(taxa_info_path))
taxa_info <- read.csv2(taxa_info_path)

# compact key->name map for workers
taxon_name_by_key <- setNames(
  as.character(taxa_info$acceptedScientificName),
  as.character(taxa_info$acceptedTaxonKey)
)

#Select unique taxonkeys
accepted_taxonkeys <- unique(taxa_info$acceptedTaxonKey)
rm(taxa_info); gc()


#--------------------------------------------
#------ Prepare results directory upfront ----
#--------------------------------------------
results_dir <- pr("data","projects", project, "Results_CV_Boyce")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)


#--------------------------------------------
#------ Parallel plan & progress setup -------
#--------------------------------------------
options(future.globals.maxSize = 2 * 1024^3)
plan(multisession, workers = 2)
handlers(global = TRUE)

with_progress({
  p <- progressr::progressor(along = accepted_taxonkeys)
  
  #------------------------------------------
  #----------- Worker function --------------
  #------------------------------------------
  run_one_key <- function(key) {
    
    #---------------------------------------------------------
    #------ Ensure required packages are loaded in worker ----
    #---------------------------------------------------------
    if (!"sdm" %in% .packages()) suppressPackageStartupMessages(library(sdm))
    if (!"blockCV" %in% .packages()) suppressPackageStartupMessages(library(blockCV))
    if (!"dplyr" %in% .packages()) suppressPackageStartupMessages(library(dplyr))
    if (!"sf" %in% .packages()) suppressPackageStartupMessages(library(sf))
    if (!"terra" %in% .packages()) suppressPackageStartupMessages(library(terra))
    if (!"raster" %in% .packages()) suppressPackageStartupMessages(library(raster))
    if (!"ecospat" %in% .packages()) suppressPackageStartupMessages(library(ecospat))
    
    
    #--------------------------------------------
    #---------Create perworker tempdirs ---------
    #--------------------------------------------
    td_terra  <- file.path(tempdir(), paste0("terra_",  Sys.getpid()))
    td_raster <- file.path(tempdir(), paste0("raster_", Sys.getpid()))
    dir.create(td_terra,  showWarnings = FALSE, recursive = TRUE)
    dir.create(td_raster, showWarnings = FALSE, recursive = TRUE)
    
    # Force on-disk processing
    terra::terraOptions(tempdir = td_terra, memfrac = 0.6, todisk = TRUE)
    raster::rasterOptions(tmpdir = td_raster, chunksize = 1e7, maxmemory = 1e8)
    
    
    #--------------------------------------------
    #------ Disable s2 and cleanup on exit ----
    #--------------------------------------------
    old_s2 <- sf::sf_use_s2(FALSE)
    on.exit({
      try(sf::sf_use_s2(old_s2), silent = TRUE)
      try(unlink(td_terra,  recursive = TRUE, force = TRUE), silent = TRUE)
      try(unlink(td_raster, recursive = TRUE, force = TRUE), silent = TRUE)
      gc()
    }, add = TRUE)
    
    
    #--------------------------------------------
    #-------------- Load species data -----------
    #--------------------------------------------
    species <- unname(taxon_name_by_key[as.character(key)])
    if (is.na(species) || is.null(species) || !nzchar(species)) {
      warning(sprintf("No species name found for key %s; skipping.", key))
      return(list(taxonkey = key, skipped = TRUE, reason = "no_species_name"))
    }
    
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
    species_title   <- gsub("_", " ", first_two_words)
    taxonkey        <- key
    
    
    #--------------------------------------------
    #--------- Define necessary paths------------
    #--------------------------------------------
    species_folder   <- file.path("data", "projects", project, paste0(first_two_words, "_", taxonkey))
    global_model_file_qs <- file.path(species_folder,"Climate",
                                      paste0("Climate_model_", first_two_words, "_", taxonkey, ".qs"))
    EU_model_file_qs <- file.path(species_folder,"Habitat", 
                           paste0("Habitat_model_", first_two_words, "_", taxonkey, ".qs"))
    biasgrid_file <- file.path(species_folder,"Climate","Current", "Interim", 
                        paste0("Biasgrid_", first_two_words, "_", taxonkey, ".tif"))
    ensemble_file <- file.path(species_folder,"Combined", "Current","Predictions","rasters",
                        paste0(first_two_words, "_Combined_current_ensemble",".tif"))
    
    #--------------------------------------------
    #----------- Check if files exist -----------
    #--------------------------------------------
    #Check if climate model qs exists
     if (!file.exists(global_model_file_qs)) {
      warning(sprintf("Skipping %s (%s): no Global model file.", species, taxonkey))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "no_global_qs"))
    }
    
    #Check if a habitat model exists
    if (!file.exists(EU_model_file_qs)) {
      warning(sprintf("Skipping %s (%s): no EU model file.", species, taxonkey))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "no_EU_qs"))
    }
    
    #Check if biasgrid and combined model predictions exist
    if (!file.exists(biasgrid_file)) {
      warning(sprintf("Skipping %s: missing biasgrid or ensemble raster.", species))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "missing biasgrid"))
    }
    
    
    #--------------------------------------------
    #-------Load  necessary layers and files-----
    #--------------------------------------------
    # ====== CRS SINGLE SOURCE OF TRUTH ======
    target_crs <- sf::st_crs(terra::crs(habitat_stack)) 
    
    euboundary    <- sf::st_read(euboundary_path, quiet = TRUE)
    habitat_stack <- terra::rast(habitat_path)
    wwf_eco_biome <- sf::st_read(wwf_path, quiet = TRUE)%>%
      sf::st_make_valid() %>%
      sf::st_transform(target_crs)
    
    
    #--------------------------------------------
    #----Load  data stored in climate model qs---
    #--------------------------------------------
    climatemodel <- qs::qread(global_model_file_qs)
    methods  <- climatemodel$top5_models
    global_presabs <- climatemodel$global_presabs
    rm(climatemodel)
    gc()
    
    biasgrid_sub <- terra::rast(biasgrid_file)
    Global_climate_for_eu <- terra::rast(ensemble_file) 

    
    #------------- Correlation filter --------------------
    presabs_df <- terra::extract(habitat_stack, terra::vect(eu_presabs), ID = FALSE)
    cor_matrix <- stats::cor(presabs_df, use = "complete.obs")
    drop_vars  <- caret::findCorrelation(cor_matrix, cutoff = 0.7, exact = TRUE, names = TRUE)
    fullstack  <- subset(habitat_stack, !(names(habitat_stack) %in% drop_vars))
    
    occ.full.data.df <- terra::extract(fullstack, terra::vect(eu_presabs), ID = FALSE) |>
      dplyr::mutate(occ = eu_presabs$species)
    if (anyNA(occ.full.data.df)) warning("Some EU points fell in NA habitat cells")
    
    # sdm data (Raster*/sp), numeric species
    eu_presabs_num <- eu_presabs |>
      dplyr::mutate(species = ifelse(species == "present", 1, 0)) |>
      dplyr::select(-decimalLongitude, -decimalLatitude)
    eu_presabs_sp <- methods::as(eu_presabs_num, "Spatial")
    fullstack_r   <- raster::stack(fullstack)
    
    # Ensure Raster* advertises SAME WKT as sf
    raster::crs(fullstack_r) <- sf::st_crs(terra::crs(habitat_stack))$wkt
    
    rm(presabs_df, cor_matrix, drop_vars, fullstack); gc()
    
    # ==== CONFIG ====
    set.seed(42)
    hex_size_km <- 100
    hex_size_m  <- hex_size_km * 1000
    
    # ==== HELPERS ====
    count_pres <- function(spdf, col = "species") sum(spdf[[col]] == 1L, na.rm = TRUE)
    count_abs  <- function(spdf, col = "species") sum(spdf[[col]] == 0L, na.rm = TRUE)
    
    favourability_from_prob <- function(prob, prev_ratio) {
      f <- function(p) {
        odds <- p / (1 - p)
        fav  <- odds / (prev_ratio + odds)
        fav[!is.finite(fav)] <- NA
        fav[fav < 0] <- 0
        fav[fav > 1] <- 1
        fav
      }
      if (inherits(prob, "SpatRaster")) {
        terra::app(prob, f)
      } else if (inherits(prob, "Raster")) {
        raster::calc(prob, f)
      } else {
        stop("Unsupported raster class: ", class(prob)[1])
      }
    }
    
    .build_stack <- function(x_list, k) {
      keep <- !vapply(x_list, is.null, logical(1))
      x_list <- x_list[keep]
      if (!length(x_list)) return(NULL)
      cls <- unique(vapply(x_list, function(x) class(x)[1], character(1)))
      if (all(cls == "SpatRaster")) {
        st <- terra::rast(x_list)
        names(st) <- paste0("fold_", which(keep))
        outfile <- file.path(td_terra, paste0("stack_", as.integer(runif(1,1,1e9)), ".tif"))
        st <- terra::writeRaster(st, outfile, overwrite = TRUE)
        return(st)
      } else {
        st <- raster::stack(x_list)
        names(st) <- paste0("fold_", which(keep))
        outfile <- raster::rasterTmpFile()
        st <- raster::writeRaster(st, filename = outfile, overwrite = TRUE)
        return(st)
      }
    }
    
    # ==== MAIN FLOW ====
    n_pres <- count_pres(eu_presabs_sp, "species")
    k <- 0L
    use_cv <- FALSE
    
    if (n_pres >= 40L) {
      k <- min(5L, floor(n_pres / 20L))
      use_cv <- k >= 2L
    }
    
    ens_fav_median <- NULL
    sb <- NULL
    
    if (use_cv) {
      # Final CRS sanity before blockCV — enforce same CRS as the raster
      eu_presabs_sp <- methods::as(
        sf::st_transform(sf::st_as_sf(eu_presabs_sp), sf::st_crs(raster::crs(fullstack_r))),
        "Spatial"
      )
      
      # Hex, class-balanced spatial folds
      sb <- blockCV::cv_spatial(
        x         = eu_presabs_sp,
        column    = "species",
        r         = fullstack_r,
        k         = k,
        hexagon   = TRUE,
        selection = "random",
        iteration = 200,
        size      = hex_size_m
      )
      fold_ids <- sb$folds_ids
      stopifnot(length(fold_ids) == nrow(eu_presabs_sp))
      
      # Per-fold train/predict (keep favourability only)
      fold_fav  <- vector("list", k)
      
      for (i in seq_len(k)) {
        message(sprintf("Fold %d/%d: training on folds != %d", i, k, i))
        train_idx <- which(fold_ids != i)
        eu_train  <- eu_presabs_sp[train_idx, ]
        
        sdm_data <- sdm::sdmData(
          species ~ .,
          train      = eu_train,
          predictors = fullstack_r
        )
        model <- sdm::sdm(species ~ ., data = sdm_data, methods = methods)
        
        # Prevalence ratio from TRAINING data
        pres_tr <- count_pres(eu_train, "species")
        abs_tr  <- count_abs(eu_train,  "species")
        prev_ratio <- abs_tr / max(1L, pres_tr)
        
        # Build lookup: method name -> numeric modelID
        mi <- sdm::getModelInfo(model)
        col_m <- if ("methods" %in% names(mi)) "methods" else if ("method" %in% names(mi)) "method" else stop("getModelInfo: no method column")
        stopifnot("modelID" %in% names(mi))
        mi$modelID <- as.integer(mi$modelID)
        id_by_method <- split(mi[["modelID"]], mi[[col_m]])
        
        fav_i <- vector("list", length(methods)); names(fav_i) <- methods
        
        cap01 <- function(x) { x[x < 0] <- 0; x[x > 1] <- 1; x }
        
        for (m in methods) {
          id_m <- as.integer(id_by_method[[m]][1])
          if (length(id_m) == 0L || is.na(id_m)) {
            warning("No numeric model id found for method ", m, "; skipping.")
            next
          }
          message("  Predicting with ", m, " (id=", id_m, ") ...")
          
          pr <- try(predict(model, newdata = fullstack_r, id = id_m), silent = TRUE)
          if (inherits(pr, "try-error")) {
            message("    predict() failed for ", m, " (id=", id_m, "): ",
                    conditionMessage(attr(pr, "condition")))
            next
          }
          if (inherits(pr, "SpatRaster")) pr <- raster::raster(pr)
          pr <- raster::calc(pr, fun = cap01)
          
          fv <- favourability_from_prob(pr, prev_ratio)
          fv_file <- raster::rasterTmpFile()
          fv <- raster::writeRaster(fv, filename = fv_file, overwrite = TRUE)
          
          fav_i[[m]] <- fv
          rm(pr, fv); gc()
        }
        
        fold_fav[[i]] <- fav_i
        rm(model, sdm_data, eu_train); gc()
      }
      
      # Per-method stacks (layers = held-out folds)
      fav_by_method  <- setNames(vector("list", length(methods)), methods)
      for (m in methods) {
        layers_fav_m  <- lapply(seq_len(k), function(i) fold_fav[[i]][[m]])
        fav_by_method[[m]]  <- .build_stack(layers_fav_m,  k)
      }
      rm(fold_fav); gc()
      
      # Median ensemble per fold
      ensemble_median_by_fold <- function(fbm, k) {
        stopifnot(k > 0L)
        out_layers <- vector("list", k)
        for (j in seq_len(k)) {
          layers_j <- lapply(fbm, function(st) {
            if (is.null(st)) return(NULL)
            if (inherits(st, "SpatRaster")) {
              if (terra::nlyr(st) < j) return(NULL)
              st[[j]]
            } else if (inherits(st, c("RasterLayer","RasterBrick","RasterStack"))) {
              if (raster::nlayers(st) < j) return(NULL)
              terra::rast(st[[j]])
            } else {
              NULL
            }
          })
          layers_j <- layers_j[!vapply(layers_j, is.null, logical(1))]
          if (!length(layers_j)) stop(sprintf("No layers found for fold %d", j))
          s <- terra::rast(layers_j)
          tmpstack <- file.path(td_terra, paste0("fold_", j, "_stack.tif"))
          s <- terra::writeRaster(s, tmpstack, overwrite=TRUE)
          tmpmed <- file.path(td_terra, paste0("fold_", j, "_median.tif"))
          out_layers[[j]] <- terra::app(s, median, na.rm=TRUE, filename=tmpmed, overwrite=TRUE)
          rm(s); gc()
        }
        ens_file <- file.path(td_terra, "ens_fav_median.tif")
        ens <- terra::rast(out_layers)
        terra::writeRaster(ens, ens_file, overwrite=TRUE)
        ens
      }
      
      ens_fav_median <- ensemble_median_by_fold(fav_by_method, k)
      
    } else {
      # ===== Fallback: train on ALL data (no folds) =====
      message("Not enough presences for CV (n_pres=", n_pres, "). Running train-only (no folds).")
      
      sdm_data_all <- sdm::sdmData(
        species ~ .,
        train      = eu_presabs_sp,
        predictors = fullstack_r
      )
      model_all <- sdm::sdm(species ~ ., data = sdm_data_all, methods = methods)
      
      pres_all <- count_pres(eu_presabs_sp, "species")
      abs_all  <- count_abs(eu_presabs_sp,  "species")
      prev_ratio_all <- abs_all / max(1L, pres_all)
      
      mi <- sdm::getModelInfo(model_all)
      col_m <- if ("methods" %in% names(mi)) "methods" else if ("method" %in% names(mi)) "method" else stop("getModelInfo: no method column")
      stopifnot("modelID" %in% names(mi))
      mi$modelID <- as.integer(mi$modelID)
      id_by_method <- split(mi[["modelID"]], mi[[col_m]])
      
      cap01 <- function(x) { x[x < 0] <- 0; x[x > 1] <- 1; x }
      fav_list <- list()
      for (m in methods) {
        id_m <- as.integer(id_by_method[[m]][1])
        if (length(id_m) == 0L || is.na(id_m)) {
          warning("No numeric model id found for method ", m, "; skipping.")
          next
        }
        message("  Predicting (train-only) with ", m, " (id=", id_m, ") ...")
        pr <- try(predict(model_all, newdata = fullstack_r, id = id_m), silent = TRUE)
        if (inherits(pr, "try-error")) {
          message("    predict() failed for ", m, " (id=", id_m, "): ",
                  conditionMessage(attr(pr, "condition")))
          next
        }
        if (inherits(pr, "SpatRaster")) pr <- raster::raster(pr)
        pr <- raster::calc(pr, fun = cap01)
        
        fv <- favourability_from_prob(pr, prev_ratio_all)
        fav_list[[m]] <- terra::rast(fv)
        rm(pr, fv); gc()
      }
      fav_list <- fav_list[!vapply(fav_list, is.null, logical(1))]
      stopifnot(length(fav_list) > 0L)
      s <- terra::rast(fav_list)
      ens_fav_median <- terra::app(s, median, na.rm = TRUE)
      names(ens_fav_median) <- "full"
      rm(s, fav_list, model_all, sdm_data_all); gc()
      k <- 0L
    }
    
    #---- merge with climate predictions from global model
    if (!terra::compareGeom(ens_fav_median, Global_climate_for_eu, stopOnError = FALSE)) {
      Global_climate_for_eu <- terra::resample(
        Global_climate_for_eu, ens_fav_median, method = "bilinear"
      )
    }
    ens_fav_median <- sqrt(ens_fav_median * Global_climate_for_eu)
    if (k > 0L && terra::nlyr(ens_fav_median) == k) {
      names(ens_fav_median) <- paste0("fold_", seq_len(k))
    }
    rm(Global_climate_for_eu); gc()
    
    #################################################
    #####        EVALUATION STATISTICS           ####
    #################################################
    get_boyce_cor <- function(res) {
      if (is.null(res)) return(NA_real_)
      if (!is.null(res$Spearman.cor)) return(as.numeric(res$Spearman.cor)[1])
      if (!is.null(res$cor))          return(as.numeric(res$cor)[1])
      NA_real_
    }
    safe_extract_vals <- function(r, sf_points) {
      if (!nrow(sf_points)) return(numeric(0))
      out <- try(terra::extract(r, terra::vect(sf_points)), silent = TRUE)
      if (inherits(out, "try-error") || is.null(out) || ncol(out) < 2) return(numeric(0))
      as.numeric(out[[2]])
    }
    sample_fit_vals <- function(r, n = 50000, thresh = 1e6) {
      n_cells <- try(terra::global(!is.na(r), "sum", na.rm = TRUE)[1,1], silent = TRUE)
      if (!inherits(n_cells, "try-error") && is.finite(n_cells) && n_cells <= thresh) {
        v <- as.vector(terra::values(r))
        return(v[is.finite(v)])
      }
      smp <- try(terra::spatSample(r, size = n, method = "random", na.rm = TRUE,
                                   as.points = FALSE, values = TRUE), silent = TRUE)
      if (inherits(smp, "try-error") || is.null(smp) || ncol(smp) < 1) return(numeric(0))
      as.numeric(smp[[1]])
    }
    compute_boyce_robust <- function(fit_vals, obs_vals) {
      fit_vals <- fit_vals[is.finite(fit_vals)]
      obs_vals <- obs_vals[is.finite(obs_vals)]
      if (length(obs_vals) < 5 || length(fit_vals) < 200) return(NA_real_)
      if (length(unique(fit_vals)) < 3) return(NA_real_)
      res1 <- try(ecospat.boyce(fit = fit_vals, obs = obs_vals, nclass = 0, PEplot = FALSE),
                  silent = TRUE)
      sc <- get_boyce_cor(if (!inherits(res1, "try-error")) res1 else NULL)
      if (is.finite(sc)) return(sc)
      for (nc in c(10, 20)) {
        res2 <- try(ecospat.boyce(fit = fit_vals, obs = obs_vals, nclass = nc, PEplot = FALSE),
                    silent = TRUE)
        sc2 <- get_boyce_cor(if (!inherits(res2, "try-error")) res2 else NULL)
        if (is.finite(sc2)) return(sc2)
      }
      NA_real_
    }
    
    boyce_df <- NULL
    
    if (k > 0L) {
      # ====== CV-based Boyce (train/test per fold) ======
      sb_cv <- sb
      blocks_sf <- sb_cv$blocks
      pts_sf    <- sf::st_as_sf(eu_presabs_sp)
      if (!sf::st_crs(pts_sf) == sf::st_crs(blocks_sf)) {
        pts_sf <- sf::st_transform(pts_sf, sf::st_crs(blocks_sf))
      }
      pts_sf$fold <- sb_cv$folds_ids
      blocks_sf$BLOCK_ROWID <- seq_len(nrow(blocks_sf))
      
      pts_with_block <- sf::st_join(
        pts_sf[c("species","fold")],
        blocks_sf["BLOCK_ROWID"],
        join = sf::st_within,
        left = FALSE
      )
      blocks_for_folds <- function(folds_vec) {
        unique(pts_with_block$BLOCK_ROWID[pts_with_block$fold %in% folds_vec])
      }
      
      k_layers <- terra::nlyr(ens_fav_median)
      stopifnot(k_layers == k, k > 0L)
      
      boyce_rows <- vector("list", k)
      for (j in seq_len(k)) {
        message(sprintf("Boyce fold %d/%d", j, k))
        
        pres_test  <- pts_sf %>% dplyr::filter(species == 1L, fold == j)
        pres_train <- pts_sf %>% dplyr::filter(species == 1L, fold != j)
        
        blk_test_ids  <- blocks_for_folds(j)
        blk_train_ids <- blocks_for_folds(setdiff(seq_len(k), j))
        
        blocks_test   <- blocks_sf %>% dplyr::filter(BLOCK_ROWID %in% blk_test_ids)
        blocks_train  <- blocks_sf %>% dplyr::filter(BLOCK_ROWID %in% blk_train_ids)
        
        rj <- ens_fav_median[[j]]
        
        r_test_area  <- if (nrow(blocks_test))  terra::mask(rj, terra::vect(blocks_test))  else rj * NA
        r_train_area <- if (nrow(blocks_train)) terra::mask(rj, terra::vect(blocks_train)) else rj * NA
        
        fit_test  <- sample_fit_vals(r_test_area,  n = 50000, thresh = 1e6)
        fit_train <- sample_fit_vals(r_train_area, n = 50000, thresh = 1e6)
        
        obs_test  <- safe_extract_vals(rj, pres_test)
        obs_train <- safe_extract_vals(rj, pres_train)
        
        bc_test  <- compute_boyce_robust(fit_test,  obs_test)
        bc_train <- compute_boyce_robust(fit_train, obs_train)
        
        boyce_rows[[j]] <- data.frame(
          fold       = j,
          boyce_test = bc_test,
          boyce_train= bc_train,
          n_test     = length(obs_test),
          n_train    = length(obs_train),
          stringsAsFactors = FALSE
        )
        rm(r_test_area, r_train_area, fit_test, fit_train, obs_test, obs_train, rj); gc()
      }
      boyce_df <- do.call(rbind, boyce_rows)
      
    } else {
      # ===== No CV: evaluate on full training set =====
      r_full <- ens_fav_median[[1]]
      pres_all_sf <- pts_sf %>% dplyr::filter(species == 1L)
      
      obs_all <- safe_extract_vals(r_full, pres_all_sf)
      fit_all <- sample_fit_vals(r_full, n = 50000, thresh = 1e6)
      bc_all  <- compute_boyce_robust(fit_all, obs_all)
      
      boyce_df <- data.frame(
        fold        = 0,
        boyce_test  = NA_real_,
        boyce_train = bc_all,
        n_test      = 0,
        n_train     = length(obs_all),
        stringsAsFactors = FALSE
      )
      rm(r_full, pres_all_sf, obs_all, fit_all); gc()
    }
    
    # ==== RETURN ====
    p()
    return(list(
      species    = species,
      taxonkey   = taxonkey,
      boyce_df   = boyce_df,
      k          = k,
      skipped    = FALSE
    ))
  }
  
  #-------------------------------------------------------
  #----------- Execute worker across species  -----------
  #-------------------------------------------------------
  results_list <- future.apply::future_lapply(
    accepted_taxonkeys,
    run_one_key,
    future.seed = TRUE
  )
  
})

#-------------------------------------------------------
#----------- Collect & save all Boyce results ---------
#-------------------------------------------------------
boyce_summary <- lapply(seq_along(results_list), function(i) {
  res <- results_list[[i]]
  if (isTRUE(res$skipped)) {
    return(data.frame(
      species      = if (!is.null(res$species)) res$species else NA_character_,
      taxonkey     = res$taxonkey,
      fold         = NA_integer_,
      boyce_test   = NA_real_,
      boyce_train  = NA_real_,
      n_test       = NA_integer_,
      n_train      = NA_integer_,
      k            = NA_integer_,
      skipped      = TRUE,
      reason       = if (!is.null(res$reason)) res$reason else "unknown",
      stringsAsFactors = FALSE
    ))
  }
  df_b <- res$boyce_df
  df_b$species <- res$species
  df_b$taxonkey <- res$taxonkey
  df_b$k <- res$k
  df_b$skipped <- FALSE
  df_b$reason <- NA_character_
  return(df_b)
})
boyce_all <- do.call(rbind, boyce_summary)

out_csv <- file.path(results_dir, paste0(project, "_Boyce_CV_results.csv"))
write.csv(boyce_all, out_csv, row.names = FALSE)
message("Saved cross-validated Boyce results to: ", out_csv)

# Summary stats
library(dplyr)
library(readr)

boyce_species_summary <- boyce_all %>%
  group_by(species, taxonkey) %>%
  summarise(
    n_folds_total    = n(),
    boyce_test_mean  = mean(boyce_test, na.rm = TRUE),
    boyce_test_sd    = sd(boyce_test, na.rm = TRUE),
    boyce_train_mean = mean(boyce_train, na.rm = TRUE),
    boyce_train_sd   = sd(boyce_train, na.rm = TRUE),
    boyce_gap_mean   = boyce_train_mean - boyce_test_mean,
    .groups = "drop"
  )

write_csv(boyce_species_summary, "boyce_species_summary.csv")
boyce_species_summary
