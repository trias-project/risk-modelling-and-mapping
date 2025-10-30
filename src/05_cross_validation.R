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

# Installs all sdm plugin backends if missing (incl. maxent). Safe if already present.
suppressWarnings(try(sdm::installAll(), silent = TRUE))


#--------------------------------------------
#- Load helper functions and configurations -
#--------------------------------------------
source(pr("src","helper_functions.R"))
source(pr("src","Configurations.R"))


#--------------------------------------------
#---------   Load shape of Europe   ---------
#--------------------------------------------
euboundary <- sf::st_read(here("./data/external/GIS/Europe/EUROPE.shp")) 


#--------------------------------------------
#-------- Prepare habitat rasters (once) ----
#--------------------------------------------
habitat_files  <- list.files(pr("data/external/habitat"), pattern = 'tif$', full.names = TRUE)
habitat_rasters <- lapply(habitat_files, terra::rast)

# Compute common intersection extent across all rasters
common_ext <- Reduce(terra::intersect, lapply(habitat_rasters, terra::ext))

# Crop all rasters to the common (smallest) extent
habitat_rasters <- lapply(habitat_rasters, terra::crop, common_ext)

# Combine into raster stack 
habitat_stack <- terra::rast(habitat_rasters)
rm(habitat_rasters)

# Scale + mask
habitat_stack <- terra::scale(habitat_stack, center = TRUE, scale = TRUE)
ref <- habitat_stack[[1]]
habitat_stack <- terra::mask(habitat_stack, ref)

# Persist single multi-layer for workers
hab_cache_dir <- pr("data/external/habitat/_cache")
dir.create(hab_cache_dir, recursive = TRUE, showWarnings = FALSE)
habitat_path <- file.path(hab_cache_dir, "EU_habitat_scaled_masked.tif")
terra::writeRaster(habitat_stack, filename = habitat_path, overwrite = TRUE)

rm(habitat_stack, ref); gc()

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
accepted_taxonkeys <- unique(taxa_info$acceptedTaxonKey)
# accepted_taxonkeys <- accepted_taxonkeys[c(1:10)]
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
plan(multisession, workers = 10)
handlers(global = TRUE)

with_progress({
  p <- progressr::progressor(along = accepted_taxonkeys)
  
  #------------------------------------------
  #----------- Worker function --------------
  #------------------------------------------
  run_one_key <- function(key) {
    # Ensure required pkgs in worker
    if (!"sdm" %in% .packages()) suppressPackageStartupMessages(library(sdm))
    if (!"blockCV" %in% .packages()) suppressPackageStartupMessages(library(blockCV))
    if (!"dplyr" %in% .packages()) suppressPackageStartupMessages(library(dplyr))
    if (!"sf" %in% .packages()) suppressPackageStartupMessages(library(sf))
    if (!"terra" %in% .packages()) suppressPackageStartupMessages(library(terra))
    if (!"raster" %in% .packages()) suppressPackageStartupMessages(library(raster))
    if (!"ecospat" %in% .packages()) suppressPackageStartupMessages(library(ecospat))
    
    # Per-worker tempdirs
    td_terra  <- file.path(tempdir(), paste0("terra_",  Sys.getpid()))
    td_raster <- file.path(tempdir(), paste0("raster_", Sys.getpid()))
    dir.create(td_terra,  showWarnings = FALSE, recursive = TRUE)
    dir.create(td_raster, showWarnings = FALSE, recursive = TRUE)
    
    # Force on-disk processing
    terra::terraOptions(tempdir = td_terra, memfrac = 0.6, todisk = TRUE)
    raster::rasterOptions(tmpdir = td_raster, chunksize = 1e7, maxmemory = 1e8)
    
    # Disable s2 + cleanup on exit
    old_s2 <- sf::sf_use_s2(FALSE)
    on.exit({
      try(sf::sf_use_s2(old_s2), silent = TRUE)
      try(unlink(td_terra,  recursive = TRUE, force = TRUE), silent = TRUE)
      try(unlink(td_raster, recursive = TRUE, force = TRUE), silent = TRUE)
      gc()
    }, add = TRUE)
    
    #------------- Species metadata -------------
    species <- unname(taxon_name_by_key[as.character(key)])
    if (is.na(species) || is.null(species) || !nzchar(species)) {
      warning(sprintf("No species name found for key %s; skipping.", key))
      return(list(taxonkey = key, skipped = TRUE, reason = "no_species_name"))
    }
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
    species_title   <- gsub("_", " ", first_two_words)
    taxonkey        <- key
    
    #------------- Locate EU model + rasters --------
    species_folder   <- pr("data", "projects", project, paste0(first_two_words, "_", taxonkey))
    EU_model_file_qs <- pr("data", "projects", project, paste0(first_two_words, "_", taxonkey),
                           paste0("EU_model_", first_two_words, "_", taxonkey, ".qs"))
    if (!file.exists(EU_model_file_qs)) {
      warning(sprintf("Skipping %s (%s): no EU model file.", species, taxonkey))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "no_EU_qs"))
    }
    biasgrid_file <- pr("data", "projects", project, paste0(first_two_words, "_", taxonkey),
                        "Rasters", "Interim", paste0("Biasgrid_", first_two_words, "_", taxonkey, ".tif"))
    ensemble_file <- pr("data", "projects", project, paste0(first_two_words, "_", taxonkey),
                        "Rasters", "Global",  paste0("Ensemble_median_", first_two_words, "_", taxonkey, ".tif"))
    if (!file.exists(biasgrid_file) || !file.exists(ensemble_file)) {
      warning(sprintf("Skipping %s: missing biasgrid or ensemble raster.", species))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "missing_bias_or_ensemble"))
    }
    global_model_file_qs <- pr("data", "projects", project, paste0(first_two_words, "_", taxonkey),
                               paste0("Global_model_", first_two_words, "_", taxonkey, ".qs"))
    if (!file.exists(global_model_file_qs)) {
      warning(sprintf("Skipping %s (%s): no Global model file.", species, taxonkey))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "no_global_qs"))
    }
    
    #------------- Reopen static data -------------
    euboundary    <- sf::st_read(euboundary_path, quiet = TRUE)
    habitat_stack <- terra::rast(habitat_path)
    wwf_eco_biome <- sf::st_read(wwf_path, quiet = TRUE)
    
    # ====== CRS SINGLE SOURCE OF TRUTH ======
    target_crs <- sf::st_crs(terra::crs(habitat_stack))  # WKT-aware
    
    # Transform/validate vectors to target CRS
    euboundary <- euboundary |>
      sf::st_make_valid() |>
      sf::st_transform(target_crs)
    
    wwf_eco_biome <- wwf_eco_biome |>
      sf::st_make_valid() |>
      sf::st_transform(target_crs)
    
    # Load species-level objects (then free big list)
    EUmodels <- qs::qread(global_model_file_qs)
    methods  <- EUmodels$top5_models
    EU.occ.sf <- EUmodels$occurrences1km |>
      sf::st_as_sf(coords = c("decimalLongitude","decimalLatitude"), crs = 4326) |>
      sf::st_transform(target_crs)
    rm(EUmodels); gc()
    
    biasgrid_sub <- terra::rast(biasgrid_file)
    Global_climate_for_eu <- terra::rast(ensemble_file) |> terra::project(habitat_stack)
    
    #------------- Fast point-in-polygon -------------
    EU.occ.sf <- sf::st_crop(EU.occ.sf, sf::st_bbox(euboundary))
    inside <- lengths(sf::st_intersects(EU.occ.sf, euboundary)) > 0
    eu_occ <- EU.occ.sf[inside, , drop = FALSE]
    
    if (nrow(eu_occ) == 0) {
      warning(sprintf("0 EU occurrences for %s — skipping EU model.", species))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "no_eu_occ"))
    }
    
    # Unique cells + drop NA habitat
    coords_mat <- sf::st_coordinates(eu_occ)
    cells      <- terra::cellFromXY(habitat_stack[[1]], coords_mat)
    keep       <- !duplicated(cells)
    eu_occ     <- eu_occ[keep, , drop = FALSE]
    
    vals_all <- terra::extract(habitat_stack, terra::vect(eu_occ), ID = FALSE)
    non_na   <- stats::complete.cases(vals_all)
    eu_occ   <- eu_occ[non_na, , drop = FALSE]
    
    #-----------------------------------------------
    #------ Limit to 10,000 occupied grid cells ----
    #-----------------------------------------------
    if (nrow(eu_occ) > 10000) {
      if (occurrence_thinning_method == "random") {
        set.seed(101)
        eu_occ <- eu_occ[sample(nrow(eu_occ), 10000, replace = FALSE), ]
      } else if (occurrence_thinning_method == "kmeans_clustering") {
        habitat_data <- terra::extract(habitat_stack, eu_occ, ID = FALSE)
        set.seed(101)
        clust <- kmeans(habitat_data, centers = n_clusters, iter.max = 10, nstart = 1)$cluster
        occ_habitat <- cbind(eu_occ, habitat_data, clust)
        max_per_cluster <- 10000/n_clusters
        row_sample <- sapply(1:10000, function(x) {
          rowids <- which(clust == x)
          sample(rowids, min(max_per_cluster, length(rowids)), replace = FALSE)
        })
        eu_occ <- occ_habitat[row_sample,] %>%
          dplyr::select(decimalLongitude, decimalLatitude, geometry, species)
      }
    }
    
    euocc_xy <- sf::st_coordinates(eu_occ) |> as.data.frame()
    if (nrow(euocc_xy) < 20) {
      warning(sprintf("%d EU occurrences for %s — skipping EU model.", nrow(euocc_xy), species))
      return(list(species = species, taxonkey = taxonkey, skipped = TRUE, reason = "too_few_eu_occ"))
    }
    
    #------------- Align bias grid to habitat -----------
    biasgrid_aligned <- terra::project(biasgrid_sub, habitat_stack[[1]], method = "bilinear")
    biasgrid_aligned <- terra::mask(biasgrid_aligned, habitat_stack[[1]])
    
    #------------- Invaded ecoregions mask --------------
    stopifnot(sf::st_crs(wwf_eco_biome) == sf::st_crs(eu_occ))
    hit_mat <- sf::st_intersects(wwf_eco_biome, eu_occ, sparse = FALSE)
    polys_with_pts <- rowSums(hit_mat) > 0
    wwf_eco_biome_filtered <- wwf_eco_biome[polys_with_pts, , drop = FALSE]
    
    inside_mask <- terra::rasterize(
      terra::vect(wwf_eco_biome_filtered),
      biasgrid_aligned, field = 1, background = NA
    )
    biasgrid_temp <- terra::ifel(!is.na(inside_mask), biasgrid_aligned, 1)
    biasgrid_eu   <- terra::mask(biasgrid_temp, biasgrid_aligned)
    
    #------------- Pseudoabsences ------------------------
    EU_points <- terra::spatSample(
      biasgrid_eu, size = 10000, method = "weights", as.points = TRUE, na.rm = TRUE
    )
    
    # presences (sf) — add explicit XY columns and align column order
    eu_coords <- sf::st_coordinates(eu_occ)
    eu_occ_pa <- eu_occ |>
      dplyr::mutate(
        decimalLongitude = eu_coords[, 1],
        decimalLatitude  = eu_coords[, 2],
        species = "present"
      ) |>
      dplyr::select(decimalLongitude, decimalLatitude, species, geometry)
    
    # pseudoabsences (sf)
    EU_points_sf <- EU_points[, 0] |>
      sf::st_as_sf() |>
      dplyr::mutate(
        coords = sf::st_coordinates(geometry),
        decimalLongitude = coords[, 1],
        decimalLatitude  = coords[, 2],
        species = "absent"
      ) |>
      dplyr::select(decimalLongitude, decimalLatitude, species, geometry)
    
    # combine
    eu_presabs <- rbind(eu_occ_pa, EU_points_sf)
    
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
        
        boyce_test  <- compute_boyce_robust(fit_test,  obs_test)
        boyce_train <- compute_boyce_robust(fit_train, obs_train)
        
        boyce_rows[[j]] <- data.frame(
          fold          = j,
          n_pres_test   = length(obs_test),
          n_pres_train  = length(obs_train),
          n_fit_test    = length(fit_test),
          n_fit_train   = length(fit_train),
          uniq_fit_test = length(unique(fit_test)),
          uniq_fit_trn  = length(unique(fit_train)),
          boyce_test    = boyce_test,
          boyce_train   = boyce_train
        )
        rm(rj, r_test_area, r_train_area, fit_test, fit_train, obs_test, obs_train); gc()
      }
      
      boyce_df <- do.call(rbind, boyce_rows)
      
    } else {
      # ====== Train-only Boyce (no folds) ======
      pts_sf <- sf::st_as_sf(eu_presabs_sp) %>% dplyr::filter(species == 1L)
      rj <- ens_fav_median[[1]]
      fit_all   <- sample_fit_vals(rj, n = 50000, thresh = 1e6)
      obs_train <- safe_extract_vals(rj, pts_sf)
      boyce_train <- compute_boyce_robust(fit_all, obs_train)
      
      boyce_df <- data.frame(
        fold          = 0,
        n_pres_test   = 0L,
        n_pres_train  = length(obs_train),
        n_fit_test    = 0L,
        n_fit_train   = length(fit_all),
        uniq_fit_test = 0L,
        uniq_fit_trn  = length(unique(fit_all)),
        boyce_test    = NA_real_,
        boyce_train   = boyce_train
      )
      rm(rj, fit_all, obs_train); gc()
    }
    
    summarize_boyce <- function(boyce_df, species, key, k_folds) {
      data.frame(
        species          = species,
        gbif_key         = key,
        k_folds          = as.integer(k_folds),
        boyce_test_mean  = mean(boyce_df$boyce_test,  na.rm = TRUE),
        boyce_test_sd    = sd(  boyce_df$boyce_test,  na.rm = TRUE),
        boyce_train_mean = mean(boyce_df$boyce_train, na.rm = TRUE),
        boyce_train_sd   = sd(  boyce_df$boyce_train, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }
    
    boyce_summary <- summarize_boyce(boyce_df,
                                     species = species,
                                     key     = key,
                                     k_folds = k)
    
    # ===== Persist per-species results =====
    out_dir <- pr("data", "projects", project, "Results_CV_Boyce")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)   # ensure exists in worker
    message("Writing outputs to: ", out_dir)
    
    summ_file_tmp <- file.path(out_dir, sprintf("boyce_summary_%s_%s.tmp.csv", first_two_words, taxonkey))
    summ_file     <- file.path(out_dir, sprintf("boyce_summary_%s_%s.csv",      first_two_words, taxonkey))
    utils::write.csv(boyce_summary, summ_file_tmp, row.names = FALSE, na = "")
    ok1 <- file.rename(summ_file_tmp, summ_file)
    if (!ok1) utils::file.copy(summ_file_tmp, summ_file, overwrite = TRUE)
    
    fold_file_tmp <- file.path(out_dir, sprintf("boyce_folds_%s_%s.tmp.csv", first_two_words, taxonkey))
    fold_file     <- file.path(out_dir, sprintf("boyce_folds_%s_%s.csv",      first_two_words, taxonkey))
    utils::write.csv(boyce_df, fold_file_tmp, row.names = FALSE, na = "")
    ok2 <- file.rename(fold_file_tmp, fold_file)
    if (!ok2) utils::file.copy(fold_file_tmp, fold_file, overwrite = TRUE)
    
    ens_file <- file.path(out_dir, sprintf("ens_fav_median_%s_%s.tif", first_two_words, taxonkey))
    terra::writeRaster(ens_fav_median, ens_file, overwrite = TRUE)
    
    rm(ens_fav_median, boyce_df, habitat_stack,
       biasgrid_sub, biasgrid_aligned, biasgrid_temp, biasgrid_eu, inside_mask,
       wwf_eco_biome, wwf_eco_biome_filtered, occ.full.data.df, fullstack_r, eu_presabs_sp,
       EU_points, EU_points_sf, eu_occ, eu_occ_pa, eu_presabs, sb); gc()
    
    list(
      species    = species,
      taxonkey   = taxonkey,
      status     = "ok",
      k          = k,              # 0 for fallback; 1..5 for CV
      summary_csv= summ_file,
      folds_csv  = fold_file,
      ens_tif    = ens_file
    )
  }  # <-- close run_one_key()
  ## -------- end worker --------
  
  # ------- run all keys in parallel -------
  keys_to_run <- accepted_taxonkeys  # or subset for testing
  
  results <- future_lapply(
    keys_to_run,
    function(k) {
      p(sprintf("Running key %s", k))
      tryCatch(
        run_one_key(k),
        error = function(e) list(taxonkey = k, status = "error", message = conditionMessage(e))
      )
    },
    future.seed = TRUE
  )
})  # <-- closes with_progress

# ------- merge per-species summaries into a single CSV -------
out_dir <- pr("data", "projects", project, "Results_CV_Boyce")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

summ_files <- list.files(out_dir, pattern = "^boyce_summary_.*\\.csv$", full.names = TRUE)

if (length(summ_files)) {
  merged <- do.call(rbind, lapply(summ_files, utils::read.csv, stringsAsFactors = FALSE))
  merged_path <- file.path(out_dir, sprintf("boyce_summary_ALL_%s.csv", project))
  utils::write.csv(merged, merged_path, row.names = FALSE, na = "")
  message("Merged summary written to: ", merged_path)
} else {
  warning("No per-species summary files found to merge.")
}
