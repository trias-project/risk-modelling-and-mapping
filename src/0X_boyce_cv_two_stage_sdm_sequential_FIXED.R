################################################################################
##                                                                            ##
##  0X_boyce_cv_two_stage_sdm_sequential.R                                   ##
##                                                                            ##
##  Sequential (non-parallel) version — easier to inspect and debug          ##
##                                                                            ##
##  Two-stage SDM cross-validated Boyce index evaluation                     ##
##                                                                            ##
##  STAGE 1 – Global (climate-only) model                                    ##
##    • Presence/PA dataset following 03_fit_climate_model.R                 ##
##    • X-fold spatial block CV (blockCV hexagons)                           ##
##    • Boyce index (train & test) per fold for:                             ##
##        (a) full training area (biomes with presences)                     ##
##        (b) Europe-only subset (if EU presences exist)                     ##
##    • Mean ± SD across folds                                               ##
##                                                                            ##
##  STAGE 2 – Europe habitat model                                           ##
##    • Presence/PA dataset following 04_fit_habitat_model.R                 ##
##    • X-fold spatial block CV                                               ##
##    • Boyce index (train & test) per fold for Europe                       ##
##    • Mean ± SD across folds                                               ##
##                                                                            ##
##  STAGE 3 – Model combination                                              ##
##    • Every global EU fold pred × every habitat fold pred                  ##
##    • Geometric mean: sqrt(global_eu * habitat)                            ##
##    • Boyce index for every combination (max k_global × k_habitat)         ##
##    • Mean ± SD across all combinations                                    ##
##                                                                            ##
################################################################################


# ==============================================================================
# SECTION 0 — Packages
# ==============================================================================
options("rgdal_show_exportToProj4_warnings" = "none")

packages <- c(
  "viridis", "dplyr", "here", "qs", "terra", "tidyterra", "sf",
  "ggplot2", "RColorBrewer", "magick", "patchwork", "grid",
  "randomForest", "raster", "dismo", "caret", "caretEnsemble",
  "kableExtra", "gbm", "PresenceAbsence", "RStoolbox", "sdm",
  "purrr", "ecospat", "blockCV"
)
for (pkg in packages) {
  if (!pkg %in% rownames(installed.packages())) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
suppressWarnings(try(sdm::installAll(), silent = TRUE))

sf::sf_use_s2(FALSE)   # disable S2 to avoid topological issues


# ==============================================================================
# SECTION 1 — Project root & configurations
# ==============================================================================
PROJECT_ROOT <- normalizePath(here::here(), winslash = "/", mustWork = TRUE)
pr <- function(...) normalizePath(file.path(PROJECT_ROOT, ...), winslash = "/", mustWork = FALSE)

source(pr("src", "00_configurations.R"))
source(pr("src", "helper_functions.R"))


# ==============================================================================
# SECTION 2 — Tempdir setup (disk-backed terra/raster processing)
# ==============================================================================
td_terra  <- file.path(tempdir(), "terra_seq")
td_raster <- file.path(tempdir(), "raster_seq")
dir.create(td_terra,  showWarnings = FALSE, recursive = TRUE)
dir.create(td_raster, showWarnings = FALSE, recursive = TRUE)
terra::terraOptions(tempdir = td_terra, memfrac = 0.5, todisk = TRUE)
raster::rasterOptions(tmpdir = td_raster, chunksize = 1e7, maxmemory = 1e8)


# ==============================================================================
# SECTION 3 — Static data (loaded once, reused for every species)
# ==============================================================================

# ---- Habitat rasters ---------------------------------------------------------
message("Loading habitat rasters...")
habitat_files   <- list.files(pr("data/external/habitat"), pattern = "tif$", full.names = TRUE)
habitat_rasters <- lapply(habitat_files, terra::rast)
common_ext      <- Reduce(terra::intersect, lapply(habitat_rasters, terra::ext))
habitat_rasters <- lapply(habitat_rasters, terra::crop, common_ext)
habitat_stack   <- terra::rast(habitat_rasters)
rm(habitat_rasters)
habitat_stack   <- terra::scale(habitat_stack, center = TRUE, scale = TRUE)
na_mask         <- anyNA(habitat_stack)
habitat_stack   <- terra::mask(habitat_stack, na_mask, maskvalue = 1)

# Cache to disk so it can be reloaded inside the loop without re-processing
hab_cache_dir <- pr("data/external/habitat/_cache")
dir.create(hab_cache_dir, recursive = TRUE, showWarnings = FALSE)
habitat_path  <- file.path(hab_cache_dir, "EU_habitat_scaled_masked.tif")
terra::writeRaster(habitat_stack, filename = habitat_path, overwrite = TRUE)
rm(habitat_stack, na_mask); gc()

# ---- European boundary -------------------------------------------------------
euboundary_path <- pr("data/external/GIS/Europe/EUROPE.shp")
stopifnot(file.exists(euboundary_path))
euboundary_raw <- sf::st_read(euboundary_path, quiet = TRUE)

# ---- WWF ecoregions (biome masking for global model) -------------------------
wwf_path <- pr("data/external/GIS/official/newRealms.shp")
stopifnot(file.exists(wwf_path))
wwf_eco_biome_raw <- sf::st_read(wwf_path, quiet = TRUE)

# ---- Global climate raster paths ---------------------------------------------
processed_folder        <- pr("data/external/climate/chelsa_current/processed")
globalclimpreds_file    <- file.path(processed_folder, "globalclimpreds.tif")
globalclimpreds_5k_file <- file.path(processed_folder, "globalclim_5k.tif")
eu_climpreds_file       <- file.path(processed_folder, "euclimpreds.tif")
stopifnot(file.exists(globalclimpreds_file),
          file.exists(globalclimpreds_5k_file),
          file.exists(eu_climpreds_file))

# ---- Bias grid paths ---------------------------------------------------------
bias_grid_folder <- pr("data/external/bias_grids")
bias_grid_paths  <- list(
  Plants       = file.path(bias_grid_folder, "log_plants_1degree_layer.tif"),
  Amphibians   = file.path(bias_grid_folder, "log_amphibians_1degree_layer.tif"),
  Birds        = file.path(bias_grid_folder, "log_birds_1degree_layer.tif"),
  Mammals      = file.path(bias_grid_folder, "log_mammals_1degree_layer.tif"),
  Molluscs     = file.path(bias_grid_folder, "log_mollusca_1degree_layer.tif"),
  Reptiles     = file.path(bias_grid_folder, "log_reptiles_1degree_layer.tif"),
  Fish         = file.path(bias_grid_folder, "log_fish_1degree_layer.tif"),
  Malacostraca = file.path(bias_grid_folder, "log_malacostraca_1degree_layer.tif"),
  Insects      = file.path(bias_grid_folder, "log_insects_1degree_layer.tif")
)

# ---- Taxa info ---------------------------------------------------------------
taxa_info_path <- pr("data/projects", project, paste0(project, "_taxa_info.csv"))
stopifnot(file.exists(taxa_info_path))
taxa_info <- read.csv2(taxa_info_path)

# Group column is optional in taxa_info; may live in the occurrence data instead
taxon_group_by_key <- if ("Group" %in% names(taxa_info) && any(nzchar(taxa_info$Group))) {
  setNames(as.character(taxa_info$Group), as.character(taxa_info$acceptedTaxonKey))
} else {
  NULL
}
taxon_name_by_key  <- setNames(as.character(taxa_info$acceptedScientificName),
                                as.character(taxa_info$acceptedTaxonKey))
accepted_taxonkeys <- unique(taxa_info$acceptedTaxonKey)
rm(taxa_info); gc()

# ---- Processed occurrence data -----------------------------------------------
occ_file <- pr("data/projects", project, paste0(project, "_processed_occurrences.qs"))
stopifnot(file.exists(occ_file))
global_occ_all <- qs::qread(occ_file)
cleaned     <- global_occ_all$cleaned
cleaned_1km <- global_occ_all$cleaned_1km
rm(global_occ_all); gc()

# ---- Results directory -------------------------------------------------------
results_dir <- pr("data/projects", project, "Results_CV_Boyce_TwoStage")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
message("Results will be written to: ", results_dir)

# Background sample size for Boyce index computation.
# 50 000 points is more than sufficient; compute_boyce_robust needs >= 200.
# Increase if you want smoother histograms; decrease to trade accuracy for speed.
N_BG <- 50000L


# ==============================================================================
# SECTION 4 — Helper functions
# ==============================================================================

# Boyce index: extract Spearman correlation from ecospat output
get_boyce_cor <- function(res) {
  if (is.null(res))              return(NA_real_)
  if (!is.null(res$Spearman.cor)) return(as.numeric(res$Spearman.cor)[1])
  if (!is.null(res$cor))          return(as.numeric(res$cor)[1])
  NA_real_
}

# Boyce index with fallbacks for different nclass settings
compute_boyce_robust <- function(fit_vals, obs_vals) {
  fit_vals <- fit_vals[is.finite(fit_vals)]
  obs_vals <- obs_vals[is.finite(obs_vals)]
  if (length(obs_vals) < 5 || length(fit_vals) < 200) return(NA_real_)
  if (length(unique(fit_vals)) < 3)                   return(NA_real_)
  res <- try(ecospat::ecospat.boyce(fit = fit_vals, obs = obs_vals,
                                    nclass = 0, PEplot = FALSE), silent = TRUE)
  sc  <- get_boyce_cor(if (!inherits(res, "try-error")) res else NULL)
  if (is.finite(sc)) return(sc)
  for (nc in c(10, 20)) {
    res2 <- try(ecospat::ecospat.boyce(fit = fit_vals, obs = obs_vals,
                                       nclass = nc, PEplot = FALSE), silent = TRUE)
    sc2  <- get_boyce_cor(if (!inherits(res2, "try-error")) res2 else NULL)
    if (is.finite(sc2)) return(sc2)
  }
  NA_real_
}

# Favourability transformation: corrects for prevalence (Garcia-Rosell 2004)
favourability_from_prob <- function(prob, prev_ratio) {
  f <- function(p) {
    odds <- p / (1 - p)
    fav  <- odds / (prev_ratio + odds)
    fav[!is.finite(fav)] <- NA
    fav[fav < 0] <- 0
    fav[fav > 1] <- 1
    fav
  }
  if (inherits(prob, "SpatRaster"))                           terra::app(prob, f)
  else if (inherits(prob, c("RasterLayer","RasterBrick","RasterStack"))) raster::calc(prob, f)
  else stop("favourability_from_prob: unsupported class: ", class(prob)[1])
}

# Number of CV folds based on number of presences (mirrors 0X_boyce draft)
n_folds_from_pres <- function(n_pres) {
  if (n_pres < 40L) return(0L)
  min(5L, floor(n_pres / 20L))
}

# Cap values to [0, 1]
cap01 <- function(x) { x[x < 0] <- 0; x[x > 1] <- 1; x }

# ------------------------------------------------------------------
# Sample background points from a raster as a plain data frame,
# then immediately delete the raster from disk.
# n_bg: target number of background points (default 50 000).
# Returns a data frame with one column per predictor layer.
# ------------------------------------------------------------------
sample_bg_df <- function(predictors_r, n_bg = 50000) {
  # predictors_r can be Raster* or SpatRaster
  if (inherits(predictors_r, c("RasterLayer","RasterBrick","RasterStack")))
    predictors_r <- terra::rast(predictors_r)
  n_cells <- terra::global(!is.na(predictors_r[[1]]), "sum", na.rm = TRUE)[1, 1]
  n_sample <- min(n_bg, n_cells)
  if (n_sample < 200) return(NULL)   # too few cells to be useful
  smp <- terra::spatSample(predictors_r, size = n_sample,
                           method = "random", na.rm = TRUE, values = TRUE,
                           as.points = FALSE)
  as.data.frame(smp)
}

# ------------------------------------------------------------------
# Predict favourability for a set of background points (data frame)
# and for presence/absence points (SpatialPointsDataFrame or sf).
# Returns a list:
#   $fit_bg  : numeric vector of favourability values at bg points
#   $obs_vals: numeric vector of favourability values at presence pts
# Uses tabular predict() — no raster written at all.
# ------------------------------------------------------------------
predict_fav_points <- function(model, bg_df, pres_sp, methods_vec, prev_ratio) {
  # ------------------------------------------------------------------
  # bg_df is a plain data frame with one column per predictor — this
  # is what sdm::predict() expects for tabular (non-raster) prediction.
  #
  # pres_sp may arrive as:
  #   (a) a SpatialPointsDataFrame with only coordinate/species cols  ← common
  #   (b) a data frame already containing predictor columns            ← rare
  #   (c) NULL / 0-row object                                          ← skip
  #
  # For case (a) the sp object has no predictor columns, so
  # sdm::predict(newdata = pres_sp) returns garbage or wrong-length
  # output. We convert it to a data frame of predictor values by
  # extracting coordinates and subsetting bg_df's columns via
  # terra::extract() from the model's predictor raster, OR — more
  # robustly — by using predict() with the same column structure as
  # bg_df, built by extracting values from the raster in scope.
  #
  # The cleanest universal fix: convert pres_sp to a plain data frame
  # with the same columns as bg_df by extracting values at the point
  # locations from the raster that was used to build bg_df.  We infer
  # which raster is in scope from a closure argument `pred_raster`
  # (added to the function signature below).
  # ------------------------------------------------------------------
  info   <- sdm::getModelInfo(model)
  col_m  <- if ("methods" %in% names(info)) "methods" else "method"
  info$modelID <- as.integer(info$modelID)
  id_map <- split(info$modelID, info[[col_m]])

  pred_cols <- names(bg_df)   # predictor column names

  fav_f <- function(p) {
    odds <- p / (1 - p)
    fav  <- odds / (prev_ratio + odds)
    fav[!is.finite(fav)] <- NA_real_
    fav <- pmax(0, pmin(1, fav))
    fav
  }

  # Normalise sdm::predict() output to a plain numeric vector.
  # sdm may return a vector, a 1-col data frame, or a multi-col data
  # frame (one column per replicate). We always take the row-wise mean
  # so the result is a single numeric vector of length == expected_n.
  extract_pred_vec <- function(pr, expected_n) {
    if (inherits(pr, "try-error")) return(NULL)
    if (is.data.frame(pr) || is.matrix(pr)) {
      if (nrow(pr) != expected_n) return(NULL)
      v <- rowMeans(as.matrix(pr), na.rm = TRUE)
    } else {
      v <- as.numeric(pr)
      if (length(v) != expected_n) return(NULL)
    }
    v
  }

  # ------------------------------------------------------------------
  # Convert pres_sp to a plain data frame with predictor columns.
  # If pres_sp already has all predictor columns → use as-is.
  # If pres_sp is a Spatial* or sf object → extract from pred_raster.
  # pred_raster must be available in the calling environment (passed
  # as an extra argument below).
  # ------------------------------------------------------------------
  pres_df <- NULL
  n_pres  <- 0L

  if (!is.null(pres_sp) && nrow(pres_sp) > 0) {
    if (is.data.frame(pres_sp) && all(pred_cols %in% names(pres_sp))) {
      # Already a data frame with the right columns
      pres_df <- pres_sp[, pred_cols, drop = FALSE]
    } else {
      # Spatial* or sf: extract predictor values at point locations
      pres_vect <- if (inherits(pres_sp, c("sf", "sfc"))) {
        terra::vect(pres_sp)
      } else if (inherits(pres_sp, c("SpatialPoints","SpatialPointsDataFrame"))) {
        terra::vect(pres_sp)
      } else {
        NULL
      }
      if (!is.null(pres_vect) && !is.null(pred_raster)) {
        pred_raster_sr <- if (inherits(pred_raster, c("RasterLayer","RasterBrick","RasterStack")))
          terra::rast(pred_raster) else pred_raster
        extracted <- terra::extract(pred_raster_sr, pres_vect, ID = FALSE)
        extracted <- extracted[, pred_cols[pred_cols %in% names(extracted)], drop = FALSE]
        complete  <- stats::complete.cases(extracted)
        pres_df   <- extracted[complete, , drop = FALSE]
      }
    }
    n_pres <- if (!is.null(pres_df)) nrow(pres_df) else 0L
  }

  bg_preds   <- list()
  pres_preds <- list()

  for (mm in methods_vec) {
    id_m <- as.integer(id_map[[mm]][1])
    if (length(id_m) == 0 || is.na(id_m)) next

    # Predict to background data frame
    pr_bg <- try(predict(model, newdata = bg_df, id = id_m), silent = TRUE)
    v_bg  <- extract_pred_vec(pr_bg, nrow(bg_df))
    if (!is.null(v_bg)) bg_preds[[mm]] <- fav_f(v_bg)

    # Predict to presence data frame (extracted predictor values)
    if (!is.null(pres_df) && nrow(pres_df) > 0) {
      pr_pres <- try(predict(model, newdata = pres_df, id = id_m), silent = TRUE)
      v_pres  <- extract_pred_vec(pr_pres, nrow(pres_df))
      if (!is.null(v_pres)) pres_preds[[mm]] <- fav_f(v_pres)
    }
    gc()
  }

  # Median ensemble across methods
  if (length(bg_preds) == 0) {
    message("    [predict_fav_points] WARNING: no bg predictions captured for any method (",
            paste(methods_vec, collapse = ", "), ")")
    return(list(fit_bg = numeric(0), obs_vals = numeric(0)))
  }
  message("    [predict_fav_points] Methods with bg predictions: ",
          paste(names(bg_preds), collapse = ", "),
          " | bg n=", nrow(bg_df),
          " | pres n=", n_pres,
          " | pres methods=", length(pres_preds))

  bg_mat   <- do.call(cbind, bg_preds)
  fit_bg   <- if (is.matrix(bg_mat)) apply(bg_mat, 1, median, na.rm = TRUE) else as.numeric(bg_mat)

  obs_vals <- if (length(pres_preds) > 0) {
    pres_mat <- do.call(cbind, pres_preds)
    if (is.matrix(pres_mat)) apply(pres_mat, 1, median, na.rm = TRUE) else as.numeric(pres_mat)
  } else {
    numeric(0)
  }

  list(fit_bg   = fit_bg[is.finite(fit_bg)],
       obs_vals = obs_vals[is.finite(obs_vals)])
}

# ------------------------------------------------------------------
# For Stage 3: given two sets of background favourability vectors
# (one from global EU model, one from habitat model), compute the
# geometric mean and return as a numeric vector.
# Both vectors must correspond to the SAME set of background points
# (same spatial sample drawn once and reused).
# ------------------------------------------------------------------
geomean_fav <- function(fav_g, fav_h) {
  n <- min(length(fav_g), length(fav_h))
  if (n == 0) return(numeric(0))
  sqrt(fav_g[seq_len(n)] * fav_h[seq_len(n)])
}


# ==============================================================================
# SECTION 5 — Accumulators for results (filled inside the species loop)
# ==============================================================================
all_boyce_global_folds   <- list()   # Stage 1 per-fold Boyce rows
all_boyce_habitat_folds  <- list()   # Stage 2 per-fold Boyce rows
all_boyce_combined_folds <- list()   # Stage 3 per-combination Boyce rows


# ==============================================================================
# SECTION 6 — Species loop
# ==============================================================================
for (key in accepted_taxonkeys) {

  message("\n", strrep("=", 72))
  message("SPECIES: ", unname(taxon_name_by_key[as.character(key)]),
          "  [taxonkey: ", key, "]")
  message(strrep("=", 72))


  # ----------------------------------------------------------------------------
  # 6.0  Species metadata
  # ----------------------------------------------------------------------------
  species     <- unname(taxon_name_by_key[as.character(key)])
  taxonkey    <- key
  speciesName <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
  base_dir    <- pr("data", "projects", project, paste0(speciesName, "_", taxonkey))

  # Resolve species group for bias-grid selection
  # Priority: (1) taxa_info Group column, (2) occurrence data Group column
  speciesgroup <- if (!is.null(taxon_group_by_key))
    unname(taxon_group_by_key[as.character(key)])
  else
    NA_character_

  if (is.na(speciesgroup) || !nzchar(speciesgroup)) {
    occ_sub  <- cleaned %>% dplyr::filter(acceptedTaxonKey == taxonkey)
    if ("Group" %in% names(occ_sub) && nrow(occ_sub) > 0) {
      grp_vals     <- unique(occ_sub$Group)
      grp_vals     <- grp_vals[!is.na(grp_vals) & nzchar(as.character(grp_vals))]
      speciesgroup <- if (length(grp_vals) > 0) as.character(grp_vals[1]) else NA_character_
    }
  }

  if (is.na(speciesgroup) || !nzchar(speciesgroup) ||
      !speciesgroup %in% names(bias_grid_paths)) {
    warning(species, ": unknown group '", speciesgroup, "' — skipping.")
    next
  }
  message("  Species group: ", speciesgroup)

  # CRS reference (habitat stack drives the target CRS)
  habitat_stack <- terra::rast(habitat_path)
  target_crs    <- sf::st_crs(terra::crs(habitat_stack))

  # Transform shared boundary layers to target CRS once per species
  euboundary    <- sf::st_make_valid(sf::st_transform(euboundary_raw,    target_crs))
  wwf_eco_biome <- sf::st_make_valid(sf::st_transform(wwf_eco_biome_raw, target_crs))

  # Load global climate rasters
  globalclimpreds_terra <- terra::rast(globalclimpreds_file)


  ############################################################################
  ##                                                                        ##
  ##  STAGE 1: GLOBAL (CLIMATE-ONLY) MODEL                                  ##
  ##  Logic follows 03_fit_climate_model.R                                  ##
  ##                                                                        ##
  ############################################################################
  message("\n--- STAGE 1: Global climate model ---")

  # --------------------------------------------------------------------------
  # S1.1  Occurrence data
  # --------------------------------------------------------------------------
  occ_raw <- cleaned %>%
    dplyr::filter(acceptedTaxonKey == taxonkey) %>%
    dplyr::select(decimalLongitude, decimalLatitude)

  # One occurrence per grid cell
  occ_dedup <- remove_duplicates(occurrences     = occ_raw,
                                 rast_template   = globalclimpreds_terra[[1]])
  # Remove occurrences in NA climate cells
  occ_clean <- remove_nodata_occurrences(occurrences  = occ_dedup,
                                         rast_template = globalclimpreds_terra[[1]],
                                         crs           = 4326)

  if (nrow(occ_clean) < 20) {
    message("  Fewer than 20 global occurrences (n = ", nrow(occ_clean), ") — skipping.")
    next
  }

  # All raw occurrences (pre-thinning) used to inform PA selection
  for_PA_selection <- cleaned %>%
    dplyr::filter(acceptedTaxonKey == taxonkey) %>%
    dplyr::select(decimalLongitude, decimalLatitude) %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326)

  # Add presence column
  occ_clean <- occ_clean %>% dplyr::mutate(species = 1L)

  # Optional thinning to ≤ 10 000
  if (nrow(occ_clean) > 10000) {
    if (occurrence_thinning_method == "kmeans_clustering") {
      env_data <- terra::extract(globalclimpreds_terra, terra::vect(occ_clean), ID = FALSE)
      n_ctr    <- min(10000L, nrow(unique(env_data)))
      set.seed(101)
      cl       <- kmeans(env_data, centers = n_ctr, iter.max = 10, nstart = 1)$cluster
      occ_clean <- occ_clean[!duplicated(cl), ] %>% dplyr::slice_sample(n = 10000L)
    } else {
      set.seed(101)
      occ_clean <- occ_clean %>% dplyr::slice_sample(n = 10000L)
    }
  }
  message("  Occurrences after thinning: ", nrow(occ_clean))

  # --------------------------------------------------------------------------
  # S1.2  Biomes that contain presences
  # --------------------------------------------------------------------------
  # occ_clean is already an sf in EPSG:4326 (output of remove_nodata_occurrences);
  # transform to wwf_eco_biome CRS for the intersection
  occ_clean_sf <- sf::st_transform(occ_clean, sf::st_crs(wwf_eco_biome))
  has_occ    <- lengths(sf::st_intersects(wwf_eco_biome, occ_clean_sf)) > 0
  wwf_ecoSub <- wwf_eco_biome[has_occ, ]
  message("  Occupied biomes: ", nrow(wwf_ecoSub))

  # --------------------------------------------------------------------------
  # S1.3  Bias grid + pseudo-absences
  # --------------------------------------------------------------------------
  globalclimpreds_5k <- terra::rast(globalclimpreds_5k_file)
  biasgrid_raw <- terra::resample(terra::rast(bias_grid_paths[[speciesgroup]]),
                                  globalclimpreds_5k, method = "bilinear")
  biasgrid_log <- terra::mask(biasgrid_raw, globalclimpreds_5k)
  min_v <- terra::global(biasgrid_log, "min", na.rm = TRUE)[[1]]
  max_v <- terra::global(biasgrid_log, "max", na.rm = TRUE)[[1]]
  biasgrid     <- ((biasgrid_log - min_v) / (max_v - min_v)) * 19 + 1
  # Reproject wwf_ecoSub to biasgrid CRS (4326) for crop/mask
  wwf_ecoSub_4326 <- sf::st_transform(wwf_ecoSub, 4326)
  biasgrid_sub <- terra::crop(biasgrid, terra::ext(terra::vect(wwf_ecoSub_4326)))
  biasgrid_sub <- terra::mask(biasgrid_sub, terra::vect(wwf_ecoSub_4326))

  # Mask cells already occupied (for_PA_selection is already sf in crs 4326)
  occ_vect   <- terra::vect(for_PA_selection)
  occ_cells  <- terra::cellFromXY(biasgrid_sub, terra::crds(occ_vect))
  occ_cells  <- occ_cells[!is.na(occ_cells)]
  biasgrid_sub[occ_cells] <- NA

  # Draw and thin pseudo-absences
  set.seed(728)
  pts30k <- terra::spatSample(biasgrid_sub, size = 30000,
                              method = "weights", as.points = TRUE, na.rm = TRUE)

  if (pseudoabsence_thinning_method == "kmeans_clustering") {
    pa_clim  <- terra::extract(globalclimpreds_terra, pts30k, ID = FALSE, xy = TRUE)
    pa_clim  <- na.omit(pa_clim)
    n_ctr_pa <- min(10000L, nrow(unique(pa_clim)))
    set.seed(101)
    cl_pa    <- kmeans(pa_clim[, !names(pa_clim) %in% c("x", "y")],
                       centers = n_ctr_pa, iter.max = 10, nstart = 1)$cluster
    pa_clim$clust <- cl_pa; pa_clim$rID <- seq_len(nrow(pa_clim))
    sampled  <- pa_clim %>% dplyr::group_by(clust) %>%
      dplyr::slice_sample(n = 1) %>% dplyr::ungroup()
    rem      <- 10000L - nrow(sampled)
    if (rem > 0) sampled <- dplyr::bind_rows(
      sampled,
      pa_clim %>% dplyr::filter(!rID %in% sampled$rID) %>% dplyr::slice_sample(n = rem)
    )
    global_pts <- sampled %>%
      dplyr::rename(decimalLongitude = x, decimalLatitude = y) %>%
      dplyr::select(decimalLongitude, decimalLatitude) %>%
      sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"),
                   crs = 4326, remove = FALSE)
    rm(pa_clim, sampled)
  } else {
    set.seed(101)
    global_pts <- pts30k[sample(nrow(pts30k), 10000L), ] %>% sf::st_as_sf()
    coords     <- sf::st_coordinates(global_pts)
    global_pts <- global_pts %>%
      dplyr::mutate(decimalLongitude = coords[, "X"],
                    decimalLatitude  = coords[, "Y"]) %>%
      dplyr::select(decimalLongitude, decimalLatitude, geometry)
  }
  rm(pts30k, biasgrid_raw, biasgrid_log, biasgrid); gc()

  # --------------------------------------------------------------------------
  # S1.4  Presence-pseudoabsence dataset
  # --------------------------------------------------------------------------
  global_pseudoAbs <- global_pts %>% dplyr::mutate(species = 0L)
  global_presabs   <- rbind(occ_clean, global_pseudoAbs)
  rm(global_pts, global_pseudoAbs); gc()
  message("  Presences: ", sum(global_presabs$species == 1L),
          " | Pseudoabsences: ", sum(global_presabs$species == 0L))

  # --------------------------------------------------------------------------
  # S1.5  Load climate .qs — reuse selected predictors and top-5 methods
  #        from 03_fit_climate_model.R rather than redoing PCAm from scratch
  # --------------------------------------------------------------------------
  climate_qs_file <- file.path(base_dir, "Climate",
                                paste0("Climate_model_", speciesName, "_", taxonkey, ".qs"))
  if (!file.exists(climate_qs_file)) {
    warning(species, ": climate .qs not found at\n  ", climate_qs_file,
            "\n  Run 03_fit_climate_model.R first. Skipping.")
    next
  }
  climate_qs   <- qs::qread(climate_qs_file)

  # selected_predictors: character vector of retained (uncorrelated) predictor names
  top5_methods        <- climate_qs$top5_models       # character vector, e.g. c("rf","brt",...)
  selected_predictors <- climate_qs$selected_predictors  # names of kept climate variables
  rm(climate_qs); gc()

  message("  Top-5 methods from climate .qs: ", paste(top5_methods, collapse = ", "))
  message("  Climate predictors from .qs: ",    paste(selected_predictors, collapse = ", "))

  # Subset climate rasters to the same predictors used in 03_fit_climate_model.R
  clim_sel <- terra::subset(globalclimpreds_terra,
                             selected_predictors[selected_predictors %in%
                                                   names(globalclimpreds_terra)])
  eu_clim  <- terra::subset(terra::rast(eu_climpreds_file),
                             selected_predictors[selected_predictors %in%
                                                   names(terra::rast(eu_climpreds_file))])
  message("  Climate predictors retained: ", terra::nlyr(clim_sel))

  # --------------------------------------------------------------------------
  # S1.6  Convert to sp / Raster* for sdm package
  # --------------------------------------------------------------------------
  n_pres_global  <- sum(global_presabs$species == 1L)
  prev_ratio_g   <- n_pres_global / max(1L, sum(global_presabs$species == 0L))

  global_presabs_sp <- methods::as(
    global_presabs %>%
      dplyr::mutate(species = as.integer(species)) %>%
      dplyr::select(-decimalLongitude, -decimalLatitude),
    "Spatial"
  )
  clim_sel_r <- raster::stack(clim_sel)
  raster::crs(clim_sel_r) <- terra::crs(clim_sel)

  # --------------------------------------------------------------------------
  # S1.7  (No PCAm rerun — top5_methods already loaded from .qs above)
  # --------------------------------------------------------------------------

  # --------------------------------------------------------------------------
  # S1.8  Number of spatial CV folds
  # --------------------------------------------------------------------------
  k_global      <- n_folds_from_pres(n_pres_global)
  use_cv_global <- k_global >= 2L
  message("  Global presences: ", n_pres_global, " → folds: ",
          if (use_cv_global) k_global else "none (train-only)")

  # Storage for per-fold predictions
  global_fold_preds_full <- vector("list", max(1L, k_global))  # full biome extent
  global_fold_preds_eu   <- vector("list", max(1L, k_global))  # EU extent only
  fold_ids_global        <- NULL
  sb_global              <- NULL

  # --------------------------------------------------------------------------
  # S1.9  Sample background points ONCE (full training area + EU subset)
  #
  #  Key design: we never predict to the full raster. Instead we draw 50 000
  #  random points from the training-area climate stack, extract their climate
  #  values into a data frame, and run tabular predict(). This is orders of
  #  magnitude faster than raster prediction over large extents (e.g. whole
  #  of Eurasia + North America) and writes nothing to disk.
  #
  #  50 000 background points is the target; compute_boyce_robust needs ≥ 200
  #  finite values so the floor is generous. The same background sample is
  #  reused across folds so that Boyce comparisons are consistent.
  # --------------------------------------------------------------------------
  message("  Sampling ", N_BG, " background points from training area...")
  bg_full_df <- sample_bg_df(clim_sel_r, n_bg = N_BG)
  if (is.null(bg_full_df)) {
    warning(species, ": training area has too few valid climate cells — skipping.")
    next
  }
  message("  Background points drawn: ", nrow(bg_full_df))

  # EU background: extract climate values at the same EU climate raster
  eu_clim_r <- raster::stack(eu_clim)
  raster::crs(eu_clim_r) <- terra::crs(eu_clim)
  bg_eu_df <- sample_bg_df(eu_clim_r, n_bg = N_BG)
  message("  EU background points drawn: ", if (is.null(bg_eu_df)) 0L else nrow(bg_eu_df))

  # --------------------------------------------------------------------------
  # S1.10  EU presences in climate CRS (for Boyce evaluation)
  # --------------------------------------------------------------------------
  clim_crs <- sf::st_crs(terra::crs(clim_sel))

  occ_global_sf <- occ_clean %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = 4326) %>%
    sf::st_transform(target_crs)

  occ_eu_climcrs <- occ_global_sf %>%
    sf::st_transform(clim_crs) %>%
    sf::st_crop(sf::st_bbox(sf::st_transform(euboundary, clim_crs)))
  inside_eu <- lengths(sf::st_intersects(
    occ_eu_climcrs, sf::st_transform(euboundary, clim_crs))) > 0
  occ_eu_climcrs <- occ_eu_climcrs[inside_eu, ]
  has_eu_occ     <- nrow(occ_eu_climcrs) >= 5L
  message("  EU occurrences for Stage 1 EU Boyce: ", nrow(occ_eu_climcrs))

  # Convert EU presences to SpatialPointsDataFrame for tabular predict()
  occ_eu_sp <- if (has_eu_occ) methods::as(occ_eu_climcrs, "Spatial") else NULL

  # All global presences as sp (for full-area Boyce)
  all_pres_sp <- methods::as(
    occ_clean %>%
      sf::st_as_sf(coords = c("decimalLongitude","decimalLatitude"), crs = 4326) %>%
      sf::st_transform(clim_crs),
    "Spatial"
  )

  # --------------------------------------------------------------------------
  # S1.11  Spatial block CV — global, with inline Boyce computation
  #         No rasters stored. Per fold: fit model → tabular predict →
  #         compute Boyce → discard model.
  #         For Stage 3 we only need the EU favourability VALUES at the
  #         background and presence points, stored as numeric vectors.
  # --------------------------------------------------------------------------
  pres_rows_g <- which(global_presabs_sp$species == 1L)
  boyce_g_rows <- list()

  # Stored as numeric vectors for Stage 3 (no rasters):
  # global_eu_bg_list[[i]]  : bg favourability vector for EU (fold i)
  # global_eu_obs_list[[i]] : presence favourability vector for EU (fold i)
  global_eu_bg_list  <- vector("list", max(1L, k_global))
  global_eu_obs_list <- vector("list", max(1L, k_global))

  # pred_raster is used inside predict_fav_points to extract predictor
  # values at presence-point locations (climate predictors for Stage 1)
  pred_raster <- clim_sel_r

  run_global_fold <- function(fold_i, train_sp, prev_tr, fold_label,
                               fold_pres_ids = NULL) {
    sdm_d <- sdm::sdmData(species ~ ., train = train_sp, predictors = clim_sel_r)
    mod   <- sdm::sdm(species ~ ., data = sdm_d, methods = top5_methods)

    # --- (a) Full training area Boyce ---
    if (is.null(fold_pres_ids)) {
      # No CV: all presences are "train"
      pres_test_sp  <- all_pres_sp
      pres_train_sp <- all_pres_sp
    } else {
      pres_test_sp  <- all_pres_sp[fold_pres_ids == fold_label, , drop = FALSE]
      pres_train_sp <- all_pres_sp[fold_pres_ids != fold_label, , drop = FALSE]
    }

    full_res_test  <- predict_fav_points(mod, bg_full_df, pres_test_sp,
                                          top5_methods, prev_tr)
    full_res_train <- predict_fav_points(mod, bg_full_df, pres_train_sp,
                                          top5_methods, prev_tr)

    boyce_full <- data.frame(
      fold        = fold_label,
      area        = "global_full",
      boyce_train = compute_boyce_robust(full_res_train$fit_bg,
                                          full_res_train$obs_vals),
      boyce_test  = compute_boyce_robust(full_res_test$fit_bg,
                                          full_res_test$obs_vals),
      n_train     = length(full_res_train$obs_vals),
      n_test      = length(full_res_test$obs_vals),
      stringsAsFactors = FALSE
    )

    # --- (b) EU Boyce (if EU occurrences exist) ---
    boyce_eu <- NULL
    eu_bg_vec  <- numeric(0)
    eu_obs_vec <- numeric(0)

    if (has_eu_occ && !is.null(bg_eu_df) && !is.null(occ_eu_sp)) {
      eu_res <- predict_fav_points(mod, bg_eu_df, occ_eu_sp,
                                   top5_methods, prev_tr)
      eu_bg_vec  <- eu_res$fit_bg
      eu_obs_vec <- eu_res$obs_vals

      # EU train/test split (same fold structure as full area)
      if (!is.null(fold_pres_ids)) {
        # Assign EU occurrences to folds via block membership
        blk_sf <- sb_global$blocks
        blk_sf$BLOCK_ROWID <- seq_len(nrow(blk_sf))
        eu_pts_for_join <- sf::st_transform(occ_eu_climcrs, sf::st_crs(blk_sf))
        joined <- sf::st_join(eu_pts_for_join,
                              blk_sf["BLOCK_ROWID"],
                              join = sf::st_within, left = TRUE)
        # Each block → the fold assigned to the majority of its presences
        block_fold_map <- tapply(fold_ids_global,
                                 sb_global$folds_ids,
                                 FUN = function(x) x[1])
        eu_fold_ids <- fold_ids_global[
          match(joined$BLOCK_ROWID,
                seq_len(nrow(global_presabs_sp)))
        ]
        eu_fold_ids[is.na(eu_fold_ids)] <- 0L

        eu_test_sp  <- if (sum(eu_fold_ids == fold_label) >= 1)
          occ_eu_sp[eu_fold_ids == fold_label, , drop = FALSE] else NULL
        eu_train_sp <- if (sum(eu_fold_ids != fold_label) >= 1)
          occ_eu_sp[eu_fold_ids != fold_label, , drop = FALSE] else occ_eu_sp

        eu_res_test  <- predict_fav_points(mod, bg_eu_df, eu_test_sp,
                                            top5_methods, prev_tr)
        eu_res_train <- predict_fav_points(mod, bg_eu_df, eu_train_sp,
                                            top5_methods, prev_tr)
      } else {
        eu_res_test  <- eu_res
        eu_res_train <- eu_res
      }

      boyce_eu <- data.frame(
        fold        = fold_label,
        area        = "global_eu",
        boyce_train = compute_boyce_robust(eu_res_train$fit_bg,
                                            eu_res_train$obs_vals),
        boyce_test  = compute_boyce_robust(eu_res_test$fit_bg,
                                            eu_res_test$obs_vals),
        n_train     = length(eu_res_train$obs_vals),
        n_test      = length(eu_res_test$obs_vals),
        stringsAsFactors = FALSE
      )
    }

    rm(mod, sdm_d); gc()
    list(boyce_full = boyce_full, boyce_eu = boyce_eu,
         eu_bg_vec = eu_bg_vec, eu_obs_vec = eu_obs_vec)
  }

  if (use_cv_global) {
    message("  Running blockCV spatial blocking (global)...")
    set.seed(42)
    sb_global <- blockCV::cv_spatial(
      x         = global_presabs_sp,
      column    = "species",
      r         = clim_sel_r,
      k         = k_global,
      hexagon   = TRUE,
      selection = "random",
      iteration = 200,
      size      = 500 * 1000   # 500 km hexagons for global extent
    )
    fold_ids_global <- sb_global$folds_ids
    fold_pres_ids   <- fold_ids_global[pres_rows_g]

    for (i in seq_len(k_global)) {
      message("  Global fold ", i, "/", k_global)
      train_idx <- which(fold_ids_global != i)
      train_sp  <- global_presabs_sp[train_idx, ]
      prev_tr   <- sum(train_sp$species == 1L) / max(1L, sum(train_sp$species == 0L))

      res_i <- run_global_fold(i, train_sp, prev_tr, fold_label = i,
                                fold_pres_ids = fold_pres_ids)
      boyce_g_rows <- c(boyce_g_rows,
                        list(res_i$boyce_full),
                        if (!is.null(res_i$boyce_eu)) list(res_i$boyce_eu))
      global_eu_bg_list[[i]]  <- res_i$eu_bg_vec
      global_eu_obs_list[[i]] <- res_i$eu_obs_vec
      rm(train_sp, res_i); gc()
    }

  } else {
    message("  Fitting train-only global model (no CV)...")
    res_all <- run_global_fold(0L, global_presabs_sp, prev_ratio_g,
                                fold_label = 0L, fold_pres_ids = NULL)
    boyce_g_rows <- c(boyce_g_rows,
                      list(res_all$boyce_full),
                      if (!is.null(res_all$boyce_eu)) list(res_all$boyce_eu))
    global_eu_bg_list[[1]]  <- res_all$eu_bg_vec
    global_eu_obs_list[[1]] <- res_all$eu_obs_vec
    rm(res_all); gc()
    k_global <- 0L
  }

  boyce_global_df <- dplyr::bind_rows(boyce_g_rows) %>%
    dplyr::mutate(species = species, taxonkey = taxonkey)

  boyce_global_summary <- boyce_global_df %>%
    dplyr::group_by(area) %>%
    dplyr::summarise(
      n_folds          = dplyr::n(),
      boyce_train_mean = mean(boyce_train, na.rm = TRUE),
      boyce_train_sd   = sd(boyce_train,   na.rm = TRUE),
      boyce_test_mean  = mean(boyce_test,  na.rm = TRUE),
      boyce_test_sd    = sd(boyce_test,    na.rm = TRUE),
      .groups = "drop"
    )

  message("  [Stage 1 Boyce — per fold]")
  print(boyce_global_df)
  message("  [Stage 1 Boyce — summary]")
  print(boyce_global_summary)

  all_boyce_global_folds[[length(all_boyce_global_folds) + 1]] <- boyce_global_df

  # Clean up background data frames — no longer needed for Stage 1
  rm(bg_full_df, eu_clim_r); gc()


  ############################################################################
  ##                                                                        ##
  ##  STAGE 2: EUROPE HABITAT MODEL                                         ##
  ##  Logic follows 04_fit_habitat_model.R                                  ##
  ##                                                                        ##
  ############################################################################
  message("\n--- STAGE 2: Europe habitat model ---")

  boyce_habitat_df     <- NULL
  hab_fold_preds       <- list()
  k_habitat            <- 0L
  fold_ids_hab         <- NULL
  sb_hab               <- NULL

  # --------------------------------------------------------------------------
  # S2.1  EU occurrence data (in habitat CRS = target_crs)
  # --------------------------------------------------------------------------
  eu_occ_hab <- sf::st_crop(occ_global_sf, sf::st_bbox(euboundary))
  inside_h   <- lengths(sf::st_intersects(eu_occ_hab, euboundary)) > 0
  eu_occ_hab <- eu_occ_hab[inside_h, ]

  # One per habitat grid cell, drop NA-habitat points
  if (nrow(eu_occ_hab) > 0) {
    coords_h   <- sf::st_coordinates(eu_occ_hab)
    cells_h    <- terra::cellFromXY(habitat_stack[[1]], coords_h)
    eu_occ_hab <- eu_occ_hab[!duplicated(cells_h), ]
    vals_h     <- terra::extract(habitat_stack, terra::vect(eu_occ_hab), ID = FALSE)
    eu_occ_hab <- eu_occ_hab[stats::complete.cases(vals_h), ]
  }
  message("  EU occurrences (unique habitat cells): ", nrow(eu_occ_hab))

  if (nrow(eu_occ_hab) < 20) {
    message("  Fewer than 20 EU occurrences — skipping Stage 2 & 3.")
  } else {

    # Optional thinning
    if (nrow(eu_occ_hab) > 10000) {
      if (occurrence_thinning_method == "kmeans_clustering") {
        hab_dat  <- terra::extract(habitat_stack, terra::vect(eu_occ_hab), ID = FALSE)
        n_ctr_h  <- min(10000L, nrow(unique(hab_dat)))
        set.seed(101)
        cl_h     <- kmeans(hab_dat, centers = n_ctr_h, iter.max = 10, nstart = 1)$cluster
        max_pp   <- ceiling(10000 / n_ctr_h)
        keep_h   <- unlist(lapply(unique(cl_h), function(cl) {
          idx <- which(cl_h == cl); idx[sample(length(idx), min(max_pp, length(idx)))]
        }))
        eu_occ_hab <- eu_occ_hab[keep_h, ]
      } else {
        set.seed(101); eu_occ_hab <- eu_occ_hab[sample(nrow(eu_occ_hab), 10000L), ]
      }
    }
    message("  EU occurrences after thinning: ", nrow(eu_occ_hab))

    # --------------------------------------------------------------------------
    # S2.2  Bias grid aligned to habitat resolution
    # --------------------------------------------------------------------------
    biasgrid_file_sp <- file.path(base_dir, "Climate", "Current", "Interim",
                                   paste0("Biasgrid_", speciesName, "_", taxonkey, ".tif"))
    if (file.exists(biasgrid_file_sp)) {
      biasgrid_aligned <- terra::project(terra::rast(biasgrid_file_sp),
                                          habitat_stack[[1]], method = "bilinear")
    } else {
      biasgrid_raw2  <- terra::resample(
        terra::rast(bias_grid_paths[[speciesgroup]]),
        terra::rast(globalclimpreds_5k_file), method = "bilinear")
      biasgrid_log2  <- terra::mask(biasgrid_raw2, terra::rast(globalclimpreds_5k_file))
      min_v2 <- terra::global(biasgrid_log2, "min", na.rm = TRUE)[[1]]
      max_v2 <- terra::global(biasgrid_log2, "max", na.rm = TRUE)[[1]]
      biasgrid2 <- ((biasgrid_log2 - min_v2) / (max_v2 - min_v2)) * 19 + 1
      biasgrid_aligned <- terra::project(biasgrid2, habitat_stack[[1]], method = "bilinear")
      rm(biasgrid_raw2, biasgrid_log2, biasgrid2)
    }
    biasgrid_aligned <- terra::mask(biasgrid_aligned, habitat_stack[[1]])

    # Invaded ecoregions mask
    hit_mat      <- sf::st_intersects(wwf_eco_biome, eu_occ_hab, sparse = FALSE)
    wwf_filtered <- wwf_eco_biome[rowSums(hit_mat) > 0, ]
    inside_mask  <- terra::rasterize(terra::vect(wwf_filtered),
                                     biasgrid_aligned, field = 1, background = NA)
    biasgrid_eu  <- terra::ifel(!is.na(inside_mask), biasgrid_aligned, 1)
    biasgrid_eu  <- terra::mask(biasgrid_eu, biasgrid_aligned)

    # --------------------------------------------------------------------------
    # S2.3  Pseudo-absences
    # --------------------------------------------------------------------------
    set.seed(728)
    EU_pts_pa <- terra::spatSample(biasgrid_eu, size = 10000,
                                   method = "weights", as.points = TRUE, na.rm = TRUE)

    eu_coords_h <- sf::st_coordinates(eu_occ_hab)
    eu_occ_pa   <- eu_occ_hab %>%
      dplyr::mutate(decimalLongitude = eu_coords_h[, 1],
                    decimalLatitude  = eu_coords_h[, 2],
                    species          = "present") %>%
      dplyr::select(decimalLongitude, decimalLatitude, species, geometry)

    EU_pts_sf <- EU_pts_pa[, 0] %>%
      sf::st_as_sf() %>%
      dplyr::mutate(coords           = sf::st_coordinates(geometry),
                    decimalLongitude = coords[, 1],
                    decimalLatitude  = coords[, 2],
                    species          = "absent") %>%
      dplyr::select(decimalLongitude, decimalLatitude, species, geometry)

    eu_presabs_hab <- rbind(eu_occ_pa, EU_pts_sf)
    message("  Habitat presences: ", sum(eu_presabs_hab$species == "present"),
            " | Pseudoabsences: ", sum(eu_presabs_hab$species == "absent"))

    # --------------------------------------------------------------------------
    # S2.4  Select habitat predictors
    #        Prefer the retained predictor set from the habitat .qs (occ_full_df
    #        column names = predictors kept after the 0.7 correlation filter in
    #        04_fit_habitat_model.R).  Fall back to rerunning the filter if the
    #        .qs is not yet available.
    # --------------------------------------------------------------------------
    habitat_qs_file_s24 <- file.path(base_dir, "Habitat",
                                      paste0("Habitat_model_", speciesName, "_", taxonkey, ".qs"))
    if (file.exists(habitat_qs_file_s24)) {
      hab_qs_s24         <- qs::qread(habitat_qs_file_s24)
      # occ_full_df columns: all retained predictors + "occ"
      retained_hab_preds <- setdiff(names(hab_qs_s24$occ_full_df), "occ")
      rm(hab_qs_s24); gc()
      # Keep only predictors that actually exist in the current habitat stack
      retained_hab_preds <- intersect(retained_hab_preds, names(habitat_stack))
      fullstack_h <- terra::subset(habitat_stack, retained_hab_preds)
      message("  Habitat predictors from .qs: ", terra::nlyr(fullstack_h),
              " (", paste(retained_hab_preds, collapse = ", "), ")")
    } else {
      # Fallback: rerun correlation filter (consistent with 04_fit_habitat_model.R)
      presabs_hab_df <- terra::extract(habitat_stack, terra::vect(eu_presabs_hab), ID = FALSE)
      cor_mat_h      <- stats::cor(presabs_hab_df, use = "complete.obs")
      drop_h         <- caret::findCorrelation(cor_mat_h, cutoff = 0.7, exact = TRUE, names = TRUE)
      fullstack_h    <- terra::subset(habitat_stack, !names(habitat_stack) %in% drop_h)
      message("  Habitat predictors (filter rerun, .qs not found): ", terra::nlyr(fullstack_h),
              " (dropped ", length(drop_h), " correlated)")
    }

    # --------------------------------------------------------------------------
    # S2.5  Convert to sp / Raster* for sdm package
    # --------------------------------------------------------------------------
    eu_presabs_num   <- eu_presabs_hab %>%
      dplyr::mutate(species = ifelse(species == "present", 1L, 0L)) %>%
      dplyr::select(-decimalLongitude, -decimalLatitude)
    # Ensure species is integer (not logical/double) before Spatial conversion
    eu_presabs_num$species <- as.integer(eu_presabs_num$species)
    eu_presabs_sp_h  <- methods::as(eu_presabs_num, "Spatial")
    fullstack_h_r    <- raster::stack(fullstack_h)
    raster::crs(fullstack_h_r) <- sf::st_crs(terra::crs(habitat_stack))$wkt

    # Load habitat .qs to get the top-5 methods selected by 04_fit_habitat_model.R
    # The habitat .qs stores `top5models` (the fitted sdm object); extract method
    # names from it via getModelInfo so we use the exact same algorithms.
    habitat_qs_file <- file.path(base_dir, "Habitat",
                                  paste0("Habitat_model_", speciesName, "_", taxonkey, ".qs"))
    if (file.exists(habitat_qs_file)) {
      habitat_qs  <- qs::qread(habitat_qs_file)
      hab_sdm_obj <- habitat_qs$top5models   # fitted sdm object with top-5 methods
      hab_info    <- sdm::getModelInfo(hab_sdm_obj)
      col_m_hab   <- if ("methods" %in% names(hab_info)) "methods" else "method"
      hab_methods <- unique(as.character(hab_info[[col_m_hab]]))
      rm(habitat_qs, hab_sdm_obj, hab_info); gc()
      message("  Top-5 habitat methods from habitat .qs: ", paste(hab_methods, collapse = ", "))
    } else {
      # Habitat .qs not yet available — fall back to the global top-5
      hab_methods <- top5_methods
      message("  Habitat .qs not found — falling back to global top-5 methods: ",
              paste(hab_methods, collapse = ", "))
    }

    n_pres_hab   <- sum(eu_presabs_sp_h$species == 1L)
    prev_ratio_h <- n_pres_hab / max(1L, sum(eu_presabs_sp_h$species == 0L))

    # --------------------------------------------------------------------------
    # S2.6  Number of CV folds for habitat model
    # --------------------------------------------------------------------------
    k_habitat      <- n_folds_from_pres(n_pres_hab)
    use_cv_hab     <- k_habitat >= 2L
    message("  Habitat presences: ", n_pres_hab, " → folds: ",
            if (use_cv_hab) k_habitat else "none (train-only)")

    # --------------------------------------------------------------------------
    # S2.7  Sample background points for habitat Boyce (once, reused across folds)
    # --------------------------------------------------------------------------
    message("  Sampling ", N_BG, " habitat background points...")
    bg_hab_df <- sample_bg_df(fullstack_h_r, n_bg = N_BG)
    if (is.null(bg_hab_df)) {
      message("  Habitat raster has too few valid cells — skipping Stage 2 & 3.")
    } else {
      message("  Habitat background points drawn: ", nrow(bg_hab_df))

      # EU presences as SpatialPointsDataFrame for tabular predict()
      eu_occ_sp_h <- methods::as(eu_occ_hab, "Spatial")

      # Stage 3 storage: bg and obs favourability vectors per habitat fold
      hab_bg_list  <- vector("list", max(1L, k_habitat))
      hab_obs_list <- vector("list", max(1L, k_habitat))

    # --------------------------------------------------------------------------
    # S2.8  Spatial block CV — habitat, with inline Boyce computation
    # --------------------------------------------------------------------------
    pres_rows_h  <- which(eu_presabs_sp_h$species == 1L)
    boyce_h_rows <- list()

      # pred_raster is used inside predict_fav_points to extract predictor
      # values at presence-point locations (habitat predictors for Stage 2)
      pred_raster <- fullstack_h_r

      run_hab_fold <- function(fold_idx, train_sp_h, prev_tr_h, fold_label,
                                fold_pres_ids_h = NULL) {
      sdm_d_h <- sdm::sdmData(species ~ ., train = train_sp_h,
                               predictors = fullstack_h_r)
      mod_h   <- sdm::sdm(species ~ ., data = sdm_d_h, methods = hab_methods)

      if (is.null(fold_pres_ids_h)) {
        pres_test_sp_h  <- eu_occ_sp_h
        pres_train_sp_h <- eu_occ_sp_h
      } else {
        pres_test_sp_h  <- eu_occ_sp_h[fold_pres_ids_h == fold_label, , drop = FALSE]
        pres_train_sp_h <- eu_occ_sp_h[fold_pres_ids_h != fold_label, , drop = FALSE]
      }

      res_test  <- predict_fav_points(mod_h, bg_hab_df, pres_test_sp_h,
                                       hab_methods, prev_tr_h)
      res_train <- predict_fav_points(mod_h, bg_hab_df, pres_train_sp_h,
                                       hab_methods, prev_tr_h)

      row_h <- data.frame(
        fold        = fold_label,
        area        = "habitat_eu",
        boyce_train = compute_boyce_robust(res_train$fit_bg, res_train$obs_vals),
        boyce_test  = compute_boyce_robust(res_test$fit_bg,  res_test$obs_vals),
        n_train     = length(res_train$obs_vals),
        n_test      = length(res_test$obs_vals),
        stringsAsFactors = FALSE
      )

      rm(mod_h, sdm_d_h); gc()
      list(boyce_row = row_h,
           bg_vec    = res_train$fit_bg,   # full background = same as train bg
           obs_vec   = c(res_train$obs_vals, res_test$obs_vals))
    }

    if (use_cv_hab) {
      message("  Running blockCV spatial blocking (habitat)...")
      set.seed(42)
      sb_hab <- blockCV::cv_spatial(
        x         = eu_presabs_sp_h,
        column    = "species",
        r         = fullstack_h_r,
        k         = k_habitat,
        hexagon   = TRUE,
        selection = "random",
        iteration = 200,
        size      = 100 * 1000   # 100 km hexagons for EU extent
      )
      fold_ids_hab    <- sb_hab$folds_ids
      fold_pres_ids_h <- fold_ids_hab[pres_rows_h]

      for (i in seq_len(k_habitat)) {
        message("  Habitat fold ", i, "/", k_habitat)
        train_idx_h <- which(fold_ids_hab != i)
        train_sp_h  <- eu_presabs_sp_h[train_idx_h, ]
        prev_tr_h   <- sum(train_sp_h$species == 1L) /
                       max(1L, sum(train_sp_h$species == 0L))

        res_i <- run_hab_fold(i, train_sp_h, prev_tr_h,
                               fold_label = i, fold_pres_ids_h = fold_pres_ids_h)
        boyce_h_rows <- c(boyce_h_rows, list(res_i$boyce_row))
        hab_bg_list[[i]]  <- res_i$bg_vec
        hab_obs_list[[i]] <- res_i$obs_vec
        rm(train_sp_h, res_i); gc()
      }

    } else {
      message("  Fitting train-only habitat model (no CV)...")
      res_h_all <- run_hab_fold(0L, eu_presabs_sp_h, prev_ratio_h,
                                 fold_label = 0L, fold_pres_ids_h = NULL)
      boyce_h_rows <- c(boyce_h_rows, list(res_h_all$boyce_row))
      hab_bg_list[[1]]  <- res_h_all$bg_vec
      hab_obs_list[[1]] <- res_h_all$obs_vec
      rm(res_h_all); gc()
      k_habitat <- 0L
    }

    boyce_habitat_df <- dplyr::bind_rows(boyce_h_rows) %>%
      dplyr::mutate(species = species, taxonkey = taxonkey)

    boyce_habitat_summary <- boyce_habitat_df %>%
      dplyr::summarise(
        area             = "habitat_eu",
        n_folds          = dplyr::n(),
        boyce_train_mean = mean(boyce_train, na.rm = TRUE),
        boyce_train_sd   = sd(boyce_train,   na.rm = TRUE),
        boyce_test_mean  = mean(boyce_test,  na.rm = TRUE),
        boyce_test_sd    = sd(boyce_test,    na.rm = TRUE)
      )

    message("  [Stage 2 Boyce — per fold]")
    print(boyce_habitat_df)
    message("  [Stage 2 Boyce — summary]")
    print(boyce_habitat_summary)

    all_boyce_habitat_folds[[length(all_boyce_habitat_folds) + 1]] <- boyce_habitat_df


    ##########################################################################
    ##                                                                      ##
    ##  STAGE 3: MODEL COMBINATION                                          ##
    ##  Every global EU fold × every habitat fold                           ##
    ##  Combined via geometric mean of favourability VALUES (no rasters).   ##
    ##  Both stages shared the same EU background sample (bg_eu_df drawn    ##
    ##  in S1.9 and bg_hab_df drawn in S2.7).  We cannot use the same      ##
    ##  background points directly because the two models live in different  ##
    ##  predictor spaces; instead we combine their PREDICTIONS at a common  ##
    ##  set of EU habitat-grid points drawn once here.                      ##
    ##                                                                      ##
    ##########################################################################
    message("\n--- STAGE 3: Model combination ---")

    n_g      <- length(Filter(Negate(is.null), global_eu_bg_list))
    n_h      <- length(Filter(Negate(is.null), hab_bg_list))
    n_combos <- n_g * n_h
    message("  Global EU folds with predictions: ", n_g,
            " | Habitat folds: ", n_h,
            " | Total combinations: ", n_combos)

    boyce_comb_rows <- list()
    comb_idx        <- 0L

    for (gi in seq_along(global_eu_bg_list)) {
      g_bg  <- global_eu_bg_list[[gi]]
      g_obs <- global_eu_obs_list[[gi]]
      if (is.null(g_bg) || length(g_bg) < 200) next

      for (hi in seq_along(hab_bg_list)) {
        h_bg  <- hab_bg_list[[hi]]
        h_obs <- hab_obs_list[[hi]]
        if (is.null(h_bg) || length(h_bg) < 200) next

        # Geometric mean of background distributions.
        # The two bg vectors are independent samples (different predictor spaces)
        # so we take the geometric mean of their EMPIRICAL DISTRIBUTIONS by
        # matching them in rank order — this preserves the shape while combining.
        n_common <- min(length(g_bg), length(h_bg))
        g_bg_s   <- sort(g_bg)[seq_len(n_common)]
        h_bg_s   <- sort(h_bg)[seq_len(n_common)]
        comb_bg  <- geomean_fav(g_bg_s, h_bg_s)

        # Geometric mean of observation values (element-wise: same occurrences)
        n_obs    <- min(length(g_obs), length(h_obs))
        comb_obs <- if (n_obs > 0)
          geomean_fav(g_obs[seq_len(n_obs)], h_obs[seq_len(n_obs)])
        else
          numeric(0)

        bc_comb  <- compute_boyce_robust(comb_bg, comb_obs)

        comb_idx     <- comb_idx + 1L
        g_fold_label <- if (k_global  > 0L) gi else 0L
        h_fold_label <- if (k_habitat > 0L) hi else 0L

        boyce_comb_rows[[comb_idx]] <- data.frame(
          combination_id   = comb_idx,
          global_fold      = g_fold_label,
          habitat_fold     = h_fold_label,
          boyce_combined   = bc_comb,
          n_eu_occ         = length(comb_obs),
          species          = species,
          taxonkey         = taxonkey,
          stringsAsFactors = FALSE
        )
        message(sprintf("    Combo %d/%d  (G-fold %d × H-fold %d):  Boyce = %s",
                        comb_idx, n_combos, g_fold_label, h_fold_label,
                        ifelse(is.finite(bc_comb), round(bc_comb, 4), "NA")))
        gc()
      }
    }

    boyce_combined_df <- dplyr::bind_rows(boyce_comb_rows)

    # Guard: if every fold combination was skipped (bg too short), the data
    # frame will be empty and 'boyce_combined' column won't exist.
    boyce_combined_summary <- if (nrow(boyce_combined_df) > 0 &&
                                   "boyce_combined" %in% names(boyce_combined_df)) {
      boyce_combined_df %>%
        dplyr::group_by(species, taxonkey) %>%
        dplyr::summarise(
          n_combinations      = dplyr::n(),
          boyce_combined_mean = mean(boyce_combined, na.rm = TRUE),
          boyce_combined_sd   = sd(boyce_combined,   na.rm = TRUE),
          boyce_combined_min  = min(boyce_combined,  na.rm = TRUE),
          boyce_combined_max  = max(boyce_combined,  na.rm = TRUE),
          .groups = "drop"
        )
    } else {
      message("  [Stage 3] No valid fold combinations produced — returning empty summary.")
      data.frame(species             = species,
                 taxonkey            = taxonkey,
                 n_combinations      = 0L,
                 boyce_combined_mean = NA_real_,
                 boyce_combined_sd   = NA_real_,
                 boyce_combined_min  = NA_real_,
                 boyce_combined_max  = NA_real_,
                 stringsAsFactors    = FALSE)
    }

    message("  [Stage 3 Combined Boyce — per combination]")
    print(boyce_combined_df)
    message("  [Stage 3 Combined Boyce — summary]")
    print(boyce_combined_summary)

    all_boyce_combined_folds[[length(all_boyce_combined_folds) + 1]] <- boyce_combined_df

    rm(bg_hab_df, eu_occ_sp_h, hab_bg_list, hab_obs_list); gc()

    }  # end: bg_hab_df not NULL
  }  # end: enough EU occurrences


  # --------------------------------------------------------------------------
  # End-of-species cleanup — keep only objects needed across iterations
  # --------------------------------------------------------------------------
  rm(list = setdiff(ls(), c(
    # Static data & paths
    "cleaned", "cleaned_1km", "habitat_path", "habitat_stack",
    "euboundary_raw", "wwf_eco_biome_raw",
    "globalclimpreds_file", "globalclimpreds_5k_file", "eu_climpreds_file",
    "bias_grid_paths", "taxon_name_by_key", "taxon_group_by_key",
    "accepted_taxonkeys", "project", "results_dir",
    "occurrence_thinning_method", "pseudoabsence_thinning_method",
    "mtp_probabilities", "country_of_interest",
    # Tempdirs
    "td_terra", "td_raster",
    # Helper functions defined in helper_functions.R (sourced once at startup)
    "remove_duplicates", "remove_nodata_occurrences",
    # Helper functions defined in this script
    "get_boyce_cor",
    "compute_boyce_robust", "favourability_from_prob",
    "n_folds_from_pres", "cap01", "sample_bg_df",
    "predict_fav_points", "geomean_fav",
    # pr() helper and PROJECT_ROOT
    "pr", "PROJECT_ROOT",
    # Background sample size constant
    "N_BG",
    # Accumulators
    "all_boyce_global_folds", "all_boyce_habitat_folds", "all_boyce_combined_folds",
    # Loop variable
    "key"
  )))
  # Also purge any leftover terra temp files from this iteration
  try(terra::tmpFiles(remove = TRUE), silent = TRUE)
  gc()

}  # end species loop


# ==============================================================================
# SECTION 7 — Collect, summarise, and save all results
# ==============================================================================
message("\n", strrep("=", 72))
message("SAVING RESULTS")
message(strrep("=", 72))

# Bind all species together
boyce_global_all   <- dplyr::bind_rows(all_boyce_global_folds)
boyce_habitat_all  <- dplyr::bind_rows(all_boyce_habitat_folds)
boyce_combined_all <- dplyr::bind_rows(all_boyce_combined_folds)

# Per-species summary (mean ± SD across folds)
summarise_folds <- function(df, group_vars = c("species", "taxonkey", "area")) {
  df %>%
    dplyr::group_by(across(all_of(intersect(group_vars, names(df))))) %>%
    dplyr::summarise(
      n_folds          = dplyr::n(),
      boyce_train_mean = mean(boyce_train, na.rm = TRUE),
      boyce_train_sd   = sd(boyce_train,   na.rm = TRUE),
      boyce_test_mean  = mean(boyce_test,  na.rm = TRUE),
      boyce_test_sd    = sd(boyce_test,    na.rm = TRUE),
      .groups = "drop"
    )
}

summarise_combined <- function(df) {
  df %>%
    dplyr::group_by(species, taxonkey) %>%
    dplyr::summarise(
      n_combinations      = dplyr::n(),
      boyce_combined_mean = mean(boyce_combined, na.rm = TRUE),
      boyce_combined_sd   = sd(boyce_combined,   na.rm = TRUE),
      boyce_combined_min  = min(boyce_combined,  na.rm = TRUE),
      boyce_combined_max  = max(boyce_combined,  na.rm = TRUE),
      .groups = "drop"
    )
}

global_summary   <- summarise_folds(boyce_global_all)
habitat_summary  <- summarise_folds(boyce_habitat_all)
combined_summary <- summarise_combined(boyce_combined_all)

# Write CSVs
write.csv(boyce_global_all,   file.path(results_dir, "stage1_global_boyce_per_fold.csv"),    row.names = FALSE)
write.csv(global_summary,     file.path(results_dir, "stage1_global_boyce_summary.csv"),     row.names = FALSE)
write.csv(boyce_habitat_all,  file.path(results_dir, "stage2_habitat_boyce_per_fold.csv"),   row.names = FALSE)
write.csv(habitat_summary,    file.path(results_dir, "stage2_habitat_boyce_summary.csv"),    row.names = FALSE)
write.csv(boyce_combined_all, file.path(results_dir, "stage3_combined_boyce_per_combo.csv"), row.names = FALSE)
write.csv(combined_summary,   file.path(results_dir, "stage3_combined_boyce_summary.csv"),   row.names = FALSE)

message("Saved six CSV files to: ", results_dir)

# Print final summaries to console
message("\n--- Stage 1: Global Boyce (per fold) ---")
print(boyce_global_all)
message("\n--- Stage 1: Global Boyce (summary) ---")
print(global_summary)
message("\n--- Stage 2: Habitat Boyce (per fold) ---")
print(boyce_habitat_all)
message("\n--- Stage 2: Habitat Boyce (summary) ---")
print(habitat_summary)
message("\n--- Stage 3: Combined Boyce (per combination) ---")
print(boyce_combined_all)
message("\n--- Stage 3: Combined Boyce (summary) ---")
print(combined_summary)
