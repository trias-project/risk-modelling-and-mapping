#--------------------------------------------
#-------  Define wiSDM configurations -------
#--------------------------------------------
project <- "test_ias_02"  #TIP: don't use a long project name to avoid Windows errors

species_to_model <- c("Arthurdendyus triangulatus")

occurrence_thinning_method <- "kmeans_clustering" #either "random" or "kmeans_clustering"

pseudoabsence_thinning_method <- "kmeans_clustering" #either "random" or "kmeans_clustering"

mtp_probabilities <- c(0.01, 0.05, 0.1) #Define MTP thresholds (0.01 = 1%; 0.05 = 5%,...)

custom_eu_boundary_path <- "./data/external/gadm/europe_selected_countries_wgs84cea_v3-1.gpkg"  #either NULL or a path to a custom Europe boundary vector layer

custom_country_boundary_path <- NULL  #either NULL or a path to a shapefile, if not NULL this overwrites country_of_interest

country_of_interest <- "Europe"

update_files <- "no" #either "ask", "yes", "no"

workflow <-"single_step" #either single_step or two_step

boyce_background_size <- 50000 #Number of non NA pixels in Europe to be selected for Boyce index calculation

user_specific_climate_data <- "./data/external/file_paths_custom_data.csv" #either NULL or a path to a CSV manifest with user-supplied climate rasters

user_specific_landcover_data <- "./data/external/file_paths_custom_landcover_data.csv" #either NULL or a path to a CSV manifest with user-supplied land-cover rasters

habitat_filter_near_zero_variance_predictors <- FALSE #either TRUE or FALSE; if TRUE, near-zero-variance habitat predictors are removed before model fitting



#------------------------------------------------------------
#!ONLY ONCE: Add your GBIF credentials to your ./Renviron file
#------------------------------------------------------------

required_vars <- c("GBIF_USER", "GBIF_PWD", "GBIF_EMAIL")
missing_vars <- required_vars[!nzchar(Sys.getenv(required_vars))]

if (length(missing_vars) > 0) {
  
  example_lines <- paste0(
    if ("GBIF_USER" %in% missing_vars) 'GBIF_USER="your_gbif_username"\n' else "",
    if ("GBIF_PWD" %in% missing_vars)  'GBIF_PWD="your_gbif_password"\n' else "",
    if ("GBIF_EMAIL" %in% missing_vars) 'GBIF_EMAIL="your_email"\n' else ""
  )
  
  message(
    "The following environment variables are missing from your ~/.Renviron file: ",
    paste(missing_vars, collapse = ", "), "\n\n",
    example_lines, "\n",
    "Please add them to your ~/.Renviron file.\n\n",
    "You can open this file by running:\n",
    "  usethis::edit_r_environ()\n"
  )
}

update_files <- tolower(update_files)
