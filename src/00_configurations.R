#--------------------------------------------
#-------  Define wiSDM configurations -------
#--------------------------------------------
project <- "PA prob & Alternative Treshold & Ensemble Boyce"  #TIP: don't use a long project name to avoid Windows errors

species_to_model <- c("Verbena bonariensis", "Erigeron karvinskianus", "Leycesteria formosa", "Asclepias syriaca L.")

occurrence_thinning_method <- "kmeans_clustering" #either "random" or "kmeans_clustering"

mtp_probabilities <- c(0.01, 0.05) #Define MTP thresholds (0.01 = 1%; 0.05 = 5%,...)

country_of_interest <- "Belgium"


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
