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

#Install usethis if not installed yet
if (!"usethis" %in% rownames(installed.packages()) ) { install.packages("usethis")}
  
# run once in R to open ~/.Renviron for editing
usethis::edit_r_environ()

#Add the following lines to the ./Renviron file, adjust the fields to your specific credentials
GBIF_USER = "your_gbif_username"
GBIF_PWD = "your_gbif_password"
GBIF_EMAIL = "your_email"

#Restart R for changes to take effect