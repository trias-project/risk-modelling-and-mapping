#--------------------------------------------
#-------  Define wiSDM configurations -------
#--------------------------------------------
project <- "PA prob & Alternative Treshold & Ensemble Boyce"

species_to_model <- c("Verbena bonariensis", "Erigeron karvinskianus", "Leycesteria formosa", "Asclepias syriaca L.")

country_of_interest <- "Belgium"

occurrence_thinning_method <- "random" #either "random" or "kmeans_clustering"

n_clusters <- 10000 #Only relevant when occurrence_thinning_method is "kmeans_clustering"
