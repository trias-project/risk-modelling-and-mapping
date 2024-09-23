results<-list()
# Use foreach and %dopar% to run the loop in parallel
gc()
n_cores <- detectCores()
cluster<-makeCluster(n_cores-1)
registerDoParallel(cluster) 
results <- foreach(i = seq_along(split_df), .packages = c("dplyr", "terra", "dismo", "sf", "raster", "rnaturalearth", "ggplot2", "mapview", "sdm", "kableExtra", "caret", "caretEnsemble", "here","tidyterra" )) %dopar% { 
  
  # Retrieve species name and data
  species <- names(split_df)[i]
  global.occ.LL.cleaned <- split_df[[i]]
  taxonkey <- unique(global.occ.LL.cleaned$acceptedTaxonKey)
  
  # (Your data processing and modeling code here)
  
  # Store results in a list for this species
 results[[species]]<-list(
    species = species,
    taxonkey = taxonkey
  )
}


