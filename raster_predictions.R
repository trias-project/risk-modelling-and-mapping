#Option 1
#Original: 268.02 
system.time({
  global_model<-raster::predict(eu_climpreds.10,global_stack,type="prob")
})

#Option 2: 326.14, 2nd time: 479.00 

eu_terra<-rast(eu_climpreds.10)
system.time({
  global_model_terra<-terra::predict(eu_terra,global_stack,type="prob", cores=2, na.rm=TRUE)
})


#Option 3: SO FAR ALWAYS THE FASTEST BUT NOT BY A LOT (30sec-10sec)
gc()
library(terra)
library(future)
library(future.apply)
library(progressr)
handlers(global = TRUE)  # Set up global progress handler
system.time({
ncores<-floor(parallel::detectCores()/2)
eu_terra<-rast(eu_climpreds.10)
# Define chunk size
chunk_size <- ceiling(nrow(eu_terra) / ncores)

# Create a list of row indices for each chunk
chunk_indices <- split(seq_len(nrow(eu_terra)), ceiling(seq_along(seq_len(nrow(eu_terra))) / chunk_size))

# Extract raster chunks
# Put raster chunks in list
r_list<- vector(mode = "list", length = ncores)
for (core in 1:ncores) {
  r_list[[core]]<- wrap(eu_terra[min(chunk_indices[[core]]):max(chunk_indices[[core]]), ,drop=FALSE])
} #SpatRasters need to be wrapped before sending out to different cores


predict_parallel <- function(chunk_raster, model, ...) {
  # Unwrap raster for processing
  unwrapped_raster <- unwrap(chunk_raster)
  
  # Perform prediction
  predicted_raster<- predict( unwrapped_raster, model, type = "prob", na.rm=TRUE)
  
  #Wrap the raster again
  wrap(predicted_raster)
}

#In Parallel
plan(strategy = "multisession", workers=ncores) #Set up parallel

with_progress({
  p <- progressor(along = r_list) # Set up a progressor object for the r_list
out_list<- future_lapply(r_list, FUN = function(chunk) {
  on.exit(p(), add = TRUE)  # Ensure progress is updated even if an error occurs and update progress for each chunk
  predict_parallel(chunk, model = global_stack)})
})

plan(strategy = "sequential")
rm(r_list)
out_list<- lapply(out_list, unwrap) #unwrap chunks
gc() # Clean up memory after processing
# Merge the chunks 
global_model_parallel<- do.call(terra::merge, out_list)
rm(out_list)
gc()
})


#Option 4: 271.78 
gc()
system.time({
  global_model2<-raster::predict(eu_climpreds.10,global_stack,type="prob", cores=2)
})
