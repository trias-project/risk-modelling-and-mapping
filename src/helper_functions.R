
#-----------------------------------------------------------------------------------
#This function calculates the number of decimal places in any given numeric value 
# eg., 15.21 has 2 decimal places, 15.2569 has 4 decimal places, 15.25690 also has 4, as 0 in the end doesn't count
#-----------------------------------------------------------------------------------
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}


#-----------------------------------------------------------------------------------
#Divide a numerical value by 10
#-----------------------------------------------------------------------------------
divide10<-function(x){
  value<-x/10
  return(value)
}

#-----------------------------------------------------------------------------------
#Divide occurrence column with either y=0 (absences) or y=1 (presences)
#-----------------------------------------------------------------------------------
add.occ<-function(x,y){
  occ<-rep(y,nrow(x))
  cbind(x,occ)
}

#-----------------------------------------------------------------------------------
#Function to return threshold where sens=spec from caret results 
#-----------------------------------------------------------------------------------
findThresh<-function(df){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    mutate(observed= ifelse(obs == "present",1,0)) %>%
    select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::optimal.thresholds(df,opt.methods = 2)
  return(result)
}


#-----------------------------------------------------------------------------------
#Recalculate accuracy for a given model with the threshold that has been optimized
#-----------------------------------------------------------------------------------
accuracyStats<-function(df,y){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    mutate(observed= ifelse(obs == "present",1,0)) %>%
    select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::presence.absence.accuracy(df,threshold = y,st.dev=FALSE)
  return(result)
}


#-----------------------------------------------------------------------------------
# Model predictions for a large raster in a more efficient way using parallellization 
#-----------------------------------------------------------------------------------
predict_large_raster<-function(rasterstack, model, type) {
  
  # Ensure that connections are closed even in case of an error
  on.exit({
    plan(strategy = "sequential")  # Ensure that the parallel plan is returned to sequential
    gc()  # Trigger garbage collection
    closeAllConnections()  # Close any open file connections
  }, add = TRUE)
  
  gc() #Free up memory
  
  ncores<-min(4, availableCores()-2)  #Set up number of cores
  
  if(class(rasterstack)!="SpatRaster"){
    raster_terra<-rast(rasterstack)  #Convert raster to terra raster format if not already
  }else{
    raster_terra<-rasterstack
  }
  
  chunk_size <- ceiling(nrow(raster_terra) / ncores)   # Define chunk size
  
  # Create a list of row indices for each chunk
  chunk_indices <- split(seq_len(nrow(raster_terra)), ceiling(seq_along(seq_len(nrow(raster_terra))) / chunk_size))
  
  # Extract raster chunks and put raster chunks in list
  r_list<- vector(mode = "list", length = ncores)
  for (core in 1:ncores) {
    r_list[[core]]<- wrap(raster_terra[min(chunk_indices[[core]]):max(chunk_indices[[core]]), ,drop=FALSE])
  } #SpatRasters need to be wrapped before sending out to different cores
  
  # Save model to disk if it’s large
  saveRDS(model, "model.rds")
  options(future.globals.maxSize = 4.5 * 1024^3)
  plan(strategy = "multisession", workers=ncores) #Set up parallel
  
  out_list <- future_lapply(r_list,  function(chunk) {
    model <- readRDS("model.rds")  # Load model from disk
    unwrapped_raster <- unwrap(chunk)  # Unwrap raster for processing
    predicted_raster <- predict(unwrapped_raster, model, type = type, na.rm = TRUE)
    rm(unwrapped_raster)
    wrap(predicted_raster)  # Wrap the raster again
  }, future.seed = TRUE)
  
  
  plan(strategy = "sequential")   #Close parallel processing
  file.remove("model.rds")
  rm(r_list) #Remove large objects we don't need anymore
  out_list<- lapply(out_list, unwrap) #unwrap chunks
  gc() # Clean up memory after processing
  model_parallel<- do.call(terra::merge, out_list)  # Merge the chunks 
  rm(out_list) #Remove large objects we don't need anymore
  gc()  #Final garbage collect
  options(future.globals.maxSize = 500 * 1024^2)  # Reset to 500 MB
  return(model_parallel)
}


#-----------------------------------------------------------------------------------
# PDF export function
#-----------------------------------------------------------------------------------
exportPDF<-function(rst,taxonkey,taxonName,nameextension,is.diff="FALSE"){
  filename=file.path(pdfOutput,paste("be_",taxonkey, "_",nameextension,sep=""))
  pdf(file=filename,width=10,height=8,paper="a4r")
  par(bty="n")#to turn off box around plot
  ifelse(is.diff=="TRUE", brks<-seq(-1, 1, by=0.2), brks <- seq(0, 1, by=0.1)) 
  nb <- length(brks)-1 
  pal <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  cols<-pal(nb)
  maintitle<-paste(taxonName,taxonkey,"_",nameextension, sep= " ")
  plot(rst, breaks=brks, col=cols,main=maintitle, lab.breaks=brks,axes=FALSE)
  dev.off() 
} 


#-----------------------------------------------------------------------------------
# Export PNG function
#-----------------------------------------------------------------------------------
exportPNG<-function(rst,taxonkey,taxonName,nameextension,is.diff="FALSE"){
  filename=file.path(pdfOutput,paste("be_",taxonkey, "_",nameextension,sep=""))
  png(file=filename)
  par(bty="n")#to turn off box around plot
  ifelse(is.diff=="TRUE", brks<-seq(-1, 1, by=0.25), brks <- seq(0, 1, by=0.1)) 
  nb <- length(brks)-1 
  pal <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  cols<-pal(nb)
  maintitle<-paste(taxonName,taxonkey,"_",nameextension, sep= " ")
  plot(rst, breaks=brks, col=cols,main=maintitle, lab.breaks=brks,axes=FALSE)
  dev.off() 
} 


#-----------------------------------------------------------------------------------
# Generate pseudoabsences
#-----------------------------------------------------------------------------------
generate_pseudoabs <- function(index = NULL,mask, alternative_mask, n, p) {
  tryf_values <- c(50,100, 150)  # tryf values to attempt in each stage
  current_raster <- mask  # Start with the initial raster layer
  
  # Attempt to generate points
  for (tryf in tryf_values) {
    # Generate random points
    suppressWarnings(pseudoabs <- as.data.frame(
      randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = FALSE, 
        cellnumbers = FALSE, 
        tryf = tryf, 
        warn = 2, 
        lonlatCorrection = TRUE
      )
    )
    )
    # Check if the number of pseudoabsences reaches required amount
    if (nrow(pseudoabs) == n) {
      # If index is provided, include it in the message (only for lists)
      if (!is.null(index)) {
        message(paste0(n, " out of ", n, " pseudoabsences generated while accounting for observer bias in set ", index))
      } else {
        message(paste0(n, " out of ", n, " pseudoabsences generated while accounting for observer bias."))
      }
      return(pseudoabs)  # Return dataset if the required amount of pseudoabsences are generated
    }
  }
  
  # If unsuccessful with biasgrid ecoregions raster, switch to the full ecoregions raster and retry
  current_raster <- alternative_mask
  
  for (tryf in tryf_values) {
    pseudoabs <- as.data.frame(
      randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = FALSE, 
        cellnumbers = FALSE, 
        tryf = tryf, 
        warn = 2, 
        lonlatCorrection = TRUE
      )
    )
    
    # Check if the number of rows meets the desired count
    if (nrow(pseudoabs) == n) {
      # If index is provided, include it in the warning (only for lists)
      if (!is.null(index)) {
        warning(paste0(n, " out of ", n, " pseudoabsences generated without accounting for observer bias in set ", index))
      } else {
        warning(paste0(n, " out of ", n, " pseudoabsences generated without accounting for observer bias."))
      }
      return(pseudoabs)  # Return dataset if enough pseudoabsences were generated
    }
  }
  
  # If all attempts fail, return the last generated dataframe with fewer pseudoabsences than requested
  # If index is provided, include it in the warning
  if (!is.null(index)) {
    warning(paste0("Could not generate the required number of pseudoabsences: ", n, " out of ", n, " pseudoabsences generated without accounting for observer bias in set ", index))
  } else {
    warning(paste0("Could not generate the required number of pseudoabsences: ", n, " out of ", n, " pseudoabsences generated without accounting for observer bias."))
  }
  
  return(pseudoabs)  # Return the pseudoabs data, even if incomplete
}


#-----------------------------------------------------------------------------------
