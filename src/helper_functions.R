
#-----------------------------------------------------------------------------------
#This function calculates the number of decimal places in any given numeric value 
# eg., 15.21 has 2 decimal places, 15.2569 has 4 decimal places, 15.25690 also has 4, as 0 in the end doesn't count
#-----------------------------------------------------------------------------------
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    # Remove trailing zeros and split at the decimal point
    split_result <- strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]]
    # Check if there are any decimals
    if (length(split_result) > 1) {
      nchar(split_result[[2]]) # Count characters in the decimal part
    } else {
      return(0) # No decimal part
    }
  } else {
    return(0) # No decimals for whole numbers
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
# Only keep one occurrence point per grid cell of a spatRaster
#-----------------------------------------------------------------------------------
remove_duplicates <- function(occurrences, rast_template){
  
  #Initial dataset
  initial_occurrences<-nrow(occurrences)
  
  #Indicate for each occurrence point in which cell of the raster it falls
  occurrences$cell <- terra::cellFromXY( rast_template, occurrences) 
  
  #Remove occurrences that don't fall in any cell of the raster and duplicate occurrences
  occurrences <- occurrences %>%
    dplyr::filter(!is.na(cell)) %>% 
    dplyr::distinct(cell, .keep_all=TRUE) %>% 
    dplyr::select(1:2)
  
  #Print how many occurrences were removed
  print(paste(initial_occurrences - nrow(occurrences), "duplicate occurrence records removed."))
  
  return(occurrences)
  
}


#-----------------------------------------------------------------------------------
# Remove occurrences that fall in grid cells with NA values
#-----------------------------------------------------------------------------------
remove_nodata_occurrences <- function(occurrences, rast_template, crs){
  
  #Store number of initial occurrences
  initial_occurrences<-nrow(occurrences)
  
  #Remove occurrences in NA cells and convert to sf
  env<- terra::extract(rast_template, occurrences, xy = F, ID = F)
  
  occurrences <- cbind(env, occurrences)%>%
    dplyr::filter(!is.na(.[[1]])) %>%   # keep only rows where first column is not NA
    dplyr::select(c(decimalLongitude, decimalLatitude)) %>%
    sf::st_as_sf(coords = c("decimalLongitude", "decimalLatitude"), crs = crs, remove = FALSE)
  
  #Print how many occurrences were removed
  print(paste(initial_occurrences - nrow(occurrences), "occurrence records in grid cells with NAs removed."))
  
  return(occurrences)
}


#-----------------------------------------------------------------------------------
#Divide occurrence column with either y=0 (absences) or y=1 (presences)
#-----------------------------------------------------------------------------------
add.occ<-function(x,y){
  occ<-rep(y,nrow(x))
  cbind(x,occ)
}


#-----------------------------------------------------------------------------------
#Define the favourability transformation function
#-----------------------------------------------------------------------------------
favourability_from_prob <- function(prob_raster, prev_ratio) {
  odds <- prob_raster / (1 - prob_raster)
  fav <- odds / (prev_ratio + odds)
  fav[is.infinite(fav)] <- NA
  fav[fav < 0] <- 0
  fav[fav > 1] <- 1
  return(fav)
}


#-----------------------------------------------------------------------------------
#Function to return threshold where sens=spec from caret results 
#-----------------------------------------------------------------------------------
findThresh<-function(df){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    dplyr::mutate(observed= ifelse(obs == "present",1,0)) %>%
    dplyr::select(rowIndex,observed,predicted=present)
  result<-PresenceAbsence::optimal.thresholds(df,opt.methods = 2)
  return(result)
}


#-----------------------------------------------------------------------------------
#Recalculate accuracy for a given model with the threshold that has been optimized
#-----------------------------------------------------------------------------------
accuracyStats<-function(df,y){
  df<-df[c("rowIndex","obs","present")]
  df<-df %>%
    dplyr::mutate(observed= ifelse(obs == "present",1,0)) %>%
    dplyr::select(rowIndex,observed,predicted=present)
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
    raster_terra<-terra::rast(rasterstack)  #Convert raster to terra raster format if not already
  }else{
    raster_terra<-rasterstack
  }
  
  chunk_size <- ceiling(nrow(raster_terra) / ncores)   # Define chunk size
  
  # Create a list of row indices for each chunk
  chunk_indices <- split(seq_len(nrow(raster_terra)), ceiling(seq_along(seq_len(nrow(raster_terra))) / chunk_size))
  
  # Extract raster chunks and put raster chunks in list
  r_list<- vector(mode = "list", length = ncores)
  for (core in 1:ncores) {
    r_list[[core]]<- terra::wrap(raster_terra[min(chunk_indices[[core]]):max(chunk_indices[[core]]), ,drop=FALSE])
  } #SpatRasters need to be wrapped before sending out to different cores
  
  # Save model to disk if it’s large
  saveRDS(model, "model.rds")
  options(future.globals.maxSize = 4.5 * 1024^3)
  plan(strategy = "multisession", workers=ncores) #Set up parallel
  
  out_list <- future_lapply(r_list,  function(chunk) {
    model <- readRDS("model.rds")  # Load model from disk
    unwrapped_raster <- terra::unwrap(chunk)  # Unwrap raster for processing
    predicted_raster <- terra::predict(unwrapped_raster, model, type = type, na.rm = TRUE)
    rm(unwrapped_raster)
    terra::wrap(predicted_raster)  # Wrap the raster again
  }, future.seed = TRUE)
  
  
  plan(strategy = "sequential")   #Close parallel processing
  file.remove("model.rds")
  rm(r_list) #Remove large objects we don't need anymore
  out_list<- lapply(out_list, terra::unwrap) #unwrap chunks
  gc() # Clean up memory after processing
  model_parallel<- do.call(terra::merge, out_list)  # Merge the chunks 
  rm(out_list) #Remove large objects we don't need anymore
  gc()  #Final garbage collect
  options(future.globals.maxSize = 500 * 1024^2)  # Reset to 500 MB
  return(model_parallel)
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
  pal <- grDevices::colorRampPalette(rev(brewer.pal(11, 'Spectral')))
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
      dismo::randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = TRUE, 
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
      dismo::randomPoints(
        current_raster, 
        n, 
        p, 
        ext = NULL, 
        extf = 1.1, 
        excludep = TRUE, 
        prob = TRUE, 
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
# Recode factor levels to absent (0) and present(1), and set present as the reference level
#-----------------------------------------------------------------------------------
factorVars<-function(df,var){
  df[,c(var)]<-as.factor(df[,c(var)])
  levels(df[,c(var)])<-c("absent","present")
  df[,c(var)]<-relevel(df[,c(var)], ref = "present")
  return(df)
}


#-----------------------------------------------------------------------------------
#----------------Create folders when they don't exist yet---------------------------
#-----------------------------------------------------------------------------------
create_folder <- function(path, name) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message(paste0("Folder '", name, "' created at path: '", path, "' 🎉"))
  }
}


#-----------------------------------------------------------------------------------
#-------------------------Export PDF new function-----------------------------------
#-----------------------------------------------------------------------------------
exportPDF <- function(predictions=NULL, period=NULL, scenario, occ_data=NULL, dataType, returnPredictions=FALSE,returnPNG=FALSE, providedPNG=NULL, exportPNG=FALSE, LabelValue=NULL, LabelName=NULL, Label2Value=NULL, Label2Name=NULL, PDF_title, PNG_folder=NULL, PDF_folder, filename){
  
  #Set scenario to "" if period is Current
  if(period=="Current") scenario<-""
  
  #Define scenario title
  scenarioTitle<- switch(paste0(period,scenario),
                         "Current" = "Current",
                         "all" = "all",
                         "2041-2070ssp126" = "2041-2070: SSP1-2.6",
                         "2041-2070ssp370" = "2041-2070: SSP3-7.0",
                         "2041-2070ssp585" = "2041-2070: SSP5-8.5",
                         "2071-2100ssp126" = "2071-2100: SSP1-2.6",
                         "2071-2100ssp370" = "2071-2100: SSP3-7.0",
                         "2071-2100ssp585" = "2071-2100: SSP5-8.5"
  )
  
  PNG_filename <- paste0(filename, ".png")
  PDF_filename <- paste0(filename, ".pdf")
  
  # Define file paths
  if(!is.null(PNG_folder)){
    plot_png_path <- file.path(PNG_folder, PNG_filename)# If not EU or Global
  }
  plot_pdf_path <- file.path(PDF_folder, PDF_filename)
  
  
  #If png is not provided, create a PNG based on the input predictions
  if(is.null(providedPNG)){
    
    #Get extent
    exten<-as.vector(terra::ext(predictions))
    
    #Settings for plot
    if (dataType == "Diff") {
      brks <- seq(-1, 1, by = 0.25)
    } else if (dataType %in% c("Suit", "Conf", "Masked_Suit", "Stdev")) {
      brks <- seq(0, 1, by = 0.2)
    }
    
    if (dataType != "Binary"){
      nb <- length(brks) - 1
      viridis_palette <- viridis::viridis(nb)
      
      #Create dummy raster with all values
      template <- predictions
      values(template) <- rep(brks, length.out = ncell(template))
      template <- mask(template, predictions)
      
      #Create plot
      suppressMessages(
        country_plot <- ggplot() + 
          geom_spatraster(data=template)+
          geom_spatraster(data = predictions) +
          scale_fill_gradientn(colors = viridis_palette, 
                               breaks = brks, 
                               labels = brks, 
                               na.value = "transparent") +
          theme_bw() +
          theme(axis.title = element_blank())+
          theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
          coord_sf(xlim = c(exten[1] - (exten[1] * 0.12), exten[2]- (exten[2] * 0.05)), 
                   ylim = c(exten[3], exten[4] + (exten[4] * 0.04)))
      )
      
    }else{
      #Create plot
      suppressMessages(
        country_plot <- ggplot() + 
          geom_spatraster(data = predictions) +
          scale_fill_manual(values = c("Absent" = "lightgrey", "Present" = "#085099"),
                            na.value = "transparent",
                            na.translate=FALSE)+
          theme_bw() +
          theme(axis.title = element_blank())+
          theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
          coord_sf(xlim = c(exten[1] - (exten[1] * 0.12), exten[2]- (exten[2] * 0.05)), 
                   ylim = c(exten[3], exten[4] + (exten[4] * 0.04)))
      )
    }
    # Define text label, fill label, and hjust based on dataType
    text_label <- ifelse(dataType=="Diff",paste(scenarioTitle, "- current"), scenarioTitle)
    
    fill_label <- switch(dataType,
                         "Suit" = "Suitability",
                         "Diff" = "Suitability difference",
                         "Conf" = "Confidence",
                         "Masked_Suit" = "Suitability",
                         "Binary" = "Suitability",
                         "Stdev" = "Standard deviation")
    
    hjust_value <- ifelse(dataType == "Diff", -0.264 , -0.24) 
    x_value<-ifelse(exten[1]>180, 800000, -35)
    # Update the plot
    country_plot <- country_plot +
      labs(fill = fill_label) +
      annotate("text",
               x = x_value, y = Inf,       # Position at top-right
               label = text_label,     # Text to display
               hjust = 0,
               vjust = 2.4,            # Adjust text alignment to the right and above
               size = 4.8,
               color = "#636363",
               fontface = "bold")+
      theme(aspect.ratio=NULL)
    
    if(!is.null(occ_data)){
      crs_value<-st_crs(occ_data)
      #Only show occurrences that fall within raster cells
      occ_data<- terra::extract(predictions, occ_data, xy = T, ID=F)%>%
        dplyr::filter(!is.na(.[, 1]))%>% #Keep rows that do not have any NA values in value column
        dplyr::select(c(x,y))%>%
        dplyr::rename(decimalLongitude=x,
                      decimalLatitude=y)%>%
        st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), crs=crs_value)
      
      suppressMessages(
        country_plot<-country_plot +
          geom_sf(data = occ_data, color = "black", fill = "red", 
                  size = 1.5, shape = 21)+
          coord_sf(xlim = c(exten[1] - (exten[1] * 0.12), exten[2]- (exten[2] * 0.05)), 
                   ylim = c(exten[3], exten[4] + (exten[4] * 0.04)))
      )
    }
    
    if(!is.null(LabelValue)){
      
      assertthat::assert_that(!is.null(LabelName), msg = "LabelValue is provided but LabelName is not.")
      
      country_plot<-country_plot +
        annotate("text",
                 x = x_value, y = Inf,       # Position at top-right
                 label = paste(LabelName,"=",LabelValue),     # Text to display
                 hjust = 0,
                 vjust = 4.5,            # Adjust text alignment to the right and above
                 size = 4.4,
                 color = "#919191",
                 fontface = "bold")
    }
    
    if(!is.null(Label2Value)){
      
      assertthat::assert_that(!is.null(Label2Name), msg = "Label2Value is provided but Label2Name is not.")
      
      country_plot<-country_plot +
        annotate("text",
                 x = x_value, y = Inf,       # Position at top-right
                 label = paste(Label2Name,"=",Label2Value),     # Text to display
                 hjust = 0,
                 vjust = 6.5,            # Adjust text alignment to the right and above
                 size = 4.4,
                 color = "#919191",
                 fontface = "bold")
    }
    
    #Create an empty plot to fill PDF
    empty_plot <- ggplot() + 
      theme_void() + 
      theme(plot.background = element_blank()) 
    
    #Create final plot
    plot_final<-country_plot 
    
    # Save plot temporarily as a PNG file
    ggplot2::ggsave(filename = PNG_filename, plot = plot_final, 
                    device = "png", width =7.7 , height = 6.94, units = "in", dpi= 300, path= PDF_folder)
    
  }else{
    # Save plot temporarily as a PNG file
    ggplot2::ggsave(filename = PNG_filename, plot =providedPNG, 
                    device = "png", width =7.7 , height = 6.94, units = "in", dpi= 300, path= PDF_folder)
  }
  
  # Read the PNG image back in
  img <- magick::image_read(here::here(PDF_folder, PNG_filename))
  
  # Start a PDF device for output (A4 portrait in inches)
  pdf(plot_pdf_path, width = 8.27, height = 11.69)
  
  grid.newpage()
  
  # Add title at the top of the PDF
  grid.text(
    label = PDF_title,
    x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 14, fontface = "bold")
  )
  
  # Insert the PNG centered on the page
  # - x = 0.5 centers horizontally
  # - y = 0.5 centers vertically
  # - just = "center" anchors the image by its center point
  # - width/height use npc units so the image scales to the page
  grid::grid.raster(
    img,
    x = unit(0.5, "npc"),
    y = unit(0.627, "npc"),  # Put plot under title
    width = unit(0.95, "npc"),
    height = unit(0.6, "npc"),  # roughly matches PNG aspect ratio
    just = "center"
  )
  
  # Close the PDF device (single close)
  dev.off()
  
  # Print confirmation
  print(paste(PDF_filename," has been created.", sep=""))
  
  # Remove the temporary PNG file
  file.remove(here::here(PDF_folder, PNG_filename))
  gc()
  
  #Store PNG file in PNG folder if exportPNG is TRUE
  if(exportPNG){
    ggplot2::ggsave(filename = PNG_filename, plot = country_plot, 
                    device = "png", width = 7.7, height = 6.94, units = "in", dpi = 300,
                    path = PNG_folder)
    print(paste(PNG_filename," has been created.", sep="")) 
  }
  
  
  #Store plots or models
  if (returnPredictions & returnPNG) {
    
    return(list("model" = predictions,
                "png"=country_plot,
                "scenario"=scenario))
  }
  
  if (returnPredictions & !returnPNG) {
    # return(setNames(list("model" = predictions), scenario))
    return(list("model" = predictions,
                "scenario"=scenario))
  }
  
  if (!returnPredictions & returnPNG) {
    # return(setNames(list("model" = predictions), scenario))
    return(list("png" = country_plot,
                "scenario"=scenario))
  }
}


#-----------------------------------------------------------------------------------
#------------------------- Standardize residuals -----------------------------------
#-----------------------------------------------------------------------------------
stdres<-function(obs.numeric, yhat){
  num<-obs.numeric-yhat #Obtain residuals (the difference between observed and predicted values)
  denom<-sqrt(yhat * (1 - yhat) + 1e-10) #Approximates the residual variance in logistic regression,  + 1e-10 is added in case predicted values are 0
  return(num/denom)#Standardize
}


#-----------------------------------------------------------------------------------
#--------Return number of elements that are equal or less than a threshold ---------
#-----------------------------------------------------------------------------------

GetLength <- function(x, y) {
  sum(x <= y)
}

#-----------------------------------------------------------------------------------
#--------- classify based on probabilities compared to a confidence level ----------
#-----------------------------------------------------------------------------------
CPconf<-function(pA,pB,confidence){
  if(pA > confidence && pB< confidence){
    predClass<-"classA"
  }else if(pA < confidence && pB> confidence){
    predClass<-"classB"
  }else if(pA< confidence && pB< confidence){
    predClass<-"noClass"
  }else{
    predClass<-"bothClasses"
  }
  return(predClass)
}


#-----------------------------------------------------------------------------------
#---------------------- calculate confidence of each prediction --------------------
#-----------------------------------------------------------------------------------
get.confidence<-function(pvalA,pvalB){
  secondHighest<-ifelse(pvalA>pvalB,pvalB,pvalA)
  conf<-(1-secondHighest)
  return(conf)
}



#-----------------------------------------------------------------------------------
#-------------- Return presence/absence based on values a and b --------------------
#-----------------------------------------------------------------------------------
forcedCp<-function(pvalA,pvalB){
  ifelse(pvalA>pvalB,"presence","absence")
}


#-----------------------------------------------------------------------------------
#--- Extract probability of presence and absence from prediction raster ------------
#-----------------------------------------------------------------------------------
extractVals<-function(predras){
  vals <-  as.numeric(terra::values(predras))
  vals[is.nan(vals)] <- NA
  coord <-  terra::xyFromCell(predras,1:terra::ncell(predras))
  raster_fitted <- cbind(coord,vals)
  raster_fitted.df<-as.data.frame(raster_fitted)
  raster_fitted.df1<-na.omit(raster_fitted.df)
  raster_fitted.df1$presence<-raster_fitted.df1$vals
  raster_fitted.df1$absence<- (1-raster_fitted.df1$presence)
  return(raster_fitted.df1)
}


#-----------------------------------------------------------------------------------
#-------------- Class conformal prediction --------------------
#-----------------------------------------------------------------------------------
classConformalPrediction<-function(x,y){
  #Extract model results
  ens_results <- x
  ens_calib<-ens_results$ens_model$pred
  
  # Filter and extract calibration data for presence and absence
  calibPresence<-ens_calib %>%
    dplyr::filter(obs=='present')%>%
    dplyr::select(present)
  calibPresence<-unname(unlist(calibPresence[c("present")]))
  
  calibAbsence<-ens_calib %>%
    dplyr::filter(obs=='absent')%>%
    dplyr::select(absent)
  calibAbsence<-unname(unlist(calibAbsence[c("absent")]))
  
  #Extract predicted values
  predicted.values<-extractVals(y)
  testPresence<-predicted.values$presence
  testAbsence<-predicted.values$absence
  
  #derive p.Values for class A
  smallrA<-lapply(testPresence,function(x) GetLength(calibPresence,x))#For each value in testPresence, you calculate the number of values in calibPresence that are smaller or equal to the testPresence value
  smallrA_1<- unlist (smallrA)+1 #Create a vector of resulting values and add 1 to each value
  nCalibSet<-length(calibPresence)+1
  pvalA<-smallrA_1+1/nCalibSet
  
  # derive p.Values for Class B
  smallrB<-lapply(testAbsence,function(x) GetLength(calibAbsence,x))
  smallrB_1<- unlist (smallrB)+1
  nCalibSetB<-length(calibAbsence)
  pvalB<-smallrB_1/nCalibSetB
  
  pvalsdf<-as.data.frame(cbind(pvalA,pvalB,0.20))
  #raster_cp_20<-mapply(CPconf,pvalsdf$pvalA,pvalsdf$pvalB,pvalsdf[3])
  #table(raster_cp_20)
  
  pvalsdf$conf<-get.confidence(pvalsdf$pvalA,pvalsdf$pvalB)
  pvalsdf_1<-cbind(pvalsdf,predicted.values)
  
  return(pvalsdf_1)  
}


#-----------------------------------------------------------------------------------
#--------------------------- Create confidence maps --------------------------------
#-----------------------------------------------------------------------------------
confidenceMaps<-function(x,original_raster,taxonName, taxonNameTitle, nameExtension, taxonKey ,scenario, regionName, scenarioTitle, dataType, folder, GlobalModel=FALSE, resampling_rast=NULL, country_sf=NULL){
  # Create a SpatVector from the data.xyz
  data.xyz <- x[c("x","y","conf")]
  points <- terra::vect(data.xyz, geom = c("x", "y"), crs = terra::crs(original_raster))
  
  # Rasterize the points using the original SpatRaster as a template
  rst <- terra::rasterize(points, original_raster, field = "conf")
  
  #If global model is used, resample map
  if(GlobalModel){
    rst_to_export<- rst%>%
      terra::project(terra::crs(country_sf))%>%
      terra::resample(resampling_rast, method="bilinear") 
  }else{
    
    rst_to_export<-rst
  }
  
  #Export raster
  raster_file<-paste(taxonName, "_", taxonKey, "_", scenario, "_confidence_", regionName, ".tif", sep="")
  terra::writeRaster(rst_to_export,
                     filename=file.path(folder, raster_file),
                     overwrite=TRUE)
  #Print
  print(paste(raster_file," has been created.", sep=""))
  
  exportPDF(predictions=rst_to_export,
            taxonName,
            nameExtension, 
            taxonNameTitle,
            taxonKey=taxonKey, 
            scenario, 
            regionName,
            returnPredictions=FALSE,
            returnPNG=FALSE, 
            dataType="Conf")
  return(rst)
  
}


#----------------------------------------------------------
#---------------- Assess response curves-------------------
#----------------------------------------------------------
#evaluate predictions while varying only the selected variable (x) and keeping all other variables at their observed values
partial_gbm<-function(x){
  m.gbm<-pdp::partial(bestModel$models$gbm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,n.trees= bestModel$models$gbm$finalModel$n.trees, which.class = 1,grid.resolution=nrow(bestModel.train))
}

partial_glm<-function(x){
  m.glm<-pdp::partial(bestModel$models$glm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

partial_mars<-function(x){
  m.mars<-pdp::partial(bestModel$models$earth$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                       prob=TRUE,which.class = 2,grid.resolution=nrow(bestModel.train)) # class=2 because in earth pkg, absense is the first class
}

partial_rf<-function(x){
  pdp::partial(bestModel$models$rf$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
               prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

#----------------------------------------------------------
#---------------- Plot response curves-------------------
#----------------------------------------------------------

responseCurves<-function(x,y) {
  colors <- c("GLM" = "gray", "GBM"="red","RF"="blueviolet","MARS"= "hotpink") 
  ggplot(all_dfs,(aes(x=.data[[x]],y=.data[[y]]))) +
    geom_line(aes(color = data), size =1.2, position=position_dodge(width=0.2))+
    theme_bw()+
    labs(y="Partial probability", x= gsub("//..*","",x),color="Legend") +
    scale_color_manual(values = colors)
}  



#----------------------------------------------------------
#------------- Evaluate model predictions------------------
#----------------------------------------------------------
eu_eval<-function (ras,y){
  indep.bil<-terra::extract(ras,y,method="bilinear")
  indep.bil.df<-as.data.frame(indep.bil)
  indep.bil.df<-indep.bil.df %>%
    dplyr::mutate(predicted= ifelse(indep.bil >= 0.5,"present","absent")) 
  indep.bil.df$observed<-rep("present",nrow(indep.bil.df))
  indep.bil.df$predicted<-as.factor(indep.bil.df$predicted)
  indep.bil.df$observed<-as.factor(indep.bil.df$observed)
  xtab<-table(indep.bil.df$predicted,indep.bil.df$observed)
  return(xtab)
}

#----------------------------------------------------------
#-------------   update_files_logic   ----------------------
#----------------------------------------------------------
#' @param dest_file output file to test existance of
#' @param update_files whether to ask, or ignore existance of input environmental layers
update_files_logic <- function(dest_file,
                               dest_folder = NULL,
                               update_files) {
  
  # normalize update_files
  if (is.factor(update_files)) update_files <- as.character(update_files)
  update_files <- trimws(tolower(update_files))  # remove spaces, make lowercase
  
  if (!update_files %in% c("yes","no","ask")) {
    stop("update_files must be 'yes', 'no', or 'ask'")
  }
  
  # ---------- SINGLE FILE ----------
  if (!is.data.frame(dest_file)) {
    
    if (update_files == "yes") {
      
      update_files_final <- TRUE
      
    } else if (update_files == "no") {
      
      update_files_final <- !file.exists(dest_file)
      
    } else {  # "ask"
      
      if (file.exists(dest_file)) {
        update_files_final <- askYesNo(
          paste("Download\n", basename(dest_file), "\n again?")
        )
      } else {
        update_files_final <- TRUE
      }
    }
    
    # ---------- MULTIPLE FILES ----------
  } else {
    
    if (update_files == "yes") {
      
      dest_file$update_file <- TRUE
      
    } else if (update_files == "no") {
      
      dest_file$update_file <- !file.exists(
        file.path(dest_folder, dest_file$file)
      )
      
    } else {  # "ask"
      
      dest_file$update_file <- FALSE
      
      for (i in seq_len(nrow(dest_file))) {
        f <- file.path(dest_folder, dest_file$file[i])
        
        if (!file.exists(f)) {
          dest_file$update_file[i] <- TRUE
        } else {
          dest_file$update_file[i] <- askYesNo(
            paste("Download\n", basename(f), "\n again?")
          )
        }
      }
    }
    
    update_files_final <- dplyr::filter(dest_file, update_file)
  }
  
  return(update_files_final)
}

#----------------------------------------------------------
#-------------   safe_download_zenodo   -------------------
#----------------------------------------------------------

#' a wrapper for zen4R::download_zenodo to delete the failed file and thus trigger 
#' a redownload with update_files == FALSE

safe_download_zenodo <- function(doi, path, files, timeout = 600, quiet = FALSE) {
  tryCatch(
    {
      zen4R::download_zenodo(
        doi    = doi,
        path   = path,
        files  = files,
        timeout = timeout,
        quiet  = quiet
      )
    },
    error = function(e) {
      msg <- conditionMessage(e)
      
      # Try to extract the dest_file name from the error message
      # Adapt the regex to match the actual message format
      m <- regexpr("dest_file ['\"]?([^'\" ]+)['\"]?", msg)
      if (m[1] != -1) {
        fname <- regmatches(msg, m)
        # fname now contains something like \"dest_file 'myfile.ext'\"
        # extract just the filename part
        fname_only <- sub(".*dest_file ['\"]?([^'\" ]+)['\"]?.*", "\\1", fname)
        
        file_to_remove <- file.path(path, fname_only)
        if (file.exists(file_to_remove)) {
          unlink(file_to_remove)
        }
      }
      
      # Re-display the original error
      stop(e)
    }
  )
}


#----------------------------------------------------------
#---- Check, and if necessary, redownloaded tif files -----
#----------------------------------------------------------
read_or_redownload <- function(file, folder, doi, max_attempts = 3) {
  
  file_path <- file.path(folder, file)
  attempt <- 1
  
  while (attempt <= max_attempts) {
    
    r <- tryCatch({
      r <- terra::rast(file_path)
      terra::ncell(r)
      
    }, error = function(e) {
      NULL
    })
    
    if (!is.null(r)) {
      return(r)  # success
    }
    
    message(paste("Corrupt file detected:", file, "- redownloading (attempt", attempt, ")"))
    
    # Remove corrupt file if it exists
    if (file.exists(file_path)) {
      file.remove(file_path)
    }
    
    # Redownload
    zen4R::download_zenodo(
      doi = doi,
      path = folder,
      files = file,
      timeout = 600,
      quiet = FALSE
    )
    
    attempt <- attempt + 1
  }
  
  stop(paste("Failed to obtain valid raster after", max_attempts, "attempts:", file))
}



#----------------------------------------------------------------------
#- Make predictions per model algorithm and dataset and obtain median -
#----------------------------------------------------------------------

compute_median_favourability <- function(model,
                                         datasets,
                                         top5_methods,
                                         prev_ratio) {
  
  #---------------------------
  #---- Make predictions -----
  #---------------------------
  
  env_favourability <- list()
  for(modelmethod in top5_methods){
    
    message("Predicting for method: ", modelmethod,".")
    
    for(dataset_name in names(datasets)) {
      
      #Load datasets
      dataset <- datasets[[dataset_name]]
      IDs <-dataset$ID
      dataset<-dplyr::select(dataset, -ID)
      
      #Predict for dataset
      dataset_suit <- predict(model,
                              newdata = dataset,
                              method = modelmethod)
      
      #Convert suitability to favourability
      dataset_fav<- favourability_from_prob(dataset_suit[[1]], prev_ratio)
      
      #Store in list
      env_favourability[[modelmethod]][[dataset_name]] <- data.frame(ID = IDs,
                                                                     fav = dataset_fav)
      
      #Clean up
      rm(dataset_suit, dataset_fav, IDs, dataset)
      
    }
  }  
  
  
  #-----------------------------------------
  #---- Calculate median favourability  ----
  #-----------------------------------------
  median_favourability<-lapply(
    names(datasets),
    function(dataset_name) {
      
      fav_matrix <- do.call(
        cbind,
        lapply(env_favourability, function(x) x[[dataset_name]]$fav)
      )
      
      data.frame(ID = env_favourability[[1]][[dataset_name]]$ID,
                 median_favourability = matrixStats::rowMedians(fav_matrix,na.rm = TRUE))
    }
  )
  
  names(median_favourability) <- names(datasets)
  
  return(median_favourability)
}


#------------------------------------------
#----- Define Boyce helper functions -----
#------------------------------------------
#fit_vals are suitability or favourability values of a background sample or all available pixels in the target region
#obs_vals are suitability or favourability values at occurrence locations

compute_boyce_robust <- function(fit_vals, obs_vals) {
  
  #------------------------------------------
  #----------- Basic checks -----------------
  #------------------------------------------
  #Define conditions for number of occurrences and total (or background) points
  if (length(obs_vals) < 5 || length(fit_vals) < 200) return(NA_real_)
  if (length(unique(fit_vals)) < 3)                   return(NA_real_)
  
  
  #------------------------------------------
  #---- Try moving window boyce ---
  #-----------------------------------------
  boyce_result <- try(ecospat::ecospat.boyce(fit = fit_vals, 
                                             obs = obs_vals,
                                             nclass = 0, #moving window boyce index
                                             PEplot = FALSE), 
                      silent = TRUE)
  
  
  #------------------------------------------
  #--- Return boyce index if all went well -----
  #------------------------------------------
  if (!inherits(boyce_result, "try-error") && !is.null(boyce_result$cor) && is.finite(boyce_result$cor)){
    return(boyce_result$cor)
    
  }else{
    
    #------------------------------------------
    #------------ Try binned Boyce ------------
    #------------------------------------------
    for (nc in c(10L, 20L)) {
      boyce_result2 <- try(
        ecospat::ecospat.boyce(
          fit    = fit_vals,
          obs    = obs_vals,
          nclass = nc,
          PEplot = FALSE
        ),
        silent = TRUE
      )
      
      if (!inherits(boyce_result2, "try-error") && !is.null(boyce_result2$cor) && is.finite(boyce_result2$cor)){
        return(boyce_result2$cor)
      }
    }
  }  
  #------------------------------------------
  #--Return NA if binned boyce also fails ---
  #------------------------------------------
  return(NA_real_)
  
}


#------------------------------------------
#---Calculate model validation metrics ----
#------------------------------------------
compute_validation_metrics <- function(species, type, region, fold, all_suit_vals, occ_suit_vals, abs_suit_vals) {
  
  #Define number of presences and pseudoabsences
  n_pres   <- length(occ_suit_vals)
  n_abs    <- length(abs_suit_vals)
  
  #Define minimum number of presences and pseudoabsences needed
  if (n_pres < 2L || n_abs < 2L)   return(c(fold = fold,
                                            auc = NA_real_,
                                            boyce = NA_real_,
                                            tss = NA_real_))
  
  #---------------
  #----- AUC -----
  #---------------
  # AUC via Wilcoxon statistic (exact, no ties correction needed for our use)
  auc_val <- tryCatch({
    wt <- wilcox.test(occ_suit_vals, abs_suit_vals, alternative = "greater", exact = FALSE)
    as.numeric(wt$statistic) / (n_pres * n_abs)
  }, error = function(e) NA_real_)
  
  
  #---------------
  #- Boyce index -
  #---------------
  boyce_val<-compute_boyce_robust(all_suit_vals, occ_suit_vals)
  
  
  #---------------
  #----- tss -----
  #---------------
  # max TSS: sweep candidate thresholds = unique sorted predicted values
  tss_val <- tryCatch({
    all_preds   <- c(occ_suit_vals, abs_suit_vals)
    all_labels  <- c(rep(1L, n_pres), rep(0L, n_abs))
    ord         <- order(all_preds, decreasing = TRUE)
    preds_s     <- all_preds[ord]
    labels_s    <- all_labels[ord]
    
    # Cumulative TP and FP as threshold descends
    tp_cum <- cumsum(labels_s == 1L)
    fp_cum <- cumsum(labels_s == 0L)
    sens   <- tp_cum / n_pres      # sensitivity 
    spec   <- 1 - fp_cum / n_abs   # specificity
    tss_v  <- sens + spec - 1
    
    # Identify the optimal threshold
    best_idx <- which.max(tss_v)
    
    list(TSS        = tss_v[best_idx],
         sensitivity = sens[best_idx],
         specificity = spec[best_idx],
         threshold   = preds_s[best_idx]
    )
  }, error = function(e) {
    list(TSS = NA_real_,
         sensitivity = NA_real_,
         specificity = NA_real_,
         threshold = NA_real_
    )
  })
  
  
  #------------------
  #- Return metrics -
  #------------------
  return(data.frame(Species = species,
                    Type = type,
                    Region = region,
                    test_fold = fold, 
                    n_pres = as.numeric(n_pres),
                    n_abs = as.numeric(n_abs),
                    auc = as.numeric(auc_val),
                    boyce = as.numeric(boyce_val),
                    tss = as.numeric(tss_val$TSS),
                    sens = as.numeric(tss_val$sensitivity),
                    spec = as.numeric(tss_val$specificity)))
 
}


#-----------------------------------------------------------------
#--Extract raster values at presence/absence points and return----
#-----------------------------------------------------------------
extract_env <- function(pres_abs_points, raster) {
  
  env_values <- terra::extract(raster,
                               terra::vect(pres_abs_points),
                               ID = FALSE,
                               xy = FALSE)
  
  # Combine with species label 
  df <- cbind(species = pres_abs_points$species, 
              ID = pres_abs_points$ID,
              env_values)

  #remove NA rows
  df <- df[complete.cases(df), ]
  
  list(
    presences = df[df$species == 1, -1],
    absences  = df[df$species == 0, -1]
  )
}


#-----------------------------------------------------------------
#--Calculate the geometric mean for ensemble validation------------
#-----------------------------------------------------------------
ensemble_geom_mean <- function(hab_df, clim_df, type,value_col_hab = "median_favourability", value_col_clim = "median_favourability") {
  
  merged <- dplyr::inner_join(hab_df, clim_df, by = "ID", suffix = c("_hab", "_clim"))
  
  geom_mean<-sqrt(merged[[paste0(value_col_hab, "_hab")]] *
         merged[[paste0(value_col_clim, "_clim")]])
  
  #Create a message for printing
  n_all  <- max(nrow(hab_df), nrow(clim_df))
  n_merged <- nrow(merged)
  perc_retained <- 100 * n_merged / n_all
  n_lost <- n_all  - n_merged

  
  message(sprintf(paste("Ensemble validation:", "%d European",type ,"points retained (%.0f%%).",
    "%d point(s) excluded due to missing predictions in the climate or habitat model."
  ), n_merged, perc_retained, n_lost))
  
  return(geom_mean)
}


#-----------------------------------------------------------------
#--Put NO CV validation data in right format ----------
#-----------------------------------------------------------------
summarise_validation<- function(df, validation){
  
  #Check that all necessary columns are present
  required_cols <- c("Species", "Type", "Region", "test_fold","auc", "boyce", "tss", "sens", "spec")
  missing_cols <- setdiff(required_cols, names(df))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  if(validation=="Cross-validation"){
    
    df<-df %>%
      dplyr::group_by(Species, Type, Region)%>%
      dplyr::summarise(n_folds   = n_distinct(test_fold),
                       mean_auc  = mean(auc, na.rm = TRUE),
                       sd_auc    = sd(auc, na.rm = TRUE),
                       mean_boyce = mean(boyce, na.rm = TRUE),
                       sd_boyce   = sd(boyce, na.rm = TRUE),
                       mean_tss  = mean(tss, na.rm = TRUE),
                       sd_tss    = sd(tss, na.rm = TRUE),
                       mean_sens  = mean(sens, na.rm = TRUE),
                       sd_sens    = sd(sens, na.rm = TRUE),
                       mean_spec  = mean(spec, na.rm = TRUE),
                       sd_spec = sd(spec, na.rm = TRUE))%>%
      dplyr::mutate(validation = "Cross-validation")
  }else{
  #Prepare data in right format
  df<- df %>%
  dplyr::transmute(Species,
                   Type,
                   Region,
                   n_folds = NA_real_,
                   mean_auc = auc,
                   sd_auc = NA_real_,
                   mean_boyce = boyce,
                   sd_boyce = NA_real_,
                   mean_tss = tss,
                   sd_tss = NA_real_,
                   mean_sens = sens,
                   sd_sens = NA_real_,
                   mean_spec = spec,
                   sd_spec= NA_real_,
                   validation = "No cross-validation")
 
  }
  return(df)
}
