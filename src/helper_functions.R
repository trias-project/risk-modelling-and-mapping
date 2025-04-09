
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
      dismo::randomPoints(
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
exportPDF <- function(predictions=NULL,taxonName, nameExtension,taxonNameTitle, taxonKey, scenario, regionName, occ_data=NULL,occ=FALSE,dataType, returnPredictions=FALSE,returnPNG=FALSE, providePNG=FALSE, PNGprovided, keep_PNG=FALSE){
 
  #Define scenario title
  scenario_titles <- c(
    "hist" = "historical",
    "historical" = "historical",
    "rcp26" = "RCP 2.6",
    "rcp45" = "RCP 4.5",
    "rcp85" = "RCP 8.5",
    "all" = "all"
  )
  
  scenarioTitle <- scenario_titles[scenario]
  
  # Define file name suffix and paths based on dataType
  suffix <- switch(dataType,
                   "Suit" = "",
                   "Diff" = "_hist_diff",
                   "Conf" = "_confidence",
                   "Masked_Suit"= "_masked"
  )
  
  # Construct file names
  PNG_file <- paste(taxonName, "_", taxonKey, "_", scenario, suffix, "_", regionName, ".png", sep = "")
  PDF_file <- paste(taxonName, "_", taxonKey, "_", scenario, suffix, "_", regionName, ".pdf", sep = "")
  
  # Define folder paths
  PNG_folder_path <- switch(regionName,
                   "EU"= file.path(PNG_folder,"Europe"),
                   "Global"= file.path(PNG_folder,"Global"),
                   file.path(PNG_folder, regionName))# If not EU or Global
  
  PDF_folder_path <- switch(regionName,
                          "EU"= file.path(PDF_folder,"Europe"),
                          "Global"=file.path(PDF_folder,"Global"),
                          file.path(PDF_folder, regionName))
  
  # Define file paths
  plot_png_path <- file.path(PNG_folder_path, PNG_file)# If not EU or Global
  
  plot_pdf_path <- file.path(PDF_folder_path, PDF_file)

  
  #If png is not provided, create a PNG based on the input predictions
  if(!providePNG){
  
  #Get extent
  exten<-as.vector(terra::ext(predictions))
    
  #Settings for plot
  ifelse(dataType=="Diff", brks<-seq(-1, 1, by=0.25), brks <- seq(0, 1, by=0.1))
  nb <- length(brks) - 1
  viridis_palette <- viridis::viridis(nb)
  
  #Create plot
  country_plot<-ggplot() + 
    geom_spatraster(data = predictions) +
    scale_fill_gradientn(colors = viridis_palette, 
                         breaks = brks, 
                         labels = brks, 
                         na.value = NA) +
    theme_bw() +
    theme(axis.title = element_blank())+
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "cm"))+
    coord_sf(xlim = c(exten[1], exten[2]), 
             ylim = c(exten[3], exten[4] + 15000))
  
  # Define text label, fill label, and hjust based on dataType
  text_label <- ifelse(dataType=="Diff",paste(scenarioTitle, "- historical"), scenarioTitle)
  
  fill_label <- switch(dataType,
                       "Suit" = "Suitability",
                       "Diff" = "Suitability difference",
                       "Conf" = "Confidence",
                       "Masked_Suit" = "Suitability")
  
  hjust_value <- ifelse(dataType == "Diff", 1.1, 1.2)
  
  # Update the plot
  country_plot <- country_plot +
    labs(fill = fill_label) +
    annotate("text",
             x = Inf, y = Inf,       # Position at top-right
             label = text_label,     # Text to display
             hjust = hjust_value,
             vjust = 2.5,            # Adjust text alignment to the right and above
             size = 4.8,
             color = "#636363",
             fontface = "bold")
  
  if(occ){
    country_plot<-country_plot +
      geom_sf(data = occ_data, color = "black", fill = "red", 
              size = 1.5, shape = 21)
  }
  
  #Create an empty plot to fill PDF
  empty_plot <- ggplot() + 
    theme_void() + 
    theme(plot.background = element_blank()) 
  
  #Create final plot
  plot_final<-country_plot /empty_plot 
  
  # Save plot as a PNG file
  ggplot2::ggsave(filename = PNG_file, plot = plot_final, 
         device = "png", width =8.27 , height = 11.69, path= PNG_folder_path)
  
   }else{
     plot_final=PNGprovided
     # Save plot as a PDF file
     ggplot2::ggsave(filename = PNG_file, plot = plot_final, 
            device = "png", width =8.27 , height = 11.69, path= PNG_folder_path)
   }
  
  # Read the PNG image back in
  img <- magick::image_read(plot_png_path)
  
  # Start a PDF device for output
  pdf(plot_pdf_path, width = 8.27, height = 11.69)
  
  # Create a layout for title and image
  grid.newpage()
  
  # Add title at the top of the PDF
  grid.text(
    label = bquote(italic(.(taxonNameTitle)) ~ .(nameExtension) ~ "(" * .(taxonKey) * ")"),
    x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  # Add the PNG image below the title
    grid::grid.raster(img, width = unit(0.9, "npc"), height = unit(0.9, "npc"), y = 0.47)

  # Close the PDF device
  while (dev.cur() > 1) dev.off()
  
  #Print
  print(paste(PDF_file," has been created.", sep=""))
  
  #Remove PNG file if you don't want to keep it (default)
  if(!keep_PNG){
  file.remove(plot_png_path)
  }else{
    #Print
    print(paste(PNG_file," has been created.", sep="")) 
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