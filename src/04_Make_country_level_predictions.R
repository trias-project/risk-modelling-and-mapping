#--------------------------------------------
#-----------To do: specify project-----------
#--------------------------------------------
#specify project name
projectname<-"Project_Frédérique"


#--------------------------------------------
#-----------  Load packages  ----------------
#--------------------------------------------
packages <- c("viridis", "geoR","ape","dplyr", "grid", "here", "qs","terra", "sf", "purrr", "progressr", "caret",
              "ggplot2","RColorBrewer","magick","patchwork","tidyterra", "gbm","kableExtra", "rnaturalearth"
)

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#-- Load or install moranfast from github  --
#--------------------------------------------
if (!"moranfast" %in% rownames(installed.packages())) {
  if (!"devtools" %in% rownames(installed.packages())) {
    install.packages("devtools")
  }
  library(devtools)  
  install_github("mcooper/moranfast")  # Install moranfast from GitHub
}
library(moranfast)


#--------------------------------------------
#---  Load right version of caretEnsemble  --
#--------------------------------------------
desired_version <- "2.0.3"

# Check if caretEnsemble is installed
if ("caretEnsemble" %in% rownames(installed.packages())) {
  # Get the current version of caretEnsemble
  current_version <- packageVersion("caretEnsemble")
  # Compare current version with the desired version
  if (as.character(current_version) != desired_version) {
    # Uninstall the current version if it's not the desired version
    remove.packages("caretEnsemble")
    # Install the specific version
    devtools::install_github("zachmayer/caretEnsemble@2.0.3")
    # Load 
    library(caretEnsemble)
  } else {
    library(caretEnsemble)
  }
  
} else {
  # If caretEnsemble is not installed, install the specific version
  devtools::install_github("zachmayer/caretEnsemble@2.0.3")
  # Load 
  library(caretEnsemble)
  rm(current_version, desired_version)
}


#--------------------------------------------
#----------  To do: specify country  --------
#--------------------------------------------
#If you'd like to predict for another country, change the shapefile
country_name<-"Belgium"
country<-sf::st_read(here("./data/external/GIS/Belgium/belgium_boundary.shp"))
country_ext<-terra::ext(country) 
country_vector <- terra::vect(country) #Convert to a SpatVector, used for masking


#--------------------------------------------
#------------ Load country shape ------------
#--------------------------------------------
#Only used when no EU model could be fitted
belgium<-rnaturalearth::ne_countries(country="Belgium", scale=10)[1]
belgium_ext<-terra::ext(belgium) 
belgium_vector <- terra::vect(belgium) #Convert to a SpatVector, used for masking


#--------------------------------------------
#----------- Load taxa info  ----------------
#--------------------------------------------
taxa_info<-read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  dplyr::pull(speciesKey)%>%
  unique()


#-------------------------------------------------
#---------- Load habitat raster data -------------
#-------------------------------------------------
#ONE LAYER HAS A SLIGHTLY DIFFERENT EXTENT: CUT ALL OTHERS TO THIS EXTENT
habitat<-list.files((here("./data/external/habitat")),pattern='tif',full.names = T)
habitat_stack<-terra::rast(habitat[c(1:5,7)]) #Distance to water (layer 6) has another extent and we're not sure whether it is correct: leave it out!
#habitat_stack2<-rast(habitat[6])

#habitat_stack1<-crop(habitat_stack1, ext(habitat_stack2))
#habitat_stack<-c(habitat_stack1, habitat_stack2)


#-------------------------------------------------
#---------- Define resampling raster -------------
#-------------------------------------------------
#Used as a template when global model predictions need to be in same raster format as those for the EU model
resampling_raster<-rast(here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26/anngdd100.tif"))


#--------------------------------------------
#--------Source helper functions-------------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#------ Create start df for model info-------
#--------------------------------------------
model_info<-taxa_info%>%
  dplyr::select(speciesKey, acceptedScientificName)%>%
  dplyr::mutate(Final_model=NA,
         n_presences = NA,
         Threshold= NA,
         Specificity=NA,
         Sensitivity=NA,
         AUC=NA,
         PCC=NA,
         Kappa =NA,
         Morans_I_method=NA,
         Morans_I = NA,
         Pvalue_Morans_I= NA,
         correlation_glm_gbm=NA,
         correlation_glm_rf=NA,
         correlation_glm_earth=NA,
         correlation_gbm_rf=NA,
         correlation_gbm_earth=NA,
         correlation_rf_earth=NA)

#--------------------------------------------
#-----------  Start loop   ------------------
#--------------------------------------------
with_progress({
  p <- progressr::progressor(along = 1:length(accepted_taxonkeys)) 
for(key in accepted_taxonkeys){
  
  p()
  #--------------------------------------------
  #-------  Extract species data   ------------
  #--------------------------------------------
  species<-taxa_info%>%
    dplyr::filter(speciesKey==key)%>%
    dplyr::pull(acceptedScientificName)%>%
    unique()
  
  #Extract first two words of species name 
  first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
  
  #Extract rest of species name
  rest_of_name <- if (grepl("^\\S+\\s+\\S+$", species)) "" else sub("^\\S+\\s+\\S+\\s+", "", species)
  
  #Specify species for plot title
  species_title <- gsub("_", " ", first_two_words)
  
  #Define taxonkey
  taxonkey<- key
  
  
  #--------------------------------------------
  #----- Specify folder and file paths --------
  #--------------------------------------------
  raster_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Rasters")
  raster_country_folder<-file.path(raster_folder, country_name)
  PDF_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PDFs")
  PDF_country_folder<-file.path(PDF_folder, country_name)
  PNG_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PNGs")
  PNG_country_folder <- file.path(PNG_folder, country_name)
  VarImp_folder <-file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Variable_importance")
  RC_folder<-file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Response_curves")
  eu_model_file<-file.path("./data/projects",projectname,paste0(first_two_words,"_",taxonkey), paste0("EU_model_",first_two_words,"_",taxonkey,".qs"))
  global_model_file<-file.path("./data/projects",projectname, paste0(first_two_words,"_",taxonkey),paste0("Global_model_",first_two_words,"_",taxonkey,".qs"))
  
  
  #-----------------------------------------------------------------------------------
  #-Create folders for variable importance, response curves, and Belgian rasters -----
  #-----------------------------------------------------------------------------------
  # Define the folder paths
  folder_paths<-list(list("path"=VarImp_folder,
                          "name"= "Variable_importance"),
                     list("path"= RC_folder,
                          "name"= "Response_curves"),
                     list("path"= raster_country_folder,
                          "name"= paste0("Rasters/",country_name)),
                     list("path"= PDF_country_folder,
                          "name"= paste0("PDF/",country_name)),
                     list("path"= PNG_country_folder,
                          "name"= paste0("PNG/",country_name)))
  
  # Check and create each folder if necessary
  lapply(folder_paths, function(folder){
    create_folder(folder$path, folder$name)
  })
  
  
  #--------------------------------------------
  #--------- Read in model data   -------------
  #--------------------------------------------
  #Use European model if it could be fitted, otherwise use the global model
  if(file.exists(eu_model_file)){
    
    #-------Read in eu model object that was stored as part of  script 03_fit_European_model and load data-------
    eumodel<-qs::qread(eu_model_file)
    euocc<-eumodel$euocc1
    bestModel<-terra::unwrap(eumodel$bestModel)
    eu_presabs.coord<-eumodel$eu_presabs.coord
    occ.full.data.forCaret<-eumodel$occ.full.data.forCaret
    model_correlation<-eumodel$model_correlation
    available_models <- rownames(model_correlation)
    Final_model<-"EU model"
    
    #---------Load raster with EU predictions using the best model------------
    #Define file paths
    EU_predictions_file<-file.path(raster_folder,"Europe",paste(first_two_words,"_",taxonkey,"_hist_EU.tif",sep=""))
    fullstack_be_file<- file.path(raster_folder,"Interim", paste0("Fullstack_be_",first_two_words,"_",taxonkey,".tif"))
    
    #Load files
    fullstack_be<-terra::rast(fullstack_be_file)
    ens_pred_hab_eu1<-terra::rast(EU_predictions_file)
    
    #-----------Print statement----------
    print(paste("Using EU model for species",first_two_words))
    
  }else if(file.exists(global_model_file)){
    #Print warning
    warning(paste0("Using global model for ", species, " because no EU model could be fitted"))
    
    #---------Read in global model object that was stored as part of  script 02_fit_global_model and load data----------
    globalmodels<-qs::qread(global_model_file)
    euocc<-globalmodels$occurrences
    bestModel<-globalmodels$global_ensemble_model
    eu_presabs.coord<-globalmodels$global_presabs
    occ.full.data.forCaret<-globalmodels$global_data_df_uncor
    model_accuracy<-globalmodels$model_accuracy
    model_correlation<-globalmodels$model_correlation
    available_models <- rownames(model_correlation)
    Final_model<-"Global model"
    
    #----------Load rasterlayers--------------
    #Define file path
    global_predictions<-file.path(raster_folder,"Global",paste("Global_model_",first_two_words,"_",taxonkey,".tif",sep=""))
    eu_climpreds10_file <-file.path("./data/projects",projectname,paste0(first_two_words,"_",taxonkey),"Rasters","Interim",paste0("EU_climpreds10_",first_two_words,"_",taxonkey,".tif"))
    
    #Load and process rasters (no need to add habitat stack as these were not used in global model)
    #I decided to not project the data to the crs of country as the model has been fitted on WGS84 and reprojecting can create tiny changes in the pixel values
    ens_pred_hab_eu1<-terra::rast(global_predictions)
    
    fullstack_be<-terra::rast(eu_climpreds10_file)%>%
      #project(crs(country))%>%
      #resample( habitat_stack, method="bilinear")%>% #Make sure the climatic layers have the same resolution (1000 1000)and align with the habitat stack layer or the climatic layers used for the eu model (interpolation)
      terra::crop(belgium_ext)%>%
      terra::mask(belgium_vector)
    
    #---------------Print statement-----------
    print(paste("Using global model for species",first_two_words))
    
  }else{
    warning(paste0("Skipping species ", species, " because no EU model or global model could be fitted"))
    model_info[model_info$speciesKey == key, ]$Final_model<-"None fitted"
    next  # Skip the rest of the loop and move to the next iteration
    
  }
  

  #--------------------------------------------
  #------- Get model performance ---------
  #--------------------------------------------
  #identify threshold where sensitivity=specifity
  thresholds<- findThresh(bestModel$ens_model$pred)
  
  #Using thresholds identified for each model in the previous step, assess performance of each model
  # accuracy measures
  thresholds.df<-accuracyStats(bestModel$ens_model$pred,thresholds$predicted) 
  
  
  #--------------------------------------------
  #------- Subset country occurrences ---------
  #--------------------------------------------
  #occ.eu is in WGS84, convert to same projection as country level shapefile (which is the same proj used for model outputs)
  #suppressWarnings({
    #occ.country <- euocc %>%
     # st_transform(st_crs(country)) %>%
     # st_intersection(country) %>%
     # mutate(coords = st_coordinates(.)) %>%
     # transmute(geometry, 
               # decimalLongitude = coords[,1], 
                #decimalLatitude = coords[,2])
  #}) 
  
  
  #--------------------------------------------
  #-------- Plot country occurrences ---------
  #--------------------------------------------
  #ggplot()+ 
  #geom_sf(data = country,  colour = "black", fill = "white")+
  #geom_point(data=occ.country, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
  #labs(x="Longitude", y="Latitude")+
  #theme_bw()
  
  
  #--------------------------------------------
  #-Create country predictions using best model -
  #--------------------------------------------
  # creates  country level rasters using the European level models
  #ens_pred_hab_be<-terra::predict(fullstack_be,bestModel,type="prob", na.rm=TRUE)
  
  
  #--------------------------------------------
  #-------- ¨Plot predictions for country -----
  #--------------------------------------------
  #brks <- seq(0, 1, by=0.1)
  #nb <- length(brks) - 1
  #viridis_palette <- viridis(nb)
  
  #ggplot() + 
  #geom_spatraster(data = ens_pred_hab_be) +
  #scale_fill_gradientn(colors = viridis_palette, 
  #                    breaks = brks, 
  #                   labels = brks, 
  #                  na.value = NA) +
  #geom_sf(data = occ.country, color = "black", fill = "red", 
  #       size = 1.5, shape = 21) +
  #theme_bw() +
  #labs(fill = "Suitability")
  
  
  #-------------------------------------------------
  #- Prepare layers for fitting --
  #-------------------------------------------------
  
  if(file.exists(eu_model_file)){
  
  #- Clip habitat raster stack to extent of country --
  habitat_only_stack<-terra::crop(habitat_stack,country)
  habitat_only_stack_be<-terra::mask(habitat_only_stack,country)
  
  
  #-Create individual RCP climate raster stacks for country --
  be26 <- list.files((here::here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26")),pattern='tif',full.names = T)
  belgium_stack26 <- terra::rast(be26)
  
  be45 <- list.files((here::here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45")),pattern='tif',full.names = T)
  belgium_stack45 <- terra::rast(be45)
  
  be85 <- list.files((here::here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85")),pattern='tif',full.names = T)
  belgium_stack85 <- terra::rast(be85)
  
  
  #-Combine habitat stacks with climate stacks for each RCP scenario --
  fullstack26_list <- list(belgium_stack26,habitat_only_stack_be)
  fullstack26 <- terra::rast(fullstack26_list) 
  
  fullstack45_list <- list(belgium_stack45,habitat_only_stack_be)
  fullstack45 <- terra::rast(fullstack45_list) 
  
  fullstack85_list <- list(belgium_stack85,habitat_only_stack_be)
  fullstack85 <- terra::rast(fullstack85_list) 
  
  
  }else if(file.exists(global_model_file)){
    
    #-Only use future climate layers--
    be26 <- list.files((here("./data/external/climate/Global_finalRCP/belgium_rcps/rcp26")),pattern='tif',full.names = T)
    fullstack26 <- terra::rast(be26)
    
    be45 <- list.files((here("./data/external/climate/Global_finalRCP/belgium_rcps/rcp45")),pattern='tif',full.names = T)
    fullstack45 <- terra::rast(be45)
    
    be85 <- list.files((here("./data/external/climate/Global_finalRCP/belgium_rcps/rcp85")),pattern='tif',full.names = T)
    fullstack85 <- terra::rast(be85)
    
  }
  
  
  #-------------------------------------------------
  #------------- Store layers in a list ------------
  #-------------------------------------------------
  country_layers<-list(
    "historical"=list("layers"=fullstack_be,
                      "scenario"= "hist",
                      "scenario_title"="historical"),
    "rcp26"=list("layers"=fullstack26,
                 "scenario"="rcp26",
                 "scenario_title"="RCP 2.6"), 
    "rcp45"=list("layers"=fullstack45,
                 "scenario"="rcp45",
                 "scenario_title"="RCP 4.5"),
    "rcp85"=list("layers"=fullstack85,
                 "scenario"="rcp85",
                 "scenario_title"="RCP 8.5")
  )
  
  
  #--------------------------------------------------------------------
  #------- Create and export risk maps for each scenario ------
  #--------------------------------------------------------------------
  pred_list<-lapply(country_layers,function(country_layer) {
    
    # Extract layers and scenario
    layers <- country_layer[[1]]  # climatic and habitat layers
    scenario <- country_layer[[2]]    # scenario name
    
    # Check model predictors
    required_vars <- bestModel$models$glm$coefnames 
    
    # Subset only the required predictors
    predictors <- layers[[required_vars]]
    
    #Make predictions
    predictions<-terra::predict(predictors,bestModel,type="prob", na.rm=TRUE)
    
    #Reproject predictions if the global model was used
    if( Final_model=="Global model"){
      
      #Resample
      predictions_to_export<-predictions%>%
        terra::project(crs(country))%>%
        terra::resample(resampling_raster, method="bilinear") #Make sure the climatic layers have the same resolution (1000 1000)and align with the habitat stack layer or the climatic layers used for the eu model (interpolation)
      
    }else{
      predictions_to_export<-predictions
    }
    
    #------------Export country predictions as raster and PDF -----------
    #Export raster
    raster_file<-paste(first_two_words, "_", taxonkey, "_", scenario, "_", country_name, ".tif", sep="")
    
    writeRaster(predictions_to_export,
                filename=file.path(raster_country_folder, raster_file),
                overwrite=TRUE)
    
    print(paste(raster_file, "has been created"))
    
    #export PDFs and PNGs
    exportPDF(predictions=predictions_to_export,
              taxonName=first_two_words,
              nameExtension=rest_of_name,
              taxonNameTitle=species_title,
              taxonKey=taxonkey,
              scenario=scenario,
              dataType="Suit",
              regionName=country_name,
              returnPredictions=TRUE,
              returnPNG=TRUE,
              keep_PNG=TRUE)
    
  }
  )
  
  
  
  #--------------------------------------------------------------------
  #-------------- Export PDF with the four predictions ----------------
  #--------------------------------------------------------------------
  #TO SOLVE:different suitability axes are used and plotted when global model is used!!
  plot_combined <- (pred_list[[1]][[2]]+ pred_list[[2]][[2]] + plot_layout(ncol = 2)) /
    (pred_list[[3]][[2]] + pred_list[[4]][[2]] + plot_layout(ncol = 2)) /
    plot_spacer() +
    plot_layout(heights = c(1, 1, 1), guides = "collect") &
    theme(legend.position = "top",
          legend.direction = "horizontal",
          legend.key.width = unit(12, "mm"),
          legend.key.height = unit(6, "mm"))
  
  
  exportPDF(taxonName=first_two_words,
            nameExtension=rest_of_name,
            taxonNameTitle=species_title,
            taxonKey=taxonkey,
            scenario="all",
            dataType="Suit",
            regionName=country_name,
            providePNG=TRUE,
            PNGprovided=plot_combined)
  
  
  #--------------------------------------------------------------------------------
  #-------If global model is used, export predictions that are not resampled ------
  #--------------------------------------------------------------------------------
  #This is done to stay as close as possible to the original model and only resample at the end
 if(Final_model=="Global model"){
   pred_list<-lapply(country_layers,function(country_layer) {
    
    # Extract layers and scenario 
    layers <- country_layer[[1]]  # climatic and habitat layers
    scenario <- country_layer[[2]]    # scenario name
    
    # Check model predictors
    required_vars <- bestModel$models$glm$coefnames 
    
    # Subset only the required predictors
    predictors <- layers[[required_vars]]
    
    #Make predictions
    predictions<-terra::predict(predictors,bestModel,type="prob", na.rm=TRUE)
    
    #Store in list
    return(list("model" = predictions,
                "scenario"=scenario))
   }
  )
 }

  
  #--------------------------------------------------------------------
  #-------------- Create and export "difference maps" ----------------
  #--------------------------------------------------------------------
  #These maps display the difference between future predictions and current (historical) predictions
  diff_list <- lapply(pred_list[2:4], function(x) {
    #Define scenario
    scenario <- x$scenario  # scenario name
    
    #Subtract historical suitabilities from future suitabilities
    pred_diff<-x$model - pred_list$historical$model
    
    #Reproject predictions if the global model was used
    if( Final_model=="Global model"){
      #Resample
      pred_diff_to_export<-pred_diff%>%
        terra::project(crs(country))%>%
        terra::resample(resampling_raster, method="bilinear") #Make sure the climatic layers have the same resolution (1000 1000)and align with the habitat stack layer or the climatic layers used for the eu model (interpolation)
      
    }else{
      pred_diff_to_export<-pred_diff
    }
    
    #------------Export difference as raster and PDF -----------
    #Export raster
    raster_file<-paste(first_two_words, "_", taxonkey, "_", scenario, "_hist_diff_", country_name, ".tif", sep="")
    
    terra::writeRaster(pred_diff_to_export,
                filename=file.path(raster_country_folder, raster_file),
                overwrite=TRUE)
    
    print(paste(raster_file, "has been created"))
    
    #export PDFs and store PNGs in plot_list
    exportPDF(predictions=pred_diff_to_export,
              taxonName=first_two_words,
              nameExtension=rest_of_name,
              taxonNameTitle=species_title,
              taxonKey=taxonkey,
              scenario=scenario,
              regionName=country_name,
              dataType = "Diff",
              returnPredictions=TRUE,
              returnPNG=TRUE)
  }
  )
  
  
  #--------------------------------------------------------------------
  #-------------- Check spatial correlation of residuals----------------
  #--------------------------------------------------------------------
  ### This is done to assess whether occurrence data should be thinned
  # Load dataframe with predictions and observations of best model 
  predEns1<-bestModel$ens_model$pred
  
  #Create vector of observed values (0-1 format)
  obs.numeric<-ifelse(predEns1$obs == "absent",0,1)
  
  #Create a vector with standardized residuals
  hab.res<-stdres(obs.numeric,predEns1$present)
  
  # Load coordinates of presences and absences used to fit the best model (stored in eu_presabs.coord) and combine with explanatory data
  res.best.coords1<-cbind(eu_presabs.coord,occ.full.data.forCaret)
  
  #Remove rows with NAs (in explanatory data)
  removedNAs.coords<-na.omit(res.best.coords1)
  
  #Combine coordinate data with residuals
  res.best.df<-cbind(removedNAs.coords[1:2],hab.res) #Standardized residuals
  res.best.df$res<-obs.numeric-predEns1$present #Regular residuals
  
  #Clean up
  rm(predEns1, obs.numeric, hab.res, res.best.coords1,removedNAs.coords)
  gc()
  
  #---------- Check Morans I.- THIS CODE FAILS FOR LARGE DATASETS: NEEDS TOO MUCH RAM
  #If Moran's I is very low (<0.10) or not significant, do not need to thin occurrences.
  #If the observed value is significantly greater than the expected value then there is positive autocorrelation
  #https://stats.oarc.ucla.edu/r/faq/how-can-i-calculate-morans-i-in-r/
  
  if(nrow(res.best.df)<10000){
    occ.dists <- as.matrix(dist(res.best.df[1:2]))
    occ.dists.inv <- 1/occ.dists
    diag(occ.dists.inv) <- 0
    autocor_result<-ape::Moran.I(res.best.df$hab.res,occ.dists.inv,scaled=TRUE,alternative="greater") #‘greater’ evaluates whether the data exhibit more spatial autocorrelation than expected
    MoransI_method<-"Original (ape)"
    observedmoransI<-autocor_result$observed #Observed Moran's I
    pvalue_moransI<-autocor_result$p.value # p-value
    
  }else{
    ##ALTERNATIVE CODE from: https://github.com/mcooper/moranfast
    #Note that it calculates it a bit differently from the code in ape
    autocor_result<-moranfast::moranfast(res.best.df$hab.res, res.best.df$x, res.best.df$y, alternative="greater")
    observedmoransI<-autocor_result$observed #Observed Moran's I
    pvalue_moransI<-autocor_result$p.value # p-value
    MoransI_method<-"C++ alternative (moranfast)"
  }
  
  #-----------------------------------------------------------------------------
  #--Quantify confidence of predicted values using class conformal prediction --
  #-----------------------------------------------------------------------------
  
  # quantify confidence for country-level predictions based on historical climate and under RCP scenarios of climate change
  set.seed(1609)
  pvalsdf_hist<-classConformalPrediction(bestModel,pred_list$historical$model)
  set.seed(447)
  pvalsdf_rcp26<-classConformalPrediction(bestModel,pred_list$rcp26$model)
  set.seed(568)
  pvalsdf_rcp45<-classConformalPrediction(bestModel,pred_list$rcp45$model)
  set.seed(988)
  pvalsdf_rcp85<-classConformalPrediction(bestModel,pred_list$rcp85$model)
  
  # option to export confidence and pvals as csv 
  # write.csv(pvalsdf_hist,file=paste(genOutput,"confidence_",taxonkey, "_hist.csv",sep=""))
  
  
  #-----------------------------------------------------------------------------
  #--------------- Create and export confidence maps --------------------------
  #-----------------------------------------------------------------------------
  if(Final_model!="Global model"){
  hist.conf.map<-confidenceMaps(x=pvalsdf_hist,original_raster=pred_list$historical$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="hist", regionName=country_name, folder=raster_country_folder)
  rcp26.conf.map<-confidenceMaps(x=pvalsdf_rcp26,original_raster=pred_list$rcp26$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="rcp26", regionName=country_name, folder=raster_country_folder)
  rcp45.conf.map<-confidenceMaps(x=pvalsdf_rcp45,original_raster=pred_list$rcp45$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="rcp45", regionName=country_name, folder=raster_country_folder)
  rcp85.conf.map<-confidenceMaps(x=pvalsdf_rcp85,original_raster=pred_list$rcp85$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="rcp85", regionName=country_name, folder=raster_country_folder)
  }else{
    hist.conf.map<-confidenceMaps(x=pvalsdf_hist,original_raster=pred_list$historical$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="hist", regionName=country_name, folder=raster_country_folder, GlobalModel=TRUE, resampling_rast=resampling_raster, country_sf=country)
    rcp26.conf.map<-confidenceMaps(x=pvalsdf_rcp26,original_raster=pred_list$rcp26$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="rcp26", regionName=country_name, folder=raster_country_folder,GlobalModel=TRUE, resampling_rast=resampling_raster, country_sf=country)
    rcp45.conf.map<-confidenceMaps(x=pvalsdf_rcp45,original_raster=pred_list$rcp45$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="rcp45", regionName=country_name, folder=raster_country_folder,GlobalModel=TRUE, resampling_rast=resampling_raster, country_sf=country)
    rcp85.conf.map<-confidenceMaps(x=pvalsdf_rcp85,original_raster=pred_list$rcp85$model, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="rcp85", regionName=country_name, folder=raster_country_folder, GlobalModel=TRUE,resampling_rast=resampling_raster, country_sf=country) 
    
  }
  
  #-----------------------------------------------------------------------------
  #-------------- Mask areas below a set confidence level ----------------------
  #-----------------------------------------------------------------------------
  
  # Cutoff for "high" confidence can be modified below. Cutoff should be a value between 0 and 1. Values that are less than the cutoff are shown in gray.
  cutoff<-0.70
  
  m1<-hist.conf.map < cutoff
  hist_masked<-terra::mask(pred_list$historical$model,m1,maskvalue=TRUE)
  
  m2<-rcp26.conf.map < cutoff
  rcp26_masked<-terra::mask(pred_list$rcp26$model,m2,maskvalue=TRUE)
  
  m3<-rcp45.conf.map < cutoff
  rcp45_masked<-terra::mask(pred_list$rcp45$model,m3,maskvalue=TRUE)
  
  m4<-rcp85.conf.map < cutoff
  rcp85_masked<-terra::mask(pred_list$rcp85$model,m4,maskvalue=TRUE)
  
  
  #-----------------------------------------------------------------------------
  #-------------------- Export rasterfiles of masked predictions ------------------
  #-----------------------------------------------------------------------------
  
  masked_maps<-list(
    "historical"=list("layers"=hist_masked,
                      "scenario"= "hist",
                      "scenario_title"="historical"),
    "rcp26"=list("layers"=rcp26_masked,
                 "scenario"="rcp26",
                 "scenario_title"="RCP 2.6"), 
    "rcp45"=list("layers"=rcp45_masked,
                 "scenario"="rcp45",
                 "scenario_title"="RCP 4.5"),
    "rcp85"=list("layers"=rcp85_masked,
                 "scenario"="rcp85",
                 "scenario_title"="RCP 8.5")
  )
  
  walk(masked_maps, function(x) {
    scenario<-x$scenario
    masked_confidencemap<- x$layers
    
    #Reproject predictions if the global model was used
    if( Final_model=="Global model"){
      masked_confidencemap<- masked_confidencemap%>%
        terra::project(terra::crs(country))%>%
        terra::resample(resampling_raster, method="bilinear") #Make sure the climatic layers have the same resolution (1000 1000)and align with the habitat stack layer or the climatic layers used for the eu model (interpolation)
      
    }
    
    #export raster
    raster_file<-paste(first_two_words, "_", taxonkey, "_", scenario, "_",cutoff,"_masked_predictions_", country_name, ".tif", sep="")
    terra::writeRaster(masked_confidencemap,
                filename=file.path(raster_country_folder, raster_file),
                overwrite=TRUE)
    
    print(paste(raster_file, "has been created"))
    
    #export PDFs and store PNGs in plot_list
    exportPDF(predictions=masked_confidencemap,
              taxonName=first_two_words,
              nameExtension=rest_of_name,
              taxonNameTitle=species_title,
              taxonKey=taxonkey,
              scenario=scenario,
              regionName=country_name,
              dataType = "Masked_Suit",
              returnPredictions=FALSE,
              returnPNG=FALSE)
  }
  )
  
  
  #-----------------------------------------------------------------------------
  #-------------- Get variable importance of best european model ---------------
  #-----------------------------------------------------------------------------
  #Variable importance for each model is calculated and then averaged by the weight of the overall model in the ensembled object.
  variableImportance<-caret::varImp(bestModel)
  #kable(variableImportance,digits=2,caption="Variable Importance") %>%
  #kable_styling(bootstrap_options = c("striped"))
  
  #Store variable importance
  Varimp_path<-file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PDFs")
  write.csv(variableImportance,file = file.path(VarImp_folder, paste0(taxonkey,"_varImp_EU_model.csv")))
  
  
  #-----------------------------------------------------------------------------
  #-------------- Store model info  ---------------
  #-----------------------------------------------------------------------------
  model_info[model_info$speciesKey == key, ]$Final_model<-Final_model
  model_info[model_info$speciesKey == key, ]$Threshold<-thresholds.df$threshold
  model_info[model_info$speciesKey == key, ]$AUC<-thresholds.df$AUC
  model_info[model_info$speciesKey == key, ]$PCC<-thresholds.df$PCC
  model_info[model_info$speciesKey == key, ]$Sensitivity<-thresholds.df$sensitivity
  model_info[model_info$speciesKey == key, ]$Specificity<-thresholds.df$specificity
  model_info[model_info$speciesKey == key, ]$Kappa<-thresholds.df$Kappa
  model_info[model_info$speciesKey==key,]$Morans_I_method<-MoransI_method
  model_info[model_info$speciesKey==key,]$Morans_I<-observedmoransI
  model_info[model_info$speciesKey==key,]$Pvalue_Morans_I<-pvalue_moransI
  model_info[model_info$speciesKey==key,]$n_presences<-nrow(euocc)
  model_info[model_info$speciesKey==key,]$correlation_glm_gbm<-if(all(c("glm", "gbm") %in% available_models)) model_correlation["glm", "gbm"] else NA
  model_info[model_info$speciesKey==key,]$correlation_glm_rf<-if(all(c("glm", "rf") %in% available_models)) model_correlation["glm", "rf"] else NA
  model_info[model_info$speciesKey==key,]$correlation_glm_earth<-if(all(c("glm", "earth") %in% available_models)) model_correlation["glm", "earth"] else NA
  model_info[model_info$speciesKey==key,]$correlation_gbm_rf<-if(all(c("gbm", "rf") %in% available_models)) model_correlation["gbm", "rf"] else NA
  model_info[model_info$speciesKey==key,]$correlation_gbm_earth<-if(all(c("gbm", "earth") %in% available_models)) model_correlation["gbm", "earth"] else NA
  model_info[model_info$speciesKey==key,]$correlation_rf_earth<-if(all(c("rf", "earth") %in% available_models)) model_correlation["rf", "earth"] else NA
  

  #--------------------------------------------
  #-------- End of loop -----------------------
  #--------------------------------------------
  print(paste("Country-level predictions were created for", species))
  rm(list = setdiff(ls(), c("p","get.confidence","extractVals","GetLength","resampling_raster","classConformalPrediction","projectname", "stdres","exportPDF", "create_folder", "country","accuracyStats","country_name", "country_ext", "country_vector", "habitat_stack", "model_info", "accepted_taxonkeys", "taxa_info", "key", "findThresh", "confidenceMaps", "belgium_ext", "belgium_vector")))

}
})
write.csv(model_info,file = file.path("./data/projects",projectname,"Overview_model_performance.csv"))





# #-----------------------------------------------------------------------------
# #------Generate and export response curves in order of variable importance ---
# #-----------------------------------------------------------------------------
# topPreds <- variableImportance[with(variableImportance,order(-overall)),]
# varNames<-rownames(topPreds)
# 
# 
# ## combine predictions from each model for each variable
# ## train data needs to be the training data used in the individual models used to build the ensemble model. This info can be extracted from the best ensemble model (ie. bestModel)
# bestModel.train<-bestModel$models[[1]]$trainingData
# 
# 
# gbm.partial.list<-lapply(varNames,partial_gbm)
# glm.partial.list<-lapply(varNames,partial_glm)
# rf.partial.list<-lapply(varNames,partial_rf)
# mars.partial.list<-lapply(varNames,partial_mars)
# 
# names(glm.partial.list)<-varNames
# names(gbm.partial.list)<-varNames
# names(rf.partial.list)<-varNames
# names(mars.partial.list)<-varNames
# 
# glm.partial.df<-as.data.frame(glm.partial.list)
# gbm.partial.df<-as.data.frame(gbm.partial.list)
# rf.partial.df<-as.data.frame(rf.partial.list)
# mars.partial.df<-as.data.frame(mars.partial.list)
# 
# predx<-data.frame()
# predy<-data.frame()
# 
# for (i in varNames){
#   predx <- rbind(predx, as.data.frame(paste(i,i,sep=".")))
#   predy<- rbind(predy,as.data.frame(paste(i,"yhat",sep=".")))
# }
# names(predx)<-""
# names(predy)<-""
# 
# predx1<-t(predx)
# predy1<-t(predy)
# 
# 
# glm.partial.df$data<-'GLM'
# gbm.partial.df$data<-'GBM'
# rf.partial.df$data<-'RF'
# mars.partial.df$data<-'MARS'
# 
# all_dfs<-rbind.data.frame(glm.partial.df,gbm.partial.df,rf.partial.df,mars.partial.df)
# 
# allplots<-map2(predx1,predy1, ~responseCurves(.x,.y))
# 
# #export plots as PNGs
# for(i in seq_along(allplots)){
#   png(file.path(folder_paths[[2]]$path,paste0(taxonkey,"_",i,".png")),width = 5, height = 5, units = "in",res=300)
#   print(allplots[[i]])
#   dev.off()
# }

