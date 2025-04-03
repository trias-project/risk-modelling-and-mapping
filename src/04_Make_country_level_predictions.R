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
  pull(speciesKey)%>%
  unique()


#-------------------------------------------------
#---------- Load habitat raster data -------------
#-------------------------------------------------
#ONE LAYER HAS A SLIGHTLY DIFFERENT EXTENT: CUT ALL OTHERS TO THIS EXTENT
habitat<-list.files((here("./data/external/habitat")),pattern='tif',full.names = T)
habitat_stack<-rast(habitat[c(1:5,7)]) #Distance to water (layer 6) has another extent and we're not sure whether it is correct: leave it out!
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
  select(speciesKey, acceptedScientificName)%>%
  mutate(Final_model=NA,
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
  p <- progressor(along = 1:length(accepted_taxonkeys)) 
for(key in accepted_taxonkeys){
  
  p()
  #--------------------------------------------
  #-------  Extract species data   ------------
  #--------------------------------------------
  species<-taxa_info%>%
    filter(speciesKey==key)%>%
    pull(acceptedScientificName)%>%
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
    eumodel<-qread(eu_model_file)
    euocc<-eumodel$euocc1
    bestModel<-unwrap(eumodel$bestModel)
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
    fullstack_be<-rast(fullstack_be_file)
    ens_pred_hab_eu1<-rast(EU_predictions_file)
    
    #-----------Print statement----------
    print(paste("Using EU model for species",first_two_words))
    
  }else if(file.exists(global_model_file)){
    #Print warning
    warning(paste0("Using global model for ", species, " because no EU model could be fitted"))
    
    #---------Read in global model object that was stored as part of  script 02_fit_global_model and load data----------
    globalmodels<-qread(global_model_file)
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
    ens_pred_hab_eu1<-rast(global_predictions)
    
    fullstack_be<-rast(eu_climpreds10_file)%>%
      #project(crs(country))%>%
      #resample( habitat_stack, method="bilinear")%>% #Make sure the climatic layers have the same resolution (1000 1000)and align with the habitat stack layer or the climatic layers used for the eu model (interpolation)
      crop(belgium_ext)%>%
      mask(belgium_vector)
    
    #---------------Print statement-----------
    print(paste("Using global model for species",first_two_words))
    
  }else{
    warning(paste0("Skipping species ", species, " because no EU model or global model could be fitted"))
    model_info[model_info$speciesKey == key, ]$Final_model<-"None fitted"
    next  # Skip the rest of the loop and move to the next iteration
    
  }
  
#--------------------------------------------
#-------- Plot country occurrences ---------
#--------------------------------------------
ggplot()+ 
    geom_sf(data = country,  colour = "black", fill = NA)+
    geom_point(data=occ.country, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
    labs(x="Longitude", y="Latitude")+
    theme_bw()


  #--------------------------------------------
  #-Create country predictions using best model -
  #--------------------------------------------
  # creates  country level rasters using the European level models
  system.time({
    ens_pred_hab_be<-terra::predict(fullstack_be,bestModel,type="prob", na.rm=TRUE)
  })
  
  
  #--------------------------------------------
  #-------- ¨Plot predictions for country -----
  #--------------------------------------------
  brks <- seq(0, 1, by=0.1)
  nb <- length(brks) - 1
  viridis_palette <- viridis(nb)
  
  country_plot<-ggplot() + 
    geom_spatraster(data = ens_pred_hab_be) +
    scale_fill_gradientn(colors = viridis_palette, 
                         breaks = brks, 
                         labels = brks, 
                         na.value = NA) +
    geom_sf(data = occ.country, color = "black", fill = "red", 
            size = 1.5, shape = 21) +
    theme_bw() +
    labs(fill = "Suitability")
  
  #Create an empty plot to fill PDF
  empty_plot <- ggplot() + 
    theme_void() + 
    theme(plot.background = element_blank()) 
  
  #Create final plot
  plot_final<-country_plot /empty_plot 
  
  
  #-------------------------------------------------
  #- Export country predictions as raster and PDF --
  #-------------------------------------------------
  #---------Specify folder paths------------
  raster_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Rasters")
  PDF_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PDFs")
  
  #---------------Export raster-------------
  writeRaster(ens_pred_hab_be,
              filename=file.path(raster_folder,paste(first_two_words,"_",taxonkey,"_hist_",country_name,".tif",sep="")),
              overwrite=TRUE)
  
  #---------------Export PDF----------------
  #Define the file paths
  plot_png_path <- file.path(PDF_folder,paste(first_two_words,"_",taxonkey,"_hist_",country_name,".png",sep=""))
  plot_pdf_path <- file.path(PDF_folder,paste(first_two_words,"_",taxonkey,"_hist_",country_name,".pdf",sep=""))
  
  # Save each plot as a PDF file
  ggsave(filename = paste0(first_two_words,"_",taxonkey,"_hist_",country_name,".png"), plot = plot_final, 
         device = "png", width =8.27 , height = 11.69, path= PDF_folder)
  
  # Read the PNG image back in
  img <- image_read(plot_png_path)
  
  # Start a PDF device for output
  pdf(plot_pdf_path, width = 8.27, height = 11.69)
  
  # Create a layout for title and image
  grid.newpage()
  
  # Add title at the top of the PDF
  grid.text(
    label = bquote(italic(.(first_two_words)) ~ .(rest_of_name) ~ "(" * .(taxonkey) * ")"),
    x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 12, fontface = "bold")
  )
  
  # Add the PNG image below the title
  grid.raster(img, width = unit(0.9, "npc"), height = unit(0.9, "npc"), y = 0.47)
  
  # Close the PDF device
  while (dev.cur() > 1) dev.off()
  
  # Remove the PNG file from the local directory
  file.remove(plot_png_path)
  
  
  #-------------------------------------------------
  #- Clip habitat raster stack to extent of country --
  #-------------------------------------------------
habitat_only_stack<-terra::crop(habitat_stack,country)
habitat_only_stack_be<-terra::mask(habitat_only_stack,country)


#-----------------------------------------------------------
#-Create individual RCP climate raster stacks for country --
#-----------------------------------------------------------
be26 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26")),pattern='tif',full.names = T)
belgium_stack26 <- rast(be26)

be45 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45")),pattern='tif',full.names = T)
belgium_stack45 <- rast(be45)

be85 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85")),pattern='tif',full.names = T)
belgium_stack85 <- rast(be85)


#--------------------------------------------------------------------
#-Combine habitat stacks with climate stacks for each RCP scenario --
#--------------------------------------------------------------------
fullstack26_list <- list(belgium_stack26,habitat_only_stack_be)
fullstack26 <- rast(fullstack26_list) 

fullstack45_list <- list(belgium_stack45,habitat_only_stack_be)
fullstack45 <- rast(fullstack45_list) 

fullstack85_list <- list(belgium_stack85,habitat_only_stack_be)
fullstack85 <- rast(fullstack85_list) 

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

### Create and export RCP risk maps for each RCP scenario
ens_pred_hist<-raster::predict(fullstack_be,bestModel,type="prob")
ens_pred_hab26<-raster::predict(fullstack26,bestModel,type="prob")
crs(ens_pred_hab26)<-laea_grs80
writeRaster(ens_pred_hab26, filename=file.path(rasterOutput,paste("be_",taxonkey, "_rcp26.tif",sep="")), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab26,taxonkey,taxonName=taxonName,"rcp26.pdf")
ens_pred_hab45<-raster::predict(fullstack45,bestModel,type="prob")
crs(ens_pred_hab45)<-laea_grs80
writeRaster(ens_pred_hab45, filename=file.path(rasterOutput,paste("be_",taxonkey, "_rcp45.tif",sep="")), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab45,taxonkey,taxonName=taxonName,"rcp45.pdf")
ens_pred_hab85<-raster::predict(fullstack85,bestModel,type="prob")
crs(ens_pred_hab85)<-laea_grs80
writeRaster(ens_pred_hab85, filename=file.path(rasterOutput,paste("be_",taxonkey, "_rcp85.tif",sep="")), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab85,taxonkey,taxonName=taxonName,"rcp85.pdf")



### Create and export RCP risk maps for each RCP scenario

par(mfrow=c(2,2), mar= c(2,3,0.8,0.8))
plot(ens_pred_hist,breaks=brks, col=cols,lab.breaks=brks)
plot(ens_pred_hab26,breaks=brks, col=cols,lab.breaks=brks)
plot(ens_pred_hab45,breaks=brks, col=cols,lab.breaks=brks)
plot(ens_pred_hab85,breaks=brks, col=cols,lab.breaks=brks)



### Create and export "difference maps": the difference between predicted risk by each RCP scenario and historical climate
hist26_diff_hab<-overlay(ens_pred_hab26, ens_pred_hist, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist26_diff_hab,filename=file.path(rasterOutput,paste("be_",taxonkey, "_rcp26_diff.tif",sep="")) , format="GTiff",overwrite=TRUE) 
exportPDF(hist26_diff_hab,taxonkey,taxonName=taxonName,"rcp26_diff.pdf","TRUE")


hist45_diff_hab<-overlay(ens_pred_hab45, ens_pred_hist, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist45_diff_hab,filename=file.path(rasterOutput,paste("be_",taxonkey, "_rcp45_diff.tif",sep="")), format="GTiff",overwrite=TRUE) 
exportPDF(hist45_diff_hab,taxonkey,taxonName=taxonName,"rcp45_diff.pdf","TRUE")


hist85_diff_hab<-overlay(ens_pred_hab85, ens_pred_hist, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist85_diff_hab, filename=file.path(rasterOutput,paste("be_",taxonkey, "_rcp_85_diff.tif",sep="")), format="GTiff",overwrite=TRUE) 
exportPDF(hist85_diff_hab,taxonkey,taxonName=taxonName,"rcp85_diff.pdf","TRUE")

par(mfrow=c(2,2), mar= c(2,3,0.8,0.8))
plot(hist26_diff_hab)
plot(hist45_diff_hab)
plot(hist85_diff_hab)




### Check spatial autocorrelation of residuals to assess whether occurrence data should be thinned
#### derive residuals from best model
predEns1<-bestModel$ens_model$pred
obs.numeric<-ifelse(predEns1$obs == "absent",0,1)


#### standardize residuals
stdres<-function(obs.numeric, yhat){
  num<-obs.numeric-yhat
  denom<-sqrt(yhat*(1-yhat))
  return(num/denom)
}
hab.res<-stdres(obs.numeric,predEns1$present)

# specify corresponding model number from eu_presabs.coord datafile to join data with xy locations. If best model is "X1", join with eu_presabs.coord$X1


res.best.coords1<-cbind(coordinates(eu_presabs.coord$X1),occ.full.data.forCaret$X1)
removedNAs.coords<-na.omit(res.best.coords1)
res.best.coords<-cbind(removedNAs.coords,hab.res)
res.best.geo<-as.geodata(res.best.coords,coords.col=1:2,data.col = 3)
summary(res.best.geo) #note distance is in meters


### Check Morans I.

#If Moran's I is very low (<0.10), or not significant, do not need to thin occurrences.
library(ape)
res.best.df<-as.data.frame(res.best.coords)
occ.dists <- as.matrix(dist(cbind(res.best.df[1], res.best.df[2])))
occ.dists.inv <- 1/occ.dists
diag(occ.dists.inv) <- 0
Moran.I(res.best.df$hab.res,occ.dists.inv,scaled=TRUE,alternative="greater")


### Code for Mondrian conformal prediction functions

# functions needed for conformal prediction function


GetLength<-function(x,y){
  length(x[which(x<= y)])
}




CPconf<-function(pA,pB,confidence){
  if(pA > confidence && pB< confidence){
    predClass<-"classA"
  }else if(pA < confidence && pB> confidence){
    predClass<-"classB"
  }else if(pA< confidence && pB< confidence){
    predClass<-"noClass"
  }else{
    predClass<-"bothClasses"
    
    return(predClass)
  }}


#function to calculate confidence of each prediction

get.confidence<-function(pvalA,pvalB){
  secondHighest<-ifelse(pvalA>pvalB,pvalB,pvalA)
  conf<-(1-secondHighest)
  return(conf)
}

forcedCp<-function(pvalA,pvalB){
  ifelse(pvalA>pvalB,"presence","absence")
}

extractVals<-function(predras){
  library(raster)
  vals <-  raster::values(predras)
  coord <-  raster::xyFromCell(predras,1:ncell(predras))
  raster_fitted <- cbind(coord,vals)
  raster_fitted.df<-as.data.frame(raster_fitted)
  raster_fitted.df1<-na.omit(raster_fitted.df)
  raster_fitted.df1$presence<-raster_fitted.df1$vals
  raster_fitted.df1$absence<- (1-raster_fitted.df1$presence)
  return(raster_fitted.df1)
}


classConformalPrediction<-function(x,y){
  ens_results<- get("x")
  ens_calib<-ens_results$ens_model$pred
  calibPresence<-ens_calib %>%
    filter(obs=='present')%>%
    select(present)
  calibPresence<-unname(unlist(calibPresence[c("present")]))
  calibAbsence<-ens_calib %>%
    filter(obs=='absent')%>%
    select(absent)
  calibAbsence<-unname(unlist(calibAbsence[c("absent")]))
  predicted.values<-extractVals(y)
  
  
  testPresence<-predicted.values$presence
  testAbsence<-predicted.values$absence
  
  #derive p.Values for class A
  smallrA<-lapply(testPresence,function(x) GetLength(calibPresence,x))
  smallrA_1<- unlist (smallrA)+1
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
}


### Quantify confidence of predicted values using class conformal prediction


# quantify confidence for country level predictions based on historical climate and under RCP scenarios of climate change

set.seed(1609)
pvalsdf_hist<-classConformalPrediction(bestModel,ens_pred_hist)
set.seed(447)
pvalsdf_rcp26<-classConformalPrediction(bestModel,ens_pred_hab26)
set.seed(568)
pvalsdf_rcp45<-classConformalPrediction(bestModel,ens_pred_hab45)
set.seed(988)
pvalsdf_rcp85<-classConformalPrediction(bestModel,ens_pred_hab85)

# option to export confidence and pvals as csv 
# write.csv(pvalsdf_hist,file=paste(genOutput,"confidence_",taxonkey, "_hist.csv",sep=""))


### Create confidence maps
brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
pal <- colorRampPalette(rev(brewer.pal(4, 'Spectral')))
cols<-pal(nb)


confidenceMaps<-function(x,taxonkey,taxonName,maptype){
  pvals_dataframe<-get("x")
  data.xyz <- pvals_dataframe[c("x","y","conf")]
  rst <- rasterFromXYZ(data.xyz)
  crs(rst)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") 
  plot(rst,breaks=brks, col=cols,lab.breaks=brks)
  writeRaster(rst, filename=file.path(rasterOutput,paste("be_",taxonkey, "_",maptype,".tif",sep="")), format="GTiff",overwrite=TRUE)
  exportPDF(rst,taxonkey,taxonName=taxonName,nameextension= paste(maptype,".pdf",sep=""))
  return(rst)
}

par(mfrow=c(2,2), mar= c(2,3,0.8,0.8))
hist.conf.map<-confidenceMaps(pvalsdf_hist,taxonkey,taxonName,maptype="hist_conf")
rcp26.conf.map<-confidenceMaps(pvalsdf_rcp26,taxonkey,taxonName,maptype="rcp26_conf")
rcp45.conf.map<-confidenceMaps(pvalsdf_rcp45,taxonkey,taxonName,maptype="rcp45_conf")
rcp85.conf.map<-confidenceMaps(pvalsdf_rcp85,taxonkey,taxonName,maptype="rcp85_conf")



### Mask areas of below a set confidence level  

# Cutoff for "high" confidence can be modified below. Cutoff should be a value between 0 and 1. Values that are less than the cutoff are shown in gray.
cutoff<-0.70

conf.brks <- seq(0,1, by=0.1) 
nb <- length(conf.brks) 
pal <- colorRampPalette(rev(brewer.pal(4, 'Spectral')))
cols<-pal(nb)

par(mfrow=c(2,2), mar= c(2,3,0.9,0.8))
m1<-hist.conf.map < cutoff
hist_masked<-mask(ens_pred_hist,m1,maskvalue=TRUE)
plot(hist_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)
plot(country,add=TRUE,border="dark gray")

m2<-rcp26.conf.map < cutoff
rcp26_masked<-mask(ens_pred_hab26,m2,maskvalue=TRUE)
plot(rcp26_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)
plot(country,add=TRUE,border="dark gray")

m3<-rcp45.conf.map < cutoff
rcp45_masked<-mask(ens_pred_hab45,m3,maskvalue=TRUE)
plot(rcp45_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)
plot(country,add=TRUE,border="dark gray")

m4<-rcp85.conf.map < cutoff
rcp85_masked<-mask(ens_pred_hab85,m4,maskvalue=TRUE)
plot(rcp85_masked,breaks=conf.brks, col=cols,lab.breaks=conf.brks)
plot(country,add=TRUE,border="dark gray")

### confidence map of best model at EU level
brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
pal <- colorRampPalette(rev(brewer.pal(4, 'Spectral')))
set.seed(792)  
pvalsdf_hist_eu<-classConformalPrediction(bestModel,ens_pred_hab_eu1$X6)
hist.conf.map.eu<-confidenceMaps(pvalsdf_hist_eu,taxonkey,taxonName,maptype="hist_conf_eu")


### Get variable importance of best european model

variableImportance<-varImp(bestModel)
kable(variableImportance,digits=2,caption="Variable Importance") %>%
  kable_styling(bootstrap_options = c("striped"))
write.csv(variableImportance,file = paste0(genOutput,taxonkey,"_varImp_EU_model.csv"))


### Generate and export response curves in order of variable importance
topPreds <- variableImportance[with(variableImportance,order(-overall)),]
varNames<-rownames(topPreds)
## combine predictions from each model for each variable
## train data needs to be the training data used in the individual models used to build the ensemble model. This info can be extracted from the best ensemble model (ie. bestModel)
bestModel.train<-bestModel$models[[1]]$trainingData

partial_gbm<-function(x){
  m.gbm<-pdp::partial(bestModel$models$gbm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,n.trees= bestModel$models$gbm$finalModel$n.trees, which.class = 1,grid.resolution=nrow(bestModel.train))
}



gbm.partial.list<-lapply(varNames,partial_gbm)

partial_glm<-function(x){
  m.glm<-pdp::partial(bestModel$models$glm$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                      prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

glm.partial.list<-lapply(varNames,partial_glm)

partial_rf<-function(x){
  pdp::partial(bestModel$models$rf$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
               prob=TRUE,which.class = 1,grid.resolution=nrow(bestModel.train))
}

rf.partial.list<-lapply(varNames,partial_rf)


partial_mars<-function(x){
  m.mars<-pdp::partial(bestModel$models$earth$finalModel,pred.var=paste(x),train = bestModel.train,type="classification",
                       prob=TRUE,which.class = 2,grid.resolution=nrow(bestModel.train)) # class=2 because in earth pkg, absense is the first class
}

mars.partial.list<-lapply(varNames,partial_mars)


names(glm.partial.list)<-varNames
names(gbm.partial.list)<-varNames
names(rf.partial.list)<-varNames
names(mars.partial.list)<-varNames

glm.partial.df<-as.data.frame(glm.partial.list)
gbm.partial.df<-as.data.frame(gbm.partial.list)
rf.partial.df<-as.data.frame(rf.partial.list)
mars.partial.df<-as.data.frame(mars.partial.list)

predx<-data.frame()
predy<-data.frame()

for (i in varNames){
  predx <- rbind(predx, as.data.frame(paste(i,i,sep=".")))
  predy<- rbind(predy,as.data.frame(paste(i,"yhat",sep=".")))
}
names(predx)<-""
names(predy)<-""

predx1<-t(predx)
predy1<-t(predy)


glm.partial.df$data<-'GLM'
gbm.partial.df$data<-'GBM'
rf.partial.df$data<-'RF'
mars.partial.df$data<-'MARS'

all_dfs<-rbind.data.frame(glm.partial.df,gbm.partial.df,rf.partial.df,mars.partial.df)


responseCurves<-function(x,y) {
  colors <- c("GLM" = "gray", "GBM"="red","RF"="blueviolet","MARS"= "hotpink") 
  ggplot(all_dfs,(aes(x=.data[[x]],y=.data[[y]]))) +
    geom_line(aes(color = data), size =1.2, position=position_dodge(width=0.2))+
    theme_bw()+
    labs(y="Partial probability", x= gsub("//..*","",x),color="Legend") +
    scale_color_manual(values = colors)
}  

allplots<-map2(predx1,predy1, ~responseCurves(.x,.y))

#export plots as PNGs
for(i in seq_along(allplots)){
  png(paste0(genOutput,taxonkey,"_",i,".png"),width = 5, height = 5, units = "in",res=300)
  print(allplots[[i]])
  dev.off()
}



### Plot response curves

par(mfrow=c(3,4))
for(i in seq_along(allplots)){
  print(allplots[[i]])
}


###  Evaluate the performance of each the EU level ensemble models using independent data set from the future 
#####################################################################

# read in and prepare independent data
#2011-2021
eval.data<-read.csv("C:/Users/amyjs/Documents/projects/xps15/xps15/wiSDM/data/external/0001753-230828120925497/0001753-230828120925497.csv",header=TRUE,sep ="\t",quote="")

#enter value for max coordinate uncertainty in meters.

eval.data.occ<-eval.data %>%
  filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters< 1000) 

eval.data.occ$lon_dplaces<-sapply(na.omit(eval.data.occ$decimalLongitude), function(x) decimalplaces(x))
eval.data.occ$lat_dplaces<-sapply(eval.data.occ$decimalLatitude, function(x) decimalplaces(x))
eval.data.occ[eval.data.occ$lon_dplaces < 4& eval.data.occ$lat_dplaces < 4 , ]<-NA
eval.data.occ<-eval.data.occ[ which(!is.na(eval.data.occ$lon_dplaces)),]
eval.data.occ<-within(eval.data.occ,rm("lon_dplaces","lat_dplaces"))

eval.data.occ<-eval.data.occ[c("decimalLongitude", "decimalLatitude")]
coordinates(eval.data.occ)<- c("decimalLongitude", "decimalLatitude")
proj4string(eval.data.occ)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")#specify here the existing coord.sys of the data
eval.data.occ.proj<-spTransform(eval.data.occ,rmiproj)

########################################################################
#Convert predicted probabilities of EU level risk maps into binary risk maps (present/absence) using thresholds from earlier step

# Eu level
binary_eu_rasters<-sapply(names(thresholds), function(x) raster::reclassify(ens_pred_hab_eu1[[x]],c(0,thresholds[[x]]$predicted,0, thresholds[[x]]$predicted,1,1)),simplify=FALSE)




eu_eval<-function (ras,y){
  indep.bil<-raster::extract(ras,y,method="bilinear")
  indep.bil.df<-as.data.frame(indep.bil)
  indep.bil.df<-indep.bil.df %>%
    mutate(predicted= ifelse(indep.bil >= 0.5,"present","absent")) 
  indep.bil.df$observed<-rep("present",nrow(indep.bil.df))
  indep.bil.df$predicted<-as.factor(indep.bil.df$predicted)
  indep.bil.df$observed<-as.factor(indep.bil.df$observed)
  xtab<-table(indep.bil.df$predicted,indep.bil.df$observed)
  return(xtab)
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


testeval.eu.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_eu_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.eu.bin.rast

binary_be_rasters<-sapply(names(thresholds), function(x) raster::reclassify(ens_pred_hab_be[[x]],c(0,thresholds[[x]]$predicted,0, thresholds[[x]]$predicted,1,1)),simplify=FALSE)
testeval.be.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_be_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.be.bin.rast
