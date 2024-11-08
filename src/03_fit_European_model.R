
  ##### To do: Specify project and country for which you aim to make predictions
#specify project name
projectname<-"Test_Frédérique"

##### 0. Load packages
options("rgdal_show_exportToProj4_warnings"="none")
library(sf)
library(here)
library(ggplot2)
library(RColorBrewer)
library(terra)
library(dplyr)
library(qs)
library(magick)
library(patchwork)
library(grid)

# Define the desired version
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


#r specify country
#If you'd like to predict for another country, change the shapefile
country<-st_read(here("./data/external/GIS/Belgium/belgium_boundary.shp"))
country_ext<-terra::ext(country) 
country_vector <- vect(country) #Convert country to a SpatVector that can be used for masking


##### 1. Source helper functions
source("./src/helper_functions.R")


#### Load shapefiles that will be used
euboundary<-st_read(here("./data/external/GIS/Europe/EUROPE.shp")) 


#### Create RasterStack of European climate variables from RMI
rmiclimrasters <- list.files((here("./data/external/climate/rmi_corrected")),pattern='tif',full.names = T) 
rmiclimrasters #shows all available climate data
rmiclimpreds<-rast(rmiclimrasters) #includes all available climate data


#### Load species that should be modelled
taxa_info<-read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  pull(accepted_taxonkeys)%>%
  unique()


### Start loop

system.time({
  for(key in accepted_taxonkeys){
    #Extract species name
    species<-taxa_info%>%
      filter(accepted_taxonkeys==key)%>%
      pull(scientificName)%>%
      unique()
    
    #Extract first two words of species name
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
    
    #Extract rest of species name
    rest_of_name <- sub("^\\w+\\s+\\w+\\s+(.*)", "\\1", species)  
    
    #Define taxonkey
    taxonkey<- key
    
    #Read in globalmodels object that was stored as part of  script 02_fit_global_model
    globalmodels<-qread( paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey,"/Global_model_",first_two_words,"_",taxonkey,".qs"
    ))
    
    #Read in different data objects stored in globalmodels
    global.occ.sf<-globalmodels$occurrences
    biasgrid_sub<-terra::unwrap(globalmodels$biasgrid)
    global_model<-terra::unwrap(globalmodels$global_model_predictions)
    model_accuracy<-globalmodels$model_accuracy

    
    ### Create European subset of occurrences
    eu_occ<-st_join(euboundary, global.occ.sf)%>%
      select(decimalLatitude, decimalLongitude, species)%>%
      filter(!is.na(decimalLatitude))%>%
      st_drop_geometry()
    
    # Check if there are any occurrences in Europe
    if (nrow(eu_occ) == 0) {
      warning(paste("No occurrences in Europe for species:", species, 
                    "\n- European model cannot be constructed, skipping to the next species."))
      next  # Skip to the next species in the loop
    }
  
    
    #Plot European subset
    ggplot()+ 
      geom_sf(data = euboundary,  colour = "black", fill = NA)+
      geom_point(data=eu_occ, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
      labs(x="Longitude", y="Latitude")+
      theme_bw()
    
    
    ### Transform eu occurrence dataset with unique presences back to an sf dataframe 

    euocc <- st_as_sf(eu_occ, coords=c("decimalLongitude", "decimalLatitude"), crs= 4326, remove=FALSE)  %>%
      st_transform(crs=st_crs(rmiclimpreds))%>%
      st_coordinates()
  
    
    ### Clip bias grid to European extent
    
    ecoregions_eu<-terra::crop(biasgrid_sub, euboundary)
    biasgrid_eu <- project(ecoregions_eu, rmiclimpreds) #reproject the ecoregions raster to match the spatial properties of rmi climpreds
    
    
    #Plot biasgrid
    ggplot()+ 
      geom_sf(data = euboundary,  colour = "black", fill = NA)+
      geom_spatraster(data=biasgrid_eu)+
      scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange")+
      labs(x="Longitude", y="Latitude")+
      theme_bw()
    
    
    ### Mask areas of high habitat suitability from global climate model
    #Create a mask of the global_model rasterlayer, containing only areas that are predicted to contain occurrences
    m<- global_model >= model_accuracy$threshold
    
    #Mask the global_model layer with this occurrence layer (i.e., convert all areas where occurrences are predicted to NA)
    global_mask<-mask(global_model,m,maskvalue=TRUE)
    global_masked_proj<-project(global_mask,biasgrid_eu)
    
    
    
    
    # Restrict pseudoSamplingArea to environmental layers
    #New: mask pseudoSamplingArea with one of the environmental layers to make sure no pseudoabsences are generated outside the environmental layers
    global_masked_proj<-terra::mask(global_masked_proj, rmiclimpreds[[1]]) 
    
    
    # Plot global_masked_proj
    brks <- seq(0, 1, by=0.1) 
    nb <- length(brks)-1 
    # Generate Viridis palette
    viridis_palette <- viridis(nb)
    
    ggplot()+ 
      geom_spatraster(data=global_masked_proj)+
      geom_sf(data = euboundary,  colour = "black", fill = NA)+
      scale_fill_gradientn(colors = viridis_palette, breaks = brks, labels = brks , na.value = NA) +
      labs(x="Longitude", y="Latitude")+
      theme_bw()
    
    
    ### Combine areas of low predicted habitat suitability with bias grid to exclude low sampled areas and areas of high suitability
    
    #create_pseudoSamplingArea
    pseudoSamplingArea<-mask(global_masked_proj,biasgrid_eu)
    
    
    
    ### Randomly locate pseudo absences within "pseudoSamplingArea" 
    
    # set number of pseudoabsences equal to the number of presences
    numb.eu.pseudoabs<-nrow(euocc)
    
    # takes 10 draws of random pseudoabsences, returns as dataframes and names them X1-X10
    setlist<-seq(1,10,1)
    set.seed(120)
    
    # Generate pseudoabsences 10 times, store in a list with 10 datasets
    pseudoabs_pts <- lapply(setlist, generate_pseudoabs, mask = raster(pseudoSamplingArea), alternative_mask = raster(global_masked_proj), n = numb.eu.pseudoabs, p = euocc)
    names(pseudoabs_pts) <- paste0("X", setlist)
    
    
    
    ###  Prepare occurrence (presence-pseudoabsence) datasets for modelling
    
    # extract data from predictors for absences
    pseudoabs_pts1<-lapply(pseudoabs_pts, function(x) terra::extract(rmiclimpreds,x, ID=FALSE))
    
    # add occ column with value 0 (indicating absences)
    pseudoabs_pts2<-lapply(pseudoabs_pts1, function(x) add.occ(x,0))
    
    # extract eu presences and add presence indicator 
    presence<-as.data.frame(euocc)
    names(presence)<- c("x","y")
    presence1<-terra::extract(rmiclimpreds,presence, ID=FALSE)
    occ<-rep(1,nrow(presence1))
    presence1<-cbind(presence1,occ)
    
    # join each pseudoabsence set with presences 
    eu_presabs.pts<-lapply(pseudoabs_pts2,  function(x) rbind(x, presence1))
    eu_presabs.coord<-lapply(pseudoabs_pts, function(x) rbind(x,presence))
    
    
    
    ### Identify highly correlated climate predictors from training data

    # convert eu data to dataframe
    eu_presabs.pts.df<-lapply(eu_presabs.pts,function(x) as.data.frame(x))
    
    # find attributes that are highly corrected 
    highlyCorrelated_climate <-lapply(names(eu_presabs.pts.df),function(x) findCorrelation(cor(eu_presabs.pts.df[[x]],use = 'complete.obs'), cutoff=0.7,exact=TRUE,names=TRUE))
    
    highlyCorrelated_climate 
    eupreds<-as.data.frame(highlyCorrelated_climate[1])
    kable(eupreds) %>%
      kable_styling(bootstrap_options = c("striped"))
 
    
    ### Remove highly correlated climate predictors from training data
    drop_climate<-highlyCorrelated_climate[[1]]
    keep_layers <- !(names(rmiclimpreds) %in% drop_climate)
    rmiclimpreds_uncor <- subset(rmiclimpreds, keep_layers)
    
    
    ### Add habitat and anthropogenic predictors
    habitat<-list.files((here("./data/external/habitat")),pattern='tif',full.names = T)
    habitat_stack<-rast(habitat[c(1:5,7)]) #Distance to water has another extent
    fullstack<-c(rmiclimpreds_uncor,habitat_stack) #combine uncorrelated climate variable selected earlier with habitat
    
    # clip fullstack to belgium extent (if using another country, replace with country boundary shapefile)
    fullstack_crop<-crop(fullstack,country_ext)
    fullstack_be<-mask(fullstack_crop,country_vector)
    
    occ.full.data <-lapply(eu_presabs.coord, function(x) extract(fullstack,x, ID=FALSE))

    
    ### Identify highly correlated predictors from the habitat/anthropogenic/climate stack (full stack)
    
    # find attributes that are highly correlated
    highlyCorrelated_full <-lapply(names(occ.full.data),function(x)
      findCorrelation(cor(occ.full.data[[x]],use = 'complete.obs'), cutoff=0.7,exact=TRUE,names=TRUE))
    
    highlyCorrelated_vec<-unlist(highlyCorrelated_full)
    eupreds1<-as.data.frame(highlyCorrelated_vec)
    kable(eupreds1) %>%
      kable_styling(bootstrap_options = c("striped"))
    
    ### Remove highly correlated predictors from occ.full.data
    occ.full.data<-sapply(names(occ.full.data),function (x) occ.full.data[[x]][,!(colnames(occ.full.data[[x]]) %in% highlyCorrelated_vec)],simplify=FALSE)

    
    ### Remove highly correlated predictors from fullstack
    keep_layers <- !(names(fullstack) %in% highlyCorrelated_vec)
    fullstack <- subset(fullstack, keep_layers)

    
    ###  Identify and remove near zero variance predictors
    # identify low variance predictors
    nzv_preds<-lapply(names(occ.full.data),function(x) caret::nearZeroVar(occ.full.data[[x]],names=TRUE))
    nzv_preds
    nzv_preds.vec<-unique(unlist(nzv_preds))
    nzv_preds.vec
    # remove near zero variance predictors. They don't contribute to the model.
    occ.full.data<-sapply(names(occ.full.data),function (x) occ.full.data[[x]][,!(colnames(occ.full.data[[x]]) %in% nzv_preds.vec)],simplify=FALSE)
  
    
    
    ### Build models with climate and habitat data
    
    # prepare data for modeling
    
    occ.full.data.df<-lapply(occ.full.data, function(x) as.data.frame(x))
    
    occ.full.data.df<- sapply(names(occ.full.data.df), function (x) cbind(occ.full.data.df[[x]],occ=eu_presabs.pts.df[[x]]$occ, deparse.level=0),simplify=FALSE)
    
    #Recode factor levels of column 'occ' to absent (0) and present(1), and set present as the reference level
    occ.full.data.factor<-sapply(names(occ.full.data.df), function (x) factorVars(occ.full.data.df[[x]], "occ"),simplify=FALSE)
    
    #Replace NA values with 0 in all columns 
    occ.full.data.forCaret<-sapply(names(occ.full.data.factor), function (x) replace(occ.full.data.factor[[x]], is.na(occ.full.data.factor[[x]]),0),simplify=FALSE)
    
    # method = LOOCV (aka "jacknife" ) should be used when occurrences are smaller than n=10 for each predictor in the model)
    if(nrow(presence)<10){
      control<-trainControl(method="LOOCV",
                            savePredictions="final", 
                            preProc=c("center","scale"),
                            classProbs=TRUE) 
    }else{
      control <- trainControl(method="cv",
                              number=4,
                              savePredictions="final", 
                              preProc=c("center","scale"),
                              classProbs=TRUE)
    }
    
    mylist<-list(
      glm =caretModelSpec(method = "glm",maxit=100),
      gbm= caretModelSpec(method = "gbm"),
      rf = caretModelSpec(method = "rf", importance = TRUE),
      earth= caretModelSpec(method = "earth"))
    
    # set.seed(167)
    eu_models<-sapply(names(occ.full.data.forCaret), function(x) model_train_habitat <- caretList(
      occ~., 
      data= occ.full.data.forCaret[[x]],
      trControl=control,
      tuneList=mylist), 
      simplify=FALSE)
    
    
    ### Display model evaluation statistics
    EU_ModelResults1<-sapply(names(eu_models), function(x) resamples(eu_models[[x]]),simplify=FALSE)
    Results.summary<-sapply(names(EU_ModelResults1), function(x) summary(EU_ModelResults1[[x]]),simplify=FALSE)
    Results.summary
  
    
    #show_euModel_correlation
    Model.cor<-sapply(names(eu_models), function(x) modelCor(resamples(eu_models[[x]])),simplify=FALSE)
    Model.cor
 
    
    
    ### Create ensemble model 
    set.seed(458)
    lm_ens_hab<-sapply(names(eu_models), function (x) caretEnsemble(eu_models[[x]], trControl=trainControl(method="cv",                                                               number=10,savePredictions= "final",classProbs = TRUE)),simplify=FALSE)

    
    
    ###  Evaluate the performance of each the EU level ensemble models based on results from CV
    # identify threshold where sensitivity=specifity
    thresholds<-sapply(names(lm_ens_hab), function(x) findThresh(lm_ens_hab[[x]]$ens_model$pred),simplify=FALSE)
  
    
    ### Using thresholds identified for each model in the previous step, assess performance of each model
    # accuracy measures
    thresholds.df<-sapply(names(thresholds), function(x) accuracyStats(lm_ens_hab[[x]]$ens_model$pred,thresholds[[x]]$predicted),simplify=FALSE)
    thresholds.comb<-do.call(rbind,thresholds.df)
    kable(thresholds.comb,digits=2)
    
    
    ### Select the best model
    # specify best model below based on highes PCC, and, in case there are multiple rows with the same PCC, the highest AUC
    bestmodelname <- thresholds.comb %>%
      filter(PCC == max(PCC)) %>%      
      slice_max(AUC, n = 1)%>%
      rownames()
    
    bestModel<-lm_ens_hab[[bestmodelname]]
  
    
    ###  Use best EU level ensemble models to predict at European level 

    system.time({
      ens_pred_hab_eu1<-terra::predict(fullstack,lm_ens_hab[[bestmodelname]],type="prob", na.rm = TRUE)
    }) 
  
    
    
    #Create sf with occurrence points
    euocc1<-st_as_sf(as.data.frame(euocc), coords=c("X","Y"),crs=st_crs(rmiclimpreds))
  
    
    ### plot the best EU level ensemble model
    
    #Plot
    brks <- seq(0, 1, by=0.1)
    nb <- length(brks) - 1
    viridis_palette <- viridis(nb)
    
    eu_plot<-ggplot() + 
      geom_spatraster(data = ens_pred_hab_eu1) +
      scale_fill_gradientn(colors = viridis_palette, 
                           breaks = brks, 
                           labels = brks, 
                           na.value = NA) +
      geom_sf(data = euocc1, color = "black", fill = "red", 
              size = 1.5, shape = 21) +
      theme_bw() +
      labs(fill = "Suitability")+
      coord_sf(xlim = c(2254476, 6005897), 
               ylim = c(1363659, 5469923))
    
    #Create an empty plot to fill PDF
    empty_plot <- ggplot() + 
      theme_void() + 
      theme(plot.background = element_blank()) 
    
    #Create final plot
    plot_final<-eu_plot /empty_plot 

    
    ### Export European raster and PDF
    #Create folders
    raster_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Rasters")
    PDF_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PDFs")
    dir.create(raster_folder, recursive = TRUE, showWarnings = FALSE)
    dir.create(PDF_folder, recursive = TRUE, showWarnings = FALSE)
    
    #Export raster
    writeRaster(ens_pred_hab_eu1,
                filename=file.path(raster_folder,paste(first_two_words,"_",taxonkey,"_hist_EU.tif",sep="")),
                overwrite=TRUE)
    
    #Export PDF 
    #Define the file paths
    plot_png_path <- file.path(PDF_folder,paste(first_two_words,"_",taxonkey,"_hist_EU.png",sep=""))
    plot_pdf_path <- file.path(PDF_folder,paste(first_two_words,"_",taxonkey,"_hist_EU.pdf",sep=""))
    
    # Save each plot as a PDF file
    ggsave(filename = paste0(first_two_words,"_",taxonkey,"_hist_EU.png"), plot = plot_final, 
           device = "png", width =8.27 , height = 11.69, path= PDF_folder)
    
    # Read the PNG image back in
    img <- image_read(plot_png_path)
    
    # Start a PDF device for output
    pdf(plot_pdf_path, width = 8.27, height = 11.69)
    
    # Create a layout for title and image
    grid.newpage()
    
    # Add title at the top of the PDF
    grid.text(
      label = bquote(italic(.(first_two_words)) ~ .(rest_of_name) ~ "( " * .(taxonkey) * ")"),
      x = 0.5, y = 0.95, just = "center", gp = gpar(fontsize = 12, fontface = "bold")
    )
    
    # Add the PNG image below the title
    grid.raster(img, width = unit(0.9, "npc"), height = unit(0.9, "npc"), y = 0.47)
    
    # Close the PDF device
    while (dev.cur() > 1) dev.off()
    
    # Remove the PNG file from the local directory
    file.remove(plot_png_path)
    
    
    ### Use best EU level ensemble model to predict for Belgium only  
    
    # creates  country level rasters using the European level models
    system.time({
      ens_pred_hab_be<-terra::predict(fullstack_be,lm_ens_hab[[bestmodelname]],type="prob", na.rm=TRUE)
    })
    
    
    # export rasters as GeoTiffs
    #lapply(names(ens_pred_hab_be), function(x) writeRaster(ens_pred_hab_be[[x]], filename=file.path(rasterOutput,paste(x,"_",taxonkey,"_hist.tif",sep="")),  format="GTiff",overwrite=TRUE))

    
    ### Plot Belgian model
    brks <- seq(0, 1, by=0.1)
    nb <- length(brks) - 1
    viridis_palette <- viridis(nb)
    
    ggplot() + 
      geom_spatraster(data = ens_pred_hab_be) +
      scale_fill_gradientn(colors = viridis_palette, 
                           breaks = brks, 
                           labels = brks, 
                           na.value = NA) +
      theme_bw() +
      labs(fill = "Suitability")
    
    
    ### Close loop

    eumodel <-list(species = species,
                   taxonkey = taxonkey,
                   euocc1 = euocc1,
                   bestModel=bestModel,
                   fullstack_be=terra::wrap(fullstack_be)
    )
    
    qsave(eumodel, paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey,"/EU_model_",first_two_words,"_",taxonkey,".qs"))
    
    print(paste("European model has been created for", species))
    
  }
})


# Clean R environment
rm(list = ls())

