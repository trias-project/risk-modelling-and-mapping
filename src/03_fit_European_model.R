#--------------------------------------------
#-----------To do: specify project ----------
#--------------------------------------------
#specify project name
projectname<-"Project_Frédérique"


#--------------------------------------------
#-----------    Load packages      ----------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")

packages <- c( "viridis","dplyr", "here", "qs","terra", "tidyterra","sf", "ggplot2","RColorBrewer","magick","patchwork","grid", "randomForest", "progressr", "raster", "dismo", "caret", "caretEnsemble", "kableExtra","gbm", "PresenceAbsence")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#---  Load right version of caretEnsemble  --
#--------------------------------------------
# Define the required version
desired_version <- "2.0.3"

if ("caretEnsemble" %in% rownames(installed.packages())) {
  current_version <- packageVersion("caretEnsemble")
  if (as.character(current_version) != desired_version) {
    remove.packages("caretEnsemble")
    devtools::install_github("zachmayer/caretEnsemble@2.0.3")
    library(caretEnsemble)
    rm(current_version, desired_version)
  } else {
    library(caretEnsemble)
    rm(current_version, desired_version)
  }
  
} else {
  devtools::install_github("zachmayer/caretEnsemble@2.0.3")
  library(caretEnsemble)
  rm(desired_version)
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
#------- Source helper fucntions     --------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#---------   Load shape of Europe   ---------
#--------------------------------------------
euboundary<-sf::st_read(here("./data/external/GIS/Europe/EUROPE.shp")) 


#--------------------------------------------
#-------- Load European habitat rasters -----
#--------------------------------------------
habitat<-list.files((here("./data/external/habitat")),pattern='tif',full.names = T)
habitat_stack<-terra::rast(habitat[c(1:5,7)]) #Distance to water (layer 6) has another extent and we're not sure if this layer is correct: leave it out


#--------------------------------------------
#-------- Load European climate rasters -----
#--------------------------------------------
rmiclimrasters <- list.files((here("./data/external/climate/rmi_corrected")),pattern='tif',full.names = T) 
rmiclimrasters 
rmiclimpreds<-terra::rast(rmiclimrasters) 
rmiclimpreds<-terra::crop(rmiclimpreds, ext(habitat_stack))


#---------------------------------------------
#----- Remove NA pixels from predictors ------
#---------------------------------------------
#This is to avoid that some layers have NA while others have values in certain pixels
#First mask pixels in the rasterstack where at least one layer has NA
na_mask_rmiclimpreds <- app(rmiclimpreds, function(x) any(is.na(x)))
na_mask_habitat_stack<- app(habitat_stack, function(x) any(is.na(x)))
rmiclimpreds<- terra::mask(rmiclimpreds, na_mask_rmiclimpreds, maskvalue=1)
habitat_stack<- terra::mask(habitat_stack, na_mask_habitat_stack, maskvalue=1)

#Second mask rmiclimpreds with habitat_stack and vice versa
rmiclimpreds<-terra::mask(rmiclimpreds, habitat_stack[[1]])
habitat_stack<-terra::mask(habitat_stack, rmiclimpreds[[1]])


#--------------------------------------------
#------------- Load species data -----------
#--------------------------------------------
taxa_info<-read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  dplyr::pull(speciesKey)%>%
  unique()


#--------------------------------------------
#----------- Start modelling loop  ----------
#--------------------------------------------

with_progress({
    p <- progressr::progressor(along = 1:length(accepted_taxonkeys)) 
  for(key in accepted_taxonkeys){ #Approx. 13 min per species
    
    #--------------------------------------------
    #---------------Map progress  ---------------
    #--------------------------------------------
    p()
    
    #--------------------------------------------
    #--------Extract species-specific data  -----
    #--------------------------------------------
    #Extract species name
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
    #-- Define file path of global model file  --
    #--------------------------------------------  
    global_model_file<-file.path("./data/projects",projectname, paste0(first_two_words,"_",taxonkey),paste0("Global_model_",first_two_words,"_",taxonkey,".qs"))
   
   
    #--------------------------------------------
    #-Check if global model exists, if not, skip-
    #--------------------------------------------
    if(file.exists(global_model_file)){
      
      #This was stored as part of  script 02_fit_global_model
      globalmodels<-qs::qread( paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey,"/Global_model_",first_two_words,"_",taxonkey,".qs"))
      
      #Extract different data objects stored in globalmodels
      global.occ.sf<-globalmodels$occurrences
      model_accuracy<-globalmodels$model_accuracy
      
    }else{
      warning(paste0("Skipping species ", species, " because no global model could be fitted"))
      next  # Skip the rest of the loop and move to the next iteration
    }
    
    
    #--------------------------------------------
    #------------ Import raster layers ----------
    #--------------------------------------------
    #Define file paths
    biasgrid_file<- file.path("./data/projects",projectname,paste0(first_two_words,"_",taxonkey),"Rasters","Interim",paste0("Biasgrid_",first_two_words,"_",taxonkey,".tif"))
    global_model_file<- file.path("./data/projects",projectname,paste0(first_two_words,"_",taxonkey),"Rasters","Global",paste0("Global_model_",first_two_words,"_",taxonkey,".tif"))
    
    #Load rasterlayers
    biasgrid_sub<-terra::rast(biasgrid_file)
    global_model<-terra::rast(global_model_file)
    
    
    #--------------------------------------------
    #------------ Define folder paths -----------
    #--------------------------------------------
    raster_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "Rasters")
    raster_EU_folder <- file.path(raster_folder, "Europe")
    PDF_folder <- file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PDFs")
    PDF_EU_folder <- file.path(PDF_folder, "Europe")
    PNG_folder<-file.path("./data/projects", projectname, paste0(first_two_words, "_", taxonkey), "PNGs")
    PNG_EU_folder <- file.path(PNG_folder, "Europe")
    
    
    #-------------------------------------------------
    #--------------- Create EU folders ---------------
    #-------------------------------------------------
    # Define the folder paths
    folder_paths<-list(list("path"= raster_EU_folder,
                            "name"= "Rasters/Europe"),
                       list("path"= PDF_EU_folder,
                            "name"= "PDF/Europe"),
                       list("path"= PNG_EU_folder,
                            "name"= "PNG/Europe"))
    
    # Check and create each folder if necessary
    lapply(folder_paths, function(folder){
      create_folder(folder$path, folder$name)
    })
    
    
    #-----------------------------------------------
    #----- Create subset of European records -------
    #-----------------------------------------------
    #Check for occurrences that fall within Europe
    eu_occ <- global.occ.sf[sf::st_intersects(global.occ.sf, euboundary, sparse = FALSE), ] %>%
      dplyr::select(decimalLatitude, decimalLongitude, species) %>%
      dplyr::filter(!is.na(decimalLatitude))%>%
      sf::st_transform(crs=st_crs(rmiclimpreds)) 
    
    # Convert to crs of rmiclimpreds
    eu_occ<-eu_occ%>%
      sf::st_coordinates()%>%
      cbind(., eu_occ)%>%
      dplyr::select(-c(decimalLatitude, decimalLongitude))
    
    #Only keep occurrences in pixels that have predictor data (not NA's)
    extracted_value <- terra::extract(rmiclimpreds[[1]], vect(eu_occ))
    eu_occ$extracted_value<-extracted_value[,2]
    eu_occ <- eu_occ[!is.na(eu_occ$extracted_value), ]
    
    # Keep XY coordinates
    euocc<-eu_occ%>%
      sf::st_coordinates()
    
    
    #------------------------------------------------
    #----- Check if at least 20 European records ----
    #------------------------------------------------
    if (nrow(euocc) < 20) {
      warning(paste(nrow(euocc)," occurrences in Europe for species:", species, 
                    "\n- European model cannot be constructed, skipping to the next species."))
      next  # Skip to the next species in the loop
    }
  
    
    #--------------------------------------------
    #-------- Plot European occurrences ---------
    #--------------------------------------------
    #ggplot()+ 
     # geom_sf(data = euboundary,  colour = "black", fill = NA)+
      #geom_point(data=eu_occ, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
      #labs(x="Longitude", y="Latitude")+
      #theme_bw()
  
    
    #--------------------------------------------
    #----- Clip biasgrid to European extent -----
    #--------------------------------------------
    ecoregions_eu<-terra::crop(biasgrid_sub, euboundary)
    biasgrid_eu <- terra::project(ecoregions_eu, rmiclimpreds) #reproject the ecoregions raster to match the spatial properties of rmi climpreds
    
    
    #--------------------------------------------
    #----------- Visualize biasgrid  ------------
    #--------------------------------------------
    #ggplot()+ 
     # geom_sf(data = euboundary,  colour = "black", fill = NA)+
      #geom_spatraster(data=biasgrid_eu)+
      #scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange")+
      #labs(x="Longitude", y="Latitude")+
      #theme_bw()
    
    
    #----------------------------------------------
    #------- Mask areas of high suitability -------
    #----------------------------------------------
    #Create a mask of the global_model rasterlayer, containing only areas that are predicted to contain occurrences
    m<- global_model >= model_accuracy$threshold
    
    #Mask the global_model layer with this occurrence layer (i.e., only keep pixels where absences are predicted, rest becomes NA)
    global_mask<-terra::mask(global_model,m,maskvalue=TRUE)
    global_masked_proj<-terra::project(global_mask,biasgrid_eu)
    
    #New: mask with one of the environmental layers to make sure no pseudoabsences are generated outside the environmental layers
    global_masked_proj<-terra::mask(global_masked_proj, rmiclimpreds[[1]]) 
    
    
    #--------------------------------------------
    #--Create sampling area for pseudoabsences --
    #--------------------------------------------
    # Combine areas of low predicted habitat suitability with bias grid to exclude pixels with no sampling effort or falling outside ecoregions with global occurrences! (NA in biasgrid_eu)
    pseudoSamplingArea<-terra::mask(global_masked_proj,biasgrid_eu)
    
    
    #--------------------------------------------
    #-------- Plot pseudosampling area ----------
    #--------------------------------------------
    #brks <- seq(0, 1, by=0.1) 
    #nb <- length(brks)-1 
    # Generate Viridis palette
    #viridis_palette <- viridis(nb)
    
    #ggplot()+ 
     # geom_spatraster(data=global_masked_proj)+
      #geom_sf(data = euboundary,  colour = "black", fill = NA)+
      #scale_fill_gradientn(colors = viridis_palette, breaks = brks, labels = brks , na.value = NA) +
      #labs(x="Longitude", y="Latitude")+
      #theme_bw()
    
    
    #--------------------------------------------
    #Randomly generate pseudoabsences in pseudoSamplingArea
    #--------------------------------------------
    # Set number of pseudoabsences equal to the number of presences
    numb.eu.pseudoabs<-nrow(euocc)
    
    # Generate pseudoabsences 10 times, store in a list with 10 datasets and names them X1-X10
    setlist<-seq(1,10,1)
    set.seed(120)
    pseudoabs_pts <- lapply(setlist, generate_pseudoabs, mask = raster(pseudoSamplingArea), alternative_mask = raster(global_masked_proj), n = numb.eu.pseudoabs, p = euocc)
    names(pseudoabs_pts) <- paste0("X", setlist)
    
    
    #--------------------------------------------
    #Prepare presence-absence dataset for modelling
    #--------------------------------------------
    # extract data from environmental predictors for each list of absences
    pseudoabs_pts1<-lapply(pseudoabs_pts, function(x) terra::extract(rmiclimpreds,x, ID=FALSE))
    
    # add occ column with value 0 (indicating absences)
    pseudoabs_pts2<-lapply(pseudoabs_pts1, function(x) add.occ(x,0))
    
    # extract environmental data for eu presences and add presence indicator (1)
    presence<-as.data.frame(euocc)
    names(presence)<- c("x","y")
    presence1<-terra::extract(rmiclimpreds,presence, ID=FALSE)
    presence1$occ<-1
    
    # join each pseudoabsence set with presences 
    eu_presabs.pts<-lapply(pseudoabs_pts2,  function(x) rbind(x, presence1))
    eu_presabs.coord<-lapply(pseudoabs_pts, function(x) rbind(x,presence))
    
    
    #--------------------------------------------
    #--Remove highly correlated predictors from training data --
    #--------------------------------------------
    # convert eu data to dataframe
    eu_presabs.pts.df<-lapply(eu_presabs.pts,function(x) as.data.frame(x))
    
    # find attributes that are highly corrected 
    highlyCorrelated_climate <-lapply(eu_presabs.pts.df, function(df) as.data.frame(cor(df[, 1:13], use = "complete.obs")))
    
    #Calculate the mean correlation over the 10 datsets and identify highly correlated variables
    mean_correlation_matrix <- Reduce("+", highlyCorrelated_climate) / length(highlyCorrelated_climate)
    drop_climate<-caret::findCorrelation(as.matrix(mean_correlation_matrix), cutoff=0.7,exact=TRUE,names=TRUE)
    
    #Only keep layers that are not highly correlated
    rmiclimpreds_uncor <- subset(rmiclimpreds, !(names(rmiclimpreds) %in% drop_climate))
    
    
    #--------------------------------------------
    #- Add habitat and anthropogenic predictors -
    #--------------------------------------------
    #combine uncorrelated climate variable selected earlier with habitat layers
    fullstack<-c(rmiclimpreds_uncor,habitat_stack) 
    
    
    #-----------------------------------------------------------
    #- Extract predictor values for presences and pseudoabsences
    #-----------------------------------------------------------
    occ.full.data <-lapply(eu_presabs.coord, function(x) terra::extract(fullstack,x, ID=FALSE))

    
    #--------------------------------------------
    #--- Remove highly correlated predictors ----
    #--------------------------------------------
    # find attributes that are highly correlated in at least one of the models and remove them from all
    highlyCorrelated_full <-lapply(names(occ.full.data),function(x) caret::findCorrelation(cor(occ.full.data[[x]],use = 'complete.obs'), cutoff=0.7,exact=TRUE,names=TRUE))
    highlyCorrelated_vec<-unique(unlist(highlyCorrelated_full))
    
    #highlyCorrelated_full <-lapply(occ.full.data, function(df) as.data.frame(cor(df, use = "complete.obs")))
    #Calculate the mean correlation over the 10 datsets and identify highly correlated variables
    #mean_correlation_full_matrix <- Reduce("+", highlyCorrelated_full) / length(highlyCorrelated_full)
    #highlyCorrelated_vec<-findCorrelation(as.matrix(mean_correlation_full_matrix), cutoff=0.7,exact=TRUE,names=TRUE)
    
    # Remove highly correlated predictors from dataset holding occurrences
    occ.full.data<-sapply(names(occ.full.data),function (x) occ.full.data[[x]][,!(colnames(occ.full.data[[x]]) %in% highlyCorrelated_vec)],simplify=FALSE)
    
    # Remove highly correlated predictors from rasterlayers
    fullstack <- subset(fullstack, !(names(fullstack) %in% highlyCorrelated_vec))

    
    #--------------------------------------------
    #--- Remove near-zero variance predictors ---
    #--------------------------------------------
    # identify low variance predictors
    nzv_preds<-lapply(names(occ.full.data),function(x) caret::nearZeroVar(occ.full.data[[x]],names=TRUE))
    nzv_preds.vec<-unique(unlist(nzv_preds))
    
    # remove near zero variance predictors. They don't contribute to the model.
    occ.full.data<-sapply(names(occ.full.data),function (x) occ.full.data[[x]][,!(colnames(occ.full.data[[x]]) %in% nzv_preds.vec)],simplify=FALSE)
  
    #remove them from fullstack
    fullstack <- fullstack%>%
      subset(!names(fullstack) %in% nzv_preds.vec)
    
    
    #--------------------------------------------
    #-------- Prepare data for modelling --------
    #--------------------------------------------
    #Convert to dataframe
    occ.full.data.df<-lapply(occ.full.data, function(x) as.data.frame(x))
    
    #Add column with occurrence data (occ)
    occ.full.data.df<- sapply(names(occ.full.data.df), function (x) cbind(occ.full.data.df[[x]],occ=eu_presabs.pts.df[[x]]$occ, deparse.level=0),simplify=FALSE)
    
    #Recode factor levels of column 'occ' to absent (0) and present(1), and set present as the reference level
    occ.full.data.forCaret<-sapply(names(occ.full.data.df), function (x) factorVars(occ.full.data.df[[x]], "occ"),simplify=FALSE)
    
    
    #--------------------------------------------
    #- Run models with climate and habitat data -
    #--------------------------------------------
    control <- caret::trainControl(method="cv",
                              number=4,
                              savePredictions="final", 
                              preProc=c("center","scale"),
                              classProbs=TRUE)
   
    mylist<-list(
      glm =caretEnsemble::caretModelSpec(method = "glm",maxit=100),
      gbm= caretEnsemble::caretModelSpec(method = "gbm"),
      rf = caretEnsemble::caretModelSpec(method = "rf", importance = TRUE),
      earth= caretEnsemble::caretModelSpec(method = "earth"))
    
    set.seed(167)
    eu_models<-sapply(names(occ.full.data.forCaret), function(x) model_train_habitat <- caretEnsemble::caretList(
      occ~., 
      data= occ.full.data.forCaret[[x]],
      trControl=control,
      tuneList=mylist), 
      simplify=FALSE)
    
    
    #--------------------------------------------
    #---- Display model evaluation statistics----
    #--------------------------------------------
    EU_ModelResults1<-sapply(names(eu_models), function(x) caret::resamples(eu_models[[x]]),simplify=FALSE)
    Results.summary<-sapply(names(EU_ModelResults1), function(x) summary(EU_ModelResults1[[x]]),simplify=FALSE)
  
    #show_euModel_correlation
    Model.cor<-sapply(names(eu_models), function(x) caret::modelCor(resamples(eu_models[[x]])),simplify=FALSE)
 
    
    #--------------------------------------------
    #---------- Create ensemble model -----------
    #--------------------------------------------
    #9 folds (e.g., 90 records) will be used for training, and 1 fold (e.g., 10 records) will be used for validation.
    control <- caret::trainControl(method="cv",
                              number=10,
                              savePredictions="final", 
                              classProbs=TRUE)
    
    set.seed(458)
    lm_ens_hab<-sapply(names(eu_models), function (x) caretEnsemble::caretEnsemble(eu_models[[x]],trControl= control), simplify=FALSE)
    
    
    #--------------------------------------------
    #- Evaluate each ensemble model's performance -
    #-------------------------------------------- 
    #based on results from CV, 
    #identify threshold where sensitivity=specifity
    thresholds<-sapply(names(lm_ens_hab), function(x) findThresh(lm_ens_hab[[x]]$ens_model$pred),simplify=FALSE)
  
    #Using thresholds identified for each model in the previous step, assess performance of each model
    # accuracy measures
    thresholds.df<-sapply(names(thresholds), function(x) accuracyStats(lm_ens_hab[[x]]$ens_model$pred,thresholds[[x]]$predicted),simplify=FALSE)
    thresholds.comb<-do.call(rbind,thresholds.df)
    
    
    #--------------------------------------------
    #---------Select best ensemble model --------
    #--------------------------------------------
    # select best model based on highest PCC, and, in case there are multiple rows with the same PCC, the highest AUC
    bestmodelname <- thresholds.comb %>%
      dplyr::filter(PCC == max(PCC)) %>%      
      dplyr::slice_max(AUC, n = 1)%>%
      rownames()
    
    bestModel<-lm_ens_hab[[bestmodelname]]
  
    #Get performance of best model
    model_performance<-thresholds.comb%>%
      dplyr::filter(rownames(.) == bestmodelname)
    
    #Get correlation between separate models of best model
    Model.cor<-Model.cor[[bestmodelname]]
     
   
    #--------------------------------------------
    #--Store variable importance of best model --
    #--------------------------------------------
    variableImportance_bestmodel<-caret::varImp(bestModel)
    
    
    #--------------------------------------------
    #-Create European predictions using best model -
    #--------------------------------------------
    ens_pred_hab_eu1<-terra::predict(fullstack,lm_ens_hab[[bestmodelname]],type="prob", na.rm = TRUE)
    
  
    #-----------------------------------------------------------------------------
    #---------- Create a confidence map of the best model at EU level ------------
    #-----------------------------------------------------------------------------
    set.seed(792)  
    pvalsdf_hist_eu<-classConformalPrediction(bestModel,ens_pred_hab_eu1)
    hist.conf.map.eu<-confidenceMaps(x=pvalsdf_hist_eu,original_raster=ens_pred_hab_eu1, taxonKey=taxonkey,taxonName=first_two_words,taxonNameTitle=species_title, nameExtension=rest_of_name, scenario="hist", regionName="Europe", folder=raster_EU_folder)
    
    
    #--------------------------------------------
    #- Create sf df with occurrences for plotting -
    #--------------------------------------------
    euocc1<-sf::st_as_sf(as.data.frame(euocc), coords=c("X","Y"),crs=st_crs(rmiclimpreds))
  
    
    #--------------------------------------------
    #- Export Eu predictions as raster and PDF --
    #-------------------------------------------
    #---------------Export raster-------------
    terra::writeRaster(ens_pred_hab_eu1,
                filename=file.path(raster_EU_folder,paste(first_two_words,"_",taxonkey,"_hist_EU.tif",sep="")),
                overwrite=TRUE)
    
    #---------------Export PDF----------------
    exportPDF(predictions=ens_pred_hab_eu1,
              taxonName=first_two_words,
              nameExtension=rest_of_name,
              dataType="",
              taxonNameTitle=species_title,
              taxonKey=taxonkey,
              scenario="hist",
              regionName="EU",
              returnPredictions=FALSE,
              returnPNG=FALSE)
    
    
    #--------------------------------------------
    #- Store datasets linked to best Model
    #--------------------------------------------
    eu_presabs.coord<-eu_presabs.coord[[bestmodelname]]
    occ.full.data.forCaret<-occ.full.data.forCaret[[bestmodelname]]
    
    
    #--------------------------------------------
    #- Store layers used to create best model for country
    #--------------------------------------------
    # create a rasterstack for Belgium in this case
    fullstack_crop<-terra::crop(fullstack,country_ext)
    fullstack_be<-terra::mask(fullstack_crop,country_vector)
    
    
    #--------------------------------------------
    #- Save best model, european occurrences, and layers for Belgium -
    #--------------------------------------------
    eumodel <-list(species = species,
                   taxonkey = taxonkey, 
                   euocc1 = euocc1, #sf dataset with coordinates of presences (geometry format)
                   bestModel=bestModel, #Best ensemble model
                   variable_importance = variableImportance_bestmodel, #Variable importance in each separate model and overall
                   occ.full.data.forCaret = occ.full.data.forCaret, #Data used to fit best model
                   eu_presabs.coord = eu_presabs.coord, #XY coordinates of presences and pseudoabsences
                   model_performance = model_performance, #Performance of best model
                   model_correlation = Model.cor #Correlation between separate models underlying best ensemble model
    )
   
    #Save eumodel as .qs file
    qs::qsave(eumodel, paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey,"/EU_model_",first_two_words,"_",taxonkey,".qs")) 
    
    
    #--------------------------------------------
    #- ------ Save rasters-----------------------
    #--------------------------------------------
    #Define raster path
    fullstack_be_file<- file.path(raster_folder,"Interim", paste0("Fullstack_be_",first_two_words,"_",taxonkey,".tif"))
    
    #Save locally because if in .qs file parts of the metadata will be stored in a Temp file that will be erased over time
    terra::writeRaster(fullstack_be, filename =fullstack_be_file, overwrite = TRUE)

    
    #--------------------------------------------
    #-------- End of loop -----------------------
    #--------------------------------------------
    print(paste("European model has been created for", species))
    rm(list = setdiff(ls(), c("p", "factorVars","accuracyStats", "findThresh", "projectname", "generate_pseudoabs", "create_folder", "country_name", "country_ext", "country_vector", "euboundary", "habitat_stack", "rmiclimpreds", "accepted_taxonkeys", "taxa_info", "key", "add.occ", "confidenceMaps", "classConformalPrediction","extractVals","GetLength","get.confidence", "exportPDF")))
  }
})


#--------------------------------------------
#----------- Clean R environment ------------
#--------------------------------------------
rm(list = ls())

