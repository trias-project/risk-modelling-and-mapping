#--------------------------------------------
#-----------To do: specify project ----------
#--------------------------------------------
#specify project name
projectname<-"Test_Frédérique"


#--------------------------------------------
#-------------- Load packages ---------------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")

packages <- c( "dplyr", "stringr", "here", "qs","CoordinateCleaner","terra", "raster", "sf", "rnaturalearth", 
               "ggplot2","tidyterra","mapview", "dismo", "sdm", "caret", "viridisLite", "kableExtra","future", "future.apply"
               )

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#--- Load right version of caretEnsemble-----
#--------------------------------------------
# Define the required version
desired_version <- "2.0.3"

if ("caretEnsemble" %in% rownames(installed.packages())) {
  current_version <- packageVersion("caretEnsemble")
  if (as.character(current_version) != desired_version) {
    remove.packages("caretEnsemble")
    devtools::install_github("zachmayer/caretEnsemble@2.0.3")
    library(caretEnsemble)
  } else {
    library(caretEnsemble)
  }
  
} else {
  devtools::install_github("zachmayer/caretEnsemble@2.0.3")
  library(caretEnsemble)
  rm(current_version, desired_version)
}


#--------------------------------------------
#-------- Source helper functions------------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#--- Load global occurrences and taxa info---
#--------------------------------------------
global<-qread( paste0("./data/projects/",projectname,"/",projectname,"_occurrences.qs"))
taxa_info<-read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  pull(speciesKey)%>%
  unique()


#--------------------------------------------
#-------- Filter global occurrences----------
#--------------------------------------------
#remove unverified records
identificationVerificationStatus_to_discard <- c( "unverified",
                                                  "unvalidated",
                                                  "not validated",
                                                  "under validation",
                                                  "not able to validate",
                                                  "control could not be conclusive due to insufficient knowledge",
                                                  "1",
                                                  "uncertain",
                                                  "unconfirmed",
                                                  "Douteux",
                                                  "Invalide",
                                                  "Non r\u00E9alisable",
                                                  "verification needed" ,
                                                  "Probable",
                                                  "unconfirmed - not reviewed",
                                                  "validation requested")

#enter value for max coordinate uncertainty in meters, default = 1000
global.occ<-global %>%
  dplyr::filter(speciesKey%in%accepted_taxonkeys) %>%   
  dplyr::filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters<= 1000) %>%
  dplyr::filter(!str_to_lower(identificationVerificationStatus) %in% identificationVerificationStatus_to_discard)

#Remove coordinates that for both lon and lat values, have less than 4 decimal places
global.occ$lon_dplaces<-sapply(global.occ$decimalLongitude, function(x) decimalplaces(x))
global.occ$lat_dplaces<-sapply(global.occ$decimalLatitude, function(x) decimalplaces(x))
global.occ[global.occ$lon_dplaces < 4 & global.occ$lat_dplaces < 4 , ]<-NA
global.occ<-global.occ[ which(!is.na(global.occ$lon_dplaces)),]
global.occ<-within(global.occ,rm("lon_dplaces","lat_dplaces")) # n= 1758


#--------------------------------------------
#------------ Define species group-----------
#--------------------------------------------
global.occ <- global.occ%>%
  dplyr::mutate(Group = case_when(kingdom == "Plantae" ~ "Plants",
                                  class == "Aves" ~ "Birds",
                                  phylum == "Mollusca" ~ "Molluscs",
                                  class == "Amphibia" ~ "Amphibians",
                                  class == "Mammalia" ~ "Mammals",
                                  class == "Crocodylia" ~ "Reptiles",
                                  class == "Testudines" ~ "Reptiles",
                                  class == "Sphenodontia" ~ "Reptiles",
                                  class == "Squamata" ~ "Reptiles",
                                  TRUE ~ NA_character_))


#--------------------------------------------
#-------Prepare occurrence dataset-----------
#--------------------------------------------
global.occ.LL<-data.frame(global.occ)[c(4,3,2,1,10)] #decimalLon, decimalLat, species, acceptedtaxonkey, Group, n= 1758
rm(global.occ, global)


#--------------------------------------------
#-----------Do coordinate cleaning-----------
#--------------------------------------------
# OPTIONAL: Coordinates are tested for several things: whether they are in capitals, whether ... . For each coordinate a column per test is added indicating wether the result is potentially problematic (FALSE) or a clean coordinate (TRUE)
#flags_report<-clean_coordinates(x = global.occ.LL, lon= "decimalLongitude", lat= "decimalLatitude",
#  tests = c("capitals", 
# "centroids","gbif", "institutions", 
# "seas", "zeros"))

# Only keep coordinates that are not flagged as potentially problematic
cleaned<-clean_coordinates(x = global.occ.LL, lon= "decimalLongitude", lat= "decimalLatitude",
                           tests = c("capitals", 
                                     "centroids","gbif", "institutions", 
                                     "zeros"),value="clean")


#--------------------------------------------
#--------Load global climate rasters --------
#--------------------------------------------
globalclimrasters <- list.files((here("./data/external/climate/trias_CHELSA")),pattern='tif',full.names = T) #import CHELSA data
globalclimpreds_terra <- terra::rast(globalclimrasters)


#--------------------------------------------
#--------Load European climate rasters-------
#--------------------------------------------
euclimrasters <- list.files((here("./data/external/climate/chelsa_eu_clips")),pattern='tif',full.names = T)
eu_climpreds<-rast(euclimrasters)
eu_climpreds.10<-divide10(eu_climpreds)  # correct for integer format of Chelsa preds


#--------------------------------------------
#--------- Load shape of the world ----------
#--------------------------------------------
world<-ne_countries(scale=50)


#--------------------------------------------
#-------Load file paths to bias grids -------
#--------------------------------------------
bias_grid_paths <- list(
  Plants = here("./data/external/bias_grids/final/trias/plants_1deg_min5.tif"),
  Amphibians = here("./data/external/bias_grids/final/trias/amphib_1deg_min5.tif"),
  Birds = here("./data/external/bias_grids/final/trias/birds_1deg_min5.tif"),
  Mammals = here("./data/external/bias_grids/final/trias/mammals_1deg_min5.tif"),
  Molluscs = here("./data/external/bias_grids/final/trias/molluscs_1deg_min5.tif"),
  Reptiles = here("./data/external/bias_grids/final/trias/reptiles_1deg_min5.tif")
)


#--------------------------------------------
#------- Split dataframe by taxonkey --------
#--------------------------------------------
sort(unique(cleaned$species))
split_df<-split(cleaned,cleaned$species) 


#--------------------------------------------
#---------------- Clean up ------------------
#--------------------------------------------
rm(global.occ.LL,cleaned)
gc()


#--------------------------------------------
#-------Start loop for SDM modelling --------
#--------------------------------------------
system.time({ # 5 species (43 min)
  for(i in seq_along (split_df)){ 
    globalmodels<-list()
    species<-names(split_df)[i]
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)  # Extract first two words of species name
    global.occ.LL.cleaned<-split_df[[i]]
    taxonkey<-unique(global.occ.LL.cleaned$speciesKey)
    speciesgroup<-unique(global.occ.LL.cleaned$Group)
    global.occ.LL.cleaned<-global.occ.LL.cleaned %>%
      dplyr::select(c(decimalLongitude,decimalLatitude))
    
    
    #--------------------------------------------
    #--------- Create species folder ------------
    #--------------------------------------------
    species_folder<-paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey)
    if (!dir.exists(species_folder)) {
      dir.create(species_folder, recursive = TRUE)
    }
    
    #--------------------------------------------
    #------ Remove duplicates per grid cell -----
    #--------------------------------------------
    global.occ.LL.cleaned$cell<-terra::cellFromXY( globalclimpreds_terra, global.occ.LL.cleaned) #Indicate for each occurrence point in which cell of the raster it falls
    global.occ.LL.cleaned<-global.occ.LL.cleaned[!is.na(global.occ.LL.cleaned$cell),]
    unique_occurrences <- !duplicated(global.occ.LL.cleaned$cell)# Identify unique occurrences
    global.occ.LL.cleaned <- global.occ.LL.cleaned[unique_occurrences, 1:2] # Subset the occurrence points to keep only one occurrence per raster cell 
    
    global.occ.LL.cleaned<- terra::extract(globalclimpreds_terra, global.occ.LL.cleaned, xy = T, ID=F)%>%
      dplyr::filter(rowSums(is.na(.[, 1:(ncol(.) - 2)])) == 0)%>% #Keep rows that do not have any NA values in column 1- 3rd last 
      dplyr::select(c(x,y))%>%
      dplyr::rename(decimalLongitude=x,
                    decimalLatitude=y) #Extract climatic values of occurrence points from each raster layer and remove occurrence points that fall in cells with NA values in at least one rasterlayer
    
    #Convert to sf dataframe
    global.occ.LL.cleaned$species<- rep(1,length(global.occ.LL.cleaned$decimalLongitude)) #adds columns indicating species presence (1) needed for modeling
    global.occ.sf<-st_as_sf(global.occ.LL.cleaned, coords=c("decimalLongitude", "decimalLatitude"), crs=4326, remove= FALSE)
    
    
    #--------------------------------------------
    #------ Define number of pseudoabsences -----
    #--------------------------------------------
    numb.global.pseudoabs <-length(global.occ.sf$decimalLongitude) #sets the number of pseudoabsences equal to number of unique presences
    rm(global.occ.LL.cleaned)
    
    
    #--------------------------------------------
    #------ Plot distribution of occurrences ----
    #--------------------------------------------
    ggplot()+ 
      geom_sf(data = world,  colour = "black", fill = NA)+
      geom_point(data=global.occ.sf, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
      labs(x="Longitude", y="Latitude")+
      theme_bw()
    
    
    #--------------------------------------------
    #- Select ecoregions containing occurrences -
    #--------------------------------------------
    wwf_eco<-sf::st_read(here("./data/external/GIS/official/wwf_terr_ecos.shp"))
    wwf_eco<-sf::st_transform(wwf_eco, 4326) %>%
      sf::st_make_valid()
    occ_ecoIntersect <- sf::st_intersects(wwf_eco,global.occ.sf) 
    wwf_ecoSub1<-wwf_eco[lengths(occ_ecoIntersect) > 0,1]
    
    
    #--------------------------------------------
    #------------- Plot ecoregions --------------
    #--------------------------------------------
    ggplot()+ 
      geom_sf(data = world,  colour = "black", fill = NA)+
      geom_sf(data=wwf_ecoSub1, fill="#f7786f")+
      labs(x="Longitude", y="Latitude")+
      theme_bw()
    
    
    #--------------------------------------------
    #------ Import right bias grid --------------
    #--------------------------------------------
    if (speciesgroup %in% names(bias_grid_paths)) {
      biasgrid <- terra::rast(bias_grid_paths[[speciesgroup]])
    } else {
      stop("No bias grid available for this species. Species has to be one of the following: Plants, Amphibians, Birds, Mammals, Molluscs, or Reptiles.")
    }
    
    
    #--------------------------------------------
    #Mask biasgrid by ecoregions with occurrences 
    #--------------------------------------------
    wwf_ecoSub1_ext<-terra::ext(wwf_ecoSub1) 
    wwf_ecoSub1_vector <- vect(wwf_ecoSub1) #Convert wwf_ecoSub1 to a SpatVector that can be used for masking
    biasgrid_crop <- terra::crop(biasgrid, wwf_ecoSub1_ext) #Crop biasgrid to extent wwf_ecoSub1
    biasgrid_sub <- terra::mask(biasgrid_crop, wwf_ecoSub1_vector)#Mask cropped biasgrid with SpatVector
    
    #Mask biasgrid with one of the climatic layers, to make sure it doesn't extend beyond them
    climategrid_rast<-terra::crop(globalclimpreds_terra[[1]], wwf_ecoSub1_ext)
    biasgrid_sub<-terra::mask(biasgrid_sub, climategrid_rast) 
    
    biasgrid_sub_raster <- raster(biasgrid_sub) #Convert SpatRaster back to normal raster object
    
    
    #--------------------------------------------
    #---------------Visualize biasgrid-----------
    #--------------------------------------------
    ggplot()+ 
      geom_sf(data = world,  colour = "black", fill = NA)+
      geom_spatraster(data=biasgrid_sub)+
      scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange")+
      labs(x="Longitude", y="Latitude")+
      theme_bw()
    
    
    #--------------------------------------------
    #---------- Generate pseudoabsences ---------
    #--------------------------------------------
    # generates pseudo absences equal to the number of presences, from cells that are not NA in biasgrid_sub and not from cells that have occurrence points as indicated in global.occ.sf. 
    #extf 1.1 increases the size of extent with 5% at each side of the extent (default value) but when ext is NULL it won't do anything.
    #excludep = TRUE indicates that presence points are excluded from the background, prob: if TRUE the values in mask are interpreted as probability weights (only works for rasters with a modest size that can be loaded into RAM)
    
    #Create alternative raster consisting of only ecoregions without biasgrid mask, used only when not enough pseudoabsence points can be generated using biasgrid_sub as mask layer
    ecoregions_crop<-terra::crop(globalclimpreds_terra[[1]], wwf_ecoSub1_ext) #Crop one of the climate rasters to extent ecoregions
    ecoregions_raster<-mask(ecoregions_crop,wwf_ecoSub1_vector) #Mask with ecoregions vector
    
    #Generate pseudoabsences
    set.seed(728)
    global_points <- generate_pseudoabs( mask = biasgrid_sub_raster, alternative_mask = raster(ecoregions_raster) , n = numb.global.pseudoabs, p =  st_drop_geometry(global.occ.sf))
    
    
    #--------------------------------------------
    #--- Create presence-pseudoabsence dataset---
    #--------------------------------------------
    global_pseudoAbs<-global_points %>%
      as.data.frame()%>%
      mutate(species = rep(0,length(x)))%>%
      st_as_sf(coords=c("x", "y"), crs=4326, remove=FALSE)%>%
      dplyr::rename(decimalLongitude=x,
                    decimalLatitude=y)
    
    global_presabs<- rbind(global.occ.sf,global_pseudoAbs)# join pseudoabsences with presences 
    rm(global_points)
    
    
    #--------------------------------------------
    #--Visualize presence-pseudoabsence dataset--
    #--------------------------------------------
    mapview(biasgrid_sub_raster, 
            col.regions = colorRampPalette(c("blue", "orange")),
            alpha=1, 
            na.color = "transparent", 
            layer.name = "Bias Grid") +
      mapview(global_presabs, zcol = "species", 
              col.regions = c("red", "yellow"),
              layer.name = "Species distribution")
    
    
    #--------------------------------------------
    #---- Extract climate data for modelling-----
    #--------------------------------------------
    global.data <- sdm::sdmData(species~.,train=vect(global_presabs),predictors=globalclimpreds_terra) 
    global.data.df<-as.data.frame(global.data)
    
    
    #--------------------------------------------
    #--- Remove highly correlated predictors---
    #--------------------------------------------
    # Identify highly correlated predictors
    correlationMatrix<-cor(global.data.df[,-c(1,2)]) #Calculate pearson correlation among environmental values
    highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7,exact=TRUE,names=TRUE)#Returns names of environmental variables to be removed because they are correlated more than 0.7 with other variables.  If two variables have a high correlation, the function removes the variable with the largest mean absolute correlation.
    preds<-as.data.frame(highlyCorrelated)
    kable(preds) %>%
      kable_styling(bootstrap_options = c("striped"))
    
    # Remove highly correlated predictors from dataframe 
    global.data.df.subset<- global.data.df %>%
      select (-all_of(highlyCorrelated), -rID) %>% 
      mutate(species = as.factor(species)) %>%
      mutate(species = recode_factor(species, 
                                     '0' = "absent",
                                     '1' = "present")) %>%  # Later steps require non numeric dependent variable
      mutate(species = relevel(species, ref = "present")) 
    
    
    #--------------------------------------------
    #--Correct climate data from integer format--
    #--------------------------------------------
    global.data.df.uncor<-cbind("species"=  global.data.df.subset$species,divide10(global.data.df.subset[,-c(1)]))
    
    
    #--------------------------------------------
    #--- Run multiple machine learning models ---
    #--------------------------------------------
    #preProc: preprocessing of predictors (environmental data). method = "center" subtracts the mean of the predictor's data (again from the data in x) from the predictor values while method = "scale" divides by the standard deviation.
    control <- caret::trainControl(method="cv",
                                   number=10,
                                   savePredictions="final", 
                                   preProc=c("center","scale"),
                                   classProbs=TRUE)
    classList1 <- c("glm","gbm","rf","earth")
    set.seed(457)
    global_train <- caretEnsemble::caretList(species~., 
                                             data= global.data.df.uncor,
                                             trControl=control,
                                             methodList=classList1,
                                             metric="Accuracy")
    
    
    #--------------------------------------------
    #--Return accurracy, kappa and correlation --
    #--------------------------------------------
    GlobalModelResults<-resamples(global_train) #Returns the results for each model 
    
    # Display accuracy of each model
    Global.Mod.Accuracy<-summary(GlobalModelResults)
    kable(Global.Mod.Accuracy$statistics$Accuracy,digits=2) %>%
      kable_styling(bootstrap_options = c("striped"))
    
    # Display kappa of each model
    kable(Global.Mod.Accuracy$statistics$Kappa,digits=2) %>%
      kable_styling(bootstrap_options = c("striped"))
    
    #Display correlation among models.
    #Weakly correlated algorithms are persuasive for stacking them in ensemble.
    Global.Mod.Cor<-modelCor(resamples(global_train))
    kable(Global.Mod.Cor,digits=2)%>%
      kable_styling(bootstrap_options = c("striped"))
    
    
    #--------------------------------------------
    #---------- Create ensemble model -----------
    #--------------------------------------------
    #combine individual models into one
    set.seed(478)
    global_stack <- caretEnsemble(
      global_train, 
      metric="Accuracy",
      trControl=trainControl(method="cv",
                             number=10,
                             savePredictions= "final",
                             classProbs=TRUE))
    print(global_stack)
    
    
    #--------------------------------------------
    #Identify best threshold and get accurracy
    #--------------------------------------------
    #Identify threshold that maximizes spec=sens
    global.ens.thresh<-findThresh(global_stack$ens_model$pred)
    #Return accurracy
    ensemble_accurracy<-accuracyStats(global_stack$ens_model$pred,global.ens.thresh$predicted)
    
    
    #--------------------------------------------
    #-- Get variable importance of global model--
    #--------------------------------------------
    variableImportance_global<-varImp(global_stack)
    kable(variableImportance_global,digits=2,caption="Variable Importance") %>%
      kable_styling(bootstrap_options = c("striped"))
    
    
    #--------------------------------------------
    #------------------ Clean up-----------------
    #--------------------------------------------
    rm(list = setdiff(ls(), c("eu_climpreds.10", "global_stack", "split_df", "taxonkey", "species", "ensemble_accurracy", "accuracyStats", "decimalplaces", "divide10", "findThresh", "predict_large_raster", "globalclimpreds_terra","bias_grid_paths", "i", "globalmodels","global.occ.sf", "biasgrid_sub", "world", "projectname", "first_two_words", "generate_pseudoabs", "variableImportance_global")))
    
    
    #--------------------------------------------
    #-------- Make predictions for Europe--------
    #--------------------------------------------
    system.time({
      global_model <- predict(eu_climpreds.10,global_stack,type="prob", na.rm = TRUE) #235.05
    })
    
    
    #--------------------------------------------
    #-------------- Plot predictions-------------
    #--------------------------------------------
    brks <- seq(0, 1, by=0.1) 
    nb <- length(brks)-1 
    # Generate Viridis palette
    viridis_palette <- viridis(nb)
    
    ggplot() + 
      #geom_sf(data = world,  colour = "grey", fill = NA)+
      geom_spatraster(data = global_model) +
      scale_fill_gradientn(colors = viridis_palette, breaks = brks, labels = brks, na.value = NA) +
      geom_sf(data = global.occ.sf, color = "black", fill = "red", size =1.5, shape = 21) +
      coord_sf(xlim = c(-10, 40), ylim = c(35, 72)) + 
      labs(fill = "Suitability")+
      theme_bw()
    
    
    #--------------------------------------------
    #-- Export results as .qs list
    #--------------------------------------------
    globalmodels <-list(species = species,
                        taxonkey = taxonkey,
                        global_ensemble_model = global_stack, 
                        model_accuracy = ensemble_accurracy,
                        variable_importance = variableImportance_global,
                        global_model_predictions = terra::wrap(global_model), #Needs to be wrapped to export it or will return a null pointer error
                        occurrences=global.occ.sf,
                        biasgrid=terra::wrap(biasgrid_sub)
    )
    
    qsave(globalmodels, paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey,"/Global_model_",first_two_words,"_",taxonkey,".qs"))
    
    print(paste("Global model has been created for", species))
    
  }
})


#--------------------------------------------
#---------- Clean R environment--------------
#--------------------------------------------
rm(list = ls())

