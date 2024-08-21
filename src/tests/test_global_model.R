# Load necessary packages
library(testthat)








#---------------------------------------------------------------------------------------------------------------------------------------------
#--------Function to run the workflow of Amy (https://github.com/amyjsdavis/wiSDM/blob/main/src/trias_sdm.Rmd) and return metrics-------------
#---------------------------------------------------------------------------------------------------------------------------------------------
run_old_workflow <- function(downloadkey) {
  
  #Retrieve downloaded records
  temp_dir<- tempdir()
  gbif_download_key<-downloadkey
  global <- rgbif::occ_download_get(downloadkey, path=temp_dir, overwrite = TRUE) %>%
    rgbif::occ_download_import()
  taxonkey<-unique(global$acceptedTaxonKey)
  
  ### 2. Create a global SDM 
  ##### 2. Specify paths for output (defaults to file structure in ReadMe)
  rasterOutput<-here("data/processed/geotiffs/")
  pdfOutput<-here("data/processed/pdf/")
  genOutput<-here("data/processed/general//")

  ####3. Filter global occurrence data 

  decimalplaces <- function(x) {
    if (abs(x - round(x)) > .Machine$double.eps^0.5) {
      nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
    } else {
      return(0)
    }
  }

  

  
  #remove unverified records
  identificationVerificationStatus_to_discard <- c("unverified", "unvalidated","not able to validate","control could not be conclusive due to insufficient knowledge")
  
  #enter value for max coordinate uncertainty in meters.
  
  global.occ<-global %>%
    filter(speciesKey==taxonkey) %>%   #using taxonKey filters out accepted synonyms
    filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters<= 1000) %>%
    filter(!str_to_lower(identificationVerificationStatus) %in% identificationVerificationStatus_to_discard)
  
  global.occ$lon_dplaces<-sapply(global.occ$decimalLongitude, function(x) decimalplaces(x))
  global.occ$lat_dplaces<-sapply(global.occ$decimalLatitude, function(x) decimalplaces(x))
  global.occ[global.occ$lon_dplaces < 4& global.occ$lat_dplaces < 4 , ]<-NA
  global.occ<-global.occ[ which(!is.na(global.occ$lon_dplaces)),]
  global.occ<-within(global.occ,rm("lon_dplaces","lat_dplaces"))
  global.occ<-global.occ[which( global.occ$year > 1970 & global.occ$year < 2011),]

  
  #### Convert global occurrences to spatial points needed for modelling
  
  global.occ<-global.occ[c("decimalLongitude", "decimalLatitude")]
  coordinates(global.occ)<- c("decimalLongitude", "decimalLatitude")
  global.occ.LL<-data.frame(global.occ)[c(1:2)] #extract long and lat 
  
  #### Flag and remove centroids and invalid georeferenced points
  global.occ.LL$species<-rep("Vaccinium corymbosum",nrow(global.occ.LL))  
  
  flags_report<-clean_coordinates(x = global.occ.LL, lon= "decimalLongitude", lat= "decimalLatitude",
                                  tests = c("capitals", 
                                            "centroids","gbif", "institutions", 
                                            "seas", "zeros"))
  
  cleaned<-clean_coordinates(x = global.occ.LL, lon= "decimalLongitude", lat= "decimalLatitude",
                             tests = c("capitals", 
                                       "centroids","gbif", "institutions", 
                                       "seas", "zeros"),value="clean")
  global.occ.LL.cleaned<-subset(cleaned,select= -c(species))
  
  #### Create global rasterstack using CHELSA data for model building
  globalclimrasters <- list.files((here("./data/external/climate/trias_CHELSA")),pattern='tif',full.names = T) #import CHELSA data
  globalclimpreds <- stack(globalclimrasters)

  
  #### Use SDMtab command from the SDMPlay package to remove duplicates per grid cell
  
  global.SDMtable<- SDMPlay:::SDMtab(global.occ.LL.cleaned, globalclimpreds, unique.data = TRUE,background.nb= 0) #
  numb.global.pseudoabs <-length(global.SDMtable$id) #sets the number of pseudoabsences equal to number of unique presences
  
  
  global.occ.sp<-global.SDMtable[c("longitude", "latitude")]
  coordinates(global.occ.sp)<- c("longitude", "latitude")
  global.occ.sp$species<- rep(1,length(global.occ.sp$latitude)) #adds columns indicating species presence needed for modeling
  
  ### Select wwf ecoregions that contain global occurrence points
  wwf_eco<-sf::st_read(here("./data/external/GIS/official/wwf_terr_ecos.shp"))
  wwf_eco<-as_Spatial(wwf_eco)
  crs(global.occ.sp)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  crs(wwf_eco)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  occ_ecoIntersect <- over(wwf_eco,global.occ.sp) 
  wwf_ecoSub1 <- wwf_eco[!is.na(occ_ecoIntersect$species),]

  ### Specify and import bias grids for relevant taxonomic group (e.g vascular plants) 
  biasgrid<-raster(here("./data/external/bias_grids/final/trias/plants_1deg_min5.tif"))### specify appropriate bias grid here

  ### Subset bias grid by ecoregions containing occurrence points
  ext_wwf_ecoSub<-extent(wwf_ecoSub1)
  biasgrid_crop<-crop(biasgrid,ext_wwf_ecoSub)
  biasgrid_sub<-mask(biasgrid_crop,wwf_ecoSub1)

  
  ###  Use randomPoints function from dismo package to locate pseduobasences within the bias grid subset 
  # generates pseudo absences equal to (or close to) the number of presences.
  set.seed(728)
  global_points<-randomPoints(biasgrid_sub,numb.global.pseudoabs, global.occ.sp, ext=NULL, extf=1.1, excludep=TRUE, prob=FALSE, cellnumbers=FALSE, tryf=70, warn=2, lonlatCorrection=TRUE) 
  # will throw a warning if randomPoints generated is less than numb.pseudoabs. If this happens, increase the number of tryf or ignore bias grid and sample from ecoregion only.

  ### Extract generated pseudo absences and create presence-pseudobasence dataset 
  global_pseudoAbs<-as.data.frame(global_points)
  coordinates(global_pseudoAbs)<-c("x","y")
  crs(global_pseudoAbs)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  global_pseudoAbs$species<-rep(0,length(global_pseudoAbs$x))
  global_presabs<- spRbind(global.occ.sp,global_pseudoAbs) # join pseudoabsences with presences (occurrences)

  
  ### Extract climate data for global scale modelling
  global.data <- sdmData(species~.,train=global_presabs, predictors=globalclimpreds)
  global.data.df<-as.data.frame(global.data)

  ### Identify highly correlated predictors
  correlationMatrix<-cor(global.data.df[,-c(1)])
  highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7,exact=TRUE,names=TRUE)
  preds<-as.data.frame(highlyCorrelated)
  kable(preds) %>%
    kable_styling(bootstrap_options = c("striped"))

  
  ### Remove highly correlated predictors from dataframe 
  global.data.df.subset<-select (global.data.df,-c(highlyCorrelated))
  global.data.df.subset<-within(global.data.df.subset,rm("rID"))
  global.data.df.subset$species<-as.factor(global.data.df.subset$species) #later steps require non numeric dependent variable
  levels(global.data.df.subset$species)<-c("absent","present")
  global.data.df.subset$species <- relevel(global.data.df.subset$species, ref = "present") 

  ### Correct global clim preds values from integer format
  divide10<-function(x){
    value<-x/10
    return(value)
  }
  
  global.data.df.uncor<-cbind("species"=  global.data.df.subset$species,divide10(global.data.df.subset[,-c(1)]))

  
  
  
  ### Use caretList from Caret package to run multiple machine learning models
  control <- trainControl(method="cv",number=10,savePredictions="final", preProc=c("center","scale"),classProbs=TRUE)
  classList1 <- c("glm","gbm","rf","earth")
  set.seed(457)
  global_train <- caretList(
    species~., data= global.data.df.uncor,
    trControl=control,
    methodList=classList1)

  GlobalModelResults<-resamples(global_train)
  Global.Mod.Accuracy<-summary(GlobalModelResults)# displays accuracy of each model
  kable(Global.Mod.Accuracy$statistics$Accuracy,digits=2) %>%
    kable_styling(bootstrap_options = c("striped"))

  GlobalModelResults<-resamples(global_train)
  kable(Global.Mod.Accuracy$statistics$Kappa,digits=2) %>%
    kable_styling(bootstrap_options = c("striped"))

  Global.Mod.Cor<-modelCor(resamples(global_train))# shows correlation among models.Weakly correlated algorithms are persuasive for stacking them in ensemble.
  kable(Global.Mod.Cor,digits=2)%>%
    kable_styling(bootstrap_options = c("striped"))

  ### Create ensemble model (combine individual models into one) 
  set.seed(478)
  global_stack <- caretEnsemble(
    global_train, 
    trControl=trainControl(method="cv",
                           number=10,
                           savePredictions= "final",classProbs=TRUE ))
  print(global_stack)

  
  ### Function to return threshold where sens=spec from caret results 
  findThresh<-function(df){
    df[c("rowIndex","obs","present")]
    df<-df %>%
      mutate(observed= ifelse(obs == "present",1,0)) %>%
      select(rowIndex,observed,predicted=present)
    result<-PresenceAbsence::optimal.thresholds(df,opt.methods = 2)
    return(result)
  }
  
  #accuracy measures
  accuracyStats<-function(df,y){
    df[c("rowIndex","obs","present")]
    df<-df %>%
      mutate(observed= ifelse(obs == "present",1,0)) %>%
      select(rowIndex,observed,predicted=present)
    result<-PresenceAbsence::presence.absence.accuracy(df,threshold = y,st.dev=FALSE)
    return(result)
  }

  
  ### Identify threshold and performance of global ensemble model

  global.ens.thresh<-findThresh(global_stack$ens_model$pred)
  results<-accuracyStats(global_stack$ens_model$pred,global.ens.thresh$predicted)

  return(results)
  }





#------------------------------------------------------------------------
#--------Function to run the new workflow and return metrics-------------
#------------------------------------------------------------------------
run_new_workflow <- function() {
  
  #Retrieve downloaded records
  temp_dir<- tempdir()
  gbif_download_key<-downloadkey
  global <- rgbif::occ_download_get(downloadkey, path=temp_dir, overwrite = TRUE) %>%
    rgbif::occ_download_import()
  taxonkey<-unique(global$acceptedTaxonKey)
  
  
  
  
  
  return(results)
}




#------------------------------------------------------------------------
#--------Function to run the new workflow and return metrics-------------
#------------------------------------------------------------------------

# Run the workflows
old_results <- run_old_workflow(downloadkey="0076914-240626123714530")
new_results <- run_new_workflow(downloadkey="0076914-240626123714530")

# Run test
test_that("New workflow produces the same metrics as the old workflow", {
  expect_identical(old_results$AUC, new_results$AUC)
  expect_identical(old_results$specificity, new_results$specificity)
  expect_identical(old_results$sensitivity, new_results$sensitivity)
  expect_identical(old_results$Kappa, new_results$Kappa)
  expect_identical(old_results$PCC, new_results$PCC)
  expect_identical(old_results$threshold, new_results$threshold)
})

# If there are differences, print them 
if (!isTRUE(all.equal(old_results, new_results))) {
  cat("Differences detected in the workflow outputs:\n")
  print(all.equal(old_results, new_results))
} else {
  cat("No differences detected. The workflows produce the same metrics 🎉.\n")
}
