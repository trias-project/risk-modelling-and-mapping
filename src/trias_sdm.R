#R script to use for generating maps from GBIF global download of occurrence data

library(sdm)
library(rgdal)
library(raster)
library(spocc)
library(sp)
library(rgeos)
library(SDMPlay)
library(tidyverse)
library(caret)
library(caretEnsemble)
library(trias)
library(rgbif)
library(RColorBrewer)
library(maptools)
library(dismo)
library(sf)
library(geoR)
######GLOBAL SDM FOR WEIGHTING PSEUDOABSENCES IN EUROPEAN-LEVEL SDMs
#all paths are relative to the risk_modelling folder. Must have a fresh session for this to work.
setwd("C:../risk_modelling")


###read in global download for species from GBIF  
gbif_filename<- ".csv" #add name of species being modeled
global<-read.csv(file=paste("./data/external/PRA_mammals/gbif_speciesFiles/",gbif_filename, sep=""))


 taxonkey<-"" #add GBIF taxonKey for the species being modeled.
 taxonName<-" " #add name of species being modeled
 

##filter data to keep only those points with acceptable spatial accuracy (for us this 4 decimal places for either lat or lon) 
#script to count number of decimal places
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}

##this option is for retrieving synomyms when working with subspecies                                                                                                    
#taxonlist<- c("8421432", "3663378","3663277")

#remove unverified records
identificationVerificationStatus_to_discard <- c( "unverified", "unvalidated","not able to validate",
                                                  "control could not be conclusive due to insufficient knowledge")
global.occ<-global %>%
  #filter(taxonKey==taxonkey) %>%   #using taxonKey filters out accepted synonyms
  #filter(taxonKey%in% taxonlist) %>%  this option is for retrieving synomyms when working with subspecies  
  filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters< 1000) %>%
  filter(!str_to_lower(identificationVerificationStatus) %in% identificationVerificationStatus_to_discard)

global.occ$lon_dplaces<-sapply(global.occ$decimalLongitude, function(x) decimalplaces(x))
global.occ$lat_dplaces<-sapply(global.occ$decimalLatitude, function(x) decimalplaces(x))
global.occ[global.occ$lon_dplaces < 4& global.occ$lat_dplaces < 4 , ]<-NA
global.occ<-global.occ[ which(!is.na(global.occ$lon_dplaces)),]
global.occ<-within(global.occ,rm("lon_dplaces","lat_dplaces"))
global.occ<-global.occ[which( global.occ$year > 1975 & global.occ$year < 2006),]



###Create SpatialPoints dataframe
global.occ<-global.occ[c("decimalLongitude", "decimalLatitude")]
coordinates(global.occ)<- c("decimalLongitude", "decimalLatitude")
plot(global.occ)


###PREDICTOR VARIABLES USING GLOBAL CLIMATE DATA (CHELSA)####
#create global rasterstack using CHELSA data for model building
globalclimrasters <- list.files("C:./data/external/climate/trias_CHELSA",pattern='tif',full.names = T) #insert path where all climate rasters are located
globalclimpreds <- stack(globalclimrasters)


##use SDMtab command from the SDMPlay package to remove duplicates per grid cell  

global.occ.LL<-data.frame(global.occ)[c(1:2)] #extract long and lat neededfor SDMtable
global.SDMtable<- SDMPlay:::SDMtab(global.occ.LL, globalclimpreds, unique.data = TRUE,background.nb= 0) #
numb.pseudoabs <- length(global.SDMtable$id) #sets the number of pseudoabsences equal to number of unique presences


##transform filtered occurrence dataset with unique presences back to a SpatialPoints dataframe. 
global.occ.sp<-global.SDMtable[c("longitude", "latitude")]
coordinates(global.occ.sp)<- c("longitude", "latitude")
global.occ.sp$species<- rep(1,length(global.occ.sp$latitude)) #adds columns indicating species presence needed for modeling


#select wwf ecoregions that contain occurrence points
wwf_eco<-shapefile("C:./data/external/GIS/official/wwf_terr_ecos.shp")
crs(global.occ)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
occ_ecoIntersect <- over(wwf_eco,global.occ) 
wwf_ecoSub1 <- wwf_eco[!is.na(occ_ecoIntersect),]


## import bias grids for relevant taxonomic group (e.g vascular plants) counts of less than 5 removed
biasgrid<-raster("C:./data/external/bias_grids/final/trias/mammals_1deg_min5.tif")

# subset bias grid by ecoregions containing occurrence points
ext_wwf_ecoSub<-extent(wwf_ecoSub1)
biasgrid_crop<-crop(biasgrid,ext_wwf_ecoSub)
biasgrid_sub<-mask(biasgrid_crop,wwf_ecoSub1)
plot(biasgrid_sub) #not run


# use randomPoints function from dismo package to locate pseduobasences within the bias grid subset
set.seed(728)
global_points<-randomPoints(biasgrid_sub, numb.pseudoabs, global.occ.sp, ext=NULL, extf=1.1, excludep=TRUE, prob=FALSE,
                                cellnumbers=FALSE, tryf=50, warn=2, lonlatCorrection=TRUE) #will throw a warning if randomPoints generated is less than numb.pseudoabs. If this happens, increase the number of tryf.

# extract generated pseudo absences and create presence-pseudobasence dataset 
global_pseudoAbs<-as.data.frame(global_points)
coordinates(global_pseudoAbs)<-c("x","y")
global_pseudoAbs$species<-rep(0,length(global_pseudoAbs$x))
global_presabs<- spRbind(global.occ.sp,global_pseudoAbs) # join pseudoabsences with presences (occurrences)
writeOGR(global_presabs, dsn="C:./data/results/sdm_occ_pts/global", layer=paste(taxonName,"_globalSDM"), driver = "ESRI Shapefile", overwrite_layer=TRUE)


# extract data for modelling
#use  sdmData to extract climate data for each occurrence
global.data <- sdmData(species~.,train=global_presabs, predictors=globalclimpreds)
global.data.df<-as.data.frame(global.data)

# identify highly correlated predictors
correlationMatrix<-cor(global.data.df[,-1])
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7,exact=TRUE,names=TRUE)
print(highlyCorrelated)

# remove highly correlated predictors from dataframe 
global.data.df.subset<-select (global.data.df,-c(highlyCorrelated))
global.data.df.subset<-within(global.data.df.subset,rm("rID"))
global.data.df.subset$species<-as.factor(global.data.df.subset$species)#later steps require non numeric dependent variable
levels(global.data.df.subset$species)<-c("absent","present")
global.data.df.subset$species <- relevel(global.data.df.subset$species, ref = "present")

# correct global clim preds values from integer format
divide10<-function(x){
  value<-x/10
  return(value)
}

global.data.df.uncor<-cbind("species"= global.data.df$species,divide10(global.data.df[,-c(1)]))

#####################  
###START MODELLING###
#####################
# use Caret package to run multiple machine learning models
control <- trainControl(method="repeatedcv",number=5, repeats=5, savePredictions="final", preProc=c("center","scale"),classProbs=TRUE)
classList1 <- c("glm","gbm","rf","adaboost","knn")
set.seed(457)
global_train <- caretList(
  species~., data= global.data.df.uncor,
  trControl=control,
  methodList=classList1
)

modelResults<-resamples(global_train)
summary(modelResults)
modelCor(resamples(global_train))

# create ensemble model  
set.seed(478)
global_stack <- caretStack(
  global_train, 
  method="glm",
  trControl=trainControl(method="cv",
                         number=10,
                         savePredictions= "final"  ))
print(global_stack)


# create rasterstack of CHELSA data clipped to European modeling extent for prediction
euclimrasters <- list.files("C:./data/external/climate/chelsa_eu_clips",pattern='tif',full.names = T)
eu_climpreds<-stack(euclimrasters)
eu_climpreds.10<-divide10(eu_climpreds) # correct for integer format of Chelsa preds 

# prediction to the extent of Europe
global_model<-raster::predict(eu_climpreds.10,global_stack,type="prob")
writeRaster(global_model, filename=paste("GlobalEnsEU_",taxonkey, ".tif",sep=""), format="GTiff",overwrite=TRUE) 

##########################END GLOBAL MODEL####################################################################################################################

############################################
#########EUROPEAN LEVEL MODEL##############
###########################################

#obtain european presences from the occurrence cube 
 cube_europe <- read.csv("C:./data/external/data_cube/eu_modellingtaxa_cube.csv")
 
#extract occurrence data for the target species meeting our criteria
 
 occ.eu<-cube_europe %>% 
   filter(taxonKey==taxonkey) %>% 
   filter(year >"1975") %>% #select presences 1976 or later to be consistent with climate data
  # filter(year <"2006") %>%
   filter(min_coord_uncertainty <= 1000) 
 
 occ.eu$eea_cell_code<-str_remove(occ.eu$eea_cell_code,"1km")#so that cell code matches with chelsa data
 
# join lat lons from EEA 1km grid centroids to occurrence cube
centroids<-st_read("./data/external/GIS/EEA_fullgrid_1km_final_centroids.shp")
 merged<-merge(occ.eu,centroids,by.x="eea_cell_code",by.y="EEA_fullgr") 
occ.eu1<-cbind(merged[c("eea_cell_code","latitude","longitude")])  

###Create SpatialPoints dataframe needed for SDMtab 
euocc<-occ.eu1[c("longitude", "latitude")]#extract long and lat for SDMtable
coordinates(euocc)<- c("longitude", "latitude")
euocc1<-data.frame(euocc)[c(1:2)] 



###PREDICTOR VARIABLES USING EUROPEAN LEVEL CLIMATE DATA FROM RMI ####
# create RasterStack of european climate variables
rmiclimrasters <- list.files("C:./data/external/climate/rmi_corrected",pattern='tif',full.names = T) #insert path where all climate rasters are located
rmiclimpreds <- stack(rmiclimrasters)


# use SDMtab command from the SDMPlay package to remove duplicates per grid cell  
euocc.SDMtable<- SDMPlay:::SDMtab(euocc1, rmiclimpreds, unique.data = TRUE,background.nb= 0)
numb.pseudoabs <- length(euocc.SDMtable$id) 



#transform eu occurrence dataset with unique presences back to a SpatialPoints dataframe. 
euocc<-euocc.SDMtable[c("longitude", "latitude")]
coordinates(euocc)<- c("longitude", "latitude")
euocc$occ<- rep(1,length(euocc$latitude))#adds columns indicating species presence needed for modeling
crs(euocc)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs ")



#clip bias grid to European extent
euboundary<-shapefile("C:/Users/amyjs/Documents/projects/Trias/modeling/GIS/EUROPE.shp") 
studyextent<-euboundary
ecoregions_eu<-crop(biasgrid_sub,studyextent)
biasgrid_eu<-projectRaster(ecoregions_eu,rmiclimpreds)

###############################################################################
# Mask areas of high habitat suitability from global climate model
#Read in global raster from earlier step if needed
#globalfilename<-paste("C:./GlobalEnsEU_", taxonkey, ".tif",sep="")
#global_model<-raster(globalfilename)

wgs84_gcs<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
crs(global_model)<-wgs84_gcs
m<-global_model >.5
global_mask<-raster:::mask(global_model,m,maskvalue=TRUE)
global_masked_proj<-projectRaster(global_mask,biasgrid_eu)

# overlay low predicted habitat suitability on bias grid to exclude low sampled areas 
pseudoSamplingArea<-mask(biasgrid_eu,global_masked_proj)
plot(pseudoSamplingArea)

# ranndomly locate pseudo absences within "pseudoSamplingArea" 
euocc_points<-randomPoints(pseudoSamplingArea, numb.pseudoabs, euocc, ext=NULL, extf=1.1, excludep=TRUE, prob=FALSE,
                              cellnumbers=FALSE, tryf=50, warn=2, lonlatCorrection=TRUE)


euocc_pseudoAbs<-as.data.frame(euocc_points)
coordinates(euocc_pseudoAbs)<-c("x","y")
euocc_pseudoAbs$occ<-rep(0,length(euocc_pseudoAbs$x))

crs(euocc_pseudoAbs)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs ")

#join filtered cube occurrences with pseudo absences 

eu_presabs<- spRbind(euocc,euocc_pseudoAbs)
crs(eu_presabs)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs ")
#confirm equal number of presences and absences
table(eu_presabs@data$occ)

#export occ points as shapefile 
writeOGR(obj=eu_presabs, dsn="C:./data/results/sdm_occ_pts/europe", layer=paste(taxonName,"_EUpresabs",sep=""), driver="ESRI Shapefile", overwrite_layer = TRUE) 


#create SDM data object for European data 
occeu.sdmdata <- sdmData(occ~.,train=eu_presabs, predictors=rmiclimpreds) 
occeu.sdmdata

##########################################################################################

#identify and filter out highly correlated predictors########

#convert eu data to dataframe
occeu.sdmdata.df<-as.data.frame(occeu.sdmdata)

#identify highly correlated predictors drop first two columns which always be rID and occ
correlationMatrix <- cor(occeu.sdmdata.df[,-c(1:2)])

# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected 
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.70,exact=TRUE,names=TRUE)

# remove highly correlated predictors from RMI rasterstack

rmiclimpreds_uncor<-dropLayer(rmiclimpreds,highlyCorrelated)

##########################################################################################################################

# add habitat and anthropogenic predictors
habitat<-list.files("C:./data/external/habitat/landcover",pattern='tif',full.names = T)
habitat_stack<-stack(habitat)
fullstack<-stack(rmiclimpreds_uncor,habitat_stack) #combine uncorrelated climate variable selected earlier with habitat

# uncomment below for vertebrates and invertebrates to include distance to water
dist2water<-raster("C:./data/external/habitat/distance2water_EEA_1km.tif")
dist2water<-extend(dist2water,habitat_stack)
fullstack<-stack(rmiclimpreds_uncor,habitat_stack,dist2water)

# clip fullstack to belgium extent
belgie<-shapefile("C:./data/external/GIS/belgium_boundary.shp")
fullstack_crop<-crop(fullstack,belgie)
fullstack_be<-mask(fullstack_crop,belgie)

occ.full.data <- sdmData(occ~.,train=eu_presabs, predictors=fullstack) 
occ.full.data

# convert eu data to dataframe
occeu.full.data.df<-as.data.frame(occ.full.data)
occeu.full.data.df<-occeu.full.data.df[,-1]


#identify and remove highly correlated predictors from the full stack (combined stack of habitat/anthropogenic and climate)
correlationMatrix <- cor(occeu.full.data.df[,-1])

# find attributes that are highly corrected 
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7,exact=TRUE,names=TRUE)

# print names of highly correlated attributes
print(highlyCorrelated)

#remve highly correlated preds
occeu.full.data.df1<-select (occeu.full.data.df,-c(highlyCorrelated))
nzv_preds<-nearZeroVar(occeu.full.data.df1,names=TRUE)
occeu.full.data.df1<-select (occeu.full.data.df1,-c(nzv_preds)) #

# build models with climate and habitat data

occeu.full.data.df1$occ<-as.factor(occeu.full.data.df1$occ)
levels(occeu.full.data.df1$occ)<-c("absent","present")

control <- trainControl(method="cv",number=10,savePredictions="final", preProc=c("center","scale"),classProbs=TRUE)
mylist<-list(
  glm =caretModelSpec(method = "glm",maxit=50),
  gbm= caretModelSpec(method = "gbm"),
  rf = caretModelSpec(method = "rf", importance = TRUE),
  adaboost =caretModelSpec(method = "adaboost"),
knn = caretModelSpec(method = "knn"))
set.seed(457)
model_train_habitat <- caretList(
  occ~., data=occeu.full.data.df1,
  trControl=control,
  tuneList=mylist
)

## display model evaluation statistics
modelResults1<-resamples(model_train_habitat)
summary(modelResults1)
modelCor(resamples(model_train_habitat))

#create ensemble
set.seed(478)
lm_ens_hab<-caretEnsemble(model_train_habitat, trControl=trainControl(method="cv",
                                                                  number=10,
                                                                  savePredictions= "final",classProbs = TRUE))

lm_ens_hab
variableImportance<-varImp(lm_ens_hab)
write.csv(variableImportance,file = paste(taxonkey,"_varImp.csv"))

#################################################################
#function to export pdfs 
exportPDF<-function(rst,taxonkey,taxonName,nameextension,is.diff="FALSE"){
  filename=paste("be_",taxonkey, "_",nameextension,sep="")
  pdf(file=filename,width=10,height=8,paper="a4r")
  par(bty="n")#to turn off box aroudn plot
  ifelse(is.diff=="TRUE", brks<-seq(-1, 1, by=0.25), brks <- seq(0, 1, by=0.1)) 
  nb <- length(brks)-1 
  pal <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
  cols<-pal(nb)
  maintitle<-paste(taxonName,filename, sep= " ")
  plot(rst, breaks=brks, col=cols,main=maintitle, lab.breaks=brks,axes=FALSE)
  dev.off() 
}  



#################################
##  use eu level model to predict for europe (uncomment)
# ens_pred_hab_eu<-raster::predict(fullstack,lm_ens_hab,type="prob")
# ens_pred_hab_eu1<-1-ens_pred_hab_eu
# writeRaster(ens_pred_hab_eu1,filename=paste("eu_",taxonkey, "_hist.tif",sep="") , format="GTiff",overwrite=TRUE) 

#######################
#predict for Belgium only using eu-level ensemble model to forecast risk under historical climate conditions
#etrs89<-CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs" )
laea_grs80<-crs(belgie)
ens_pred_hab<-raster::predict(fullstack_be,lm_ens_hab,type="prob")
ens_pred_hab1<-1-ens_pred_hab
crs(ens_pred_hab1)<-laea_grs80
writeRaster(ens_pred_hab1, filename=paste("be_",taxonkey, "_hist.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab1,taxonkey,taxonName=taxonName,"hist.pdf")

# use eu-level ensemble model to forecast risk under RCP climate change scenarios
## clip habitat stack to belgium while removing historical climate layers to combine with rcp climate scenarios
habitat_stack<-stack(habitat,dist2water)
habitat_only_stack<-crop(habitat_stack,belgie)
habitat_only_stack_be<-crop(habitat_only_stack,belgie)

## create individual RCP climate data stacks for Belgium
be26 <- list.files("C:./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26",pattern='tif',full.names = T)
belgium_stack26 <- stack(be26)

be45 <- list.files("C:./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45",pattern='tif',full.names = T)
belgium_stack45 <- stack(be45)

be85 <- list.files("C:./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85",pattern='tif',full.names = T)
belgium_stack85 <- stack(be85)

## combine habitat stacks with climate stacks for each RCP scenario
fullstack26<-stack(be26,habitat_only_stack_be)
fullstack45<-stack(be45,habitat_only_stack_be)
fullstack85<-stack(be85,habitat_only_stack_be)

## export RCP risk maps for each RCP scenario
ens_pred_hab26<-raster::predict(fullstack26,lm_ens_hab,type="prob")
ens_pred_hab26_1<-1-ens_pred_hab26
crs(ens_pred_hab26_1)<-laea_grs80
writeRaster(ens_pred_hab26_1, filename=paste("be_",taxonkey, "_rcp26.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab26_1,taxonkey,taxonName=taxonName,"rcp26.pdf")

ens_pred_hab45<-raster::predict(fullstack45,lm_ens_hab,type="prob")
ens_pred_hab45_1<-1-ens_pred_hab45
crs(ens_pred_hab45_1)<-laea_grs80
writeRaster(ens_pred_hab45_1, filename=paste("be_",taxonkey, "_rcp45.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab45_1,taxonkey,taxonName=taxonName,"rcp45.pdf")

ens_pred_hab85<-raster::predict(fullstack85,lm_ens_hab,type="prob")
ens_pred_hab85_1<-1-ens_pred_hab85
crs(ens_pred_hab85_1)<-laea_grs80
writeRaster(ens_pred_hab85_1, filename=paste("be_",taxonkey, "_rcp85.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab85_1,taxonkey,taxonName=taxonName,"rcp85.pdf")

# create and export "difference maps": the difference between predicted risk by each RCP scenario and historical climate

hist26_diff_hab<-overlay(ens_pred_hab26_1, ens_pred_hab1, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist26_diff_hab,filename=paste("be_",taxonkey, "_rcp26_diff.tif",sep="") , format="GTiff",overwrite=TRUE) 
exportPDF(hist26_diff_hab,taxonkey,taxonName=taxonName,"rcp26_diff.pdf","TRUE")

hist45_diff_hab<-overlay(ens_pred_hab45_1, ens_pred_hab1, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist45_diff_hab,filename=paste("be_",taxonkey, "_rcp45_diff.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(hist45_diff_hab,taxonkey,taxonName=taxonName,"rcp45_diff.pdf","TRUE")

hist85_diff_hab<-overlay(ens_pred_hab85_1, ens_pred_hab1, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist85_diff_hab, filename=paste("be_",taxonkey, "_rcp_85_diff.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(hist85_diff_hab,taxonkey,taxonName=taxonName,"rcp85_diff.pdf","TRUE")

# output data used for SDM
df4export<-cbind(occeu.full.data.df1,coordinates(occ.full.data))
write.csv(df4export,file=paste(taxonkey,"_sdmdata.csv",sep=""))


# Check spatial autocorrelation of residuals to assess whether occurrence data should be thinned

# derive residuals from unthinned model
predEns1<-lm_ens_hab$ens_model$pred
obs.numeric<-ifelse(predEns1$obs == "absent",0,1)

# # derive sensitivity and specificity
# table(predEns1$pred,predEns1$obs)
# sensitivity(predEns1$pred,predEns1$obs)
# specificity(predEns1$pred,predEns1$obs)


#######################################################################################################
## standardize the residuals: generalized residuals/sq root of the variance according to Davis et al 2015
stdres<-function(obs.numeric, yhat){
  num<-obs.numeric-yhat
  denom<-sqrt(yhat*(1-yhat))
  return(num/denom)
}
hab.res<-stdres(obs.numeric,predEns1$present)


res.best.coords<-cbind(coordinates(occ.full.data),hab.res)
res.best.geo<-as.geodata(res.best.coords,coords.col=1:2,data.col = 3)
summary(res.best.geo) #note distance is in meters

# optional: run variogram
#res.best.vario<-variog(res.best.geo,coords=res.best.geo$coords, data=res.best.geo$data,max.dist=50000,option = "bin")
#plot(res.best.vario)


# check MoransI
##if Moran's I is very low, or not significant, skip thinning.
library(ape)
res.best.df<-as.data.frame(res.best.coords)
occ.dists <- as.matrix(dist(cbind(res.best.df[1], res.best.df[2])))
occ.dists.inv <- 1/occ.dists
diag(occ.dists.inv) <- 0
Moran.I(res.best.df$hab.res,occ.dists.inv,scaled=TRUE,alternative="greater")


##########conformal prediction######################### 
#####for non-thinned  historical climate data with habitat variables
######################################################################
#Run the code in the section to generate maps of confidence for belgian risk maps

GetLength<-function(x,y){
  length(x[which(x >= y)])
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

#extract values from raster stack for use in conformal prediction algorithm in a way that they can be linked back to the belgium raster (ie include spatial info)

extractVals<-function(predras){
  vals <- values(predras)
  coord <- xyFromCell(predras,1:ncell(predras))
  raster_fitted <- cbind(coord,vals)
  raster_fitted.df<-as.data.frame(raster_fitted)
  raster_fitted.df1<-na.omit(raster_fitted.df)
  raster_fitted.df1$absence<-raster_fitted.df1$vals
  raster_fitted.df1$presence<- (1-raster_fitted.df1$absence)
  return(raster_fitted.df1)
}


####combine into single function
# Arguments
# output from call to caretEnsemble (e.g. lm_ens_hab) 
# y predicted raster output from applying caret ensemble model (e.g ens_pred_hab) to rasterstack


classConformalPrediction<-function(x,y){
ens_results<- get("x")
ens_calib<-ens_results$ens_model$pred
calibPresence<-subset(ens_calib$present,ens_calib$obs=='present',)
calibAbsence<-subset(ens_calib$absent,ens_calib$obs=='absent',)
predicted.values<-extractVals(y)

testPresence<-predicted.values$presence
testAbsence<-predicted.values$absence

#derive p.Values for class A
smallrA<-lapply(testPresence,function(x) GetLength(calibPresence,x))
smallrA_1<- unlist (smallrA)
nCalibSet<-length(calibPresence)+1
pvalA<-smallrA_1/nCalibSet

# derive p.Values for Class B
smallrB<-lapply(testAbsence,function(x) GetLength(calibAbsence,x))
smallrB_1<- unlist (smallrB)
nCalibSetB<-length(calibAbsence)+1
pvalB<-smallrB_1/nCalibSetB

pvalsdf<-as.data.frame(cbind(pvalA,pvalB,0.20))
raster_cp_20<-mapply(CPconf,pvalsdf$pvalA,pvalsdf$pvalB,pvalsdf[3])
table(raster_cp_20)

pvalsdf$conf<-get.confidence(pvalsdf$pvalA,pvalsdf$pvalB)
pvalsdf_1<-cbind(pvalsdf,predicted.values)
}

pvalsdf_hist<-classConformalPrediction(lm_ens_hab,ens_pred_hab)

# export csv 
write.csv(pvalsdf_hist,file=paste("confidence_",taxonkey, "_hist.csv",sep=""))


##############################
###Create and export confidence maps
confidenceMaps<-function(x,taxonkey,taxonName,maptype){
pvals_dataframe<-get("x")
data.xyz <- pvals_dataframe[c("x","y","conf")]
rst <- rasterFromXYZ(data.xyz)
crs(rst)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs") 
plot(rst)
writeRaster(rst, filename=paste("be_",taxonkey, "_",maptype,".tif",sep=""), format="GTiff",overwrite=TRUE)
exportPDF(rst,taxonkey,taxonName=taxonName,nameextension= paste(maptype,".pdf",sep=""))
}



hist.conf.map<-confidenceMaps(pvalsdf_hist,taxonkey,taxonName,maptype="hist_conf")
rcp26.conf.map<-confidenceMaps(pvalsdf_rcp26,taxonkey,taxonName,maptype="rcp26_conf")
rcp45.conf.map<-confidenceMaps(pvalsdf_rcp45,taxonkey,taxonName,maptype="rcp45_conf")
rcp85.conf.map<-confidenceMaps(pvalsdf_rcp85,taxonkey,taxonName,maptype="rcp85_conf")

#End conformal prediction
###############################################################







