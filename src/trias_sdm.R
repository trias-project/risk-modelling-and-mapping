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
######GLOBAL SDM FOR WEIGHTING PSEUDOABSENCES IN EUROPEAN-LEVEL SDMs

setwd("C:../risk_modelling")
#all paths are relative to the risk_modelling folder


###read in global download for species from GBIF  
gbif_filename<- "Cyperus eragrostis.csv"
global<-read.csv(file=paste("./data/external/PRA_Plants/gbif_speciesFiles/",gbif_filename, sep=""))


 taxonkey<-"2715482" #GBIF taxonKey for the species being modeled.
 taxonName<-"Cyperus eragrostis" #names of species being modeling


##filter data to keep only those points with acceptable spatial accuracy (for us this 4 decimal places for either lat or lon) 
#script to count number of decimal places
decimalplaces <- function(x) {
  if (abs(x - round(x)) > .Machine$double.eps^0.5) {
    nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed = TRUE)[[1]][[2]])
  } else {
    return(0)
  }
}
global.occ<-global %>%
 #filter(taxonKey==taxonkey) %>%   #using taxonKey filters out accepted synonyms
  filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters< 709)

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


##create rasterstack of CHELSA data clipped to European modeling extent for prediction

#import eu_climpreds (global chelsa rasters clipped to European modeling extent for prediction
euclimrasters <- list.files("C:./data/external/climate/chelsa_eu_clips",pattern='tif',full.names = T)
eu_climpreds<-stack(euclimrasters)

#function for correcting global clim preds values
divide10<-function(x){
  value<-x/10
  return(value)
}

eu_chelsapreds<-divide10(eu_climpreds) #it is better to correct them on the fly than to store them corrected

#import WWF ecoregion layer clipped to distribution of target species for restricting pseudoabsence selection
wwf_eco<-shapefile("C:./data/external/GIS/official/wwf_terr_ecos.shp")

#select wwf ecoregions that contain occurrence points
crs(global.occ)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
occ_ecoIntersect <- over(wwf_eco,global.occ) 
wwf_ecoSub1 <- wwf_eco[!is.na(occ_ecoIntersect),]
ext_wwf_ecoSub<-extent(wwf_ecoSub1)

##import global ecoregions raster that has vascular plants counts of less than 6 removed
ecoregions_filtered<-raster("C:./data/external/GIS/official/wwf_terr_ecos_plants.tif")
#then subset ecoregions containing occurrence points
eco_filt_crop<-crop(ecoregions_filtered,ext_wwf_ecoSub)
ecoregions_filtered_sub<-mask(eco_filt_crop,wwf_ecoSub1)
#plot(ecoregions_filtered_sub) #not run



#extract long and lat for SDMtable
global.occ.LL<-data.frame(global.occ)[c(1:2)] 
##use SDMtab command from the SDMPlay package to remove duplicates per grid cell  
global.SDMtable<- SDMPlay:::SDMtab(global.occ.LL, globalclimpreds, unique.data = TRUE,background.nb= 0) #
numb.pseudoabs <- length(global.SDMtable$id) #select the same number of pseudoabsences as unique presences


##transform filtered occurrence dataset with unique presences back to a SpatialPoints dataframe. 
global.occ.sp<-global.SDMtable[c("longitude", "latitude")]
coordinates(global.occ.sp)<- c("longitude", "latitude")
global.occ.sp$species<- rep(1,length(global.occ.sp$latitude)) #adds columns indicating species presence needed for modeling

#use randomPoints function from dismo package to locate pseduobasences within the ecoregions
library(dismo)
set.seed(728)
global_points<-randomPoints(ecoregions_filtered_sub, numb.pseudoabs, global.occ.sp, ext=NULL, extf=1.1, excludep=TRUE, prob=FALSE,
                          cellnumbers=FALSE, tryf=50, warn=2, lonlatCorrection=TRUE) #will throw a warning if randomPoints generated is less than numb.pseudoabs. If this happens, increase the number of tryf.

global_pseudoAbs<-as.data.frame(global_points)
coordinates(global_pseudoAbs)<-c("x","y")
global_pseudoAbs$species<-rep(0,length(global_pseudoAbs$x))

#join filtered global presences with pseudo absences

library(maptools)
global_presabs<- spRbind(global.occ.sp,global_pseudoAbs)
writeOGR(global_presabs, dsn="C:./data/results/sdm_occ_pts/global", layer=paste(taxonName,"_globalSDM"), driver = "ESRI Shapefile", overwrite_layer=TRUE)

#####################  
###START MODELLING###
#####################

#use the sdmData to extract data for each occurrence across all the climate variables
global.data <- sdmData(species~.,train=global_presabs, predictors=globalclimpreds)
                     
global.data

global.data1<-as.data.frame(global.data)

####automatic predictor selection########

#Chelsa predictors need to be divided by 10 per README
# rescale predictors in dataframe and drop ID and response to calculate correlation matrix
global.data2<-(sapply(global.data1[,-c(1:2),], function(x) divide10(x)))
correlationMatrix <- cor(global.data2)

# summarize the correlation matrix
#print(correlationMatrix)
# find attributes that are highly corrected 
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.7,exact=TRUE,names=TRUE)
# print names of highly correlated attributes
print(highlyCorrelated)


#scale global clim preds to actual values to match model
globalclimpreds10<-divide10(globalclimpreds)
globalclimpreds_sub<-dropLayer(globalclimpreds10,highlyCorrelated)

#refit model with selected subset of uncorrelated climate variables

global.data.sub <- sdmData(species~.,train=global_presabs, predictors=globalclimpreds_sub)
global.data.sub

global.data.sub.df<-as.data.frame(global.data.sub)
global.data.sub.df<-within(global.data.sub.df,rm("rID"))
global.data.sub.df$species<-as.factor(global.data.sub.df$species)
levels(global.data.sub.df$species)<-c("absent","present")
global.data.sub.df$species <- relevel(global.data.sub.df$species, ref = "present")

control <- trainControl(method="repeatedcv",number=5, repeats=10, savePredictions="final", preProc=c("center","scale"),classProbs=TRUE)
classList1 <- c("glm","gbm","rf","adaboost","knn")
set.seed(457)
global_train <- caretList(
  species~., data=global.data.sub.df,
  trControl=control,
  methodList=classList1
)

modelResults<-resamples(global_train)
summary(modelResults)
modelCor(resamples(global_train))

set.seed(478)
global_stack <- caretStack(
  global_train, 
  method="glm",
  trControl=trainControl(method="cv",
                         number=10,
                         savePredictions= "final"  ))
print(global_stack)

#use global model to predict only for the extent of Europe
GlobalEns_pred_eu<-raster::predict(eu_chelsapreds,global_stack,type="prob")
writeRaster(GlobalEns_pred_eu, filename=paste("GlobalEnsEU_",taxonkey, ".tif",sep=""), format="GTiff",overwrite=TRUE) 

##########################END GLOBAL MODEL####################################################################################################################

############################################
#########EUROPEAN LEVEL MODEL##############
###########################################

#obtain european presences from the data cube if available, otherwise skip to ##Create European subset
#Read in european data cube file
# cube_europe <- read.csv("C:./data/external/data_cube/cube_europe.csv")
# 
# #extract occurrence data for the target species meeting our criteria
# 
# occ.eu<-cube_europe %>% 
#   filter(taxonKey==taxonkey) %>% #look up taxonKey by species in "cube_europe_taxa.txt" %>%
#   filter(year >="1975") %>% #select presences 1975 or later to be consistent with climate data
#   filter(min_coord_uncertainty <= 708) 
# 
# occ.eu$eea_cell_code<-str_remove(occ.eu$eea_cell_code,"1km")#so that cell code matches with chelsa data

#####create European subset if data cube not available
euboundary<-shapefile("C:/Users/amyjs/Documents/projects/Trias/modeling/GIS/EUROPE.shp")
ext<-extent(euboundary)


crs(global.occ)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
occ_euIntersect <- over(global.occ,euboundary,returnList=TRUE) 
occ.eu  <- global.occ[!is.na(occ_euIntersect),]
#remove commenting to see plots
#plot(euboundary)
#plot(occ.eu,add=TRUE)






###join lat lons from centroids to eu level occurrence data file
##skip if Trias EU occurrence cube data not available
#centroids <- shapefile("./data/external/GIS/EEA_fullgrid_1km_final_centroids.shp")
#centroids1<-cbind(centroids[1:3])


# merged<-merge(occ.eu,centroids1,by.x="eea_cell_code",by.y="EEA_fullgr") #1000 obs
# str(merged)
# 
# occ.eu1<-cbind(merged[1],merged[7:8])

# use this to create spatial points dataframe to extract climate data from eu climate rasters (any location within the same 1km eea grid cell of the 
# data will have the same climate data)



###PREDICTOR VARIABLES USING EUROPEAN LEVEL CLIMATE DATA (RMI)####
#1. create RasterStack of european climate variables
rmiclimrasters <- list.files("C:./data/external/climate/rmi_corrected",pattern='tif',full.names = T) #insert path where all climate rasters are located

rmiclimpreds <- stack(rmiclimrasters)
###Create SpatialPoints dataframe needed for SDMtab command if using occurrence cube
# euocc<-occ.eu1[c("longitude", "latitude")]#extract long and lat for SDMtable
# coordinates(euocc)<- c("longitude", "latitude")
# euocc1<-data.frame(euocc)[c(1:2)] 
#########################################################


###Create SpatialPoints dataframe needed for SDMtab command if NOT using occurrence cube

occ.eu1<-spTransform(occ.eu,crs(rmiclimpreds))
euocc<-occ.eu1@coords
euocc1<-data.frame(euocc)[c(1:2)] 
names(euocc1)<-c("longitude", "latitude")

#########EU SDM modeling
#use SDMtab command from the SDMPlay package to remove duplicates per grid cell  
euocc.SDMtable<- SDMPlay:::SDMtab(euocc1, rmiclimpreds, unique.data = TRUE,background.nb= 0)
numb.pseudoabs <- length(euocc.SDMtable$id) 



#transform occurrence dataset with unique presences back to a SpatialPoints dataframe. 
euocc<-euocc.SDMtable[c("longitude", "latitude")]
coordinates(euocc)<- c("longitude", "latitude")
euocc$occ<- rep(1,length(euocc$latitude))#adds columns indicating species presence needed for modeling
crs(euocc)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs ")

set.seed(728)
#######SELECT RANDOM PSEUDO ABSENCES from low sampled areas
#clip ecoregions file to european extent
studyextent<-euboundary
ecoregions_eu<-crop(ecoregions_filtered_sub,studyextent)
euregions_eu1<-projectRaster(ecoregions_eu,rmiclimpreds)

###############################################################################
###MASK OUT AREAS OF LOW HABITAT SUITABILITY FROM GLOBAL CLIMATE MODEL#####
##Read in global raster from earlier step if needed
#globalfilename<-paste("C:./GlobalEnsEU_", taxonkey, ".tif",sep="")
#global_model<-raster(globalfilename)
global_model<-GlobalEns_pred_eu

wgs84_gcs<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
crs(global_model)<-wgs84_gcs

m<-global_model >.5
global_mask<-raster:::mask(global_model,m,maskvalue=TRUE)
global_masked_proj<-projectRaster(global_mask,euregions_eu1)

##OVERLAY GLOBAL HABITAT SUITABILITY WITH FILTERED ECOREGIONS
pseudoSamplingArea<-mask(euregions_eu1,global_masked_proj)
plot(pseudoSamplingArea)

library(dismo)
euocc_points_res<-randomPoints(pseudoSamplingArea, numb.pseudoabs, euocc, ext=NULL, extf=1.1, excludep=TRUE, prob=FALSE,
                               cellnumbers=FALSE, tryf=25, warn=2, lonlatCorrection=TRUE) 

euocc_pseudoAbs_res<-as.data.frame(euocc_points_res)
coordinates(euocc_pseudoAbs_res)<-c("x","y")
euocc_pseudoAbs_res$occ<-rep(0,length(euocc_pseudoAbs_res$x))

crs(euocc_pseudoAbs_res)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs ")
#join filtered occ presences with pseudo absences 
library(maptools)
eu_presabs_res<- spRbind(euocc,euocc_pseudoAbs_res)

crs(eu_presabs_res)<-CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,-0,-0,-0,0 +units=m +no_defs ")
#confirm equal number of presences and absences
table(eu_presabs_res@data$occ)

#export occ points as shapefile for later viewing
writeOGR(obj=eu_presabs_res, dsn="C:./data/results/sdm_occ_pts/europe", layer=paste(taxonName,"_EUpresabs",sep=""), driver="ESRI Shapefile", overwrite_layer = TRUE) 


#create SDM data object 
occeu.data.res <- sdmData(occ~.,train=eu_presabs_res, predictors=rmiclimpreds) 
occeu.data.res

##########################################################################################
#################
####identify and filter out highly correlated predictors########
library(caret)

#convert eu data to dataframe
occeu.data.res.df<-as.data.frame(occeu.data.res)

#identify highly correlated predictors drop first two columns which always be rID and occ
correlationMatrix <- cor(occeu.data.res.df[,-c(1:2)])

# summarize the correlation matrix
print(correlationMatrix)
# find attributes that are highly corrected 
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.70,exact=TRUE,names=TRUE)
# print indexes of highly correlated attributes
print(highlyCorrelated)
#random forests
 control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# # run the RFE algorithm to select predictors
 results <- rfe(occeu.data.res.df[,3:15], as.factor(occeu.data.res.df[,2]), sizes=c(1:13), rfeControl=control)
# # summarize the results
 print(results)
 # list the chosen features
 predictors(results)
## plot the results
# plot(results, type=c("g", "o")) # not run

#remove highly correlated predictors from dataframe and rasterstack
occeu.data.df1.uncor<-select (occeu.data.res.df[,-1],-c(highlyCorrelated))
rmiclimpreds_uncor<-dropLayer(rmiclimpreds,highlyCorrelated)

##########################################################################################################################

#investigate role of habitat and anthropogenic factors on predicted risk
habitat<-list.files("C:./data/external/habitat",pattern='tif',full.names = T)
habitat_stack<-stack(habitat)


fullstack<-stack(rmiclimpreds_uncor,habitat_stack) #combine uncorrelated climate variable selected earlier with habitat

##clip fullstack to belgium extent
belgie<-shapefile("C:./data/external/GIS/belgium_boundary.shp")
habitatstack1_be<-crop(fullstack,belgie)
fullstack_be<-mask(habitatstack1_be,belgie)

occ.full.data <- sdmData(occ~.,train=eu_presabs_res, predictors=fullstack) 
occ.full.data



###convert eu data to dataframe
occeu.full.data.df<-as.data.frame(occ.full.data)
occeu.full.data.df<-occeu.full.data.df[,-1]


###identify highly correlated predictors from the combined stack of habitat/anthropogenic and climate

correlationMatrix <- cor(occeu.full.data.df[,-1])
# summarize the correlation matrix
print(head(correlationMatrix))
# find attributes that are highly corrected 
highlyCorrelated <- findCorrelation(correlationMatrix, cutoff=0.8,exact=TRUE,names=TRUE)
# print names of highly correlated attributes
print(highlyCorrelated)


####Identify best subset of predictors for each algorithm
# #random forests
# control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# # run the RFE algorithm to select predictors
# endColumn<-length(occeu.full.data.df)
# results <- rfe(occeu.full.data.df[,-1], as.factor(occeu.full.data.df[,1]), sizes=c(1:endColumn), rfeControl=control)

# # summarize the results
# print(results)
# # list the chosen features
# predictors(results)
# # plot the results
# plot(results, type=c("g", "o"))

# rfe_rank<-function(object, x, y) {
#   vimp <- varImp(object)
#   vimp <- vimp[order(vimp$Overall,decreasing = TRUE),,drop = FALSE]
#   vimp$var <- rownames(vimp)                  
#   vimp
# }
# rfe_rank(results,occeu.full.data.df[,-1],as.factor(occeu.full.data.df[,1]))


occeu.full.data.df1<-select (occeu.full.data.df,-c(highlyCorrelated))
nzv_preds<-nearZeroVar(occeu.full.data.df1,names=TRUE)
occeu.full.data.df1<-select (occeu.full.data.df1,-c(nzv_preds))

##now build models with climate and habitat data
fulldata<-occeu.full.data.df1
fulldata$occ<-as.factor(fulldata$occ)
levels(fulldata$occ)<-c("absent","present")

control <- trainControl(method="cv",number=10,savePredictions="final", preProc=c("center","scale"),classProbs=TRUE)
mylist<-list(
  glm =caretModelSpec(method = "glm",maxit=50),
  gbm= caretModelSpec(method = "gbm"),
  rf = caretModelSpec(method = "rf", importance = TRUE),
  adaboost =caretModelSpec(method = "adaboost"),
knn = caretModelSpec(method = "knn"))
set.seed(457)
model_train_habitat <- caretList(
  occ~., data=fulldata,
  trControl=control,
  tuneList=mylist
)


modelResults1<-resamples(model_train_habitat)
summary(modelResults1)
modelCor(resamples(model_train_habitat))

set.seed(478)
lm_stack <- caretStack(
  model_train_habitat, 
  method="glm",
  na.action=na.pass,
  trControl=trainControl(method="cv",
                         number=10,
                         savePredictions= "final",classProbs = TRUE))
lm_stack

##does caretEnsemble give the same results as caretStack?
set.seed(478)
lm_ens_hab<-caretEnsemble(model_train_habitat, trControl=trainControl(method="cv",
                                                                  number=10,
                                                                  savePredictions= "final",classProbs = TRUE))

lm_ens_hab
variableImportance<-varImp(lm_ens_hab)
write.csv(variableImportance,file = paste(taxonkey,".csv"))
#the answer is yes
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
##predict for europe

ens_pred_hab_eu<-raster::predict(fullstack,lm_ens_hab,type="prob")
ens_pred_hab_eu1<-1-ens_pred_hab_eu
writeRaster(ens_pred_hab_eu1,filename=paste("eu_",taxonkey, "_hist.tif",sep="") , format="GTiff",overwrite=TRUE) 

#######################
##predict for Belgium only
#etrs89<-CRS("+proj=longlat +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +no_defs" )
laea_grs80<-crs(belgie)
ens_pred_hab<-raster::predict(fullstack_be,lm_ens_hab,type="prob")
ens_pred_hab1<-1-ens_pred_hab
crs(ens_pred_hab1)<-laea_grs80
writeRaster(ens_pred_hab1, filename=paste("be_",taxonkey, "_hist.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(ens_pred_hab1,taxonkey,taxonName=taxonName,"hist.pdf")

##create habitat-climate stacks for each rcp scenario
habitat_only_stack<-crop(habitat_stack,belgie)
habitat_only_stack_be<-crop(habitat_only_stack,belgie)

be26 <- list.files("C:./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26",pattern='tif',full.names = T)
belgium_stack26 <- stack(be26)


be45 <- list.files("C:./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45",pattern='tif',full.names = T)
belgium_stack45 <- stack(be45)

be85 <- list.files("C:./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85",pattern='tif',full.names = T)
belgium_stack85 <- stack(be85)


fullstack26<-stack(be26,habitat_only_stack_be)
fullstack45<-stack(be45,habitat_only_stack_be)
fullstack85<-stack(be85,habitat_only_stack_be)

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

##calculate difference between predicted risk by each RCP scenario and historical climate

hist26_diff_hab<-overlay(ens_pred_hab26_1, ens_pred_hab1, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist26_diff_hab,filename=paste("be_",taxonkey, "_rcp26_diff.tif",sep="") , format="GTiff",overwrite=TRUE) 
exportPDF(hist26_diff_hab,taxonkey,taxonName=taxonName,"rcp26_diff.pdf","TRUE")


hist45_diff_hab<-overlay(ens_pred_hab45_1, ens_pred_hab1, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist45_diff_hab,filename=paste("be_",taxonkey, "_rcp45_diff.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(hist45_diff_hab,taxonkey,taxonName=taxonName,"rcp45_diff.pdf","TRUE")

hist85_diff_hab<-overlay(ens_pred_hab85_1, ens_pred_hab1, fun=function(r1,r2){return(r1-r2)})
writeRaster(hist85_diff_hab, filename=paste("be_",taxonkey, "_rcp_85_diff.tif",sep=""), format="GTiff",overwrite=TRUE) 
exportPDF(hist85_diff_hab,taxonkey,taxonName=taxonName,"rcp85_diff.pdf","TRUE")

#output data used for SDM
df4export<-cbind(fulldata,coordinates(occ.full.data))
write.csv(df4export,file=paste(taxonkey,"_sdmdata.csv",sep=""))


##check spatial autocorrelation of residuals
library(geoR)
#derive residuals from unthinned model
predEns1<-lm_ens_hab$ens_model$pred
obs.numeric<-ifelse(predEns1$obs == "absent",0,1)

#derive sensitivity and specificity
table(predEns1$pred,predEns1$obs)
sensitivity(predEns1$pred,predEns1$obs)
specificity(predEns1$pred,predEns1$obs)


#######################################################################################################
#standardize the residuals: generalized residuals/sq root of the variance according to Davis et al 2015
stdres<-function(obs.numeric, yhat){
  num<-obs.numeric-yhat
  denom<-sqrt(yhat*(1-yhat))
  return(num/denom)
}
hab.res<-stdres(obs.numeric,predEns1$present)

library(geoR)
res.best.coords<-cbind(coordinates(occ.full.data),hab.res)
res.best.geo<-as.geodata(res.best.coords,coords.col=1:2,data.col = 3)
summary(res.best.geo) #note distance is in meters

res.best.vario<-variog(res.best.geo,coords=res.best.geo$coords, data=res.best.geo$data,max.dist=50000,option = "bin")
plot(res.best.vario)


#check MoransI
library(ape)
res.best.df<-as.data.frame(res.best.coords)
occ.dists <- as.matrix(dist(cbind(res.best.df[1], res.best.df[2])))
occ.dists.inv <- 1/occ.dists
diag(occ.dists.inv) <- 0
Moran.I(res.best.df$hab.res,occ.dists.inv,scaled=TRUE,alternative="greater")
##if Moran's I is very low, or not significant, skip thinning.

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
ens_calib<-lm_ens_hab$ens_model$pred

calibPresence<-subset(ens_calib$present,ens_calib$obs=='present',)
calibAbsence<-subset(ens_calib$absent,ens_calib$obs=='absent',)

#extract values from raster stack for use in conformal prediction algorithm in a way that they can be linked back to the belgium raster (ie include spatial info)

conf.hist<-extractVals(ens_pred_hab)

testPresence<-conf.hist$presence
testAbsence<-conf.hist$absence


# run CP scripts
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
pvalsdf_hist<-cbind(pvalsdf,conf.hist)

# #export csv 
write.csv(pvalsdf_hist,file=paste("confidence_",taxonkey, "_hist.csv",sep=""))


##conformal prediction for RCP scenarios

conf.hist<-extractVals(ens_pred_hab26)

testPresence<-conf.hist$presence
testAbsence<-conf.hist$absence

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
pvalsdf_rcp26<-cbind(pvalsdf,conf.hist)

# #export csv 
write.csv(pvalsdf_rcp26,file=paste("confidence_",taxonkey, "_rcp26.csv",sep=""))


conf.hist<-extractVals(ens_pred_hab45)

testPresence<-conf.hist$presence
testAbsence<-conf.hist$absence

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
pvalsdf_rcp45<-cbind(pvalsdf,conf.hist)

# #export csv 

write.csv(pvalsdf_rcp45,file=paste("confidence",taxonkey,"_rcp45.csv",sep=""))

conf.hist<-extractVals(ens_pred_hab85)

testPresence<-conf.hist$presence
testAbsence<-conf.hist$absence

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
pvalsdf_rcp85<-cbind(pvalsdf,conf.hist)

# #export csv 
write.csv(pvalsdf_rcp85,file=paste("confidence",taxonkey,"_rcp85.csv",sep=""))

###Create and export confidence maps

data.xyz.hist <- pvalsdf_hist[c(5,6,4)]
rst.hist <- rasterFromXYZ(data.xyz.hist)
crs(rst.hist)<-laea_grs80
plot(rst.hist)
writeRaster(rst.hist, filename=paste("be_",taxonkey, "_hist_conf.tif",sep=""), format="GTiff",overwrite=TRUE)
exportPDF(rst.hist,taxonkey,taxonName=taxonName,"hist_conf.pdf")

data.xyz.26 <- pvalsdf_rcp26[c(5,6,4)]
rst.26 <- rasterFromXYZ(data.xyz.26)
crs(rst.26)<-laea_grs80
plot(rst.26)
writeRaster(rst.26, filename=paste("be_",taxonkey, "_rcp26_conf.tif",sep=""), format="GTiff",overwrite=TRUE)
exportPDF(rst.26,taxonkey,taxonName=taxonName,"rcp26_conf.pdf")


data.xyz.45 <- pvalsdf_rcp45[c(5,6,4)]
rst.45 <- rasterFromXYZ(data.xyz.45)
crs(rst.45)<-laea_grs80
writeRaster(rst.45, filename=paste("be_",taxonkey, "_rcp45_conf.tif",sep=""), format="GTiff",overwrite=TRUE)
exportPDF(rst.45,taxonkey,taxonName=taxonName,"rcp45_conf.pdf")

data.xyz.85 <- pvalsdf_rcp85[c(5,6,4)]
rst.85 <- rasterFromXYZ(data.xyz.85)
crs(rst.85)<-laea_grs80
writeRaster(rst.85, filename=paste("be_",taxonkey, "_rcp85_conf.tif",sep=""), format="GTiff",overwrite=TRUE)
exportPDF(rst.85,taxonkey,taxonName=taxonName,"rcp85_conf.pdf")
#End conformal prediction
###############################################################


rm(list = ls())
gc()

removeTmpFiles()






