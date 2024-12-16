#--------------------------------------------
#-----------To do: specify project-----------
#--------------------------------------------
#specify project name
projectname<-"Test_Frédérique"


#--------------------------------------------
#-----------  Load packages  ----------------
#--------------------------------------------
packages <- c("viridis", "dplyr", "grid", "here", "qs","terra", "sf", "ggplot2","RColorBrewer","magick","patchwork"
)

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


#--------------------------------------------
#----------- Load taxa info  ----------------
#--------------------------------------------
taxa_info<-read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))
accepted_taxonkeys<-taxa_info%>%
  pull(speciesKey)%>%
  unique()


#--------------------------------------------
#-----------Load country data----------------
#--------------------------------------------
#If you'd like to predict for another country, change the shapefile
country<-st_read(here("./data/external/GIS/Belgium/belgium_boundary.shp"))
country_ext<-terra::ext(country) 
country_vector <- terra::vect(country) #Convert country to a SpatVector that can be used for masking


#--------------------------------------------
#--------Source helper functions-------------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#-----------  Start loop   ----------------
#--------------------------------------------
for(key in accepted_taxonkeys){
  #Extract species name
  species<-taxa_info%>%
    filter(accepted_taxonkeys==key)%>%
    pull(scientificName)%>%
    unique()
  
  #Extract first two words of species name
  first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
  
  #Define taxonkey
  taxonkey<- key
  
  #Read in globalmodels object that was stored as part of  script 03_fit_European_model
  eumodel<-qread( paste0("./data/projects/",projectname,"/",first_two_words,"_",taxonkey,"/EU_model_",first_two_words,"_",taxonkey,".qs"))
  
  #Read in different data objects stored in globalmodels
  euocc<-eumodel$euocc1
  bestModel<-unwrap(eumodel$bestModel)
  fullstack_be<-unwrap(eumodel$fullstack_be)
  
  
### Subset Belgium occurrences 
#occ.eu is in WGS84, convert to same projection as country level shapefile (which is the same proj used for model outputs)
occ.eu.proj  <- spTransform(occ.eu,crs(country))
occ.country <- occ.eu.proj[country,]
plot(country)
plot(occ.country,pch=21,bg="green",cex=1,add=TRUE)

### plot the best EU level ensemble model showing only Belgium

brks <- seq(0, 1, by=0.1) 
nb <- length(brks)-1 
pal <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
cols<-pal(nb)
plot(ens_pred_hab_be$X6, breaks=brks, col=cols,lab.breaks=brks) # specify best model
plot(occ.country,pch=21,cex=1,add=TRUE)




### Clip habitat raster stack to Belgium 
habitat_stack<-stack(habitat)
habitat_only_stack<-crop(habitat_stack,country)
habitat_only_stack_be<-crop(habitat_only_stack,country)

### Create individual RCP (2.6, 4.5, 8.5) climate raster stacks for Belgium
be26 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp26")),pattern='tif',full.names = T)
belgium_stack26 <- stack(be26)

be45 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp45")),pattern='tif',full.names = T)
belgium_stack45 <- stack(be45)

be85 <- list.files((here("./data/external/climate/byEEA_finalRCP/belgium_rcps/rcp85")),pattern='tif',full.names = T)
belgium_stack85 <- stack(be85)


### Combine habitat stacks with climate stacks for each RCP scenario
fullstack26<-stack(be26,habitat_only_stack_be)
fullstack45<-stack(be45,habitat_only_stack_be)
fullstack85<-stack(be85,habitat_only_stack_be)


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


testeval.eu.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_eu_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.eu.bin.rast

binary_be_rasters<-sapply(names(thresholds), function(x) raster::reclassify(ens_pred_hab_be[[x]],c(0,thresholds[[x]]$predicted,0, thresholds[[x]]$predicted,1,1)),simplify=FALSE)
testeval.be.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_be_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.be.bin.rast
