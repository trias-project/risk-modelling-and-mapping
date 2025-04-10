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

testeval.eu.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_eu_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.eu.bin.rast

binary_be_rasters<-sapply(names(thresholds), function(x) raster::reclassify(ens_pred_hab_be[[x]],c(0,thresholds[[x]]$predicted,0, thresholds[[x]]$predicted,1,1)),simplify=FALSE)
testeval.be.bin.rast<-sapply(names(binary_eu_rasters), function(x) eu_eval(binary_be_rasters[[x]],eval.data.occ.proj),simplify=FALSE)
testeval.be.bin.rast
