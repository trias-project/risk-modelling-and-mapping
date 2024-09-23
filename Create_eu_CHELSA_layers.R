#Crop global CHELSA climate layers (terra format) to extent of Europe (as in paper)
europe_ext<-ext(-12, 31.6, 35, 71) #xmin, xmax, ymin, ymax
euclimpreds<-crop(globalclimpreds_terra, europe_ext)
plot(euclimpreds[[1]])

#Divide all raster values by 10
eu_climpreds_terra.10<-divide10(euclimpreds)

#Predict model 
system.time({
  global_model<-predict_large_raster(eu_climpreds_terra.10,global_stack,type="prob")
})
