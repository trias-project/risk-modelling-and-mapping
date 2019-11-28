# Alien species risk modeling and mapping

This respository contains the framework and R code in development for risk modelling and mapping of alien species throughout Belgium and greater Europe at 1 km2 resolution as part of the TrIAS project.


### Repository Structure

    ├── LICENSE            <- MIT.
    ├── README.md          <- Description of this repository.
    ├── src                <- Source code for use in this project.
    ├── data               <- data for modelling
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    ├── maps               <- Risk maps for Belgium (PDF).
    
      
 ## R Workflow  
  ### A.) Develop global scale climate-only species distribution models (SDMs)
  These SDMs use all acceptable species occurrence data available at the date of download from GBIF to create climate suitability maps with global coverage for weighing pseudo absences. 
   1. Download georeferenced occurrence data for the target species from GBIF with as many filters as possible to decrease the load on GBIF
    2. Further filter data, by extracting points that match the time period being modeled and with the minimal acceptable geographic accuracy.
    3. Using the raster package, create a rasterstack of the global climate layers (e.g. Annual Temperature, Annual Precipitation etc.) We used data from CHELSA ()
    use SDMtab command from the SDMPlay package to remove duplicates per grid cell. 
    4. Next,import WWF ecoregions layer (https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world) clipped to the distribution of target species. This step restricting pseudoabsence selection to those areas that theoretically reachable by the target species. 
    5.Use randomPoints function from dismo package to randomly locate pseduobasences within the ecoregions; join these pseuboabsences to the presence data created in step 2.
    6. Run multiple SDMs (using different algorithms) using the sdm package
    7. Combine the SDMs by creating a weighted ensemble using a model evaluation statistic (we used TSS). Only predict to the limit of your study extent to reduce computation time.(Our study extent is Europe)
    8. Evaluate accuracy of the model using 10-fold cross validation
  
  ### B) Generate European level SDMs
  These SDMs predict the risk of invasion by alien species at 1km2 spatial resolution (using the EEA 1km2 grid) and the European data cube. It is from these models that the risk maps for Belgium are extracted.
  
    Data Prep
  1. Extract occurrence data that is Acer negundo and meets our criteria from the European Data Cube
  2. Using the raster package, create a rasterstack of the european historical climate layers provided by RMI (these will be publicly available for download. Also create a rasterstack clipped to Belgium for prediction.
  
   Restrict pseudoabsences to areas of low habitat suitability that have been adequately sampled
   
   3. Using the raster package, mask areas of high habitability of the global climate raster to create a low habitat suitability global layer
   4. Using 1 degree grid of record counts by target species group (e.g if modelling a tree species, use the Vascular Plants target group grid; the target group grid summarizes all records in GBIF off a broad taxonomic group such as Vascular Plants, Birds,etc. by 1 degree grid cell. These will be publicly available for download). Use this grid to exclude areas that have received little to no sampling effort from the global climate layer to reduce effects of sampling bias. Use this raster to randomly locate pseudoabsences.
   5.Join pseudoabsence to acer presence data created above.
   6. Run multiple SDMs(using different algorithms),create weighted ensembles that limit prediction to the study extent (e.g Belgium). 
   7. Evaluate accuracy of model. Investigate addition of biophysical,landcover, and anthropogenic variables to model accuracy.
   
  ###  C) Forecast species distributions under climate change scenarios
   8. Build raster stacks of the same variables used to calibrate the European level model, replacing historical climate predictors with the future climate predictors by RCP scenario provided by RMI.Clip rasters to study extent.
   9. Choose the best performing European level model(s) from step 7 and use it to forecast species distributions using the raster stack created in 8.
  
 
      


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
