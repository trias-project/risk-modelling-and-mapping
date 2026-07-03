#--------------------------------------------
#-----------    Load packages      ----------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")
terra::setGDALconfig("GDAL_PAM_ENABLED", "FALSE")#Prevent terra from writing aux.xml files

packages <- c( "viridis","dplyr", "here", "qs2", "tidyterra","sf", "ggplot2","RColorBrewer","magick","patchwork","grid", "randomForest", "progressr", "raster", "dismo", "caret", "caretEnsemble", "kableExtra","gbm", "PresenceAbsence", "RStoolbox", "sdm", "purrr", "terra")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}

sdm::installAll()


#--------------------------------------------
#- Load helper functions and configurations -
#--------------------------------------------
source(file.path("src", "helper_functions.R"))
source(file.path("src", "00_configurations.R"))


#--------------------------------------------
#---- Define habitat raster file paths ------
#--------------------------------------------
processed_folder<-file.path("data", "external", "habitat", "processed")
habitatstack_file <- file.path(processed_folder, "habitat_stack.tif")


#--------------------------------------------
#---------   Load euboundary  ---------
#--------------------------------------------
euboundary_path<-file.path("data", "external", "GIS", "Europe", "EUboundary.shp")
euboundary<-sf::st_read(euboundary_path)


#---------------------------------------------
#----- Load country boundary ------
#---------------------------------------------
habitat_stack<-terra::rast(habitatstack_file)
if(tolower(country_of_interest)!="europe"||!is.null(custom_country_boundary_path)){
  country_boundary <- sf::read_sf(here::here("data","external","GIS","Country","country.shp"))%>%
    sf::st_transform(crs(habitat_stack[[1]]))%>%
    terra::vect()
}else{
  country_boundary<-euboundary%>%
    terra::vect()
}


#--------------------------------------------
#------------- Load species data -----------
#--------------------------------------------
taxa_info <- read.csv2(paste0("./data/projects/",project,"/",project,"_taxa_info.csv"))
accepted_taxonkeys <- taxa_info %>%
  dplyr::pull(acceptedTaxonKey) %>%
  unique()


#--------------------------------------------
#------------- Clean up -----------
#--------------------------------------------
rm(habitat_stack)
gc()


#--------------------------------------------
#----------- Start modelling loop  ----------
#--------------------------------------------

with_progress({
  p <- progressr::progressor(along = 1:length(accepted_taxonkeys)) 
  
  for(key in accepted_taxonkeys){ 
    #--------------------------------------------
    #---------------Map progress  ---------------
    #--------------------------------------------
    p()
    start_time <- Sys.time()
    
    #--------------------------------------------
    #--------Extract species-specific data  -----
    #--------------------------------------------
    #Extract species name
    species <- taxa_info%>%
      dplyr::filter(acceptedTaxonKey==key)%>%
      dplyr::pull(acceptedScientificName)%>%
      unique()
    
    #Extract first two words of species name
    speciesName <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
    
    #Extract rest of species name
    nameExtension <- if (grepl("^\\S+\\s+\\S+$", species)) "" else sub("^\\S+\\s+\\S+\\s+", "", species)
    
    #Specify species for plot title
    species_title <- gsub("_", " ", speciesName)
    
    #Define taxonkey
    taxonkey<- key
    
    
    #---------------------------------------------
    #-- Prepare filenames and titles for export --
    #---------------------------------------------
    #Prepare PDF title 
    PDF_title<-bquote(italic(.(gsub("_", " ", speciesName))) ~ .(nameExtension) ~ "(" * .(taxonkey) * ")")
    
    #Prepare current and future basefile
    basefile<-  paste0(speciesName, "_Habitat_")
    combined_basefile<-  paste0(speciesName, "_Combined_")
    global_basefile<-  paste0(speciesName, "_Climate_")
    
    
    #--------------------------------------------
    #-- Define file path of global model file  --
    #--------------------------------------------  
    base_dir <- file.path("data", "projects", project, paste0(speciesName,"_",taxonkey))
    global_model_file <- file.path(base_dir,"Climate",
                                   paste0("Climate_model_",speciesName,"_",taxonkey,".qs"))
    
    
    #--------------------------------------------
    #-Check if global model exists, if not, skip-
    #--------------------------------------------
    if(file.exists(global_model_file)){
      
      #This was stored as part of  script 02_fit_global_model
      globalmodels <- qs2::qs_read(global_model_file)
      
      #Extract different data objects stored in globalmodels
      global.occ.sf <- globalmodels$occurrences1km %>% # FULL occurrence with coordinateUncertainty <= 1km
        sf::st_as_sf(.,coords = c("decimalLongitude", "decimalLatitude"), crs=4326)
      
    }else{
      warning(paste0("Skipping species ", species, " because no global model could be fitted"))
      next  # Skip the rest of the loop and move to the next iteration
    }
    
    
    #---------------------------------------------
    #------------- Define folders ----------------
    #---------------------------------------------
    raster_folder <- file.path(base_dir, "Habitat","Current", "Predictions", "Rasters")
    climate_raster_folder <- file.path( base_dir,"Climate", "Current", "Predictions", "Rasters")
    
    
    #--------------------------------------------
    #------------ Import raster layers ----------
    #--------------------------------------------
    #Define file paths
    biasgrid_file <- file.path(base_dir,"Climate", "Current", "Interim", paste0("Biasgrid_",speciesName,"_",taxonkey,".tif"))
    #Load rasterlayers
    habitat_stack <- terra::rast(habitatstack_file)
    biasgrid_sub <- terra::rast(biasgrid_file)
    
    
    #-------------------------------------------------
    #--------------- Create EU folders ---------------
    #-------------------------------------------------
    # Define outputs, periods, and scenarios
    periods   <- c("Current","2041-2070", "2071-2100")
    scenarios <- c("ssp126", "ssp370", "ssp585")
    outputs   <- c("Rasters", "PDFs", "PNGs")
    
    #Create folders for each combination
    scenario_folders <- list()
    for(period in periods){
      for(output in outputs){
        if(period=="Current"){
          loop_list <- list(list(path = file.path(base_dir, "Habitat", period,"Predictions",output),
                                 name = paste("Habitat", period, "Predictions", output,  sep = "/")),
                            list(path = file.path(base_dir, "Combined", period,"Predictions",output),
                                 name = paste("Combined", period,"Predictions", output,  sep = "/")),
                            list(path = file.path(base_dir, "Habitat", period,"Diagnostics", "Variable_importance"),
                                 name = paste("Habitat", period, "Diagnostics", "Variable_importance", output,  sep = "/")),
                            list(path = file.path(base_dir, "Habitat", period,"Diagnostics", "Response_curves"),
                                 name = paste("Habitat", period, "Diagnostics", "Response_curves", output,  sep = "/")),
                            list(path = file.path(base_dir, "Habitat", period,"Diagnostics", "Confidence_maps",output),
                                 name = paste("Habitat", period, "Diagnostics", "Confidence_maps", output,  sep = "/")),
                            list(path = file.path(base_dir, "Combined", period,"Diagnostics", "Confidence_maps",output),
                                 name = paste("Combined", period, "Diagnostics", "Confidence_maps", output,  sep = "/")))
          scenario_folders <- c(scenario_folders, loop_list)  
          
        }else{
          for(scenario in scenarios){
            loop_list <- list(list(path = file.path(base_dir, "Combined", period, scenario, "Predictions", output),
                                   name = paste("Combined", period, scenario, output, sep = "/")),
                              list(path = file.path(base_dir, "Combined", period, scenario, "Diagnostics", "Confidence_maps", output),
                                   name = paste("Combined", period, scenario,"Diagnostics", "Confidence_maps",  output, sep = "/")))
            scenario_folders <- c(scenario_folders, loop_list)
          }
        }
      }
    }
    
    # Check and create each folder if necessary
    lapply(scenario_folders, function(folder){
      create_folder(folder$path, folder$name)
    })
    
    
    #-----------------------------------------------
    #----- Create subset of European records -------
    #-----------------------------------------------
    #Check for occurrences that fall within Europe
    eu_occ <- global.occ.sf%>%
      st_transform(crs = st_crs(habitat_stack)) %>%
      sf::st_filter(euboundary) %>%
      sf::st_coordinates() %>%
      as.data.frame()
    rm(global.occ.sf)
    gc()
    
    
    #-----------------------------------------------
    #----------- Process occurrences ---------------
    #-----------------------------------------------
    #Keep only one occurrence per grid cell
    eu_occ <- remove_duplicates(occurrences =  eu_occ, rast_template = habitat_stack[[1]])%>%
      dplyr::rename("decimalLongitude" = X,
                    "decimalLatitude" = Y)
    
    #Remove occurrences within grid cells with NA values
    eu_occ <- remove_nodata_occurrences(occurrences = eu_occ, rast_template= habitat_stack[[1]], st_crs(habitat_stack))
    
    
    #-----------------------------------------------
    #------ Limit to 10,000 occupied grid cells ----
    #-----------------------------------------------
    if(nrow(eu_occ) > 10000){
      if(occurrence_thinning_method == "random"){
        print("Thinning occurrences randomly")
        set.seed(101) 
        eu_occ <- eu_occ[sample(nrow(eu_occ), 10000, replace=FALSE), ]
      }else if (occurrence_thinning_method == "kmeans_clustering"){
        print("Thinning occurrences based on k-means clustering")
        #Extract environmental data in each occurrence grid cell
        habitat_data <- terra::extract(habitat_stack, eu_occ, ID = FALSE)
        
        #Check how many unique rows there are and set centers to lowest of either 10000 or #unique rows
        unique_centers<-nrow(unique(habitat_data))
        center_number<-min(unique_centers, 10000)
        
        # K-means clustering
        set.seed(101)
        clust <- kmeans(habitat_data, centers = center_number,iter.max = 10, nstart = 1)$cluster
        occ_habitat <- cbind(eu_occ, habitat_data, clust)%>%
          dplyr::mutate(rID =row_number())
        
        # Keep 1 occurrence per cluster
        sampled <- occ_habitat %>%
          dplyr::group_by(clust) %>%
          dplyr::slice_sample(n = 1) %>%
          dplyr::ungroup()
        
        # How many presences do we still need
        remaining <- 10000 - nrow(sampled)
        
        # sample extra occurrences if fewer than 10000
        if (remaining > 0) {
          # Randomly sample additional presences excluding already chosen ones
          extra_occ <- occ_habitat %>%
            dplyr::filter(!rID %in% sampled$rID)%>%
            dplyr::slice_sample(n = remaining) 
          
          eu_occ <- bind_rows(sampled, extra_occ)
          rm(extra_occ)
          
        } else {
          eu_occ <- sampled
        }
        
        # Keep only occurrence columns
        eu_occ <- eu_occ %>%
          dplyr::select(decimalLongitude, decimalLatitude, geometry)
        
        rm(habitat_data, occ_habitat, sampled, remaining, unique_centers, center_number, clust)
        
      }
    }
    
    
    #------------------------------------------------
    #----- Check if at least 20 European records ----
    #------------------------------------------------
    if (nrow(eu_occ) < 20) {
      warning(paste(nrow(eu_occ)," occurrences in Europe for species:", species, 
                    "\n- European model cannot be constructed, skipping to the next species."))
      next  # Skip to the next species in the loop
    }
    
    
    #--------------------------------------------
    #-------- Plot European occurrences ---------
    #--------------------------------------------
    # ggplot()+
    # geom_sf(data = euboundary,  colour = "black", fill = NA)+
    # geom_point(data=eu_occ, aes(x=decimalLongitude, y= decimalLatitude),  fill="green", shape = 22, colour = "black", size=3)+
    # labs(x="Longitude", y="Latitude")+
    # theme_bw()
    
    
    #--------------------------------------------
    #----- Clip biasgrid to European extent -----
    #--------------------------------------------
    # Reproject biasgrid_sub to match CRS, extent, and resolution of habitat_stack
    biasgrid_aligned <- terra::project(
      biasgrid_sub,
      habitat_stack[[1]],
      method = "bilinear" 
    )
    
    #Mask biasgrid with habitat raster (so no PA can be selected in cells that are NA in habitat rasters)
    biasgrid_aligned <- terra::mask(biasgrid_aligned, habitat_stack[[1]])
    
    
    #-------------------------------------------
    #------- Select invaded WWF ecoregions------
    #-------------------------------------------
    #Load WWF ecoregions
    wwf_ecoregions <- sf::st_read(here("./data/external/GIS/official/wwf_terr_ecos.shp")) %>%
      sf::st_transform(crs = st_crs(habitat_stack))
    
    # Identify which polygons contain at least one occurrence
    polygons_with_points <- lengths(sf::st_intersects(wwf_ecoregions, eu_occ)) > 0
    
    # Subset only those polygons
    wwf_ecoregions_filtered <- wwf_ecoregions[polygons_with_points, ]
    # plot(wwf_ecoregions_filtered[4], key.pos = NULL)
    
    #Clean up ecoregions
    rm(wwf_ecoregions)
    gc()
    
    
    #----------------------------------------------------------------------------------------
    #---- biasgrid: keep values inside invaded ecoregions, set outside to 1 (lowest value)---
    #----------------------------------------------------------------------------------------
    # Step 1: Rasterize WWF polygons to match biasgrid_aligned
    inside_mask <- terra::rasterize(vect(wwf_ecoregions_filtered), biasgrid_aligned, field = 1, background = NA)
    
    # Step 2: Apply logic — keep original where inside_mask, else 1
    biasgrid_temp <- terra::ifel(!is.na(inside_mask), biasgrid_aligned, 1)
    
    # Step 3: Restore NA values from the original biasgrid
    biasgrid_eu <- mask(biasgrid_temp, biasgrid_aligned)
    
    
    #--------------------------------------------
    #----------- Visualize biasgrid  ------------
    #--------------------------------------------
    # ggplot()+
    # geom_sf(data = sf::st_transform(euboundary, crs=st_crs(habitat_stack)),  colour = "black", fill = NA)+
    # geom_spatraster(data=biasgrid_eu)+
    # scale_fill_continuous(na.value = "transparent",low = "blue", high = "orange")+
    # labs(x="Longitude", y="Latitude")+
    # theme_bw()
    
    
    #--------------------------------------------------------------
    #Generate pseudoabsences weighted by sampling bias---------
    #--------------------------------------------------------------
    set.seed(728)
    global_points <- terra::spatSample(
      biasgrid_eu,
      size = 30000, #three times the number we need
      method = "weights",     # weighted random sampling
      as.points = TRUE,       # return SpatVector of points
      na.rm = TRUE            # ignore NA pixels
    )
    
    #Select 10000 pseudoabsences
    if(nrow(global_points) > 10000){
      if(pseudoabsence_thinning_method == "random"){
        print("Thinning pseudoabsences randomly")
        set.seed(101) 
        global_points <- global_points[sample(nrow(global_points), 10000, replace=FALSE), ]%>%
          sf::st_as_sf()
        
        coords <- sf::st_coordinates(global_points)
        
        global_points<-global_points%>%
          dplyr::mutate(decimalLongitude = coords[, "X"],
                        decimalLatitude  = coords[, "Y"])%>%
          dplyr::select(decimalLongitude, decimalLatitude, geometry)
        
      }else if (pseudoabsence_thinning_method == "kmeans_clustering"){
        print("Thinning pseudoabsences based on k-means clustering")
        
        #Extract environmental data from pseudoabsences
        pa_habitat_data <- terra::extract(habitat_stack, global_points, ID = FALSE, xy = TRUE)
        
        #Check how many unique rows there are and set centers to lowest of either 10000 or #unique rows
        unique_centers<-nrow(unique(pa_habitat_data))
        center_number<-min(unique_centers, 10000)
        
        # K-means clustering
        set.seed(101)
        clust <- kmeans(pa_habitat_data[, !names(pa_habitat_data) %in% c("x", "y")], centers = center_number,iter.max = 10, nstart = 1)$cluster
        pa_habitat <- cbind(pa_habitat_data, clust)%>%
          mutate(rID =row_number())
        
        # Keep 1 pseudoabsence per cluster
        sampled <- pa_habitat %>%
          dplyr::group_by(clust) %>%
          dplyr::slice_sample(n = 1) %>%
          dplyr::ungroup()
        
        # How many pseudoabsences do we still need
        remaining <- 10000 - nrow(sampled)
        
        # Randomly sample extra pseudoabsences if fewer than 10000
        if (remaining > 0) {
          extra_pa <- pa_habitat %>%
            dplyr::filter(!rID %in% sampled$rID)%>%
            dplyr::slice_sample(n = remaining) 
          
          global_points <- bind_rows(sampled, extra_pa)
          rm(extra_pa)
          
        } else {
          global_points<- sampled
        }
        
        # Keep only three columns
        global_points <- global_points %>%
          dplyr::rename("decimalLongitude" = x,
                        "decimalLatitude" = y)%>%
          dplyr::select(decimalLongitude, decimalLatitude)%>%
          sf::st_as_sf(coords=c("decimalLongitude", "decimalLatitude"), crs=crs(biasgrid_eu), remove=FALSE)
        
        rm(pa_habitat_data, pa_habitat, sampled, remaining, unique_centers, center_number, clust)
        
      }
    }
    
    
    #--------------------------------------------
    #--- Create presence-pseudoabsence dataset ---
    #--------------------------------------------
    # Format presence data (eu_occ)
    eu_occ <- eu_occ %>%
      dplyr::mutate(species = "present") %>%
      dplyr::relocate(decimalLongitude, decimalLatitude, species, geometry)
    
    #Format pseudoabsence data (global_points) 
    global_points_sf <- global_points %>% #Keep only geometry
      dplyr::mutate(species = "absent") %>%
      dplyr::relocate(decimalLongitude, decimalLatitude, species, geometry) #Reorder columns
    
    #Combine presence and pseudoabsence data
    eu_presabs <- rbind(eu_occ, global_points_sf)
    
    
    #-----------------------------------------------------------
    #--Remove highly correlated predictors from training data --
    #-----------------------------------------------------------
    # Extract raster values at eu_presabs points
    presabs_df <- terra::extract(habitat_stack, terra::vect(eu_presabs), ID = FALSE)
    
    # Compute correlation matrix
    cor_matrix <- cor(presabs_df, use = "complete.obs")
    
    # Identify highly correlated variables
    drop_vars <- caret::findCorrelation(cor_matrix, cutoff = 0.7, exact = TRUE, names = TRUE)
    
    # Subset fullstack to keep only uncorrelated predictors
    fullstack <- subset(habitat_stack, !(names(habitat_stack) %in% drop_vars))
    
    
    #-----------------------------------------------------------
    #- Extract predictor values for presences and pseudoabsences
    #-----------------------------------------------------------
    # Extract raster values from fullstack
    occ.full.data.df <- terra::extract(fullstack, terra::vect(eu_presabs), ID = FALSE) %>%
      dplyr::mutate(occ = eu_presabs$species) 
    
    if (anyNA(occ.full.data.df)) warning("Some pseudoabsence points or occurrences fall within NA grid cells")
    
    
    #--------------------------------------------
    #------- Run models with habitat data -------
    #--------------------------------------------
    #Convert present absent to 1 0
    eu_presabs <- eu_presabs %>%
      dplyr::mutate(species = ifelse(species == "present", 1, 0))
    
    #Define SDM data and methods
    sdm_data <- sdm::sdmData(species~.,train=vect(eu_presabs), predictors= fullstack ) 
    methods <- c("glm", "gam", "bioclim", "brt", "rf", "glmpoly", "mars", "maxent", "fda","cart")
    
    #run model
    set.seed(2025)
    model <- sdm(
      species ~ ., 
      data = sdm_data,
      methods = methods  # 10 models
    )
    
    
    #--------------------------------------------
    #---  Make predictions using each model  ---
    #-------------------------------------------- 
    #Define prevalence ratio
    n1 <- nrow(eu_occ)  # presences
    n0 <- nrow(global_points_sf) # pseudoabsences 
    prev_ratio <- n1 / n0
    
    # Get model info
    info <- sdm::getModelInfo(model)
    
    #Define extents to cut fullstack into 4 latitudinal blocks to make predictions more efficient
    nblocks <- 4
    e <- terra::ext(fullstack)
    ybreaks <- seq(e$ymin, e$ymax, length.out = nblocks + 1)
    exts <- lapply(1:nblocks, function(i) ext(e$xmin, e$xmax, ybreaks[i], ybreaks[i+1]))
    pred_blocks <- vector("list", nblocks)
    
    #Create empty list to store models in
    modeloutput<-list()
    
    for(modelmethod in methods){
      
      print(modelmethod)
      
      pred_raster <- try({
        
        for(rasterblock in seq_along(exts)) {
          block_r <- crop(fullstack, exts[[rasterblock]])
          
          # Make predictions for each block
          pred_blocks[[rasterblock]] <- predict(model,
                                                newdata = block_r,
                                                method = modelmethod)
        }
        
        # Merge blocks only if all succeed
        do.call(terra::merge, pred_blocks)
        
      }, silent = TRUE)
      
      
      # If prediction failed entirely (full raster + blocks), skip to next method
      if(inherits(pred_raster, "try-error")) {
        message("Skipping method ", modelmethod, " due to prediction failure.")
        next
      } else{
        message("Predictions successfully completed for method '", modelmethod, "'.")
      }
      
      # Get model IDs
      model_ids <- info$modelID[info$method == modelmethod]
      
      # Subset using those IDs
      method_model <- model[[model_ids]]  
      
      #Apply the transformation to the raster
      fav_raster <- favourability_from_prob(pred_raster, prev_ratio)
      
      #Store
      modeloutput[[modelmethod]]<-list(fav_raster=fav_raster,
                                       model=method_model)
      
      rm(fav_raster, method_model)
    }
    
    
    #------------------------------------------------
    #----- Check if at least 8 algorithms worked ----
    #------------------------------------------------
    # if (length(modeloutput) < 8) {
    #   warning(paste0("Prediction skipped for species '", species, "': Only ",
    #                  length(modeloutput), " out of 10 algorithms successfully produced predictions.\n",
    #                  "At least 9 are required to continue: moving on to the next species."))
    #   next  # Skip to the next species 
    # }
    
    
    #---------------------------------------------
    #--------CREATE ENSEMBLE USING PCAm method----
    #---------------------------------------------
    # List favourability rasters
    fav_rasters_list <- lapply(modeloutput, function(x) x$fav_raster)
    
    # Combine into a SpatRaster stack
    fav_stack <- terra::rast(fav_rasters_list)
    
    # Assign layer names based on model methods
    names(fav_stack) <- names(modeloutput)
    
    #Clean up
    rm(modeloutput)
    gc()
    
    #make PCA: sample a representative number of pixels (in total less than 6000000 non NA cells)
    set.seed(100)
    pca_result <- rasterPCA(fav_stack, nSamples = 100000, spca = FALSE, maskCheck = TRUE)
    
    
    #-----------------GET TOP 5 variance models----------------
    # Step 1: Extract PC1 loadings from princomp object
    loadings <- pca_result$model$loadings[, 1]  # Comp.1 = PC1
    names(loadings) <- rownames(pca_result$model$loadings)
    
    # Step 2: Convert raster stack to matrix (rows = pixels, cols = models)
    fav_matrix <- as.matrix(fav_stack)
    
    # Step 3: Calculate variance along PC1 for each model
    model_variances <- setNames(numeric(nlyr(fav_stack)), names(fav_stack))
    
    for (lyr in 1:nlyr(fav_stack)) {
      model_vals <- fav_matrix[, lyr]
      centered <- model_vals - mean(model_vals, na.rm = TRUE)
      projection <- centered * loadings[lyr]
      model_variances[lyr] <- var(projection, na.rm = TRUE)
    }
    
    # Step 4: Select top 5 models with highest variance on PC1
    top5_models <- names(sort(model_variances, decreasing = TRUE))[1:5]
    cat("Top 5 models by variance along PC1:\n", paste(top5_models, collapse = ", "), "\n")
    
    # Get model IDs
    top_ids <- info$modelID[info$method %in% top5_models]
    
    # Subset using those IDs
    top5models <- model[[top_ids]]  
    
    #Clean up
    gc()
    
    #--------Create ensemble predictions using those 5 models------
    # Step 5: Subset fav_stack to top 5 layers
    top5_stack <- subset(fav_stack, top5_models)
    
    # Step 6: Compute pixel-wise median = consensus model
    consensus_habitat <- app(top5_stack, median)
    
    # Step 7: Compute pixel-wise mean
    consensus_habitat_mean <- mean(top5_stack, na.rm=TRUE)
    
    # Step 8: Compute pixel-wise population standard deviation
    consensus_habitat_sd <- stdev(top5_stack, pop=TRUE)
    
    # Step 9: Crop to extent of country if relevant
    if(tolower(country_of_interest)!="europe"||!is.null(custom_country_boundary_path)){
      ensemble_habitat_suitability<- consensus_habitat%>%
        terra::crop(country_boundary)%>%
        terra::mask(country_boundary)
      
      ensemble_habitat_sd <- consensus_habitat_sd%>%
        terra::crop(country_boundary)%>%
        terra::mask(country_boundary)
      
      ensemble_habitat_mean <- consensus_habitat_mean%>%
        terra::crop(country_boundary)%>%
        terra::mask(country_boundary)
    }else{
      ensemble_habitat_suitability<-consensus_habitat
      ensemble_habitat_sd <- consensus_habitat_sd
      ensemble_habitat_mean<- consensus_habitat_mean 
    }
    
    
    #--------------------------------------------------
    #-- Create map with ensemble habitat suitability --
    #--------------------------------------------------
    #Define name of files
    base_file <- paste0(basefile, "current_ensemble")
    
    #Export PDFs with and without occurrences plotted
    for (occs in list(NULL, eu_occ)){
      filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
      
      exportPDF(predictions = ensemble_habitat_suitability,
                dataType = "Suit",
                period = "Current",
                returnPredictions = FALSE,
                returnPNG = FALSE,
                occ_data=occs,
                exportPNG=TRUE,
                PDF_title = PDF_title,
                PNG_folder=file.path(base_dir, "Habitat", "Current", "Predictions", "PNGs"),
                PDF_folder=file.path(base_dir, "Habitat", "Current", "Predictions","PDFs"),
                filename = filename)
    }
    
    # Export ensemble raster (favorability) 
    current_habitat_folder <- file.path(base_dir, "Habitat", "Current", "Predictions", "Rasters")
    habitat_ensemble_file <- file.path(current_habitat_folder, paste0(base_file,".tif"))
    terra::writeRaster(ensemble_habitat_suitability, filename = habitat_ensemble_file, overwrite = TRUE)
    
    
    #--------------------------------------------------
    #---------- Create map with ensemble SD -----------
    #--------------------------------------------------
    #Define name of files
    filename <- paste0(basefile, "current_ensemble_SD")
    
    #Export PDFs with and without occurrences plotted
    exportPDF(predictions = ensemble_habitat_sd,
              dataType = "Stdev",
              period = "Current",
              returnPredictions = FALSE,
              returnPNG = FALSE,
              occ_data=NULL,
              exportPNG=TRUE,
              PDF_title = PDF_title,
              PNG_folder=file.path(base_dir, "Habitat", "Current", "Diagnostics","Confidence_maps", "PNGs"),
              PDF_folder=file.path(base_dir, "Habitat", "Current", "Diagnostics","Confidence_maps", "PDFs"),
              filename = filename)
    
    # Export ensemble raster (favorability) 
    current_sd_habitat_folder <- file.path(base_dir, "Habitat", "Current", "Diagnostics", "Confidence_maps", "Rasters")
    habitat_sd_ensemble_file <- file.path(current_sd_habitat_folder, paste0(filename,".tif"))
    terra::writeRaster(ensemble_habitat_sd, filename = habitat_sd_ensemble_file, overwrite = TRUE)
    
    
    #------------------------------------------
    #------------ Create binary map -----------
    #------------------------------------------
    #Create lists for storing thresholds
    habitat_thresholds<-list()
    climhab_thresholds<-list()
    
    # Get predictor values at occurrence points
    predictors_only <- occ.full.data.df%>%
      dplyr::filter(occ=="present")%>%
      dplyr::select(-occ)
    
    # Predict for top 5 models
    pred_vals <- list()
    for (method in top5_models) {
      pred_vals[[method]] <- predict(model, newdata = predictors_only, method = tolower(method))
    }
    
    # Favourability transformation
    fav_vals <- lapply(pred_vals, function(p) favourability_from_prob(p[[1]], prev_ratio))
    
    #Create one df with the median favorability value for each occurrence
    fav_vals <- fav_vals %>%
      do.call(cbind, .) %>%
      as.data.frame() %>%
      dplyr::mutate(median = apply(., 1, median, na.rm = TRUE)) %>% #1 = apply to rows
      dplyr::select(median)
    
    # Create binary maps
    for (probs in mtp_probabilities){
      
      #Define mtp_pct and mtp_value
      mtp_value<- probs*100
      mtp_pct <- paste0(mtp_value, "%")
      
      # Obtain threshold
      to_omit <- floor(probs * nrow(fav_vals)) #Define how many of lowest ranked occs to omit based on mtp threshold
      thr <- sort(fav_vals$median)[to_omit + 1]
      cat(paste0("Mean ",mtp_pct," minimum training presence threshold habitat model: ", round(thr, 4), "\n"))
      
      # Create binary raster using MTP threshold
      binary_map_pct <- ensemble_habitat_suitability >= thr  
      binary_map_pct <- as.factor( binary_map_pct*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
      levels( binary_map_pct) <- data.frame(ID = c(0, 1),
                                            class = c("Absent", "Present"))
      
      #Store raster
      binary_file <- file.path (raster_folder, paste0(basefile,"current_binary",mtp_value,"pct.tif"))
      terra::writeRaster(binary_map_pct, filename = binary_file, overwrite = TRUE)
      
      # export as PDF and PNG with and without occurrences plotted 
      base_file<- paste0(basefile, "current_binary",mtp_value,"pct")
      for (occs in list(NULL, eu_occ)){
        filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
        exportPDF(predictions = binary_map_pct,
                  dataType = "Binary",
                  period = "Current",
                  returnPredictions = FALSE,
                  returnPNG = TRUE,
                  occ_data=occs,
                  exportPNG=TRUE,
                  LabelValue= round(thr,3),
                  LabelName=paste0(mtp_pct, " MTP threshold"),
                  PDF_title = PDF_title,
                  PNG_folder=file.path(base_dir, "Habitat","Current", "Predictions", "PNGs"),
                  PDF_folder=file.path(base_dir,"Habitat" ,"Current", "Predictions", "PDFs"),
                  filename = filename)
      }
      
      habitat_thresholds[[mtp_pct]] <- thr
      rm(binary_map_pct, binary_file, thr)
    }
    
    
    #---------------------------------------------
    #-- Get response curves of 5 selected models -
    #---------------------------------------------
    response_list<-list()
    varimp_list<-list()
    
    for(topmethod in top5_models){
      # Get model id
      id <- info$modelID[info$method == topmethod]
      
      #Get response curve
      response_curves<-sdm::getResponseCurve(model,id)@response
      
      #Get variable importance
      varimp<-sdm::getVarImp(model,id)@varImportance
      
      #Store
      response_list[[topmethod]]<-response_curves
      varimp_list[[topmethod]]<-varimp
    }
    
    # Convert list to a dataframe
    response_df <- purrr::imap_dfr(response_list, function(model_list, model_name) {
      imap_dfr(model_list, function(df, var_name) {
        response_df <- df %>%
          setNames(c("Predictor_value", "Response"))%>%
          mutate( Algorithm = model_name,
                  Predictor = var_name)})}) %>%
      dplyr::select(Algorithm,Predictor, Predictor_value, Response)
    
    
    varimp_df <- imap_dfr(varimp_list, function(df, model_name) {
      df %>%
        setNames(c("Predictor", "corTest" , "AUCtest"))%>%
        dplyr::mutate(Algorithm = model_name)})%>%
      dplyr::select(Algorithm,Predictor, corTest, AUCtest)
    
    
    # Plot response curves
    response_plot <- ggplot(response_df, aes(x = Predictor_value,
                                             y = Response, 
                                             color = Algorithm)) +
      geom_line(size=0.8) +
      facet_wrap(~ Predictor, scales = "free_x")+
      labs(title= "Habitat response curves" ,x= "Predictor value")+
      theme_bw()
    
    # Plot variable importance 
    varimp_plot <- ggplot(varimp_df, aes(x = Predictor, y = corTest)) +
      geom_col(fill = "steelblue") +
      coord_flip() +  #horizontal bars
      facet_wrap(~ Algorithm) +  
      geom_hline(yintercept = 0, color = "black") + 
      labs( x = "Variable",
            y = "Importance",
            title = "Variable importance per model") +
      theme_bw()
    
    #Save plot
    PNG_folder <- file.path(base_dir, "Habitat", "Current", "Diagnostics")
    
    ggplot2::ggsave(filename = paste(basefile, "variable_importance.png"), plot = varimp_plot ,  device = "png", width =8.27 , height = 5.845, path= file.path(PNG_folder, "Variable_importance"))
    ggplot2::ggsave(filename = paste(basefile, "response_curves.png"), plot = response_plot,  device = "png", width =8.27 , height = 5.845, path=  file.path(PNG_folder, "Response_curves"))
    
    
    #------------------------------------------------------------    
    #-- Create final predictions combining habitat and climate --
    #------------------------------------------------------------
    if(tolower(country_of_interest)!="europe"||!is.null(custom_country_boundary_path)){
      consensus_climate<-terra::rast(file.path( climate_raster_folder,
                                                paste0(speciesName, "_Climate_current_ensemble_Europe.tif")))%>%
        terra::project(consensus_habitat)
    }else{
      consensus_climate<-terra::rast( file.path(climate_raster_folder,
                                                paste0(speciesName,"_Climate_current_ensemble.tif")))%>%
        terra::project(consensus_habitat)
    }
    
    #Combine suitability predictions by global model (climate) and EU habitat model
    clim_hab <- sqrt(consensus_habitat * consensus_climate)
    
    
    #--------------------------------------------------
    #--Export maps with final suitability predictions -
    #--------------------------------------------------
    # Crop to extent of country if relevant
    if(tolower(country_of_interest)!="europe"||!is.null(custom_country_boundary_path)){
      ensemble_combined_suitability<- clim_hab%>%
        terra::crop(country_boundary)%>%
        terra::mask(country_boundary)
    }else{
      ensemble_combined_suitability<-clim_hab
    }
    
    #Define name of files
    base_file <- paste0(combined_basefile, "current_ensemble")
    
    #Export raster file
    # Export continuous suitability raster
    clim_hab_file <- file.path(base_dir, "Combined", "Current", "Predictions", "Rasters",
                               paste0(base_file,".tif"))
    terra::writeRaster(ensemble_combined_suitability, filename = clim_hab_file, overwrite = T)
    
    #Export PDFs with and without occurrences plotted
    for (occs in list(NULL, eu_occ)){
      filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
      
      exportPDF(predictions = ensemble_combined_suitability,
                dataType = "Suit",
                period = "Current",
                returnPredictions = FALSE,
                returnPNG = FALSE,
                occ_data=occs,
                exportPNG=TRUE,
                PDF_title = PDF_title,
                PNG_folder=file.path(base_dir, "Combined", "Current", "Predictions", "PNGs"),
                PDF_folder=file.path(base_dir, "Combined", "Current", "Predictions","PDFs"),
                filename = filename)
    }
    
    
    #--------------------------------------------------
    #----- Create maps with final SD predictions ------
    #--------------------------------------------------
    #Load climate layers
    mean_climate_path<- file.path( base_dir,"Climate", "Current", "Interim",
                                   paste0(global_basefile, "current_ensemble_mean.tif"))
    sd_climate_path <- file.path( base_dir,"Climate", "Current","Diagnostics", "Confidence_maps", "Rasters",
                                  paste0(global_basefile, "current_ensemble_SD.tif"))
    consensus_climate_mean <- terra::rast(mean_climate_path)
    consensus_climate_sd <- terra::rast(sd_climate_path)
    
    #reproject mean climate to mean habitat crs
    consensus_climate_mean <- terra::project(consensus_climate_mean,
                                             ensemble_habitat_mean,
                                             method = "bilinear")
    consensus_climate_sd <- terra::project(consensus_climate_sd,
                                           ensemble_habitat_mean,
                                           method = "bilinear")
    
    # small floor to avoid division by zero
    eps <- 1e-6    
    
    # compute geometric mean
    S <- sqrt(consensus_climate_mean * ensemble_habitat_mean)
    
    # compute relative SDs 
    sd_climate <- consensus_climate_sd / (consensus_climate_mean + eps)
    sd_habitat <- ensemble_habitat_sd / (ensemble_habitat_mean + eps)
    
    # combined relative uncertainty 
    sd_comb <- sqrt(sd_climate^2 + sd_habitat^2)
    
    # final sd of geometric mean
    Final_SD <- 0.5 * S * sd_comb
    
    names(Final_SD) <- "sd_geometric_mean"
    
    #Define name of files
    filename <- paste0(combined_basefile, "current_ensemble_SD")
    
    #Export raster file
    clim_comb_sd_file <- file.path(base_dir, "Combined", "Current", "Diagnostics", "Confidence_maps", "Rasters",
                                   paste0(filename,".tif"))
    terra::writeRaster(Final_SD, filename = clim_comb_sd_file, overwrite = T)
    
    #Export PDFs and PNGs
    exportPDF(predictions = Final_SD,
              dataType = "Stdev",
              period = "Current",
              returnPredictions = FALSE,
              returnPNG = FALSE,
              occ_data=NULL,
              exportPNG=TRUE,
              PDF_title = PDF_title,
              PNG_folder=file.path(base_dir, "Combined", "Current", "Diagnostics", "Confidence_maps", "PNGs"),
              PDF_folder=file.path(base_dir, "Combined", "Current", "Diagnostics", "Confidence_maps", "PDFs"),
              filename = filename)
    
    
    #------------------------------------------
    #------------ Create binary map -----------
    #------------------------------------------
    # Get predicted values at occurrence points
    vals_occ <- terra::extract(clim_hab, terra::vect(eu_occ), ID=FALSE)
    
    # Create binary maps
    for (probs in mtp_probabilities){
      
      #Define mtp_pct and mtp_value
      mtp_value<- probs*100
      mtp_pct <- paste0(mtp_value, "%")
      
      # Obtain threshold
      to_omit <- floor(probs * nrow(vals_occ)) #Define how many of lowest ranked occs to omit based on mtp threshold
      thr <- sort(vals_occ[[1]])[to_omit + 1]
      cat(paste0("Mean ",mtp_pct," minimum training presence threshold: ", round(thr, 4), "\n"))
      
      # Create binary raster using MTP threshold
      binary_map_pct <- ensemble_combined_suitability >= thr 
      binary_map_pct <- as.factor( binary_map_pct*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
      levels( binary_map_pct) <- data.frame(ID = c(0, 1),
                                            class = c("Absent", "Present"))
      
      #Store raster
      raster_folder <- file.path(base_dir, "Combined","Current", "Predictions", "Rasters")
      binary_file <- file.path (raster_folder, paste0(combined_basefile,"current_binary",mtp_value,"pct.tif"))
      terra::writeRaster(binary_map_pct, filename = binary_file, overwrite = TRUE)
      
      # export as PDF and PNG with and without occurrences plotted 
      base_file<- paste0(combined_basefile, "current_binary",mtp_value,"pct")
      for (occs in list(NULL, eu_occ)){
        filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
        exportPDF(predictions = binary_map_pct,
                  dataType = "Binary",
                  period = "Current",
                  returnPredictions = FALSE,
                  returnPNG = TRUE,
                  occ_data=occs,
                  exportPNG=TRUE,
                  LabelValue= round(thr,3),
                  LabelName=paste0(mtp_pct, " MTP threshold"),
                  PDF_title = PDF_title,
                  PNG_folder=file.path(base_dir, "Combined","Current", "Predictions", "PNGs"),
                  PDF_folder=file.path(base_dir,"Combined" ,"Current", "Predictions", "PDFs"),
                  filename = filename)
      }
      
      climhab_thresholds[[mtp_pct]]<-thr
      rm(binary_map_pct, binary_file, thr)
    }
    
    
    #--------------------------------------------
    #--- Create maps with future projections ----
    #--------------------------------------------
    for (period in c("2041-2070","2071-2100")){
      for(scenario in c("ssp126", "ssp370", "ssp585")){
        
        #--------------------------------
        #--- Create suitability maps ----
        #--------------------------------
        print(paste("[FUTURE] Projecting:", period,scenario))
        
        #Get climate data for specific period and scenario
        future_folder <- file.path(base_dir, "Climate", period, scenario, "Predictions", "Rasters")
        ensemble_file <- file.path(future_folder, paste0(global_basefile, period,"_",scenario,"_ensemble.tif"))
        future_climate <- terra::rast(ensemble_file)%>%
          terra::project(ensemble_habitat_suitability)
        
        #Final ensemble predictions
        final_ensemble<-sqrt(ensemble_habitat_suitability * future_climate)
        
        # Export future ensemble raster (favorability) 
        future_folder <- file.path(base_dir, "Combined", period, scenario, "Predictions", "Rasters")
        ensemble_file <- file.path(future_folder, paste0(combined_basefile, period,"_",scenario,"_ensemble.tif"))
        terra::writeRaster(final_ensemble, filename = ensemble_file, overwrite = TRUE)
        
        # Export ensemble predictions as PDF and PNG with and without occurrences
        base_file <- paste0(combined_basefile, scenario,"_", period,"_ensemble")
        
        for (occs in list(NULL, eu_occ)){
          filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
          
          exportPDF(predictions = final_ensemble,
                    dataType = "Suit",
                    period = period,
                    scenario = scenario,
                    returnPredictions = FALSE,
                    returnPNG = TRUE,
                    occ_data=occs,
                    exportPNG=TRUE,
                    PDF_title=PDF_title,
                    PNG_folder=file.path(base_dir, "Combined", period, scenario, "Predictions", "PNGs"),
                    PDF_folder=file.path(base_dir, "Combined", period, scenario, "Predictions", "PDFs"),
                    filename = filename)
        }
        
        
        # Create binarized ensemble predictions for future
        for(probs in mtp_probabilities){
          
          #Define mtp_pct and mtp_value
          mtp_value<- probs*100
          mtp_pct <- paste0(mtp_value, "%")
          mtp_thr <- paste0(mtp_value, "pct")
          
          #Get threshold value and apply to consensus predictions
          threshold<- climhab_thresholds[[mtp_pct]]
          binary_map_future <- final_ensemble  >= threshold
          binary_map_future <- as.factor( binary_map_future*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
          levels( binary_map_future) <- data.frame(ID = c(0, 1),
                                                   class = c("Absent", "Present"))
          
          #Store raster
          future_combined_folder <- file.path(base_dir, "Combined", period, scenario, "Predictions", "Rasters")
          binary_file <- file.path(future_combined_folder, paste0(combined_basefile, period,"_",scenario,"_binary",mtp_thr,".tif"))
          terra::writeRaster(binary_map_future, filename = binary_file, overwrite = TRUE)
          
          # Export binarized ensemble predictions as PDF and PNG with and without occurrences 
          base_file <- paste0(combined_basefile, period,"_", scenario, "_binary",mtp_thr)
          
          for (occs in list(NULL, eu_occ)){
            
            filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
            exportPDF(predictions = binary_map_future,
                      dataType = "Binary",
                      period = period,
                      scenario = scenario,
                      occ_data = occs,
                      returnPredictions = FALSE,
                      returnPNG = FALSE,
                      exportPNG = TRUE,
                      LabelValue= round(threshold,2),
                      LabelName=paste0(mtp_pct, " MTP threshold"),
                      PDF_title=PDF_title,
                      PNG_folder=file.path(base_dir, "Combined", period, scenario, "Predictions", "PNGs"),
                      PDF_folder=file.path(base_dir, "Combined", period, scenario, "Predictions", "PDFs"),
                      filename=filename)
          }
        }
        
        #--------------------------------
        #---- Create confidence maps ----
        #--------------------------------
        #Define file paths for future climate SD and mean files
        future_sd_folder <- file.path(base_dir, "Climate", period, scenario, "Diagnostics", "Confidence_maps", "Rasters")
        sd_future_climate_path <-  file.path(future_sd_folder, paste0(global_basefile, period,"_",scenario,"_ensemble_SD.tif"))
        future_mean_folder <- file.path(base_dir, "Climate", "Current", "Interim")
        mean_future_climate_path <- file.path(future_mean_folder, paste0(global_basefile, period,"_",scenario,"_ensemble_mean.tif"))
        
        #Load future climate SD and mean files
        consensus_future_climate_mean <- terra::rast(mean_future_climate_path)
        consensus_future_climate_sd <- terra::rast(sd_future_climate_path)
        
        #reproject mean future_climate to mean habitat crs
        consensus_future_climate_mean <- terra::project(consensus_future_climate_mean,
                                                        ensemble_habitat_mean,
                                                        method = "bilinear")
        consensus_future_climate_sd <- terra::project(consensus_future_climate_sd,
                                                      ensemble_habitat_mean,
                                                      method = "bilinear")
        
        # small floor to avoid division by zero; adjust if needed
        eps <- 1e-6    
        
        # compute geometric mean
        S <- sqrt(consensus_future_climate_mean * ensemble_habitat_mean)
        
        # compute relative SDs safely
        sd_future_climate <- consensus_future_climate_sd / (consensus_future_climate_mean + eps)
        sd_habitat <- ensemble_habitat_sd / (ensemble_habitat_mean + eps)
        
        # combined relative uncertainty (root-sum-of-squares)
        sd_comb <- sqrt(sd_future_climate^2 + sd_habitat^2)
        
        # final sd of geometric mean
        Final_future_SD <- 0.5 * S * sd_comb
        
        names(Final_future_SD) <- "sd_geometric_mean"
        
        #Define name of files
        filename <- paste0(combined_basefile, period, "_",scenario,"_ensemble_SD")
        
        #Export raster file
        future_sd_file <- file.path(base_dir, "Combined", period, scenario, "Diagnostics", "Confidence_maps", "Rasters",
                                    paste0(filename,".tif"))
        terra::writeRaster(Final_future_SD, filename = future_sd_file, overwrite = T)
        
        #Export PDFs and PNGs
        exportPDF(predictions = Final_future_SD,
                  dataType = "Stdev",
                  period = period,
                  scenario = scenario,
                  returnPredictions = FALSE,
                  returnPNG = FALSE,
                  occ_data=NULL,
                  exportPNG=TRUE,
                  PDF_title = PDF_title,
                  PNG_folder=file.path(base_dir, "Combined", period, scenario, "Diagnostics", "Confidence_maps", "PNGs"),
                  PDF_folder=file.path(base_dir, "Combined", period, scenario, "Diagnostics", "Confidence_maps", "PDFs"),
                  filename = filename)
        
        rm(S, consensus_future_climate_mean, consensus_future_climate_sd, sd_future_climate, sd_comb, Final_future_SD)
      }
    }
    
    
    #--------------------------------------------
    #- Save best model, european occurrences, and layers for Belgium -
    #--------------------------------------------
    habitatmodel <- list(species = species,
                         taxonkey = taxonkey,
                         eu_occ = eu_occ,  # sf of filtered EU occurrences
                         eu_presabs = eu_presabs,    # sf of presence + pseudoabsence data
                         occ_full_df = occ.full.data.df, # presabs data and their habitat values
                         prevalence_ratio = prev_ratio, # used for favourability scaling
                         habitat_thresholds = habitat_thresholds,
                         climhab_thresholds = climhab_thresholds,
                         response_df = response_df,
                         varimp_df = varimp_df,
                         selected_predictors = names(fullstack),
                         top5models = top5models, #model object holding selected models
                         top5_models = top5_models
    )
    
    #Save eumodel as .qs file
    qs2::qs_save(habitatmodel, 
              file.path(base_dir,"Habitat", paste0("Habitat_model_",speciesName,"_",taxonkey,".qs"))) 
    
    
    #--------------------------------------------
    #-------- End of loop -----------------------
    #--------------------------------------------
    end_time <- Sys.time()
    elapsed<-difftime(end_time, start_time, units="mins")
    cat("Habitat and ensemble model have been created for", species_title, "in", round(elapsed, 2), "minutes\n\n")
    
    rm(list = setdiff(ls(), c("p", "project",  "habitatstack_file","create_folder", "custom_country_boundary_path","country_boundary", "split_df","euboundary", "habitat_stack",  "accepted_taxonkeys", "taxa_info", "key", "exportPDF", "remove_duplicates", "remove_nodata_occurrences", "favourability_from_prob", "mtp_probabilities", "occurrence_thinning_method", "mtp_probabilities", "pseudoabsence_thinning_method", "country_of_interest")))
    
    #Clean terra tempfiles
    terra::tmpFiles(remove = TRUE)
    
  }
})


#--------------------------------------------
#----------- Clean R environment ------------
#--------------------------------------------
rm(list = ls())

