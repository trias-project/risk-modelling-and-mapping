#--------------------------------------------
#-----------    Load packages      ----------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")
setGDALconfig("GDAL_PAM_ENABLED", "FALSE")#Prevent terra from writing aux.xml files

packages <- c( "viridis","dplyr", "here", "qs","terra", "tidyterra","sf", "ggplot2","RColorBrewer","magick","patchwork","grid", "randomForest", "progressr", "raster", "dismo", "caret", "caretEnsemble", "kableExtra","gbm", "PresenceAbsence", "RStoolbox", "sdm")

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
#---------   Load shape of Europe   ---------
#--------------------------------------------
euboundary<-sf::st_read(here("./data/external/GIS/Europe/EUROPE.shp")) 


#--------------------------------------------
#-------- Load European habitat rasters -----
#--------------------------------------------
# Load all habitat rasters
habitat_files <- list.files(file.path("./data/external/habitat"), pattern = 'tif$', full.names = TRUE)
habitat_rasters <- lapply(habitat_files, terra::rast)

# compute common intersection extent across all rasters
common_ext <- Reduce(intersect, lapply(habitat_rasters, ext))

# Crop all rasters to the common (smallest) extent
habitat_rasters <- lapply(habitat_rasters, terra::crop, common_ext)

# Combine into raster stack 
habitat_stack <- terra::rast(habitat_rasters)
rm(habitat_rasters)

#Scale habitat rasters
habitat_stack <- terra::scale(habitat_stack, center = TRUE, scale = TRUE)


#---------------------------------------------
#----- Remove NA pixels from predictors ------
#---------------------------------------------
#This is to avoid that some layers have NA while others have values in certain pixels
na_mask_habitat_stack <- anyNA(habitat_stack)
habitat_stack <- terra::mask(habitat_stack, na_mask_habitat_stack, maskvalue=1)


#---------------------------------------------
#--------- Load WWF ecoregions file ----------
#---------------------------------------------
wwf_eco_biome <- sf::st_read(here("./data/external/GIS/official/wwf_terr_ecos.shp")) %>%
  sf::st_transform(crs = st_crs(habitat_stack))


#--------------------------------------------
#------------- Load species data -----------
#--------------------------------------------
taxa_info <- read.csv2(paste0("./data/projects/",project,"/",project,"_taxa_info.csv"))
accepted_taxonkeys <- taxa_info %>%
  dplyr::pull(acceptedTaxonKey) %>%
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
    basefile<-  paste0(speciesName, "_Europe_")
    global_basefile<-  paste0(speciesName, "_", taxonkey, "_Global_")
    
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
      globalmodels <- qs::qread(global_model_file)
      
      #Extract different data objects stored in globalmodels
      global.occ.sf <- globalmodels$occurrences1km %>% # FULL occurrence per km²
        sf::st_as_sf(.,coords = c("decimalLongitude", "decimalLatitude"), crs=4326)
      
    }else{
      warning(paste0("Skipping species ", species, " because no global model could be fitted"))
      next  # Skip the rest of the loop and move to the next iteration
    }
    
    
    #--------------------------------------------
    #------------ Import raster layers ----------
    #--------------------------------------------
    #Define file paths
    biasgrid_file <- file.path(base_dir,"Climate", "Current", "Interim", paste0("Biasgrid_",speciesName,"_",taxonkey,".tif"))
    global_model_file <- file.path(base_dir,"Climate", "Current","Predictions","Rasters",
                                   paste0(speciesName,"_Climate_current_ensemble.tif"))
    
    #Load rasterlayers
    biasgrid_sub <- terra::rast(biasgrid_file)
    global_climate_for_eu <- terra::rast(global_model_file)%>%
      terra::project( habitat_stack)
        
    
    #-------------------------------------------------
    #--------------- Create EU folders ---------------
    #-------------------------------------------------
    # Define outputs, periods, and scenarios
    outputs   <- c("Rasters", "PDFs", "PNGs")
    periods   <- c("Current","2041-2070", "2071-2100")
    scenarios <- c("ssp126", "ssp370", "ssp585")
    
    #Create folders for each combination
    scenario_folders <- list()
    
    for(output in outputs){
      for(period in periods){
        if(period=="Current"){
          loop_list <- list(list(path = file.path(base_dir, output, "Europe", period),
                                 name = paste(output, "Europe", period,  sep = "/")))
          scenario_folders <- c(scenario_folders, loop_list)  
          
        }else{
          for(scenario in scenarios){
            loop_list <- list(list(path = file.path(base_dir, output, "Europe", period, scenario),
                                   name = paste(output, "Europe", period, scenario, sep = "/")))
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
    eu_occ <- sf::st_filter(global.occ.sf, euboundary) %>%
      st_transform(crs = st_crs(habitat_stack)) %>%
      sf::st_coordinates() %>%
      as.data.frame()
 
    
    #-----------------------------------------------
    #----------- Process occurrences ---------------
    #-----------------------------------------------
    #Keep only one occurrence per grid cell
    eu_occ <- remove_duplicates(occurrences =  eu_occ, rast_template = habitat_stack[[1]])
    
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
        
        # K-means clustering
        set.seed(101)
        clust <- kmeans(habitat_data, centers = 10000,iter.max = 10, nstart = 1)$cluster
        occ_habitat <- cbind(eu_occ, habitat_data, clust)
        
        # Keep 1 occurrence per cluster
        eu_occ <- occ_habitat %>%
          dplyr::group_by(clust) %>%
          dplyr::slice_sample(n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::select(decimalLongitude, decimalLatitude, geometry, species)
        
        rm(habitat_data, occ_habitat)
        
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
    # geom_sf(data = sf::st_transform(euboundary, crs=st_crs(habitat_stack)),  colour = "black", fill = NA)+
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
    # Identify which polygons contain at least one occurrence
    polygons_with_points <- lengths(sf::st_intersects(wwf_eco_biome, eu_occ)) > 0
    
    # Subset only those polygons
    wwf_eco_biome_filtered <- wwf_eco_biome[polygons_with_points, ]
   # plot(wwf_eco_biome_filtered[4], key.pos = NULL)
    
    
    #----------------------------------------------------------------------------------------
    #---- biasgrid: keep values inside invaded ecoregions, set outside to 1 (lowest value)---
    #----------------------------------------------------------------------------------------
    # Step 1: Rasterize WWF polygons to match biasgrid_aligned
    inside_mask <- terra::rasterize(vect(wwf_eco_biome_filtered), biasgrid_aligned, field = 1, background = NA)
    
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
    #Generate pseudoabsences in according to sampling bias---------
    #--------------------------------------------------------------
    set.seed(728)
    global_points <- terra::spatSample(
      biasgrid_eu,
      size = 10000,
      method = "weights",     # weighted random sampling
      as.points = TRUE,       # return SpatVector of points
      na.rm = TRUE            # ignore NA pixels
    )
    
    
    #--------------------------------------------
    #--- Create presence-pseudoabsence dataset ---
    #--------------------------------------------
    # Format presence data (eu_occ)
    eu_occ <- eu_occ %>%
      dplyr::mutate(species = "present") %>%
      dplyr::relocate(decimalLongitude, decimalLatitude, species, geometry)
  
    #Format pseudoabsence data (global_points) 
    global_points_sf <- global_points[, 0] %>% #Keep only geometry
      sf::st_as_sf() %>% #Convert to sf
      dplyr::mutate(coords = st_coordinates(geometry),#Extract coordinates as matrix
             decimalLongitude = coords[,1], #Add coordinates to designated column
             decimalLatitude  = coords[,2], #Add coordinates to designated column
             species = "absent") %>%
      dplyr::select(-coords) %>%  # drop helper column
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
    sdm_data <- sdm::sdmData(species~.,train=vect(eu_presabs),predictors= fullstack ) 
    methods <- c("glm", "gam", "bioclim", "brt", "rf", "glmpoly", "mars", "maxent", "fda","cart")
    
    #run model
    model <- sdm(
      species ~ ., data = sdm_data,
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
    
    #Create empty list to store models in
    modeloutput<-list()
    
    #Around 22 min for one species
      for(modelmethod in methods){
        
        print(modelmethod)
        
        #Create raster with predictions for Europe
        pred_raster <- predict(model,
                               newdata = fullstack,
                               method = modelmethod)
        
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
    
    
    #---------------------------------------------
    #--------CREATE ENSEMBLE USING PCAm method----
    #---------------------------------------------
    # List favourability rasters
    fav_rasters_list <- lapply(modeloutput, function(x) x$fav_raster)
    
    # Combine into a SpatRaster stack
    fav_stack <- terra::rast(fav_rasters_list)
    
    # Assign layer names based on model methods
    names(fav_stack) <- names(modeloutput)
    
    #make PCA
    pca_result <- rasterPCA(fav_stack, nSamples = NULL, spca = FALSE, maskCheck = TRUE)
    
    
    #-----------------GET TOP 5 variance models----------------
    # Step 1: Recover original raster stack used in rasterPCA
    fav_stack <- eval(pca_result$call$img)
    
    # Step 2: Extract PC1 loadings from princomp object
    loadings <- pca_result$model$loadings[, 1]  # Comp.1 = PC1
    names(loadings) <- rownames(pca_result$model$loadings)
    
    # Step 3: Convert raster stack to matrix (rows = pixels, cols = models)
    fav_matrix <- as.matrix(fav_stack)
    
    # Step 4: Calculate variance along PC1 for each model
    model_variances <- setNames(numeric(nlyr(fav_stack)), names(fav_stack))
    
    for (lyr in 1:nlyr(fav_stack)) {
      model_vals <- fav_matrix[, lyr]
      centered <- model_vals - mean(model_vals, na.rm = TRUE)
      projection <- centered * loadings[lyr]
      model_variances[lyr] <- var(projection, na.rm = TRUE)
    }
    
    # Step 5: Select top 5 models with highest variance on PC1
    top5_models <- names(sort(model_variances, decreasing = TRUE))[1:5]
    cat(paste("Top 5 models by variance along PC1:\n", top5_models))
    
    # Get model IDs
    top_ids <- info$modelID[info$method %in% top5_models]
    
    # Subset using those IDs
    top5models <- model[[top_ids]]  
    
    # Step 6: Subset fav_stack to top 5 layers
    top5_stack <- subset(fav_stack, top5_models)
    
    # Step 7: Compute pixel-wise median = consensus model
    consensus_habitat <- app(top5_stack, median)

   
    #------------------------------------------------------------    
    #-- Create final predictions combining habitat and climate --
    #------------------------------------------------------------
    #Combine suitability predictions by global model (climate) and EU habitat model
    clim_hab <- sqrt(consensus_habitat * global_climate_for_eu)
    
    #Extract suitability predictions of EU occurrences
    vals_occ <- terra::extract(clim_hab, terra::vect(eu_occ), ID=FALSE)
    
    #Define 5% and 1% reclassification thresholds
    threshold_5pct <- quantile(vals_occ, probs = 0.05, na.rm = TRUE)
    threshold_1pct <- quantile(vals_occ, probs = 0.01, na.rm = TRUE)
    
    # Create binary map: 1 if ≥ threshold, 0 otherwise
    clim_hab_binary_1pct <- clim_hab >= threshold_1pct
    clim_hab_binary_5pct <- clim_hab >= threshold_5pct
    
    # Convert TRUE FALSE to present Absent
    clim_hab_binary_1pct <- as.factor( clim_hab_binary_1pct*1) 
    levels(clim_hab_binary_1pct) <- data.frame(ID = c(0, 1),
                                          class = c("Absent", "Present"))

    clim_hab_binary_5pct <- as.factor( clim_hab_binary_5pct*1) 
    levels(clim_hab_binary_5pct) <- data.frame(ID = c(0, 1),
                                               class = c("Absent", "Present"))
    
    
    #------------------------------------------------------------    
    #---------- Calculate sensitivity and Boyce Index -----------
    #------------------------------------------------------------
    
    #SENSITIVITY
    vals_1pct <- terra::extract(clim_hab_binary_1pct, eu_occ, ID = FALSE)
    FN_1pct <- sum(vals_1pct$class == "Absent")#Counts false negatives
    TP_1pct <- sum(vals_1pct$class == "Present") #Counts true positives
    final_sensitivity_1pct <- TP_1pct / (TP_1pct + FN_1pct)
    
    vals_5pct <- terra::extract(clim_hab_binary_5pct, eu_occ, ID = FALSE)
    FN_5pct <- sum(vals_5pct$class == "Absent") #Counts false negatives
    TP_5pct  <- sum(vals_5pct$class == "Present") #Counts true positives
    final_sensitivity_5pct <- TP_5pct / (TP_5pct + FN_5pct)
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
    response_df <- imap_dfr(response_list, function(model_list, model_name) {
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
    #Combine suitability predictions by global model (climate) and EU habitat model
    clim_hab <- sqrt(consensus_habitat * global_climate_for_eu)
    
    
    #------------------------------------------
    #-Calculate Boyce index to display on map--
    #------------------------------------------
    #Extract all raster values (excluding NAs)
    pred_vals <- values(clim_hab)
    pred_vals <- pred_vals[!is.na(pred_vals)]
    
    #Extract suitability values at occurrence locations
    obs_vals <- terra::extract(clim_hab, eu_occ, ID = FALSE) 
    obs_vals <- obs_vals[!is.na(obs_vals)]
    
    #Calculate Boyce index
    boyce_result <- ecospat::ecospat.boyce(
      fit = pred_vals,
      obs = obs_vals,
      nclass = 0  # continuous Boyce Index
    )
    
    if (exists("boyce_result") && !is.null(boyce_result$cor)) {
      boyce_val <- round(boyce_result$cor, 3)
    } else {
      boyce_val <- NA
    }
    
    
    #------------------------------------------------------------    
    #---------- Create map with ensemble suitability ------------
    #------------------------------------------------------------
    #Define name of files
    base_file <- paste0(basefile, "_hist_ensemble")
    
    #Export PDFs with and without occurrences plotted
    for (occs in list(NULL, eu_occ)){
      filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
      exportPDF(predictions = clim_hab,
                dataType = "Suit",
                scenario = "hist",
                returnPredictions = FALSE,
                returnPNG = FALSE,
                occ_data=occs,
                exportPNG=TRUE,
                LabelValue= boyce_val,
                LabelName="Boyce index",
                PDF_title = PDF_title,
                PNG_folder=file.path(base_dir, "Combined", "Current", "Predictions", "PNGs"),
                PDF_folder=file.path(base_dir, "Combined", "Current", "Predictions","PDFs"),
                filename = filename)
    }
    
    
    #------------------------------------------
    #------------ Create binary map -----------
    #------------------------------------------
    binary_maps <- list( "0.05" = list( binary_raster = clim_hab_binary_5pct,
                                        EU_sensitivity = final_sensitivity_5pct,
                                        boyce = boyce_val,
                                        mtp_pct = "5%",
                                        mtp_value = 5),
                         "0.01" = list(binary_raster = clim_hab_binary_1pct,
                                       EU_sensitivity = final_sensitivity_1pct,
                                       boyce = boyce_val,
                                       mtp_pct = "1%",
                                       mtp_value = 1))
    
    for (pct in seq_along(binary_maps)){
      
      binary_map_pct <- binary_maps[[pct]]$binary_raster
      EU_sensitivity <- binary_maps[[pct]]$EU_sensitivity
      mtp_pct <- binary_maps[[pct]]$mtp_pct
      mtp_value <- binary_maps[[pct]]$mtp_value
      #boyce_ind <- binary_maps[[pct]]$boyce
      
      # export as PDF and PNG with and without occurrences plotted and return as PNG
      base_file<- paste0(basefile, "_hist_ensemble_binary",mtp_value,"pct")
      for (occs in list(NULL, eu_occ)){
        filename <- ifelse(is.null(occs), base_file, paste0(base_file, "_occ"))
        exportPDF(predictions = binary_map_pct,
                  dataType = "Binary",
                  scenario = "hist",
                  returnPredictions = FALSE,
                  returnPNG = FALSE,
                  occ_data=occs,
                  exportPNG=TRUE,
                  LabelValue= mtp_pct,
                  LabelName="MTP threshold",
                  Label2Value=round(EU_sensitivity,3),
                  Label2Name="Sensitivity",
                  PDF_title = PDF_title,
                  PNG_folder=here::here(base_dir, "PNGs", "Europe","Current"),
                  PDF_folder=here::here(base_dir, "PDFs","Europe" ,"Current"),
                  filename = filename)
      }
    }

    
    #--------------------------------------------
    #--- Create maps with future projections ----
    #--------------------------------------------
    for (period in c("2041-2070","2071-2100")){
      for(scenario in c("ssp126", "ssp370", "ssp585")){
        
        print(paste("[FUTURE] Projecting:", period,scenario))
        
        #Get climate data for specific period and scenario
        future_folder <- file.path(base_dir, "Rasters", "Global", period, scenario)
        ensemble_file <- file.path(future_folder, paste0(global_basefile, period,"_",scenario,"_ensemble.tif"))
        future_climate <- terra::rast(ensemble_file)%>%
          terra::project(consensus_habitat)
        
        #Final ensemble predictions
        final_ensemble<-sqrt(consensus_habitat * future_climate)
        
        # Export ensemble predictions as PDF and PNG with and without occurrences
        base_file <- paste0(basefile, scenario,"_", period,"_final_ensemble")
        
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
                    PNG_folder=file.path(base_dir, "PNGs", "Europe", period, scenario),
                    PDF_folder=file.path(base_dir, "PDFs", "Europe", period, scenario),
                    filename = filename)
        }
        
        
        # Create binarized ensemble predictions for future
        for(MTP_threshold in c("threshold_1pct","threshold_5pct")){
          
          mtp_label <- switch(MTP_threshold,
                              "threshold_5pct" = "5%",
                              "threshold_1pct" = "1%")  
          mtp_thr <- switch(MTP_threshold,
                              "threshold_5pct" = "5pct",
                              "threshold_1pct" = "1pct")  
          
          #Get threshold value and apply to consensus predictions
          threshold<-get(MTP_threshold)
          binary_map_future <- final_ensemble  >= threshold
          binary_map_future <- as.factor( binary_map_future*1) #Convert TRUE/FALSE to 1/0 and then to Present/Absent
          levels( binary_map_future) <- data.frame(ID = c(0, 1),
                                                   class = c("Absent", "Present"))
          
          #Store raster
          future_europe_folder <- file.path(base_dir, "Rasters", "Europe", period, scenario)
          binary_file <- file.path(future_europe_folder, paste0(basefile, period,"_",scenario,"_final_binary",mtp_thr,".tif"))
          terra::writeRaster(binary_map_future, filename = binary_file, overwrite = TRUE)
          
          # Export binarized ensemble predictions as PDF and PNG with and without occurrences 
          base_file <- paste0(basefile, period,"_", scenario, "_final_binary",mtp_thr)
          
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
                      LabelValue= mtp_label,
                      LabelName="MTP threshold",
                      PDF_title=PDF_title,
                      PNG_folder=file.path(base_dir, "PNGs","Europe", period, scenario),
                      PDF_folder=file.path(base_dir, "PDFs", "Europe",period, scenario),
                      filename=filename)
          }
        }
      }
    }
    
  
    #--------------------------------------------
    #- Save best model, european occurrences, and layers for Belgium -
    #--------------------------------------------
    eumodel <- list(species = species,
                    taxonkey = taxonkey,
                    eu_occ = eu_occ,  # sf of filtered EU occurrences
                    eu_presabs = eu_presabs,    # sf of presence + pseudoabsence data
                    occ_full_df = occ.full.data.df, # presabs data and their habitat values
                    prevalence_ratio = prev_ratio,  # used for favourability scaling
                    eu_5pct_threshold = threshold_5pct,# 5% min training presence threshold from EU model
                    eu_1pct_threshold = threshold_1pct, # 1% min training presence threshold from EU model
                    final_boyce = boyce_val,
                    boyce_result = boyce_result,# full object from ecospat.boyce()
                    response_df = response_df,
                    varimp_df = varimp_df,
                    top5models = top5models #model object holding selected models
    )
   
    #Save eumodel as .qs file
    qs::qsave(eumodel, 
              file.path(base_dir, paste0("EU_model_",speciesName,"_",taxonkey,".qs"))) 
    
    
    #--------------------------------------------
    #- ------ Save rasters-----------------------
    #--------------------------------------------
    # Export continuous suitability raster
    clim_hab_file <- file.path(base_dir, "Rasters","Europe", "Current",
                               paste0(paste0(basefile, "final_hist_ensemble.tif")))
    
    terra::writeRaster(clim_hab, filename = clim_hab_file, overwrite = T)
    
    # Export binary raster — 1% threshold
    clim_hab_binary_1pct_file <- file.path(raster_EU_folder,
                                           paste0("Climate_Habitat_binary_1pct_",
                                                  speciesName,
                                                  "_", 
                                                  taxonkey,
                                                  ".tif"))
    
    terra::writeRaster(clim_hab_binary_1pct, filename = clim_hab_binary_1pct_file, overwrite = T)
    
    # Export binary raster — 5% threshold
    clim_hab_binary_5pct_file <- file.path(raster_EU_folder, 
                                           paste0("Climate_Habitat_binary_5pct_",
                                                  speciesName,
                                                  "_",
                                                  taxonkey,
                                                  ".tif"))
    
    terra::writeRaster(clim_hab_binary_5pct, filename = clim_hab_binary_5pct_file, overwrite = T)
    
    
   
    #--------------------------------------------
    #-------- End of loop -----------------------
    #--------------------------------------------
    print(paste("European model has been created for", species_title))
    rm(list = setdiff(ls(), c("p", "project",  "create_folder",  "euboundary", "habitat_stack",  "accepted_taxonkeys", "taxa_info", "key", "exportPDF", "remove_duplicates", "wwf_eco_biome", "remove_nodata_occurrences", "favourability_from_prob", "mtp_probabilities")))
  }
})


#--------------------------------------------
#----------- Clean R environment ------------
#--------------------------------------------
rm(list = ls())

