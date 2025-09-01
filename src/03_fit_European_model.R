#--------------------------------------------
#-----------To do: specify project ----------
#--------------------------------------------
#specify project name
projectname <-"PA prob & Alternative Treshold & Ensemble Boyce"


#--------------------------------------------
#-----------    Load packages      ----------
#--------------------------------------------
options("rgdal_show_exportToProj4_warnings"="none")

packages <- c( "viridis","dplyr", "here", "qs","terra", "tidyterra","sf", "ggplot2","RColorBrewer","magick","patchwork","grid", "randomForest", "progressr", "raster", "dismo", "caret", "caretEnsemble", "kableExtra","gbm", "PresenceAbsence", "RStoolbox", "sdm")

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}

sdm::installAll()


#--------------------------------------------
#------- Source helper fucntions     --------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#---------   Load shape of Europe   ---------
#--------------------------------------------
euboundary<-sf::st_read(here("./data/external/GIS/Europe/EUROPE.shp")) 


#--------------------------------------------
#-------- Load European habitat rasters -----
#--------------------------------------------
# Load all habitat rasters
habitat_files <- list.files(here("./data/external/habitat"), pattern = 'tif$', full.names = TRUE)
habitat_stack <- terra::rast(habitat_files)

# Assign meaningful and unique names (based on file order)
names(habitat_stack) <- c(
  "corine_perAgriculture",
  "corine_perConiferous",
  "corine_perdeciduous",
  "corine_pergrass",
  "dist_to_water_log1p",         # <- from distance_to_water_masked_log1p.tif
  "esm_index",
  "water_line_log1p",            # <- from total_water_length_log1p.tif
  "water_polygon_cover"          # <- from total_water_polygon_cover_proportion.tif
)

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
taxa_info <- read.csv2(paste0("./data/projects/",projectname,"/",projectname,"_taxa_info.csv"))
accepted_taxonkeys <- taxa_info %>%
  dplyr::pull(speciesKey) %>%
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
    species<-taxa_info%>%
      dplyr::filter(speciesKey==key)%>%
      dplyr::pull(acceptedScientificName)%>%
      unique()
    
    #Extract first two words of species name
    first_two_words <- sub("^(\\w+)\\s+(\\w+).*", "\\1_\\2", species)
    
    #Extract rest of species name
    rest_of_name <- if (grepl("^\\S+\\s+\\S+$", species)) "" else sub("^\\S+\\s+\\S+\\s+", "", species)
    
    #Specify species for plot title
    species_title <- gsub("_", " ", first_two_words)
    
    #Define taxonkey
    taxonkey<- key
    
    
    #--------------------------------------------
    #-- Define file path of global model file  --
    #--------------------------------------------  
    species_folder <- here::here("data", "projects", projectname, paste0(first_two_words,"_",taxonkey))
    global_model_file <- here::here(species_folder,
                                  paste0("Global_model_",first_two_words,"_",taxonkey,".qs"))
    
    
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
    biasgrid_file <- here::here("data", "projects",projectname, paste0(first_two_words,"_",taxonkey),"Rasters","Interim",paste0("Biasgrid_",first_two_words,"_",taxonkey,".tif"))
    global_model_file <- here::here("data","projects",projectname,paste0(first_two_words,"_",taxonkey),"Rasters","Global",paste0("Ensemble_median_",first_two_words,"_",taxonkey,".tif"))
    
    #Load rasterlayers
    biasgrid_sub <- terra::rast(biasgrid_file)
    global_climate_for_eu <- terra::rast(global_model_file)%>%
      terra::project( habitat_stack)
        
    #--------------------------------------------
    #------------ Define folder paths -----------
    #--------------------------------------------
    PDF_folder <- here::here("data", "projects", projectname, paste0(first_two_words, "_", taxonkey), "PDFs")
    PNG_folder <- here::here("data","projects", projectname, paste0(first_two_words, "_", taxonkey), "PNGs")
    raster_EU_folder <- here::here("data", "projects", projectname, paste0(first_two_words, "_", taxonkey), "Rasters", "Europe")


    
    #-------------------------------------------------
    #--------------- Create EU folders ---------------
    #-------------------------------------------------
    # Define the folder paths
    folder_paths<-list(list("path"= raster_EU_folder,
                            "name"= "Rasters/Europe"),
                       list("path"= here::here(PDF_folder, "Europe"),
                            "name"= here::here(PNG_folder, "Europe"),
                            "name"= "PNG/Europe"))
    
    # Check and create each folder if necessary
    lapply(folder_paths, function(folder){
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
    eu_occ <- remove_nodata_occurrences(occurrences = eu_occ, rast_template= habitat_stack, st_crs(habitat_stack))
    
    # Keep XY coordinates
    euocc<-eu_occ%>%
      sf::st_coordinates()
    
    
    #------------------------------------------------
    #----- Check if at least 20 European records ----
    #------------------------------------------------
    if (nrow(euocc) < 20) {
      warning(paste(nrow(euocc)," occurrences in Europe for species:", species, 
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
    # Ensure CRS match
    if (!st_crs(eu_occ) == st_crs(euboundary)) {
      euboundary <- st_transform(euboundary, st_crs(eu_occ))
      print("transforming CRS of euboundary")
    }
    
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
    
    if (anyNA(occ.full.data.df)) warning("Some pseudoabsence points fall within NA grid cells")
    
    occ_counts <- table(occ.full.data.df$occ)
    
    
    #--------------------------------------------
    #-- Run models with climate and habitat data -
    #--------------------------------------------
    #Define prevalence ratio
    n1 <- occ_counts["present"]   # presences
    n0 <- occ_counts["absent"]    # pseudoabsences 
    prev_ratio <- n1 / n0
    
    #define methods and data
    eu_presabs <- eu_presabs %>%
      dplyr::mutate(species = ifelse(species == "present", 1, 0))
    
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
    
    # Get model info
    info <- sdm::getModelInfo(model)
    
    #Get presence data and their associated climatological values
    pres_features <- occ.full.data.df %>%
      dplyr::filter(occ == "present") %>%
      dplyr::select(-occ)

    #Create empty list to store models in
    modeloutput<-list()
    
    #Around 22min for one species
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
        
        #Get threshold
        fitted_model <- predict(method_model, newdata = pres_features, type = "response")
        fitted_favourability <- favourability_from_prob(fitted_model[[1]], prev_ratio)
        
        # 1% minimum training presence threshold
        threshold_1pct <- quantile(fitted_favourability, probs = 0.01, na.rm = TRUE)
        
        #5% minimum training presence threshold
        threshold_5pct <- quantile(fitted_favourability, probs = 0.05, na.rm = TRUE)
        
        # Binarize rasters using the thresholds
        binary_1pct <- fav_raster >= threshold_1pct 
        binary_5pct <- fav_raster >= threshold_5pct
        
        # Plot
        plot(fav_raster, main = paste0("Favourability (", modelmethod,")"))
        plot(binary_1pct, main = paste0("Binary map ",modelmethod," (1% threshold)"))
        plot(binary_5pct, main = paste0("Binary map ",modelmethod," (5% threshold)"))
        
        #Store
        modeloutput[[modelmethod]]<-list(fav_raster=fav_raster,
                                         binary1pct=binary_1pct,
                                         binary5pct=binary_5pct,
                                         model=method_model)
        
        rm(fav_raster, binary_1pct, binary_5pct, method_model)
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
    cat("Top 5 models by variance along PC1:\n")
    print(top5_models)
    
    # Step 6: Subset fav_stack to top 5 layers
    top5_stack <- subset(fav_stack, top5_models)
    
    # Step 7: Compute pixel-wise median = consensus model
    consensus_habitat <- app(top5_stack, median)
    
    # Step 8: Plot result
    plot(consensus_habitat, main = "Consensus (Median of Top 5 Models by Variance on PC1)")

    
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
    
    # Plot (optional)
    plot(clim_hab_binary_1pct, main = "Binary Map (1% Threshold)")
    plot(clim_hab_binary_5pct, main = "Binary Map (5% Threshold)")
    
    
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
    
    
    #BOYCE
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
    
    #Prepare title and species extension
    species_title <- gsub("_", " ", first_two_words)
    rest_of_name <- if (grepl("^\\S+\\s+\\S+$", species)) "" else sub("^\\S+\\s+\\S+\\s+", "", species)
    
    # Export PDFs with and without occurrences plotted
    for (occs in list(NULL, eu_occ)){
      exportPDF(predictions = clim_hab,
                taxonName = first_two_words,
                nameExtension= rest_of_name,
                dataType = "Suit",
                taxonNameTitle = species_title,
                taxonKey = taxonkey,
                scenario = "hist",
                regionName = "EU",
                returnPredictions = FALSE,
                returnPNG = FALSE,
                occ_data=occs,
                exportPNG=TRUE,
                LabelValue= boyce_val,
                LabelName="Boyce index")
    }
    

    #------------------------------------------
    #------------ Create binary map -----------
    #------------------------------------------
    binary_maps <- list(
      "0.05" = list(
        binary_raster = clim_hab_binary_5pct,
        EU_sensitivity = final_sensitivity_5pct,
        boyce = boyce_val,
        mtp_pct = "5%"
      ),
      "0.01" = list(
        binary_raster = clim_hab_binary_1pct,
        EU_sensitivity = final_sensitivity_1pct,
        boyce = boyce_val,
        mtp_pct = "1%"
      )
    )
    
    for (pct in seq_along(binary_maps)){
      
      binary_map_pct <- binary_maps[[pct]]$binary_raster
      EU_sensitivity <- binary_maps[[pct]]$EU_sensitivity
      mtp_pct <- binary_maps[[pct]]$mtp_pct
      #boyce_ind <- binary_maps[[pct]]$boyce
      
      # export as PDF and PNG with and without occurrences plotted and return as PNG
      for (occs in list(NULL, eu_occ)){
        exportPDF(predictions = binary_map_pct,
                  taxonName = first_two_words,
                  nameExtension= rest_of_name,
                  dataType = "Binary",
                  taxonNameTitle = species_title,
                  taxonKey = taxonkey,
                  scenario = "hist",
                  regionName = "EU",
                  returnPredictions = FALSE,
                  returnPNG = FALSE,
                  occ_data=occs,
                  exportPNG=TRUE,
                  LabelValue= mtp_pct,
                  LabelName="MTP threshold",
                  Label2Value=round(EU_sensitivity,3),
                  Label2Name="Sensitivity")
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
                    boyce_result = boyce_result           # full object from ecospat.boyce()
    )
   
    #Save eumodel as .qs file
    qs::qsave(eumodel, 
              here::here(species_folder, paste0("EU_model_",first_two_words,"_",taxonkey,".qs"))) 
    
    
    #--------------------------------------------
    #- ------ Save rasters-----------------------
    #--------------------------------------------
    #Define raster path
    fullstack_be_file<- file.path(raster_folder,"Interim", paste0("Fullstack_be_",first_two_words,"_",taxonkey,".tif"))
    
    #Save locally because if in .qs file parts of the metadata will be stored in a Temp file that will be erased over time
    terra::writeRaster(fullstack_be, filename =fullstack_be_file, overwrite = TRUE)

    
    #--------------------------------------------
    #-------- End of loop -----------------------
    #--------------------------------------------
    print(paste("European model has been created for", species))
    rm(list = setdiff(ls(), c("p", "factorVars","accuracyStats", "findThresh", "projectname", "generate_pseudoabs", "create_folder", "country_name", "country_ext", "country_vector", "euboundary", "habitat_stack", "rmiclimpreds", "accepted_taxonkeys", "taxa_info", "key", "add.occ", "confidenceMaps", "classConformalPrediction","extractVals","GetLength","get.confidence", "exportPDF")))
  }
})


#--------------------------------------------
#----------- Clean R environment ------------
#--------------------------------------------
rm(list = ls())

