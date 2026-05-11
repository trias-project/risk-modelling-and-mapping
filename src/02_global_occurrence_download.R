#--------------------------------------------
#-----------------Load packages--------------
#--------------------------------------------
packages <- c("rgbif", "dplyr", "purrr", "assertthat", "readr", "here", "retry", 
              "CoordinateCleaner", "remotes", "stringr")

installed_packages <- installed.packages() |>
  as.data.frame()

for (package in packages) {
  print(package)
  if (!package %in% rownames(installed.packages()) ) { 
    install.packages( package ) 
    }
  library(package, character.only = TRUE)
}

# Install correct version of qs
req_qs_version <- "0.27.3"

if (!"qs" %in% installed_packages$Package){
  warning("qs is not installed => installing")
  remotes::install_version("qs", version = req_qs_version)
}else{
  qs_version <- installed_packages |>
    dplyr::filter(Package == "qs") |>
    dplyr::pull(Version)
  
  if(qs_version != req_qs_version){
    warning(paste("Version", qs_version, "of qs is installed, while", req_qs_version, 
                  "is required => installing correct version"))
    remotes::install_version("qs", version = req_qs_version)
  }else{
    print("Correct version of qs installed")
  }
}

library(qs)


#--------------------------------------------
#- Load helper functions and configurations -
#--------------------------------------------
source(here::here("src", "helper_functions.R"))
source(here::here("src", "00_configurations.R"))


#--------------------------------------------
#--------- Create project folder ------------
#--------------------------------------------
create_folder(here::here("data", "projects", project), project)


#--------------------------------------------
#-----------Retrieve GBIF taxonkeys----------
#--------------------------------------------
# Match species names with the GBIF backbone, retrieve taxon keys from GBIF when a match is found
taxon_df <- as.data.frame(species_to_model)

mapped_taxa <- purrr::map_dfr(
  taxon_df$species_to_model,
  ~ {
    tryCatch(
      {
        # Add a small delay to avoid API misses
        Sys.sleep(0.2)
        
        data <- rgbif::name_backbone(name = .x,
                                     curlopts=list(http_version=2))
        if (length(data) == 0) {
          stop("No match with the GBIF backbone found")
        }
        data
      },
      error = function(e) {
        NULL
      }
    )
  }
)

#Make sure that only species info is stored as it is possible that genus information is captured when the species part of the name is not clear
# mapped_taxa <- mapped_taxa %>%
#   dplyr::filter(rank == "SPECIES")

#Make sure that all species were mapped to the GBIF backbone, if not an error will appear indicating which species are missing
assertthat::assert_that(
  nrow(mapped_taxa) == length(species_to_model),
  msg = paste0("The following species could not be found in the GBIF backbone taxonomy: ",
               paste(species_to_model[!sapply(species_to_model, function(x) any(grepl(x, mapped_taxa$scientificName)))], collapse = ", "))
)

not_accepted <- mapped_taxa %>%
  dplyr::filter(status != "ACCEPTED")

if (nrow(not_accepted) != 0) {
  warning(paste0("The following species do not have an accepted taxonomic status in the GBIF backbone: ",paste(unique(not_accepted$scientificName), collapse = ", "),". Their corresponding accepted species names will be used for downloading occurrence data.")
  )
} else {
  paste0("All species are accepted taxa in the GBIF backbone 🎉")
}

#Extract taxonkeys of each species, for synonyms the acceptedUsageKey is stored
accepted_taxonkeys <- mapped_taxa %>%
  dplyr::filter(status == "ACCEPTED") %>%
  dplyr::pull(usageKey)

if (nrow(not_accepted > 0)) {
  synonym_taxonkeys <- mapped_taxa %>%
    dplyr::filter(status != "ACCEPTED") %>%
    dplyr::pull(acceptedUsageKey)
  
  accepted_taxonkeys <- c(accepted_taxonkeys, synonym_taxonkeys)
}

#Keep unique accepted taxonkeys
accepted_taxonkeys <- unique(accepted_taxonkeys)


#--------------------------------------------
#-----------Define download settings---------
#--------------------------------------------
#All basis of record types, except `FOSSIL SPECIMEN` and `LIVING SPECIMEN`, which can have misleading location information (e.g. location of captive animal).
basis_of_record <- c(
  "OBSERVATION", 
  "HUMAN_OBSERVATION",
  "MATERIAL_SAMPLE",
  "PRESERVED_SPECIMEN", 
  "UNKNOWN", 
  "MACHINE_OBSERVATION",
  "OCCURRENCE"
)

#Time period
year_begin <- 1971
year_end <- 2025

#Only georeferenced points
hasCoordinate <- TRUE


#--------------------------------------------
#---------------Perform download-------------
#--------------------------------------------
gbif_user  <- Sys.getenv("GBIF_USER",   unset = NA)
gbif_pwd   <- Sys.getenv("GBIF_PWD",    unset = NA)
gbif_email <- Sys.getenv("GBIF_EMAIL",  unset = NA)

gbif_download_key <-  rgbif::occ_download(
  pred_in("taxonKey", accepted_taxonkeys),
  pred_in("basisOfRecord", basis_of_record),
  pred_gte("year", year_begin),
  pred_lte("year", year_end),
  pred("hasCoordinate", hasCoordinate),
  pred("occurrenceStatus", "PRESENT"),
  pred("hasGeospatialIssue", FALSE), #Remove default geospatial issues
  pred_or(
    pred_lte("coordinateUncertaintyInMeters",5000),
    pred_isnull("coordinateUncertaintyInMeters")),
  user  =  gbif_user,
  pwd   = gbif_pwd,
  email = gbif_email,
  curlopts = list(http_version = 2)
)


rgbif::occ_download_wait(gbif_download_key)#Check download status


#--------------------------------------------
#--------------Retrieve download-------------
#--------------------------------------------
rgbif::occ_download_get(gbif_download_key, path = here::here("data"), overwrite = TRUE)
metadata <- rgbif::occ_download_meta(key = gbif_download_key)
gbif_download_key <- metadata$key

#extract_GBIF_occurrence
data.path <- here::here("data", gbif_download_key)
unzip(paste0(data.path,".zip"),exdir = data.path, overwrite = TRUE)
global <- as.data.frame(data.table::fread(paste0(data.path,"/occurrence.txt"),
                                          select=c("acceptedTaxonKey","acceptedScientificName", "decimalLatitude", "decimalLongitude", "kingdom", "phylum", "class","order", "genus", "coordinateUncertaintyInMeters", "identificationVerificationStatus"),
                                          header = TRUE))


#--------------------------------------------
# Process ambiguous synonyms
#--------------------------------------------

#Get unique taxonkeys that are not part of the accepted taxonkeys (ambiguous keys)
ambiguous<-global%>%
  dplyr::filter(!acceptedTaxonKey %in% accepted_taxonkeys)%>%
  dplyr::select(acceptedTaxonKey, acceptedScientificName)%>%
  dplyr::distinct()

if (nrow(ambiguous) > 0) {
  #Map these with the GBIF backbone
  mapped_ambiguous<- purrr::map_dfr(
    ambiguous$acceptedScientificName,
    ~ {
      tryCatch(
        {
          # Add a small delay to avoid API misses
          Sys.sleep(0.2)
          
          data <- rgbif::name_backbone(name = .x)
          if (length(data) == 0) {
            stop("No match with the GBIF backbone found")
          }
          data
        },
        error = function(e) {
          NULL
        }
      )
    }
  )
  
  #Keep original acceptedScientificName and the species it was mapped to
  mapped_ambiguous <- mapped_ambiguous %>% 
    dplyr::select(verbatim_name, species)
  
  #Map the species-level against the GBIF backbone
  mapped_ambiguous_species <- purrr::map_dfr(
    mapped_ambiguous$species,
    ~ {
      tryCatch(
        {
          # Add a small delay to avoid API misses
          Sys.sleep(0.2)
          
          data <- rgbif::name_backbone(name = .x)
          if (length(data) == 0) {
            stop("No match with the GBIF backbone found")
          }
          data
        },
        error = function(e) {
          NULL
        }
      )
    }
  )
  
  #Create a df with the following columns:
  #verbatim_name= original acceptedScientificName in df 'global'
  #usageKey = taxonKey of the mapped species
  #scientificName = scientific name of the mapped species
  mapped<- mapped_ambiguous %>%
    dplyr::select(species, verbatim_name) %>%
    left_join(mapped_ambiguous_species, by = c("species" = "verbatim_name"))%>%
    dplyr::select(verbatim_name, usageKey, scientificName)
  
  #Add this info to df 'global'
  global<-left_join(global, mapped, by = c("acceptedScientificName" = "verbatim_name"))
  
  #Overwrite acceptedScientificName and acceptedTaxonKey if necessary
  global<-global %>%
    dplyr::mutate(
      acceptedScientificName = ifelse(!acceptedTaxonKey %in% accepted_taxonkeys,
                                      scientificName, 
                                      acceptedScientificName),
      acceptedTaxonKey = ifelse(!acceptedTaxonKey %in% accepted_taxonkeys, 
                                usageKey,
                                acceptedTaxonKey))%>%
    dplyr::select(-c(usageKey, scientificName))%>%
    dplyr::filter(acceptedTaxonKey %in% accepted_taxonkeys) 
}

#--------------------------------------------
#-------- Filter global occurrences----------
#--------------------------------------------
#remove unverified records
identificationVerificationStatus_to_discard <- c( "unverified",
                                                  "unvalidated",
                                                  "not validated",
                                                  "under validation",
                                                  "not able to validate",
                                                  "control could not be conclusive due to insufficient knowledge",
                                                  "1",
                                                  "uncertain",
                                                  "unconfirmed",
                                                  "douteux",
                                                  "invalide",
                                                  "non r\u00E9alisable",
                                                  "verification needed" ,
                                                  "probable",
                                                  "unconfirmed - not reviewed",
                                                  "validation requested",
                                                  "unconfirmed - plausible")

#enter value for max coordinate uncertainty in meters, default = 5000
#Reasoning: there are several invasive species datasets with default uncertainty of 5km
global.occ<-global %>%
  dplyr::filter(acceptedTaxonKey%in%accepted_taxonkeys) %>%   
  dplyr::filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters <= 5000) %>% 
  dplyr::filter(!str_to_lower(identificationVerificationStatus) %in% identificationVerificationStatus_to_discard)

#Remove coordinates that for both lon and lat values, have less than 4 decimal places
global.occ$lon_dplaces<-sapply(global.occ$decimalLongitude, function(x) decimalplaces(x))
global.occ$lat_dplaces<-sapply(global.occ$decimalLatitude, function(x) decimalplaces(x))
rows_to_na <- which(global.occ$lon_dplaces < 4 & global.occ$lat_dplaces < 4)
global.occ[rows_to_na, ] <- NA
global.occ<-global.occ[ which(!is.na(global.occ$lon_dplaces)),]
global.occ<-within(global.occ,rm("lon_dplaces","lat_dplaces")) 


#--------------------------------------------
#------------ Define species group-----------
#--------------------------------------------
fish_orders <- c(
  "Acipenseriformes", "Albuliformes", "Amiiformes", "Anguilliformes", "Atheriniformes",
  "Batrachoidiformes", "Beloniformes", "Beryciformes", "Blenniiformes", "Carangiformes",
  "Characiformes", "Chimaeriformes", "Clupeiformes", "Cypriniformes", "Cyprinodontiformes",
  "Elopiformes", "Esociformes", "Gadiformes", "Gasterosteiformes", "Gonorynchiformes",
  "Gymnotiformes", "Holocentriformes", "Istiophoriformes", "Lampriformes", "Lophiiformes",
  "Mugiliformes", "Myliobatiformes", "Myxiniformes", "Notacanthiformes", "Ophidiiformes",
  "Osmeriformes", "Osteoglossiformes", "Perciformes", "Percopsiformes", "Petromyzontiformes",
  "Pleuronectiformes", "Polypteriformes", "Pristiformes", "Rajiformes", "Salmoniformes",
  "Scorpaeniformes", "Siluriformes", "Squaliformes", "Stomiiformes", "Syngnathiformes",
  "Tetraodontiformes", "Torpediniformes", "Trachiniformes", "Zeiformes"
)

global.occ <- global.occ %>%
  dplyr::mutate(Group = case_when(
    kingdom == "Plantae" ~ "Plants",
    class == "Aves" ~ "Birds",
    phylum == "Mollusca" ~ "Molluscs",
    class == "Amphibia" ~ "Amphibians",
    class == "Mammalia" ~ "Mammals",
    class %in% c("Crocodylia", "Testudines", "Sphenodontia", "Squamata") ~ "Reptiles",
    class == "Malacostraca" ~ "Malacostraca",
    class == "Insecta" ~ "Insects",
    order %in% fish_orders ~ "Fish",
    TRUE ~ NA_character_
  ))


#--------------------------------------------
#-------Prepare occurrence dataset-----------
#--------------------------------------------
global.occ.LL<-global.occ%>%
  dplyr::rename(species= acceptedScientificName)%>%
  dplyr::select(decimalLongitude, decimalLatitude, species, acceptedTaxonKey, Group, coordinateUncertaintyInMeters) 
rm(global.occ, global)


#--------------------------------------------
#-----------Do coordinate cleaning-----------
#--------------------------------------------
# Clean coordinates based on their proximity to country centroids, capitals, biodiversity institutions, GBIF headquarters, and the 0/0 point
cleaned<-global.occ.LL%>%
  CoordinateCleaner::cc_cen(buffer=100) %>% # remove points within a buffer of 100m around country centroids, default 1km
  CoordinateCleaner::cc_cap(buffer=100) %>% # remove capitals centroids (buffer 100m), default 10km
  CoordinateCleaner::cc_inst(buffer=100) %>% # remove zoo and herbaria records buffer of 100 m around biodiversity institutes, default 100m
  CoordinateCleaner::cc_gbif(buffer=100)%>% #remove around GBIF headquarters in Copenhagen (buffer 100m), default 100m
  CoordinateCleaner::cc_zero() #Remove around the 0/0 point (buffer 0.5 degrees)


#---------------------------------------------
#--- Prepare data at 1km coord uncertainty ---
#---------------------------------------------
#These data will be used for the European model
cleaned_1km<-cleaned%>%
  dplyr::filter(is.na(coordinateUncertaintyInMeters)| coordinateUncertaintyInMeters <= 1000)


#--------------------------------------------
#------------------Save data-----------------
#--------------------------------------------
# Number of unique taxa
unique_keys <- unique(cleaned$acceptedTaxonKey)
unique_names <- unique(cleaned$species)
n <- length(unique_keys)

#Create dataset taxa_info containing scientific name, canonical name, taxonkeys, gbif download key,...
taxa_info <- data.frame(acceptedTaxonKey = unique_keys,
                        acceptedScientificName = unique_names,
                        year_begin = rep_len(metadata[["request"]][["predicate"]][["predicates"]][[3]][["value"]], n),
                        year_end = rep_len(metadata[["request"]][["predicate"]][["predicates"]][[4]][["value"]], n),
                        gbif_download_key = rep_len(gbif_download_key, n),
                        gbif_download_created = rep_len(format(
                          strptime(metadata$created, "%Y-%m-%dT%H:%M:%S"),
                          "%Y-%m-%d %H:%M:%S"
                        ), n),
                        project = rep_len(project, n)
)

#Save occurrence data as .qs file and taxa info as .csv
occurrences <-list(cleaned_1km = cleaned_1km,
                   cleaned = cleaned)

qs::qsave(occurrences, paste0("./data/projects/",project,"/",project,"_processed_occurrences.qs"))
write.csv2(taxa_info, paste0("./data/projects/",project,"/",project,"_taxa_info.csv"), row.names = FALSE)


#--------------------------------------------
#---- Clean up environment and local disk----
#--------------------------------------------
# Remove the zipped folder
suppressWarnings(file.remove(paste0(data.path, ".zip"), full.names = TRUE))

# Remove the unzipped folder 
unlink(data.path, recursive = TRUE)

# Clean R environment
rm(list = ls())
