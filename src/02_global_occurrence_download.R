#--------------------------------------------
#-----------------Load packages--------------
#--------------------------------------------
packages <- c("rgbif", "dplyr", "purrr", "assertthat", "readr", "here", "qs", "retry")

for (package in packages) {
  print(package)
  if (!package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}


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
mapped_taxa <- mapped_taxa %>%
  dplyr::filter(rank == "SPECIES")

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
gbif_download_key <-  rgbif::occ_download(
  pred_in("taxonKey", accepted_taxonkeys),
  pred_in("basisOfRecord", basis_of_record),
  pred_gte("year", year_begin),
  pred_lte("year", year_end),
  pred("hasCoordinate", hasCoordinate),
  user  =  rstudioapi::askForPassword("GBIF username"),
  pwd   = rstudioapi::askForPassword("GBIF password"),
  email = rstudioapi::askForPassword("Email address for notification"),
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
global <- as.data.frame(data.table::fread(paste0(data.path,"/occurrence.txt"),header = TRUE))

global <- dplyr::select(global, c(acceptedTaxonKey,acceptedScientificName, decimalLatitude, decimalLongitude, kingdom, phylum, class,order, genus, coordinateUncertaintyInMeters, identificationVerificationStatus))


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
#------------------Save data-----------------
#--------------------------------------------
# Number of unique taxa
unique_keys <- unique(global$acceptedTaxonKey)
unique_names <- unique(global$acceptedScientificName)
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
qs::qsave(global, paste0("./data/projects/",project,"/",project,"_occurrences.qs"))
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