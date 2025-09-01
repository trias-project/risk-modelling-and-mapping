#--------------------------------------------
#----To do: specify project and species------
#--------------------------------------------
#specify project name
project <- "PA prob & Alternative Treshold & Ensemble Boyce"

# specify the scientific name of the species to be modelled; see an example below
species <- c("Cabomba caroliniana",
           "Cipangopaludina chinensis",
           "Salvinia ×molesta",
           "Xenopus laevis",
           "Apalone ferox",
           "Muntiacus reevesi",
           "Myriophyllum aquaticum",
           "Vaccinium corymbosum",
           "Alternanthera philoxeroides",
           "Mustela vison",
           "Procyon lotor",
           "Faxonius virilis",
           "Faxonius immunis",
           "Obama nungara",
           "Procambarus clarkii",
           "Pacifastacus leniusculus",
           "Crassula helmsii",
           "Neogobius melanostomus",
           "Pseudorasbora parva",
           "Elaphe taeniura",
           "Heracleum mantegazzianum",
           "Ambrosia artemisiifolia",
           "Vespa velutina",
           "Vespa mandarinia")


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
#---------Load helper functions--------------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
#-------------Create folders-----------------
#--------------------------------------------
# Define the folder paths
folder_paths <- list(list("path" = file.path("./data/projects",project),
                        "name" = project),
                   list("path" = file.path("./data/raw"),
                        "name" = "raw")
)

# Check and create each folder if necessary
lapply(folder_paths, function(folder){
  create_folder(folder$path, folder$name)
})

#--------------------------------------------
#-----------Retrieve GBIF taxonkeys----------
#--------------------------------------------
# Match species names with the GBIF backbone, retrieve taxon keys from GBIF when a match is found
taxon_df <- as.data.frame(species)

mapped_taxa <- purrr::map_dfr(
  taxon_df$species,
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

#Make sure that only species info is stored as it is possible that genus information is captured when the species part of the name is not clear
mapped_taxa <- mapped_taxa %>%
  dplyr::filter(rank == "SPECIES")

#Make sure that all species were mapped to the GBIF backbone, if not an error will appear indicating which species are missing
assertthat::assert_that(
  nrow(mapped_taxa) == length(species),
  msg = paste0("The following species could not be found in the GBIF backbone taxonomy: ",
               paste(species[!sapply(species, function(x) any(grepl(x, mapped_taxa$scientificName)))], collapse = ", "))
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
year_end <- 2024

#Only georeferenced points
hasCoordinate <- TRUE


#--------------------------------------------
#---------------Perform download-------------
#--------------------------------------------
gbif_user  <- rstudioapi::askForPassword("GBIF username")
gbif_pwd   <- rstudioapi::askForPassword("GBIF password")
gbif_email <- rstudioapi::askForPassword("Email address for notification")

gbif_download_key <- retry::retry(
  {
    rgbif::occ_download(
      pred_in("taxonKey", accepted_taxonkeys),
      pred_in("basisOfRecord", basis_of_record),
      pred_gte("year", year_begin),
      pred_lte("year", year_end),
      pred("hasCoordinate", hasCoordinate),
      user  = gbif_user,
      pwd   = gbif_pwd,
      email = gbif_email
    )
  },
  when = "error",
  max_tries = 3,
  interval = 5
)

rgbif::occ_download_wait(gbif_download_key)#Check download status


#--------------------------------------------
#--------------Retrieve download-------------
#--------------------------------------------
rgbif::occ_download_get(gbif_download_key, path = here::here("data","raw"), overwrite = TRUE)
metadata <- rgbif::occ_download_meta(key = gbif_download_key)
gbif_download_key <- metadata$key

#extract_GBIF_occurrence
raw.path <- here::here("data", "raw", gbif_download_key)
unzip(paste0(raw.path,".zip"),exdir = raw.path, overwrite = TRUE)
global <- as.data.frame(data.table::fread(paste0(raw.path,"/occurrence.txt"),header = TRUE))
  
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
#Create dataset taxa_info containing scientific name, canonical name, taxonkeys, gbif download key,...
taxa_info<-data.frame(speciesKey=unique(global$speciesKey),
                      acceptedScientificName=unique(global$species),
                      year_begin=metadata[["request"]][["predicate"]][["predicates"]][[3]][["value"]],
                      year_end=metadata[["request"]][["predicate"]][["predicates"]][[4]][["value"]],
                      gbif_download_key = gbif_download_key,
                      gbif_download_created = format(strptime(metadata$created, "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d %H:%M:%S"),
                      project = project)

#Save occurrence data as .qs file and taxa info as .csv
qs::qsave(global, paste0("./data/projects/",project,"/",project,"_occurrences.qs"))
write.csv2(taxa_info, paste0("./data/projects/",project,"/",project,"_taxa_info.csv"), row.names = FALSE)


#--------------------------------------------
#---- Clean up environment and local disk----
#--------------------------------------------
# Remove the zipped folder
suppressWarnings(file.remove(paste0(raw.path, ".zip"), full.names = TRUE))

# Remove the unzipped folder 
unlink(raw.path, recursive = TRUE)

# Clean R environment
rm(list = ls())


