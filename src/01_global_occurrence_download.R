#--------------------------------------------
#----To do: specify project and species------
#--------------------------------------------
#specify project name
project<-"Test_Frédérique"

# specify the scientific name of the species to be modelled
species<-c("Elodea densa","Koenigia polystachya", "Hydrocharis laevigata")


#--------------------------------------------
#-----------------Load packages--------------
#--------------------------------------------
library(purrr)
library(rgbif)
library(dplyr)
library(assertthat)
library(readr)
library(here)
library(trias)
library(qs)


#--------------------------------------------
#--Create a project folder and a raw folder--
#--------------------------------------------
# Define the folder paths
project_path <- file.path("./data/projects",project)
raw_path<-file.path("./data/raw")

# Check and create each folder if necessary
create_folder(project_path, project)
create_folder(raw_path, "raw")

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
mapped_taxa<-mapped_taxa %>%
  dplyr::filter(rank =="SPECIES")

#Make sure that all species were mapped to the GBIF backbone, if not an error will appear indicating which species are missing
assertthat::assert_that(nrow(mapped_taxa)==length(species),
                        msg=paste0("The following species could not be found in the GBIF backbone taxonomy: "
                                   ,species[!sapply(species, function(x) any(grepl(x,mapped_taxa$scientificName)))])
)

not_accepted <- mapped_taxa %>%
  dplyr::filter(status !="ACCEPTED")

if (nrow(not_accepted)!=0) {
  warning(paste0("The following species do not have an accepted taxonomic status in the GBIF backbone: ",paste(unique(not_accepted$scientificName), collapse=", "),". Their corresponding accepted species names will be used for downloading occurrence data.")
  )
} else {
  paste0("All species are accepted taxa in the GBIF backbone 🎉")
}

#Extract taxonkeys of each species, for synonyms the acceptedUsageKey is stored
accepted_taxonkeys<-mapped_taxa %>%
  dplyr::filter(status =="ACCEPTED")%>%
  pull(usageKey)

if(nrow(not_accepted!=0)){
  synonym_taxonkeys<-mapped_taxa %>%
    dplyr::filter(status !="ACCEPTED")%>%
    pull(acceptedUsageKey)
  
  accepted_taxonkeys<-c(accepted_taxonkeys, synonym_taxonkeys)
}

#Create taxa_input file
taxa_input_file<-mapped_taxa[,2:3]


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
year_end <-2010

#Only georeferenced points
hasCoordinate <- TRUE


#--------------------------------------------
#---------------Perform download-------------
#--------------------------------------------
#Note that GBIF credentials are required
gbif_download_key <- occ_download(
  pred_in("taxonKey", accepted_taxonkeys),
  pred_in("basisOfRecord", basis_of_record),
  pred_gte("year", year_begin),
  pred_lte("year", year_end),
  pred("hasCoordinate", hasCoordinate),
  user = rstudioapi::askForPassword("GBIF username"),
  pwd = rstudioapi::askForPassword("GBIF password"),
  email = rstudioapi::askForPassword("Email address for notification")
)

occ_download_wait(gbif_download_key)#Check download status


#--------------------------------------------
#--------------Retrieve download-------------
#--------------------------------------------
#gbif_download_key<-"0076914-240626123714530"
#gbif_download_key<-"0064066-240626123714530" #4 species
occ_download_get(gbif_download_key, path = here("data","raw"), overwrite=TRUE)
metadata <- occ_download_meta(key = gbif_download_key)
gbif_download_key<-metadata$key

#extract_GBIF_occurrence
raw.path<- here("data", "raw", gbif_download_key)
unzip(paste0(raw.path,".zip"),exdir=raw.path, overwrite=TRUE)
global<-as.data.frame(data.table::fread(paste0(raw.path,"/occurrence.txt"),header=TRUE))
global<-select(global, c(speciesKey,acceptedScientificName, decimalLatitude, decimalLongitude, kingdom, phylum, class, coordinateUncertaintyInMeters, identificationVerificationStatus))


#--------------------------------------------
#------------------Save data-----------------
#--------------------------------------------
#Create dataset taxa_info containing scientific name, canonical name, taxonkeys, gbif download key,...
taxa_info<-taxa_input_file %>%
  cbind(accepted_taxonkeys) %>%
  as.data.frame() %>%  
  mutate(gbif_download_key = gbif_download_key,
         gbif_download_created = format(strptime(metadata$created, "%Y-%m-%dT%H:%M:%S"), "%Y-%m-%d %H:%M:%S"),
         project = project)

#Save occurrence data as .qs file and taxa info as .csv
qsave(global, paste0("./data/projects/",project,"/",project,"_occurrences.qs"))
write.csv2(taxa_info, paste0("./data/projects/",project,"/",project,"_taxa_info.csv"), row.names=FALSE)


#--------------------------------------------
#---- Clean up environment and local disk----
#--------------------------------------------
# Remove the zipped folder
suppressWarnings(file.remove(paste0(raw.path, ".zip"), full.names = TRUE))

# Remove the unzipped folder 
unlink(raw.path, recursive = TRUE)

# Clean R environment
rm(list = ls())


