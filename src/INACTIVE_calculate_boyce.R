#------------------------------------------
#-Calculate Boyce index to display on map- climate model --
#------------------------------------------
#Extract all raster values (excluding NAs)
# all_suit_vals <- values(consensus_median)
# all_suit_vals <- all_suit_vals[!is.na(all_suit_vals)]

#Extract suitability values at occurrence locations
# occ_suit_vals <- terra::extract(consensus_median, vect(global.occ.sf))[,2]
# occ_suit_vals <- occ_suit_vals[!is.na(occ_suit_vals)]

#Compute Boyce only if there are enough occurrences
# if (length(occ_suit_vals) > 0) {
#   boyce_result <- ecospat.boyce(
#     fit = all_suit_vals,
#     obs = occ_suit_vals,
#     nclass = 0
#   )
#   boyce_val <- round(boyce_result$cor, 3)
# } else {
#   warning(paste("No EU occurrences available to calculate Boyce index for", species))
#   boyce_val <- "NA (no EU data)"
# }


#------------------------------------------
#-Calculate Boyce index to display on map- habitat model--
#------------------------------------------
#Extract all raster values (excluding NAs)
all_suit_vals <- values(consensus_habitat)
all_suit_vals <- all_suit_vals[!is.na(all_suit_vals)]

#Extract suitability values at occurrence locations
occ_suit_vals <- terra::extract(consensus_habitat, vect(eu_occ))[,2]
occ_suit_vals <- occ_suit_vals[!is.na(occ_suit_vals)]

#Compute Boyce only if there are enough occurrences
if (length(occ_suit_vals) > 0) {
  boyce_result <- ecospat.boyce(
    fit = all_suit_vals,
    obs = occ_suit_vals,
    nclass = 0
  )
  boyce_val <- round(boyce_result$cor, 3)
} else {
  warning(paste("No EU occurrences available to calculate Boyce index for", species))
  boyce_val <- "NA (no EU data)"
}



#------------------------------------------
#-Calculate Boyce index to display on map- ensemble model
#------------------------------------------
#Extract all raster values (excluding NAs)
all_suit_vals <- values(clim_hab)
all_suit_vals <- all_suit_vals[!is.na(all_suit_vals)]

#Extract suitability values at occurrence locations
occ_suit_vals <- terra::extract(clim_hab, vect(eu_occ))[,2]
occ_suit_vals <- occ_suit_vals[!is.na(occ_suit_vals)]

#Compute Boyce only if there are enough occurrences
if (length(occ_suit_vals) > 0) {
  boyce_result <- ecospat.boyce(
    fit = all_suit_vals,
    obs = occ_suit_vals,
    nclass = 0
  )
  boyce_val <- round(boyce_result$cor, 3)
} else {
  warning(paste("No EU occurrences available to calculate Boyce index for", species))
  boyce_val <- "NA (no EU data)"
}