#--------------------------------------------
#-----------  Load packages  ----------------
#--------------------------------------------
packages <- c("curl"
)

for(package in packages) {
  print(package)
  if( ! package %in% rownames(installed.packages()) ) { install.packages( package ) }
  library(package, character.only = TRUE)
}



#--------------------------------------------
#--------- Source helper functions ----------
#--------------------------------------------
source("./src/helper_functions.R")


#--------------------------------------------
