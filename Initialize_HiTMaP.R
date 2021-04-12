install.packages("remotes") # allows HiTMaP to be installed from GitHub
install.packages("config")  # allows HiTMaP to be run from configuration file. TODO: add to DESCRIPTION
install.packages("enviPat") # TODO: add to DESCRIPTION (imports)
install.packages("fftwtools") # TODO: add to DESCRIPTION
BiocManager::install("EBImage") # TODO: add to DESCRIPTION
BiocManager::install(c("Cardinal", "EBImage", "ChemmineR")) # TODO: add to DESCRIPTION
library(remotes)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github(
  "MASHUOA/HiTMaP",
  auth_token ="cf6d877f8b6ada1865987b13f9e9996c1883014a", # TODO: hide token
  force=TRUE,
  upgrade="always" # suppresses prompt
)
warnings()
#Update all dependencies
#BiocManager::install(ask = F)
library(HiTMaP)
