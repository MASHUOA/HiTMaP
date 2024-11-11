library(devtools)
library(remotes)
remotes::install_github("sneumann/Rdisop")
BiocManager::install(c("EBImage","ChemmineR","Cardinal"))
BiocManager::install(c("XVector", "Biostrings", "KEGGREST","cleaver"))
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github(
  "MASHUOA/HiTMaP",
  force=TRUE,
  upgrade="always",
  verbose=TRUE
)
