


#install the git package
install.packages("remotes")
#library(devtools)
library(remotes)
install.packages("enviPat")
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
install.packages("fftwtools")
install.packages("")
BiocManager::install("EBImage")
BiocManager::install(c("Cardinal", "EBImage", "ChemmineR"))
remotes::install_github("MASHUOA/HiTMaP",auth_token ="cf6d877f8b6ada1865987b13f9e9996c1883014a",force=T)
3
no
#Update all dependencies
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
BiocManager::install(ask = F)
yes
library(HiTMaP)