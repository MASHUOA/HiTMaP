#install the git package
install.packages("remotes")
install.packages("config")
install.packages("textshaping")
install.packages("ragg")
install.packages("pkgdown")
install.packages("devtools")
BiocManager::install(c("EBImage","ChemmineR","Cardinal"))
library(remotes)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github(
  "MASHUOA/HiTMaP",
  force=TRUE,
  upgrade="always",
  verbose=TRUE,
  build=F
)
warnings()
#Update all dependencies
#BiocManager::install(ask = F)
library(HiTMaP)
