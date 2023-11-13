install.packages("BiocManager")
install.packages("remotes")
install.packages("config")
install.packages("textshaping")
install.packages("ragg")
install.packages("pkgdown")
install.packages("devtools")
install.packages("promises")
install.packages("magrittr")
if (!require(plotly)) devtools::install_github("ropensci/plotly@async")
if (!require(htmlwidgets)) devtools::install_github("ramnathv/htmlwidgets")
install.packages("tidyverse")
install.packages("reshape2")
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
