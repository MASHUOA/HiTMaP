install.packages("remotes")
if (!require(plotly)) devtools::install_github("ropensci/plotly@async")
if (!require(htmlwidgets)) devtools::install_github("ramnathv/htmlwidgets")
install.packages("tidyverse")
install.packages("reshape2")
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
