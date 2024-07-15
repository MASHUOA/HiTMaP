
install.packages("BiocManager")
install.packages("systemfonts")
install.packages("remotes")
install.packages("devtools")
library(devtools)
if (!require(plotly)) remotes::install_github("ropensci/plotly@async")
if (!require(htmlwidgets)) remotes::install_github("ramnathv/htmlwidgets")
install.packages("tidyverse")
install.packages("reshape2")
BiocManager::install(c("EBImage","ChemmineR","matter","Cardinal"),ask=F)
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
remotes::install_github(
  "MASHUOA/HiTMaP",
  force=TRUE,
  upgrade="always",
  verbose=TRUE
)
