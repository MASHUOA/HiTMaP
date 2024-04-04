
install.packages("BiocManager")
install.packages("devtools")
if (!require(plotly)) devtools::install_github("ropensci/plotly@async")
if (!require(htmlwidgets)) devtools::install_github("ramnathv/htmlwidgets")
install.packages("tidyverse")
install.packages("reshape2")
library(devtools)
BiocManager::install(c("EBImage","ChemmineR","matter","Cardinal"))
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS" = "true")
devtools::install_github(
  "kuwisdelu/Cardinal",
  force=TRUE,
  upgrade="always",
  verbose=TRUE
)
devtools::install_github(
  "MASHUOA/HiTMaP",
  force=TRUE,
  upgrade="always",
  verbose=TRUE
)
