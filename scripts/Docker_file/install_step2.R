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
  verbose=TRUE
)
"https://bioconductor.statistik.tu-dortmund.de/packages/3.17/bioc/src/contrib/Archive/Cardinal/Cardinal_3.2.0.tar.gz"->packageurl
install.packages(packageurl, repos=NULL, type="source")