#install the git package
install.packages("remotes")
install.packages("config")
install.packages("textshaping")
install.packages("ragg")
install.packages("pkgdown")
install.packages("devtools")
install.packages("BiocManager")
install.packages("systemfonts")
if (!require(plotly)) remotes::install_github("ropensci/plotly@async")
if (!require(htmlwidgets)) remotes::install_github("ramnathv/htmlwidgets")
install.packages("tidyverse")
install.packages("reshape2")

