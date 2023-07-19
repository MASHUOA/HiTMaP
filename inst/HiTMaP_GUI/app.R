
if (!require(shiny)) install.packages("shiny")
library(shiny)

if (!require(future)) install.packages("future")
library(future)

if (!require(DT)) devtools::install_github("rstudio/DT")
library(DT)

library(HiTMaP)

if (!require(data.table)) install.packages("data.table")
library(data.table)

if (!require(shinyFiles)) install.packages("shinyFiles")
library(shinyFiles)

if (!require(plotly)) devtools::install_github("ropensci/plotly@async")
library(plotly)

if (!require(htmlwidgets)) devtools::install_github("ramnathv/htmlwidgets")
library(htmlwidgets)

if (!require(shinythemes)) devtools::install_github("rstudio/shinythemes")
library(shinythemes)

HiTMaP:::Peptide_modification(retrive_ID=NULL,update_unimod=F)
plan(multiprocess)
options(shiny.maxRequestSize = 3000*1024^2)







enableBookmarking("server")
shinyApp(ui = ui, server = server )