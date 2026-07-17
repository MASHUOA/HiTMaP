
if (!requireNamespace("shiny", quietly = TRUE)) stop("Missing GUI dependency: shiny")
library(shiny)

if (!requireNamespace("future", quietly = TRUE)) stop("Missing GUI dependency: future")
library(future)

if (!requireNamespace("DT", quietly = TRUE)) stop("Missing GUI dependency: DT")
library(DT)

library(HiTMaP)

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing GUI dependency: data.table")
library(data.table)

if (!requireNamespace("shinyFiles", quietly = TRUE)) stop("Missing GUI dependency: shinyFiles")
library(shinyFiles)

if (!requireNamespace("plotly", quietly = TRUE)) stop("Missing GUI dependency: plotly")
library(plotly)

if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Missing GUI dependency: htmlwidgets")
library(htmlwidgets)


HiTMaP:::Peptide_modification(retrive_ID=NULL,update_unimod=F)
plan(multicore)
options(shiny.maxRequestSize = 3000*1024^2)



if (!exists("WorkingDir_global", envir=globalenv()))  WorkingDir_global<<-"~/expdata"


enableBookmarking("server")
shinyApp(ui = ui, server = server )
