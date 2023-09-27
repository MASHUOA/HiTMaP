if (!require(shiny)) install.packages("shiny")
library(shiny)

if (!require(future)) install.packages("future")
library(future)

if (!require(DT)) devtools::install_github("rstudio/DT")
library(DT)

if (!require(data.table)) install.packages("data.table")
library(data.table)

if (!require(shinyFiles)) install.packages("shinyFiles")
library(shinyFiles)

if (!require(plotly)) devtools::install_github("ropensci/plotly@async")
library(plotly)

if (!require(htmlwidgets)) devtools::install_github("ramnathv/htmlwidgets")
library(htmlwidgets)

if (!require(promises)) install.packages("promises")
library(promises)

if (!require(shinythemes)) install.packages("shinythemes")
library(shinythemes)
