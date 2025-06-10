#' Enhanced HiTMaP Application with Pipeline Integration
#' 
#' Main application file that launches the enhanced GUI with pipeline functionality

# Load required packages
if (!require(shiny)) install.packages("shiny")
library(shiny)

if (!require(shinydashboard)) install.packages("shinydashboard")
library(shinydashboard)

if (!require(shinyjs)) install.packages("shinyjs")
library(shinyjs)

if (!require(DT)) install.packages("DT")
library(DT)

if (!require(plotly)) install.packages("plotly")
library(plotly)

if (!require(ggplot2)) install.packages("ggplot2")
library(ggplot2)

if (!require(future)) install.packages("future")
library(future)

if (!require(data.table)) install.packages("data.table")
library(data.table)

if (!require(magrittr)) install.packages("magrittr")
library(magrittr)

if (!require(shinyFiles)) install.packages("shinyFiles")
library(shinyFiles)

if (!require(jsonlite)) install.packages("jsonlite")
library(jsonlite)

if (!require(zip)) install.packages("zip")
library(zip)

# Load HiTMaP package
library(HiTMaP)

# Source the enhanced UI and server components
source("enhanced_ui.R")
source("enhanced_server.R")

# Initialize HiTMaP components
HiTMaP:::Peptide_modification(retrive_ID = NULL, update_unimod = FALSE)

# Set up parallel processing
plan(multisession)

# Configure Shiny options
options(shiny.maxRequestSize = 3000*1024^2)  # 3GB max file size

# Set up global working directory
if (!exists("WorkingDir_global", envir = globalenv())) {
  WorkingDir_global <<- "~/expdata"
  if (!dir.exists(WorkingDir_global)) {
    dir.create(WorkingDir_global, recursive = TRUE)
  }
}

# Enable bookmarking
enableBookmarking("server")

# Create the Shiny app
shinyApp(
  ui = enhanced_ui,
  server = enhanced_server,
  options = list(
    launch.browser = TRUE,
    host = "0.0.0.0",
    port = 3838
  )
)