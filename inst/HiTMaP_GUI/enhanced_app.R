#' Enhanced HiTMaP Application with Pipeline Integration
#' 
#' Main application file that launches the enhanced GUI with pipeline functionality

# Load required packages
if (!requireNamespace("shiny", quietly = TRUE)) stop("Missing GUI dependency: shiny")
library(shiny)

if (!requireNamespace("shinydashboard", quietly = TRUE)) stop("Missing GUI dependency: shinydashboard")
library(shinydashboard)

if (!requireNamespace("shinyjs", quietly = TRUE)) stop("Missing GUI dependency: shinyjs")
library(shinyjs)

if (!requireNamespace("DT", quietly = TRUE)) stop("Missing GUI dependency: DT")
library(DT)

if (!requireNamespace("plotly", quietly = TRUE)) stop("Missing GUI dependency: plotly")
library(plotly)

if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Missing GUI dependency: ggplot2")
library(ggplot2)

if (!requireNamespace("future", quietly = TRUE)) stop("Missing GUI dependency: future")
library(future)

if (!requireNamespace("data.table", quietly = TRUE)) stop("Missing GUI dependency: data.table")
library(data.table)

if (!requireNamespace("magrittr", quietly = TRUE)) stop("Missing GUI dependency: magrittr")
library(magrittr)

if (!requireNamespace("shinyFiles", quietly = TRUE)) stop("Missing GUI dependency: shinyFiles")
library(shinyFiles)

if (!requireNamespace("jsonlite", quietly = TRUE)) stop("Missing GUI dependency: jsonlite")
library(jsonlite)

if (!requireNamespace("zip", quietly = TRUE)) stop("Missing GUI dependency: zip")
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
