#' Launch Enhanced HiTMaP GUI
#' 
#' Convenience function to launch the enhanced GUI with pipeline functionality

#' Launch Enhanced HiTMaP GUI
#' 
#' @param port Port number for the Shiny application (default: 3838)
#' @param host Host address (default: "0.0.0.0" for all interfaces)
#' @param launch_browser Whether to launch browser automatically (default: TRUE)
#' @param working_dir Working directory for data files (default: "~/expdata")
#' @return None (launches Shiny app)
#' 
#' @export
launch_enhanced_hitmap_gui <- function(port = 3838, 
                                      host = "0.0.0.0", 
                                      launch_browser = TRUE,
                                      working_dir = "~/expdata") {
  
  # Set up working directory
  if (!dir.exists(working_dir)) {
    dir.create(working_dir, recursive = TRUE)
    message(paste("Created working directory:", working_dir))
  }
  
  # Set global working directory
  assign("WorkingDir_global", working_dir, envir = globalenv())
  
  message("=== Starting Enhanced HiTMaP GUI ===")
  message(paste("Working directory:", working_dir))
  message(paste("URL: http://localhost:", port, sep = ""))
  message("Features available:")
  message("  - Pipeline-based workflow")
  message("  - Real-time monitoring")
  message("  - Advanced parameter configuration")
  message("  - Interactive results visualization")
  message("  - Checkpoint and resume functionality")
  
  # Check for required packages
  required_packages <- c("shiny", "shinydashboard", "shinyjs", "DT", "plotly", 
                        "ggplot2", "future", "data.table", "magrittr", 
                        "shinyFiles", "jsonlite", "zip")
  
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
  
  if (length(missing_packages) > 0) {
    message("Installing missing packages:", paste(missing_packages, collapse = ", "))
    install.packages(missing_packages)
  }
  
  # Set the working directory and launch
  setwd(working_dir)
  
  # Source the enhanced app
  gui_dir <- system.file("HiTMaP_GUI", package = "HiTMaP")
  
  if (gui_dir == "" || !file.exists(file.path(gui_dir, "enhanced_app.R"))) {
    # If package installation, use local files
    gui_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  }
  
  if (file.exists(file.path(gui_dir, "enhanced_app.R"))) {
    source(file.path(gui_dir, "enhanced_app.R"))
  } else {
    message("Enhanced GUI files not found. Using basic launcher...")
    
    # Basic launcher if enhanced files not available
    ui <- fluidPage(
      titlePanel("HiTMaP Enhanced GUI"),
      wellPanel(
        h3("GUI Loading..."),
        p("The enhanced HiTMaP GUI is being prepared."),
        p("If this takes too long, please check the installation.")
      )
    )
    
    server <- function(input, output, session) {
      showNotification("Enhanced GUI files not found. Please check installation.", 
                      type = "error", duration = NULL)
    }
    
    shinyApp(ui = ui, server = server, options = list(
      host = host,
      port = port,
      launch.browser = launch_browser
    ))
  }
}

#' Quick Launch Function
#' 
#' @export
hitmap_gui <- function() {
  launch_enhanced_hitmap_gui()
}

#' Launch with Development Settings
#' 
#' @export
launch_hitmap_dev <- function() {
  launch_enhanced_hitmap_gui(
    port = 3939,
    host = "127.0.0.1",
    launch_browser = TRUE,
    working_dir = "./dev_data"
  )
}

# Auto-launch if run directly
if (sys.nframe() == 0) {
  launch_enhanced_hitmap_gui()
}