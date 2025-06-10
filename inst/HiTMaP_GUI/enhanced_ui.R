#' Enhanced HiTMaP GUI with Pipeline Integration
#' 
#' Updated UI that integrates the new pipeline system with the existing GUI

# Source the pipeline UI components
source("pipeline_ui.R")

# Enhanced UI with pipeline integration
enhanced_ui <- function(request) {
  if (!require(shiny)) install.packages("shiny")
  library(shiny)
  if (!require(shinyjs)) install.packages("shinyjs")
  library(shinyjs)
  if (!require(shinydashboard)) install.packages("shinydashboard")
  library(shinydashboard)
  
  dashboardPage(
    
    # Dashboard Header
    dashboardHeader(
      title = "HiTMaP - Enhanced Pipeline Interface",
      titleWidth = 350,
      
      # Header menu items
      dropdownMenu(
        type = "notifications",
        icon = icon("bell"),
        headerText = "Pipeline Status",
        notificationItem(
          text = div(id = "header_status", "Ready"),
          icon = icon("info"),
          status = "info"
        )
      ),
      
      dropdownMenu(
        type = "tasks",
        icon = icon("tasks"),
        headerText = "Quick Actions",
        taskItem(
          text = "Quick Start",
          value = 0,
          color = "green"
        )
      )
    ),
    
    # Dashboard Sidebar
    dashboardSidebar(
      width = 250,
      
      sidebarMenu(
        id = "main_sidebar",
        
        menuItem("Pipeline Workflow", 
                tabName = "pipeline_main", 
                icon = icon("play-circle"),
                selected = TRUE),
        
        menuItem("Project Setup", 
                tabName = "project_setup", 
                icon = icon("folder")),
        
        menuItem("Configuration", 
                tabName = "configuration", 
                icon = icon("cogs")),
        
        menuItem("Monitoring", 
                tabName = "monitoring", 
                icon = icon("chart-line")),
        
        menuItem("Results", 
                tabName = "results", 
                icon = icon("table")),
        
        menuItem("Legacy Interface", 
                tabName = "legacy", 
                icon = icon("history")),
        
        hr(),
        
        # Quick status in sidebar
        div(style = "padding: 10px;",
          h5("Quick Status", style = "color: white;"),
          div(id = "sidebar_status", 
              style = "font-size: 10px; color: #ecf0f1;",
              "Project: Not loaded")
        )
      )
    ),
    
    # Dashboard Body
    dashboardBody(
      
      # Include CSS for custom styling
      tags$head(
        tags$style(HTML("
          .pipeline-step-complete { background-color: #d4edda !important; }
          .pipeline-step-active { background-color: #fff3cd !important; }
          .pipeline-step-pending { background-color: #f8f9fa !important; }
          .info-box { 
            padding: 10px; 
            margin: 5px; 
            border-radius: 5px; 
            text-align: center; 
            color: white; 
          }
          .bg-blue { background-color: #3498db; }
          .bg-green { background-color: #2ecc71; }
          .bg-orange { background-color: #e67e22; }
          .bg-purple { background-color: #9b59b6; }
          .progress { height: 25px; }
          .pipeline-control-panel {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 20px;
          }
          .step-indicator {
            display: inline-block;
            width: 30px;
            height: 30px;
            border-radius: 50%;
            line-height: 30px;
            text-align: center;
            margin-right: 10px;
            font-weight: bold;
          }
          .step-complete { background-color: #28a745; color: white; }
          .step-active { background-color: #ffc107; color: black; }
          .step-pending { background-color: #6c757d; color: white; }
        "))
      ),
      
      # Include shinyjs
      useShinyjs(),
      
      tabItems(
        
        # ===== MAIN PIPELINE TAB =====
        tabItem(
          tabName = "pipeline_main",
          
          # Pipeline Control Panel
          div(class = "pipeline-control-panel",
            fluidRow(
              column(8,
                h2("HiTMaP Pipeline Control Center", icon("rocket"))
              ),
              column(4,
                div(style = "text-align: right; padding-top: 10px;",
                  actionButton("emergency_stop", "Emergency Stop", 
                              icon = icon("stop"), 
                              class = "btn-danger"),
                  actionButton("help_guide", "Help", 
                              icon = icon("question-circle"), 
                              class = "btn-info")
                )
              )
            )
          ),
          
          # Pipeline status and control
          pipeline_control_ui(),
          
          # Main pipeline steps
          pipeline_step_ui(),
          
          # Monitoring section
          pipeline_monitoring_ui()
        ),
        
        # ===== PROJECT SETUP TAB =====
        tabItem(
          tabName = "project_setup",
          
          fluidRow(
            column(12,
              wellPanel(
                h3("Project Setup", icon("folder-open")),
                
                fluidRow(
                  column(6,
                    h4("Data Files"),
                    fileInput("setup_data_files", "Select imzML Files:",
                             multiple = TRUE, accept = c(".imzML")),
                    
                    fileInput("setup_database", "Select Database:",
                             accept = c(".fasta", ".fa")),
                    
                    textInput("setup_project_name", "Project Name:",
                             placeholder = "Enter project name")
                  ),
                  
                  column(6,
                    h4("Project Directory"),
                    shinyDirButton("setup_project_dir", "Choose Project Directory", 
                                  "Select project directory", 
                                  icon = icon("folder")),
                    br(), br(),
                    verbatimTextOutput("setup_project_path"),
                    
                    br(),
                    
                    h4("Configuration"),
                    fileInput("setup_config_file", "Load Configuration:",
                             accept = c(".json", ".rds")),
                    
                    actionButton("create_new_project", "Create New Project",
                                icon = icon("plus"), class = "btn-success")
                  )
                )
              )
            )
          ),
          
          fluidRow(
            column(12,
              wellPanel(
                h3("Recent Projects"),
                DT::dataTableOutput("recent_projects_table")
              )
            )
          )
        ),
        
        # ===== CONFIGURATION TAB =====
        tabItem(
          tabName = "configuration",
          
          fluidRow(
            column(6,
              wellPanel(
                h3("Parameter Presets", icon("magic")),
                
                selectInput("config_preset", "Choose Preset:",
                           choices = list(
                             "Custom" = "custom",
                             "High Precision" = "high_precision", 
                             "High Throughput" = "high_throughput",
                             "Discovery Mode" = "discovery",
                             "Targeted Analysis" = "targeted",
                             "Metabolomics" = "metabolomics"
                           )),
                
                conditionalPanel(
                  condition = "input.config_preset != 'custom'",
                  div(class = "alert alert-info",
                      textOutput("preset_description")
                  )
                ),
                
                actionButton("apply_preset", "Apply Preset", 
                            icon = icon("check"), class = "btn-primary"),
                
                hr(),
                
                h4("Save/Load Configuration"),
                textInput("config_name", "Configuration Name:"),
                actionButton("save_config", "Save Configuration", 
                            icon = icon("save"), class = "btn-success"),
                
                br(), br(),
                
                fileInput("load_config", "Load Configuration File:",
                         accept = c(".json", ".rds"))
              )
            ),
            
            column(6,
              wellPanel(
                h3("Advanced Parameters", icon("sliders-h")),
                
                tabsetPanel(
                  tabPanel("Search",
                    br(),
                    sliderInput("adv_threshold", "Intensity Threshold:",
                               min = 0.00001, max = 0.1, value = 0.001, step = 0.00001),
                    sliderInput("adv_fdr", "FDR Cutoff:",
                               min = 0.001, max = 0.5, value = 0.05, step = 0.001),
                    numericInput("adv_peptide_filter", "Min Peptides/Protein:",
                                value = 2, min = 1, max = 10)
                  ),
                  
                  tabPanel("Processing", 
                    br(),
                    numericInput("adv_threads", "Processing Threads:",
                                value = 4, min = 1, max = 32),
                    selectInput("adv_parallel_method", "Parallel Method:",
                               choices = c("multicore", "multisession", "sequential")),
                    checkboxInput("adv_debug_mode", "Debug Mode", FALSE)
                  ),
                  
                  tabPanel("Export",
                    br(),
                    checkboxGroupInput("adv_export_formats", "Export Formats:",
                                      choices = list("CSV" = "csv", "Excel" = "xlsx", "RDS" = "rds"),
                                      selected = "csv"),
                    checkboxInput("adv_auto_export", "Auto-export results", TRUE),
                    textInput("adv_export_prefix", "Export prefix:", value = "hitmap")
                  )
                )
              )
            )
          )
        ),
        
        # ===== MONITORING TAB =====
        tabItem(
          tabName = "monitoring",
          
          fluidRow(
            column(8,
              wellPanel(
                h3("Real-time Monitoring", icon("chart-line")),
                
                # Performance metrics
                fluidRow(
                  column(3,
                    div(class = "info-box bg-blue",
                        h4(textOutput("monitor_cpu_usage")),
                        p("CPU Usage")
                    )
                  ),
                  column(3,
                    div(class = "info-box bg-green",
                        h4(textOutput("monitor_memory_usage")),
                        p("Memory Usage")
                    )
                  ),
                  column(3,
                    div(class = "info-box bg-orange",
                        h4(textOutput("monitor_elapsed_time")),
                        p("Elapsed Time")
                    )
                  ),
                  column(3,
                    div(class = "info-box bg-purple",
                        h4(textOutput("monitor_est_remaining")),
                        p("Est. Remaining")
                    )
                  )
                ),
                
                br(),
                
                # Progress visualization
                h4("Processing Progress"),
                plotOutput("progress_visualization", height = "300px"),
                
                br(),
                
                # Resource usage plot
                h4("Resource Usage"),
                plotOutput("resource_usage_plot", height = "200px")
              )
            ),
            
            column(4,
              wellPanel(
                h3("Process Control", icon("control")),
                
                div(id = "process_controls",
                  actionButton("pause_processing", "Pause", 
                              icon = icon("pause"), class = "btn-warning"),
                  actionButton("resume_processing", "Resume", 
                              icon = icon("play"), class = "btn-success"),
                  actionButton("abort_processing", "Abort", 
                              icon = icon("stop"), class = "btn-danger")
                ),
                
                br(), br(),
                
                h4("Live Log"),
                div(style = "height: 300px; overflow-y: scroll; background: black; color: green; padding: 10px; font-family: monospace;",
                    verbatimTextOutput("live_log")
                ),
                
                br(),
                
                h4("Error Console"),
                div(style = "height: 150px; overflow-y: scroll; background: #f8d7da; color: #721c24; padding: 10px;",
                    verbatimTextOutput("error_console")
                )
              )
            )
          )
        ),
        
        # ===== RESULTS TAB =====
        tabItem(
          tabName = "results",
          pipeline_results_ui()
        ),
        
        # ===== LEGACY INTERFACE TAB =====
        tabItem(
          tabName = "legacy",
          
          wellPanel(
            h3("Legacy HiTMaP Interface", icon("history")),
            div(class = "alert alert-warning",
                strong("Note: "), "This is the original HiTMaP interface. ",
                "For new analyses, we recommend using the Pipeline Workflow tab."
            )
          ),
          
          # Include original UI components here
          tabsetPanel(
            tabPanel("Project", 
              # Original project UI
              verbatimTextOutput("legacy_project_placeholder")
            ),
            
            tabPanel("Pre-processing",
              # Original preprocessing UI  
              verbatimTextOutput("legacy_preprocess_placeholder")
            ),
            
            tabPanel("Proteomics annotation",
              # Original proteomics UI
              verbatimTextOutput("legacy_proteomics_placeholder")
            )
          )
        )
      )
    )
  )
}

# Helper function to create step indicators
create_step_indicator <- function(step_number, status = "pending") {
  class_name <- paste("step-indicator", paste0("step-", status))
  span(class = class_name, step_number)
}

# Helper function for preset descriptions
get_preset_description <- function(preset) {
  descriptions <- list(
    "high_precision" = "Optimized for high-resolution data with strict FDR controls and comprehensive modifications. Best for detailed protein characterization.",
    "high_throughput" = "Fast processing with relaxed parameters. Suitable for screening large datasets or initial exploration.",
    "discovery" = "Balanced parameters for discovering new proteins and peptides. Moderate stringency with comprehensive search options.",
    "targeted" = "Focused analysis with high precision. Ideal for targeted protein analysis and validation studies.",
    "metabolomics" = "Specialized parameters for small molecule analysis. Lower mass range and modified scoring."
  )
  
  return(descriptions[[preset]] %||% "Custom configuration")
}