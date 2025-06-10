#' HiTMaP Pipeline GUI Components
#' 
#' Enhanced UI components that integrate with the new pipeline system

#' Pipeline Control Panel UI
pipeline_control_ui <- function() {
  fluidRow(
    column(12,
      wellPanel(
        h3("HiTMaP Pipeline Control", icon("play-circle")),
        
        fluidRow(
          column(3,
            h4("Pipeline Status"),
            div(id = "pipeline_status_display", 
                style = "font-family: monospace; background: #f8f9fa; padding: 10px; border-radius: 5px;",
                verbatimTextOutput("pipeline_status_text")
            ),
            br(),
            actionButton("refresh_status", "Refresh Status", 
                        icon = icon("refresh"), class = "btn-info btn-sm")
          ),
          
          column(3,
            h4("Quick Actions"),
            actionButton("pipeline_quick_start", "Quick Start", 
                        icon = icon("rocket"), class = "btn-success"),
            br(), br(),
            actionButton("pipeline_step_by_step", "Step-by-Step", 
                        icon = icon("list"), class = "btn-primary"),
            br(), br(),
            actionButton("pipeline_resume", "Resume Project", 
                        icon = icon("play"), class = "btn-warning")
          ),
          
          column(3,
            h4("Pipeline Templates"),
            selectInput("pipeline_template", "Choose Template:",
                       choices = list(
                         "Custom" = "custom",
                         "High Precision" = "high_precision",
                         "High Throughput" = "high_throughput",
                         "Discovery Mode" = "discovery",
                         "Targeted Analysis" = "targeted"
                       )),
            actionButton("apply_template", "Apply Template", 
                        icon = icon("magic"), class = "btn-info")
          ),
          
          column(3,
            h4("Project Management"),
            actionButton("save_checkpoint", "Save Checkpoint", 
                        icon = icon("save"), class = "btn-secondary"),
            br(), br(),
            actionButton("load_checkpoint", "Load Checkpoint", 
                        icon = icon("folder-open"), class = "btn-secondary"),
            br(), br(),
            downloadButton("export_pipeline_results", "Export Results", 
                          class = "btn-success")
          )
        )
      )
    )
  )
}

#' Pipeline Step Configuration UI
pipeline_step_ui <- function() {
  tabsetPanel(
    id = "pipeline_steps_tabs",
    
    # Step 1: Project Initialization
    tabPanel("1. Initialize", 
      icon = icon("cog"),
      wellPanel(
        h4("Project Initialization"),
        
        fluidRow(
          column(6,
            fileInput("pipeline_data_files", "Select Data Files (.imzML)",
                     multiple = TRUE, accept = c(".imzML")),
            
            textInput("pipeline_project_folder", "Project Folder:", 
                     placeholder = "Leave empty for auto-detection"),
            
            numericInput("pipeline_threads", "Number of Threads:", 
                        value = 4, min = 1, max = 32, step = 1)
          ),
          
          column(6,
            div(style = "margin-top: 25px;",
              actionButton("run_init", "Initialize Project", 
                          icon = icon("play"), class = "btn-primary"),
              br(), br(),
              
              conditionalPanel(
                condition = "output.pipeline_initialized == true",
                div(class = "alert alert-success",
                    icon("check"), " Project initialized successfully!")
              )
            )
          )
        )
      )
    ),
    
    # Step 2: Candidate Generation
    tabPanel("2. Candidates", 
      icon = icon("database"),
      wellPanel(
        h4("Candidate Generation"),
        
        fluidRow(
          column(4,
            fileInput("pipeline_database", "Database File (.fasta)",
                     accept = c(".fasta", ".fa")),
            
            selectInput("pipeline_digestion", "Digestion Site:",
                       choices = c("trypsin", "chymotrypsin", "pepsin", "custom"),
                       selected = "trypsin"),
            
            conditionalPanel(
              condition = "input.pipeline_digestion == 'custom'",
              textInput("pipeline_custom_digestion", "Custom Digestion Pattern:")
            ),
            
            sliderInput("pipeline_missed_cleavages", "Missed Cleavages:",
                       min = 0, max = 5, value = c(0, 1))
          ),
          
          column(4,
            selectizeInput("pipeline_adducts", "Adducts:",
                          choices = c("M+H", "M+Na", "M+K", "M+NH4"),
                          selected = "M+H", multiple = TRUE),
            
            selectInput("pipeline_modifications", "Modification Profile:",
                       choices = list(
                         "None" = "none",
                         "Basic" = "basic", 
                         "Comprehensive" = "comprehensive",
                         "Custom" = "custom"
                       )),
            
            conditionalPanel(
              condition = "input.pipeline_modifications == 'custom'",
              textAreaInput("pipeline_custom_mods", "Custom Modifications (JSON):",
                           rows = 3)
            )
          ),
          
          column(4,
            checkboxInput("pipeline_decoy_search", "Enable Decoy Search", TRUE),
            
            conditionalPanel(
              condition = "input.pipeline_decoy_search",
              selectInput("pipeline_decoy_mode", "Decoy Mode:",
                         choices = c("isotope", "element", "adduct"))
            ),
            
            sliderInput("pipeline_mz_range", "m/z Range:",
                       min = 0, max = 5000, value = c(700, 4000)),
            
            div(style = "margin-top: 20px;",
              actionButton("run_candidates", "Generate Candidates", 
                          icon = icon("play"), class = "btn-primary")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.candidates_generated == true",
          br(),
          div(class = "alert alert-info",
              h5("Candidates Summary:"),
              textOutput("candidates_summary")
          )
        )
      )
    ),
    
    # Step 3: Preprocessing
    tabPanel("3. Preprocess", 
      icon = icon("cogs"),
      wellPanel(
        h4("Data Preprocessing"),
        
        fluidRow(
          column(4,
            selectInput("pipeline_preprocess_profile", "Preprocessing Profile:",
                       choices = list(
                         "Minimal" = "minimal",
                         "Standard" = "standard",
                         "Comprehensive" = "comprehensive",
                         "Custom" = "custom"
                       )),
            
            numericInput("pipeline_ppm", "Mass Tolerance (ppm):",
                        value = 5, min = 0.1, max = 50, step = 0.1),
            
            sliderInput("pipeline_preprocess_mz_range", "Processing m/z Range:",
                       min = 0, max = 5000, value = c(700, 4000))
          ),
          
          column(4,
            conditionalPanel(
              condition = "input.pipeline_preprocess_profile == 'custom'",
              h5("Custom Preprocessing Parameters:"),
              
              selectInput("pipeline_smooth_method", "Smoothing:",
                         choices = list("Disable" = "disable", "Gaussian" = "gaussian", 
                                      "Moving Average" = "ma", "Savitzky-Golay" = "sgolay")),
              
              selectInput("pipeline_baseline_method", "Baseline Reduction:",
                         choices = list("Disable" = "disable", "Local Minimum" = "locmin", 
                                      "Median" = "median")),
              
              selectInput("pipeline_peak_method", "Peak Picking:",
                         choices = list("Adaptive" = "adaptive", "Simple" = "simple", 
                                      "MAD" = "mad", "Disable" = "disable"))
            )
          ),
          
          column(4,
            checkboxInput("pipeline_force_preprocess", "Force Reprocessing", FALSE),
            checkboxInput("pipeline_use_rds", "Use Preprocessed RDS", TRUE),
            
            div(style = "margin-top: 40px;",
              actionButton("run_preprocess", "Start Preprocessing", 
                          icon = icon("play"), class = "btn-primary")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.preprocessing_complete == true",
          br(),
          div(class = "alert alert-success",
              icon("check"), " Preprocessing completed for all files"
          )
        )
      )
    ),
    
    # Step 4: Segmentation
    tabPanel("4. Segment", 
      icon = icon("th"),
      wellPanel(
        h4("Spatial Segmentation"),
        
        fluidRow(
          column(4,
            numericInput("pipeline_segment_count", "Number of Segments:",
                        value = 4, min = 2, max = 20, step = 1),
            
            selectInput("pipeline_segment_method", "Segmentation Method:",
                       choices = list(
                         "Spatial k-Means" = "spatialKMeans",
                         "Spatial Shrunken Centroids" = "spatialShrunkenCentroids",
                         "Manual (File-based)" = "def_file",
                         "None" = "none"
                       )),
            
            conditionalPanel(
              condition = "input.pipeline_segment_method == 'def_file'",
              fileInput("pipeline_segment_def", "Segmentation Definition File:",
                       accept = c(".csv", ".txt"))
            )
          ),
          
          column(4,
            sliderInput("pipeline_variance_coverage", "Variance Coverage:",
                       min = 0.1, max = 1.0, value = 0.8, step = 0.05),
            
            numericInput("pipeline_smooth_range", "Smoothing Range:",
                        value = 1, min = 0, max = 5, step = 1),
            
            selectInput("pipeline_segment_ncomp", "Components:",
                       choices = list("Auto-detect" = "auto-detect", "3" = 3, "6" = 6, "9" = 9))
          ),
          
          column(4,
            div(style = "margin-top: 40px;",
              actionButton("run_segmentation", "Start Segmentation", 
                          icon = icon("play"), class = "btn-primary")
            ),
            
            br(), br(),
            
            conditionalPanel(
              condition = "output.segmentation_complete == true",
              actionButton("preview_segmentation", "Preview Results", 
                          icon = icon("eye"), class = "btn-info")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.segmentation_complete == true",
          br(),
          div(class = "alert alert-success",
              textOutput("segmentation_summary")
          )
        )
      )
    ),
    
    # Step 5: PMF Search
    tabPanel("5. Search", 
      icon = icon("search"),
      wellPanel(
        h4("PMF Search"),
        
        fluidRow(
          column(4,
            sliderInput("pipeline_search_threshold", "Intensity Threshold:",
                       min = 0.0001, max = 0.1, value = 0.001, step = 0.0001),
            
            numericInput("pipeline_search_ppm", "Search Tolerance (ppm):",
                        value = 5, min = 0.1, max = 50, step = 0.1),
            
            selectInput("pipeline_score_method", "Scoring Method:",
                       choices = c("SQRTP", "Cosine", "Correlation"))
          ),
          
          column(4,
            sliderInput("pipeline_fdr_cutoff", "FDR Cutoff:",
                       min = 0.001, max = 0.2, value = 0.05, step = 0.001),
            
            numericInput("pipeline_peptide_filter", "Min. Peptides per Protein:",
                        value = 2, min = 1, max = 10, step = 1),
            
            checkboxInput("pipeline_plot_scores", "Generate Score Plots", FALSE)
          ),
          
          column(4,
            h5("Region Selection:"),
            radioButtons("pipeline_region_selection", "",
                        choices = list(
                          "All Regions" = "all",
                          "Specific Regions" = "specific",
                          "Filter by Size" = "filter"
                        ), selected = "all"),
            
            conditionalPanel(
              condition = "input.pipeline_region_selection == 'specific'",
              textInput("pipeline_specific_regions", "Region Names (comma-separated):")
            ),
            
            div(style = "margin-top: 20px;",
              actionButton("run_pmf_search", "Start PMF Search", 
                          icon = icon("play"), class = "btn-primary")
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.pmf_search_complete == true",
          br(),
          div(class = "alert alert-success",
              h5("PMF Search Results:"),
              textOutput("pmf_search_summary")
          )
        )
      )
    ),
    
    # Step 6: Summarization
    tabPanel("6. Summarize", 
      icon = icon("chart-bar"),
      wellPanel(
        h4("Generate Summaries"),
        
        fluidRow(
          column(6,
            h5("Summary Options:"),
            checkboxInput("pipeline_protein_summary", "Protein Summary", TRUE),
            checkboxInput("pipeline_peptide_summary", "Peptide Summary", TRUE),
            checkboxInput("pipeline_region_summary", "Region Feature Summary", FALSE)
          ),
          
          column(6,
            div(style = "margin-top: 20px;",
              actionButton("run_summarize", "Generate Summaries", 
                          icon = icon("play"), class = "btn-primary")
            ),
            
            br(), br(),
            
            conditionalPanel(
              condition = "output.summarization_complete == true",
              div(class = "alert alert-success",
                  icon("check"), " Summaries generated successfully!"
              )
            )
          )
        ),
        
        conditionalPanel(
          condition = "output.summarization_complete == true",
          br(),
          hr(),
          h5("Final Results Summary:"),
          div(id = "final_results_summary",
              verbatimTextOutput("final_results_text")
          )
        )
      )
    )
  )
}

#' Pipeline Monitoring UI
pipeline_monitoring_ui <- function() {
  fluidRow(
    column(6,
      wellPanel(
        h4("Pipeline Progress"),
        
        div(id = "progress_container",
            htmlOutput("pipeline_progress_bar")
        ),
        
        br(),
        
        h5("Processing Log:"),
        div(style = "height: 300px; overflow-y: scroll; background: #f8f9fa; padding: 10px; border-radius: 5px;",
            verbatimTextOutput("pipeline_log")
        )
      )
    ),
    
    column(6,
      wellPanel(
        h4("Performance Monitoring"),
        
        conditionalPanel(
          condition = "output.benchmarks_available == true",
          plotOutput("benchmark_plot", height = "200px"),
          br()
        ),
        
        h5("Current Status:"),
        tableOutput("pipeline_status_table"),
        
        br(),
        
        h5("Quick Actions:"),
        fluidRow(
          column(6,
            actionButton("pipeline_pause", "Pause", 
                        icon = icon("pause"), class = "btn-warning btn-sm"),
            actionButton("pipeline_stop", "Stop", 
                        icon = icon("stop"), class = "btn-danger btn-sm")
          ),
          column(6,
            actionButton("pipeline_debug", "Debug Mode", 
                        icon = icon("bug"), class = "btn-info btn-sm"),
            actionButton("view_detailed_status", "Detailed Status", 
                        icon = icon("info"), class = "btn-secondary btn-sm")
          )
        )
      )
    )
  )
}

#' Results Visualization UI
pipeline_results_ui <- function() {
  tabsetPanel(
    id = "results_tabs",
    
    tabPanel("Summary", 
      icon = icon("table"),
      fluidRow(
        column(12,
          h4("Results Overview"),
          
          conditionalPanel(
            condition = "output.has_results == true",
            
            fluidRow(
              column(3,
                div(class = "info-box bg-blue",
                    h3(textOutput("total_peptides_count")),
                    p("Total Peptides")
                )
              ),
              column(3,
                div(class = "info-box bg-green", 
                    h3(textOutput("total_proteins_count")),
                    p("Total Proteins")
                )
              ),
              column(3,
                div(class = "info-box bg-orange",
                    h3(textOutput("total_regions_count")),
                    p("Total Regions")
                )
              ),
              column(3,
                div(class = "info-box bg-purple",
                    h3(textOutput("total_files_count")),
                    p("Files Processed")
                )
              )
            ),
            
            br(),
            
            tabsetPanel(
              tabPanel("Peptides", DT::dataTableOutput("peptides_table")),
              tabPanel("Proteins", DT::dataTableOutput("proteins_table")),
              tabPanel("By Region", DT::dataTableOutput("regions_table"))
            )
          ),
          
          conditionalPanel(
            condition = "output.has_results == false",
            div(class = "alert alert-info",
                h4("No Results Available"),
                p("Complete the pipeline to view results here.")
            )
          )
        )
      )
    ),
    
    tabPanel("Visualizations", 
      icon = icon("chart-line"),
      fluidRow(
        column(6,
          wellPanel(
            h5("Result Distribution"),
            selectInput("viz_plot_type", "Plot Type:",
                       choices = list(
                         "Peptides by File" = "peptides_file",
                         "Peptides by Region" = "peptides_region", 
                         "Proteins by File" = "proteins_file",
                         "Mass Distribution" = "mass_dist",
                         "Score Distribution" = "score_dist"
                       )),
            plotOutput("results_plot")
          )
        ),
        
        column(6,
          wellPanel(
            h5("Interactive Exploration"),
            conditionalPanel(
              condition = "output.has_results == true",
              plotlyOutput("interactive_plot")
            ),
            conditionalPanel(
              condition = "output.has_results == false",
              div(class = "alert alert-info", "Complete analysis to view interactive plots")
            )
          )
        )
      )
    ),
    
    tabPanel("Export", 
      icon = icon("download"),
      wellPanel(
        h4("Export Results"),
        
        fluidRow(
          column(4,
            h5("Export Format:"),
            radioButtons("export_format", "",
                        choices = list("CSV" = "csv", "Excel" = "excel", "RDS" = "rds"),
                        selected = "csv"),
            
            h5("Export Types:"),
            checkboxGroupInput("export_types", "",
                              choices = list(
                                "Summaries" = "summaries",
                                "Candidates" = "candidates", 
                                "Raw Results" = "raw_results",
                                "Project Object" = "project"
                              ),
                              selected = c("summaries"))
          ),
          
          column(4,
            h5("Export Options:"),
            textInput("export_prefix", "File Prefix:", value = "hitmap_results"),
            checkboxInput("export_timestamp", "Add Timestamp", TRUE),
            
            br(),
            
            downloadButton("download_results", "Download Results", 
                          class = "btn-success", icon = icon("download"))
          ),
          
          column(4,
            h5("Quick Exports:"),
            downloadButton("quick_export_peptides", "Peptides (CSV)", 
                          class = "btn-info btn-sm"),
            br(), br(),
            downloadButton("quick_export_proteins", "Proteins (CSV)", 
                          class = "btn-info btn-sm"),
            br(), br(),
            downloadButton("quick_export_project", "Full Project (RDS)", 
                          class = "btn-warning btn-sm")
          )
        )
      )
    )
  )
}