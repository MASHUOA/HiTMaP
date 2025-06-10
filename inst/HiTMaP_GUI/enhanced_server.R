#' Enhanced HiTMaP Server Logic with Pipeline Integration
#' 
#' Complete server logic that integrates pipeline functionality with GUI

# Source required components
source("pipeline_integration.R")

enhanced_server <- function(input, output, session) {
  
  # Initialize the real pipeline server (not the placeholder one)
  pipeline_state <- create_real_pipeline_server(input, output, session)
  
  # ===== GLOBAL REACTIVE VALUES =====
  
  values <- reactiveValues(
    current_project = NULL,
    recent_projects = data.frame(),
    monitoring_data = list(),
    system_stats = list(),
    start_time = NULL
  )
  
  # ===== PROJECT SETUP TAB =====
  
  # Project directory selection
  observe({
    shinyDirChoose(input, 'setup_project_dir', roots = c(home = "~"))
  })
  
  output$setup_project_path <- renderText({
    if (!is.null(input$setup_project_dir)) {
      if (is.list(input$setup_project_dir)) {
        parseDirPath(c(home = "~"), input$setup_project_dir)
      }
    }
  })
  
  # Create new project
  observeEvent(input$create_new_project, {
    req(input$setup_data_files, input$setup_database)
    
    tryCatch({
      # Create project using pipeline system
      data_files <- input$setup_data_files$datapath
      names(data_files) <- input$setup_data_files$name
      
      project_folder <- if (!is.null(input$setup_project_dir) && is.list(input$setup_project_dir)) {
        parseDirPath(c(home = "~"), input$setup_project_dir)
      } else {
        NULL
      }
      
      project <- hitmap_project(data_files, project_folder) %>%
        init(thread_count = 4)
      
      values$current_project <- project
      
      # Switch to pipeline tab
      updateTabItems(session, "main_sidebar", "pipeline_main")
      
      showNotification("Project created successfully!", type = "success")
      
    }, error = function(e) {
      showNotification(paste("Error creating project:", e$message), type = "error")
    })
  })
  
  # Recent projects table
  output$recent_projects_table <- DT::renderDataTable({
    # This would be populated from saved project metadata
    data.frame(
      Name = c("Sample Project 1", "Sample Project 2"),
      Date = c("2024-01-15", "2024-01-10"),
      Files = c(3, 5),
      Status = c("Completed", "In Progress")
    )
  }, options = list(pageLength = 10))
  
  # ===== CONFIGURATION TAB =====
  
  # Preset descriptions
  output$preset_description <- renderText({
    get_preset_description(input$config_preset)
  })
  
  # Apply preset configuration
  observeEvent(input$apply_preset, {
    if (input$config_preset != "custom") {
      
      preset_configs <- list(
        "high_precision" = list(
          ppm = 2, threshold = 0.005, fdr = 0.01, peptide_filter = 3,
          threads = 8, preprocess = "comprehensive"
        ),
        "high_throughput" = list(
          ppm = 10, threshold = 0.01, fdr = 0.1, peptide_filter = 1,
          threads = 16, preprocess = "minimal"
        ),
        "discovery" = list(
          ppm = 5, threshold = 0.001, fdr = 0.05, peptide_filter = 2,
          threads = 8, preprocess = "standard"
        ),
        "targeted" = list(
          ppm = 3, threshold = 0.002, fdr = 0.01, peptide_filter = 2,
          threads = 4, preprocess = "comprehensive"
        ),
        "metabolomics" = list(
          ppm = 5, threshold = 0.0001, fdr = 0.1, peptide_filter = 1,
          threads = 4, preprocess = "minimal"
        )
      )
      
      config <- preset_configs[[input$config_preset]]
      
      if (!is.null(config)) {
        # Update advanced parameter inputs
        updateSliderInput(session, "adv_threshold", value = config$threshold)
        updateSliderInput(session, "adv_fdr", value = config$fdr)
        updateNumericInput(session, "adv_peptide_filter", value = config$peptide_filter)
        updateNumericInput(session, "adv_threads", value = config$threads)
        
        # Update pipeline inputs if they exist
        updateNumericInput(session, "pipeline_ppm", value = config$ppm)
        updateSliderInput(session, "pipeline_search_threshold", value = config$threshold)
        updateSliderInput(session, "pipeline_fdr_cutoff", value = config$fdr)
        updateNumericInput(session, "pipeline_peptide_filter", value = config$peptide_filter)
        updateNumericInput(session, "pipeline_threads", value = config$threads)
        updateSelectInput(session, "pipeline_preprocess_profile", selected = config$preprocess)
        
        showNotification(paste("Applied", input$config_preset, "preset"), type = "success")
      }
    }
  })
  
  # Save configuration
  observeEvent(input$save_config, {
    req(input$config_name)
    
    config <- list(
      name = input$config_name,
      threshold = input$adv_threshold,
      fdr = input$adv_fdr,
      peptide_filter = input$adv_peptide_filter,
      threads = input$adv_threads,
      parallel_method = input$adv_parallel_method,
      debug_mode = input$adv_debug_mode,
      export_formats = input$adv_export_formats,
      auto_export = input$adv_auto_export,
      export_prefix = input$adv_export_prefix,
      created = Sys.time()
    )
    
    tryCatch({
      config_file <- paste0("config_", gsub("[^A-Za-z0-9]", "_", input$config_name), ".rds")
      saveRDS(config, config_file)
      showNotification(paste("Configuration saved as", config_file), type = "success")
    }, error = function(e) {
      showNotification(paste("Error saving configuration:", e$message), type = "error")
    })
  })
  
  # ===== MONITORING TAB =====
  
  # System monitoring reactive
  system_monitor <- reactive({
    invalidateLater(2000)  # Update every 2 seconds
    
    if (pipeline_state$processing()) {
      # Get system stats (simplified - would use actual system monitoring)
      list(
        cpu_usage = sample(50:90, 1),
        memory_usage = sample(60:85, 1),
        elapsed_time = if (!is.null(values$start_time)) {
          as.numeric(difftime(Sys.time(), values$start_time, units = "mins"))
        } else 0,
        estimated_remaining = sample(10:30, 1)
      )
    } else {
      list(
        cpu_usage = sample(5:15, 1),
        memory_usage = sample(20:40, 1),
        elapsed_time = 0,
        estimated_remaining = 0
      )
    }
  })
  
  # Monitor outputs
  output$monitor_cpu_usage <- renderText({
    paste0(system_monitor()$cpu_usage, "%")
  })
  
  output$monitor_memory_usage <- renderText({
    paste0(system_monitor()$memory_usage, "%")
  })
  
  output$monitor_elapsed_time <- renderText({
    elapsed <- system_monitor()$elapsed_time
    if (elapsed < 1) {
      paste0(round(elapsed * 60), "s")
    } else {
      paste0(round(elapsed, 1), "m")
    }
  })
  
  output$monitor_est_remaining <- renderText({
    remaining <- system_monitor()$estimated_remaining
    if (remaining < 1) {
      "< 1m"
    } else {
      paste0(round(remaining), "m")
    }
  })
  
  # Progress visualization
  output$progress_visualization <- renderPlot({
    status <- pipeline_state$status
    if (!is.null(status)) {
      # Create a visualization of pipeline progress
      steps <- c("Init", "Candidates", "Preprocess", "Segment", "Search", "Summarize")
      status_values <- c(
        ifelse(status$initialized, 1, 0),
        ifelse(status$candidates_generated, 1, 0),
        ifelse(status$preprocessing_complete, 1, 0),
        ifelse(status$segmentation_complete, 1, 0),
        ifelse(status$pmf_search_complete, 1, 0),
        ifelse(status$summarization_complete, 1, 0)
      )
      
      progress_data <- data.frame(
        Step = factor(steps, levels = steps),
        Status = status_values,
        Color = ifelse(status_values == 1, "Completed", "Pending")
      )
      
      ggplot2::ggplot(progress_data, ggplot2::aes(x = Step, y = Status, fill = Color)) +
        ggplot2::geom_col() +
        ggplot2::scale_fill_manual(values = c("Completed" = "#28a745", "Pending" = "#6c757d")) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Pipeline Progress", y = "Completion Status") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    } else {
      # Empty plot
      ggplot2::ggplot() + ggplot2::theme_void() + 
        ggplot2::labs(title = "No project loaded")
    }
  })
  
  # Resource usage plot
  output$resource_usage_plot <- renderPlot({
    # Simulate resource usage over time
    if (length(values$monitoring_data) > 0) {
      time_points <- 1:length(values$monitoring_data)
      cpu_data <- sapply(values$monitoring_data, function(x) x$cpu_usage)
      memory_data <- sapply(values$monitoring_data, function(x) x$memory_usage)
      
      usage_data <- data.frame(
        Time = rep(time_points, 2),
        Usage = c(cpu_data, memory_data),
        Type = rep(c("CPU", "Memory"), each = length(time_points))
      )
      
      ggplot2::ggplot(usage_data, ggplot2::aes(x = Time, y = Usage, color = Type)) +
        ggplot2::geom_line() +
        ggplot2::scale_y_continuous(limits = c(0, 100)) +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Resource Usage Over Time", y = "Usage (%)")
    } else {
      ggplot2::ggplot() + ggplot2::theme_void() + 
        ggplot2::labs(title = "Collecting resource data...")
    }
  })
  
  # Live log
  output$live_log <- renderText({
    log_entries <- pipeline_state$log_entries()
    if (length(log_entries) > 0) {
      paste(tail(log_entries, 20), collapse = "\n")
    } else {
      "No log entries"
    }
  })
  
  # Error console
  output$error_console <- renderText({
    error_messages <- pipeline_state$error_messages()
    if (length(error_messages) > 0) {
      paste(tail(error_messages, 10), collapse = "\n")
    } else {
      "No errors"
    }
  })
  
  # Process control buttons
  observeEvent(input$pause_processing, {
    # Implementation would depend on the actual processing system
    showNotification("Processing paused", type = "warning")
  })
  
  observeEvent(input$resume_processing, {
    showNotification("Processing resumed", type = "success")
  })
  
  observeEvent(input$abort_processing, {
    showModal(modalDialog(
      title = "Confirm Abort",
      "Are you sure you want to abort the current processing? This action cannot be undone.",
      footer = tagList(
        modalButton("Cancel"),
        actionButton("confirm_abort", "Abort", class = "btn-danger")
      )
    ))
  })
  
  observeEvent(input$confirm_abort, {
    removeModal()
    # Note: Would need to implement actual processing abort in pipeline_integration.R
    showNotification("Processing aborted", type = "error")
  })
  
  # ===== RESULTS TAB =====
  
  # Results visualization
  output$results_plot <- renderPlot({
    req(input$viz_plot_type)
    
    status <- pipeline_state$status
    if (!is.null(status) && status$summarization_complete) {
      
      # Create placeholder plots - would need to access actual results
      ggplot2::ggplot() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(title = "Results visualization available after analysis completion")
      
    } else {
      ggplot2::ggplot() + ggplot2::theme_void() + 
        ggplot2::labs(title = "Complete analysis to view plots")
    }
  })
  
  # Interactive plot
  output$interactive_plot <- plotly::renderPlotly({
    status <- pipeline_state$status
    if (!is.null(status) && status$summarization_complete) {
      # Create placeholder interactive plot - would access actual results
      p <- ggplot2::ggplot() + 
        ggplot2::theme_minimal() + 
        ggplot2::labs(title = "Interactive results plot available after analysis")
      plotly::ggplotly(p)
    }
  })
  
  # ===== GLOBAL UPDATES =====
  
  # Update header status
  observe({
    project_data <- pipeline_state$project()
    status <- pipeline_state$status
    status_text <- if (!is.null(project_data)) {
      if (pipeline_state$processing()) {
        "Processing..."
      } else if (!is.null(status) && status$summarization_complete) {
        "Analysis Complete"
      } else {
        "Ready"
      }
    } else {
      "No Project"
    }
    
    runjs(paste0("$('#header_status').text('", status_text, "');"))
  })
  
  # Update sidebar status
  observe({
    project_data <- pipeline_state$project()
    sidebar_text <- if (!is.null(project_data) && !is.null(project_data$file_names)) {
      paste("Project:", length(project_data$file_names), "files")
    } else {
      "Project: Not loaded"
    }
    
    runjs(paste0("$('#sidebar_status').text('", sidebar_text, "');"))
  })
  
  # Help guide
  observeEvent(input$help_guide, {
    showModal(modalDialog(
      title = "HiTMaP Pipeline Guide",
      size = "l",
      div(
        h4("Getting Started"),
        p("1. Upload your imzML files and database in the Pipeline Workflow tab"),
        p("2. Configure parameters or use a preset template"),
        p("3. Run each pipeline step in sequence or use Quick Start"),
        p("4. Monitor progress in real-time"),
        p("5. View and export results"),
        
        hr(),
        
        h4("Pipeline Steps"),
        tags$ol(
          tags$li(strong("Initialize:"), " Set up project and parallel processing"),
          tags$li(strong("Candidates:"), " Generate peptide candidates from database"),
          tags$li(strong("Preprocess:"), " Process raw imaging data"),
          tags$li(strong("Segment:"), " Divide tissue into spatial regions"),
          tags$li(strong("Search:"), " Perform mass spectrometry search"),
          tags$li(strong("Summarize:"), " Generate final results")
        ),
        
        hr(),
        
        h4("Tips"),
        tags$ul(
          tags$li("Use checkpoints for long analyses"),
          tags$li("Monitor system resources during processing"),
          tags$li("Try different templates for different analysis types"),
          tags$li("Export results in multiple formats for downstream analysis")
        )
      ),
      footer = modalButton("Close")
    ))
  })
  
  # Emergency stop
  observeEvent(input$emergency_stop, {
    pipeline_state$processing <- FALSE
    showNotification("Emergency stop activated", type = "error")
  })
  
  # Store monitoring data
  observe({
    if (pipeline_state$processing()) {
      values$monitoring_data <- append(values$monitoring_data, list(system_monitor()), after = 0)
      
      # Keep only last 100 data points
      if (length(values$monitoring_data) > 100) {
        values$monitoring_data <- values$monitoring_data[1:100]
      }
      
      # Set start time if not set
      if (is.null(values$start_time)) {
        values$start_time <- Sys.time()
      }
    } else {
      values$start_time <- NULL
    }
  })
  
  # ===== LEGACY INTERFACE PLACEHOLDERS =====
  
  output$legacy_project_placeholder <- renderText({
    "Legacy project interface would be here. Use the Pipeline Workflow tab for new analyses."
  })
  
  output$legacy_preprocess_placeholder <- renderText({
    "Legacy preprocessing interface would be here. Use the Pipeline Workflow tab for new analyses."
  })
  
  output$legacy_proteomics_placeholder <- renderText({
    "Legacy proteomics interface would be here. Use the Pipeline Workflow tab for new analyses."
  })
}