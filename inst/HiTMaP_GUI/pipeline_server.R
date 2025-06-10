#' HiTMaP Pipeline Server Logic
#' 
#' Server-side logic for the pipeline-based GUI

# Initialize pipeline-related reactive values
init_pipeline_reactives <- function() {
  return(reactiveValues(
    project = NULL,
    current_step = 0,
    status = list(
      initialized = FALSE,
      candidates_generated = FALSE,
      preprocessing_complete = FALSE,
      segmentation_complete = FALSE,
      pmf_search_complete = FALSE,
      summarization_complete = FALSE
    ),
    processing = FALSE,
    error_messages = character(),
    log_entries = character(),
    benchmarks = list()
  ))
}

#' Pipeline Server Logic
pipeline_server <- function(input, output, session) {
  
  # Initialize reactive values
  pipeline_state <- init_pipeline_reactives()
  
  # Helper function to add log entry
  add_log_entry <- function(message) {
    timestamp <- format(Sys.time(), "%H:%M:%S")
    entry <- paste0("[", timestamp, "] ", message)
    pipeline_state$log_entries <- c(entry, pipeline_state$log_entries)
    
    # Keep only last 100 entries
    if (length(pipeline_state$log_entries) > 100) {
      pipeline_state$log_entries <- pipeline_state$log_entries[1:100]
    }
  }
  
  # Helper function to update status
  update_pipeline_status <- function() {
    if (!is.null(pipeline_state$project)) {
      status <- hitmap_project_status(pipeline_state$project)
      pipeline_state$status <- list(
        initialized = status$processing_progress$initialized,
        candidates_generated = status$processing_progress$candidates_generated,
        preprocessing_complete = status$processing_progress$preprocessing_complete,
        segmentation_complete = status$processing_progress$segmentation_complete,
        pmf_search_complete = status$processing_progress$pmf_search_complete,
        summarization_complete = length(status$results_summary$summaries_available) > 0
      )
    }
  }
  
  # Helper function to run pipeline step safely
  run_pipeline_step <- function(step_name, step_function) {
    tryCatch({
      pipeline_state$processing <- TRUE
      add_log_entry(paste("Starting", step_name, "..."))
      
      start_time <- Sys.time()
      result <- step_function()
      end_time <- Sys.time()
      
      # Store benchmark
      pipeline_state$benchmarks[[step_name]] <- as.numeric(difftime(end_time, start_time, units = "secs"))
      
      pipeline_state$project <- result
      update_pipeline_status()
      
      add_log_entry(paste(step_name, "completed successfully"))
      pipeline_state$processing <- FALSE
      
      return(TRUE)
      
    }, error = function(e) {
      pipeline_state$processing <- FALSE
      error_msg <- paste("Error in", step_name, ":", e$message)
      pipeline_state$error_messages <- c(error_msg, pipeline_state$error_messages)
      add_log_entry(error_msg)
      
      showNotification(error_msg, type = "error", duration = 10)
      return(FALSE)
    })
  }
  
  # ===== STEP 1: INITIALIZATION =====
  
  observeEvent(input$run_init, {
    req(input$pipeline_data_files)
    
    run_pipeline_step("Project Initialization", function() {
      data_files <- input$pipeline_data_files$datapath
      names(data_files) <- input$pipeline_data_files$name
      
      project_folder <- if (input$pipeline_project_folder != "") {
        input$pipeline_project_folder
      } else {
        NULL
      }
      
      project <- hitmap_project(data_files, project_folder) %>%
        init(thread_count = input$pipeline_threads)
      
      return(project)
    })
  })
  
  # ===== STEP 2: CANDIDATE GENERATION =====
  
  observeEvent(input$run_candidates, {
    req(pipeline_state$status$initialized)
    
    run_pipeline_step("Candidate Generation", function() {
      
      # Build modifications list
      if (input$pipeline_modifications == "custom") {
        modifications <- tryCatch({
          jsonlite::fromJSON(input$pipeline_custom_mods)
        }, error = function(e) {
          hitmap_create_modification_config("basic")
        })
      } else {
        modifications <- hitmap_create_modification_config(input$pipeline_modifications)
      }
      
      # Handle custom digestion
      digestion_site <- if (input$pipeline_digestion == "custom") {
        input$pipeline_custom_digestion
      } else {
        input$pipeline_digestion
      }
      
      project <- pipeline_state$project %>%
        generate_candidates(
          database = input$pipeline_database$datapath,
          digestion_site = digestion_site,
          missed_cleavages = input$pipeline_missed_cleavages[1]:input$pipeline_missed_cleavages[2],
          adducts = input$pipeline_adducts,
          modifications = modifications,
          decoy_search = input$pipeline_decoy_search,
          decoy_mode = if (input$pipeline_decoy_search) input$pipeline_decoy_mode else "isotope",
          mz_range = input$pipeline_mz_range
        )
      
      return(project)
    })
  })
  
  # ===== STEP 3: PREPROCESSING =====
  
  observeEvent(input$run_preprocess, {
    req(pipeline_state$status$candidates_generated)
    
    run_pipeline_step("Data Preprocessing", function() {
      
      # Build preprocessing parameters
      if (input$pipeline_preprocess_profile == "custom") {
        preprocess_params <- hitmap_create_preprocess_config("standard", input$pipeline_ppm)
        
        # Override with custom parameters
        preprocess_params$smooth_signal <- list(method = input$pipeline_smooth_method)
        preprocess_params$reduce_baseline <- list(method = input$pipeline_baseline_method)
        preprocess_params$peak_pick <- list(method = input$pipeline_peak_method)
      } else {
        preprocess_params <- hitmap_create_preprocess_config(input$pipeline_preprocess_profile, input$pipeline_ppm)
      }
      
      # Override force settings
      preprocess_params$force_preprocess <- input$pipeline_force_preprocess
      preprocess_params$use_preprocessrds <- input$pipeline_use_rds
      
      project <- pipeline_state$project %>%
        preprocess(
          ppm = input$pipeline_ppm,
          mz_range = input$pipeline_preprocess_mz_range,
          preprocess_params = preprocess_params
        )
      
      return(project)
    })
  })
  
  # ===== STEP 4: SEGMENTATION =====
  
  observeEvent(input$run_segmentation, {
    req(pipeline_state$status$preprocessing_complete)
    
    run_pipeline_step("Spatial Segmentation", function() {
      
      segmentation_params <- list(
        variance_coverage = input$pipeline_variance_coverage,
        smooth_range = input$pipeline_smooth_range,
        segmentation_ncomp = input$pipeline_segment_ncomp
      )
      
      # Handle segmentation definition file
      if (input$pipeline_segment_method == "def_file" && !is.null(input$pipeline_segment_def)) {
        segmentation_params$segmentation_def <- input$pipeline_segment_def$datapath
      }
      
      project <- pipeline_state$project %>%
        segment(
          segment_count = input$pipeline_segment_count,
          segmentation_method = input$pipeline_segment_method,
          segmentation_params = segmentation_params
        )
      
      return(project)
    })
  })
  
  # ===== STEP 5: PMF SEARCH =====
  
  observeEvent(input$run_pmf_search, {
    req(pipeline_state$status$segmentation_complete)
    
    run_pipeline_step("PMF Search", function() {
      
      # Handle region selection
      region_filter <- switch(input$pipeline_region_selection,
        "all" = NULL,
        "specific" = {
          if (input$pipeline_specific_regions != "") {
            strsplit(input$pipeline_specific_regions, ",")[[1]]
          } else {
            NULL
          }
        },
        "filter" = function(regions) regions[1:min(3, length(regions))]  # Example filter
      )
      
      search_params <- list(
        plot_scores = input$pipeline_plot_scores
      )
      
      project <- pipeline_state$project %>%
        search_pmf(
          region_filter = region_filter,
          threshold = input$pipeline_search_threshold,
          ppm = input$pipeline_search_ppm,
          score_method = input$pipeline_score_method,
          fdr_cutoff = input$pipeline_fdr_cutoff,
          peptide_filter = input$pipeline_peptide_filter,
          search_params = search_params
        )
      
      return(project)
    })
  })
  
  # ===== STEP 6: SUMMARIZATION =====
  
  observeEvent(input$run_summarize, {
    req(pipeline_state$status$pmf_search_complete)
    
    run_pipeline_step("Generate Summaries", function() {
      project <- pipeline_state$project %>%
        summarize(
          protein_summary = input$pipeline_protein_summary,
          peptide_summary = input$pipeline_peptide_summary,
          region_summary = input$pipeline_region_summary
        )
      
      return(project)
    })
  })
  
  # ===== QUICK START =====
  
  observeEvent(input$pipeline_quick_start, {
    req(input$pipeline_data_files, input$pipeline_database)
    
    run_pipeline_step("Quick Start Pipeline", function() {
      data_files <- input$pipeline_data_files$datapath
      names(data_files) <- input$pipeline_data_files$name
      
      project <- hitmap_quick_start(
        data_files = data_files,
        database = input$pipeline_database$datapath,
        thread_count = input$pipeline_threads %||% 4
      )
      
      return(project)
    })
  })
  
  # ===== TEMPLATE APPLICATION =====
  
  observeEvent(input$apply_template, {
    req(pipeline_state$status$initialized)
    
    template_config <- switch(input$pipeline_template,
      "high_precision" = list(
        ppm = 2,
        threshold = 0.005,
        segment_count = 8,
        preprocess_profile = "comprehensive",
        modifications = "comprehensive",
        fdr_cutoff = 0.01,
        peptide_filter = 3
      ),
      "high_throughput" = list(
        ppm = 10,
        threshold = 0.01,
        segment_count = 4,
        preprocess_profile = "minimal",
        modifications = "basic",
        fdr_cutoff = 0.1,
        peptide_filter = 1
      ),
      "discovery" = list(
        ppm = 8,
        threshold = 0.001,
        segment_count = 6,
        preprocess_profile = "standard",
        modifications = "comprehensive",
        fdr_cutoff = 0.1,
        peptide_filter = 1
      ),
      "targeted" = list(
        ppm = 3,
        threshold = 0.002,
        segment_count = 10,
        preprocess_profile = "comprehensive",
        modifications = "basic",
        fdr_cutoff = 0.05,
        peptide_filter = 2
      ),
      list()  # custom
    )
    
    if (length(template_config) > 0) {
      # Update UI elements with template values
      updateNumericInput(session, "pipeline_ppm", value = template_config$ppm)
      updateSliderInput(session, "pipeline_search_threshold", value = template_config$threshold)
      updateNumericInput(session, "pipeline_segment_count", value = template_config$segment_count)
      updateSelectInput(session, "pipeline_preprocess_profile", selected = template_config$preprocess_profile)
      updateSelectInput(session, "pipeline_modifications", selected = template_config$modifications)
      updateSliderInput(session, "pipeline_fdr_cutoff", value = template_config$fdr_cutoff)
      updateNumericInput(session, "pipeline_peptide_filter", value = template_config$peptide_filter)
      
      showNotification(paste("Applied", input$pipeline_template, "template"), type = "success")
    }
  })
  
  # ===== OUTPUTS =====
  
  # Pipeline status outputs
  output$pipeline_initialized <- reactive({ pipeline_state$status$initialized })
  output$candidates_generated <- reactive({ pipeline_state$status$candidates_generated })
  output$preprocessing_complete <- reactive({ pipeline_state$status$preprocessing_complete })
  output$segmentation_complete <- reactive({ pipeline_state$status$segmentation_complete })
  output$pmf_search_complete <- reactive({ pipeline_state$status$pmf_search_complete })
  output$summarization_complete <- reactive({ pipeline_state$status$summarization_complete })
  
  outputOptions(output, "pipeline_initialized", suspendWhenHidden = FALSE)
  outputOptions(output, "candidates_generated", suspendWhenHidden = FALSE)
  outputOptions(output, "preprocessing_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "segmentation_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "pmf_search_complete", suspendWhenHidden = FALSE)
  outputOptions(output, "summarization_complete", suspendWhenHidden = FALSE)
  
  # Pipeline status text
  output$pipeline_status_text <- renderText({
    if (is.null(pipeline_state$project)) {
      "No project loaded"
    } else {
      progress_bar(pipeline_state$project, width = 30)
    }
  })
  
  # Processing log
  output$pipeline_log <- renderText({
    paste(pipeline_state$log_entries, collapse = "\n")
  })
  
  # Candidates summary
  output$candidates_summary <- renderText({
    if (pipeline_state$status$candidates_generated && !is.null(pipeline_state$project)) {
      candidates <- pipeline_state$project$candidates
      if (!is.null(candidates)) {
        paste("Generated", nrow(candidates), "candidates |",
              "Mass range:", round(min(candidates$mz), 2), "-", round(max(candidates$mz), 2), "m/z")
      } else {
        "No candidates available"
      }
    } else {
      ""
    }
  })
  
  # Segmentation summary
  output$segmentation_summary <- renderText({
    if (pipeline_state$status$segmentation_complete && !is.null(pipeline_state$project)) {
      status <- hitmap_project_status(pipeline_state$project)
      paste("Segmentation completed for all files |",
            "Total regions:", status$results_summary$total_regions)
    } else {
      ""
    }
  })
  
  # PMF search summary
  output$pmf_search_summary <- renderText({
    if (pipeline_state$status$pmf_search_complete && !is.null(pipeline_state$project)) {
      status <- hitmap_project_status(pipeline_state$project)
      paste("PMF search completed |",
            "Peptides:", status$results_summary$total_peptides, "|",
            "Proteins:", status$results_summary$total_proteins)
    } else {
      ""
    }
  })
  
  # Final results text
  output$final_results_text <- renderText({
    if (pipeline_state$status$summarization_complete && !is.null(pipeline_state$project)) {
      capture.output(print(pipeline_state$project))
    } else {
      ""
    }
  })
  
  # Progress bar HTML
  output$pipeline_progress_bar <- renderUI({
    if (!is.null(pipeline_state$project)) {
      progress_text <- progress_bar(pipeline_state$project, width = 50)
      
      # Extract percentage
      pct_match <- regmatches(progress_text, regexpr("\\d+%", progress_text))
      percentage <- if (length(pct_match) > 0) as.numeric(gsub("%", "", pct_match)) else 0
      
      div(
        div(class = "progress",
            div(class = paste("progress-bar", if (percentage == 100) "bg-success" else "bg-primary"),
                style = paste0("width: ", percentage, "%"),
                paste0(percentage, "%")
            )
        ),
        br(),
        tags$code(progress_text)
      )
    } else {
      div(
        div(class = "progress",
            div(class = "progress-bar bg-light", style = "width: 0%", "0%")
        ),
        br(),
        tags$code("No project loaded")
      )
    }
  })
  
  # Results availability
  output$has_results <- reactive({
    pipeline_state$status$summarization_complete
  })
  outputOptions(output, "has_results", suspendWhenHidden = FALSE)
  
  # Results summaries
  output$total_peptides_count <- renderText({
    if (pipeline_state$status$summarization_complete && !is.null(pipeline_state$project$peptides)) {
      nrow(pipeline_state$project$peptides)
    } else {
      "0"
    }
  })
  
  output$total_proteins_count <- renderText({
    if (pipeline_state$status$summarization_complete && !is.null(pipeline_state$project$proteins)) {
      nrow(pipeline_state$project$proteins)
    } else {
      "0"
    }
  })
  
  output$total_regions_count <- renderText({
    if (!is.null(pipeline_state$project)) {
      status <- hitmap_project_status(pipeline_state$project)
      status$results_summary$total_regions
    } else {
      "0"
    }
  })
  
  output$total_files_count <- renderText({
    if (!is.null(pipeline_state$project)) {
      length(pipeline_state$project$metadata$data_files)
    } else {
      "0"
    }
  })
  
  # Data tables
  output$peptides_table <- DT::renderDataTable({
    if (pipeline_state$status$summarization_complete && !is.null(pipeline_state$project$peptides)) {
      DT::datatable(pipeline_state$project$peptides, 
                    options = list(scrollX = TRUE, pageLength = 25))
    }
  })
  
  output$proteins_table <- DT::renderDataTable({
    if (pipeline_state$status$summarization_complete && !is.null(pipeline_state$project$proteins)) {
      DT::datatable(pipeline_state$project$proteins, 
                    options = list(scrollX = TRUE, pageLength = 25))
    }
  })
  
  # Benchmark plot
  output$benchmarks_available <- reactive({
    length(pipeline_state$benchmarks) > 0
  })
  outputOptions(output, "benchmarks_available", suspendWhenHidden = FALSE)
  
  output$benchmark_plot <- renderPlot({
    if (length(pipeline_state$benchmarks) > 0) {
      benchmark_data <- data.frame(
        Step = names(pipeline_state$benchmarks),
        Time = unlist(pipeline_state$benchmarks)
      )
      
      ggplot2::ggplot(benchmark_data, ggplot2::aes(x = Step, y = Time)) +
        ggplot2::geom_col(fill = "steelblue") +
        ggplot2::theme_minimal() +
        ggplot2::labs(title = "Processing Times", y = "Time (seconds)") +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    }
  })
  
  # ===== CHECKPOINT AND EXPORT HANDLERS =====
  
  observeEvent(input$save_checkpoint, {
    if (!is.null(pipeline_state$project)) {
      timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
      filename <- paste0("hitmap_checkpoint_", timestamp, ".rds")
      
      tryCatch({
        pipeline_state$project %>% checkpoint(filename, "Manual checkpoint from GUI")
        showNotification(paste("Checkpoint saved:", filename), type = "success")
      }, error = function(e) {
        showNotification(paste("Error saving checkpoint:", e$message), type = "error")
      })
    } else {
      showNotification("No project to save", type = "warning")
    }
  })
  
  # Export handlers
  output$download_results <- downloadHandler(
    filename = function() {
      timestamp <- if (input$export_timestamp) {
        format(Sys.time(), "_%Y%m%d_%H%M%S")
      } else {
        ""
      }
      paste0(input$export_prefix, timestamp, ".zip")
    },
    content = function(file) {
      if (!is.null(pipeline_state$project)) {
        temp_dir <- tempdir()
        
        pipeline_state$project %>%
          export_results(
            output_folder = temp_dir,
            export_types = input$export_types,
            format = input$export_format
          )
        
        # Create zip file
        zip::zip(file, files = list.files(temp_dir, full.names = TRUE))
      }
    }
  )
  
  # Refresh status
  observeEvent(input$refresh_status, {
    update_pipeline_status()
    add_log_entry("Status refreshed")
  })
  
  return(pipeline_state)
}