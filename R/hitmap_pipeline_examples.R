#' HiTMaP Pipeline Usage Examples
#' 
#' Comprehensive examples demonstrating the pipeline-style workflow
#' using %>% operator, similar to Cardinal and dplyr

# =============================================================================
# Example 1: Basic Pipeline Workflows
# =============================================================================

#' Example 1A: Simple linear pipeline
run_example_basic_pipeline <- function() {
  
  # Basic pipeline - most common usage pattern
  result <- hitmap_project("mouse_brain.imzML") %>%
    init(thread_count = 8) %>%
    generate_candidates(database = "mouse_proteome.fasta") %>%
    preprocess() %>%
    segment(segment_count = 6) %>%
    search_pmf() %>%
    summarize()
  
  # Check results
  print(result)
  
  return(result)
}

#' Example 1B: Pipeline with custom configuration
run_example_configured_pipeline <- function() {
  
  # Pipeline with inline configuration
  result <- hitmap_project("sample.imzML") %>%
    init(thread_count = 4) %>%
    apply_config(config_strict()) %>%
    generate_candidates(
      database = "human_proteome.fasta",
      modifications = hitmap_create_modification_config("comprehensive"),
      adducts = c("M+H", "M+Na", "M+K")
    ) %>%
    preprocess(
      ppm = 3,
      preprocess_profile = "comprehensive"
    ) %>%
    segment(
      segment_count = 8,
      segmentation_method = "spatialShrunkenCentroids"
    ) %>%
    search_pmf(
      threshold = 0.005,
      fdr_cutoff = 0.01,
      peptide_filter = 3
    ) %>%
    summarize()
  
  return(result)
}

#' Example 1C: Pipeline with checkpointing
run_example_pipeline_with_checkpoints <- function() {
  
  # Pipeline with regular checkpointing
  result <- hitmap_project("large_dataset.imzML") %>%
    init(thread_count = 12) %>%
    checkpoint("checkpoint_01_init.rds", "Project initialized") %>%
    
    generate_candidates(database = "large_database.fasta") %>%
    checkpoint("checkpoint_02_candidates.rds", "Candidates generated") %>%
    
    preprocess() %>%
    checkpoint("checkpoint_03_preprocess.rds", "Preprocessing complete") %>%
    
    segment(segment_count = 10) %>%
    checkpoint("checkpoint_04_segment.rds", "Segmentation complete") %>%
    
    search_pmf() %>%
    checkpoint("checkpoint_05_pmf.rds", "PMF search complete") %>%
    
    summarize() %>%
    checkpoint("checkpoint_06_final.rds", "Analysis complete")
  
  return(result)
}

# =============================================================================
# Example 2: Conditional and Branching Pipelines
# =============================================================================

#' Example 2A: Conditional pipeline execution
run_example_conditional_pipeline <- function() {
  
  # Pipeline with conditional steps
  result <- hitmap_project("data.imzML") %>%
    init() %>%
    generate_candidates(database = "database.fasta") %>%
    
    # Conditional preprocessing based on data size
    if_then(
      condition = function(p) length(p$metadata$data_files) > 5,
      true_action = function(p) p %>% preprocess(preprocess_profile = "minimal"),
      false_action = function(p) p %>% preprocess(preprocess_profile = "comprehensive")
    ) %>%
    
    segment() %>%
    
    # Conditional search parameters based on number of candidates
    if_then(
      condition = function(p) nrow(p$results$candidates) > 10000,
      true_action = function(p) p %>% search_pmf(threshold = 0.005, fdr_cutoff = 0.01),
      false_action = function(p) p %>% search_pmf(threshold = 0.001, fdr_cutoff = 0.05)
    ) %>%
    
    summarize()
  
  return(result)
}

#' Example 2B: Pipeline branching for parameter testing
run_example_pipeline_branching <- function() {
  
  # Create base pipeline
  base_project <- hitmap_project("test_sample.imzML") %>%
    init() %>%
    generate_candidates(database = "test_db.fasta") %>%
    preprocess() %>%
    segment()
  
  # Create branches for different search parameters
  strict_branch <- branch_pipeline(
    base_project, 
    "strict_search",
    function(p) p %>% search_pmf(threshold = 0.01, fdr_cutoff = 0.01, peptide_filter = 3),
    function(p) p %>% summarize()
  )
  
  permissive_branch <- branch_pipeline(
    base_project,
    "permissive_search", 
    function(p) p %>% search_pmf(threshold = 0.0001, fdr_cutoff = 0.1, peptide_filter = 1),
    function(p) p %>% summarize()
  )
  
  # Compare results
  strict_peptides <- nrow(strict_branch$branch$peptides)
  permissive_peptides <- nrow(permissive_branch$branch$peptides)
  
  message(paste("Strict search:", strict_peptides, "peptides"))
  message(paste("Permissive search:", permissive_peptides, "peptides"))
  
  # Choose best branch and merge
  if (strict_peptides > 100) {
    final_result <- merge_branches(base_project, strict_branch$branch, "replace")
  } else {
    final_result <- merge_branches(base_project, permissive_branch$branch, "replace")
  }
  
  return(final_result)
}

# =============================================================================
# Example 3: Pipeline Templates and Reusable Workflows
# =============================================================================

#' Example 3A: Create and use pipeline templates
run_example_pipeline_templates <- function() {
  
  # Define a high-throughput template
  high_throughput_template <- create_template(
    "high_throughput",
    function(p) p %>% init(thread_count = 16),
    function(p) p %>% generate_candidates(database = "large_database.fasta"),
    function(p) p %>% preprocess(preprocess_profile = "minimal"),
    function(p) p %>% segment(segment_count = 4),
    function(p) p %>% search_pmf(threshold = 0.01),
    function(p) p %>% summarize()
  )
  
  # Define a high-precision template
  high_precision_template <- create_template(
    "high_precision",
    function(p) p %>% init(thread_count = 8),
    function(p) p %>% generate_candidates(
      database = "curated_database.fasta",
      modifications = hitmap_create_modification_config("comprehensive")
    ),
    function(p) p %>% preprocess(
      ppm = 2, 
      preprocess_profile = "comprehensive"
    ),
    function(p) p %>% segment(segment_count = 12),
    function(p) p %>% search_pmf(
      threshold = 0.005, 
      ppm = 2, 
      fdr_cutoff = 0.01
    ),
    function(p) p %>% summarize()
  )
  
  # Apply templates to different datasets
  htp_result <- hitmap_project("screening_sample.imzML") %>%
    apply_template(high_throughput_template)
  
  precision_result <- hitmap_project("discovery_sample.imzML") %>%
    apply_template(high_precision_template)
  
  return(list(
    high_throughput = htp_result,
    high_precision = precision_result
  ))
}

#' Example 3B: Reusable workflow components
run_example_reusable_components <- function() {
  
  # Define reusable workflow components
  basic_setup <- function(project, db_name) {
    project %>%
      init(thread_count = 8) %>%
      generate_candidates(database = db_name) %>%
      preprocess()
  }
  
  detailed_segmentation <- function(project, n_segments = 8) {
    project %>%
      segment(
        segment_count = n_segments,
        segmentation_method = "spatialKMeans",
        segmentation_params = list(variance_coverage = 0.9)
      )
  }
  
  strict_search <- function(project) {
    project %>%
      search_pmf(
        threshold = 0.005,
        fdr_cutoff = 0.01,
        peptide_filter = 3
      )
  }
  
  # Use components in pipeline
  result <- hitmap_project("tissue_sample.imzML") %>%
    basic_setup("tissue_proteome.fasta") %>%
    detailed_segmentation(n_segments = 10) %>%
    strict_search() %>%
    summarize()
  
  return(result)
}

# =============================================================================
# Example 4: Error Handling and Debugging
# =============================================================================

#' Example 4A: Pipeline with error handling
run_example_error_handling_pipeline <- function() {
  
  # Pipeline with comprehensive error handling
  result <- hitmap_project("problematic_data.imzML") %>%
    try_step(
      action = function(p) p %>% init(thread_count = 8),
      on_error = function(p, e) {
        message("Init failed, trying with fewer threads")
        p %>% init(thread_count = 2)
      }
    ) %>%
    
    try_step(
      action = function(p) p %>% generate_candidates(database = "primary_db.fasta"),
      on_error = function(p, e) {
        message("Primary database failed, using backup")
        p %>% generate_candidates(database = "backup_db.fasta")
      }
    ) %>%
    
    validate_state(c("initialized", "candidates"), action = "warn") %>%
    
    try_step(
      action = function(p) p %>% preprocess(preprocess_profile = "comprehensive"),
      on_error = function(p, e) {
        message("Comprehensive preprocessing failed, using minimal")
        p %>% preprocess(preprocess_profile = "minimal")
      }
    ) %>%
    
    segment() %>%
    search_pmf() %>%
    summarize()
  
  return(result)
}

#' Example 4B: Pipeline with debugging and benchmarking
run_example_debug_pipeline <- function() {
  
  # Pipeline with debugging and performance monitoring
  result <- hitmap_project("performance_test.imzML") %>%
    set_options(debug_mode = TRUE, verbose = TRUE) %>%
    
    debug_step("Starting analysis pipeline") %>%
    
    benchmark_step(
      "initialization",
      function(p) p %>% init(thread_count = 8)
    ) %>%
    
    debug_step("Initialization complete, generating candidates") %>%
    
    benchmark_step(
      "candidate_generation",
      function(p) p %>% generate_candidates(database = "large_db.fasta")
    ) %>%
    
    debug_step("Candidates generated, starting preprocessing") %>%
    
    benchmark_step(
      "preprocessing",
      function(p) p %>% preprocess()
    ) %>%
    
    debug_step("Preprocessing complete, performing segmentation") %>%
    
    benchmark_step(
      "segmentation", 
      function(p) p %>% segment(segment_count = 8)
    ) %>%
    
    debug_step("Segmentation complete, starting PMF search") %>%
    
    benchmark_step(
      "pmf_search",
      function(p) p %>% search_pmf()
    ) %>%
    
    debug_step("PMF search complete, generating summaries") %>%
    
    benchmark_step(
      "summarization",
      function(p) p %>% summarize()
    ) %>%
    
    debug_step("Analysis pipeline complete")
  
  # Print benchmark results
  if (!is.null(result$metadata$benchmarks)) {
    message("Performance Summary:")
    for (step_name in names(result$metadata$benchmarks)) {
      duration <- result$metadata$benchmarks[[step_name]]$duration_seconds
      message(paste("  ", step_name, ":", round(duration, 2), "seconds"))
    }
  }
  
  return(result)
}

# =============================================================================
# Example 5: Data Filtering and Selection
# =============================================================================

#' Example 5A: Pipeline with data filtering
run_example_filtering_pipeline <- function() {
  
  # Pipeline with multiple data filtering steps
  result <- hitmap_project(c("sample1.imzML", "sample2.imzML", "sample3.imzML")) %>%
    init() %>%
    
    # Filter to specific files if needed
    select_data(files = c("sample1", "sample2")) %>%
    
    generate_candidates(database = "database.fasta") %>%
    
    # Filter candidates by mass range
    filter_data(
      filter_type = "candidates",
      filter_function = function(candidates) candidates$mz >= 800 & candidates$mz <= 3000
    ) %>%
    
    preprocess() %>%
    segment(segment_count = 6) %>%
    
    # Filter to specific regions of interest
    filter_data(
      filter_type = "regions", 
      filter_function = function(regions) regions %in% c("region_1", "region_3", "region_5")
    ) %>%
    
    search_pmf() %>%
    summarize()
  
  return(result)
}

#' Example 5B: Pipeline with dynamic region selection
run_example_dynamic_region_pipeline <- function() {
  
  # Pipeline that adapts region selection based on segmentation results
  result <- hitmap_project("heterogeneous_tissue.imzML") %>%
    init() %>%
    generate_candidates(database = "tissue_db.fasta") %>%
    preprocess() %>%
    segment(segment_count = 12) %>%
    
    # Dynamically select regions based on size
    filter_data(
      filter_type = "regions",
      filter_function = function(regions) {
        # This would normally check actual region sizes
        # For demo, select every other region
        seq(1, length(regions), by = 2)
      }
    ) %>%
    
    search_pmf(
      region_filter = function(regions) {
        # Additional filtering at search time
        # Select regions with good signal quality (hypothetical)
        regions[c(1, 3, 5)]  # Example selection
      }
    ) %>%
    
    summarize()
  
  return(result)
}

# =============================================================================
# Example 6: Multiple File and Batch Processing
# =============================================================================

#' Example 6A: Pipeline for multiple files
run_example_multiple_files_pipeline <- function() {
  
  # Process multiple files in a single pipeline
  result <- hitmap_project(c("control_1.imzML", "control_2.imzML", 
                           "treatment_1.imzML", "treatment_2.imzML")) %>%
    init(thread_count = 16) %>%
    generate_candidates(database = "experiment_db.fasta") %>%
    
    # Process all files
    preprocess(ppm = 5) %>%
    segment(segment_count = 6) %>%
    search_pmf() %>%
    summarize()
  
  # Access results by file
  for (file_name in names(result)) {
    file_peptides <- sum(sapply(result$results$pmf_search[[file_name]], 
                               function(x) nrow(x$peptides)))
    message(paste("File", file_name, ":", file_peptides, "peptides"))
  }
  
  return(result)
}

#' Example 6B: Comparative pipeline analysis
run_example_comparative_pipeline <- function() {
  
  # Create separate pipelines for comparison
  control_pipeline <- hitmap_project(c("ctrl_1.imzML", "ctrl_2.imzML")) %>%
    init() %>%
    generate_candidates(database = "mouse_db.fasta") %>%
    preprocess() %>%
    segment() %>%
    search_pmf() %>%
    summarize()
  
  treatment_pipeline <- hitmap_project(c("treat_1.imzML", "treat_2.imzML")) %>%
    init() %>%
    generate_candidates(database = "mouse_db.fasta") %>%
    preprocess() %>%
    segment() %>%
    search_pmf() %>%
    summarize()
  
  # Combine for comparative analysis
  combined_analysis <- c(control_pipeline, treatment_pipeline) %>%
    summarize()
  
  # Export comparative results
  combined_analysis %>%
    export_results(
      output_folder = "comparative_analysis",
      export_types = c("summaries", "raw_results")
    )
  
  return(combined_analysis)
}

# =============================================================================
# Example 7: Advanced Pipeline Patterns
# =============================================================================

#' Example 7A: Iterative optimization pipeline
run_example_iterative_pipeline <- function() {
  
  # Pipeline that iteratively optimizes parameters
  result <- hitmap_project("optimization_sample.imzML") %>%
    init() %>%
    generate_candidates(database = "test_db.fasta") %>%
    preprocess() %>%
    segment() %>%
    
    # Iterative parameter optimization
    repeat_until(
      action = function(p) {
        # Try search with current threshold
        current_threshold <- get_option(p, "search_threshold", 0.01)
        
        p %>%
          search_pmf(threshold = current_threshold) %>%
          mutate_config(search_threshold = current_threshold * 0.8)  # Reduce threshold
      },
      condition = function(p) {
        # Stop when we have enough results or threshold is too low
        total_peptides <- sum(sapply(p$results$pmf_search, 
                                   function(file_results) {
                                     sum(sapply(file_results, function(x) nrow(x$peptides)))
                                   }))
        threshold <- get_option(p, "search_threshold", 0.01)
        
        return(total_peptides >= 100 || threshold < 0.0001)
      },
      max_iterations = 5
    ) %>%
    
    summarize()
  
  return(result)
}

#' Example 7B: Progressive complexity pipeline
run_example_progressive_pipeline <- function() {
  
  # Pipeline that progressively increases analysis complexity
  result <- hitmap_project("complex_sample.imzML") %>%
    init() %>%
    
    # Start with basic candidates
    generate_candidates(
      database = "core_proteome.fasta",
      modifications = hitmap_create_modification_config("none")
    ) %>%
    
    preprocess(preprocess_profile = "minimal") %>%
    segment(segment_count = 4) %>%
    search_pmf() %>%
    
    # Check if we need more complexity
    if_then(
      condition = function(p) {
        total_peptides <- nrow(p$peptides %||% data.frame())
        return(total_peptides < 50)  # Need more candidates
      },
      true_action = function(p) {
        message("Adding complexity: comprehensive modifications")
        p %>%
          generate_candidates(
            database = "extended_proteome.fasta",
            modifications = hitmap_create_modification_config("comprehensive"),
            adducts = c("M+H", "M+Na", "M+K")
          ) %>%
          preprocess(preprocess_profile = "comprehensive") %>%
          segment(segment_count = 8) %>%
          search_pmf(threshold = 0.0005)
      }
    ) %>%
    
    summarize()
  
  return(result)
}

# =============================================================================
# Example 8: Pipeline Utilities and Monitoring
# =============================================================================

#' Example 8A: Pipeline with progress monitoring
run_example_monitored_pipeline <- function() {
  
  # Pipeline with progress tracking
  result <- hitmap_project("monitored_sample.imzML") %>%
    init() %>%
    view_status() %>%  # Show initial status
    
    generate_candidates(database = "database.fasta") %>%
    debug_step(paste("Progress:", progress_bar(.))) %>%
    
    preprocess() %>%
    debug_step(paste("Progress:", progress_bar(.))) %>%
    
    segment() %>%
    debug_step(paste("Progress:", progress_bar(.))) %>%
    
    search_pmf() %>%
    debug_step(paste("Progress:", progress_bar(.))) %>%
    
    summarize() %>%
    debug_step(paste("Progress:", progress_bar(.))) %>%
    view_status(detailed = TRUE)  # Show final detailed status
  
  return(result)
}

#' Example 8B: Pipeline with comprehensive logging
run_example_logged_pipeline <- function() {
  
  # Pipeline with detailed logging and export
  result <- hitmap_project("logged_analysis.imzML") %>%
    set_options(
      debug_mode = TRUE,
      log_level = "detailed",
      export_logs = TRUE
    ) %>%
    
    debug_step("=== Starting HiTMaP Analysis ===", level = "info") %>%
    init() %>%
    debug_step("Project initialized successfully") %>%
    
    generate_candidates(database = "proteome.fasta") %>%
    debug_step(paste("Generated", nrow(.$candidates), "candidates")) %>%
    
    preprocess() %>%
    debug_step("Data preprocessing completed") %>%
    
    segment(segment_count = 6) %>%
    debug_step("Tissue segmentation completed") %>%
    
    search_pmf() %>%
    debug_step("PMF search completed") %>%
    
    summarize() %>%
    debug_step("Results summarized") %>%
    
    export_results(
      output_folder = "logged_results",
      export_types = c("summaries", "candidates", "raw_results")
    ) %>%
    
    debug_step("=== Analysis Complete ===", level = "info")
  
  # Export processing log
  writeLines(result$log, "processing_log.txt")
  
  return(result)
}

# =============================================================================
# Usage Summary and Best Practices
# =============================================================================

#' Best practices for HiTMaP pipeline usage:
#' 
#' 1. **Pipeline Structure:**
#'    - Always start with hitmap_project() or load_project()
#'    - Use init() as the first pipeline step
#'    - Follow logical order: init -> candidates -> preprocess -> segment -> search_pmf -> summarize
#' 
#' 2. **Error Handling:**
#'    - Use try_step() for potentially problematic operations
#'    - Use validate_state() to check prerequisites
#'    - Use checkpoint() for long-running analyses
#' 
#' 3. **Performance:**
#'    - Use benchmark_step() to identify bottlenecks
#'    - Use appropriate thread counts in init()
#'    - Consider data filtering to reduce processing time
#' 
#' 4. **Debugging:**
#'    - Use debug_step() for progress monitoring
#'    - Use view_status() to check pipeline state
#'    - Use progress_bar() for visual progress tracking
#' 
#' 5. **Reusability:**
#'    - Create templates for common workflows
#'    - Use configuration functions for parameter sets
#'    - Save intermediate results with checkpoint()
#' 
#' 6. **Data Management:**
#'    - Use filter_data() and select_data() for subsetting
#'    - Use export_results() for standardized output
#'    - Use proper file naming for multiple analyses