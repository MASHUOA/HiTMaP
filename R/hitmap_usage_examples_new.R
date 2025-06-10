#' Usage Examples for New HiTMaP Modular System
#' 
#' Comprehensive examples using the new unified project object system
#' with consistent lowercase naming and flexible processing options

# =============================================================================
# Example 1: Quick Start Workflows
# =============================================================================

#' Example 1A: Simple quick start
run_example_simple_quick_start <- function() {
  
  # Most basic usage - process everything with defaults
  project <- hitmap_quick_start(
    data_files = "mouse_brain.imzML",
    database = "mouse_proteome.fasta"
  )
  
  # Check results
  print(project)
  status <- hitmap_project_status(project)
  print(status)
  
  return(project)
}

#' Example 1B: Quick start with custom parameters
run_example_custom_quick_start <- function() {
  
  # Quick start with some customization
  project <- hitmap_quick_start(
    data_files = c("sample1.imzML", "sample2.imzML"),
    database = "human_database.fasta",
    ppm = 3,                    # Higher precision
    threshold = 0.002,          # Higher threshold
    segment_count = 8,          # More detailed segmentation
    thread_count = 12,          # More threads
    output_summaries = TRUE
  )
  
  return(project)
}

# =============================================================================
# Example 2: Individual Processing Steps (Maximum Flexibility)
# =============================================================================

#' Example 2A: Candidate generation only
run_example_candidates_only <- function() {
  
  # Only generate candidates - useful for database preparation
  project <- hitmap_candidates_only(
    data_files = "data.imzML",
    database = "uniprot_human.fasta",
    digestion_site = "trypsin",
    missed_cleavages = 0:2,
    adducts = c("M+H", "M+Na", "M+K"),
    modifications = hitmap_create_modification_config("comprehensive"),
    mz_range = c(500, 5000),
    output_file = "comprehensive_candidates.csv"
  )
  
  # Inspect candidate statistics
  candidates <- project$results$candidates
  message(paste("Generated", nrow(candidates), "candidates"))
  message(paste("Mass range:", min(candidates$mz), "-", max(candidates$mz)))
  message(paste("Unique proteins:", length(unique(candidates$Protein))))
  
  return(project)
}

#' Example 2B: Preprocessing only workflow
run_example_preprocessing_only <- function() {
  
  # Only preprocess data - useful for data quality assessment
  project <- hitmap_preprocess_only(
    data_files = c("file1.imzML", "file2.imzML"),
    preprocess_profile = "comprehensive",  # More thorough preprocessing
    ppm = 2,                              # High precision
    mz_range = c(700, 4000),
    thread_count = 8
  )
  
  # Check preprocessing results
  for (file_name in project$status$files_preprocessed) {
    result <- project$results$preprocessing[[file_name]]
    message(paste("File:", file_name, "- Status:", result$status))
  }
  
  return(project)
}

#' Example 2C: Segmentation only workflow
run_example_segmentation_only <- function() {
  
  # Start with a project that has preprocessing done
  project <- hitmap_preprocess_only(
    data_files = "tissue_section.imzML",
    preprocess_profile = "standard"
  )
  
  # Now do segmentation only
  project <- hitmap_segment_only(
    project = project,
    segment_count = 12,                    # Fine spatial resolution
    segmentation_method = "spatialShrunkenCentroids",
    segmentation_params = list(
      variance_coverage = 0.9,             # Higher variance coverage
      smooth_range = 2                     # More smoothing
    )
  )
  
  # Analyze segmentation results
  for (file_name in project$status$files_segmented) {
    seg_result <- project$results$segmentation[[file_name]]
    regions <- names(seg_result$segmentation_labels)
    message(paste("File:", file_name))
    message(paste("  Regions:", paste(regions, collapse = ", ")))
    
    # Check region sizes
    for (region in regions) {
      pixel_count <- length(seg_result$segmentation_labels[[region]])
      message(paste("  ", region, ":", pixel_count, "pixels"))
    }
  }
  
  return(project)
}

#' Example 2D: PMF search only workflow
run_example_pmf_only <- function() {
  
  # Start with a complete preprocessing pipeline
  project <- hitmap_create_project("data.imzML")
  project <- hitmap_init(project, thread_count = 8)
  
  # Generate candidates
  project <- hitmap_generate_candidates(
    project = project,
    database = "database.fasta",
    modifications = hitmap_create_modification_config("basic")
  )
  
  # Preprocess and segment (assuming these are done)
  project <- hitmap_preprocess_file(project, "data", ppm = 5)
  project <- hitmap_segment_file(project, "data", segment_count = 6)
  
  # Now do PMF search only with custom parameters
  project <- hitmap_pmf_only(
    project = project,
    region_filter = function(regions) regions[1:3],  # Only search first 3 regions
    threshold = 0.0005,     # Lower threshold
    ppm = 3,                # Higher precision
    fdr_cutoff = 0.01,      # Stricter FDR
    peptide_filter = 3,     # Require more peptides per protein
    search_params = list(
      top_n_features = 2000,  # More features
      plot_scores = TRUE      # Generate score plots
    )
  )
  
  # Analyze PMF results
  for (file_name in project$status$files_pmf_searched) {
    pmf_results <- project$results$pmf_search[[file_name]]
    
    for (region in names(pmf_results)) {
      result <- pmf_results[[region]]
      message(paste("File:", file_name, "Region:", region))
      message(paste("  Peptides:", nrow(result$peptides)))
      message(paste("  Proteins:", nrow(result$proteins)))
    }
  }
  
  return(project)
}

# =============================================================================
# Example 3: Flexible Workflow Patterns
# =============================================================================

#' Example 3A: Step-by-step with custom configurations
run_example_flexible_steps <- function() {
  
  # Candidates only first
  project <- hitmap_flexible_workflow(
    data_files = "data.imzML",
    steps = c("candidates"),
    config = list(
      database = "custom_database.fasta",
      digestion_site = "[KR]|{P}",  # Custom enzyme
      modifications = hitmap_create_modification_config("comprehensive")
    )
  )
  
  # Inspect candidates before proceeding
  message(paste("Generated", nrow(project$results$candidates), "candidates"))
  
  # Continue with preprocessing and segmentation
  project <- hitmap_flexible_workflow(
    data_files = "data.imzML",
    steps = c("preprocess", "segment"),
    config = list(
      ppm = 3,
      preprocess_profile = "comprehensive",
      segment_count = 8,
      segmentation_method = "spatialKMeans"
    )
  )
  
  # Finally do PMF search with optimized parameters
  project <- hitmap_flexible_workflow(
    data_files = "data.imzML",
    steps = c("pmf", "summarize"),
    config = list(
      threshold = 0.001,
      fdr_cutoff = 0.05,
      peptide_filter = 2
    )
  )
  
  return(project)
}

#' Example 3B: Resume workflow from saved project
run_example_resume_workflow <- function() {
  
  # Start a workflow and save it
  project <- hitmap_flexible_workflow(
    data_files = "large_dataset.imzML",
    steps = c("candidates", "preprocess", "segment"),
    config = list(
      thread_count = 16,
      segment_count = 10
    )
  )
  
  # Save project
  saveRDS(project, "partial_analysis.rds")
  
  # Later... load and continue
  project <- hitmap_load_project("partial_analysis.rds")
  
  # Resume with PMF search and summary
  project <- hitmap_resume_workflow(
    project = project,
    additional_steps = c("pmf", "summarize"),
    config = list(
      threshold = 0.002,  # Updated parameters
      fdr_cutoff = 0.01
    )
  )
  
  return(project)
}

# =============================================================================
# Example 4: Advanced Configuration Patterns
# =============================================================================

#' Example 4A: High-precision proteomics analysis
run_example_high_precision <- function() {
  
  # Configuration for high-precision proteomics
  high_precision_config <- list(
    # High precision settings
    ppm = 2,
    import_ppm = 1,
    threshold = 0.005,
    fdr_cutoff = 0.01,
    peptide_filter = 3,
    
    # Comprehensive processing
    preprocess_profile = "comprehensive",
    segment_count = 12,
    segmentation_method = "spatialShrunkenCentroids",
    
    # Comprehensive modifications
    modifications = hitmap_create_modification_config("comprehensive"),
    missed_cleavages = 0:2,
    adducts = c("M+H", "M+Na", "M+K"),
    
    # Enhanced search
    score_method = "SQRTP",
    decoy_search = TRUE,
    mz_range = c(800, 4000),
    
    # Full summaries
    protein_summary = TRUE,
    peptide_summary = TRUE,
    region_summary = TRUE
  )
  
  project <- hitmap_advanced(
    data_files = "high_resolution_data.imzML",
    config = high_precision_config
  )
  
  return(project)
}

#' Example 4B: Discovery proteomics (permissive settings)
run_example_discovery_mode <- function() {
  
  # Configuration for discovery proteomics
  discovery_config <- list(
    # Permissive settings
    ppm = 10,
    threshold = 0.0001,
    fdr_cutoff = 0.1,
    peptide_filter = 1,
    
    # Less stringent processing
    preprocess_profile = "minimal",
    segment_count = 6,
    
    # Basic modifications only
    modifications = hitmap_create_modification_config("basic"),
    missed_cleavages = 0:3,  # Allow more missed cleavages
    
    # Wide mass range
    mz_range = c(400, 6000)
  )
  
  project <- hitmap_advanced(
    data_files = "discovery_sample.imzML",
    config = discovery_config
  )
  
  return(project)
}

#' Example 4C: Targeted analysis with specific regions
run_example_targeted_analysis <- function() {
  
  # Create project with custom segmentation
  project <- hitmap_create_project("tumor_sample.imzML")
  project <- hitmap_init(project, thread_count = 8)
  
  # Generate comprehensive candidates
  project <- hitmap_generate_candidates(
    project = project,
    database = "human_cancer_markers.fasta",
    modifications = hitmap_create_modification_config("comprehensive")
  )
  
  # Preprocess with high precision
  project <- hitmap_preprocess_file(
    project = project,
    file_name = "tumor_sample",
    ppm = 2,
    preprocess_params = hitmap_create_preprocess_config("comprehensive", 2)
  )
  
  # Segment into many regions for detailed analysis
  project <- hitmap_segment_file(
    project = project,
    file_name = "tumor_sample",
    segment_count = 15,
    segmentation_method = "spatialKMeans"
  )
  
  # Search only tumor core regions (hypothetical region selection)
  tumor_regions <- c("region_3", "region_7", "region_11")
  
  project <- hitmap_pmf_search(
    project = project,
    file_name = "tumor_sample",
    region_names = tumor_regions,
    threshold = 0.001,
    ppm = 2,
    fdr_cutoff = 0.01,
    peptide_filter = 3
  )
  
  # Generate summaries
  project <- hitmap_generate_summaries(project)
  
  return(project)
}

# =============================================================================
# Example 5: Batch Processing and Comparison
# =============================================================================

#' Example 5A: Batch process multiple conditions
run_example_batch_processing <- function() {
  
  # Define datasets for batch processing
  datasets <- list(
    control = list(
      data_files = c("control_1.imzML", "control_2.imzML", "control_3.imzML"),
      config = list(
        database = "mouse_proteome.fasta",
        threshold = 0.001,
        segment_count = 6
      )
    ),
    
    treatment_low = list(
      data_files = c("treat_low_1.imzML", "treat_low_2.imzML", "treat_low_3.imzML"),
      config = list(
        database = "mouse_proteome.fasta",
        threshold = 0.001,
        segment_count = 6
      )
    ),
    
    treatment_high = list(
      data_files = c("treat_high_1.imzML", "treat_high_2.imzML", "treat_high_3.imzML"),
      config = list(
        database = "mouse_proteome.fasta",
        threshold = 0.001,
        segment_count = 6
      )
    )
  )
  
  # Base configuration applied to all
  base_config <- list(
    ppm = 5,
    preprocess_profile = "standard",
    modifications = hitmap_create_modification_config("basic"),
    fdr_cutoff = 0.05,
    peptide_filter = 2
  )
  
  # Batch process
  batch_results <- hitmap_batch_process(
    dataset_list = datasets,
    base_config = base_config,
    parallel_datasets = TRUE,
    output_folder = "batch_analysis_results"
  )
  
  # Compare results
  comparison <- hitmap_compare_results(
    result_list = batch_results,
    comparison_type = "peptides",
    output_file = "condition_comparison"
  )
  
  return(list(batch_results = batch_results, comparison = comparison))
}

#' Example 5B: Parameter optimization study
run_example_parameter_optimization <- function() {
  
  # Test different parameter combinations
  param_combinations <- list(
    strict = list(
      ppm = 3,
      threshold = 0.005,
      fdr_cutoff = 0.01,
      peptide_filter = 3
    ),
    
    standard = list(
      ppm = 5,
      threshold = 0.001,
      fdr_cutoff = 0.05,
      peptide_filter = 2
    ),
    
    permissive = list(
      ppm = 10,
      threshold = 0.0005,
      fdr_cutoff = 0.1,
      peptide_filter = 1
    )
  )
  
  # Create datasets for each parameter set
  datasets <- list()
  for (param_name in names(param_combinations)) {
    datasets[[param_name]] <- list(
      data_files = "test_sample.imzML",
      config = param_combinations[[param_name]]
    )
  }
  
  # Process with different parameters
  optimization_results <- hitmap_batch_process(
    dataset_list = datasets,
    base_config = list(database = "test_database.fasta"),
    output_folder = "parameter_optimization"
  )
  
  # Analyze results
  for (param_name in names(optimization_results)) {
    if (!is.null(optimization_results[[param_name]])) {
      status <- hitmap_project_status(optimization_results[[param_name]])
      message(paste("Parameters:", param_name))
      message(paste("  Peptides:", status$results_summary$total_peptides))
      message(paste("  Proteins:", status$results_summary$total_proteins))
    }
  }
  
  return(optimization_results)
}

# =============================================================================
# Example 6: Quality Control and Validation
# =============================================================================

#' Example 6A: Quality control workflow
run_example_quality_control <- function() {
  
  # Start analysis with QC checks
  project <- hitmap_create_project("qc_sample.imzML")
  project <- hitmap_init(project, thread_count = 4)
  
  # QC Check 1: File accessibility
  if (!file.exists(file.path(project$setup$working_directory, "qc_sample.imzML"))) {
    stop("Data file not accessible!")
  }
  message("✓ Data file accessibility check passed")
  
  # Generate candidates with QC
  project <- hitmap_generate_candidates(project, database = "test_db.fasta")
  
  # QC Check 2: Candidate count
  n_candidates <- nrow(project$results$candidates)
  if (n_candidates < 100) {
    warning("Very few candidates generated. Check database and parameters.")
  } else if (n_candidates > 50000) {
    warning("Very large candidate list. Consider restricting parameters.")
  }
  message(paste("✓ Candidate count check:", n_candidates, "candidates"))
  
  # Preprocess with QC
  project <- hitmap_preprocess_file(project, "qc_sample")
  
  # Segment with QC
  project <- hitmap_segment_file(project, "qc_sample", segment_count = 6)
  
  # QC Check 3: Segmentation quality
  seg_result <- project$results$segmentation[["qc_sample"]]
  region_count <- length(seg_result$segmentation_labels)
  if (region_count < 2) {
    warning("Poor segmentation - only", region_count, "regions found")
  }
  message(paste("✓ Segmentation quality check:", region_count, "regions"))
  
  # PMF search with QC
  project <- hitmap_pmf_search(project, "qc_sample")
  
  # QC Check 4: Identification rates
  pmf_results <- project$results$pmf_search[["qc_sample"]]
  total_peptides <- sum(sapply(pmf_results, function(x) nrow(x$peptides)))
  
  if (total_peptides == 0) {
    stop("No peptides identified! Check parameters.")
  } else if (total_peptides < 10) {
    warning("Very few peptides identified. Consider relaxing parameters.")
  }
  message(paste("✓ Identification rate check:", total_peptides, "peptides"))
  
  # Generate summaries
  project <- hitmap_generate_summaries(project)
  
  message("✓ All QC checks completed successfully!")
  
  return(project)
}

# =============================================================================
# Example 7: Export and Reporting
# =============================================================================

#' Example 7A: Comprehensive result export
run_example_result_export <- function() {
  
  # Complete analysis
  project <- hitmap_quick_start(
    data_files = "sample.imzML",
    database = "database.fasta"
  )
  
  # Export all results in multiple formats
  
  # CSV exports
  csv_files <- hitmap_export_results(
    project = project,
    output_folder = "results_csv",
    export_types = c("summaries", "candidates"),
    format = "csv"
  )
  
  # Excel exports
  excel_files <- hitmap_export_results(
    project = project,
    output_folder = "results_excel",
    export_types = c("summaries", "candidates"),
    format = "excel"
  )
  
  # Complete project backup
  backup_files <- hitmap_export_results(
    project = project,
    output_folder = "project_backup",
    export_types = c("summaries", "candidates", "raw_results"),
    format = "rds"
  )
  
  return(list(
    csv_files = csv_files,
    excel_files = excel_files,
    backup_files = backup_files
  ))
}

# =============================================================================
# Usage Guidelines and Best Practices
# =============================================================================

#' Best practices for the new modular system:
#' 
#' 1. **Project Management:**
#'    - Always save projects with saveRDS() at key checkpoints
#'    - Use hitmap_project_status() to check progress
#'    - Use consistent lowercase naming for all parameters
#' 
#' 2. **Flexible Processing:**
#'    - Use hitmap_quick_start() for initial exploration
#'    - Use individual step functions for maximum control
#'    - Use hitmap_flexible_workflow() for custom step combinations
#'    - Use batch processing for multiple samples
#' 
#' 3. **Configuration Management:**
#'    - Use configuration helper functions for common setups
#'    - Store configurations in the project object
#'    - Test parameters with small datasets first
#' 
#' 4. **Quality Control:**
#'    - Always check intermediate results
#'    - Use appropriate QC thresholds
#'    - Validate segmentation quality before PMF search
#' 
#' 5. **Performance Optimization:**
#'    - Use appropriate thread counts
#'    - Consider memory usage for large datasets
#'    - Use batch processing for multiple samples
#' 
#' Parameter Guidelines:
#' - ppm: 2-3 for high-resolution, 5-10 for standard resolution
#' - threshold: 0.001-0.01 depending on data quality
#' - segment_count: 4-8 for most applications, more for detailed analysis
#' - fdr_cutoff: 0.01-0.05 for strict analysis, 0.1 for discovery
#' - peptide_filter: 2-3 for reliable protein identification