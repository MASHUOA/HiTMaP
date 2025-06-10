#' Usage Examples for Modular Imaging Identification Workflow
#' 
#' This file contains comprehensive examples of how to use the new modular
#' imaging identification functions for step-by-step control.

# =============================================================================
# Example 1: Quick Start - Simple Analysis
# =============================================================================

#' Example 1A: Basic quick start with minimal configuration
run_example_quick_start_basic <- function() {
  
  # Simple analysis with default parameters
  results <- hitmap_quick_start(
    datafile = "path/to/your/data.imzML",
    database = "uniprot-bovin.fasta"
  )
  
  return(results)
}

#' Example 1B: Quick start with custom parameters
run_example_quick_start_custom <- function() {
  
  # Quick start with some customization
  results <- hitmap_quick_start(
    datafile = c("file1.imzML", "file2.imzML"),  # Multiple files
    database = "human_database.fasta",
    ppm = 10,                    # Higher mass tolerance
    threshold = 0.005,           # Higher intensity threshold
    segments = 6,                # More spatial segments
    threads = 8                  # More processing threads
  )
  
  return(results)
}

# =============================================================================
# Example 2: Advanced Configuration Workflow
# =============================================================================

#' Example 2A: Advanced workflow with comprehensive modifications
run_example_advanced_comprehensive <- function() {
  
  # Define custom configuration
  config <- list(
    database = "mouse_proteome.fasta",
    threshold = 0.002,
    ppm = 3,                     # High precision
    Thread = 12,
    segmentation_num = 8,        # More detailed segmentation
    Segmentation = "spatialShrunkenCentroids",
    preprocess_profile = "comprehensive",  # More thorough preprocessing
    modification_profile = "comprehensive", # More modifications
    Digestion_site = "trypsin",
    missedCleavages = 0:2,       # Allow more missed cleavages
    adducts = c("M+H", "M+Na", "M+K"),  # Multiple adducts
    FDR_cutoff = 0.01,          # Stricter FDR
    peptide_ID_filter = 3,       # Require more peptides per protein
    score_method = "SQRTP",
    mzrange = c(500, 5000),     # Wider mass range
    Region_feature_summary = TRUE  # Include region summary
  )
  
  results <- hitmap_advanced_workflow(
    datafile = "high_resolution_data.imzML",
    config_list = config
  )
  
  return(results)
}

#' Example 2B: Metabolomics-focused configuration
run_example_metabolomics_config <- function() {
  
  config <- list(
    threshold = 0.0001,          # Lower threshold for metabolites
    ppm = 2,                     # High precision for small molecules
    segmentation_num = 12,       # Fine spatial resolution
    preprocess_profile = "minimal",  # Less aggressive preprocessing
    modification_profile = "none",   # No protein modifications
    mzrange = c(50, 1000),      # Lower mass range
    adducts = c("M+H", "M+Na", "M+K", "M+NH4", "M-H"),  # Common metabolite adducts
    FDR_cutoff = 0.1            # More permissive for discovery
  )
  
  results <- hitmap_advanced_workflow(
    datafile = "metabolomics_data.imzML",
    config_list = config
  )
  
  return(results)
}

# =============================================================================
# Example 3: Step-by-Step Interactive Workflow
# =============================================================================

#' Example 3A: Basic step-by-step workflow
run_example_step_by_step_basic <- function() {
  
  # Create workflow controller
  controller <- hitmap_step_by_step(
    datafile = "data.imzML",
    config_list = list(
      database = "uniprot-bovin.fasta",
      ppm = 5,
      threshold = 0.001
    )
  )
  
  # Step 1: Initialize
  controller <- run_step_init(controller)
  print(controller)  # Check status
  
  # Step 2: Generate candidates
  controller <- run_step_candidates(controller)
  print(paste("Generated", nrow(controller$results$candidate_list), "candidates"))
  
  # Step 3: Preprocess and segment
  controller <- run_step_preprocess(controller, file_index = 1)
  print(paste("Found regions:", paste(names(controller$results$segmented_data[[1]]$segmentation_label), collapse = ", ")))
  
  # Step 4: PMF search
  controller <- run_step_pmf(controller, file_index = 1)
  print(paste("Identified", nrow(controller$results$file_results[[1]]$combined_peptides), "peptides"))
  
  # Step 5: Summarize
  controller <- run_step_summarize(controller)
  print("Workflow completed!")
  
  return(controller)
}

#' Example 3B: Step-by-step with parameter inspection and modification
run_example_step_by_step_interactive <- function() {
  
  # Initial configuration
  config <- list(
    database = "database.fasta",
    ppm = 5,
    threshold = 0.001,
    segmentation_num = 4
  )
  
  controller <- hitmap_step_by_step("data.imzML", config)
  
  # Step 1: Initialize
  controller <- run_step_init(controller)
  
  # Step 2: Generate candidates
  controller <- run_step_candidates(controller)
  
  # Inspect candidate list
  candidates <- controller$results$candidate_list
  print(paste("Candidate statistics:"))
  print(paste("Total candidates:", nrow(candidates)))
  print(paste("Mass range:", min(candidates$mz), "-", max(candidates$mz)))
  
  # Step 3: Preprocess - maybe adjust based on data quality
  controller <- run_step_preprocess(controller)
  
  # Inspect segmentation results
  seg_data <- controller$results$segmented_data[[1]]
  print(paste("Segmentation results:"))
  print(paste("Regions found:", length(seg_data$segmentation_label)))
  for (region in names(seg_data$segmentation_label)) {
    print(paste("Region", region, "- pixels:", length(seg_data$segmentation_label[[region]])))
  }
  
  # Decide if segmentation looks good, maybe adjust threshold
  # controller$config$threshold <- 0.002  # Increase if needed
  
  # Step 4: PMF search
  controller <- run_step_pmf(controller)
  
  # Inspect PMF results
  pmf_results <- controller$results$file_results[[1]]
  print("PMF results by region:")
  for (region in names(pmf_results$regions)) {
    region_result <- pmf_results$regions[[region]]
    print(paste("Region", region, "- peptides:", nrow(region_result$peptides), 
                "proteins:", nrow(region_result$proteins)))
  }
  
  # Step 5: Summarize
  controller <- run_step_summarize(controller)
  
  return(controller)
}

# =============================================================================
# Example 4: Individual Step Functions - Maximum Control
# =============================================================================

#' Example 4A: Manual control of each individual step
run_example_manual_control <- function() {
  
  # Step 1: Initialize project manually
  project_setup <- hitmap_init_project(
    datafile = "data.imzML",
    projectfolder = "/path/to/project",
    Thread = 8
  )
  
  # Step 2: Generate candidates with specific parameters
  candidate_list <- hitmap_generate_candidates(
    workdir = project_setup$workdir,
    database = "custom_database.fasta",
    Digestion_site = "[KR]|{P}",  # Custom digestion pattern
    missedCleavages = 0:1,
    adducts = c("M+H", "M+Na"),
    Modifications = list(
      fixed = c("Carbamidomethyl (C)"),
      fixmod_position = c("any"),
      variable = c("Oxidation (M)", "Acetyl (Protein N-term)"),
      varmod_position = c("any", "any")
    ),
    mzrange = c(800, 3500),
    BPPARAM = project_setup$BPPARAM
  )
  
  # Step 3: Preprocess with custom parameters
  preprocess_params <- list(
    force_preprocess = TRUE,
    use_preprocessRDS = FALSE,
    smoothSignal = list(method = "gaussian", sd = 0.5),
    reduceBaseline = list(method = "locmin", blocks = 200),
    peakPick = list(method = "adaptive", SNR = 2),
    peakAlign = list(tolerance = 2, units = "ppm"),
    peakFilter = list(freq.min = 0.1),
    normalize = list(method = "tic", mz = 1)
  )
  
  segmented_data <- hitmap_preprocess_segment(
    datafile = project_setup$datafile[1],
    workdir = project_setup$workdir,
    ppm = 3,
    segmentation_num = 6,
    Segmentation = "spatialKMeans",
    preprocess = preprocess_params,
    BPPARAM = project_setup$BPPARAM
  )
  
  # Step 4: Process individual regions with different parameters
  region_results <- list()
  
  for (region in names(segmented_data$segmentation_label)) {
    # Customize parameters per region if needed
    region_threshold <- if (region == "region_1") 0.001 else 0.002
    
    result <- hitmap_pmf_search_region(
      segmented_data = segmented_data,
      region_name = region,
      candidate_list = candidate_list,
      threshold = region_threshold,
      ppm = 3,
      score_method = "SQRTP",
      FDR_cutoff = 0.05,
      peptide_ID_filter = 2,
      BPPARAM = project_setup$BPPARAM
    )
    
    region_results[[region]] <- result
  }
  
  # Manual summary
  all_peptides <- do.call(rbind, lapply(region_results, function(x) x$peptides))
  all_proteins <- do.call(rbind, lapply(region_results, function(x) x$proteins))
  
  return(list(
    project_setup = project_setup,
    candidate_list = candidate_list,
    segmented_data = segmented_data,
    region_results = region_results,
    summary = list(
      total_peptides = nrow(all_peptides),
      total_proteins = nrow(all_proteins),
      regions_processed = length(region_results)
    )
  ))
}

# =============================================================================
# Example 5: Parameter Configuration Helpers
# =============================================================================

#' Example 5A: Using preprocessing configuration helpers
run_example_preprocess_configs <- function() {
  
  # Compare different preprocessing approaches
  
  # Minimal preprocessing
  minimal_params <- hitmap_create_preprocess_params("minimal", ppm = 5)
  
  # Standard preprocessing
  standard_params <- hitmap_create_preprocess_params("standard", ppm = 5)
  
  # Comprehensive preprocessing
  comprehensive_params <- hitmap_create_preprocess_params("comprehensive", ppm = 5)
  
  print("Preprocessing configurations:")
  print("Minimal:")
  print(minimal_params)
  print("Standard:")
  print(standard_params)
  print("Comprehensive:")
  print(comprehensive_params)
  
  # Use in workflow
  results_minimal <- hitmap_modular_workflow(
    datafile = "data.imzML",
    preprocess = minimal_params
  )
  
  return(results_minimal)
}

#' Example 5B: Using modification configuration helpers
run_example_modification_configs <- function() {
  
  # Different modification profiles
  no_mods <- hitmap_create_modification_params("none")
  basic_mods <- hitmap_create_modification_params("basic")
  comprehensive_mods <- hitmap_create_modification_params("comprehensive")
  
  print("Modification configurations:")
  print("None:")
  print(no_mods)
  print("Basic:")
  print(basic_mods)
  print("Comprehensive:")
  print(comprehensive_mods)
  
  # Custom modifications
  custom_mods <- list(
    fixed = c("Carbamidomethyl (C)", "TMT6plex (N-term)", "TMT6plex (K)"),
    fixmod_position = c("any", "any", "any"),
    variable = c("Oxidation (M)"),
    varmod_position = c("any")
  )
  
  # Use in workflow
  results <- hitmap_modular_workflow(
    datafile = "tmt_labeled_data.imzML",
    Modifications = custom_mods
  )
  
  return(results)
}

# =============================================================================
# Example 6: Multiple File Processing Strategies
# =============================================================================

#' Example 6A: Sequential processing of multiple files
run_example_multiple_files_sequential <- function() {
  
  files <- c("sample1.imzML", "sample2.imzML", "sample3.imzML")
  
  # Process all files in one workflow
  results <- hitmap_modular_workflow(
    datafile = files,
    database = "database.fasta",
    threshold = 0.001,
    ppm = 5
  )
  
  # Check results for each file
  for (file in files) {
    file_result <- results$file_results[[gsub(".imzML", "", file)]]
    print(paste("File:", file))
    print(paste("  Peptides:", nrow(file_result$combined_peptides)))
    print(paste("  Regions:", length(file_result$regions)))
  }
  
  return(results)
}

#' Example 6B: Individual file processing with comparison
run_example_multiple_files_individual <- function() {
  
  files <- c("control.imzML", "treatment.imzML")
  all_results <- list()
  
  # Process each file individually for comparison
  for (file in files) {
    results <- hitmap_modular_workflow(
      datafile = file,
      database = "database.fasta",
      threshold = 0.001,
      ppm = 5
    )
    
    all_results[[file]] <- results
  }
  
  # Compare results
  control_peptides <- all_results[["control.imzML"]]$summary_results$Peptide_Summary
  treatment_peptides <- all_results[["treatment.imzML"]]$summary_results$Peptide_Summary
  
  print("Comparison results:")
  print(paste("Control peptides:", nrow(control_peptides)))
  print(paste("Treatment peptides:", nrow(treatment_peptides)))
  
  # Find unique and common peptides
  control_mz <- control_peptides$mz
  treatment_mz <- treatment_peptides$mz
  
  common_mz <- intersect(control_mz, treatment_mz)
  control_unique <- setdiff(control_mz, treatment_mz)
  treatment_unique <- setdiff(treatment_mz, control_mz)
  
  print(paste("Common peptides:", length(common_mz)))
  print(paste("Control-specific:", length(control_unique)))
  print(paste("Treatment-specific:", length(treatment_unique)))
  
  return(all_results)
}

# =============================================================================
# Example 7: Quality Control and Validation
# =============================================================================

#' Example 7A: Workflow with quality control checks
run_example_with_qc <- function() {
  
  # Start step-by-step workflow for QC
  controller <- hitmap_step_by_step(
    datafile = "data.imzML",
    config_list = list(database = "database.fasta")
  )
  
  # Step 1: Initialize
  controller <- run_step_init(controller)
  
  # QC Check 1: File accessibility
  workdir <- controller$results$project_setup$workdir
  datafile <- controller$results$project_setup$datafile_imzML[1]
  if (!file.exists(file.path(workdir, datafile))) {
    stop("Data file not accessible!")
  }
  print("✓ Data file accessibility check passed")
  
  # Step 2: Generate candidates
  controller <- run_step_candidates(controller)
  
  # QC Check 2: Candidate list size
  n_candidates <- nrow(controller$results$candidate_list)
  if (n_candidates < 100) {
    warning("Very few candidates generated. Check database and parameters.")
  } else if (n_candidates > 100000) {
    warning("Very large candidate list. Consider restricting mass range or modifications.")
  }
  print(paste("✓ Candidate list size check:", n_candidates, "candidates"))
  
  # Step 3: Preprocess
  controller <- run_step_preprocess(controller)
  
  # QC Check 3: Segmentation quality
  seg_labels <- controller$results$segmented_data[[1]]$segmentation_label
  if (length(seg_labels) < 2) {
    warning("Poor segmentation - only", length(seg_labels), "regions found")
  }
  
  # Check region sizes
  region_sizes <- sapply(seg_labels, length)
  if (any(region_sizes < 10)) {
    warning("Some regions have very few pixels")
  }
  print(paste("✓ Segmentation quality check:", length(seg_labels), "regions"))
  
  # Step 4: PMF search
  controller <- run_step_pmf(controller)
  
  # QC Check 4: Identification rates
  file_result <- controller$results$file_results[[1]]
  total_peptides <- nrow(file_result$combined_peptides)
  
  if (total_peptides == 0) {
    stop("No peptides identified! Check parameters.")
  } else if (total_peptides < 10) {
    warning("Very few peptides identified. Consider relaxing parameters.")
  }
  print(paste("✓ Identification rate check:", total_peptides, "peptides"))
  
  # Step 5: Summarize
  controller <- run_step_summarize(controller)
  
  print("✓ All QC checks completed successfully!")
  
  return(controller)
}

# =============================================================================
# Usage Notes and Best Practices
# =============================================================================

#' Best practices for using the modular workflow:
#' 
#' 1. Start with hitmap_quick_start() for initial exploration
#' 2. Use hitmap_step_by_step() when you need to inspect intermediate results
#' 3. Use hitmap_advanced_workflow() for production runs with specific parameters
#' 4. Use individual step functions when you need maximum control
#' 
#' Parameter guidelines:
#' - ppm: 2-3 for high-resolution, 5-10 for standard resolution
#' - threshold: 0.001-0.01 depending on data quality and noise
#' - segmentation_num: 4-8 for most applications, more for detailed spatial analysis
#' - FDR_cutoff: 0.01-0.05 for strict analysis, 0.1 for discovery
#' - peptide_ID_filter: 2-3 for reliable protein identification
#' 
#' Troubleshooting:
#' - If no candidates: Check database path and modification parameters
#' - If poor segmentation: Adjust segmentation_num or try different method
#' - If no identifications: Lower threshold or increase ppm tolerance
#' - If too many false positives: Decrease FDR_cutoff or increase peptide_ID_filter