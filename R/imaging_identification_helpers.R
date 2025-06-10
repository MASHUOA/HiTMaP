#' Helper Functions for Modular Imaging Identification Workflow
#' 
#' This file contains convenience functions and parameter configuration helpers
#' for the modular imaging identification workflow.

#' Create Default Preprocessing Parameters
#' 
#' @param method_profile Preprocessing profile: "minimal", "standard", "comprehensive"
#' @param ppm Mass tolerance for peak alignment
#' @param force_preprocess Force preprocessing even if already processed
#' @param use_preprocessRDS Use existing preprocessed RDS files
#' @return List of preprocessing parameters
#' 
#' @export
hitmap_create_preprocess_params <- function(method_profile = "standard",
                                          ppm = 5,
                                          force_preprocess = FALSE,
                                          use_preprocessRDS = TRUE) {
  
  base_params <- list(
    force_preprocess = force_preprocess,
    use_preprocessRDS = use_preprocessRDS
  )
  
  if (method_profile == "minimal") {
    params <- list(
      smoothSignal = list(method = "disable"),
      reduceBaseline = list(method = "disable"),
      peakPick = list(method = "adaptive"),
      peakAlign = list(tolerance = ppm/2, units = "ppm"),
      peakFilter = list(freq.min = 0.01),
      normalize = list(method = "rms", mz = 1)
    )
  } else if (method_profile == "comprehensive") {
    params <- list(
      smoothSignal = list(method = "gaussian", sd = 1),
      reduceBaseline = list(method = "locmin", blocks = 100),
      peakPick = list(method = "adaptive", SNR = 3),
      peakAlign = list(tolerance = ppm/2, units = "ppm"),
      peakFilter = list(freq.min = 0.1),
      normalize = list(method = "tic", mz = 1)
    )
  } else { # standard
    params <- list(
      smoothSignal = list(method = "disable"),
      reduceBaseline = list(method = "locmin"),
      peakPick = list(method = "adaptive"),
      peakAlign = list(tolerance = ppm/2, units = "ppm"),
      peakFilter = list(freq.min = 0.05),
      normalize = list(method = "rms", mz = 1)
    )
  }
  
  return(c(base_params, params))
}

#' Create Default Modification Parameters
#' 
#' @param modification_profile Modification profile: "none", "basic", "comprehensive"
#' @return List of modification parameters
#' 
#' @export
hitmap_create_modification_params <- function(modification_profile = "basic") {
  
  if (modification_profile == "none") {
    return(list(
      fixed = NULL,
      fixmod_position = NULL,
      variable = NULL,
      varmod_position = NULL
    ))
  } else if (modification_profile == "basic") {
    return(list(
      fixed = c("Carbamidomethyl (C)"),
      fixmod_position = c("any"),
      variable = c("Oxidation (M)"),
      varmod_position = c("any")
    ))
  } else if (modification_profile == "comprehensive") {
    return(list(
      fixed = c("Carbamidomethyl (C)"),
      fixmod_position = c("any"),
      variable = c("Oxidation (M)", "Acetyl (Protein N-term)", "Gln->pyro-Glu (N-term Q)"),
      varmod_position = c("any", "any", "any")
    ))
  }
}

#' Quick Start Function - Simple Workflow
#' 
#' @param datafile Path to data file(s)
#' @param database Database filename
#' @param projectfolder Optional project folder
#' @param ppm Mass tolerance in ppm
#' @param threshold Intensity threshold
#' @param segments Number of segments for spatial analysis
#' @param threads Number of processing threads
#' @return Complete workflow results
#' 
#' @export
hitmap_quick_start <- function(datafile,
                             database = "uniprot-bovin.fasta",
                             projectfolder = NULL,
                             ppm = 5,
                             threshold = 0.001,
                             segments = 4,
                             threads = 4) {
  
  message("Starting HiTMaP Quick Start Workflow...")
  message("=== Using default parameters for rapid analysis ===")
  
  # Use standard preprocessing
  preprocess_params <- hitmap_create_preprocess_params("standard", ppm)
  
  # Use basic modifications
  mod_params <- hitmap_create_modification_params("basic")
  
  # Run modular workflow
  results <- hitmap_modular_workflow(
    datafile = datafile,
    projectfolder = projectfolder,
    database = database,
    threshold = threshold,
    ppm = ppm,
    Thread = threads,
    segmentation_num = segments,
    Segmentation = "spatialKMeans",
    preprocess = preprocess_params,
    Modifications = mod_params,
    Digestion_site = "trypsin",
    missedCleavages = 0:1,
    adducts = c("M+H"),
    FDR_cutoff = 0.05,
    peptide_ID_filter = 2,
    score_method = "SQRTP"
  )
  
  return(results)
}

#' Advanced Configuration Function
#' 
#' @param datafile Path to data file(s)
#' @param config_list Named list of configuration parameters
#' @return Complete workflow results with custom configuration
#' 
#' @export
hitmap_advanced_workflow <- function(datafile, config_list = list()) {
  
  # Default configuration
  default_config <- list(
    database = "uniprot-bovin.fasta",
    projectfolder = NULL,
    threshold = 0.001,
    ppm = 5,
    Thread = 4,
    segmentation_num = 4,
    Segmentation = "spatialKMeans",
    preprocess_profile = "standard",
    modification_profile = "basic",
    Digestion_site = "trypsin",
    missedCleavages = 0:1,
    adducts = c("M+H"),
    FDR_cutoff = 0.05,
    peptide_ID_filter = 2,
    score_method = "SQRTP",
    Decoy_search = TRUE,
    Decoy_mode = "isotope",
    mzrange = c(700, 4000),
    Protein_feature_summary = TRUE,
    Peptide_feature_summary = TRUE,
    Region_feature_summary = FALSE
  )
  
  # Merge user config with defaults
  config <- modifyList(default_config, config_list)
  
  message("Starting HiTMaP Advanced Workflow with custom configuration...")
  
  # Create preprocessing parameters
  preprocess_params <- hitmap_create_preprocess_params(
    config$preprocess_profile, 
    config$ppm
  )
  
  # Create modification parameters
  mod_params <- hitmap_create_modification_params(config$modification_profile)
  
  # Run modular workflow
  results <- hitmap_modular_workflow(
    datafile = datafile,
    projectfolder = config$projectfolder,
    database = config$database,
    threshold = config$threshold,
    ppm = config$ppm,
    Thread = config$Thread,
    segmentation_num = config$segmentation_num,
    Segmentation = config$Segmentation,
    preprocess = preprocess_params,
    Modifications = mod_params,
    Digestion_site = config$Digestion_site,
    missedCleavages = config$missedCleavages,
    adducts = config$adducts,
    FDR_cutoff = config$FDR_cutoff,
    peptide_ID_filter = config$peptide_ID_filter,
    score_method = config$score_method,
    Decoy_search = config$Decoy_search,
    Decoy_mode = config$Decoy_mode,
    mzrange = config$mzrange,
    Protein_feature_summary = config$Protein_feature_summary,
    Peptide_feature_summary = config$Peptide_feature_summary,
    Region_feature_summary = config$Region_feature_summary
  )
  
  return(results)
}

#' Step-by-Step Interactive Workflow
#' 
#' Allows user to run each step individually and inspect results
#' 
#' @param datafile Path to data file(s)
#' @param config_list Configuration parameters
#' @return Workflow controller object
#' 
#' @export
hitmap_step_by_step <- function(datafile, config_list = list()) {
  
  # Create workflow controller
  controller <- list(
    datafile = datafile,
    config = config_list,
    step = 0,
    results = list(),
    completed_steps = character(0)
  )
  
  class(controller) <- "hitmap_controller"
  
  message("HiTMaP Step-by-Step Workflow Controller Created")
  message("Available steps:")
  message("1. run_step_init()     - Initialize project")
  message("2. run_step_candidates() - Generate candidate list")
  message("3. run_step_preprocess() - Preprocess and segment data")
  message("4. run_step_pmf()       - Perform PMF search")
  message("5. run_step_summarize() - Summarize results")
  message("Use controller$run_step_X() to execute each step")
  
  return(controller)
}

#' Run Step 1: Initialize Project
#' 
#' @param controller Workflow controller object
#' @param ... Additional parameters
#' 
#' @export
run_step_init <- function(controller, ...) {
  if (!"hitmap_controller" %in% class(controller)) {
    stop("Invalid controller object")
  }
  
  message("=== Step 1: Initializing Project ===")
  
  # Get parameters
  projectfolder <- controller$config$projectfolder
  threads <- ifelse(is.null(controller$config$Thread), 4, controller$config$Thread)
  
  # Run initialization
  project_setup <- hitmap_init_project(
    datafile = controller$datafile,
    projectfolder = projectfolder,
    Thread = threads
  )
  
  # Store results
  controller$results$project_setup <- project_setup
  controller$step <- 1
  controller$completed_steps <- c(controller$completed_steps, "init")
  
  message("Step 1 completed successfully!")
  message("Next: run_step_candidates(controller)")
  
  return(controller)
}

#' Run Step 2: Generate Candidates
#' 
#' @param controller Workflow controller object
#' @param ... Additional parameters
#' 
#' @export
run_step_candidates <- function(controller, ...) {
  if (!"init" %in% controller$completed_steps) {
    stop("Must complete initialization step first")
  }
  
  message("=== Step 2: Generating Candidate List ===")
  
  # Get parameters
  database <- ifelse(is.null(controller$config$database), "uniprot-bovin.fasta", controller$config$database)
  Digestion_site <- ifelse(is.null(controller$config$Digestion_site), "trypsin", controller$config$Digestion_site)
  missedCleavages <- if(is.null(controller$config$missedCleavages)) 0:1 else controller$config$missedCleavages
  adducts <- if(is.null(controller$config$adducts)) c("M+H") else controller$config$adducts
  
  # Generate modification parameters
  mod_profile <- ifelse(is.null(controller$config$modification_profile), "basic", controller$config$modification_profile)
  mod_params <- hitmap_create_modification_params(mod_profile)
  
  candidate_list <- hitmap_generate_candidates(
    workdir = controller$results$project_setup$workdir,
    database = database,
    Digestion_site = Digestion_site,
    missedCleavages = missedCleavages,
    adducts = adducts,
    Modifications = mod_params,
    BPPARAM = controller$results$project_setup$BPPARAM,
    ...
  )
  
  controller$results$candidate_list <- candidate_list
  controller$step <- 2
  controller$completed_steps <- c(controller$completed_steps, "candidates")
  
  message("Step 2 completed successfully!")
  message(paste("Generated", nrow(candidate_list), "candidate peptides"))
  message("Next: run_step_preprocess(controller)")
  
  return(controller)
}

#' Run Step 3: Preprocess and Segment
#' 
#' @param controller Workflow controller object
#' @param file_index Index of file to process (for multiple files)
#' @param ... Additional parameters
#' 
#' @export
run_step_preprocess <- function(controller, file_index = 1, ...) {
  if (!"candidates" %in% controller$completed_steps) {
    stop("Must complete candidate generation step first")
  }
  
  if (file_index > length(controller$results$project_setup$datafile)) {
    stop("File index out of range")
  }
  
  current_file <- controller$results$project_setup$datafile[file_index]
  message(paste("=== Step 3: Preprocessing and Segmentation for file:", current_file, "==="))
  
  # Get parameters
  ppm <- ifelse(is.null(controller$config$ppm), 5, controller$config$ppm)
  segmentation_num <- ifelse(is.null(controller$config$segmentation_num), 4, controller$config$segmentation_num)
  Segmentation <- ifelse(is.null(controller$config$Segmentation), "spatialKMeans", controller$config$Segmentation)
  
  # Create preprocessing parameters
  preprocess_profile <- ifelse(is.null(controller$config$preprocess_profile), "standard", controller$config$preprocess_profile)
  preprocess_params <- hitmap_create_preprocess_params(preprocess_profile, ppm)
  
  segmented_data <- hitmap_preprocess_segment(
    datafile = current_file,
    workdir = controller$results$project_setup$workdir,
    ppm = ppm,
    segmentation_num = segmentation_num,
    Segmentation = Segmentation,
    preprocess = preprocess_params,
    BPPARAM = controller$results$project_setup$BPPARAM,
    ...
  )
  
  # Store results
  if (is.null(controller$results$segmented_data)) {
    controller$results$segmented_data <- list()
  }
  controller$results$segmented_data[[current_file]] <- segmented_data
  
  controller$step <- 3
  if (!"preprocess" %in% controller$completed_steps) {
    controller$completed_steps <- c(controller$completed_steps, "preprocess")
  }
  
  message("Step 3 completed successfully!")
  message(paste("Found", length(segmented_data$segmentation_label), "regions:", paste(names(segmented_data$segmentation_label), collapse = ", ")))
  message("Next: run_step_pmf(controller)")
  
  return(controller)
}

#' Run Step 4: PMF Search
#' 
#' @param controller Workflow controller object
#' @param file_index Index of file to process
#' @param ... Additional parameters
#' 
#' @export
run_step_pmf <- function(controller, file_index = 1, ...) {
  if (!"preprocess" %in% controller$completed_steps) {
    stop("Must complete preprocessing step first")
  }
  
  current_file <- controller$results$project_setup$datafile[file_index]
  
  if (is.null(controller$results$segmented_data[[current_file]])) {
    stop(paste("No segmented data found for file:", current_file, ". Run preprocessing first."))
  }
  
  message(paste("=== Step 4: PMF Search for file:", current_file, "==="))
  
  # Get parameters
  threshold <- ifelse(is.null(controller$config$threshold), 0.001, controller$config$threshold)
  ppm <- ifelse(is.null(controller$config$ppm), 5, controller$config$ppm)
  score_method <- ifelse(is.null(controller$config$score_method), "SQRTP", controller$config$score_method)
  FDR_cutoff <- ifelse(is.null(controller$config$FDR_cutoff), 0.05, controller$config$FDR_cutoff)
  peptide_ID_filter <- ifelse(is.null(controller$config$peptide_ID_filter), 2, controller$config$peptide_ID_filter)
  
  file_results <- hitmap_process_all_regions(
    segmented_data = controller$results$segmented_data[[current_file]],
    candidate_list = controller$results$candidate_list,
    threshold = threshold,
    ppm = ppm,
    score_method = score_method,
    FDR_cutoff = FDR_cutoff,
    peptide_ID_filter = peptide_ID_filter,
    BPPARAM = controller$results$project_setup$BPPARAM,
    ...
  )
  
  # Store results
  if (is.null(controller$results$file_results)) {
    controller$results$file_results <- list()
  }
  controller$results$file_results[[current_file]] <- file_results
  
  controller$step <- 4
  if (!"pmf" %in% controller$completed_steps) {
    controller$completed_steps <- c(controller$completed_steps, "pmf")
  }
  
  message("Step 4 completed successfully!")
  message(paste("Identified", nrow(file_results$combined_peptides), "peptides across all regions"))
  message("Next: run_step_summarize(controller)")
  
  return(controller)
}

#' Run Step 5: Summarize Results
#' 
#' @param controller Workflow controller object
#' @param ... Additional parameters
#' 
#' @export
run_step_summarize <- function(controller, ...) {
  if (!"pmf" %in% controller$completed_steps) {
    stop("Must complete PMF search step first")
  }
  
  message("=== Step 5: Summarizing Results ===")
  
  # Get parameters
  Protein_feature_summary <- ifelse(is.null(controller$config$Protein_feature_summary), TRUE, controller$config$Protein_feature_summary)
  Peptide_feature_summary <- ifelse(is.null(controller$config$Peptide_feature_summary), TRUE, controller$config$Peptide_feature_summary)
  Region_feature_summary <- ifelse(is.null(controller$config$Region_feature_summary), FALSE, controller$config$Region_feature_summary)
  
  summary_results <- hitmap_summarize_results(
    project_setup = controller$results$project_setup,
    Protein_feature_summary = Protein_feature_summary,
    Peptide_feature_summary = Peptide_feature_summary,
    Region_feature_summary = Region_feature_summary
  )
  
  controller$results$summary_results <- summary_results
  controller$step <- 5
  controller$completed_steps <- c(controller$completed_steps, "summarize")
  
  message("Step 5 completed successfully!")
  message("Workflow completed! Check controller$results for all outputs.")
  
  return(controller)
}

#' Print method for workflow controller
#' 
#' @param x Workflow controller object
#' @param ... Additional arguments
#' 
#' @export
print.hitmap_controller <- function(x, ...) {
  cat("HiTMaP Workflow Controller\n")
  cat("=========================\n")
  cat("Current step:", x$step, "\n")
  cat("Completed steps:", paste(x$completed_steps, collapse = ", "), "\n")
  cat("Data files:", length(x$datafile), "\n")
  
  if (x$step >= 1) {
    cat("Working directory:", x$results$project_setup$workdir, "\n")
  }
  
  if (x$step >= 2) {
    cat("Candidate peptides:", nrow(x$results$candidate_list), "\n")
  }
  
  if (x$step >= 3) {
    total_regions <- sum(sapply(x$results$segmented_data, function(x) length(x$segmentation_label)))
    cat("Total regions segmented:", total_regions, "\n")
  }
  
  if (x$step >= 4) {
    total_peptides <- sum(sapply(x$results$file_results, function(x) nrow(x$combined_peptides)))
    cat("Total peptides identified:", total_peptides, "\n")
  }
  
  cat("\nNext steps:\n")
  if (x$step == 0) cat("- run_step_init(controller)\n")
  if (x$step == 1) cat("- run_step_candidates(controller)\n")
  if (x$step == 2) cat("- run_step_preprocess(controller)\n")
  if (x$step == 3) cat("- run_step_pmf(controller)\n")
  if (x$step == 4) cat("- run_step_summarize(controller)\n")
  if (x$step >= 5) cat("- Workflow completed!\n")
}