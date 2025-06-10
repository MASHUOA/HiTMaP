#' Flexible HiTMaP Workflow Functions
#' 
#' Functions for individual processing steps that can be used independently
#' or as part of a complete workflow

#' Candidate Processing Only Workflow
#' 
#' @param data_files Character vector of data file paths
#' @param project_folder Project folder path
#' @param database Database filename
#' @param digestion_site Enzyme digestion site
#' @param missed_cleavages Missed cleavages allowed
#' @param adducts List of adducts
#' @param modifications Modification parameters
#' @param decoy_search Enable decoy search
#' @param mz_range m/z range
#' @param thread_count Number of threads
#' @param output_file Output filename for candidates
#' @return hitmap_project object with candidates only
#' 
#' @export
hitmap_candidates_only <- function(data_files,
                                 project_folder = NULL,
                                 database = "uniprot-bovin.fasta",
                                 digestion_site = "trypsin",
                                 missed_cleavages = 0:1,
                                 adducts = c("M+H"),
                                 modifications = list(fixed = NULL, fixmod_position = NULL,
                                                    variable = NULL, varmod_position = NULL),
                                 decoy_search = TRUE,
                                 mz_range = c(700, 4000),
                                 thread_count = 4,
                                 output_file = "candidate_list.csv") {
  
  message("=== HiTMaP Candidate Generation Only ===")
  
  # Create and initialize project
  project <- hitmap_create_project(data_files, project_folder)
  project <- hitmap_init(project, thread_count)
  
  # Generate candidates
  project <- hitmap_generate_candidates(
    project = project,
    database = database,
    digestion_site = digestion_site,
    missed_cleavages = missed_cleavages,
    adducts = adducts,
    modifications = modifications,
    decoy_search = decoy_search,
    mz_range = mz_range
  )
  
  # Export candidates if requested
  if (!is.null(output_file)) {
    output_path <- file.path(project$setup$working_directory, output_file)
    write.csv(project$results$candidates, output_path, row.names = FALSE)
    message(paste("Candidates exported to:", output_path))
  }
  
  message("Candidate generation completed!")
  
  return(project)
}

#' Preprocessing Only Workflow
#' 
#' @param data_files Character vector of data file paths
#' @param project_folder Project folder path
#' @param file_names Specific files to preprocess (NULL for all)
#' @param ppm Mass tolerance in ppm
#' @param import_ppm Import mass tolerance
#' @param mz_range m/z range
#' @param preprocess_profile Preprocessing profile: "minimal", "standard", "comprehensive"
#' @param custom_preprocess_params Custom preprocessing parameters
#' @param thread_count Number of threads
#' @return hitmap_project object with preprocessing results
#' 
#' @export
hitmap_preprocess_only <- function(data_files,
                                 project_folder = NULL,
                                 file_names = NULL,
                                 ppm = 5,
                                 import_ppm = 5,
                                 mz_range = "auto-detect",
                                 preprocess_profile = "standard",
                                 custom_preprocess_params = NULL,
                                 thread_count = 4) {
  
  message("=== HiTMaP Preprocessing Only ===")
  
  # Create and initialize project
  project <- hitmap_create_project(data_files, project_folder)
  project <- hitmap_init(project, thread_count)
  
  # Determine files to process
  if (is.null(file_names)) {
    files_to_process <- project$metadata$data_files
  } else {
    files_to_process <- intersect(file_names, project$metadata$data_files)
    if (length(files_to_process) == 0) {
      stop("No valid files found to preprocess")
    }
  }
  
  # Create preprocessing parameters
  if (is.null(custom_preprocess_params)) {
    preprocess_params <- hitmap_create_preprocess_config(preprocess_profile, ppm)
  } else {
    preprocess_params <- custom_preprocess_params
  }
  
  # Process each file
  for (file_name in files_to_process) {
    message(paste("Processing file:", file_name))
    
    project <- hitmap_preprocess_file(
      project = project,
      file_name = file_name,
      ppm = ppm,
      import_ppm = import_ppm,
      mz_range = mz_range,
      preprocess_params = preprocess_params
    )
  }
  
  message("Preprocessing completed!")
  
  return(project)
}

#' Segmentation Only Workflow
#' 
#' @param project hitmap_project object with preprocessing completed
#' @param file_names Specific files to segment (NULL for all preprocessed)
#' @param segment_count Number of segments
#' @param segmentation_method Segmentation method
#' @param segmentation_params Additional segmentation parameters
#' @return Updated hitmap_project object with segmentation results
#' 
#' @export
hitmap_segment_only <- function(project,
                               file_names = NULL,
                               segment_count = 4,
                               segmentation_method = "spatialKMeans",
                               segmentation_params = list()) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  message("=== HiTMaP Segmentation Only ===")
  
  # Determine files to process
  if (is.null(file_names)) {
    files_to_process <- project$status$files_preprocessed
  } else {
    files_to_process <- intersect(file_names, project$status$files_preprocessed)
    if (length(files_to_process) == 0) {
      stop("No valid preprocessed files found to segment")
    }
  }
  
  if (length(files_to_process) == 0) {
    stop("No preprocessed files available. Run preprocessing first.")
  }
  
  # Set default segmentation parameters
  default_params <- list(
    segmentation_def = "segmentation_def.csv",
    segmentation_ncomp = "auto-detect",
    variance_coverage = 0.8,
    smooth_range = 1,
    virtual_segmentation_file = NULL
  )
  
  # Merge with user parameters
  final_params <- modifyList(default_params, segmentation_params)
  
  # Process each file
  for (file_name in files_to_process) {
    message(paste("Segmenting file:", file_name))
    
    project <- hitmap_segment_file(
      project = project,
      file_name = file_name,
      segment_count = segment_count,
      segmentation_method = segmentation_method,
      segmentation_def = final_params$segmentation_def,
      segmentation_ncomp = final_params$segmentation_ncomp,
      variance_coverage = final_params$variance_coverage,
      smooth_range = final_params$smooth_range,
      virtual_segmentation_file = final_params$virtual_segmentation_file
    )
  }
  
  message("Segmentation completed!")
  
  return(project)
}

#' PMF Search Only Workflow
#' 
#' @param project hitmap_project object with candidates and segmentation completed
#' @param file_names Specific files to search (NULL for all segmented)
#' @param region_filter Function to filter regions (e.g., function(x) x %in% c("region_1", "region_2"))
#' @param threshold Intensity threshold
#' @param ppm Mass tolerance
#' @param score_method Scoring method
#' @param fdr_cutoff FDR cutoff
#' @param peptide_filter Minimum peptides per protein
#' @param search_params Additional search parameters
#' @return Updated hitmap_project object with PMF search results
#' 
#' @export
hitmap_pmf_only <- function(project,
                          file_names = NULL,
                          region_filter = NULL,
                          threshold = 0.001,
                          ppm = 5,
                          score_method = "SQRTP",
                          fdr_cutoff = 0.05,
                          peptide_filter = 2,
                          search_params = list()) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (!project$status$candidates_generated) {
    stop("Candidates must be generated first")
  }
  
  message("=== HiTMaP PMF Search Only ===")
  
  # Determine files to process
  if (is.null(file_names)) {
    files_to_process <- project$status$files_segmented
  } else {
    files_to_process <- intersect(file_names, project$status$files_segmented)
    if (length(files_to_process) == 0) {
      stop("No valid segmented files found to search")
    }
  }
  
  if (length(files_to_process) == 0) {
    stop("No segmented files available. Run segmentation first.")
  }
  
  # Set default search parameters
  default_params <- list(
    top_n_features = 1250,
    decoy_mode = "isotope",
    decoy_search = TRUE,
    adjust_score = FALSE,
    use_top_rank = NULL,
    plot_scores = FALSE
  )
  
  # Merge with user parameters
  final_params <- modifyList(default_params, search_params)
  
  # Process each file
  for (file_name in files_to_process) {
    message(paste("Searching file:", file_name))
    
    # Get available regions
    segmentation_result <- project$results$segmentation[[file_name]]
    available_regions <- names(segmentation_result$segmentation_labels)
    
    # Apply region filter if provided
    if (!is.null(region_filter)) {
      if (is.function(region_filter)) {
        regions_to_search <- available_regions[region_filter(available_regions)]
      } else if (is.character(region_filter)) {
        regions_to_search <- intersect(region_filter, available_regions)
      } else {
        stop("region_filter must be a function or character vector")
      }
    } else {
      regions_to_search <- available_regions
    }
    
    if (length(regions_to_search) == 0) {
      warning(paste("No regions to search for file:", file_name))
      next
    }
    
    project <- hitmap_pmf_search(
      project = project,
      file_name = file_name,
      region_names = regions_to_search,
      threshold = threshold,
      ppm = ppm,
      top_n_features = final_params$top_n_features,
      score_method = score_method,
      decoy_mode = final_params$decoy_mode,
      decoy_search = final_params$decoy_search,
      adjust_score = final_params$adjust_score,
      fdr_cutoff = fdr_cutoff,
      peptide_filter = peptide_filter,
      use_top_rank = final_params$use_top_rank,
      plot_scores = final_params$plot_scores
    )
  }
  
  message("PMF search completed!")
  
  return(project)
}

#' Complete Flexible Workflow
#' 
#' @param data_files Character vector of data file paths
#' @param project_folder Project folder path
#' @param steps Character vector of steps to perform: "candidates", "preprocess", "segment", "pmf", "summarize"
#' @param config Named list of configuration parameters
#' @return hitmap_project object with requested steps completed
#' 
#' @export
hitmap_flexible_workflow <- function(data_files,
                                   project_folder = NULL,
                                   steps = c("candidates", "preprocess", "segment", "pmf", "summarize"),
                                   config = list()) {
  
  message("=== HiTMaP Flexible Workflow ===")
  message(paste("Requested steps:", paste(steps, collapse = ", ")))
  
  # Validate steps
  valid_steps <- c("candidates", "preprocess", "segment", "pmf", "summarize")
  invalid_steps <- setdiff(steps, valid_steps)
  if (length(invalid_steps) > 0) {
    stop(paste("Invalid steps:", paste(invalid_steps, collapse = ", ")))
  }
  
  # Default configuration
  default_config <- list(
    # General
    thread_count = 4,
    
    # Candidates
    database = "uniprot-bovin.fasta",
    digestion_site = "trypsin",
    missed_cleavages = 0:1,
    adducts = c("M+H"),
    modifications = list(fixed = NULL, fixmod_position = NULL, variable = NULL, varmod_position = NULL),
    decoy_search = TRUE,
    mz_range = c(700, 4000),
    
    # Preprocessing
    ppm = 5,
    import_ppm = 5,
    preprocess_profile = "standard",
    
    # Segmentation
    segment_count = 4,
    segmentation_method = "spatialKMeans",
    
    # PMF search
    threshold = 0.001,
    score_method = "SQRTP",
    fdr_cutoff = 0.05,
    peptide_filter = 2,
    
    # Summary
    protein_summary = TRUE,
    peptide_summary = TRUE,
    region_summary = FALSE
  )
  
  # Merge configurations
  final_config <- modifyList(default_config, config)
  
  # Create and initialize project
  project <- hitmap_create_project(data_files, project_folder, final_config)
  project <- hitmap_init(project, final_config$thread_count)
  
  # Execute requested steps
  
  # Step 1: Candidates
  if ("candidates" %in% steps) {
    message("\n--- Step: Generating Candidates ---")
    project <- hitmap_generate_candidates(
      project = project,
      database = final_config$database,
      digestion_site = final_config$digestion_site,
      missed_cleavages = final_config$missed_cleavages,
      adducts = final_config$adducts,
      modifications = final_config$modifications,
      decoy_search = final_config$decoy_search,
      mz_range = final_config$mz_range
    )
  }
  
  # Step 2: Preprocessing
  if ("preprocess" %in% steps) {
    message("\n--- Step: Preprocessing ---")
    
    # Check if candidates are needed and available
    if (("pmf" %in% steps) && !project$status$candidates_generated) {
      stop("Candidates must be generated before PMF search. Add 'candidates' to steps or run separately.")
    }
    
    preprocess_params <- hitmap_create_preprocess_config(final_config$preprocess_profile, final_config$ppm)
    
    for (file_name in project$metadata$data_files) {
      project <- hitmap_preprocess_file(
        project = project,
        file_name = file_name,
        ppm = final_config$ppm,
        import_ppm = final_config$import_ppm,
        preprocess_params = preprocess_params
      )
    }
  }
  
  # Step 3: Segmentation
  if ("segment" %in% steps) {
    message("\n--- Step: Segmentation ---")
    
    if (length(project$status$files_preprocessed) == 0) {
      stop("No files have been preprocessed. Add 'preprocess' to steps or run separately.")
    }
    
    for (file_name in project$status$files_preprocessed) {
      project <- hitmap_segment_file(
        project = project,
        file_name = file_name,
        segment_count = final_config$segment_count,
        segmentation_method = final_config$segmentation_method
      )
    }
  }
  
  # Step 4: PMF Search
  if ("pmf" %in% steps) {
    message("\n--- Step: PMF Search ---")
    
    if (!project$status$candidates_generated) {
      stop("Candidates must be generated before PMF search. Add 'candidates' to steps or run separately.")
    }
    
    if (length(project$status$files_segmented) == 0) {
      stop("No files have been segmented. Add 'segment' to steps or run separately.")
    }
    
    for (file_name in project$status$files_segmented) {
      project <- hitmap_pmf_search(
        project = project,
        file_name = file_name,
        threshold = final_config$threshold,
        ppm = final_config$ppm,
        score_method = final_config$score_method,
        fdr_cutoff = final_config$fdr_cutoff,
        peptide_filter = final_config$peptide_filter
      )
    }
  }
  
  # Step 5: Summarize
  if ("summarize" %in% steps) {
    message("\n--- Step: Summarizing ---")
    
    if (length(project$status$files_pmf_searched) == 0) {
      stop("No files have been PMF searched. Add 'pmf' to steps or run separately.")
    }
    
    project <- hitmap_generate_summaries(
      project = project,
      protein_summary = final_config$protein_summary,
      peptide_summary = final_config$peptide_summary,
      region_summary = final_config$region_summary
    )
  }
  
  message("\n=== Flexible Workflow Completed ===")
  print(project)
  
  return(project)
}

#' Resume Workflow from Project
#' 
#' @param project hitmap_project object
#' @param additional_steps Steps to add to the workflow
#' @param config Updated configuration parameters
#' @return Updated hitmap_project object
#' 
#' @export
hitmap_resume_workflow <- function(project, 
                                 additional_steps = c("summarize"),
                                 config = list()) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  message("=== Resuming HiTMaP Workflow ===")
  
  # Get current status
  status <- hitmap_project_status(project)
  
  # Update configuration if provided
  if (length(config) > 0) {
    project$metadata$config <- modifyList(project$metadata$config, config)
  }
  
  # Execute additional steps
  for (step in additional_steps) {
    
    if (step == "candidates" && !status$processing_progress$candidates_generated) {
      message("--- Generating Candidates ---")
      
      config_params <- project$metadata$config
      project <- hitmap_generate_candidates(
        project = project,
        database = config_params$database %||% "uniprot-bovin.fasta",
        digestion_site = config_params$digestion_site %||% "trypsin",
        missed_cleavages = config_params$missed_cleavages %||% 0:1,
        adducts = config_params$adducts %||% c("M+H"),
        modifications = config_params$modifications %||% list(),
        decoy_search = config_params$decoy_search %||% TRUE,
        mz_range = config_params$mz_range %||% c(700, 4000)
      )
      
    } else if (step == "preprocess" && !status$processing_progress$preprocessing_complete) {
      message("--- Preprocessing Remaining Files ---")
      
      files_to_process <- setdiff(project$metadata$data_files, project$status$files_preprocessed)
      config_params <- project$metadata$config
      
      preprocess_params <- hitmap_create_preprocess_config(
        config_params$preprocess_profile %||% "standard", 
        config_params$ppm %||% 5
      )
      
      for (file_name in files_to_process) {
        project <- hitmap_preprocess_file(
          project = project,
          file_name = file_name,
          ppm = config_params$ppm %||% 5,
          import_ppm = config_params$import_ppm %||% 5,
          preprocess_params = preprocess_params
        )
      }
      
    } else if (step == "segment" && !status$processing_progress$segmentation_complete) {
      message("--- Segmenting Remaining Files ---")
      
      files_to_process <- setdiff(project$status$files_preprocessed, project$status$files_segmented)
      config_params <- project$metadata$config
      
      for (file_name in files_to_process) {
        project <- hitmap_segment_file(
          project = project,
          file_name = file_name,
          segment_count = config_params$segment_count %||% 4,
          segmentation_method = config_params$segmentation_method %||% "spatialKMeans"
        )
      }
      
    } else if (step == "pmf" && !status$processing_progress$pmf_search_complete) {
      message("--- PMF Search on Remaining Files ---")
      
      if (!status$processing_progress$candidates_generated) {
        stop("Candidates must be generated before PMF search")
      }
      
      files_to_process <- setdiff(project$status$files_segmented, project$status$files_pmf_searched)
      config_params <- project$metadata$config
      
      for (file_name in files_to_process) {
        project <- hitmap_pmf_search(
          project = project,
          file_name = file_name,
          threshold = config_params$threshold %||% 0.001,
          ppm = config_params$ppm %||% 5,
          score_method = config_params$score_method %||% "SQRTP",
          fdr_cutoff = config_params$fdr_cutoff %||% 0.05,
          peptide_filter = config_params$peptide_filter %||% 2
        )
      }
      
    } else if (step == "summarize" && length(status$results_summary$summaries_available) == 0) {
      message("--- Generating Summaries ---")
      
      config_params <- project$metadata$config
      project <- hitmap_generate_summaries(
        project = project,
        protein_summary = config_params$protein_summary %||% TRUE,
        peptide_summary = config_params$peptide_summary %||% TRUE,
        region_summary = config_params$region_summary %||% FALSE
      )
    }
  }
  
  message("=== Resume Workflow Completed ===")
  print(project)
  
  return(project)
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Create Preprocessing Configuration
#' 
#' @param profile Preprocessing profile
#' @param ppm Mass tolerance
#' @return List of preprocessing parameters
#' 
#' @export
hitmap_create_preprocess_config <- function(profile = "standard", ppm = 5) {
  
  base_params <- list(
    force_preprocess = FALSE,
    use_preprocessrds = TRUE
  )
  
  if (profile == "minimal") {
    params <- list(
      smooth_signal = list(method = "disable"),
      reduce_baseline = list(method = "disable"),
      peak_pick = list(method = "adaptive"),
      peak_align = list(tolerance = ppm/2, units = "ppm"),
      peak_filter = list(freq.min = 0.01),
      normalize = list(method = "rms", mz = 1)
    )
  } else if (profile == "comprehensive") {
    params <- list(
      smooth_signal = list(method = "gaussian", sd = 1),
      reduce_baseline = list(method = "locmin", blocks = 100),
      peak_pick = list(method = "adaptive", snr = 3),
      peak_align = list(tolerance = ppm/2, units = "ppm"),
      peak_filter = list(freq.min = 0.1),
      normalize = list(method = "tic", mz = 1)
    )
  } else { # standard
    params <- list(
      smooth_signal = list(method = "disable"),
      reduce_baseline = list(method = "locmin"),
      peak_pick = list(method = "adaptive"),
      peak_align = list(tolerance = ppm/2, units = "ppm"),
      peak_filter = list(freq.min = 0.05),
      normalize = list(method = "rms", mz = 1)
    )
  }
  
  return(c(base_params, params))
}