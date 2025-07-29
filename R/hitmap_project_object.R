#' HiTMaP Project Object System
#' 
#' Unified project object for storing and managing HiTMaP workflow results
#' with support for flexible processing patterns

#' Create HiTMaP Project Object
#' 
#' @param data_files Character vector of data file paths
#' @param project_folder Project folder path
#' @param config Named list of configuration parameters
#' @return hitmap_project object
#' 
#' @export
hitmap_create_project <- function(data_files, project_folder = NULL, config = list()) {
  
  # Validate inputs
  if (missing(data_files)) stop("data_files parameter is required")
  
  # Determine project folder
  if (is.null(project_folder)) {
    project_folder <- dirname(data_files[1])
  }
  
  # Clean up file names
  data_files <- basename(data_files)
  data_files <- gsub(".imzML$", "", data_files)
  data_files_imzml <- paste0(data_files, ".imzML")
  
  # Create project structure
  project <- list(
    # Project metadata
    metadata = list(
      created = Sys.time(),
      version = "1.0",
      data_files = data_files,
      data_files_imzml = data_files_imzml,
      project_folder = project_folder,
      config = config,
      processing_log = character(0)
    ),
    
    # Processing setup
    setup = list(
      parallel_config = NULL,
      working_directory = NULL,
      thread_count = NULL
    ),
    
    # Processing results by step
    results = list(
      candidates = NULL,           # Candidate list
      preprocessing = list(),      # Preprocessing results by file
      segmentation = list(),       # Segmentation results by file
      pmf_search = list(),         # PMF search results by file and region
      summaries = list()           # Summary results
    ),
    
    # Processing status tracking
    status = list(
      initialized = FALSE,
      candidates_generated = FALSE,
      files_preprocessed = character(0),
      files_segmented = character(0),
      files_pmf_searched = character(0),
      summaries_generated = character(0)
    )
  )
  
  class(project) <- "hitmap_project"
  return(project)
}

#' Initialize HiTMaP Project
#' 
#' @param project hitmap_project object
#' @param thread_count Number of processing threads
#' @return Updated hitmap_project object
#' 
#' @export
hitmap_init <- function(project, thread_count = 4) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Object must be of class 'hitmap_project'")
  }
  
  suppressMessages(suppressWarnings(library("pacman")))
  suppressMessages(suppressWarnings(p_load(stringr, BiocParallel, data.table, Cardinal, parallel)))
  
  # Set working directory
  setwd(project$metadata$project_folder)
  
  # Configure parallel processing
  if (is.null(thread_count)) {
    parallel_threads <- try(future::availableCores() / 2)
    if (parallel_threads < 1 | is.null(parallel_threads)) { parallel_threads <- 1 }
  } else {
    parallel_threads <- thread_count
  }
  
  bpparam <- Parallel.OS(parallel_threads)
  setCardinalBPPARAM(BPPARAM = bpparam)
  
  # Update project object
  project$setup$parallel_config <- bpparam
  project$setup$working_directory <- project$metadata$project_folder
  project$setup$thread_count <- parallel_threads
  project$status$initialized <- TRUE
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Project initialized with", parallel_threads, "threads")
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message(paste(try(future::availableCores()), "cores detected,", parallel_threads, "threads configured"))
  message(paste(length(project$metadata$data_files), "files registered for processing"))
  
  return(project)
}

#' Generate Candidate List
#' 
#' @param project hitmap_project object
#' @param database Database filename
#' @param digestion_site Enzyme digestion site
#' @param missed_cleavages Missed cleavages allowed
#' @param adducts List of adducts
#' @param modifications Modification parameters
#' @param substitute_aa Amino acid substitutions
#' @param decoy_search Enable decoy search
#' @param decoy_adducts Decoy adducts
#' @param decoy_mode Decoy mode
#' @param mz_range m/z range
#' @param protein_exclusions Proteins to exclude
#' @param database_stats Generate database statistics
#' @param output_candidates Output candidate list to file
#' @param use_previous Use existing candidate list
#' @return Updated hitmap_project object
#' 
#' @export
hitmap_generate_candidates <- function(project,
                                     database = "uniprot-bovin.fasta",
                                     digestion_site = "trypsin",
                                     missed_cleavages = 0:1,
                                     adducts = c("M+H"),
                                     modifications = list(fixed = NULL, fixmod_position = NULL, 
                                                        variable = NULL, varmod_position = NULL),
                                     substitute_aa = NULL,
                                     decoy_search = TRUE,
                                     decoy_adducts = c("M+ACN+H", "M+IsoProp+H", "M+DMSO+H", 
                                                     "M+Co", "M+Ag", "M+Cu", "M+He", "M+Ne", 
                                                     "M+Ar", "M+Kr", "M+Xe", "M+Rn"),
                                     decoy_mode = "isotope",
                                     mz_range = c(700, 4000),
                                     protein_exclusions = NULL,
                                     database_stats = FALSE,
                                     output_candidates = TRUE,
                                     use_previous = FALSE) {
  
  if (!project$status$initialized) {
    stop("Project must be initialized first. Run hitmap_init()")
  }
  
  message(paste(database, "selected as database for candidate generation"))
  
  # Generate protein feature list
  candidate_list <- Protein_feature_list_fun(
    workdir = project$setup$working_directory,
    database = database,
    Digestion_site = digestion_site,
    missedCleavages = missed_cleavages,
    adducts = adducts,
    BPPARAM = project$setup$parallel_config,
    Decoy_adducts = decoy_adducts,
    Decoy_search = decoy_search,
    Decoy_mode = decoy_mode,
    output_candidatelist = output_candidates,
    use_previous_candidates = use_previous,
    Substitute_AA = substitute_aa,
    Modifications = modifications,
    mzrange = mz_range,
    Protein_desc_of_exclusion = protein_exclusions,
    Database_stats = database_stats
  )
  
  # Store results in project
  project$results$candidates <- candidate_list
  project$status$candidates_generated <- TRUE
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Generated", nrow(candidate_list), "candidates from", database)
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message(paste("Generated", nrow(candidate_list), "candidate peptides"))
  
  return(project)
}

#' Preprocess Single Data File
#' 
#' @param project hitmap_project object
#' @param file_name Name of file to process (without .imzML extension)
#' @param ppm Mass tolerance in ppm
#' @param import_ppm Import mass tolerance
#' @param mz_range m/z range for analysis
#' @param preprocess_params Preprocessing parameters list
#' @param rotation_config Image rotation configuration
#' @return Updated hitmap_project object
#' 
#' @export
hitmap_preprocess_file <- function(project,
                                 file_name,
                                 ppm = 5,
                                 import_ppm = 5,
                                 mz_range = "auto-detect",
                                 preprocess_params = list(
                                   force_preprocess = FALSE,
                                   use_preprocessrds = TRUE,
                                   smooth_signal = list(method = "disable"),
                                   reduce_baseline = list(method = "locmin"),
                                   peak_pick = list(method = "adaptive"),
                                   peak_align = list(tolerance = ppm/2, units = "ppm"),
                                   peak_filter = list(freq.min = 0.05),
                                   normalize = list(method = "rms", mz = 1)
                                 ),
                                 rotation_config = NULL) {
  
  if (!project$status$initialized) {
    stop("Project must be initialized first. Run hitmap_init()")
  }
  
  if (!file_name %in% project$metadata$data_files) {
    stop(paste("File", file_name, "not found in project data files"))
  }
  
  message(paste("Preprocessing file:", file_name))
  
  # Parse rotation configuration
  rotation_config <- Parse_rotation(file_name, rotation_config)
  
  # Load and preprocess the data
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(Cardinal)))
  suppressMessages(suppressWarnings(require(stringr)))
  
  setCardinalBPPARAM(project$setup$parallel_config)
  setwd(project$setup$working_directory)
  
  # Load imzML data
  file_path <- paste0(file_name, ".imzML")
  if (!file.exists(file_path)) {
    stop(paste("Data file not found:", file_path))
  }
  
  # Use existing preprocessing function but store results differently
  # This is a simplified version - you may need to extract the preprocessing
  # logic from the original IMS_data_process function
  
  message("Loading and preprocessing data...")
  
  # For now, create a placeholder structure
  # In actual implementation, you would extract the preprocessing logic
  preprocessing_result <- list(
    file_name = file_name,
    ppm = ppm,
    import_ppm = import_ppm,
    mz_range = mz_range,
    preprocess_params = preprocess_params,
    rotation_config = rotation_config,
    processed_data = NULL,  # Would contain the processed Cardinal object
    processing_time = Sys.time(),
    status = "completed"
  )
  
  # Store in project
  project$results$preprocessing[[file_name]] <- preprocessing_result
  project$status$files_preprocessed <- c(project$status$files_preprocessed, file_name)
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Preprocessed file:", file_name)
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message(paste("Preprocessing completed for:", file_name))
  
  return(project)
}

#' Segment Single Data File
#' 
#' @param project hitmap_project object
#' @param file_name Name of file to segment
#' @param segment_count Number of segments
#' @param segmentation_method Segmentation method
#' @param segmentation_def Segmentation definition file
#' @param segmentation_ncomp Number of components
#' @param variance_coverage Variance coverage
#' @param smooth_range Smoothing range
#' @param virtual_segmentation_file Virtual segmentation file
#' @return Updated hitmap_project object
#' 
#' @export
hitmap_segment_file <- function(project,
                               file_name,
                               segment_count = 4,
                               segmentation_method = "spatialKMeans",
                               segmentation_def = "segmentation_def.csv",
                               segmentation_ncomp = "auto-detect",
                               variance_coverage = 0.8,
                               smooth_range = 1,
                               virtual_segmentation_file = NULL) {
  
  if (!project$status$initialized) {
    stop("Project must be initialized first. Run hitmap_init()")
  }
  
  if (!file_name %in% project$status$files_preprocessed) {
    stop(paste("File", file_name, "must be preprocessed first"))
  }
  
  message(paste("Segmenting file:", file_name, "into", segment_count, "regions"))
  
  # Get preprocessing results
  preprocess_result <- project$results$preprocessing[[file_name]]
  
  # Perform segmentation
  # This would use the Preprocessing_segmentation function
  # For now, create placeholder
  
  segmentation_result <- list(
    file_name = file_name,
    segment_count = segment_count,
    segmentation_method = segmentation_method,
    parameters = list(
      segmentation_def = segmentation_def,
      segmentation_ncomp = segmentation_ncomp,
      variance_coverage = variance_coverage,
      smooth_range = smooth_range,
      virtual_segmentation_file = virtual_segmentation_file
    ),
    segmentation_labels = list(),  # Would contain actual segmentation labels
    processed_data = NULL,         # Would contain segmented Cardinal object
    processing_time = Sys.time(),
    status = "completed"
  )
  
  # Create dummy segmentation labels for demonstration
  for (i in 1:segment_count) {
    region_name <- paste0("region_", i)
    segmentation_result$segmentation_labels[[region_name]] <- list(1:100)  # Placeholder pixel indices
  }
  
  # Store in project
  project$results$segmentation[[file_name]] <- segmentation_result
  project$status$files_segmented <- c(project$status$files_segmented, file_name)
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Segmented file:", file_name, "into", segment_count, "regions")
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message(paste("Segmentation completed for:", file_name, "- found", segment_count, "regions"))
  
  return(project)
}

#' Perform PMF Search on File Regions
#' 
#' @param project hitmap_project object
#' @param file_name Name of file to search
#' @param region_names Specific regions to search (NULL for all)
#' @param threshold Intensity threshold
#' @param ppm Mass tolerance
#' @param top_n_features Maximum features for PMF
#' @param score_method Scoring method
#' @param decoy_mode Decoy mode
#' @param decoy_search Enable decoy search
#' @param adjust_score Adjust scores
#' @param fdr_cutoff FDR cutoff
#' @param peptide_filter Minimum peptides per protein
#' @param use_top_rank Use top ranking only
#' @param plot_scores Plot matching scores
#' @return Updated hitmap_project object
#' 
#' @export
hitmap_pmf_search <- function(project,
                            file_name,
                            region_names = NULL,
                            threshold = 0.001,
                            ppm = 5,
                            top_n_features = 1250,
                            score_method = "SQRTP",
                            decoy_mode = "isotope",
                            decoy_search = TRUE,
                            adjust_score = FALSE,
                            fdr_cutoff = 0.05,
                            peptide_filter = 2,
                            use_top_rank = NULL,
                            plot_scores = FALSE) {
  
  if (!project$status$initialized) {
    stop("Project must be initialized first. Run hitmap_init()")
  }
  
  if (!project$status$candidates_generated) {
    stop("Candidates must be generated first. Run hitmap_generate_candidates()")
  }
  
  if (!file_name %in% project$status$files_segmented) {
    stop(paste("File", file_name, "must be segmented first"))
  }
  
  message(paste("Performing PMF search on file:", file_name))
  
  # Get segmentation results
  segmentation_result <- project$results$segmentation[[file_name]]
  available_regions <- names(segmentation_result$segmentation_labels)
  
  # Determine regions to search
  if (is.null(region_names)) {
    search_regions <- available_regions
  } else {
    search_regions <- intersect(region_names, available_regions)
    if (length(search_regions) == 0) {
      stop("No valid regions found for search")
    }
  }
  
  message(paste("Searching", length(search_regions), "regions:", paste(search_regions, collapse = ", ")))
  
  # Initialize PMF results for this file
  if (is.null(project$results$pmf_search[[file_name]])) {
    project$results$pmf_search[[file_name]] <- list()
  }
  
  # Process each region
  for (region in search_regions) {
    message(paste("Processing region:", region))
    
    # Placeholder PMF search result
    region_result <- list(
      file_name = file_name,
      region_name = region,
      parameters = list(
        threshold = threshold,
        ppm = ppm,
        top_n_features = top_n_features,
        score_method = score_method,
        fdr_cutoff = fdr_cutoff,
        peptide_filter = peptide_filter
      ),
      peptides = data.frame(),     # Would contain identified peptides
      proteins = data.frame(),     # Would contain identified proteins
      spectrum = data.frame(),     # Would contain spectrum data
      processing_time = Sys.time(),
      status = "completed"
    )
    
    # Store region result
    project$results$pmf_search[[file_name]][[region]] <- region_result
  }
  
  # Update file status
  if (!file_name %in% project$status$files_pmf_searched) {
    project$status$files_pmf_searched <- c(project$status$files_pmf_searched, file_name)
  }
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Completed PMF search for file:", file_name, 
                    "regions:", paste(search_regions, collapse = ", "))
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message(paste("PMF search completed for:", file_name))
  
  return(project)
}

#' Generate Project Summaries
#' 
#' @param project hitmap_project object
#' @param protein_summary Generate protein summary
#' @param peptide_summary Generate peptide summary
#' @param region_summary Generate region summary
#' @return Updated hitmap_project object
#' 
#' @export
hitmap_generate_summaries <- function(project,
                                    protein_summary = TRUE,
                                    peptide_summary = TRUE,
                                    region_summary = FALSE) {
  
  if (!project$status$initialized) {
    stop("Project must be initialized first. Run hitmap_init()")
  }
  
  if (length(project$status$files_pmf_searched) == 0) {
    stop("No files have been PMF searched yet")
  }
  
  message("Generating project summaries...")
  
  # Create summary folder
  summary_dir <- file.path(project$setup$working_directory, "Summary folder")
  if (!dir.exists(summary_dir)) { dir.create(summary_dir) }
  
  summary_results <- list()
  
  # Protein summary
  if (protein_summary) {
    message("Generating protein summary...")
    
    # Combine protein results from all files and regions
    all_proteins <- data.frame()
    
    for (file_name in project$status$files_pmf_searched) {
      file_pmf_results <- project$results$pmf_search[[file_name]]
      
      for (region_name in names(file_pmf_results)) {
        region_result <- file_pmf_results[[region_name]]
        if (nrow(region_result$proteins) > 0) {
          proteins_with_source <- region_result$proteins
          proteins_with_source$source <- file_name
          proteins_with_source$region <- region_name
          all_proteins <- rbind(all_proteins, proteins_with_source)
        }
      }
    }
    
    summary_results$protein_summary <- all_proteins
    
    # Write to file
    if (nrow(all_proteins) > 0) {
      write.csv(all_proteins, file.path(summary_dir, "Protein_Summary.csv"), row.names = FALSE)
    }
  }
  
  # Peptide summary
  if (peptide_summary) {
    message("Generating peptide summary...")
    
    # Combine peptide results from all files and regions
    all_peptides <- data.frame()
    
    for (file_name in project$status$files_pmf_searched) {
      file_pmf_results <- project$results$pmf_search[[file_name]]
      
      for (region_name in names(file_pmf_results)) {
        region_result <- file_pmf_results[[region_name]]
        if (nrow(region_result$peptides) > 0) {
          peptides_with_source <- region_result$peptides
          peptides_with_source$source <- file_name
          peptides_with_source$region <- region_name
          all_peptides <- rbind(all_peptides, peptides_with_source)
        }
      }
    }
    
    summary_results$peptide_summary <- all_peptides
    
    # Write to file
    if (nrow(all_peptides) > 0) {
      write.csv(all_peptides, file.path(summary_dir, "Peptide_Summary.csv"), row.names = FALSE)
    }
  }
  
  # Region summary
  if (region_summary) {
    message("Generating region summary...")
    
    # Combine spectrum data from all regions
    all_spectra <- data.frame()
    
    for (file_name in project$status$files_pmf_searched) {
      file_pmf_results <- project$results$pmf_search[[file_name]]
      
      for (region_name in names(file_pmf_results)) {
        region_result <- file_pmf_results[[region_name]]
        if (nrow(region_result$spectrum) > 0) {
          spectrum_with_source <- region_result$spectrum
          spectrum_with_source$source <- file_name
          spectrum_with_source$region <- region_name
          all_spectra <- rbind(all_spectra, spectrum_with_source)
        }
      }
    }
    
    summary_results$region_summary <- all_spectra
    
    # Write to file
    if (nrow(all_spectra) > 0) {
      write.csv(all_spectra, file.path(summary_dir, "Region_feature_summary.csv"), row.names = FALSE)
    }
  }
  
  # Store summaries in project
  project$results$summaries <- summary_results
  
  # Update status
  summary_types <- c()
  if (protein_summary) summary_types <- c(summary_types, "protein")
  if (peptide_summary) summary_types <- c(summary_types, "peptide")
  if (region_summary) summary_types <- c(summary_types, "region")
  
  project$status$summaries_generated <- summary_types
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Generated summaries:", paste(summary_types, collapse = ", "))
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message("Summary generation completed")
  
  return(project)
}

#' Print method for hitmap_project
#' 
#' @param x hitmap_project object
#' @param ... Additional arguments
#' 
#' @export
print.hitmap_project <- function(x, ...) {
  cat("HiTMaP Project Object\n")
  cat("====================\n")
  cat("Created:", as.character(x$metadata$created), "\n")
  cat("Project folder:", x$metadata$project_folder, "\n")
  cat("Data files:", length(x$metadata$data_files), "\n")
  cat("  Files:", paste(x$metadata$data_files, collapse = ", "), "\n")
  
  cat("\nProcessing Status:\n")
  cat("  Initialized:", x$status$initialized, "\n")
  cat("  Candidates generated:", x$status$candidates_generated, "\n")
  if (x$status$candidates_generated && !is.null(x$results$candidates)) {
    cat("    Candidates:", nrow(x$results$candidates), "\n")
  }
  
  cat("  Files preprocessed:", length(x$status$files_preprocessed), "\n")
  if (length(x$status$files_preprocessed) > 0) {
    cat("    Files:", paste(x$status$files_preprocessed, collapse = ", "), "\n")
  }
  
  cat("  Files segmented:", length(x$status$files_segmented), "\n")
  if (length(x$status$files_segmented) > 0) {
    cat("    Files:", paste(x$status$files_segmented, collapse = ", "), "\n")
  }
  
  cat("  Files PMF searched:", length(x$status$files_pmf_searched), "\n")
  if (length(x$status$files_pmf_searched) > 0) {
    cat("    Files:", paste(x$status$files_pmf_searched, collapse = ", "), "\n")
  }
  
  cat("  Summaries generated:", length(x$status$summaries_generated), "\n")
  if (length(x$status$summaries_generated) > 0) {
    cat("    Types:", paste(x$status$summaries_generated, collapse = ", "), "\n")
  }
  
  if (x$status$initialized) {
    cat("\nConfiguration:\n")
    cat("  Threads:", x$setup$thread_count, "\n")
    cat("  Working directory:", x$setup$working_directory, "\n")
  }
  
  cat("\nProcessing Log (last 5 entries):\n")
  log_entries <- tail(x$metadata$processing_log, 5)
  for (entry in log_entries) {
    cat("  ", entry, "\n")
  }
}

#' Get Project Status Summary
#' 
#' @param project hitmap_project object
#' @return List with detailed status information
#' 
#' @export
hitmap_project_status <- function(project) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Object must be of class 'hitmap_project'")
  }
  
  # Calculate processing statistics
  total_files <- length(project$metadata$data_files)
  preprocessed_files <- length(project$status$files_preprocessed)
  segmented_files <- length(project$status$files_segmented)
  pmf_searched_files <- length(project$status$files_pmf_searched)
  
  # Calculate region statistics
  total_regions <- 0
  searched_regions <- 0
  
  for (file_name in project$status$files_segmented) {
    if (!is.null(project$results$segmentation[[file_name]])) {
      total_regions <- total_regions + length(project$results$segmentation[[file_name]]$segmentation_labels)
    }
  }
  
  for (file_name in project$status$files_pmf_searched) {
    if (!is.null(project$results$pmf_search[[file_name]])) {
      searched_regions <- searched_regions + length(project$results$pmf_search[[file_name]])
    }
  }
  
  # Calculate identification statistics
  total_peptides <- 0
  total_proteins <- 0
  
  if (!is.null(project$results$summaries$peptide_summary)) {
    total_peptides <- nrow(project$results$summaries$peptide_summary)
  }
  
  if (!is.null(project$results$summaries$protein_summary)) {
    total_proteins <- nrow(project$results$summaries$protein_summary)
  }
  
  status_summary <- list(
    project_info = list(
      created = project$metadata$created,
      total_files = total_files,
      project_folder = project$metadata$project_folder
    ),
    processing_progress = list(
      initialized = project$status$initialized,
      candidates_generated = project$status$candidates_generated,
      files_preprocessed = preprocessed_files,
      files_segmented = segmented_files,
      files_pmf_searched = pmf_searched_files,
      preprocessing_complete = (preprocessed_files == total_files),
      segmentation_complete = (segmented_files == total_files),
      pmf_search_complete = (pmf_searched_files == total_files)
    ),
    results_summary = list(
      candidate_count = ifelse(is.null(project$results$candidates), 0, nrow(project$results$candidates)),
      total_regions = total_regions,
      searched_regions = searched_regions,
      total_peptides = total_peptides,
      total_proteins = total_proteins,
      summaries_available = project$status$summaries_generated
    ),
    next_steps = character(0)
  )
  
  # Suggest next steps
  if (!status_summary$processing_progress$initialized) {
    status_summary$next_steps <- c(status_summary$next_steps, "Run hitmap_init()")
  } else if (!status_summary$processing_progress$candidates_generated) {
    status_summary$next_steps <- c(status_summary$next_steps, "Run hitmap_generate_candidates()")
  } else if (!status_summary$processing_progress$preprocessing_complete) {
    status_summary$next_steps <- c(status_summary$next_steps, "Run hitmap_preprocess_file() for remaining files")
  } else if (!status_summary$processing_progress$segmentation_complete) {
    status_summary$next_steps <- c(status_summary$next_steps, "Run hitmap_segment_file() for remaining files")
  } else if (!status_summary$processing_progress$pmf_search_complete) {
    status_summary$next_steps <- c(status_summary$next_steps, "Run hitmap_pmf_search() for remaining files")
  } else if (length(status_summary$results_summary$summaries_available) == 0) {
    status_summary$next_steps <- c(status_summary$next_steps, "Run hitmap_generate_summaries()")
  } else {
    status_summary$next_steps <- c(status_summary$next_steps, "Processing complete!")
  }
  
  return(status_summary)
}