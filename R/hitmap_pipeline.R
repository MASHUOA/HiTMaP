#' HiTMaP Pipeline Methods
#' 
#' Pipeline-style methods for HiTMaP workflow using %>% operator
#' Similar to Cardinal and dplyr design patterns

# Import pipe operator
#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' Create HiTMaP Project (Pipeline Start)
#' 
#' @param data_files Character vector of data file paths
#' @param project_folder Project folder path (optional)
#' @param config Initial configuration list
#' @return hitmap_project object ready for pipeline
#' 
#' @export
hitmap_project <- function(data_files, project_folder = NULL, config = list()) {
  project <- hitmap_create_project(data_files, project_folder, config)
  return(project)
}

#' Initialize Project (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param thread_count Number of processing threads
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
init <- function(project, thread_count = 4, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object. Use hitmap_project() to create one.")
  }
  
  project <- hitmap_init(project, thread_count)
  return(project)
}

#' Generate Candidates (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param database Database filename
#' @param digestion_site Enzyme digestion site
#' @param missed_cleavages Missed cleavages allowed
#' @param adducts List of adducts
#' @param modifications Modification parameters
#' @param decoy_search Enable decoy search
#' @param mz_range m/z range
#' @param .validate Validate inputs (default TRUE)
#' @param ... Additional parameters
#' @return Updated hitmap_project object
#' 
#' @export
generate_candidates <- function(project,
                               database = "uniprot-bovin.fasta",
                               digestion_site = "trypsin",
                               missed_cleavages = 0:1,
                               adducts = c("M+H"),
                               modifications = list(fixed = NULL, fixmod_position = NULL,
                                                   variable = NULL, varmod_position = NULL),
                               decoy_search = TRUE,
                               mz_range = c(700, 4000),
                               .validate = TRUE,
                               ...) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (.validate && !project$status$initialized) {
    stop("Project must be initialized first. Use init() in pipeline.")
  }
  
  project <- hitmap_generate_candidates(
    project = project,
    database = database,
    digestion_site = digestion_site,
    missed_cleavages = missed_cleavages,
    adducts = adducts,
    modifications = modifications,
    decoy_search = decoy_search,
    mz_range = mz_range,
    ...
  )
  
  return(project)
}

#' Preprocess Files (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param file_names Specific files to preprocess (NULL for all)
#' @param ppm Mass tolerance in ppm
#' @param import_ppm Import mass tolerance
#' @param mz_range m/z range
#' @param preprocess_profile Preprocessing profile
#' @param preprocess_params Custom preprocessing parameters
#' @param .validate Validate inputs (default TRUE)
#' @param ... Additional parameters
#' @return Updated hitmap_project object
#' 
#' @export
preprocess <- function(project,
                      file_names = NULL,
                      ppm = 5,
                      import_ppm = 5,
                      mz_range = "auto-detect",
                      preprocess_profile = "standard",
                      preprocess_params = NULL,
                      .validate = TRUE,
                      ...) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (.validate && !project$status$initialized) {
    stop("Project must be initialized first. Use init() in pipeline.")
  }
  
  # Determine files to process
  if (is.null(file_names)) {
    files_to_process <- project$metadata$data_files
  } else {
    files_to_process <- intersect(file_names, project$metadata$data_files)
    if (length(files_to_process) == 0) {
      stop("No valid files found to preprocess")
    }
  }
  
  # Create preprocessing parameters if not provided
  if (is.null(preprocess_params)) {
    preprocess_params <- hitmap_create_preprocess_config(preprocess_profile, ppm)
  }
  
  # Process each file
  for (file_name in files_to_process) {
    project <- hitmap_preprocess_file(
      project = project,
      file_name = file_name,
      ppm = ppm,
      import_ppm = import_ppm,
      mz_range = mz_range,
      preprocess_params = preprocess_params,
      ...
    )
  }
  
  return(project)
}

#' Segment Files (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param file_names Specific files to segment (NULL for all preprocessed)
#' @param segment_count Number of segments
#' @param segmentation_method Segmentation method
#' @param segmentation_params Additional segmentation parameters
#' @param .validate Validate inputs (default TRUE)
#' @param ... Additional parameters
#' @return Updated hitmap_project object
#' 
#' @export
segment <- function(project,
                   file_names = NULL,
                   segment_count = 4,
                   segmentation_method = "spatialKMeans",
                   segmentation_params = list(),
                   .validate = TRUE,
                   ...) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Determine files to process
  if (is.null(file_names)) {
    files_to_process <- project$status$files_preprocessed
  } else {
    files_to_process <- intersect(file_names, project$status$files_preprocessed)
    if (length(files_to_process) == 0) {
      stop("No valid preprocessed files found to segment")
    }
  }
  
  if (.validate && length(files_to_process) == 0) {
    stop("No preprocessed files available. Use preprocess() in pipeline first.")
  }
  
  # Set default parameters
  default_params <- list(
    segmentation_def = "segmentation_def.csv",
    segmentation_ncomp = "auto-detect",
    variance_coverage = 0.8,
    smooth_range = 1,
    virtual_segmentation_file = NULL
  )
  
  final_params <- modifyList(default_params, segmentation_params)
  
  # Process each file
  for (file_name in files_to_process) {
    project <- hitmap_segment_file(
      project = project,
      file_name = file_name,
      segment_count = segment_count,
      segmentation_method = segmentation_method,
      segmentation_def = final_params$segmentation_def,
      segmentation_ncomp = final_params$segmentation_ncomp,
      variance_coverage = final_params$variance_coverage,
      smooth_range = final_params$smooth_range,
      virtual_segmentation_file = final_params$virtual_segmentation_file,
      ...
    )
  }
  
  return(project)
}

#' Perform PMF Search (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param file_names Specific files to search (NULL for all segmented)
#' @param region_filter Function or character vector to filter regions
#' @param threshold Intensity threshold
#' @param ppm Mass tolerance
#' @param score_method Scoring method
#' @param fdr_cutoff FDR cutoff
#' @param peptide_filter Minimum peptides per protein
#' @param search_params Additional search parameters
#' @param .validate Validate inputs (default TRUE)
#' @param ... Additional parameters
#' @return Updated hitmap_project object
#' 
#' @export
search_pmf <- function(project,
                      file_names = NULL,
                      region_filter = NULL,
                      threshold = 0.001,
                      ppm = 5,
                      score_method = "SQRTP",
                      fdr_cutoff = 0.05,
                      peptide_filter = 2,
                      search_params = list(),
                      .validate = TRUE,
                      ...) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (.validate && !project$status$candidates_generated) {
    stop("Candidates must be generated first. Use generate_candidates() in pipeline.")
  }
  
  # Determine files to process
  if (is.null(file_names)) {
    files_to_process <- project$status$files_segmented
  } else {
    files_to_process <- intersect(file_names, project$status$files_segmented)
    if (length(files_to_process) == 0) {
      stop("No valid segmented files found to search")
    }
  }
  
  if (.validate && length(files_to_process) == 0) {
    stop("No segmented files available. Use segment() in pipeline first.")
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
  
  final_params <- modifyList(default_params, search_params)
  
  # Process each file
  for (file_name in files_to_process) {
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
      plot_scores = final_params$plot_scores,
      ...
    )
  }
  
  return(project)
}

#' Generate Summaries (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param protein_summary Generate protein summary
#' @param peptide_summary Generate peptide summary
#' @param region_summary Generate region summary
#' @param .validate Validate inputs (default TRUE)
#' @param ... Additional parameters
#' @return Updated hitmap_project object
#' 
#' @export
summarize <- function(project,
                     protein_summary = TRUE,
                     peptide_summary = TRUE,
                     region_summary = FALSE,
                     .validate = TRUE,
                     ...) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (.validate && length(project$status$files_pmf_searched) == 0) {
    stop("No files have been PMF searched yet. Use search_pmf() in pipeline first.")
  }
  
  project <- hitmap_generate_summaries(
    project = project,
    protein_summary = protein_summary,
    peptide_summary = peptide_summary,
    region_summary = region_summary,
    ...
  )
  
  return(project)
}

#' Filter Project Data (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param filter_type Type of filtering: "files", "regions", "candidates"
#' @param filter_function Function to apply for filtering
#' @param .validate Validate inputs (default TRUE)
#' @return Filtered hitmap_project object
#' 
#' @export
filter_data <- function(project, filter_type, filter_function, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (filter_type == "files") {
    # Filter data files
    original_files <- project$metadata$data_files
    keep_files <- original_files[filter_function(original_files)]
    
    # Update metadata
    project$metadata$data_files <- keep_files
    project$metadata$data_files_imzml <- paste0(keep_files, ".imzML")
    
    # Update status
    project$status$files_preprocessed <- intersect(project$status$files_preprocessed, keep_files)
    project$status$files_segmented <- intersect(project$status$files_segmented, keep_files)
    project$status$files_pmf_searched <- intersect(project$status$files_pmf_searched, keep_files)
    
    # Update results
    project$results$preprocessing <- project$results$preprocessing[keep_files]
    project$results$segmentation <- project$results$segmentation[keep_files]
    project$results$pmf_search <- project$results$pmf_search[keep_files]
    
  } else if (filter_type == "candidates" && !is.null(project$results$candidates)) {
    # Filter candidates
    project$results$candidates <- project$results$candidates[filter_function(project$results$candidates), ]
    
  } else if (filter_type == "regions") {
    # Filter regions within segmentation results
    for (file_name in names(project$results$segmentation)) {
      seg_result <- project$results$segmentation[[file_name]]
      available_regions <- names(seg_result$segmentation_labels)
      keep_regions <- available_regions[filter_function(available_regions)]
      
      # Update segmentation labels
      project$results$segmentation[[file_name]]$segmentation_labels <- 
        seg_result$segmentation_labels[keep_regions]
      
      # Update PMF search results if they exist
      if (file_name %in% names(project$results$pmf_search)) {
        project$results$pmf_search[[file_name]] <- 
          project$results$pmf_search[[file_name]][keep_regions]
      }
    }
  }
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Filtered", filter_type, "using custom function")
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  return(project)
}

#' Select Specific Data (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param files Character vector of file names to select
#' @param regions Character vector of region names to select
#' @param .validate Validate inputs (default TRUE)
#' @return Filtered hitmap_project object
#' 
#' @export
select_data <- function(project, files = NULL, regions = NULL, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Select specific files
  if (!is.null(files)) {
    project <- filter_data(project, "files", function(x) x %in% files, .validate = FALSE)
  }
  
  # Select specific regions
  if (!is.null(regions)) {
    project <- filter_data(project, "regions", function(x) x %in% regions, .validate = FALSE)
  }
  
  return(project)
}

#' Mutate Project Configuration (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param ... Named configuration parameters to update
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
mutate_config <- function(project, ..., .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Get new configuration parameters
  new_config <- list(...)
  
  # Update project configuration
  project$metadata$config <- modifyList(project$metadata$config, new_config)
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Updated configuration:", 
                    paste(names(new_config), collapse = ", "))
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  return(project)
}

#' View Project Status (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param detailed Show detailed status information
#' @param .validate Validate inputs (default TRUE)
#' @return hitmap_project object (unchanged, for pipeline continuity)
#' 
#' @export
view_status <- function(project, detailed = FALSE, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (detailed) {
    status <- hitmap_project_status(project)
    print(status)
  } else {
    print(project)
  }
  
  return(project)
}

#' Checkpoint Project (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param file_path Path to save the project checkpoint
#' @param message Optional checkpoint message
#' @param .validate Validate inputs (default TRUE)
#' @return hitmap_project object (unchanged, for pipeline continuity)
#' 
#' @export
checkpoint <- function(project, file_path, message = NULL, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Save project
  saveRDS(project, file_path)
  
  # Add to processing log
  checkpoint_msg <- ifelse(is.null(message), 
                          paste("Checkpoint saved to", file_path),
                          paste("Checkpoint:", message, "- saved to", file_path))
  
  log_entry <- paste(Sys.time(), "-", checkpoint_msg)
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message(paste("Project checkpoint saved to:", file_path))
  
  return(project)
}

#' Export Pipeline Results (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param output_folder Output folder path
#' @param export_types Types of data to export
#' @param format Export format
#' @param .validate Validate inputs (default TRUE)
#' @return hitmap_project object (unchanged, for pipeline continuity)
#' 
#' @export
export_results <- function(project,
                          output_folder,
                          export_types = c("summaries", "candidates"),
                          format = "csv",
                          .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  exported_files <- hitmap_export_results(
    project = project,
    output_folder = output_folder,
    export_types = export_types,
    format = format
  )
  
  # Add to processing log
  log_entry <- paste(Sys.time(), "- Exported results to", output_folder)
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  return(project)
}

#' Pipeline Configuration Helpers
#' 
#' @param profile Configuration profile name
#' @param ... Additional parameters to override
#' @return Configuration list
#' 
#' @export
config_basic <- function(...) {
  config <- list(
    ppm = 5,
    threshold = 0.001,
    segment_count = 4,
    preprocess_profile = "standard",
    modifications = hitmap_create_modification_config("basic"),
    fdr_cutoff = 0.05,
    peptide_filter = 2
  )
  
  # Override with user parameters
  user_params <- list(...)
  if (length(user_params) > 0) {
    config <- modifyList(config, user_params)
  }
  
  return(config)
}

#' @rdname config_basic
#' @export
config_strict <- function(...) {
  config <- list(
    ppm = 3,
    threshold = 0.005,
    segment_count = 6,
    preprocess_profile = "comprehensive",
    modifications = hitmap_create_modification_config("basic"),
    fdr_cutoff = 0.01,
    peptide_filter = 3
  )
  
  user_params <- list(...)
  if (length(user_params) > 0) {
    config <- modifyList(config, user_params)
  }
  
  return(config)
}

#' @rdname config_basic
#' @export
config_permissive <- function(...) {
  config <- list(
    ppm = 10,
    threshold = 0.0005,
    segment_count = 4,
    preprocess_profile = "minimal",
    modifications = hitmap_create_modification_config("comprehensive"),
    fdr_cutoff = 0.1,
    peptide_filter = 1
  )
  
  user_params <- list(...)
  if (length(user_params) > 0) {
    config <- modifyList(config, user_params)
  }
  
  return(config)
}

#' Load Project for Pipeline (Pipeline Start Alternative)
#' 
#' @param file_path Path to saved project RDS file
#' @return hitmap_project object ready for pipeline
#' 
#' @export
load_project <- function(file_path) {
  
  if (!file.exists(file_path)) {
    stop("Project file not found:", file_path)
  }
  
  project <- readRDS(file_path)
  
  if (!inherits(project, "hitmap_project")) {
    stop("File does not contain a valid hitmap_project object")
  }
  
  message(paste("Loaded project created:", project$metadata$created))
  
  return(project)
}