#' HiTMaP Convenience Functions
#' 
#' Easy-to-use wrapper functions with consistent lowercase naming

#' Quick Start Function with Consistent Naming
#' 
#' @param data_files Path to data file(s)
#' @param database Database filename
#' @param project_folder Optional project folder
#' @param ppm Mass tolerance in ppm
#' @param threshold Intensity threshold
#' @param segment_count Number of segments for spatial analysis
#' @param thread_count Number of processing threads
#' @param output_summaries Generate summary files
#' @return Complete workflow results
#' 
#' @export
hitmap_quick_start <- function(data_files,
                             database = "uniprot-bovin.fasta",
                             project_folder = NULL,
                             ppm = 5,
                             threshold = 0.001,
                             segment_count = 4,
                             thread_count = 4,
                             output_summaries = TRUE) {
  
  message("=== HiTMaP Quick Start Workflow ===")
  message("Using default parameters for rapid analysis")
  
  # Default configuration
  config <- list(
    database = database,
    ppm = ppm,
    threshold = threshold,
    segment_count = segment_count,
    thread_count = thread_count,
    preprocess_profile = "standard",
    segmentation_method = "spatialKMeans",
    digestion_site = "trypsin",
    missed_cleavages = 0:1,
    adducts = c("M+H"),
    modifications = hitmap_create_modification_config("basic"),
    fdr_cutoff = 0.05,
    peptide_filter = 2,
    score_method = "SQRTP"
  )
  
  # Determine steps
  steps <- c("candidates", "preprocess", "segment", "pmf")
  if (output_summaries) {
    steps <- c(steps, "summarize")
  }
  
  # Run workflow
  results <- hitmap_flexible_workflow(
    data_files = data_files,
    project_folder = project_folder,
    steps = steps,
    config = config
  )
  
  return(results)
}

#' Advanced Workflow with Custom Configuration
#' 
#' @param data_files Path to data file(s)
#' @param config Named list of configuration parameters
#' @param steps Steps to execute
#' @return Complete workflow results
#' 
#' @export
hitmap_advanced <- function(data_files, config = list(), steps = c("candidates", "preprocess", "segment", "pmf", "summarize")) {
  
  message("=== HiTMaP Advanced Workflow ===")
  message("Using custom configuration")
  
  # Default configuration with lowercase naming
  default_config <- list(
    # General
    thread_count = 4,
    project_folder = NULL,
    
    # Database and candidates
    database = "uniprot-bovin.fasta",
    digestion_site = "trypsin",
    missed_cleavages = 0:1,
    adducts = c("M+H"),
    modifications = hitmap_create_modification_config("basic"),
    decoy_search = TRUE,
    decoy_mode = "isotope",
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
  
  # Merge user config with defaults
  final_config <- modifyList(default_config, config)
  
  # Run flexible workflow
  results <- hitmap_flexible_workflow(
    data_files = data_files,
    project_folder = final_config$project_folder,
    steps = steps,
    config = final_config
  )
  
  return(results)
}

#' Create Modification Configuration
#' 
#' @param profile Modification profile: "none", "basic", "comprehensive", "custom"
#' @param custom_modifications Custom modification list (used when profile = "custom")
#' @return List of modification parameters
#' 
#' @export
hitmap_create_modification_config <- function(profile = "basic", custom_modifications = NULL) {
  
  if (profile == "none") {
    return(list(
      fixed = NULL,
      fixmod_position = NULL,
      variable = NULL,
      varmod_position = NULL
    ))
  } else if (profile == "basic") {
    return(list(
      fixed = c("Carbamidomethyl (C)"),
      fixmod_position = c("any"),
      variable = c("Oxidation (M)"),
      varmod_position = c("any")
    ))
  } else if (profile == "comprehensive") {
    return(list(
      fixed = c("Carbamidomethyl (C)"),
      fixmod_position = c("any"),
      variable = c("Oxidation (M)", "Acetyl (Protein N-term)", "Gln->pyro-Glu (N-term Q)", "Deamidated (NQ)"),
      varmod_position = c("any", "any", "any", "any")
    ))
  } else if (profile == "custom") {
    if (is.null(custom_modifications)) {
      stop("custom_modifications must be provided when profile = 'custom'")
    }
    return(custom_modifications)
  } else {
    stop("Invalid modification profile. Use 'none', 'basic', 'comprehensive', or 'custom'")
  }
}

#' Create Search Configuration
#' 
#' @param sensitivity Search sensitivity: "strict", "standard", "permissive"
#' @param custom_params Custom parameter list
#' @return List of search parameters
#' 
#' @export
hitmap_create_search_config <- function(sensitivity = "standard", custom_params = NULL) {
  
  if (sensitivity == "strict") {
    config <- list(
      threshold = 0.005,
      ppm = 3,
      fdr_cutoff = 0.01,
      peptide_filter = 3,
      score_method = "SQRTP"
    )
  } else if (sensitivity == "standard") {
    config <- list(
      threshold = 0.001,
      ppm = 5,
      fdr_cutoff = 0.05,
      peptide_filter = 2,
      score_method = "SQRTP"
    )
  } else if (sensitivity == "permissive") {
    config <- list(
      threshold = 0.0005,
      ppm = 10,
      fdr_cutoff = 0.1,
      peptide_filter = 1,
      score_method = "SQRTP"
    )
  } else {
    stop("Invalid sensitivity. Use 'strict', 'standard', or 'permissive'")
  }
  
  # Merge with custom parameters if provided
  if (!is.null(custom_params)) {
    config <- modifyList(config, custom_params)
  }
  
  return(config)
}

#' Batch Process Multiple Datasets
#' 
#' @param dataset_list List of datasets, each containing data_files and optional config
#' @param base_config Base configuration applied to all datasets
#' @param parallel_datasets Process datasets in parallel
#' @param output_folder Base output folder
#' @return List of results for each dataset
#' 
#' @export
hitmap_batch_process <- function(dataset_list,
                               base_config = list(),
                               parallel_datasets = FALSE,
                               output_folder = NULL) {
  
  message("=== HiTMaP Batch Processing ===")
  message(paste("Processing", length(dataset_list), "datasets"))
  
  if (!is.null(output_folder) && !dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  # Function to process single dataset
  process_dataset <- function(dataset_info, dataset_name) {
    
    message(paste("Processing dataset:", dataset_name))
    
    # Extract dataset information
    data_files <- dataset_info$data_files
    dataset_config <- dataset_info$config %||% list()
    
    # Merge with base config
    final_config <- modifyList(base_config, dataset_config)
    
    # Set project folder if output folder specified
    if (!is.null(output_folder)) {
      final_config$project_folder <- file.path(output_folder, dataset_name)
      if (!dir.exists(final_config$project_folder)) {
        dir.create(final_config$project_folder, recursive = TRUE)
      }
    }
    
    # Process dataset
    tryCatch({
      result <- hitmap_advanced(
        data_files = data_files,
        config = final_config
      )
      
      # Add dataset name to result
      result$metadata$dataset_name <- dataset_name
      
      message(paste("Completed dataset:", dataset_name))
      return(result)
      
    }, error = function(e) {
      warning(paste("Failed to process dataset", dataset_name, ":", e$message))
      return(NULL)
    })
  }
  
  # Process datasets
  dataset_names <- names(dataset_list)
  if (is.null(dataset_names)) {
    dataset_names <- paste0("dataset_", seq_along(dataset_list))
  }
  
  if (parallel_datasets && length(dataset_list) > 1) {
    # Parallel processing
    suppressMessages(suppressWarnings(require(parallel)))
    
    results <- mclapply(seq_along(dataset_list), function(i) {
      process_dataset(dataset_list[[i]], dataset_names[i])
    }, mc.cores = min(length(dataset_list), detectCores() - 1))
    
    names(results) <- dataset_names
    
  } else {
    # Sequential processing
    results <- list()
    
    for (i in seq_along(dataset_list)) {
      results[[dataset_names[i]]] <- process_dataset(dataset_list[[i]], dataset_names[i])
    }
  }
  
  # Generate batch summary
  successful_datasets <- sum(!sapply(results, is.null))
  failed_datasets <- length(results) - successful_datasets
  
  message("=== Batch Processing Summary ===")
  message(paste("Successful:", successful_datasets))
  message(paste("Failed:", failed_datasets))
  
  # Create batch summary report
  if (!is.null(output_folder)) {
    batch_summary <- list(
      processing_time = Sys.time(),
      total_datasets = length(dataset_list),
      successful_datasets = successful_datasets,
      failed_datasets = failed_datasets,
      dataset_names = dataset_names,
      base_config = base_config
    )
    
    saveRDS(batch_summary, file.path(output_folder, "batch_summary.rds"))
    
    # Create summary CSV
    summary_df <- data.frame(
      dataset = dataset_names,
      status = ifelse(sapply(results, is.null), "failed", "success"),
      stringsAsFactors = FALSE
    )
    
    # Add result statistics
    for (i in seq_along(results)) {
      if (!is.null(results[[i]])) {
        status <- hitmap_project_status(results[[i]])
        summary_df$total_peptides[i] <- status$results_summary$total_peptides
        summary_df$total_proteins[i] <- status$results_summary$total_proteins
        summary_df$total_regions[i] <- status$results_summary$total_regions
      } else {
        summary_df$total_peptides[i] <- NA
        summary_df$total_proteins[i] <- NA
        summary_df$total_regions[i] <- NA
      }
    }
    
    write.csv(summary_df, file.path(output_folder, "batch_summary.csv"), row.names = FALSE)
  }
  
  return(results)
}

#' Compare Multiple Results
#' 
#' @param result_list List of hitmap_project objects to compare
#' @param comparison_type Type of comparison: "peptides", "proteins", "regions"
#' @param output_file Optional output file for comparison
#' @return Comparison data frame
#' 
#' @export
hitmap_compare_results <- function(result_list, 
                                 comparison_type = "peptides",
                                 output_file = NULL) {
  
  message("=== HiTMaP Results Comparison ===")
  
  # Validate inputs
  if (!all(sapply(result_list, function(x) inherits(x, "hitmap_project")))) {
    stop("All results must be hitmap_project objects")
  }
  
  # Get names for results
  result_names <- names(result_list)
  if (is.null(result_names)) {
    result_names <- paste0("result_", seq_along(result_list))
  }
  
  # Extract comparison data
  comparison_data <- list()
  
  for (i in seq_along(result_list)) {
    project <- result_list[[i]]
    name <- result_names[i]
    
    if (comparison_type == "peptides" && !is.null(project$results$summaries$peptide_summary)) {
      data <- project$results$summaries$peptide_summary
      data$result_name <- name
      comparison_data[[name]] <- data
      
    } else if (comparison_type == "proteins" && !is.null(project$results$summaries$protein_summary)) {
      data <- project$results$summaries$protein_summary
      data$result_name <- name
      comparison_data[[name]] <- data
      
    } else if (comparison_type == "regions" && !is.null(project$results$summaries$region_summary)) {
      data <- project$results$summaries$region_summary
      data$result_name <- name
      comparison_data[[name]] <- data
    }
  }
  
  if (length(comparison_data) == 0) {
    stop(paste("No", comparison_type, "data available for comparison"))
  }
  
  # Combine data
  combined_data <- do.call(rbind, comparison_data)
  
  # Generate comparison statistics
  stats_by_result <- aggregate(
    combined_data[, sapply(combined_data, is.numeric), drop = FALSE],
    by = list(result_name = combined_data$result_name),
    FUN = function(x) c(count = length(x), mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE))
  )
  
  # Create summary
  comparison_summary <- list(
    comparison_type = comparison_type,
    result_names = result_names,
    combined_data = combined_data,
    statistics = stats_by_result,
    comparison_time = Sys.time()
  )
  
  # Output to file if requested
  if (!is.null(output_file)) {
    write.csv(combined_data, paste0(output_file, "_combined_data.csv"), row.names = FALSE)
    write.csv(stats_by_result, paste0(output_file, "_statistics.csv"), row.names = FALSE)
    saveRDS(comparison_summary, paste0(output_file, "_comparison.rds"))
  }
  
  message(paste("Comparison completed for", length(result_list), "results"))
  message(paste("Total", comparison_type, "entries:", nrow(combined_data)))
  
  return(comparison_summary)
}

#' Export Project Results
#' 
#' @param project hitmap_project object
#' @param output_folder Output folder path
#' @param export_types Types of data to export: "summaries", "candidates", "raw_results"
#' @param format Export format: "csv", "excel", "rds"
#' @return List of exported file paths
#' 
#' @export
hitmap_export_results <- function(project,
                                output_folder,
                                export_types = c("summaries", "candidates"),
                                format = "csv") {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (!dir.exists(output_folder)) {
    dir.create(output_folder, recursive = TRUE)
  }
  
  message("=== Exporting HiTMaP Results ===")
  
  exported_files <- list()
  
  # Export summaries
  if ("summaries" %in% export_types && !is.null(project$results$summaries)) {
    
    summaries <- project$results$summaries
    
    for (summary_type in names(summaries)) {
      if (!is.null(summaries[[summary_type]]) && nrow(summaries[[summary_type]]) > 0) {
        
        filename <- paste0(summary_type, "_summary")
        
        if (format == "csv") {
          filepath <- file.path(output_folder, paste0(filename, ".csv"))
          write.csv(summaries[[summary_type]], filepath, row.names = FALSE)
          
        } else if (format == "excel") {
          suppressMessages(suppressWarnings(require(openxlsx)))
          filepath <- file.path(output_folder, paste0(filename, ".xlsx"))
          write.xlsx(summaries[[summary_type]], filepath)
          
        } else if (format == "rds") {
          filepath <- file.path(output_folder, paste0(filename, ".rds"))
          saveRDS(summaries[[summary_type]], filepath)
        }
        
        exported_files[[summary_type]] <- filepath
      }
    }
  }
  
  # Export candidates
  if ("candidates" %in% export_types && !is.null(project$results$candidates)) {
    
    filename <- "candidate_list"
    
    if (format == "csv") {
      filepath <- file.path(output_folder, paste0(filename, ".csv"))
      write.csv(project$results$candidates, filepath, row.names = FALSE)
      
    } else if (format == "excel") {
      suppressMessages(suppressWarnings(require(openxlsx)))
      filepath <- file.path(output_folder, paste0(filename, ".xlsx"))
      write.xlsx(project$results$candidates, filepath)
      
    } else if (format == "rds") {
      filepath <- file.path(output_folder, paste0(filename, ".rds"))
      saveRDS(project$results$candidates, filepath)
    }
    
    exported_files[["candidates"]] <- filepath
  }
  
  # Export raw results
  if ("raw_results" %in% export_types) {
    
    # Export entire project object
    filepath <- file.path(output_folder, "hitmap_project.rds")
    saveRDS(project, filepath)
    exported_files[["project_object"]] <- filepath
    
    # Export processing log
    log_filepath <- file.path(output_folder, "processing_log.txt")
    writeLines(project$metadata$processing_log, log_filepath)
    exported_files[["processing_log"]] <- log_filepath
    
    # Export project status
    status <- hitmap_project_status(project)
    status_filepath <- file.path(output_folder, "project_status.rds")
    saveRDS(status, status_filepath)
    exported_files[["project_status"]] <- status_filepath
  }
  
  message(paste("Exported", length(exported_files), "files to", output_folder))
  
  return(exported_files)
}

#' Load and Resume Project
#' 
#' @param project_file Path to saved hitmap_project RDS file
#' @param additional_steps Additional steps to perform
#' @param config Updated configuration
#' @return Loaded and potentially updated hitmap_project object
#' 
#' @export
hitmap_load_project <- function(project_file, 
                              additional_steps = NULL,
                              config = list()) {
  
  if (!file.exists(project_file)) {
    stop("Project file not found:", project_file)
  }
  
  message("=== Loading HiTMaP Project ===")
  
  # Load project
  project <- readRDS(project_file)
  
  if (!inherits(project, "hitmap_project")) {
    stop("File does not contain a valid hitmap_project object")
  }
  
  message(paste("Loaded project created:", project$metadata$created))
  message(paste("Project folder:", project$metadata$project_folder))
  
  # Print current status
  print(project)
  
  # Resume with additional steps if requested
  if (!is.null(additional_steps)) {
    project <- hitmap_resume_workflow(project, additional_steps, config)
  }
  
  return(project)
}