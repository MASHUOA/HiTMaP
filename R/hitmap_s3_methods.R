#' S3 Methods for HiTMaP Project Object
#' 
#' S3 methods to enhance the pipeline experience and provide 
#' seamless integration with R's generic functions

#' Print method for hitmap_project (Enhanced)
#' 
#' @param x hitmap_project object
#' @param ... Additional arguments (ignored)
#' @method print hitmap_project
#' @export
print.hitmap_project <- function(x, ...) {
  cat("HiTMaP Project Pipeline\n")
  cat("=======================\n")
  
  # Project info
  cat("Created:", format(x$metadata$created, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Files:", length(x$metadata$data_files), 
      paste0("(", paste(x$metadata$data_files, collapse = ", "), ")"), "\n")
  cat("Project folder:", x$metadata$project_folder, "\n")
  
  # Pipeline progress
  cat("\nPipeline Progress:\n")
  
  # Progress indicators
  init_status <- if (x$status$initialized) "✓" else "○"
  candidates_status <- if (x$status$candidates_generated) "✓" else "○"
  preprocess_status <- if (length(x$status$files_preprocessed) == length(x$metadata$data_files)) "✓" 
                      else if (length(x$status$files_preprocessed) > 0) "◐" else "○"
  segment_status <- if (length(x$status$files_segmented) == length(x$metadata$data_files)) "✓" 
                   else if (length(x$status$files_segmented) > 0) "◐" else "○"
  pmf_status <- if (length(x$status$files_pmf_searched) == length(x$metadata$data_files)) "✓" 
               else if (length(x$status$files_pmf_searched) > 0) "◐" else "○"
  summary_status <- if (length(x$status$summaries_generated) > 0) "✓" else "○"
  
  cat("  ", init_status, " ims_init()           - Initialize project\n")
  cat("  ", candidates_status, " ims_generate_candidates() - Generate candidate list\n")
  cat("  ", preprocess_status, " ims_preprocess()     - Preprocess data files\n")
  cat("  ", segment_status, " ims_segment()        - Segment into regions\n")
  cat("  ", pmf_status, " ims_search_pmf()     - Perform PMF search\n")
  cat("  ", summary_status, " ims_summarize()      - Generate summaries\n")
  
  # Results summary
  if (x$status$candidates_generated) {
    cat("\nResults Summary:\n")
    cat("  Candidates:", nrow(x$results$candidates), "\n")
  }
  
  if (length(x$status$files_preprocessed) > 0) {
    cat("  Preprocessed files:", length(x$status$files_preprocessed), "/", length(x$metadata$data_files), "\n")
  }
  
  if (length(x$status$files_segmented) > 0) {
    total_regions <- sum(sapply(x$results$segmentation, function(seg) length(seg$segmentation_labels)))
    cat("  Segmented files:", length(x$status$files_segmented), "/", length(x$metadata$data_files), 
        "(", total_regions, "regions )\n")
  }
  
  if (length(x$status$files_pmf_searched) > 0) {
    cat("  PMF searched files:", length(x$status$files_pmf_searched), "/", length(x$metadata$data_files), "\n")
  }
  
  if (!is.null(x$results$summaries$peptide_summary)) {
    cat("  Total peptides:", nrow(x$results$summaries$peptide_summary), "\n")
  }
  
  if (!is.null(x$results$summaries$protein_summary)) {
    cat("  Total proteins:", nrow(x$results$summaries$protein_summary), "\n")
  }
  
  # Next step suggestion
  cat("\nNext Pipeline Step:\n")
  if (!x$status$initialized) {
    cat("  project %>% init(thread_count = 4)\n")
  } else if (!x$status$candidates_generated) {
    cat("  project %>% ims_generate_candidates(database = \"your_database.fasta\")\n")
  } else if (length(x$status$files_preprocessed) < length(x$metadata$data_files)) {
    cat("  project %>% ims_preprocess()\n")
  } else if (length(x$status$files_segmented) < length(x$metadata$data_files)) {
    cat("  project %>% ims_segment(segment_count = 4)\n")
  } else if (length(x$status$files_pmf_searched) < length(x$metadata$data_files)) {
    cat("  project %>% ims_search_pmf()\n")
  } else if (length(x$status$summaries_generated) == 0) {
    cat("  project %>% ims_summarize()\n")
  } else {
    cat("  Pipeline complete! Use export_results() or checkpoint()\n")
  }
  
  invisible(x)
}

#' Summary method for hitmap_project
#' 
#' @param object hitmap_project object
#' @param ... Additional arguments (ignored)
#' @method summary hitmap_project
#' @export
summary.hitmap_project <- function(object, ...) {
  cat("HiTMaP Project Summary\n")
  cat("======================\n")
  
  # Basic information
  cat("Project created:", format(object$metadata$created, "%Y-%m-%d %H:%M:%S"), "\n")
  cat("Data files:", length(object$metadata$data_files), "\n")
  cat("Working directory:", object$metadata$project_folder, "\n")
  
  if (object$status$initialized) {
    cat("Threads configured:", object$setup$thread_count, "\n")
  }
  
  # Processing statistics
  cat("\nProcessing Statistics:\n")
  
  if (object$status$candidates_generated) {
    candidates <- object$results$candidates
    cat("  Candidates generated:", nrow(candidates), "\n")
    cat("    Mass range:", round(min(candidates$mz), 2), "-", round(max(candidates$mz), 2), "m/z\n")
    if ("Protein" %in% colnames(candidates)) {
      cat("    Unique proteins:", length(unique(candidates$Protein)), "\n")
    }
  }
  
  # File processing details
  for (file_name in object$metadata$data_files) {
    cat("\nFile:", file_name, "\n")
    
    if (file_name %in% object$status$files_preprocessed) {
      preprocess_result <- object$results$preprocessing[[file_name]]
      cat("  ✓ Preprocessed (", format(preprocess_result$processing_time, "%H:%M:%S"), ")\n")
    } else {
      cat("  ○ Not preprocessed\n")
    }
    
    if (file_name %in% object$status$files_segmented) {
      segment_result <- object$results$segmentation[[file_name]]
      region_count <- length(segment_result$segmentation_labels)
      cat("  ✓ Segmented into", region_count, "regions\n")
      
      # Show region sizes
      region_sizes <- sapply(segment_result$segmentation_labels, length)
      cat("    Region sizes:", paste(region_sizes, collapse = ", "), "pixels\n")
    } else {
      cat("  ○ Not segmented\n")
    }
    
    if (file_name %in% object$status$files_pmf_searched) {
      pmf_results <- object$results$pmf_search[[file_name]]
      total_peptides <- sum(sapply(pmf_results, function(x) nrow(x$peptides)))
      total_proteins <- sum(sapply(pmf_results, function(x) nrow(x$proteins)))
      cat("  ✓ PMF searched -", total_peptides, "peptides,", total_proteins, "proteins\n")
      
      # Show results by region
      for (region_name in names(pmf_results)) {
        region_result <- pmf_results[[region_name]]
        cat("    ", region_name, ":", nrow(region_result$peptides), "peptides,", 
            nrow(region_result$proteins), "proteins\n")
      }
    } else {
      cat("  ○ Not PMF searched\n")
    }
  }
  
  # Summary results
  if (length(object$status$summaries_generated) > 0) {
    cat("\nSummary Results:\n")
    
    for (summary_type in object$status$summaries_generated) {
      if (summary_type == "protein" && !is.null(object$results$summaries$protein_summary)) {
        protein_data <- object$results$summaries$protein_summary
        cat("  Protein summary:", nrow(protein_data), "entries\n")
        
        if ("source" %in% colnames(protein_data)) {
          cat("    By source:", paste(table(protein_data$source), collapse = ", "), "\n")
        }
      }
      
      if (summary_type == "peptide" && !is.null(object$results$summaries$peptide_summary)) {
        peptide_data <- object$results$summaries$peptide_summary
        cat("  Peptide summary:", nrow(peptide_data), "entries\n")
        
        if ("source" %in% colnames(peptide_data)) {
          cat("    By source:", paste(table(peptide_data$source), collapse = ", "), "\n")
        }
      }
      
      if (summary_type == "region" && !is.null(object$results$summaries$region_summary)) {
        region_data <- object$results$summaries$region_summary
        cat("  Region summary:", nrow(region_data), "entries\n")
      }
    }
  }
  
  # Processing log (last few entries)
  if (length(object$metadata$processing_log) > 0) {
    cat("\nRecent Processing Log:\n")
    recent_log <- tail(object$metadata$processing_log, 3)
    for (entry in recent_log) {
      cat("  ", entry, "\n")
    }
  }
  
  invisible(object)
}

#' Length method for hitmap_project (returns number of files)
#' 
#' @param x hitmap_project object
#' @method length hitmap_project
#' @export
length.hitmap_project <- function(x) {
  return(length(x$metadata$data_files))
}

#' Names method for hitmap_project (returns file names)
#' 
#' @param x hitmap_project object
#' @method names hitmap_project
#' @export
names.hitmap_project <- function(x) {
  return(x$metadata$data_files)
}

#' Subset method for hitmap_project (subset by file names)
#' 
#' @param x hitmap_project object
#' @param i Indices or names of files to subset
#' @param ... Additional arguments (ignored)
#' @method `[` hitmap_project
#' @export
`[.hitmap_project` <- function(x, i, ...) {
  
  # Get file names to keep
  if (is.numeric(i)) {
    files_to_keep <- x$metadata$data_files[i]
  } else if (is.character(i)) {
    files_to_keep <- intersect(i, x$metadata$data_files)
  } else if (is.logical(i)) {
    files_to_keep <- x$metadata$data_files[i]
  } else {
    stop("Invalid subset index")
  }
  
  if (length(files_to_keep) == 0) {
    stop("No valid files found in subset")
  }
  
  # Update metadata
  x$metadata$data_files <- files_to_keep
  x$metadata$data_files_imzml <- paste0(files_to_keep, ".imzML")
  
  # Update status
  x$status$files_preprocessed <- intersect(x$status$files_preprocessed, files_to_keep)
  x$status$files_segmented <- intersect(x$status$files_segmented, files_to_keep)
  x$status$files_pmf_searched <- intersect(x$status$files_pmf_searched, files_to_keep)
  
  # Update results
  x$results$preprocessing <- x$results$preprocessing[names(x$results$preprocessing) %in% files_to_keep]
  x$results$segmentation <- x$results$segmentation[names(x$results$segmentation) %in% files_to_keep]
  x$results$pmf_search <- x$results$pmf_search[names(x$results$pmf_search) %in% files_to_keep]
  
  # Clear summaries as they may no longer be valid
  x$results$summaries <- list()
  x$status$summaries_generated <- character(0)
  
  # Add to log
  log_entry <- paste(Sys.time(), "- Subsetted project to files:", paste(files_to_keep, collapse = ", "))
  x$metadata$processing_log <- c(x$metadata$processing_log, log_entry)
  
  return(x)
}

#' Extract results from hitmap_project
#' 
#' @param x hitmap_project object
#' @param name Name of result to extract
#' @method `$` hitmap_project
#' @export
`$.hitmap_project` <- function(x, name) {
  
  # Allow easy access to common results
  if (name == "candidates") {
    return(x$results$candidates)
  } else if (name == "peptides") {
    return(x$results$summaries$peptide_summary)
  } else if (name == "proteins") {
    return(x$results$summaries$protein_summary)
  } else if (name == "regions") {
    return(x$results$summaries$region_summary)
  } else if (name == "status") {
    return(hitmap_project_status(x))
  } else if (name == "log") {
    return(x$metadata$processing_log)
  } else if (name == "config") {
    return(x$metadata$config)
  } else {
    # Default to normal list extraction
    return(NextMethod())
  }
}

#' Combine hitmap_project objects
#' 
#' @param ... hitmap_project objects to combine
#' @method c hitmap_project
#' @export
c.hitmap_project <- function(...) {
  
  projects <- list(...)
  
  # Validate all inputs are hitmap_project objects
  if (!all(sapply(projects, function(x) inherits(x, "hitmap_project")))) {
    stop("All objects must be hitmap_project objects")
  }
  
  if (length(projects) == 1) {
    return(projects[[1]])
  }
  
  # Use first project as base
  combined_project <- projects[[1]]
  
  # Combine data files
  all_files <- unique(unlist(lapply(projects, function(x) x$metadata$data_files)))
  
  # Check for file name conflicts
  file_counts <- table(unlist(lapply(projects, function(x) x$metadata$data_files)))
  if (any(file_counts > 1)) {
    duplicate_files <- names(file_counts)[file_counts > 1]
    warning(paste("Duplicate file names found:", paste(duplicate_files, collapse = ", ")))
  }
  
  # Update metadata
  combined_project$metadata$data_files <- all_files
  combined_project$metadata$data_files_imzml <- paste0(all_files, ".imzML")
  combined_project$metadata$created <- Sys.time()
  
  # Combine processing status
  combined_project$status$files_preprocessed <- unique(unlist(lapply(projects, function(x) x$status$files_preprocessed)))
  combined_project$status$files_segmented <- unique(unlist(lapply(projects, function(x) x$status$files_segmented)))
  combined_project$status$files_pmf_searched <- unique(unlist(lapply(projects, function(x) x$status$files_pmf_searched)))
  
  # Combine results
  all_preprocessing <- list()
  all_segmentation <- list()
  all_pmf_search <- list()
  
  for (project in projects) {
    all_preprocessing <- c(all_preprocessing, project$results$preprocessing)
    all_segmentation <- c(all_segmentation, project$results$segmentation)
    all_pmf_search <- c(all_pmf_search, project$results$pmf_search)
  }
  
  combined_project$results$preprocessing <- all_preprocessing
  combined_project$results$segmentation <- all_segmentation
  combined_project$results$pmf_search <- all_pmf_search
  
  # Clear summaries - they need to be regenerated
  combined_project$results$summaries <- list()
  combined_project$status$summaries_generated <- character(0)
  
  # Combine processing logs
  all_logs <- unlist(lapply(projects, function(x) x$metadata$processing_log))
  combined_project$metadata$processing_log <- c(all_logs, 
    paste(Sys.time(), "- Combined", length(projects), "projects"))
  
  return(combined_project)
}

#' Check if hitmap_project is valid
#' 
#' @param x hitmap_project object
#' @return Logical indicating if project is valid
#' @export
is_valid_project <- function(x) {
  
  if (!inherits(x, "hitmap_project")) {
    return(FALSE)
  }
  
  # Check required structure
  required_elements <- c("metadata", "setup", "results", "status")
  if (!all(required_elements %in% names(x))) {
    return(FALSE)
  }
  
  # Check metadata structure
  required_metadata <- c("data_files", "project_folder", "created")
  if (!all(required_metadata %in% names(x$metadata))) {
    return(FALSE)
  }
  
  # Check status structure
  required_status <- c("initialized", "candidates_generated")
  if (!all(required_status %in% names(x$status))) {
    return(FALSE)
  }
  
  return(TRUE)
}

#' Convert hitmap_project to data.frame
#' 
#' @param x hitmap_project object
#' @param row.names Ignored
#' @param optional Ignored
#' @param result_type Type of result to convert: "peptides", "proteins", "summary"
#' @param ... Additional arguments
#' @method as.data.frame hitmap_project
#' @export
as.data.frame.hitmap_project <- function(x, row.names = NULL, optional = FALSE, 
                                        result_type = "peptides", ...) {
  
  if (result_type == "peptides" && !is.null(x$results$summaries$peptide_summary)) {
    return(as.data.frame(x$results$summaries$peptide_summary))
  } else if (result_type == "proteins" && !is.null(x$results$summaries$protein_summary)) {
    return(as.data.frame(x$results$summaries$protein_summary))
  } else if (result_type == "candidates" && !is.null(x$results$candidates)) {
    return(as.data.frame(x$results$candidates))
  } else if (result_type == "summary") {
    # Create a summary data frame
    status <- hitmap_project_status(x)
    summary_df <- data.frame(
      metric = c("total_files", "preprocessed_files", "segmented_files", 
                "pmf_searched_files", "total_regions", "total_peptides", "total_proteins"),
      value = c(status$project_info$total_files,
               status$processing_progress$files_preprocessed,
               status$processing_progress$files_segmented,
               status$processing_progress$files_pmf_searched,
               status$results_summary$total_regions,
               status$results_summary$total_peptides,
               status$results_summary$total_proteins),
      stringsAsFactors = FALSE
    )
    return(summary_df)
  } else {
    stop(paste("No data available for result_type:", result_type))
  }
}