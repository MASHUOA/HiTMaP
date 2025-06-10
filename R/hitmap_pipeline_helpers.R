#' HiTMaP Pipeline Helper Functions
#' 
#' Additional helper functions to enhance the pipeline experience
#' and provide specialized data manipulation capabilities

#' Pipeline-aware configuration builder
#' 
#' @param base_config Base configuration to start with
#' @param ... Named parameters to add or override
#' @return Configuration list suitable for pipeline use
#' 
#' @export
build_config <- function(base_config = list(), ...) {
  
  user_params <- list(...)
  
  # Merge configurations
  final_config <- modifyList(base_config, user_params)
  
  # Add class for potential future methods
  class(final_config) <- c("hitmap_config", "list")
  
  return(final_config)
}

#' Validate and complete preprocessing parameters
#' 
#' Ensures that all required preprocessing parameters are present with proper defaults
#' to prevent errors during imaging_identification calls. Based on Cardinal 3 preprocessing methods.
#' 
#' @param preprocess User-provided preprocessing list (can be partial)
#' @param ppm ppm tolerance for peakAlign (default 5)
#' @return Complete preprocessing list with all required parameters
#' 
#' @details
#' Cardinal 3 supports these preprocessing methods:
#' \itemize{
#'   \item smoothSignal: "disable", "gaussian", "bilateral", "adaptive", "diff", "guide", "pag", "sgolay", "ma"
#'   \item reduceBaseline: "locmin", "hull", "snip", "median"
#'   \item peakPick: "diff", "sd", "mad", "quantile", "filter", "cwt", "adaptive"
#'   \item normalize: "tic", "rms", "reference"
#' }
#' 
#' @export
validate_preprocess_params <- function(preprocess = list(), ppm = 5) {
  
  # Define complete default preprocessing structure based on Cardinal 3
  default_preprocess <- list(
    force_preprocess = FALSE,
    use_preprocessRDS = TRUE,
    smoothSignal = list(method = "disable"),
    reduceBaseline = list(method = "locmin"),
    peakPick = list(method = "adaptive"),
    peakAlign = list(tolerance = ppm/2, units = "ppm"),
    peakFilter = list(freq.min = 0.05),
    normalize = list(method = "rms", mz = 1)
  )
  
  # Valid method options from Cardinal 3
  valid_methods <- list(
    smoothSignal = c("disable", "gaussian", "bilateral", "adaptive", "diff", "guide", "pag", "sgolay", "ma"),
    reduceBaseline = c("locmin", "hull", "snip", "median"),
    peakPick = c("diff", "sd", "mad", "quantile", "filter", "cwt", "adaptive"),
    normalize = c("tic", "rms", "reference")
  )
  
  # Merge user parameters with defaults
  # User parameters override defaults
  complete_preprocess <- modifyList(default_preprocess, preprocess)
  
  # Validate critical parameters and methods
  if (is.null(complete_preprocess$smoothSignal$method)) {
    complete_preprocess$smoothSignal$method <- "disable"
  } else if (!complete_preprocess$smoothSignal$method %in% valid_methods$smoothSignal) {
    warning("Invalid smoothSignal method. Using 'disable'. Valid options: ", 
            paste(valid_methods$smoothSignal, collapse = ", "))
    complete_preprocess$smoothSignal$method <- "disable"
  }
  
  if (is.null(complete_preprocess$reduceBaseline$method)) {
    complete_preprocess$reduceBaseline$method <- "locmin"
  } else if (!complete_preprocess$reduceBaseline$method %in% valid_methods$reduceBaseline) {
    warning("Invalid reduceBaseline method. Using 'locmin'. Valid options: ", 
            paste(valid_methods$reduceBaseline, collapse = ", "))
    complete_preprocess$reduceBaseline$method <- "locmin"
  }
  
  if (is.null(complete_preprocess$peakPick$method)) {
    complete_preprocess$peakPick$method <- "adaptive"
  } else if (!complete_preprocess$peakPick$method %in% valid_methods$peakPick) {
    warning("Invalid peakPick method. Using 'adaptive'. Valid options: ", 
            paste(valid_methods$peakPick, collapse = ", "))
    complete_preprocess$peakPick$method <- "adaptive"
  }
  
  if (is.null(complete_preprocess$normalize$method)) {
    complete_preprocess$normalize$method <- "rms"
  } else if (!complete_preprocess$normalize$method %in% valid_methods$normalize) {
    warning("Invalid normalize method. Using 'rms'. Valid options: ", 
            paste(valid_methods$normalize, collapse = ", "))
    complete_preprocess$normalize$method <- "rms"
  }
  
  # Update peakAlign tolerance if ppm was changed
  if (!is.null(complete_preprocess$peakAlign)) {
    complete_preprocess$peakAlign$tolerance <- ppm/2
  }
  
  return(complete_preprocess)
}

#' Apply configuration to project (Pipeline Method)
#' 
#' @param project hitmap_project object
#' @param config Configuration list or function that returns configuration
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
apply_config <- function(project, config, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # If config is a function, call it
  if (is.function(config)) {
    config <- config()
  }
  
  # Apply configuration
  project$metadata$config <- modifyList(project$metadata$config, config)
  
  # Log the change
  log_entry <- paste(Sys.time(), "- Applied configuration with", length(config), "parameters")
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  return(project)
}

#' Conditional pipeline execution
#' 
#' @param project hitmap_project object
#' @param condition Logical condition or function that returns logical
#' @param true_action Function to execute if condition is TRUE
#' @param false_action Function to execute if condition is FALSE (optional)
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
if_then <- function(project, condition, true_action, false_action = NULL, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Evaluate condition
  if (is.function(condition)) {
    condition_result <- condition(project)
  } else {
    condition_result <- condition
  }
  
  # Execute appropriate action
  if (condition_result) {
    if (is.function(true_action)) {
      project <- true_action(project)
    }
  } else {
    if (!is.null(false_action) && is.function(false_action)) {
      project <- false_action(project)
    }
  }
  
  return(project)
}

#' Repeat pipeline step until condition is met
#' 
#' @param project hitmap_project object
#' @param action Function to execute repeatedly
#' @param condition Function that returns TRUE when to stop
#' @param max_iterations Maximum number of iterations (safety)
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
repeat_until <- function(project, action, condition, max_iterations = 10, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  iteration <- 0
  
  while (!condition(project) && iteration < max_iterations) {
    project <- action(project)
    iteration <- iteration + 1
  }
  
  if (iteration >= max_iterations) {
    warning("Maximum iterations reached in repeat_until")
  }
  
  return(project)
}

#' Try pipeline step with error handling
#' 
#' @param project hitmap_project object
#' @param action Function to execute
#' @param on_error Function to execute on error (optional)
#' @param silent Suppress error messages (default FALSE)
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
try_step <- function(project, action, on_error = NULL, silent = FALSE, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  tryCatch({
    project <- action(project)
    return(project)
  }, error = function(e) {
    
    if (!silent) {
      warning(paste("Pipeline step failed:", e$message))
    }
    
    # Log the error
    log_entry <- paste(Sys.time(), "- Pipeline step failed:", e$message)
    project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
    
    # Execute error handler if provided
    if (!is.null(on_error) && is.function(on_error)) {
      project <- on_error(project, e)
    }
    
    return(project)
  })
}

#' Validate pipeline state
#' 
#' @param project hitmap_project object
#' @param requirements List of requirements to check
#' @param action What to do if validation fails: "stop", "warn", "skip"
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
validate_state <- function(project, requirements, action = "stop", .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  validation_errors <- character(0)
  
  # Check each requirement
  for (req in requirements) {
    if (req == "initialized" && !project$status$initialized) {
      validation_errors <- c(validation_errors, "Project must be initialized")
    } else if (req == "candidates" && !project$status$candidates_generated) {
      validation_errors <- c(validation_errors, "Candidates must be generated")
    } else if (req == "preprocessed" && length(project$status$files_preprocessed) == 0) {
      validation_errors <- c(validation_errors, "Files must be preprocessed")
    } else if (req == "segmented" && length(project$status$files_segmented) == 0) {
      validation_errors <- c(validation_errors, "Files must be segmented")
    } else if (req == "searched" && length(project$status$files_pmf_searched) == 0) {
      validation_errors <- c(validation_errors, "PMF search must be completed")
    } else if (req == "summarized" && length(project$status$summaries_generated) == 0) {
      validation_errors <- c(validation_errors, "Summaries must be generated")
    }
  }
  
  # Handle validation errors
  if (length(validation_errors) > 0) {
    error_msg <- paste("Validation failed:", paste(validation_errors, collapse = "; "))
    
    if (action == "stop") {
      stop(error_msg)
    } else if (action == "warn") {
      warning(error_msg)
    }
    
    # Log validation failure
    log_entry <- paste(Sys.time(), "- Validation failed:", error_msg)
    project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  }
  
  return(project)
}

#' Debug pipeline step
#' 
#' @param project hitmap_project object
#' @param message Debug message to display
#' @param level Debug level: "info", "debug", "trace"
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object (unchanged)
#' 
#' @export
debug_step <- function(project, message, level = "debug", .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Check if debug mode is enabled (could be in config)
  debug_enabled <- project$metadata$config$debug_mode %||% FALSE
  
  if (debug_enabled || level == "info") {
    formatted_message <- paste("[", toupper(level), "]", message)
    message(formatted_message)
    
    # Log debug message
    log_entry <- paste(Sys.time(), "- DEBUG:", message)
    project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  }
  
  return(project)
}

#' Benchmark pipeline step
#' 
#' @param project hitmap_project object
#' @param step_name Name of the step being benchmarked
#' @param action Function to execute and benchmark
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
benchmark_step <- function(project, step_name, action, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Record start time
  start_time <- Sys.time()
  
  # Execute action
  project <- action(project)
  
  # Record end time and calculate duration
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Store benchmark information
  if (is.null(project$metadata$benchmarks)) {
    project$metadata$benchmarks <- list()
  }
  
  project$metadata$benchmarks[[step_name]] <- list(
    start_time = start_time,
    end_time = end_time,
    duration_seconds = duration
  )
  
  # Log benchmark
  log_entry <- paste(Sys.time(), "- Benchmark:", step_name, "completed in", 
                    round(duration, 2), "seconds")
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  message(paste("Step", step_name, "completed in", round(duration, 2), "seconds"))
  
  return(project)
}

#' Set pipeline options
#' 
#' @param project hitmap_project object
#' @param ... Named options to set
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
set_options <- function(project, ..., .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  options <- list(...)
  
  # Initialize options if not present
  if (is.null(project$metadata$pipeline_options)) {
    project$metadata$pipeline_options <- list()
  }
  
  # Update options
  project$metadata$pipeline_options <- modifyList(project$metadata$pipeline_options, options)
  
  # Log option changes
  log_entry <- paste(Sys.time(), "- Set pipeline options:", 
                    paste(names(options), "=", options, collapse = ", "))
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  return(project)
}

#' Get pipeline option value
#' 
#' @param project hitmap_project object
#' @param option_name Name of option to retrieve
#' @param default_value Default value if option not found
#' @return Option value
#' 
#' @export
get_option <- function(project, option_name, default_value = NULL) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (is.null(project$metadata$pipeline_options)) {
    return(default_value)
  }
  
  return(project$metadata$pipeline_options[[option_name]] %||% default_value)
}

#' Create pipeline branch
#' 
#' @param project hitmap_project object
#' @param branch_name Name for the branch
#' @param ... Pipeline steps to execute in the branch
#' @return List containing original project and branched project
#' 
#' @export
branch_pipeline <- function(project, branch_name, ...) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Create a copy for the branch
  branched_project <- project
  
  # Add branch information
  branched_project$metadata$branch_name <- branch_name
  branched_project$metadata$branched_from <- project$metadata$created
  branched_project$metadata$created <- Sys.time()
  
  # Execute pipeline steps on the branch
  pipeline_steps <- list(...)
  
  for (step in pipeline_steps) {
    if (is.function(step)) {
      branched_project <- step(branched_project)
    }
  }
  
  # Log branching
  log_entry <- paste(Sys.time(), "- Created branch:", branch_name)
  project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
  
  return(list(
    original = project,
    branch = branched_project
  ))
}

#' Merge pipeline branches
#' 
#' @param main_project Main hitmap_project object
#' @param branch_project Branched hitmap_project object
#' @param merge_strategy Strategy for merging: "replace", "combine"
#' @return Merged hitmap_project object
#' 
#' @export
merge_branches <- function(main_project, branch_project, merge_strategy = "replace") {
  
  if (!inherits(main_project, "hitmap_project") || !inherits(branch_project, "hitmap_project")) {
    stop("Both inputs must be hitmap_project objects")
  }
  
  if (merge_strategy == "replace") {
    # Replace main with branch results
    merged_project <- branch_project
    merged_project$metadata$branch_name <- NULL
    merged_project$metadata$branched_from <- NULL
    
  } else if (merge_strategy == "combine") {
    # Combine results from both branches
    merged_project <- c(main_project, branch_project)
    
  } else {
    stop("Invalid merge strategy. Use 'replace' or 'combine'")
  }
  
  # Log merge
  log_entry <- paste(Sys.time(), "- Merged branch using strategy:", merge_strategy)
  merged_project$metadata$processing_log <- c(merged_project$metadata$processing_log, log_entry)
  
  return(merged_project)
}

#' Create pipeline template
#' 
#' @param name Template name
#' @param ... Pipeline functions to include in template
#' @return Function that can be applied to a hitmap_project
#' 
#' @export
create_template <- function(name, ...) {
  
  pipeline_steps <- list(...)
  
  # Return a function that applies the template
  template_function <- function(project) {
    
    # Add template information to project
    project$metadata$template_applied <- name
    project$metadata$template_applied_time <- Sys.time()
    
    # Execute all pipeline steps
    for (step in pipeline_steps) {
      if (is.function(step)) {
        project <- step(project)
      }
    }
    
    # Log template application
    log_entry <- paste(Sys.time(), "- Applied template:", name)
    project$metadata$processing_log <- c(project$metadata$processing_log, log_entry)
    
    return(project)
  }
  
  # Add template metadata
  attr(template_function, "template_name") <- name
  attr(template_function, "template_steps") <- length(pipeline_steps)
  class(template_function) <- c("hitmap_template", "function")
  
  return(template_function)
}

#' Apply template to project
#' 
#' @param project hitmap_project object
#' @param template Template function created by create_template
#' @param .validate Validate inputs (default TRUE)
#' @return Updated hitmap_project object
#' 
#' @export
apply_template <- function(project, template, .validate = TRUE) {
  
  if (.validate && !inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  if (!inherits(template, "hitmap_template")) {
    stop("Template must be created with create_template()")
  }
  
  # Apply the template
  project <- template(project)
  
  return(project)
}

#' Null coalescing operator helper
#' 
#' @param a First value
#' @param b Second value (used if a is NULL)
#' @return a if not NULL, otherwise b
#' 
#' @export
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}

#' Check if project meets condition
#' 
#' @param project hitmap_project object
#' @param condition Condition to check
#' @return Logical value
#' 
#' @export
meets_condition <- function(project, condition) {
  
  if (!inherits(project, "hitmap_project")) {
    return(FALSE)
  }
  
  if (condition == "initialized") {
    return(project$status$initialized)
  } else if (condition == "has_candidates") {
    return(project$status$candidates_generated)
  } else if (condition == "has_preprocessing") {
    return(length(project$status$files_preprocessed) > 0)
  } else if (condition == "has_segmentation") {
    return(length(project$status$files_segmented) > 0)
  } else if (condition == "has_pmf_results") {
    return(length(project$status$files_pmf_searched) > 0)
  } else if (condition == "has_summaries") {
    return(length(project$status$summaries_generated) > 0)
  } else if (condition == "complete") {
    return(length(project$status$files_pmf_searched) == length(project$metadata$data_files) &&
           length(project$status$summaries_generated) > 0)
  } else {
    return(FALSE)
  }
}

#' Pipeline progress bar
#' 
#' @param project hitmap_project object
#' @param width Width of progress bar (characters)
#' @return Character string representing progress
#' 
#' @export
progress_bar <- function(project, width = 50) {
  
  if (!inherits(project, "hitmap_project")) {
    stop("Input must be a hitmap_project object")
  }
  
  # Calculate progress based on completed steps
  total_steps <- 5  # init, candidates, preprocess, segment, pmf, summarize
  completed_steps <- 0
  
  if (project$status$initialized) completed_steps <- completed_steps + 1
  if (project$status$candidates_generated) completed_steps <- completed_steps + 1
  if (length(project$status$files_preprocessed) == length(project$metadata$data_files)) completed_steps <- completed_steps + 1
  if (length(project$status$files_segmented) == length(project$metadata$data_files)) completed_steps <- completed_steps + 1
  if (length(project$status$files_pmf_searched) == length(project$metadata$data_files)) completed_steps <- completed_steps + 1
  
  # Calculate progress percentage
  progress_pct <- completed_steps / total_steps
  filled_width <- round(width * progress_pct)
  empty_width <- width - filled_width
  
  # Create progress bar
  filled_part <- paste(rep("█", filled_width), collapse = "")
  empty_part <- paste(rep("░", empty_width), collapse = "")
  
  progress_str <- paste0("[", filled_part, empty_part, "] ", 
                        round(progress_pct * 100), "% (", 
                        completed_steps, "/", total_steps, " steps)")
  
  return(progress_str)
}