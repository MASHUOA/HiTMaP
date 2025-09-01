#' UniMod Modification Helper Functions
#'
#' Enhanced functions for handling protein modifications using the UniMod database
#' with improved ambiguity resolution and user-friendly input parsing.

#' Retrieve and Format UniMod Modifications
#'
#' This function provides an enhanced interface to the UniMod database for
#' organizing modification arguments in the correct format for HiTMaP workflows.
#' It handles ambiguous modification names and provides multiple matching options.
#'
#' @param modifications Character vector or list of modification names/IDs
#' @param positions Character vector of positions ("any", "N-term", "C-term", etc.)
#' @param mod_type Character, either "fixed" or "variable" (default: "variable")
#' @param update_unimod Logical, whether to update the UniMod database (default: FALSE)
#' @param interactive Logical, whether to prompt user for ambiguous matches (default: TRUE)
#' @param show_matches Logical, whether to display found matches (default: TRUE)
#' 
#' @return List with properly formatted modification arguments for HiTMaP workflows:
#'   - modifications: Named list with fixed/variable modifications and positions
#'   - summary: Summary of matched modifications
#'   - ambiguous: List of ambiguous matches requiring user attention
#'   - not_found: Vector of modifications that couldn't be matched
#'
#' @examples
#' # Simple usage
#' mod_result <- format_unimod_modifications(
#'   modifications = c("Oxidation", "Acetyl"),
#'   mod_type = "variable"
#' )
#' 
#' # Use the result in imaging_identification
#' imaging_identification(
#'   datafile = "data.imzML",
#'   Modifications = mod_result$modifications,
#'   # ... other parameters
#' )
#' 
#' # Multiple modifications with positions
#' mod_result <- format_unimod_modifications(
#'   modifications = c("Carbamidomethyl", "Oxidation", "Phospho"),
#'   positions = c("C", "M", "STY"),
#'   mod_type = "variable"
#' )
#'
#' @export
format_unimod_modifications <- function(modifications = NULL,
                                       positions = NULL,
                                       mod_type = "variable",
                                       update_unimod = FALSE,
                                       interactive = TRUE,
                                       show_matches = TRUE) {
  
  # Input validation
  if (is.null(modifications)) {
    return(list(
      modifications = list(fixed = NULL, fixmod_position = NULL, fixmod_one_letter = NULL,
                          variable = NULL, varmod_position = NULL, varmod_one_letter = NULL),
      summary = "No modifications specified",
      ambiguous = list(),
      not_found = character(0)
    ))
  }
  
  # Ensure required packages are loaded
  suppressMessages(suppressWarnings(require(protViz)))
  suppressMessages(suppressWarnings(require(XML)))
  suppressMessages(suppressWarnings(require(stringr)))
  
  # Load or update UniMod database
  unimod_data <- load_unimod_database(update_unimod)
  
  # Default positions if not specified
  if (is.null(positions)) {
    positions <- rep("any", length(modifications))
  } else if (length(positions) == 1 && length(modifications) > 1) {
    positions <- rep(positions, length(modifications))
  } else if (length(positions) != length(modifications)) {
    stop("Length of positions must match length of modifications or be 1")
  }
  
  # Direct matching with immediate position filtering
  matched_mods <- data.frame()
  failed_mods <- character()
  
  for (i in seq_along(modifications)) {
    mod_name <- modifications[i]
    mod_pos <- positions[i]
    
    if (show_matches) {
      message("Searching for: ", mod_name, " at position: ", mod_pos)
    }
    
    # Use the rewritten search function
    match_result <- find_modification_matches(mod_name, mod_pos, unimod_data)
    
    if (nrow(match_result) == 1 && match_result$code_name[1] == "FAILED_MATCH") {
      failed_mods <- c(failed_mods, mod_name)
    } else {
      matched_mods <- rbind(matched_mods, match_result[1, ])
      if (show_matches) {
        message("Matched: ", mod_name, " -> ", match_result$full_name[1])
      }
    }
  }
  
  # Format results
  if (nrow(matched_mods) == 0) {
    return(list(
      modifications = list(fixed = NULL, fixmod_position = NULL, fixmod_one_letter = NULL,
                          variable = NULL, varmod_position = NULL, varmod_one_letter = NULL),
      summary = "No modifications matched",
      not_found = failed_mods
    ))
  }
  
  # Create modification list
  mod_names <- matched_mods$code_name
  mod_positions <- matched_mods$one_letter
  mod_positions_type <- matched_mods$position_key
  
  if (mod_type == "fixed") {
    modifications_list <- list(
      fixed = mod_names,
      fixmod_position = mod_positions_type,
      fixmod_one_letter = mod_positions,
      variable = NULL,
      varmod_position = NULL,
      varmod_one_letter = NULL
    )
  } else {
    modifications_list <- list(
      fixed = NULL,
      fixmod_position = NULL,
      fixmod_one_letter = NULL,
      variable = mod_names,
      varmod_position = mod_positions_type,
      varmod_one_letter = mod_positions
    )
  }
  
  summary_text <- paste("Matched", length(mod_names), mod_type, "modifications:", 
                       paste(mod_names, collapse = ", "))
  
  return(list(
    modifications = modifications_list,
    summary = summary_text,
    not_found = failed_mods
  ))
}

#' Load UniMod Database
#'
#' Internal function to load or update the UniMod database
#'
#' @param update_unimod Logical, whether to update from web
#' @return Processed UniMod data frame
load_unimod_database <- function(update_unimod = FALSE) {
  
  # Update database if requested
  if (update_unimod) {
    message("Updating UniMod database from web...")
    tryCatch({
      unimodurl <- url("http://www.unimod.org/xml/unimod_tables.xml")
      unimod.list <- XML::xmlToList(
        XML::xmlParse(
          scan(unimodurl, what = character(), quiet = TRUE)
        )
      )
      save(unimod.list, file = paste0(path.package("HiTMaP"), "/data/unimod.list.rda"))
      message("UniMod database updated successfully")
    }, error = function(e) {
      warning("Failed to update UniMod database from web. Using local version.")
      message("Error: ", e$message)
    })
  }
  
  # Load existing database
  if (!exists("unimod.df", envir = globalenv())) {
    data("unimod.list", package = "HiTMaP")
    
    unimod.df <- lapply(unimod.list, function(x) {
      x <- unname(x)
      x <- data.frame(do.call('rbind', lapply(x, function(y) {
        y[match(names(x[[1]]), names(y))]
      })), stringsAsFactors = FALSE)
      
      if (".attrs" %in% colnames(x)) {
        attrs_tb <- data.frame(do.call('rbind', lapply(x$.attrs, function(y) {
          y[match(names(x$.attrs[[1]]), names(y))]
        })), stringsAsFactors = FALSE)
        x <- cbind(x, attrs_tb)
      }
      return(x)
    })
    
    assign("unimod.df", unimod.df, envir = globalenv())
  } else {
    unimod.df <- get("unimod.df", envir = globalenv())
  }
  
  # Process and merge UniMod tables
  unimod.modification.df <- merge(
    unimod.df$modifications, 
    unimod.df$specificity,
    by.x = "record_id", 
    by.y = "mod_key", 
    all.x = TRUE
  )
  
  # Process positions
  unimod.df$positions_new <- unimod.df$positions
  unimod.df$positions_new$position_ext <- stringr::str_replace(
    unimod.df$positions_new$position,
    "Anywhere|Any N-term|Any C-term", 
    "Peptide"
  )
  unimod.df$positions_new$position_ext <- stringr::str_replace(
    unimod.df$positions_new$position_ext,
    "Protein N-term|Protein C-term", 
    "Protein"
  )
  
  unimod.modification.df <- merge(
    unimod.modification.df,
    unimod.df$positions_new,
    by.x = "position_key",
    by.y = "record_id"
  )
  
  return(unimod.modification.df)
}

#' Match Modifications with UniMod Database
#'
#' Internal function to perform fuzzy matching of modification names
#'
#' @param modifications Character vector of modification names
#' @param positions Character vector of positions
#' @param unimod_data UniMod data frame
#' @param interactive Logical, whether to prompt for ambiguous matches
#' @param show_matches Logical, whether to show matches
#' @return List with matched, ambiguous, and not found modifications
match_modifications_with_unimod <- function(modifications, positions, unimod_data, 
                                           interactive = TRUE, show_matches = TRUE) {
  
  results <- list(
    matched = data.frame(),
    ambiguous = list(),
    not_found = character()
  )
  
  for (i in seq_along(modifications)) {
    mod_name <- modifications[i]
    mod_pos <- positions[i]
    
    # Try different matching strategies
    matches <- find_modification_matches(mod_name, mod_pos, unimod_data)
    
    if (nrow(matches) == 0) {
      # No matches found
      results$not_found <- c(results$not_found, mod_name)
      if (show_matches) {
        message(paste("Warning: No matches found for:", mod_name))
      }
      
    } else if (nrow(matches) == 1) {
      # Single match found
      results$matched <- rbind(results$matched, matches[1, ])
      if (show_matches) {
        message(paste("Matched:", mod_name, "->", matches$full_name[1]))
      }
      
    } else {
      # Multiple matches found - handle ambiguity
      if (show_matches) {
        message(paste("Multiple matches found for:", mod_name))
        print(matches[, c("full_name", "code_name", "one_letter", "position")])
      }
      
      if (interactive) {
        # Interactive selection
        cat("\nSelect the correct modification:\n")
        for (j in 1:nrow(matches)) {
          cat(sprintf("%d: %s (%s) at %s\n", j, 
                     matches$full_name[j], 
                     matches$code_name[j], 
                     matches$position[j]))
        }
        cat("0: Skip this modification\n")
        
        choice <- as.numeric(readline(paste("Enter choice (0-", nrow(matches), "): ")))
        
        if (!is.na(choice) && choice > 0 && choice <= nrow(matches)) {
          results$matched <- rbind(results$matched, matches[choice, ])
          message(paste("Selected:", matches$full_name[choice]))
        } else {
          results$not_found <- c(results$not_found, mod_name)
          message(paste("Skipped:", mod_name))
        }
      } else {
        # Non-interactive - store ambiguous matches for later resolution
        results$ambiguous[[mod_name]] <- matches
        # Use first match as default
        results$matched <- rbind(results$matched, matches[1, ])
        if (show_matches) {
          message(paste("Using first match for:", mod_name, "->", matches$full_name[1]))
          message("(Multiple matches available - see ambiguous results)")
        }
      }
    }
  }
  
  return(results)
}

#' Find Modification Matches
#'
#' Internal function to find matches using various strategies
#'
#' @param mod_name Character, modification name to search
#' @param mod_pos Character, position specification
#' @param unimod_data UniMod data frame
#' @return Data frame with matches
find_modification_matches <- function(mod_name, mod_pos, unimod_data) {
  
  # Priority 1: ex_code_name + position match
  if ("ex_code_name" %in% colnames(unimod_data) && tolower(mod_pos) != "any") {
    ex_pos_matches <- unimod_data[
      tolower(stringr::str_replace_all(unimod_data$ex_code_name, "_", "-")) == 
      tolower(stringr::str_replace_all(mod_name, "_", "-")), 
    ]
    if (nrow(ex_pos_matches) > 0) {
      filtered <- filter_by_position(ex_pos_matches, mod_pos)
      if (nrow(filtered) > 0 && any(filtered$hidden == 0)) {
        best_match <- filtered[filtered$hidden == 0, ][1, , drop = FALSE]
        return(best_match)
      }
    }
  }
  
  # Priority 2: code_name + position match  
  if (tolower(mod_pos) != "any") {
    code_pos_matches <- unimod_data[
      tolower(stringr::str_replace_all(unimod_data$code_name, "_", "-")) == 
      tolower(stringr::str_replace_all(mod_name, "_", "-")), 
    ]
    if (nrow(code_pos_matches) > 0) {
      filtered <- filter_by_position(code_pos_matches, mod_pos)
      if (nrow(filtered) > 0 && any(filtered$hidden == 0)) {
        best_match <- filtered[filtered$hidden == 0, ][1, , drop = FALSE]
        return(best_match)
      }
    }
  }
  
  # Priority 3: full_name + position match
  if (tolower(mod_pos) != "any") {
    full_pos_matches <- unimod_data[
      tolower(unimod_data$full_name) == tolower(mod_name), 
    ]
    if (nrow(full_pos_matches) > 0) {
      filtered <- filter_by_position(full_pos_matches, mod_pos)
      if (nrow(filtered) > 0 && any(filtered$hidden == 0)) {
        best_match <- filtered[filtered$hidden == 0, ][1, , drop = FALSE]
        return(best_match)
      }
    }
  }
  
  # Priority 4: Suggest partial matches and return failure signal
  partial_matches <- unimod_data[
    grepl(tolower(mod_name), tolower(unimod_data$code_name)) |
    grepl(tolower(mod_name), tolower(unimod_data$full_name)) |
    (("ex_code_name" %in% colnames(unimod_data)) && 
     grepl(tolower(mod_name), tolower(unimod_data$ex_code_name))), 
  ]
  
  if (nrow(partial_matches) > 0) {
    message("No exact match for '", mod_name, "' at position '", mod_pos, "'")
    message("Suggested alternatives:")
    suggestions <- partial_matches[1:min(5, nrow(partial_matches)), 
                                  c("record_id", "code_name", "ex_code_name", "one_letter")]
    print(suggestions)
  }
  
  # Return failure signal
  return(data.frame(record_id = 1, code_name = "FAILED_MATCH", 
                    ex_code_name = mod_name, one_letter = mod_pos,
                    hidden = 0, stringsAsFactors = FALSE))
}

#' Find Synonym Matches
#'
#' Internal function to match common modification synonyms
#'
#' @param mod_name Character, modification name
#' @param unimod_data UniMod data frame
#' @return Data frame with synonym matches
find_synonym_matches <- function(mod_name, unimod_data) {
  
  # Common synonyms mapping
  synonyms <- list(
    "oxidation" = c("oxidation", "ox"),
    "phosphorylation" = c("phospho", "phos", "phosphorylation"),
    "acetylation" = c("acetyl", "ace", "acetylation"),
    "methylation" = c("methyl", "me", "methylation"),
    "ubiquitination" = c("ubiquitin", "ub", "ubiquitination"),
    "carbamidomethyl" = c("carbamidomethyl", "cam", "iodoacetamide"),
    "propionamide" = c("propionamide", "acrylamide"),
    "dimethyl" = c("dimethyl", "dime", "dma"),
    "trimethyl" = c("trimethyl", "trime", "tma"),
    "deamidation" = c("deamidated", "deamidation"),
    "citrullination" = c("citrulline", "citrullination"),
    "nitrosylation" = c("nitrosyl", "nitrosylation"),
    "hydroxylation" = c("hydroxyl", "hydroxylation")
  )
  
  matches <- data.frame()
  mod_lower <- tolower(mod_name)
  
  for (canonical in names(synonyms)) {
    if (mod_lower %in% synonyms[[canonical]]) {
      # Find matches for canonical name
      canonical_matches <- unimod_data[
        grepl(canonical, tolower(unimod_data$full_name)) |
        grepl(canonical, tolower(unimod_data$code_name)), 
      ]
      matches <- rbind(matches, canonical_matches)
    }
  }
  
  return(matches)
}

#' Rank Modification Matches
#'
#' Internal function to rank matches by quality/relevance
#'
#' @param matches Data frame of matched modifications  
#' @param mod_name Original search term
#' @return Ranked matches with best matches first
rank_modification_matches <- function(matches, mod_name, mod_pos = NULL) {
  
  if (nrow(matches) == 0) return(matches)
  
  # Create quality score for each match
  matches$quality_score <- 0
  mod_name_lower <- tolower(mod_name)
  
  # HEAVILY penalize complex/labeled modifications
  matches$quality_score[grepl("label:|\\+|->|loss", tolower(matches$code_name))] <- -1000
  
  # Highest priority: exact name + position match
  if ("ex_code_name" %in% colnames(matches)) {
    matches$quality_score[tolower(matches$ex_code_name) == mod_name_lower] <- 2000
  }
  
  # HUGE boost for position match
  if (!is.null(mod_pos) && tolower(mod_pos) != "any") {
    if (nchar(mod_pos) == 1 && grepl("^[A-Z]$", toupper(mod_pos))) {
      # Single amino acid match
      matches$quality_score[matches$one_letter == toupper(mod_pos)] <- 
        matches$quality_score[matches$one_letter == toupper(mod_pos)] + 5000
    } else if (grepl("term", tolower(mod_pos))) {
      # Terminal position match
      matches$quality_score[grepl(mod_pos, matches$one_letter, ignore.case = TRUE)] <- 
        matches$quality_score[grepl(mod_pos, matches$one_letter, ignore.case = TRUE)] + 5000
    }
  }
  
  # Exact matches get high score
  matches$quality_score[tolower(matches$code_name) == mod_name_lower] <- 1500
  matches$quality_score[tolower(matches$full_name) == mod_name_lower] <- 1000
  
  # Prefer visible (non-hidden) modifications
  matches$quality_score[matches$hidden == 0] <- matches$quality_score[matches$hidden == 0] + 500
  
  # Matches that start with the search term
  matches$quality_score[grepl(paste0("^", mod_name_lower), tolower(matches$code_name))] <- 
    matches$quality_score[grepl(paste0("^", mod_name_lower), tolower(matches$code_name))] + 100
  
  # Shorter names are better (less verbose)
  matches$quality_score <- matches$quality_score + (50 - pmin(nchar(matches$code_name), 50))
  
  # Sort by quality score (highest first)
  matches <- matches[order(matches$quality_score, decreasing = TRUE), ]
  
  return(matches)
}

#' Filter by Position
#'
#' Internal function to filter modifications by amino acid position
#'
#' @param matches Data frame of modification matches
#' @param position Character, position specification
#' @return Filtered data frame
filter_by_position <- function(matches, position) {
  
  if (nrow(matches) == 0) return(matches)
  
  # Handle different position formats
  pos_lower <- tolower(position)
  
  # Single amino acid (e.g., "C", "M", "K")
  if (nchar(position) == 1 && grepl("^[A-Z]$", toupper(position))) {
    filtered <- matches[
      `|`(matches$one_letter == toupper(position),
      tolower(matches$position_key) == pos_lower), 
    ]
    if (nrow(filtered) > 0) return(filtered)
  }
  
  # Multiple amino acids (e.g., "STY", "KR")
  if (nchar(position) > 1 && grepl("^[A-Z]+$", toupper(position))) {
    aa_list <- strsplit(toupper(position), "")[[1]]
    for (aa in aa_list) {
      filtered <- matches[
        grepl(aa, matches$one_letter, fixed = TRUE), 
      ]
      if (nrow(filtered) > 0) return(filtered)
    }
  }
  
  # Terminal positions
  if (grepl("term", pos_lower)) {
    if (grepl("n.?term", pos_lower)) {
      filtered <- matches[grepl("N-term", matches$position), ]
    } else if (grepl("c.?term", pos_lower)) {
      filtered <- matches[grepl("C-term", matches$position), ]
    }
    if (exists("filtered") && nrow(filtered) > 0) return(filtered)
  }
  
  # If no specific filtering worked, return original matches
  return(matches)
}

#' Format for HiTMaP Workflow
#'
#' Internal function to format matched modifications for HiTMaP workflow
#'
#' @param match_results List with matched modifications
#' @param mod_type Character, "fixed" or "variable"
#' @return Properly formatted modification list
format_for_hitmap_workflow <- function(match_results, mod_type = "variable") {
  
  if (nrow(match_results$matched) == 0) {
    return(list(
      modifications = list(fixed = NULL, fixmod_position = NULL, 
                          variable = NULL, varmod_position = NULL),
      summary = "No modifications matched",
      ambiguous = match_results$ambiguous,
      not_found = match_results$not_found
    ))
  }
  
  # Extract modification names and positions
  # Use code_name for shorter, cleaner names
  mod_names <- match_results$matched$code_name
  
  mod_positions <- ifelse(
    is.na(match_results$matched$one_letter) | match_results$matched$one_letter == "",
    "any",
    match_results$matched$one_letter
  )
  
  # Create properly formatted modifications list
  if (mod_type == "fixed") {
    modifications <- list(
      fixed = mod_names,
      fixmod_position = mod_positions,
      fixmod_one_letter = mod_positions,
      variable = NULL,
      varmod_position = NULL,
      varmod_one_letter = NULL
    )
  } else {
    modifications <- list(
      fixed = NULL,
      fixmod_position = NULL,
      fixmod_one_letter = NULL,
      variable = mod_names,
      varmod_position = mod_positions,
      varmod_one_letter = mod_positions
    )
  }
  
  # Create summary
  summary_text <- paste(
    "Matched", length(mod_names), mod_type, "modifications:",
    paste(mod_names, collapse = ", ")
  )
  
  if (length(match_results$not_found) > 0) {
    summary_text <- paste(
      summary_text, "\n",
      "Not found:", paste(match_results$not_found, collapse = ", ")
    )
  }
  
  if (length(match_results$ambiguous) > 0) {
    summary_text <- paste(
      summary_text, "\n",
      "Ambiguous matches for:", paste(names(match_results$ambiguous), collapse = ", ")
    )
  }
  
  return(list(
    modifications = modifications,
    summary = summary_text,
    ambiguous = match_results$ambiguous,
    not_found = match_results$not_found
  ))
}

#' Quick Modification Setup
#'
#' Convenience function for common modification scenarios
#'
#' @param scenario Character, preset scenario ("standard", "comprehensive", "minimal", "custom")
#' @param custom_mods Character vector, custom modifications (used with scenario = "custom")
#' @param custom_positions Character vector, custom positions (used with scenario = "custom")
#' @param update_unimod Logical, whether to update UniMod database
#' 
#' @return Properly formatted modification list for HiTMaP workflow
#'
#' @examples
#' # Standard proteomics modifications
#' mods <- quick_modification_setup("standard")
#' 
#' # Custom modifications
#' mods <- quick_modification_setup("custom", 
#'   custom_mods = c("Oxidation", "Phospho"), 
#'   custom_positions = c("M", "STY")
#' )
#'
#' @export
quick_modification_setup <- function(scenario = "standard", 
                                   custom_mods = NULL, 
                                   custom_positions = NULL,
                                   update_unimod = FALSE) {
  
  if (scenario == "standard") {
    # Standard MALDI imaging setup - typically no fixed Cys modification
    var_mods <- format_unimod_modifications(
      modifications = c("Oxidation", "Acetyl"),
      positions = c("M", "N-term"),
      mod_type = "variable",
      update_unimod = update_unimod,
      interactive = FALSE
    )
    
    return(list(
      modifications = list(
        fixed = NULL,
        fixmod_position = NULL,
        variable = var_mods$modifications$variable,
        varmod_position = var_mods$modifications$varmod_position
      ),
      summary = paste("Standard MALDI imaging modifications:", var_mods$summary)
    ))
    
  } else if (scenario == "comprehensive") {
    # Comprehensive MALDI imaging PTM discovery - no fixed modifications
    var_mods <- format_unimod_modifications(
      modifications = c("Oxidation", "Acetyl", "Phospho", "Deamidation", "Methylation", "Nitrosylation"),
      positions = c("M", "N-term", "STY", "NQ", "KR", "C"),
      mod_type = "variable",
      update_unimod = update_unimod,
      interactive = FALSE
    )
    
    return(list(
      modifications = list(
        fixed = NULL,
        fixmod_position = NULL,
        variable = var_mods$modifications$variable,
        varmod_position = var_mods$modifications$varmod_position
      ),
      summary = paste("Comprehensive MALDI imaging PTM discovery:", var_mods$summary)
    ))
    
  } else if (scenario == "minimal") {
    return(list(
      modifications = list(fixed = NULL, fixmod_position = NULL, 
                          variable = NULL, varmod_position = NULL),
      summary = "No modifications (minimal scenario)"
    ))
    
  } else if (scenario == "custom") {
    if (is.null(custom_mods)) {
      stop("custom_mods must be provided for custom scenario")
    }
    
    return(format_unimod_modifications(
      modifications = custom_mods,
      positions = custom_positions,
      mod_type = "variable",
      update_unimod = update_unimod,
      interactive = FALSE
    ))
    
  } else {
    stop("Scenario must be one of: 'standard', 'comprehensive', 'minimal', 'custom'")
  }
}

#' Print Modification Summary
#'
#' Print a formatted summary of modifications
#'
#' @param mod_result Result from format_unimod_modifications or quick_modification_setup
#' @export
print_modification_summary <- function(mod_result) {
  cat("=== Modification Summary ===\n")
  cat(mod_result$summary, "\n")
  
  if (length(mod_result$ambiguous) > 0) {
    cat("\n=== Ambiguous Matches ===\n")
    for (name in names(mod_result$ambiguous)) {
      cat(paste("For", name, ":\n"))
      matches <- mod_result$ambiguous[[name]]
      for (i in 1:nrow(matches)) {
        cat(sprintf("  %d: %s (%s) at %s\n", i, 
                   matches$full_name[i], 
                   matches$code_name[i], 
                   matches$position[i]))
      }
    }
  }
  
  if (length(mod_result$not_found) > 0) {
    cat("\n=== Not Found ===\n")
    cat(paste(mod_result$not_found, collapse = ", "), "\n")
  }
  
  cat("\n=== HiTMaP Format ===\n")
  cat("Use the following in your imaging_identification call:\n")
  cat("Modifications = mod_result$modifications\n")
}