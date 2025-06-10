#' Modular Imaging Identification Workflow Functions
#' 
#' This file contains modular functions for step-by-step control of the imaging identification workflow.
#' Each function represents a specific processing step with clear inputs and outputs.

#' Step 1: Initialize and Setup Project
#' 
#' @param datafile Path(s) to the imzML data file(s)
#' @param projectfolder Optional project folder path
#' @param Thread Number of threads for parallel processing
#' @return List containing workdir, datafile names, and BPPARAM
#' 
#' @export
hitmap_init_project <- function(datafile, projectfolder = NULL, Thread = 4) {
  suppressMessages(suppressWarnings(library("pacman")))
  suppressMessages(suppressWarnings(p_load(stringr, BiocParallel, data.table, Cardinal, parallel)))
  
  if (missing(datafile)) stop("Missing data file, Choose single or multiple imzml file(s) for analysis")
  
  # Retrieve/parse the working dir info, and convert the filenames
  if (is.null(projectfolder)) {
    workdir <- base::dirname(datafile[1])
  } else { 
    workdir <- projectfolder 
  }
  
  datafile <- basename(datafile)
  datafile <- gsub(".imzML$", "", datafile)
  datafile_imzML <- paste0(datafile, ".imzML")
  
  setwd(paste0(workdir[1], "/"))
  
  # Set the parallel processing parameter
  if (is.null(Thread)) {
    parallel <- try(future::availableCores() / 2)
    if (parallel < 1 | is.null(parallel)) { parallel <- 1 }
    BPPARAM <- HiTMaP:::Parallel.OS(parallel)
    setCardinalBPPARAM(BPPARAM = BPPARAM)
  } else {
    parallel <- Thread
    BPPARAM <- HiTMaP:::Parallel.OS(parallel)
    setCardinalBPPARAM(BPPARAM = BPPARAM)
  }
  
  message(paste(try(future::availableCores()), "Cores detected,", parallel, "threads will be used for computing"))
  message(paste(length(datafile), "files were selected and will be used for Searching"))
  
  return(list(
    workdir = workdir,
    datafile = datafile,
    datafile_imzML = datafile_imzML,
    parallel = parallel,
    BPPARAM = BPPARAM
  ))
}

#' Step 2: Generate Protein Feature List (Candidate Generation)
#' 
#' @param workdir Working directory path
#' @param database Fasta database filename
#' @param Digestion_site Enzyme digestion site specification
#' @param missedCleavages Number of missed cleavages allowed
#' @param adducts List of adducts to consider
#' @param Modifications List of modifications
#' @param Substitute_AA Amino acid substitutions
#' @param Decoy_search Enable decoy search
#' @param Decoy_adducts Decoy adducts list
#' @param Decoy_mode Decoy search mode
#' @param mzrange m/z range for analysis
#' @param Protein_desc_of_exclusion Proteins to exclude
#' @param Database_stats Generate database statistics
#' @param output_candidatelist Output candidate list
#' @param use_previous_candidates Use existing candidate list
#' @param BPPARAM BiocParallel parameter
#' @return Protein feature list data frame
#' 
#' @export
hitmap_generate_candidates <- function(workdir, 
                                     database = "uniprot-bovin.fasta",
                                     Digestion_site = "trypsin",
                                     missedCleavages = 0:1,
                                     adducts = c("M+H"),
                                     Modifications = list(fixed = NULL, fixmod_position = NULL, variable = NULL, varmod_position = NULL),
                                     Substitute_AA = NULL,
                                     Decoy_search = TRUE,
                                     Decoy_adducts = c("M+ACN+H", "M+IsoProp+H", "M+DMSO+H", "M+Co", "M+Ag", "M+Cu", "M+He", "M+Ne", "M+Ar", "M+Kr", "M+Xe", "M+Rn"),
                                     Decoy_mode = "isotope",
                                     mzrange = c(700, 4000),
                                     Protein_desc_of_exclusion = NULL,
                                     Database_stats = FALSE,
                                     output_candidatelist = TRUE,
                                     use_previous_candidates = FALSE,
                                     BPPARAM = bpparam()) {
  
  message(paste(database, "was selected as database. Candidates will be generated through Proteomics mode"))
  
  Protein_feature_list <- Protein_feature_list_fun(
    workdir = workdir,
    database = database,
    Digestion_site = Digestion_site,
    missedCleavages = missedCleavages,
    adducts = adducts,
    BPPARAM = BPPARAM,
    Decoy_adducts = Decoy_adducts,
    Decoy_search = Decoy_search,
    Decoy_mode = Decoy_mode,
    output_candidatelist = output_candidatelist,
    use_previous_candidates = use_previous_candidates,
    Substitute_AA = Substitute_AA,
    Modifications = Modifications,
    mzrange = mzrange,
    Protein_desc_of_exclusion = Protein_desc_of_exclusion,
    Database_stats = Database_stats
  )
  
  return(Protein_feature_list)
}

#' Step 3: Preprocess and Segment IMS Data
#' 
#' @param datafile Single data file name
#' @param workdir Working directory
#' @param ppm Mass tolerance in ppm
#' @param import_ppm Import mass tolerance
#' @param mzrange m/z range for analysis
#' @param segmentation_num Number of segments
#' @param Segmentation Segmentation method
#' @param Segmentation_def Segmentation definition file
#' @param Segmentation_ncomp Number of components for segmentation
#' @param Segmentation_variance_coverage Variance coverage for segmentation
#' @param Smooth_range Smoothing range
#' @param Virtual_segmentation_rankfile Virtual segmentation rank file
#' @param rotate Image rotation parameters
#' @param preprocess Preprocessing parameters list
#' @param BPPARAM BiocParallel parameters
#' @return List containing segmented imdata and segmentation labels
#' 
#' @export
hitmap_preprocess_segment <- function(datafile,
                                    workdir,
                                    ppm = 5,
                                    import_ppm = 5,
                                    mzrange = "auto-detect",
                                    segmentation_num = 4,
                                    Segmentation = "spatialKMeans",
                                    Segmentation_def = "segmentation_def.csv",
                                    Segmentation_ncomp = "auto-detect",
                                    Segmentation_variance_coverage = 0.8,
                                    Smooth_range = 1,
                                    Virtual_segmentation_rankfile = NULL,
                                    rotate = NULL,
                                    preprocess = list(
                                      force_preprocess = FALSE,
                                      use_preprocessRDS = TRUE,
                                      smoothSignal = list(method = "disable"),
                                      reduceBaseline = list(method = "locmin"),
                                      peakPick = list(method = "adaptive"),
                                      peakAlign = list(tolerance = ppm/2, units = "ppm"),
                                      peakFilter = list(freq.min = 0.05),
                                      normalize = list(method = c("rms", "tic", "reference")[1], mz = 1)
                                    ),
                                    BPPARAM = bpparam()) {
  
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(Cardinal)))
  suppressMessages(suppressWarnings(require(stringr)))
  
  setCardinalBPPARAM(BPPARAM)
  
  # Resolve the filename and rotation info
  rotate <- Parse_rotation(datafile, rotate)
  
  message(paste("Preprocessing and segmenting:", datafile))
  message(paste("Segmentation method:", Segmentation, "| Segments:", segmentation_num))
  
  # Perform IMS pre-processing and image segmentation
  setwd(workdir)
  segmentation_res <- Preprocessing_segmentation(
    datafile = datafile,
    workdir = workdir,
    segmentation_num = segmentation_num,
    ppm = ppm, 
    import_ppm = import_ppm,
    Bypass_Segmentation = FALSE,
    mzrange = mzrange,
    Segmentation = Segmentation,
    Segmentation_def = Segmentation_def,
    Segmentation_ncomp = Segmentation_ncomp,
    Segmentation_variance_coverage = Segmentation_variance_coverage,
    Smooth_range = Smooth_range,
    colorstyle = "Set1",
    Virtual_segmentation_rankfile = Virtual_segmentation_rankfile,
    rotate = rotate,
    BPPARAM = BPPARAM,
    preprocess = preprocess
  )
  
  return(list(
    imdata = segmentation_res$imdata,
    segmentation_label = segmentation_res$segmentation_label,
    datafile = datafile,
    workdir = workdir
  ))
}

#' Step 4: Perform PMF Search on Single Region
#' 
#' @param segmented_data Output from hitmap_preprocess_segment
#' @param region_name Name of the region to process
#' @param candidate_list Protein feature list from hitmap_generate_candidates
#' @param threshold Intensity threshold
#' @param ppm Mass tolerance in ppm
#' @param TopNFeatures Maximum number of features to use
#' @param score_method Scoring method
#' @param Decoy_mode Decoy search mode
#' @param Decoy_search Enable decoy search
#' @param adjust_score Adjust scores
#' @param FDR_cutoff FDR cutoff threshold
#' @param peptide_ID_filter Minimum peptides for protein ID
#' @param use_top_rank Use top ranking peptides only
#' @param plot_matching_score Plot matching scores
#' @param preprocess Preprocessing parameters
#' @param BPPARAM BiocParallel parameters
#' @return List containing PMF search results for the region
#' 
#' @export
hitmap_pmf_search_region <- function(segmented_data,
                                   region_name,
                                   candidate_list,
                                   threshold = 0.001,
                                   ppm = 5,
                                   TopNFeatures = 1250,
                                   score_method = "SQRTP",
                                   Decoy_mode = "isotope",
                                   Decoy_search = TRUE,
                                   adjust_score = FALSE,
                                   FDR_cutoff = 0.05,
                                   peptide_ID_filter = 2,
                                   use_top_rank = NULL,
                                   plot_matching_score = FALSE,
                                   preprocess = list(peakAlign = list(method = "Disable", tolerance = 0, units = "ppm")),
                                   BPPARAM = bpparam()) {
  
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(ggplot2)))
  
  imdata <- segmented_data$imdata
  segmentation_label <- segmented_data$segmentation_label
  datafile <- segmented_data$datafile
  workdir <- segmented_data$workdir
  
  if (!region_name %in% names(segmentation_label)) {
    stop(paste("Region", region_name, "not found in segmentation labels"))
  }
  
  message(paste("Processing PMF search for region:", region_name))
  
  # Create output directory
  output_dir <- paste0(workdir, "/", gsub(".imzML$", "", datafile), " ID/", region_name, "/")
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  # Extract region-specific data
  imdata_sb <- imdata[, unlist(segmentation_label[[region_name]])]
  imdata_ed <- imdata_sb
  
  # Apply peak alignment if specified
  if (preprocess$peakAlign$method != "Disable" && !is.null(preprocess$peakAlign$tolerance) && preprocess$peakAlign$tolerance > 0) {
    message("Applying peak alignment with tolerance:", preprocess$peakAlign$tolerance, preprocess$peakAlign$units)
    imdata_ed <- imdata_ed %>% peakAlign(tolerance = preprocess$peakAlign$tolerance, units = preprocess$peakAlign$units)
    imdata_ed <- imdata_ed %>% process()
  }
  
  # Generate spectrum for the region
  spectrum_file_table <- summarizeFeatures(imdata_ed, "mean")
  spectrum_file_table <- data.frame(
    mz = spectrum_file_table@featureData@listData[["mz"]],
    mean = spectrum_file_table@featureData@listData[["mean"]]
  )
  peaklist <- spectrum_file_table
  colnames(peaklist) <- c("m.z", "intensities")
  
  # Filter peaklist
  peaklist <- peaklist[peaklist$intensities > 0, ]
  write.csv(peaklist, paste0(output_dir, "Spectrum.csv"), row.names = FALSE)
  
  # Generate filtered processed peaklist for PMF
  deconv_peaklist <- HiTMaP:::isopattern_ppm_filter_peaklist(peaklist, ppm = ppm, threshold = 0)
  deconv_peaklist_thres_id <- deconv_peaklist$intensities >= max(deconv_peaklist$intensities) * threshold
  deconv_peaklist_topN_id <- rank(deconv_peaklist$intensities) > (max(rank(deconv_peaklist$intensities)) - TopNFeatures)
  deconv_peaklist_PMF_id <- deconv_peaklist_thres_id | deconv_peaklist_topN_id
  peaklist_pmf <- deconv_peaklist[deconv_peaklist_PMF_id, ]
  
  message(paste(nrow(peaklist_pmf), "mz features subjected to PMF search"))
  
  if (nrow(peaklist_pmf) == 0) {
    warning("No features found above threshold for PMF search")
    return(list(
      region = region_name,
      peptides = data.frame(),
      proteins = data.frame(),
      peaklist = peaklist,
      output_dir = output_dir
    ))
  }
  
  # Perform first round of peptide search
  Peptide_Summary_searchlist <- unique(candidate_list)
  mz_feature_list <- Do_PMF_search(peaklist_pmf, Peptide_Summary_searchlist, BPPARAM = BPPARAM, ppm = ppm)
  mz_feature_list <- unique(mz_feature_list)
  
  # Merge results
  Peptide_Summary_searchlist <- as.data.table(Peptide_Summary_searchlist)
  mz_feature_list <- as.data.table(mz_feature_list)
  Peptide_Summary_searchlist$Intensity <- NULL
  mz_feature_list$mz <- as.character(mz_feature_list$mz)
  Peptide_Summary_searchlist$mz <- as.character(Peptide_Summary_searchlist$mz)
  Peptide_Summary_searchlist <- merge(Peptide_Summary_searchlist, mz_feature_list, by.x = "mz", by.y = "mz", all.x = TRUE, sort = FALSE)
  Peptide_Summary_searchlist$Intensity[is.na(Peptide_Summary_searchlist$Intensity)] <- 0
  Peptide_plot_list <- Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity > 0, ]
  
  if (nrow(Peptide_plot_list) == 0) {
    warning("No peptide matches found")
    return(list(
      region = region_name,
      peptides = data.frame(),
      proteins = data.frame(),
      peaklist = peaklist,
      output_dir = output_dir
    ))
  }
  
  # Score peptides
  data(isotopes, package = "enviPat")
  unique_formula <- unique(Peptide_plot_list$formula)
  Peptide_plot_list_Score <- bplapply(unique_formula, SCORE_PMF, 
                                    peaklist = peaklist, isotopes = isotopes, 
                                    score_method = score_method, charge = 1, ppm = ppm, 
                                    BPPARAM = BPPARAM)
  
  Peptide_plot_list_Score_m <- as.data.frame(do.call(rbind, Peptide_plot_list_Score))
  names(Peptide_plot_list_Score_m) <- c("Score", "delta_ppm", "Intensity")
  formula_score <- data.frame(
    formula = unique_formula,
    Score = Peptide_plot_list_Score_m$Score,
    Delta_ppm = Peptide_plot_list_Score_m$delta_ppm,
    Intensity = Peptide_plot_list_Score_m$Intensity
  )
  
  # Merge scores back
  Peptide_plot_list$Score <- NULL
  Peptide_plot_list$Delta_ppm <- NULL
  Peptide_plot_list$Intensity <- NULL
  Peptide_plot_list <- merge(Peptide_plot_list, formula_score, by = "formula")
  Peptide_plot_list$Region <- region_name
  
  # Rank peptides
  Peptide_plot_list$mz <- as.numeric(as.character(Peptide_plot_list$mz))
  Peptide_plot_list_rank <- rank_mz_feature(Peptide_plot_list, mz_feature = deconv_peaklist, BPPARAM = BPPARAM)
  Peptide_plot_list_rank <- unique(Peptide_plot_list_rank)
  
  # Apply FDR cutoff
  Score_cutoff <- FDR_cutoff_plot(Peptide_plot_list_rank, FDR_cutoff = FDR_cutoff, plot_fdr = TRUE, 
                                 outputdir = output_dir, adjust_score = adjust_score)
  if (adjust_score == FALSE) {
    Peptide_plot_list_2nd <- Peptide_plot_list_rank[((Peptide_plot_list_rank$Score >= Score_cutoff) & (!is.na(Peptide_plot_list_rank$Intensity))), ]
  } else {
    Peptide_plot_list_2nd <- Score_cutoff[[2]]
    Score_cutoff <- Score_cutoff[[1]]
    Peptide_plot_list_2nd <- Peptide_plot_list_2nd[((Peptide_plot_list_2nd$Score >= Score_cutoff) & (!is.na(Peptide_plot_list_2nd$Intensity))), ]
  }
  
  # Save intermediate results
  write.csv(Peptide_plot_list_rank, paste0(output_dir, "Peptide_scored_ranked.csv"), row.names = FALSE)
  write.csv(Peptide_plot_list_2nd, paste0(output_dir, "Peptide_filtered.csv"), row.names = FALSE)
  
  # Protein scoring if we have results
  protein_results <- data.frame()
  if (nrow(Peptide_plot_list_2nd) > 0) {
    # Get full protein feature list for scoring
    if (exists("Protein_feature_list", envir = .GlobalEnv)) {
      Protein_feature_list <- get("Protein_feature_list", envir = .GlobalEnv)
    } else {
      Protein_feature_list <- candidate_list
    }
    
    Protein_feature_result <- protein_scoring(Protein_feature_list, Peptide_plot_list_rank, 
                                            BPPARAM = BPPARAM, scoretype = "mean", 
                                            peptide_ID_filter = peptide_ID_filter, 
                                            use_top_rank = use_top_rank)
    
    protein_scores <- Protein_feature_result[[1]]
    protein_details <- Protein_feature_result[[2]]
    
    if (nrow(protein_scores) >= 2) {
      Score_cutoff_protein <- FDR_cutoff_plot_protein(protein_scores, FDR_cutoff = FDR_cutoff, 
                                                    plot_fdr = TRUE, outputdir = output_dir, 
                                                    adjust_score = FALSE)
    } else {
      Score_cutoff_protein <- ifelse(nrow(protein_scores) > 0, protein_scores$Proscore[1], 0)
    }
    
    protein_results <- protein_scores[((protein_scores$Proscore >= Score_cutoff_protein) & 
                                     (!is.na(protein_scores$Intensity)) & 
                                     (protein_scores$isdecoy == 0)), ]
    
    # Save protein results
    write.csv(protein_scores, paste0(output_dir, "Protein_scored.csv"), row.names = FALSE)
    write.csv(protein_results, paste0(output_dir, "Protein_filtered.csv"), row.names = FALSE)
  }
  
  return(list(
    region = region_name,
    peptides = Peptide_plot_list_2nd,
    proteins = protein_results,
    peaklist = peaklist,
    output_dir = output_dir,
    score_cutoff = Score_cutoff
  ))
}

#' Step 5: Process All Regions for a Single Data File
#' 
#' @param segmented_data Output from hitmap_preprocess_segment
#' @param candidate_list Protein feature list
#' @param threshold Intensity threshold
#' @param ppm Mass tolerance in ppm
#' @param score_method Scoring method
#' @param FDR_cutoff FDR cutoff threshold
#' @param peptide_ID_filter Minimum peptides for protein ID
#' @param plot_matching_score Plot matching scores
#' @param BPPARAM BiocParallel parameters
#' @param ... Additional parameters passed to hitmap_pmf_search_region
#' @return List containing results for all regions
#' 
#' @export
hitmap_process_all_regions <- function(segmented_data,
                                     candidate_list,
                                     threshold = 0.001,
                                     ppm = 5,
                                     score_method = "SQRTP",
                                     FDR_cutoff = 0.05,
                                     peptide_ID_filter = 2,
                                     plot_matching_score = FALSE,
                                     BPPARAM = bpparam(),
                                     ...) {
  
  region_names <- names(segmented_data$segmentation_label)
  datafile <- segmented_data$datafile
  workdir <- segmented_data$workdir
  
  message(paste("Processing", length(region_names), "regions for file:", datafile))
  message(paste("Regions found:", paste(region_names, collapse = ", ")))
  
  # Clean up previous results
  output_base <- paste0(workdir, "/", gsub(".imzML$", "", datafile), " ID")
  if (dir.exists(output_base)) {
    all_files <- dir(output_base, recursive = FALSE)
    delete_files <- all_files[stringr::str_detect(all_files, "Protein_segment_PMF_RESULT_|Peptide_|PMF spectrum match.png$")]
    delete_dir <- list.dirs(path = output_base, full.names = TRUE, recursive = FALSE)
    try(suppressWarnings(file.remove(paste0(output_base, "/", delete_files))))
    unlink(delete_dir, recursive = TRUE)
  }
  
  # Process each region
  all_results <- list()
  combined_peptides <- data.frame()
  
  for (region in region_names) {
    result <- hitmap_pmf_search_region(
      segmented_data = segmented_data,
      region_name = region,
      candidate_list = candidate_list,
      threshold = threshold,
      ppm = ppm,
      score_method = score_method,
      FDR_cutoff = FDR_cutoff,
      peptide_ID_filter = peptide_ID_filter,
      plot_matching_score = plot_matching_score,
      BPPARAM = BPPARAM,
      ...
    )
    
    all_results[[region]] <- result
    
    # Combine peptide results
    if (nrow(result$peptides) > 0) {
      combined_peptides <- rbind(combined_peptides, result$peptides)
    }
    
    # Save region-specific protein results
    if (nrow(result$proteins) > 0) {
      write.csv(result$proteins, paste0(output_base, "/Protein_segment_PMF_RESULT_", region, ".csv"), row.names = FALSE)
    }
  }
  
  # Save combined peptide results
  setwd(output_base)
  write.csv(combined_peptides, "Peptide_region_file.csv", row.names = FALSE)
  
  setwd(workdir)
  
  return(list(
    datafile = datafile,
    regions = all_results,
    combined_peptides = combined_peptides,
    output_dir = output_base
  ))
}

#' Step 6: Summarize Results Across Multiple Files
#' 
#' @param project_setup Output from hitmap_init_project
#' @param Protein_feature_summary Generate protein summary
#' @param Peptide_feature_summary Generate peptide summary
#' @param Region_feature_summary Generate region feature summary
#' @return List of summary data frames
#' 
#' @export
hitmap_summarize_results <- function(project_setup,
                                   Protein_feature_summary = TRUE,
                                   Peptide_feature_summary = TRUE,
                                   Region_feature_summary = FALSE) {
  
  workdir <- project_setup$workdir
  datafile <- project_setup$datafile
  
  # Create summary folder
  summary_dir <- paste(workdir, "/Summary folder", sep = "")
  if (!dir.exists(summary_dir)) { dir.create(summary_dir) }
  
  results <- list()
  
  # Protein feature summary
  if (Protein_feature_summary) {
    message("Generating protein feature summary...")
    protein_feature_all <- NULL
    Protein_peptide_Summary_file <- NULL
    
    for (i in 1:length(datafile)) {
      datafilename <- gsub(paste(workdir, "/", sep = ""), "", gsub(".imzML", "", datafile[i]))
      currentdir <- paste0(workdir, "/", datafile[i], " ID")
      
      if (dir.exists(currentdir)) {
        setwd(currentdir)
        
        # Collect protein results
        for (protein_feature_file in dir()[stringr::str_detect(dir(), "Protein_segment_PMF_RESULT_")]) {
          protein_feature <- fread(protein_feature_file)
          
          if (nrow(protein_feature) != 0) {
            protein_feature$Source <- datafilename
            region_code <- stringr::str_replace(protein_feature_file, "Protein_segment_PMF_RESULT_", "")
            region_code <- stringr::str_replace(region_code, ".csv", "")
            protein_feature$Region <- region_code
            protein_feature_all <- rbind(protein_feature_all, protein_feature)
          }
        }
        
        # Collect peptide results
        if (file.exists("Peptide_region_file.csv")) {
          Peptide_Summary_file <- fread("Peptide_region_file.csv")
          Peptide_Summary_file$Source <- datafilename
          if (nrow(Peptide_Summary_file) != 0) {
            Protein_peptide_Summary_file <- rbind(Protein_peptide_Summary_file, Peptide_Summary_file)
          }
        }
      }
    }
    
    # Write protein summaries
    if (!is.null(protein_feature_all)) {
      write.csv(protein_feature_all, paste(summary_dir, "/Protein_Summary.csv", sep = ""), row.names = FALSE)
      results$Protein_Summary <- protein_feature_all
    }
    
    if (!is.null(Protein_peptide_Summary_file)) {
      write.csv(Protein_peptide_Summary_file, paste(summary_dir, "/Protein_peptide_Summary.csv", sep = ""), row.names = FALSE)
      results$Protein_peptide_Summary <- Protein_peptide_Summary_file
    }
    
    message("Protein feature summary completed.")
  }
  
  # Peptide feature summary
  if (Peptide_feature_summary) {
    message("Generating peptide feature summary...")
    Peptide_Summary_file_a <- NULL
    
    for (i in 1:length(datafile)) {
      currentdir <- paste0(workdir, "/", datafile[i], " ID")
      
      if (dir.exists(currentdir)) {
        setwd(currentdir)
        
        if (file.exists("Peptide_region_file.csv")) {
          Peptide_Summary_file <- fread("Peptide_region_file.csv")
          Peptide_Summary_file$Source <- gsub(".imzML", "", datafile[i])
          if (nrow(Peptide_Summary_file) != 0) {
            Peptide_Summary_file_a <- rbind(Peptide_Summary_file_a, Peptide_Summary_file)
          }
        }
      }
    }
    
    if (!is.null(Peptide_Summary_file_a)) {
      Peptide_Summary_file_a <- unique(Peptide_Summary_file_a)
      write.csv(Peptide_Summary_file_a, paste(summary_dir, "/Peptide_Summary.csv", sep = ""), row.names = FALSE)
      results$Peptide_Summary <- Peptide_Summary_file_a
    }
    
    message("Peptide feature summary completed.")
  }
  
  # Region feature summary
  if (Region_feature_summary) {
    message("Generating region feature summary...")
    Spectrum_summary <- NULL
    
    for (i in 1:length(datafile)) {
      currentdir <- paste0(workdir, "/", datafile[i], " ID")
      
      if (dir.exists(currentdir)) {
        setwd(currentdir)
        match_pattern <- "Spectrum.csv"
        spectrum_file_table_sum <- NULL
        
        for (spectrum_file in dir(recursive = TRUE)[stringr::str_detect(dir(recursive = TRUE), match_pattern)]) {
          spectrum_file_table <- fread(spectrum_file)
          if (nrow(spectrum_file_table) >= 1) {
            spectrum_file_table$Region <- gsub("/Spectrum.csv", "", spectrum_file)
            spectrum_file_table$Source <- gsub(".imzML", "", datafile[i])
            spectrum_file_table_sum[[spectrum_file]] <- spectrum_file_table
          }
        }
        
        if (!is.null(spectrum_file_table_sum)) {
          spectrum_file_table_sum <- do.call(rbind, spectrum_file_table_sum)
          Spectrum_summary[[datafile[i]]] <- spectrum_file_table_sum
        }
      }
    }
    
    if (!is.null(Spectrum_summary)) {
      Spectrum_summary <- do.call(rbind, Spectrum_summary)
      write.csv(Spectrum_summary, file = paste(summary_dir, "/Region_feature_summary.csv", sep = ""), row.names = FALSE)
      results$Region_feature_summary <- Spectrum_summary
    }
    
    message("Region feature summary completed.")
  }
  
  setwd(workdir)
  
  return(results)
}

#' Step 7: Complete Modular Workflow Function
#' 
#' @param datafile Path(s) to data files
#' @param projectfolder Project folder path
#' @param database Database filename
#' @param threshold Intensity threshold
#' @param ppm Mass tolerance
#' @param score_method Scoring method
#' @param FDR_cutoff FDR cutoff
#' @param peptide_ID_filter Peptide ID filter
#' @param Thread Number of threads
#' @param segmentation_num Number of segments
#' @param Segmentation Segmentation method
#' @param preprocess Preprocessing parameters
#' @param Protein_feature_summary Generate protein summary
#' @param Peptide_feature_summary Generate peptide summary
#' @param Region_feature_summary Generate region summary
#' @param ... Additional parameters
#' @return Complete workflow results
#' 
#' @export
hitmap_modular_workflow <- function(datafile,
                                  projectfolder = NULL,
                                  database = "uniprot-bovin.fasta",
                                  threshold = 0.001,
                                  ppm = 5,
                                  score_method = "SQRTP",
                                  FDR_cutoff = 0.05,
                                  peptide_ID_filter = 2,
                                  Thread = 4,
                                  segmentation_num = 4,
                                  Segmentation = "spatialKMeans",
                                  preprocess = list(
                                    force_preprocess = FALSE,
                                    use_preprocessRDS = TRUE,
                                    smoothSignal = list(method = "disable"),
                                    reduceBaseline = list(method = "locmin"),
                                    peakPick = list(method = "adaptive"),
                                    peakAlign = list(tolerance = ppm/2, units = "ppm"),
                                    peakFilter = list(freq.min = 0.05),
                                    normalize = list(method = "rms", mz = 1)
                                  ),
                                  Protein_feature_summary = TRUE,
                                  Peptide_feature_summary = TRUE,
                                  Region_feature_summary = FALSE,
                                  ...) {
  
  message("Starting HiTMaP modular workflow...")
  
  # Step 1: Initialize project
  project_setup <- hitmap_init_project(datafile, projectfolder, Thread)
  
  # Step 2: Generate candidates
  message("\n=== Step 2: Generating candidate list ===")
  candidate_list <- hitmap_generate_candidates(
    workdir = project_setup$workdir,
    database = database,
    BPPARAM = project_setup$BPPARAM,
    ...
  )
  
  # Step 3-5: Process each file
  all_file_results <- list()
  
  for (i in 1:length(project_setup$datafile)) {
    current_file <- project_setup$datafile[i]
    message(paste("\n=== Processing file", i, "of", length(project_setup$datafile), ":", current_file, "==="))
    
    # Step 3: Preprocess and segment
    message("Step 3: Preprocessing and segmentation...")
    segmented_data <- hitmap_preprocess_segment(
      datafile = current_file,
      workdir = project_setup$workdir,
      ppm = ppm,
      segmentation_num = segmentation_num,
      Segmentation = Segmentation,
      preprocess = preprocess,
      BPPARAM = project_setup$BPPARAM,
      ...
    )
    
    # Step 4-5: Process all regions
    message("Step 4-5: PMF search across all regions...")
    file_results <- hitmap_process_all_regions(
      segmented_data = segmented_data,
      candidate_list = candidate_list,
      threshold = threshold,
      ppm = ppm,
      score_method = score_method,
      FDR_cutoff = FDR_cutoff,
      peptide_ID_filter = peptide_ID_filter,
      BPPARAM = project_setup$BPPARAM,
      ...
    )
    
    all_file_results[[current_file]] <- file_results
  }
  
  # Step 6: Summarize results
  message("\n=== Step 6: Summarizing results across all files ===")
  summary_results <- hitmap_summarize_results(
    project_setup = project_setup,
    Protein_feature_summary = Protein_feature_summary,
    Peptide_feature_summary = Peptide_feature_summary,
    Region_feature_summary = Region_feature_summary
  )
  
  message("\nHiTMaP modular workflow completed successfully!")
  
  return(list(
    project_setup = project_setup,
    candidate_list = candidate_list,
    file_results = all_file_results,
    summary_results = summary_results
  ))
}