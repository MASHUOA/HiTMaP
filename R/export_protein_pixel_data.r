#' export_protein_pixel_data
#'
#' This function exports protein scores at pixel level and visualizes protein distributions in MALDI imaging data
#' 
#' @param datafile The data files' path for the analysis, leave it as blank to enable a GUI to select data
#' @param projectfolder Optional folder containing project data, if NULL will extract the path from datafile(s)
#' @param protein_file Path to the protein results file (CSV) with protein scores and peptide mappings
#' @param ppm Mass tolerance (in ppm) for peak integration
#' @param smooth.image Smoothing method for image processing ("gaussian", "adaptive", or "none")
#' @param contrast.enhance Contrast enhancement method ("suppression", "none")
#' @param plot_style Visualization style ("fleximaging", "ClusterOnly", "rainbow")
#' @param Component_plot_coloure Color scheme for visualization ("mono", "as.cluster")
#' @param protein_list Character vector of protein IDs to visualize. If NULL, will use top scoring proteins
#' @param max_proteins Maximum number of proteins to visualize (default: 10)
#' @param workdir Working directory to save output files
#' @param pixelwise_normalization Method for normalizing pixel intensities (default: "sum")
#' @param export_data Whether to export raw protein intensity data as CSV (default: TRUE)
#'
#' @return List containing protein distribution images and pixel-level data
#'
#' @examples
#' export_protein_pixel_data(datafile="Mouse_brain_trimmed", 
#'                           protein_file="Summary folder/Protein_peptide_Summary.csv", 
#'                           ppm=10)
#'
#' @export
export_protein_pixel_data <- function(datafile = NULL,
                                      projectfolder = NULL,
                                      protein_file = NULL,
                                      ppm = 5,
                                      smooth.image = "gaussian",
                                      contrast.enhance = "suppression",
                                      plot_style = "fleximaging",
                                      Component_plot_coloure = "mono",
                                      protein_list = NULL,
                                      max_proteins = 10,
                                      workdir = getwd(),
                                      pixelwise_normalization = "sum",
                                      export_data = TRUE) {
  
  # Load required packages
  suppressMessages(suppressWarnings(require(stringr)))
  suppressMessages(suppressWarnings(require(Cardinal)))
  suppressMessages(suppressWarnings(require(magick)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(dplyr)))
  
  # Set up working directory structure
  if (is.null(projectfolder) && !is.null(datafile)) {
    projectfolder <- dirname(datafile[1])
  }
  
  # Convert filenames
  if (!is.null(datafile)) {
    datafile <- basename(datafile)
    datafile <- gsub(".imzML$", "", datafile)
  } else {
    stop("Missing data file. Please provide imzML file for analysis.")
  }
  
  # Create output directory for protein images
  protein_output_dir <- file.path(workdir, "protein_score_images")
  if (!dir.exists(protein_output_dir)) {
    dir.create(protein_output_dir, recursive = TRUE)
  }
  
  # Load protein data
  if (is.null(protein_file)) {
    protein_file <- file.path(workdir, paste0(datafile, " ID/Protein_peptide_Summary.csv"))
    if (!file.exists(protein_file)) {
      stop("Protein file not found. Please specify the correct path.")
    }
  } else if (!file.exists(file.path(workdir, protein_file))) {
    stop("Protein file not found. Please specify the correct path.")
  }
  
  message("Loading protein data...")
  protein_data <- read.csv(file.path(workdir, protein_file))
  
  # Load imaging data
  message("Loading imaging data...")
  imdata <- tryCatch({
    Cardinal::readImzML(datafile, folder = projectfolder)
  }, error = function(e) {
    stop(paste("Error loading imaging data:", e$message))
  })
  
  # If no specific proteins provided, use top scoring
  if (is.null(protein_list)) {
    if ("Proscore" %in% colnames(protein_data)) {
      # Group by protein and get max score
      protein_summary <- protein_data %>%
        group_by(Protein) %>%
        summarize(
          Proscore = max(Proscore, na.rm = TRUE),
          desc = first(desc)
        ) %>%
        arrange(desc(Proscore)) %>%
        head(max_proteins)
      
      protein_list <- protein_summary$Protein
      protein_desc <- protein_summary$desc
    } else {
      # If no Proscore column, use unique proteins
      protein_list <- unique(protein_data$Protein)[1:min(max_proteins, length(unique(protein_data$Protein)))]
    }
  }
  
  message(paste("Processing", length(protein_list), "proteins..."))
  
  # Create an empty list to store protein pixel data
  protein_pixel_data <- list()
  protein_images <- list()
  
  # Process each protein
  for (i in seq_along(protein_list)) {
    protein_id <- protein_list[i]
    message(paste("Processing protein", protein_id))
    
    # Get peptides for this protein
    protein_peptides <- protein_data[protein_data$Protein == protein_id, ]
    
    if (nrow(protein_peptides) == 0) {
      message(paste("No peptides found for protein", protein_id))
      next
    }
    
    peptide_mzs <- unique(protein_peptides$mz)
    
    # Create protein score image by combining peptide intensities
    protein_image <- NULL
    
    for (j in seq_along(peptide_mzs)) {
      mz_value <- peptide_mzs[j]
      
      # Get peptide intensity at this m/z value for each pixel
      peptide_image <- try(
        image(imdata, mz = mz_value,
              plusminus = ppm * mz_value / 1000000,
              contrast.enhance = "none", # No enhancement for raw data
              smooth.image = "none",     # No smoothing for raw data
              normalize.image = "none",  # No normalization for raw data
              key = FALSE,
              plot = FALSE)
      )
      
      if (inherits(peptide_image, "try-error")) {
        message(paste("Error processing m/z", mz_value))
        next
      }
      
      # For first peptide, initialize protein image
      if (is.null(protein_image)) {
        protein_image <- peptide_image
      } else {
        # Add intensities (weighted by peptide score if available)
        if ("Score" %in% colnames(protein_peptides)) {
          # Get score for this peptide
          peptide_score <- mean(protein_peptides$Score[protein_peptides$mz == mz_value])
          protein_image <- protein_image + (peptide_image * peptide_score)
        } else {
          protein_image <- protein_image + peptide_image
        }
      }
    }
    
    if (is.null(protein_image)) {
      message(paste("Failed to create image for protein", protein_id))
      next
    }
    
    # Normalize protein image if requested
    if (pixelwise_normalization == "sum") {
      # Normalize to total ion current
      total_intensity <- image(imdata, 
                               contrast.enhance = "none",
                               smooth.image = "none",
                               normalize.image = "none",
                               key = FALSE,
                               plot = FALSE)
      
      protein_image <- protein_image / total_intensity
    } else if (pixelwise_normalization == "max") {
      # Normalize to maximum intensity
      protein_image <- protein_image / max(protein_image, na.rm = TRUE)
    }
    
    # Store raw pixel data
    protein_data_df <- data.frame(
      protein_id = protein_id,
      x = coordPixels(imdata)$x,
      y = coordPixels(imdata)$y,
      intensity = as.vector(protein_image)
    )
    
    protein_pixel_data[[as.character(protein_id)]] <- protein_data_df
    
    # Visualize protein distribution
    protein_name <- ifelse("desc" %in% colnames(protein_peptides), 
                         protein_peptides$desc[1], 
                         paste("Protein", protein_id))
    
    # Create an output file name
    output_name <- paste0("protein_", protein_id, "_", gsub("[^a-zA-Z0-9]", "_", substr(protein_name, 1, 30)))
    
    # Generate and save visualization
    png(file.path(protein_output_dir, paste0(output_name, ".png")), 
        width = 1200, height = 1000, res = 150)
    
    # Display visualization
    image(imdata, protein_image,
          contrast.enhance = contrast.enhance,
          smooth.image = smooth.image,
          normalize.image = "linear",
          colorscale = intensity.colors_customize1(colset = ifelse(Component_plot_coloure == "mono", 2, 1)),
          main = paste("Protein", protein_id, "-", substr(protein_name, 1, 40)),
          layout = NULL)
    
    dev.off()
    
    # Store image references
    protein_images[[as.character(protein_id)]] <- file.path(protein_output_dir, paste0(output_name, ".png"))
    
    # Optionally create grid images with underlying peptides
    if (nrow(protein_peptides) >= 2) {
      # Create a temporary file
      temp_smplist <- protein_peptides
      
      # If needed, rename columns to match cluster_image_grid expected format
      if (!"Peptide" %in% colnames(temp_smplist) && "moleculeNames" %in% colnames(temp_smplist)) {
        temp_smplist$Peptide <- temp_smplist$moleculeNames
      }
      
      # Add ClusterID column if not present
      if (!"Protein" %in% colnames(temp_smplist)) {
        temp_smplist$Protein <- protein_id
      }
      
      # Try to create cluster image grid
      tryCatch({
        cluster_image_grid(
          clusterID = protein_id,
          SMPLIST = temp_smplist,
          imdata = imdata,
          ClusterID_colname = "Protein",
          componentID_colname = "Peptide",
          Component_plot_threshold = 2,
          Component_plot_coloure = Component_plot_coloure,
          smooth.image = smooth.image,
          contrast.enhance = contrast.enhance,
          plot_style = plot_style,
          workdir = protein_output_dir,
          ppm = ppm
        )
      }, error = function(e) {
        message(paste("Error creating cluster image for protein", protein_id, ":", e$message))
      })
    }
  }
  
  # Export combined data if requested
  if (export_data) {
    combined_data <- do.call(rbind, protein_pixel_data)
    write.csv(combined_data, file.path(protein_output_dir, "protein_pixel_data.csv"), row.names = FALSE)
  }
  
  # Create a summary visualization of all proteins
  if (length(protein_pixel_data) > 1 && length(protein_pixel_data) <= 20) {
    message("Creating combined protein visualization...")
    
    # Combine all protein images into one grid
    all_images <- lapply(protein_images, image_read)
    
    # Determine grid dimensions
    n_images <- length(all_images)
    grid_width <- min(5, n_images)
    grid_height <- ceiling(n_images / grid_width)
    
    # Resize all images to the same dimensions
    all_images_resized <- lapply(all_images, function(img) {
      image_resize(img, "800x800")
    })
    
    # Create image grid
    combined_image <- image_append(
      image_join(all_images_resized),
      stack = TRUE
    )
    
    # Save combined image
    image_write(combined_image, file.path(protein_output_dir, "all_proteins.png"))
  }
  
  message("Protein score visualization complete.")
  
  # Return invisible list of results
  invisible(list(
    protein_pixel_data = protein_pixel_data,
    protein_images = protein_images,
    output_dir = protein_output_dir
  ))
}

#' Calculate Protein Scores for Each Pixel
#'
#' This function calculates protein scores for each pixel in an imaging dataset
#' 
#' @param imdata Cardinal imaging data object
#' @param protein_data Protein data frame with mappings to peptides
#' @param protein_id ID of the protein to calculate scores for
#' @param ppm Mass tolerance in parts per million
#' @param normalization Normalization method to use
#'
#' @return A matrix of protein scores for each pixel
#'
#' @examples
#' # This function is used internally by export_protein_pixel_data
#'
#' @keywords internal
calculate_protein_pixel_scores <- function(imdata, protein_data, protein_id, ppm = 5, normalization = "sum") {
  # Get peptides for this protein
  protein_peptides <- protein_data[protein_data$Protein == protein_id, ]
  
  if (nrow(protein_peptides) == 0) {
    warning(paste("No peptides found for protein", protein_id))
    return(NULL)
  }
  
  peptide_mzs <- unique(protein_peptides$mz)
  
  # Create protein score image by combining peptide intensities
  protein_image <- NULL
  
  for (j in seq_along(peptide_mzs)) {
    mz_value <- peptide_mzs[j]
    
    # Get peptide intensity at this m/z value for each pixel
    peptide_image <- try(
      image(imdata, mz = mz_value,
            plusminus = ppm * mz_value / 1000000,
            contrast.enhance = "none",
            smooth.image = "none",
            normalize.image = "none",
            key = FALSE,
            plot = FALSE)
    )
    
    if (inherits(peptide_image, "try-error")) {
      warning(paste("Error processing m/z", mz_value))
      next
    }
    
    # For first peptide, initialize protein image
    if (is.null(protein_image)) {
      protein_image <- peptide_image
    } else {
      # Add intensities (weighted by peptide score if available)
      if ("Score" %in% colnames(protein_peptides)) {
        # Get score for this peptide
        peptide_score <- mean(protein_peptides$Score[protein_peptides$mz == mz_value])
        protein_image <- protein_image + (peptide_image * peptide_score)
      } else {
        protein_image <- protein_image + peptide_image
      }
    }
  }
  
  if (is.null(protein_image)) {
    warning(paste("Failed to create image for protein", protein_id))
    return(NULL)
  }
  
  # Normalize protein image if requested
  if (normalization == "sum") {
    # Normalize to total ion current
    total_intensity <- image(imdata, 
                             contrast.enhance = "none",
                             smooth.image = "none",
                             normalize.image = "none",
                             key = FALSE,
                             plot = FALSE)
    
    protein_image <- protein_image / total_intensity
  } else if (normalization == "max") {
    # Normalize to maximum intensity
    protein_image <- protein_image / max(protein_image, na.rm = TRUE)
  }
  
  return(protein_image)
}

#' Create PCA Analysis of Protein Distributions
#'
#' This function performs PCA on protein distributions and visualizes the results
#' 
#' @param protein_pixel_data List of data frames containing protein pixel data
#' @param output_dir Directory to save output files
#' @param n_components Number of principal components to calculate
#'
#' @return List containing PCA results and visualizations
#'
#' @examples
#' # This function is used with results from export_protein_pixel_data
#'
#' @export
visualize_protein_pca <- function(protein_pixel_data, output_dir = NULL, n_components = 3) {
  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(reshape2)))
  
  if (is.null(output_dir)) {
    output_dir <- getwd()
  }
  
  # Ensure output directory exists
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Combine protein data into a wide format
  combined_data <- protein_pixel_data[[1]][, c("x", "y")]
  
  for (protein_id in names(protein_pixel_data)) {
    protein_df <- protein_pixel_data[[protein_id]]
    combined_data[[paste0("protein_", protein_id)]] <- protein_df$intensity
  }
  
  # Remove rows with any NA values
  combined_data <- na.omit(combined_data)
  
  if (nrow(combined_data) < 10) {
    warning("Not enough data points for PCA analysis after removing NAs")
    return(NULL)
  }
  
  # Extract coordinates and protein intensities
  coords <- combined_data[, c("x", "y")]
  protein_intensities <- combined_data[, grepl("protein_", colnames(combined_data))]
  
  # Perform PCA
  pca_result <- prcomp(protein_intensities, scale. = TRUE)
  
  # Get PCA scores
  pca_scores <- data.frame(pca_result$x)
  pca_scores <- cbind(coords, pca_scores)
  
  # Calculate variance explained
  var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)
  
  # Create visualizations for each principal component
  pca_images <- list()
  
  for (i in 1:min(n_components, ncol(pca_scores) - 2)) {
    pc_col <- paste0("PC", i)
    
    # Create visualization
    p <- ggplot(pca_scores, aes(x = x, y = y, fill = .data[[pc_col]])) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
      labs(
        title = paste("Principal Component", i),
        subtitle = paste0("Explains ", round(var_explained[i] * 100, 1), "% of variance"),
        x = "X coordinate",
        y = "Y coordinate"
      ) +
      theme_minimal() +
      coord_fixed()
    
    # Save plot
    output_file <- file.path(output_dir, paste0("protein_pca_component_", i, ".png"))
    ggsave(output_file, p, width = 8, height = 6, dpi = 150)
    
    pca_images[[pc_col]] <- output_file
  }
  
  # Create loadings plot
  loadings <- as.data.frame(pca_result$rotation)
  loadings$protein <- gsub("protein_", "", rownames(loadings))
  
  loadings_melted <- melt(loadings, id.vars = "protein", 
                         variable.name = "Component", 
                         value.name = "Loading")
  
  # Filter to include only the top components
  loadings_melted <- loadings_melted[loadings_melted$Component %in% paste0("PC", 1:min(n_components, ncol(loadings) - 1)), ]
  
  # Create loadings plot
  p_loadings <- ggplot(loadings_melted, aes(x = protein, y = Loading, fill = Component)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(
      title = "PCA Loadings by Protein",
      x = "Protein",
      y = "Loading"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save loadings plot
  output_file <- file.path(output_dir, "protein_pca_loadings.png")
  ggsave(output_file, p_loadings, width = 10, height = 6, dpi = 150)
  
  # Return results
  return(list(
    pca_result = pca_result,
    pca_scores = pca_scores,
    var_explained = var_explained,
    pca_images = pca_images,
    loadings_plot = output_file
  ))
}
