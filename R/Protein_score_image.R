


#' Generate_proscore_image
#'
#' This is a sccore based image function for Proteomics imaging data set
#' this function will read the Protein peptide summary file and generate protein score image by integrating the score map of the peptides.
#' The function will read the candidatelist file and pre-processed (.rds) file from the \code{IMS_analysis} function.
#' @param datafile specify the imzML data files to be used
#' @param threshold specify the intensities threshold (0 to 1 in percentage)to report a identified molecule
#' @param ppm the mz tolerance (in ppm) for peak integration
#' @param Protein_feature  \code{"IMS_analysis"} result file that contains the identified proteins/peptides information for scoring image
#'
#' @examples
#'
#' @export
#'

Generate_proscore_image <- function(projectfolder = NULL,
  datafile,
  threshold = 0.05,
  ppm = 5,
  Protein_feature = FALSE,
  plot_cluster_image = FALSE,
  plot_ion_image = FALSE,
  plot_score_image = FALSE,
  Load_candidatelist = TRUE,
  score_method = "SQRTP",
  plot_cluster_image_grid = TRUE,
  plot_ion_image_grid = TRUE,
  plot_cluster_image_maxretry = 2,
  plot_cluster_image_overwrite = F,
  smooth.image = "gaussian",
  componentID_colname = "Peptide",
  ClusterID_colname = "Protein",
  Protein_desc_of_interest = ".",
  Protein_desc_of_exclusion = NULL,
  plot_unique_component = TRUE,
  FDR_cutoff = 0.05,
  use_top_rank = NULL,
  plot_matching_score = F,
  Component_plot_coloure = "mono",
  cluster_color_scale = "blackwhite",
  plot_layout = "line",
  export_Header_table = T,
  export_footer_table = T,
  attach_summary_cluster = T,
  remove_cluster_from_grid = attach_summary_cluster,
  pixel_size_um = 50,
  img_brightness = 100,
  Thread = 4,
  cluster_rds_path = NULL,
  remove_score_outlier = F,
  Plot_score_IQR_cutoff = 0.75,
  Plot_score_abs_cutoff = -0.1,
  mzAlign_runs = "TopNfeature_mean") {
  # Load the required packages
  
  suppressMessages(suppressWarnings(library("pacman")))
  suppressMessages(suppressWarnings(
    p_load(stringr, BiocParallel, data.table, Cardinal, parallel)
  ))
  
    message("cluster image rendering...")
    workdir=paste0(projectfolder,"/")
    setwd(workdir[1])
    
    # read protein-peptide features result
    if (is.null(Protein_feature)) {
      Protein_feature_list = read.csv(
        file = paste(
          workdir[1],
          "/Summary folder/Protein_peptide_Summary.csv",
          sep = ""
        ),
        stringsAsFactors = F
      )
    } else{
      Protein_feature_list = Protein_feature
    }
    # remove peptide score outlier from result
    
    if (remove_score_outlier) {
      Protein_feature_list <- remove_pep_score_outlier(Protein_feature_list,
        abs_cutoff = Plot_score_abs_cutoff,
        IQR_LB = Plot_score_IQR_cutoff)
    }
    
    # extract the protein entries of interest
    
    if (sum(Protein_desc_of_interest != ".") >= 1) {
      Protein_feature_list_interest <- NULL
      num_of_interest <- numeric(0)
      
      for (interest_desc in Protein_desc_of_interest) {
        idx_iterest_desc <- str_detect(Protein_feature_list$desc,
          regex(interest_desc, ignore_case = T))
        if (nrow(Protein_feature_list[idx_iterest_desc, ]) != 0) {
          Protein_feature_list_interest <- rbind(Protein_feature_list_interest,
            Protein_feature_list[idx_iterest_desc, ])
        }
        num_of_interest[interest_desc] <- length(unique(Protein_feature_list[idx_iterest_desc, "Protein"]))
      }
      
      Protein_feature_list = Protein_feature_list_interest
      message(
        paste(
          num_of_interest,
          "Protein(s) found with annotations of interest:",
          Protein_desc_of_interest,
          collapse = "\n"
        )
      )
    }
    
    Protein_feature_list = as.data.frame(Protein_feature_list)
    
    if (!is.data.frame(Protein_feature_list)){
      stop("Protein feature input is not a data frame")
    }
    
    if (sum(c(componentID_colname, ClusterID_colname) %in% colnames(Protein_feature_list)) != 2) {
      stop("Protein feature input does not contain the required columns",
        c(componentID_colname, ClusterID_colname)[!(c(componentID_colname, ClusterID_colname) %in% colnames(Protein_feature_list))])
    }
    
    # generate combined IMS data for multiple files or use a link to load the pre-processed IMS data
    
    if (!is.null(cluster_rds_path)) {
      if(is.character(cluster_rds_path)){
        if (file.exists(cluster_rds_path)) {
          imdata = readRDS(cluster_rds_path)
          message("cluster imdata loaded.")
        } else{
          stop("cluster_rds_path does not exist")
        }
      } else if (!typeof(cluster_rds_path) %in% c("MSImagingExperiment","MSImagingArrays")){
        stop("cluster_rds_path is not a character")
      }
    } else{
      message("Processed image data is NULL will generate the data via",mzAlign_runs, "and target","mode.")
      cluster_rds_path <- Load_IMS_decov_combine(
        datafile = datafile,
        workdir = workdir,
        import_ppm = ppm,
        SPECTRUM_batch = "overall",
        mzAlign_runs = mzAlign_runs,
        ppm = ppm,
        threshold = 0,
        rotate = NULL,
        mzrange = "auto-detect",
        deconv_peaklist = "Target",
        preprocessRDS_rotated = T,
        target_mzlist = sort(unique(
          as.numeric(Protein_feature_list$mz)
        ), decreasing = F)
      )
      
      
      imdata = readRDS(paste0(workdir[1], "/", basename(cluster_rds_path)))
      
      message("cluster imdata generated and loaded.")
    }
    
    # test combined imdata
    if (class(imdata)[1] == "matrix") {
      do.call(Cardinal::cbind, imdata) -> imdata
      
      saveRDS(imdata,
        paste0(workdir[1], "/combinedimdata.rds"),
        compress = T)
      
    }
    
    # Setup output folder and queue the R calls for cluster image randering
    outputfolder = paste(workdir, "/Summary folder/cluster Ion images/", sep =
        "")
    if (dir.exists(outputfolder) == FALSE) {
      dir.create(outputfolder)
    }
    
    if (!(plot_unique_component)) {
      setwd(outputfolder)
      Protein_feature_list_trimmed <- Protein_feature_list
    }
    
    
    if (plot_unique_component) {
      outputfolder = paste(workdir,
        "/Summary folder/cluster Ion images/unique/",
        sep = "")
      
      if (dir.exists(outputfolder) == FALSE) {
        dir.create(outputfolder)
      }
      setwd(outputfolder)
      Protein_feature_list_unique = Protein_feature_list %>% group_by(mz) %>% dplyr::summarise(num =
          length(unique(Protein)))
      Protein_feature_list_unique_mz <- Protein_feature_list_unique$mz[Protein_feature_list_unique$num ==
          1]
      Protein_feature_list_trimmed <- Protein_feature_list[Protein_feature_list$mz %in% Protein_feature_list_unique_mz, ]
      write.csv(
        Protein_feature_list_trimmed,
        paste(
          workdir,
          "/Summary folder/Protein_feature_list_trimmed.csv",
          sep = ""
        ),
        row.names = F
      )
    }
    
    #list_of_protein_sequence[!(1:length(list_of_protein_sequence) %in% as.numeric(Protein_feature_list_trimmed$Protein))]<-""
    #unname(as.character(list_of_protein_sequence))->str_vec
    #imdata <- imdata %>% peakBin(ref=sort(unique((Protein_feature_list_trimmed$mz_align+Protein_feature_list_trimmed$mz)/2)), tolerance=ppm, units="ppm") %>% process(BPPARAM=SerialParam())
    save(
      list = c(
        "Protein_feature_list_trimmed",
        "imdata",
        "ClusterID_colname",
        "componentID_colname",
        "plot_layout",
        "export_Header_table",
        "export_footer_table",
        "attach_summary_cluster",
        "remove_cluster_from_grid",
        "smooth.image",
        "Component_plot_coloure",
        "cluster_color_scale",
        "list_of_protein_sequence",
        "outputfolder",
        "peptide_ID_filter",
        "ppm",
        "img_brightness",
        "pixel_size_um"
      ),
      file = paste0(workdir, "/cluster_img_grid.RData")
    )
    if (plot_cluster_image_grid){
      for (clusterID in unique(Protein_feature_list_trimmed$Protein)) {
      cluster_desc <- unique(Protein_feature_list_trimmed$desc[Protein_feature_list_trimmed[[ClusterID_colname]] ==
          clusterID])
      cluster_desc <- gsub(stringr::str_extract(cluster_desc, "OS=.{1,}"),
        "",
        cluster_desc)
      n_component <- nrow(unique(Protein_feature_list_trimmed[Protein_feature_list_trimmed[[ClusterID_colname]] ==
          clusterID, c(
            ClusterID_colname,
            componentID_colname,
            "moleculeNames",
            "adduct",
            "Modification"
          )]))
      if (n_component >= peptide_ID_filter) {
        if ('&'(file.exists(
          paste0(outputfolder, clusterID, "_cluster_imaging.png")
        ), !plot_cluster_image_overwrite)) {
          message(
            "Cluster image rendering Skipped file exists: No.",
            clusterID,
            " ",
            cluster_desc
          )
          next
        } else {
          fileConn <- file(paste0(workdir, "/cluster_img_scource.R"), )
          writeLines(
            c(
              "suppressMessages(suppressWarnings(require(HiTMaP)))",
              paste0("clusterID=", clusterID),
              paste0(
                "suppressMessages(suppressWarnings(load(file =\"",
                workdir,
                "/cluster_img_grid.RData\")))"
              ),
              "suppressMessages(suppressWarnings(cluster_image_grid(clusterID = clusterID,
                                      imdata=imdata,
                                      SMPLIST=Protein_feature_list_trimmed,
                                      ppm=ppm,
                                      ClusterID_colname=ClusterID_colname,
                                      componentID_colname=componentID_colname,
                                      plot_layout=plot_layout,
                                      export_Header_table=export_Header_table,
                                      export_footer_table=export_footer_table,
                                      attach_summary_cluster=attach_summary_cluster,
                                      remove_cluster_from_grid=remove_cluster_from_grid,
                                      plot_style=\"fleximaging\",
                                      smooth.image=smooth.image,
                                      Component_plot_coloure=Component_plot_coloure,
                                      cluster_color_scale=cluster_color_scale,
                                      list_of_protein_sequence=list_of_protein_sequence,
                                      workdir=outputfolder,
                                      pixel_size_um=pixel_size_um,
                                      img_brightness=img_brightness,
                                      Component_plot_threshold=peptide_ID_filter)))"
            ),
            fileConn
          )
          close(fileConn)
          
          system(paste0(
            "Rscript \"",
            paste0(workdir, "/cluster_img_scource.R\"")
          ))
          
          
          if (file.exists(paste0(outputfolder, clusterID, "_cluster_imaging.png"))) {
            message("Cluster image rendering Done: No.",
              clusterID,
              " ",
              cluster_desc)
            
          } else{
            retrytime = 1
            repeat {
              message(
                "Cluster image rendering failed and retry ",
                retrytime,
                ": No.",
                clusterID,
                " ",
                cluster_desc
              )
              system(paste0(
                "Rscript \"",
                paste0(workdir, "/cluster_img_scource.R\"")
              ))
              
              if (file.exists(paste0(
                outputfolder,
                clusterID,
                "_cluster_imaging.png"
              ))) {
                message("Cluster image rendering Done: No.",
                  clusterID,
                  " ",
                  cluster_desc)
                
                break
              } else if (retrytime >= plot_cluster_image_maxretry) {
                message(
                  "Cluster image rendering reaches maximum Retry Attempts: No.",
                  clusterID,
                  " ",
                  cluster_desc
                )
                
                break
              }
              retrytime = 1 + retrytime
            }
            
          }
          
        }
        
        
      }
    }
    }
    if (plot_score_image_grid){
      for (clusterID in unique(Protein_feature_list_trimmed$Protein)) {
      cluster_desc <- unique(Protein_feature_list_trimmed$desc[Protein_feature_list_trimmed[[ClusterID_colname]] ==
          clusterID])
      cluster_desc <- gsub(stringr::str_extract(cluster_desc, "OS=.{1,}"),
        "",
        cluster_desc)
      n_component <- nrow(unique(Protein_feature_list_trimmed[Protein_feature_list_trimmed[[ClusterID_colname]] ==
          clusterID, c(
            ClusterID_colname,
            componentID_colname,
            "moleculeNames",
            "adduct",
            "Modification"
          )]))
      if (n_component >= peptide_ID_filter) {
        if ('&'(file.exists(
          paste0(outputfolder, clusterID, "_score_imaging.png")
        ), !plot_cluster_image_overwrite)) {
          message(
            "Cluster image rendering Skipped file exists: No.",
            clusterID,
            " ",
            cluster_desc
          )
          next
        } else {
          fileConn <- file(paste0(workdir, "/cluster_img_scource.R"), )
          writeLines(
            c(
              "suppressMessages(suppressWarnings(require(HiTMaP)))",
              paste0("clusterID=", clusterID),
              paste0(
                "suppressMessages(suppressWarnings(load(file =\"",
                workdir,
                "/cluster_img_grid.RData\")))"
              ),
              "suppressMessages(suppressWarnings(cluster_score_grid(clusterID = clusterID,
                                      imdata=imdata,
                                      SMPLIST=Protein_feature_list_trimmed,
                                      ppm=ppm,
                                      ClusterID_colname=ClusterID_colname,
                                      componentID_colname=componentID_colname,
                                      plot_layout=plot_layout,
                                      export_Header_table=export_Header_table,
                                      export_footer_table=export_footer_table,
                                      attach_summary_cluster=attach_summary_cluster,
                                      remove_cluster_from_grid=remove_cluster_from_grid,
                                      plot_style=\"fleximaging\",
                                      smooth.image=smooth.image,
                                      Component_plot_coloure=Component_plot_coloure,
                                      cluster_color_scale=cluster_color_scale,
                                      list_of_protein_sequence=list_of_protein_sequence,
                                      workdir=outputfolder,
                                      pixel_size_um=pixel_size_um,
                                      img_brightness=img_brightness,
                                      Component_plot_threshold=peptide_ID_filter)))"
            ),
            fileConn
          )
          close(fileConn)
          
          system(paste0(
            "Rscript \"",
            paste0(workdir, "/cluster_img_scource.R\"")
          ))
          
          
          if (file.exists(paste0(outputfolder, clusterID, "_score_imaging.png"))) {
            message("Cluster image rendering Done: No.",
              clusterID,
              " ",
              cluster_desc)
            
          } else{
            retrytime = 1
            repeat {
              message(
                "Cluster image rendering failed and retry ",
                retrytime,
                ": No.",
                clusterID,
                " ",
                cluster_desc
              )
              system(paste0(
                "Rscript \"",
                paste0(workdir, "/cluster_img_scource.R\"")
              ))
              
              if (file.exists(paste0(
                outputfolder,
                clusterID,
                "_score_imaging.png"
              ))) {
                message("Cluster image rendering Done: No.",
                  clusterID,
                  " ",
                  cluster_desc)
                
                break
              } else if (retrytime >= plot_cluster_image_maxretry) {
                message(
                  "Cluster image rendering reaches maximum Retry Attempts: No.",
                  clusterID,
                  " ",
                  cluster_desc
                )
                
                break
              }
              retrytime = 1 + retrytime
            }
            
          }
          
        }
        
        
      }
    }
    }
    
    
  
}