#' imaging_identification
#'
#' This is a peptide mass fingerprint search function for maldi imaging data analysis
#' @param datafile the data files' path for the analysis, leave it as blank to enable a graphical user interface to select the data
#' @param projectfolder optional, if NULL script will extract the path from datafile(s), and use the first workdir as project folder
#' @param parallel the number of threads will be used in the PMF search, this option now only works for windows OS
#' @param threshold specify the intensities threshold (0 to 1 in percentage)to report a identified molecule
#' @param ppm the mz tolerance (in ppm) for peak integration
#' @param Digestion_site Set the enzyme digestion specificity by one or more regex expressions or the name of a enzyme
#' @param missedCleavages miss cleavage number allowed in this PMF search
#' @param Fastadatabase the fasta database used in this pmf search, the file should be placed in the same folder with data files
#' @param adducts the adducts list to be used for generating the PMF search candidates
#' @param Modifications set the modifications
#' @param Substitute_AA set the amino acid Substitutions
#' @param Decoy_search enable (default) or disable the decoy search
#' @param Decoy_mode select the decoy search mode between "isotope" (default), "element" and "adduct"
#' @param Decoy_adducts define the adduct list for decoy search. the decoy adducts could be "M+ACN+H","M+IsoProp+H","M+DMSO+H","M+Co","M+Ag","M+Cu","M+He","M+Ne","M+Ar","M+Kr","M+Xe" or"M+Rn".
#' @param mzrange define the mz range for the experiment, default is 700 to 4000 m/z.
#' @param use_previous_candidates set as TRUE to reload the previously generated candidate list.
#' @param IMS_analysis Set \code{"true"} if you want to perform data pre-processing and proteomics search, set \code{"false"} if you want to bypass it
#' @param FDR_cutoff set the protein FDR cutoff threshold, default is 5 percent
#' @param score_method specify the peptide spectrum scoring method, "SQRTP" is recommended.
#' @param peptide_ID_filter set the minimal count of peptides needed to identify a protein
#' @param plot_matching_score enable the spectrum matching overlay plot
#' @param Protein_feature_summary \code{"IMS_analysis"} follow-up process that will collect all the identified peptide information and associate them with possible proteins
#' @param Peptide_feature_summary \code{"IMS_analysis"} follow-up process that will summarize all datafiles identified peptides and generats a \code{"peptide shortlist"} in the result summary folder
#' @param Region_feature_summary \code{"IMS_analysis"} follow-up process that will summarize mz feature of all regions of all data files into the summary folder
#' @param plot_ion_image \code{"Peptide_feature_summarya"} follow-up process that will plot every connponents in the \code{"peptide shortlist"}. please use the cluster image grid to output the images.
#' @param preprocess a list of params that define the IMS data pre-processing procedure
#' @param spectra_segments_per_file optimal number of distinctive regions in the tissue section, a virtual segmentation will be applied to the image files with this value. To have a better PMF result you may set a value that in the sweet point of sensitivety and false discovery rate (FDR).
#' @param Segmentation set as "spatialKMeans" to enable a \code{"spatialKMeans"} Segmentation; set as "spatialShrunkenCentroids" to enable a \code{"spatialShrunkenCentroids"} Segmentation; If a region rank file was supplied, you can set this as "Virtual_segmentation" to perform a manual segmentation; Set it as "none" to bypass the segmentation.
#' @param Smooth_range \code{"Segmentation"} pixel smooth range
#' @param Virtual_segmentation_rankfile specify a region rank file contains region information for manualy region segmentation
#' @param Rotate_IMG specify a configuration file to further change the rotation of the images
#' @param plot_cluster_image_grid set as \code{"TRUE"} to enable the protein cluster image function.
#' @param plot_layout Set as \code{"line"} to plot cluster and component images for multiple data file or as \code{"grid"} to plot cluster images for single data file. In "grid" mode, Image's will be rendered into a grid with 5 columns.
#' @param ClusterID_colname Specify the cluster ID column in the result spreadsheet.
#' @param componentID_colname Specify the component ID column in the result spreadsheet.
#' @param Protein_desc_of_interest Specify a list of protein descriptions for cluster image plotting. Default setting will plot all reported proteins.
#' @param Protein_desc_of_exclusion Specify a list of protein descriptions to be excluded from cluster image plotting.
#' @param plot_unique_component Set as \code{"TRUE"} to plot only the unique components in the cluster image plotting.
#' @param Component_plot_coloure set as "mono" to use a pre-defined color scale to plot component images. Set as "as.cluster" to use the previously assigned mono color in the additive cluster binning process.
#' @param cluster_color_scale Set as "blackwhite" to use only black and white color in the cluster image plotting. using "blackwhite" in cluster_color_scale will overwrite the components' color setting.
#' @param export_Header_table Set as \code{"TRUE"} to plot the header in the cluster image plotting. Header table includes the basic information of cluster and components.
#' @param export_footer_table Set as \code{"TRUE"} to plot the footer in the cluster image plotting. Footer shows the protein coverage in the Proteomics mode.
#' @param attach_summary_cluster Set as \code{"TRUE"} to attach an enlarged cluster image to the bottom of the cluster image.
#' @param remove_cluster_from_grid Set as \code{"TRUE"} to remove the cluster image from the cluster image grid. it is recommended to set this same as the attach_summary_cluster.
#' @param plot_cluster_image_overwrite Set as true to generate the cluster images regardless the existance of previously file(s)
#' @param cluster_rds_path set as NULL if there is not preprocessed.rds available for a single file, script will load the raw data file which may reduce the signal intensities. For multiple samples, scripts will try to load the RDS file from each "ID" folder and merge the mz features via instrument resolution setting and output a combined RDS file to the project folder. For multiple files cluster images rendering user should set the attach_summary_cluster as False, and set remove_cluster_from_grid as true.
#' 
#' 
#' @return None
#'
#' @examples
#' imaging_identification(threshold=0.05, ppm=5,Digestion_site="[G]",
#'                        missedCleavages=0:1,Fastadatabase="murine_matrisome.fasta",
#'                        adducts=c("M+H","M+NH4","M+Na"),IMS_analysis=TRUE,
#'                        Protein_feature_summary=TRUE,plot_cluster_image=TRUE,
#'                        Peptide_feature_summary=TRUE,plot_ion_image=FALSE,
#'                        parallel=3,spectra_segments_per_file=5,Segmentation="spatialKMeans"
#'                        )
#'
#' @export
#'
imaging_identification<-function(
#==============Choose the imzml raw data file(s) to process  make sure the fasta file in the same folder
               datafile,
               projectfolder=NULL,
               threshold=0.005,
               ppm=5,
               mode=c("Proteomics","Metabolomics"),
               Digestion_site="trypsin",
               missedCleavages=0:1,
               Fastadatabase="uniprot-bovin.fasta",
               adducts=c("M+H"),
               Modifications=list(fixed=NULL,fixmod_position=NULL,variable=NULL,varmod_position=NULL),
               Substitute_AA=NULL,
               Decoy_search=TRUE,
               Decoy_adducts=c("M+ACN+H","M+IsoProp+H","M+DMSO+H","M+Co","M+Ag","M+Cu","M+He","M+Ne","M+Ar","M+Kr","M+Xe","M+Rn"),
               Decoy_mode = "isotope",
               mzrange=c(700,4000),
               Database_stats=F,
               adjust_score = FALSE,
               IMS_analysis=TRUE,
               PMFsearch=IMS_analysis,
               Load_candidatelist=IMS_analysis || plot_cluster_image_grid,
               Bypass_generate_spectrum=FALSE,
               peptide_ID_filter=2,
               Protein_feature_summary=TRUE,
               Peptide_feature_summary=TRUE,
               plot_ion_image=FALSE,
               parallel=detectCores(),
               spectra_segments_per_file=4,
               Segmentation=c("spatialKMeans","spatialShrunkenCentroids","Virtual_segmentation","none","def_file"),
               Segmentation_def="Segmentation_def.csv",
               Segmentation_ncomp="auto-detect",
               Segmentation_variance_coverage=0.8,
               preprocess=list(force_preprocess=FALSE,use_preprocessRDS=TRUE,smoothSignal=list(method="disable"),
                               reduceBaseline=list(method="locmin"),
                               peakPick=list(method="adaptive"),
                               peakAlign=list(tolerance=ppm/2, units="ppm"),
                               normalize=list(method=c("rms","tic","reference")[1],mz=1)),
               Smooth_range=1,
               Virtual_segmentation_rankfile=NULL,
               Rotate_IMG=NULL,
               Region_feature_summary=FALSE,
               Spectrum_validate=TRUE,
               output_candidatelist=TRUE,
               use_previous_candidates=FALSE,
               score_method="SQRTP",
               plot_cluster_image_grid=FALSE,
               plot_cluster_image_maxretry=2,
               plot_cluster_image_overwrite=F,
               smooth.image="gaussian",
               componentID_colname="Peptide",
               ClusterID_colname="Protein",
               Protein_desc_of_interest=".",
               Protein_desc_of_exclusion=NULL,
               plot_unique_component=TRUE,
               FDR_cutoff=0.05,
               use_top_rank=NULL,
               plot_matching_score=F,
               Component_plot_coloure="mono",
               cluster_color_scale="blackwhite",
               plot_layout="line",
               export_Header_table=T,
               export_footer_table=T,
               attach_summary_cluster=T,
               remove_cluster_from_grid=attach_summary_cluster,
               pixel_size_um=50,
               img_brightness=100,
               Thread=4,
               cluster_rds_path=NULL,
               remove_score_outlier=F,
               Plot_score_IQR_cutoff=0.75,
               Plot_score_abs_cutoff=-0.1,
               ...
               ){
  suppressMessages(suppressWarnings(library("pacman")))
  suppressMessages(suppressWarnings(p_load(stringr,BiocParallel,data.table,Cardinal)))

  if (missing(datafile)) stop("Missing data file, Choose single or multiple imzml file(s) for analysis")
  
# resove missing vars this is an patch during package evaluation and will be removed in next revision
  e <- environment()
  p <- parent.env(e)

  if(!exists("ppm", envir=p)) {
    pf <- parent.frame()
    pf$ppm <- ppm / 2
    message("Alignment tolerance missing, using half of identification ppm (half FWHM) to fit peak aligment algorithm")}
  
  
# retrieve/parse the working dir info, and convert the filenames
  if (is.null(projectfolder)){
    workdir<-base::dirname(datafile[1])
  }else{ workdir<-projectfolder }

  datafile <- basename(datafile)
  datafile <- gsub(".imzML$", "", datafile)
  datafile_imzML <- paste0(datafile,".imzML")
  
  setwd(paste0(workdir[1],"/"))

# Set the parallel processing parameter, multicore-fork method has been temporarily disabled due to the reduced performance in docker enviornment
  if (is.null(Thread)){
  parallel=try(detectCores()/2)
  if (parallel<1 | is.null(parallel)){parallel=1}
  BPPARAM=HiTMaP:::Parallel.OS(parallel, override_type=BiocParallel::SnowParam())
  setCardinalBPPARAM(BPPARAM = BPPARAM)
  }else{
  parallel=Thread
  BPPARAM=Parallel.OS(parallel)
  setCardinalBPPARAM(BPPARAM = BPPARAM)
  }



  message(paste(try(detectCores()), "Cores detected,",parallel, "threads will be used for computing"))

  message(paste(length(datafile), "files were selected and will be used for Searching"))

  message(paste(Fastadatabase, "was selected as database.", "Candidates will be generated through",mode[1] ,"mode" ))

  
# prepare the proteome database, use_previous_candidates=T will override the other argument and load the candidate list directly

  if(Load_candidatelist){

  Protein_feature_list<-Protein_feature_list_fun(workdir=workdir,
                                                 database=Fastadatabase,
                                                 Digestion_site=Digestion_site,
                                                 missedCleavages=missedCleavages,
                                                 adducts=adducts,
                                                 BPPARAM = BPPARAM,
                                                 Decoy_adducts=Decoy_adducts,
                                                 Decoy_search=Decoy_search,
                                                 Decoy_mode = Decoy_mode,
                                                 output_candidatelist=output_candidatelist,
                                                 use_previous_candidates=use_previous_candidates,
                                                 Substitute_AA=Substitute_AA,
                                                 Modifications=Modifications,
                                                 mzrange=mzrange,
                                                 Protein_desc_of_exclusion=Protein_desc_of_exclusion,
                                                 Database_stats=Database_stats)

  }
  
  # 

  if(IMS_analysis){
    
  message(paste(Fastadatabase,"was selected as database","\nSpectrum intensity threshold:",percent(threshold),"\nmz tolerance:",ppm,"ppm","Segmentation method:",Segmentation[1],
                "\nManual segmentation def file:",ifelse(is.null(Virtual_segmentation_rankfile),"None",Virtual_segmentation_rankfile),"\nBypass spectrum generation:",Bypass_generate_spectrum))
  
  #select candidate list for IMS annotation 
    
  Peptide_Summary_searchlist<-unique(Protein_feature_list)

  Peptide_Summary_file<-IMS_data_process(datafile=datafile, workdir=workdir,
                                                  Peptide_Summary_searchlist=Peptide_Summary_searchlist,
                                                  segmentation_num=spectra_segments_per_file,
                                                  threshold=threshold,rotate = Rotate_IMG,
                                                  ppm=ppm,mzrange=mzrange,
                                                  Segmentation=Segmentation,
                                                  Segmentation_ncomp=Segmentation_ncomp,
                                                  PMFsearch = PMFsearch,
                                                  Virtual_segmentation_rankfile = Virtual_segmentation_rankfile,
                                                  BPPARAM = BPPARAM,
                                                  Bypass_generate_spectrum=Bypass_generate_spectrum,
                                                  score_method = score_method,
                                                  Decoy_mode=Decoy_mode,
                                                  Decoy_search=Decoy_search,
                                                  adjust_score=adjust_score,
                                                  peptide_ID_filter=peptide_ID_filter,
                                                  Protein_desc_of_interest=Protein_desc_of_interest,
                                                  plot_matching_score_t=plot_matching_score,
                                                  FDR_cutoff= FDR_cutoff,
                                                  Segmentation_def=Segmentation_def,
                                                  Segmentation_variance_coverage=Segmentation_variance_coverage,
                                                  preprocess=preprocess)

  }
  
  
  
  #Summarize the protein result across the datafiles and store these summarized files into the summary folder
  
  if(Protein_feature_summary){
    message("Protein feature summary...")
    Peptide_Summary_file<-NULL
    Protein_peptide_Summary_file<-NULL
    protein_feature_all<-NULL
  for (i in 1:length(datafile)){
  datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
  currentdir<-paste0(workdir,"/",datafile[i]," ID")
  setwd(paste(currentdir,sep=""))
  protein_feature<-NULL
  for (protein_feature_file in dir()[stringr::str_detect(dir(),"Protein_segment_PMF_RESULT_")]){
    protein_feature<-fread(protein_feature_file)

    if(nrow(protein_feature)!=0){
    protein_feature$Source<-datafilename
    region_code<-str_replace(protein_feature_file,"Protein_segment_PMF_RESULT_","")
    region_code<-str_replace(region_code,".csv","")
    protein_feature$Region<-region_code
    protein_feature_all<-rbind(protein_feature_all,protein_feature)
    }

  }
  Peptide_Summary_file<-fread("Peptide_region_file.csv")
  Peptide_Summary_file$Source<-datafilename
  if(nrow(Peptide_Summary_file)!=0){
  Protein_peptide_Summary_file<-rbind(Protein_peptide_Summary_file,Peptide_Summary_file)
  }
  }
    message("Protein feature summary...Done.")
    if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}

    write.csv(protein_feature_all,paste(workdir,"/Summary folder/Protein_Summary.csv",sep=""),row.names = F)
    write.csv(Protein_peptide_Summary_file,paste(workdir,"/Summary folder/Protein_peptide_Summary.csv",sep=""),row.names = F)
  }
  
  #Summarize the protein and peptide result across the datafiles and store these summarized files into the summary folder
  
  if(Peptide_feature_summary){
    message("Peptide feature summary...")
    Peptide_Summary_file<-NULL
    Peptide_Summary_file_a<-NULL
    for (i in 1:length(datafile)){
      datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
      currentdir<-paste0(workdir,"/",datafile[i]," ID")

      setwd(paste(currentdir,sep=""))

      Peptide_Summary_file<-fread("Peptide_region_file.csv")
      Peptide_Summary_file$Source<-gsub(".imzML", "", datafile[i])
      if(nrow(Peptide_Summary_file)!=0){
      Peptide_Summary_file_a<-rbind(Peptide_Summary_file_a,Peptide_Summary_file)
      }
      }
    Peptide_Summary_file_a<-unique(Peptide_Summary_file_a)
      message("Peptide feature summary...Done.")
      if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
      write.csv(Peptide_Summary_file_a,paste(workdir,"/Summary folder/Peptide_Summary.csv",sep=""),row.names = F)
      #Peptide_feature_summary_all_files_new(datafile,workdir,threshold = threshold)
    }
  
  #Summarize the mz feature list 
  
  if(Region_feature_summary){
    message("Region feature summary...")
    Spectrum_summary<-NULL
    for (i in 1:length(datafile)){
      datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
      currentdir<-paste0(workdir,"/",datafile[i]," ID")
      setwd(currentdir)
      name <-gsub(base::dirname(datafile[i]),"",gsub(".imzML", "", datafile[i]))
      message(paste("Region_feature_summary",gsub(".imzML", "", datafile[i])))

      if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
      
      match_pattern <- "Spectrum.csv"
      spectrum_file_table_sum<-NULL
      for (spectrum_file in dir(recursive = T)[str_detect(dir(recursive = T), match_pattern)]){
        spectrum_file_table=fread(spectrum_file)
        if (nrow(spectrum_file_table)>=1){
        spectrum_file_table$Region=gsub("/Spectrum.csv","",spectrum_file)
        spectrum_file_table$Source<-gsub(".imzML", "", datafile[i])
        spectrum_file_table_sum[[spectrum_file]] <- spectrum_file_table
        }
        
      }
      spectrum_file_table_sum <- do.call(rbind,spectrum_file_table_sum)
      Spectrum_summary[[datafile[i]]]=spectrum_file_table_sum
    }
    Spectrum_summary <- do.call(rbind,Spectrum_summary)
    write.csv(Spectrum_summary,file = paste(workdir,"/Summary folder/Region_feature_summary.csv",sep=""),row.names = F)

  }

  #Protein cluster image rendering
  
  if(plot_cluster_image_grid){
    
    message("cluster image rendering...")
    
    setwd(workdir[1])
    
    # read protein-peptide features result
    
    Protein_feature_list=read.csv(file=paste(workdir[1],"/Summary folder/Protein_peptide_Summary.csv",sep=""),stringsAsFactors = F)
    
    # remove peptide score outlier from result
    
    if (remove_score_outlier){
      Protein_feature_list <- remove_pep_score_outlier(Protein_feature_list,abs_cutoff=Plot_score_abs_cutoff,IQR_LB = Plot_score_IQR_cutoff)
    }
    
    # extract the protein entries of interest
    
    if (sum(Protein_desc_of_interest!=".")>=1){
    Protein_feature_list_interest<-NULL
    num_of_interest<-numeric(0)
    
    for (interest_desc in Protein_desc_of_interest){
      idx_iterest_desc<-str_detect(Protein_feature_list$desc,regex(interest_desc,ignore_case = T))
      if(nrow(Protein_feature_list[idx_iterest_desc,])!=0){
      Protein_feature_list_interest<-rbind(Protein_feature_list_interest,Protein_feature_list[idx_iterest_desc,])
      }
      num_of_interest[interest_desc]<-length(unique(Protein_feature_list[idx_iterest_desc,"Protein"]))
      }
    
    Protein_feature_list=Protein_feature_list_interest
    message(paste(num_of_interest,"Protein(s) found with annotations of interest:",Protein_desc_of_interest,collapse = "\n"))
    }

    Protein_feature_list=as.data.frame(Protein_feature_list)
    
    # generate combined IMS data for multiple files or use a link to load the pre-processed IMS data
    
    if (!is.null(cluster_rds_path)){
    imdata=readRDS(paste0(workdir[1],"/",cluster_rds_path))
    message("cluster imdata loaded.")
    
    }else{
    cluster_rds_path<-Load_IMS_decov_combine(datafile=datafile,workdir=workdir,import_ppm=ppm,SPECTRUM_batch="overall",
                                       ppm=ppm,threshold=0,rotate=Rotate_IMG,mzrange=mzrange,
                                       deconv_peaklist="Load_exist",preprocessRDS_rotated=T)


    imdata=readRDS(paste0(workdir[1],"/",cluster_rds_path))
    message("cluster imdata generated and loaded.")
    }
    
    # Setup output folder and queue the R calls for cluster image randering
    outputfolder=paste(workdir,"/Summary folder/cluster Ion images/",sep="")
    if (dir.exists(outputfolder)==FALSE){dir.create(outputfolder)}

    if (!(plot_unique_component)){
    setwd(outputfolder)
    Protein_feature_list_trimmed<-Protein_feature_list
    }


    if (plot_unique_component){
    outputfolder=paste(workdir,"/Summary folder/cluster Ion images/unique/",sep="")

    if (dir.exists(outputfolder)==FALSE){dir.create(outputfolder)}
    setwd(outputfolder)
    Protein_feature_list_unique=Protein_feature_list %>% group_by(mz) %>% dplyr::summarise(num=length(unique(Protein)))
    Protein_feature_list_unique_mz<-Protein_feature_list_unique$mz[Protein_feature_list_unique$num==1]
    Protein_feature_list_trimmed<-Protein_feature_list[Protein_feature_list$mz %in% Protein_feature_list_unique_mz, ]
    write.csv(Protein_feature_list_trimmed,paste(workdir,"/Summary folder/Protein_feature_list_trimmed.csv",sep=""),row.names = F)
    }
    
    imdata <- imdata %>% peakBin(sort(unique(Protein_feature_list_trimmed$mz_align)), tolerance=ppm, units="ppm") %>% process(BPPARAM=SerialParam())
    
    save(list=c("Protein_feature_list_trimmed",
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
         file=paste0(workdir,"/cluster_img_grid.RData"), ascii = T,compress  = FALSE
         )

    for (clusterID in unique(Protein_feature_list_trimmed$Protein)){
      cluster_desc<-unique(Protein_feature_list_trimmed$desc[Protein_feature_list_trimmed[[ClusterID_colname]]==clusterID])
      cluster_desc<-gsub(stringr::str_extract(cluster_desc,"OS=.{1,}"),"",cluster_desc)
      n_component<-nrow(unique(Protein_feature_list_trimmed[Protein_feature_list_trimmed[[ClusterID_colname]]==clusterID,c(ClusterID_colname,componentID_colname,"moleculeNames","adduct","Modification")]))
      if (n_component>=peptide_ID_filter){
      if ('&'(file.exists(paste0(outputfolder,clusterID,"_cluster_imaging.png")),!plot_cluster_image_overwrite)){
        message("Cluster image rendering Skipped file exists: No.",clusterID," ",cluster_desc)
        next
        }else {

      fileConn<-file(paste0(workdir,"/cluster_img_scource.R"),)
      writeLines(c("suppressMessages(suppressWarnings(require(HiTMaP)))",
                   paste0("clusterID=",clusterID),
                   paste0("suppressMessages(suppressWarnings(load(file =\"", workdir,"/cluster_img_grid.RData\")))"),
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
                                      Component_plot_threshold=peptide_ID_filter)))"),
                 fileConn)
      close(fileConn)

      system(paste0("Rscript \"",paste0(workdir,"/cluster_img_scource.R\"")))


      if(file.exists(paste0(outputfolder,clusterID,"_cluster_imaging.png"))){
       message("Cluster image rendering Done: No.",clusterID," ",cluster_desc)

      }else{
        retrytime=1
        repeat{

          message("Cluster image rendering failed and retry ",retrytime,": No.",clusterID," ",cluster_desc)
          system(paste0("Rscript \"",paste0(workdir,"/cluster_img_scource.R\"")))

          if (file.exists(paste0(outputfolder,clusterID,"_cluster_imaging.png"))){
            message("Cluster image rendering Done: No.",clusterID," ",cluster_desc)

            break
          }else if(retrytime>=plot_cluster_image_maxretry){
            message("Cluster image rendering reaches maximum Retry Attempts: No.",clusterID," ",cluster_desc)

            break
          }
          retrytime=1+retrytime
        }

      }

      }


      }
      }

  
}
message("Workflow done.")
}

IMS_data_process<-function(datafile,
                                    workdir=NULL,
                                    Peptide_Summary_searchlist,
                                    segmentation_num=5,threshold=0.1,
                                    ppm,import_ppm=5,
                                    mzrange="auto-detect",
                                    Segmentation=c("spatialKMeans","spatialShrunkenCentroids","Virtual_segmentation","none","def_file"),
                                    Segmentation_def="segmentation_def.csv",
                                    Segmentation_ncomp="auto-detect",
                                    Segmentation_variance_coverage=0.8,
                                    Bypass_Segmentation=F,
                                    Smooth_range=1,
                                    colorstyle="Set1",
                                    Virtual_segmentation_rankfile=NULL,
                                    PMFsearch=TRUE,
                                    rotate=NULL,
                                    matching_validation=T,
                                    BPPARAM=bpparam(),
                                    Bypass_generate_spectrum=F,
                                    scoring=T,
                                    score_method="SQRTP",
                                    Rank=3,
                                    Decoy_mode="",
                                    Decoy_search=T,
                                    adjust_score=T,
                                    plot_matching_score_t=F,
                                    Protein_feature_list,
                                    peptide_ID_filter=2,
                                    Protein_desc_of_interest=".",
                                    FDR_cutoff=0.1,
                                    use_top_rank=NULL,
                                    preprocess=list(force_preprocess=FALSE,use_preprocessRDS=TRUE,smoothSignal=list(method="gaussian"),
                                                    reduceBaseline=list(method="locmin"),
                                                    peakPick=list(method="adaptive"),
                                                    peakAlign=list(tolerance=5, units="ppm"),
                                                    normalize=list(method=c("rms","tic","reference")[1],mz=1)),
                                    ...){
   suppressMessages(suppressWarnings(require(data.table)))
   suppressMessages(suppressWarnings(require(Cardinal)))
   suppressMessages(suppressWarnings(require(RColorBrewer)))
   suppressMessages(suppressWarnings(require(stringr)))
   getPalette = colorRampPalette(brewer.pal_n(9, colorstyle))
   setCardinalBPPARAM(BPPARAM)

   # resolve the filename, project folder and rotation info
   datafile <- paste0(workdir,"/",datafile)
   workdir <- dirname(datafile)
   datafile <- basename(datafile)
   rotate=Parse_rotation(datafile,rotate)
   datafile_imzML<-datafile

   # perform the IMS analysis
  for (z in 1:length(datafile)){
    
    #file related variables
    name <-basename(datafile[z])
    name <-gsub(".imzML$","",name)
    name <-gsub("/$","",name)
    folder<-base::dirname(datafile[z])
    if (!str_detect(datafile[z],".imzML$")){
      datafile_imzML[z]<-paste0(datafile[z],".imzML")
    }
    
    #perform the IMS pre-processing and image segmentation
         setwd(workdir[z])
         segmentation_res <- Preprocessing_segmentation(datafile=datafile[z],
                                         workdir=workdir,
                                         segmentation_num=segmentation_num,
                                         ppm=ppm,import_ppm=import_ppm,Bypass_Segmentation=Bypass_Segmentation,
                                         mzrange=mzrange,
                                         Segmentation=Segmentation,
                                         Segmentation_def=Segmentation_def,
                                         Segmentation_ncomp=Segmentation_ncomp,
                                         Segmentation_variance_coverage=Segmentation_variance_coverage,
                                         Smooth_range=Smooth_range,
                                         colorstyle=colorstyle,
                                         Virtual_segmentation_rankfile=Virtual_segmentation_rankfile,
                                         rotate=rotate,
                                         BPPARAM=BPPARAM,
                                         preprocess=preprocess)

         segmentation_label=segmentation_res$segmentation_label
         
         imdata=segmentation_res$imdata
         #imdata_ed=segmentation_res$imdata_ed
         #imdata_org=segmentation_res$imdata_org
         rm(segmentation_res)
         
   if(PMFsearch){
     if (missing(Protein_feature_list)){
       Protein_feature_list=get("Protein_feature_list", envir = .GlobalEnv)
       message("Got Protein_feature_list from global environment.")
     }
     Peptide_Summary_searchlist$Protein<-NULL
     Peptide_Summary_searchlist$pro_end<-NULL
     Peptide_Summary_searchlist$start<-NULL
     Peptide_Summary_searchlist$end<-NULL
     Peptide_Summary_searchlist<-unique(Peptide_Summary_searchlist)
     Peptide_Summary_file<-Peptide_Summary_searchlist
     Peptide_Summary_file$Intensity<-rep(0,nrow(Peptide_Summary_file))

   if (ppm>=25) {
    instrument_ppm=50
   }else{
    instrument_ppm=8
   }
   setwd(workdir[z])
   
   #remove previously genrated IMS annotation result
   if (dir.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID"))==FALSE){dir.create(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
      }else{
        all_files<-dir(paste0(gsub(".imzML$","",datafile[z])  ," ID"),recursive = F)
        delete_files<-all_files[str_detect(all_files,"Protein_segment_PMF_RESULT_|Peptide_|PMF spectrum match.png$")]
        delete_dir<-list.dirs(path = paste0(gsub(".imzML$","",datafile[z])  ," ID"), full.names = TRUE, recursive = F)
        try(suppressWarnings(file.remove(c(paste0(gsub(".imzML$","",datafile[z])  ," ID/",delete_files)))))
        unlink(delete_dir, recursive = T)
      }
   
   
   #Start annotation for each found region
    setwd(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
    Peptide_Summary_file_regions<-data.frame()
    message(paste("PMFsearch",name))
    message(paste( "region",names(segmentation_label),"Found.",sep=" ",collapse = "\n"))
    for (SPECTRUM_batch in names(segmentation_label)){
      if (dir.exists(paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/"))==FALSE){dir.create(paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/"))}
    # if (!is.null(preprocess)){
    #   if ('|'(imdata@metadata[["ibd binary type"]]!="processed",preprocess$force_preprocess)){
    #   message(paste("Using preprocessed .rda data:",datafile[z]))
    #   imdata<-imdata_ed
    #   } else if (imdata@metadata[["ibd binary type"]]=="processed"){
    #   message(paste("Using preprocessed .imzml data:",datafile[z]))
    #   imdata<-imdata_org
    #   }
    # }
      
    imdata_sb<-imdata[,unlist(segmentation_label[[SPECTRUM_batch]])]
   #generate spectrum for each found region
    spectrum_file_table<- summarizeFeatures(imdata_sb, FUN = "mean")
    spectrum_file_table<-data.frame(mz=spectrum_file_table@featureData@mz,mean=spectrum_file_table@featureData@listData[["mean"]])
    peaklist<-spectrum_file_table
    colnames(peaklist)<-c("m.z","intensities")
    savename=paste(name,SPECTRUM_batch)
    message(paste("IMS_analysis",name,"region",SPECTRUM_batch))
    
    peaklist<-peaklist[peaklist$intensities>0,]
    write.csv(peaklist,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Spectrum.csv"),row.names = F)
    
   #generate filtered processed peaklist to next step
    deconv_peaklist<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist,ppm=ppm,threshold=0)
    peaklist_pmf<-deconv_peaklist[deconv_peaklist$intensities>(max(deconv_peaklist$intensities)*threshold),]
    message(paste(nrow(peaklist_pmf),"mz features found in the spectrum") )
    
   #Do first round of peptide search to get putative result
    mz_feature_list<-Do_PMF_search(peaklist_pmf,Peptide_Summary_searchlist,BPPARAM=BPPARAM,ppm = ppm)
    mz_feature_list<-unique(mz_feature_list)
    message("Iterating peptide information...")
    Peptide_Summary_searchlist<-as.data.table(Peptide_Summary_searchlist)
    mz_feature_list<-as.data.table(mz_feature_list)
    Peptide_Summary_searchlist$Intensity<-NULL
    mz_feature_list$mz<-as.character(mz_feature_list$mz)

    Peptide_Summary_searchlist$mz<-as.character(Peptide_Summary_searchlist$mz)
    Peptide_Summary_searchlist<-merge(Peptide_Summary_searchlist,mz_feature_list,by.x="mz",by.y="mz",all.x=T,sort=F)
    Peptide_Summary_searchlist$Intensity[is.na(Peptide_Summary_searchlist$Intensity)]<-0
    Peptide_plot_list<-Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity>0,]
    Peptide_plot_list$formula<-as.character(Peptide_plot_list$formula)
  
    if(is.null(Peptide_plot_list$moleculeNames)){Peptide_plot_list$moleculeNames=Peptide_plot_list$Peptide}
    Peptide_plot_list$Region=SPECTRUM_batch
    Peptide_plot_list=Peptide_plot_list[(!is.na(Peptide_plot_list$Intensity)),]
    message(paste("1st run returns",nrow(Peptide_plot_list), "peptide candidates"))
    write.csv(Peptide_plot_list,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_1st_ID.csv"),row.names = F)
    
    
    #perform peptide-protein scoring and FDR cut off
    if (nrow(Peptide_plot_list)==0){
      next
    }
    if (scoring){
      
      suppressMessages(suppressWarnings(require(enviPat)))
      suppressMessages(suppressWarnings(require(ggplot2)))
      
      data(isotopes)

      unique_formula<-unique(Peptide_plot_list$formula)
      #Peptide_plot_list_Score=(bplapply(unique_formula,SCORE_PMF,peaklist=peaklist,isotopes=isotopes,score_method=score_method,charge = 1,ppm=ppm,BPPARAM = BPPARAM))
      Peptide_plot_list_Score=(lapply(unique_formula,SCORE_PMF,peaklist=peaklist,isotopes=isotopes,score_method=score_method,charge = 1,ppm=ppm))
      
      Peptide_plot_list_Score_m=as.data.frame(do.call(rbind, Peptide_plot_list_Score))
      names(Peptide_plot_list_Score_m)<-c("Score", "delta_ppm","Intensity")
      formula_score<-data.frame(formula=unique_formula,Score=Peptide_plot_list_Score_m$Score,Delta_ppm=Peptide_plot_list_Score_m$delta_ppm,Intensity=Peptide_plot_list_Score_m$Intensity)

      Peptide_plot_list$Score<-NULL
      Peptide_plot_list$Delta_ppm<-NULL
      Peptide_plot_list$Intensity<-NULL
      Peptide_plot_list<-merge(Peptide_plot_list,formula_score,by="formula")
    
      #generate decoy candidate list if running in "isotope" decoy mode 
      if (Decoy_search && ("isotope" %in% Decoy_mode)){
      decoy_isotopes=isotopes
      decoy_isotopes[decoy_isotopes$isotope=="13C",]=data.frame(element="C",isotope="11C",mass=10.99664516,aboudance=0.0107,ratioC=0,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="15N",]=data.frame(element="N",isotope="13N",mass=13.00603905,aboudance=0.00364000,ratioC=4,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="2H",]=data.frame(element="H",isotope="0H",mass=0.001548286,aboudance=0.00011500,ratioC=6,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="17O",]=data.frame(element="O",isotope="15O",mass=14.99069774,aboudance=0.00038000,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="18O",]=data.frame(element="O",isotope="14O",mass=13.99066884,aboudance=0.00205000,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="33S",]=data.frame(element="S",isotope="31S",mass=30.97268292,aboudance=0.0075,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="34S",]=data.frame(element="S",isotope="30S",mass=29.9762745,aboudance=0.0425,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="35S",]=data.frame(element="S",isotope="29S",mass=28.94414146,aboudance=0,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="36S",]=data.frame(element="S",isotope="28S",mass=27.97706058,aboudance=0.0001,ratioC=3,stringsAsFactors = F)
      Peptide_plot_list_decoy<-Peptide_plot_list[Peptide_plot_list$isdecoy==0,]
      Peptide_plot_list_decoy$isdecoy=1
      unique_formula<-unique(Peptide_plot_list_decoy$formula)
      #Peptide_plot_list_Score=(bplapply(unique_formula,SCORE_PMF,peaklist=peaklist,isotopes=decoy_isotopes,score_method=score_method,charge = 1,ppm=ppm,BPPARAM = BPPARAM))
      Peptide_plot_list_Score=(lapply(unique_formula,SCORE_PMF,peaklist=peaklist,isotopes=decoy_isotopes,score_method=score_method,charge = 1,ppm=ppm))
      Peptide_plot_list_Score_m=as.data.frame(do.call(rbind, Peptide_plot_list_Score))
      names(Peptide_plot_list_Score_m)<-c("Score", "delta_ppm","Intensity")
      formula_score<-data.frame(formula=unique_formula,Score=Peptide_plot_list_Score_m$Score,Delta_ppm=Peptide_plot_list_Score_m$delta_ppm,Intensity=Peptide_plot_list_Score_m$Intensity)
      Peptide_plot_list_decoy$Score<-NULL
      Peptide_plot_list_decoy$Delta_ppm<-NULL
      Peptide_plot_list_decoy$Intensity<-NULL
      Peptide_plot_list_decoy<-merge(Peptide_plot_list_decoy,formula_score,by="formula")
      Peptide_plot_list<-rbind(Peptide_plot_list,Peptide_plot_list_decoy)
      Protein_feature_list_decoy<-Protein_feature_list[Protein_feature_list$isdecoy==0,]
      Protein_feature_list_decoy$isdecoy=1
      Protein_feature_list=rbind(Protein_feature_list[Protein_feature_list$isdecoy==0,],Protein_feature_list_decoy)
      }
      Peptide_plot_list$mz<-as.numeric(as.character(Peptide_plot_list$mz))
      Peptide_plot_list_rank=rank_mz_feature(Peptide_plot_list,mz_feature=deconv_peaklist,BPPARAM = BPPARAM)
      Peptide_plot_list_rank<-unique(Peptide_plot_list_rank)
      write.csv(Peptide_plot_list_rank,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_1st_ID_score_rank_",score_method,".csv"),row.names = F)
      Peptide_plot_list_rank<-read.csv(paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_1st_ID_score_rank_",score_method,".csv"),stringsAsFactors = F)
      
      #Peptide score filtering (default=False) 
      if (adjust_score==F){
      Score_cutoff= FDR_cutoff_plot(Peptide_plot_list_rank,FDR_cutoff=FDR_cutoff,plot_fdr=T,outputdir=paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch),adjust_score = adjust_score)
      Peptide_plot_list_2nd=Peptide_plot_list_rank[((Peptide_plot_list_rank$Score>=Score_cutoff)&(!is.na(Peptide_plot_list_rank$Intensity))),]
      }else{
        Score_cutoff= FDR_cutoff_plot(Peptide_plot_list_rank,FDR_cutoff=FDR_cutoff,plot_fdr=T,outputdir=paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch),adjust_score = adjust_score)
        Peptide_plot_list_2nd=Score_cutoff[[2]]
        Score_cutoff=Score_cutoff[[1]]
        Peptide_plot_list_2nd=Peptide_plot_list_2nd[((Peptide_plot_list_2nd$Score>=Score_cutoff)&(!is.na(Peptide_plot_list_2nd$Intensity))),]
        }
      write.csv(Peptide_plot_list_2nd,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_2nd_ID_score_rank",score_method,"_Rank_above_",Rank,".csv"),row.names = F)
      
      #Plot peptide matching spectrum and score 
      if (plot_matching_score_t && nrow(Peptide_plot_list_2nd)!=0){
        try(plot_matching_score(Peptide_plot_list_2nd,peaklist,charge=1,ppm,outputdir=paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/ppm")))
      }
      Peptide_plot_list_rank$mz=as.numeric(as.character(Peptide_plot_list_rank$mz))
      Peptide_plot_list_2nd$mz=as.numeric(as.character(Peptide_plot_list_2nd$mz))

      #Protein scoring and FDR cut-off
      Index_of_protein_sequence<-get("Index_of_protein_sequence", envir = .GlobalEnv)
      Protein_feature_result<-protein_scoring(Protein_feature_list,Peptide_plot_list_rank,BPPARAM = BPPARAM,scoretype="mean",peptide_ID_filter=peptide_ID_filter,use_top_rank=use_top_rank)
      Protein_feature_list_rank<-Protein_feature_result[[2]]
      Protein_feature_list_rank_filtered<-Protein_feature_result[[3]]
      Protein_feature_list_rank_filtered_grouped<-Protein_feature_result[[4]]
      Protein_feature_result<-Protein_feature_result[[1]]


      if (nrow(Protein_feature_result)>=2){
      Score_cutoff_protein= FDR_cutoff_plot_protein(Protein_feature_result,FDR_cutoff=FDR_cutoff,plot_fdr=T,outputdir=paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch),adjust_score = F)
      }else{
            Score_cutoff_protein=Protein_feature_result$Proscore
            if (length(Score_cutoff_protein)==0) Score_cutoff_protein=0}
      Protein_feature_result_cutoff=Protein_feature_result[((Protein_feature_result$Proscore>=Score_cutoff_protein)&(!is.na(Protein_feature_result$Intensity))&(Protein_feature_result$isdecoy==0)),]
      Protein_feature_list_rank=Protein_feature_list_rank[Protein_feature_list_rank$Protein %in% Protein_feature_result_cutoff$Protein,]
      Protein_feature_list_rank$Score=round(Protein_feature_list_rank$Score,digits = 7)
      Protein_feature_list_rank$mz=round(Protein_feature_list_rank$mz,digits = 4)
      Protein_feature_list_rank<-as.data.frame(Protein_feature_list_rank)
      Protein_feature_list_rank$desc<-NULL
      Protein_feature_list_rank=merge(Protein_feature_list_rank,Index_of_protein_sequence[,c("recno","desc")],by.x="Protein",by.y="recno",all.x=T)

      Protein_feature_list_rank_cutoff<-Protein_feature_list_rank
      Protein_feature_list_rank_cutoff<-Protein_feature_list_rank_cutoff[Protein_feature_list_rank_cutoff$isdecoy==0,]
      write.csv(Protein_feature_result,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Protein_ID_score_rank_",score_method,".csv"),row.names = F)
      write.csv(Protein_feature_list_rank_filtered,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Protein_ID_score_rank_filtered_",score_method,".csv"),row.names = F)
      write.csv(Protein_feature_list_rank_filtered_grouped,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Protein_ID_score_rank_filtered_grouped_",score_method,".csv"),row.names = F)
      write.csv(Protein_feature_list_rank_cutoff,paste0(workdir[z],"/",datafile[z] ," ID/","Peptide_segment_PMF_RESULT_",SPECTRUM_batch,".csv"),row.names = F)

      #Plot Protein matching spectrum and score 
      png(paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/unique_peptide_ranking_vs_mz_feature",".png"),width = 960,height = 480)
      sp<-ggplot2::ggplot(Peptide_plot_list_rank, aes(x=mz,fill=as.factor(Rank))) +  geom_bar(stat = "bin",bins = 100) +
        labs(title="Matched peptide ranking vs. mz feature",x="mz", y = "Matched peptide ranking")+
        theme_classic() + theme(legend.title  = element_blank()) + scale_fill_brewer()
      try(print(sp))
      dev.off()
      message(paste("Plotting peptide matching number vs mz feature"))
      Peptide_plot_list_2nd=Protein_feature_list_rank_cutoff
      png(paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/unique_peptide_ranking_vs_mz_feature_2nd",".png"),width = 960,height = 480)
      sp<-ggplot2::ggplot(Peptide_plot_list_2nd, aes(x=mz,fill=as.factor(floor(Rank/max(Rank)*8)))) +  geom_bar(stat = "bin",bins = 100) +
        labs(title="Matched peptide ranking 2nd vs. mz feature",x="mz", y = "Matched peptide ranking")+
        theme_classic() + theme(legend.title  = element_blank()) + scale_fill_brewer()
      try(print(sp))
      dev.off()
      write.csv(Protein_feature_result_cutoff,paste0(workdir[z],"/",datafile[z] ," ID/","Protein_segment_PMF_RESULT_",SPECTRUM_batch,".csv"),row.names = F)
    }
    
    #Plot overall matching spectrum and score 
    Plot_PMF_all(Peptide_plot_list_2nd,peaklist=peaklist,threshold=threshold,savename)
    write.csv(peaklist,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"/Spectrum.csv"),row.names = F)
    if (nrow(Peptide_plot_list_2nd)!=0){
    Peptide_plot_list_2nd$Region<-SPECTRUM_batch
    }
    Peptide_Summary_file_regions<-rbind(Peptide_Summary_file_regions,Peptide_plot_list_2nd)

    }

  if(is.null(Peptide_Summary_file_regions$Peptide)){Peptide_Summary_file_regions$Peptide=Peptide_Summary_file_regions$moleculeNames}
  if(is.null(Peptide_Summary_file_regions$moleculeNames)){Peptide_Summary_file_regions$moleculeNames=Peptide_Summary_file_regions$Peptide}
  
  #export final result for the given data 
  write.csv(Peptide_Summary_file_regions,"Peptide_region_file.csv",row.names = F)


  }

  }

  setwd(workdir[1])
}






