#initiallization

require("pacman")
library(tcltk)

Filters <- matrix(c( "imzml file", ".imzML",
                     "Text", ".txt", "All files", "*"),
                  3, 2, byrow = TRUE)




			   
#' imaging_identification
#'
#' This is a peptide mass fingerprint search function for maldi imaging data analysis
#'
#' @param threshold specify the intensities threshold (0 to 1 in percentage)to report a identified molecule 
#' @param ppm the mz tolerance (in ppm) for peak integration
#' @param Digestion_site the enzyme digestion site specificity
#' @param missedCleavages miss cleavage number allowed in this PMF search
#' @param Fastadatabase the fasta database used in this pmf search, the file should be placed in the same folder with data files
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param PMF_analysis Set \code{"true"} if you want to have a PMF search, set \code{"false"} if you want to bypass it
#' @param Protein_feature_summary  \code{"PMF_analysis"} follow-up process that will collect all the identified peptide information and associate them with possible proteins 
#' @param plot_cluster_image  \code{"Protein_feature_summary"} follow-up process that will plot the protein cluster image 
#' @param Peptide_feature_summarya \code{"PMF_analysis"} follow-up process that will summarize all datafiles identified peptides and generats a \code{"peptide shortlist"} in the result summary folder
#' @param plot_ion_image  \code{"Peptide_feature_summarya"} follow-up process that will plot every connponents in the \code{"peptide shortlist"}
#' @param parallel the number of threads will be used in the PMF search, this option now only works for windows OS
#' @param spectra_segments_per_file optimal number of distinctive regions in the imaging, a virtual segmentation will be applied to the image files with this value. To have a better PMF result you may set a value that in the sweet point of sensitivety and false discovery rate (FDR).
#' @param spatialKMeans set true to enable a \code{"spatialKMeans"}  method for the automatic virtual segmentation. If a region rank file was supplied, you can disable this to perform a mannual segmentation.
#' @param Smooth_range \code{"spatialKMeans"} pixel smooth range 
#' @param Virtual_segmentation set \code{"TRUE"} if you want to overide the automaitic segmentation
#' @param Virtual_segmentation_rankfile specify a region rank file contains region information for manualy region segmentation
#' @return None
#'
#' @examples
#' imaging_identification(threshold=0.05, ppm=5,Digestion_site="[G]",
#'                        missedCleavages=0:1,Fastadatabase="murine_matrisome.fasta",
#'                        adducts=c("M+H","M+NH4","M+Na"),PMF_analysis=TRUE,
#'                        Protein_feature_summary=TRUE,plot_cluster_image=TRUE,
#'                        Peptide_feature_summary=TRUE,plot_ion_image=FALSE,
#'                        parallel=3,spectra_segments_per_file=5,spatialKMeans=TRUE
#'                        )
#'
#' @export
#' 
imaging_identification<-function(
#==============Choose the imzml raw data file(s) to process  make sure the fasta file in the same folder
               datafile,
               threshold=0.001, 
               ppm=5,
               mode=c("Proteomics","Metabolomics"),
               Digestion_site="[GN][LI]",
               missedCleavages=0:1,
               Fastadatabase="murine_matrisome.fasta",
               adducts=c("M+H","M+NH4","M+Na"),
               Decoy_search=T,
               Decoy_adducts=c("M+ACN+H","M+IsoProp+H","M+DMSO+H","M+Co","M+Ag","M+Cu","M+He","M+Ne","M+Ar","M+Kr","M+Xe","M+Rn"),
               Decoy_mode = "isotope",
               PMF_analysis=TRUE,
               Bypass_generate_spectrum=F,
               Protein_feature_summary=TRUE,
               plot_cluster_image=TRUE,
               Peptide_feature_summary=TRUE,
               plot_ion_image=FALSE,
               parallel=detectCores(),
               spectra_segments_per_file=5,
               spatialKMeans=TRUE,
               Smooth_range=1,
               Virtual_segmentation=FALSE,
               Virtual_segmentation_rankfile=NULL,
               rotateimg=NULL,
               Region_feature_summary=F,
               Spectrum_validate=T,
               output_candidatelist=T,
               use_previous_candidates=F,
               score_method="SQRTP",
               ...
               ){
  library("pacman")
  p_load(RColorBrewer,RCurl,bitops,magick,ggplot2,reticulate,dplyr,stringr,tcltk,
         data.table,doParallel,iterators,foreach,protViz,cleaver,MALDIquant,Biostrings,
         XVector,IRanges,Cardinal,ProtGenerics,S4Vectors,stats4,EBImage,BiocParallel,
         BiocGenerics,parallel,stats,graphics,grDevices,utils,datasets,methods)
  if (missing(datafile)) {datafile=tk_choose.files(filter = matrix(c( "imzml file", ".imzML","Text", ".txt", "All files", "*"),3, 2, byrow = TRUE),
                                        caption  = "Choose single or multiple file(s) for analysis")}
  datafile<-gsub(".imzML", "", datafile)
  workdir<-base::dirname(datafile[1])
  setwd(workdir)
  #cl <- autoStopCluster(makeCluster(parallel))
  parallel=try(detectCores()/2)
  if (parallel<1 | is.null(parallel)){parallel=1}
  BPPARAM=bpparam()
  BiocParallel::bpworkers(BPPARAM)=parallel
  bpprogressbar(BPPARAM)=TRUE
  setwd(workdir)
  message(paste(try(detectCores()), "Cores detected,",parallel, "threads will be used for computing"))

  message(paste(length(datafile), "files were selected and will be used for Searching"))
  
  message(paste(Fastadatabase, "was selected as data base", "candidates will be generated through",mode[1] ,"mode" ))
  
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
                                                 use_previous_candidates=use_previous_candidates
                                                 )

  
  if (!is.null(rotateimg)){rotateimg=read.csv(rotateimg,stringsAsFactors = F)}
  
  if(PMF_analysis){
  message(paste("Protein_feature_list was selected as database",spectra_segments_per_file,threshold,ppm,PMF_analysis,Virtual_segmentation,
                "\nVirtual_segmentation_rankfile",Virtual_segmentation_rankfile,"\nBypass_generate_spectrum",Bypass_generate_spectrum))
  Peptide_Summary_searchlist<-unique(Protein_feature_list[,c("Peptide","mz","adduct","formula","isdecoy")])
  
  Peptide_Summary_file<-PMF_Cardinal_Datafilelist(datafile, 
                                                  Peptide_Summary_searchlist,
                                                  SPECTRUM_for_average=spectra_segments_per_file,
                                                  threshold=threshold,rotate = rotateimg,
                                                  ppm=ppm,
                                                  spatialKMeans=spatialKMeans,
                                                  PMFsearch = PMF_analysis,
                                                  Virtual_segmentation=Virtual_segmentation,
                                                  Virtual_segmentation_rankfile = Virtual_segmentation_rankfile,
                                                  BPPARAM = BPPARAM,
                                                  Bypass_generate_spectrum=Bypass_generate_spectrum,
                                                  score_method = score_method,
                                                  Protein_feature_list=Protein_feature_list,
                                                  Decoy_mode=Decoy_mode,
                                                  Decoy_search=Decoy_search)
  
  }
  #Summarize the peptide list
  #Summarize the protein and peptide list across the datafiles
  if(Protein_feature_summary){
  for (i in 1:length(datafile)){
  datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
  currentdir<-paste0(datafile[i] ," ID")
  
  setwd(paste(currentdir,sep=""))
  
  Peptide_Summary_file<-fread("Peptide_Summary_file.csv",select=c("Peptide","mz","Intensity","adduct","moleculeNames"))
  Peptide_feature_list<-Peptide_Summary_file[Peptide_Summary_file$Intensity>=max(Peptide_Summary_file$Intensity)*threshold,]
  write.csv(Peptide_feature_list,"Peptide_feature_list.csv")
  uniques_intensity<-unique(Peptide_Summary_file[,c("mz","Intensity")])
  Protein_feature_list$Intensity<-NULL
  #Protein_feature_list$Intensity<-unlist(parLapply(cl=autoStopCluster(makeCluster(parallel)),Protein_feature_list$mz,intensity_sum_para,uniques_intensity))
  message("Iterating protein information")
  Protein_feature_list=merge(Protein_feature_list,uniques_intensity,by="mz")
  
  #Protein_feature_list$Intensity<-unlist(bplapply(Protein_feature_list$mz,intensity_sum_para,uniques_intensity,BPPARAM=BPPARAM))
  Protein_feature_list=Protein_feature_list[Protein_feature_list$Intensity>(threshold*max(Protein_feature_list$Intensity)),]
  write.csv(Protein_feature_list,"Cluster.csv")
  } 
    if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
    Protein_feature_summary_all_files_new(datafile,workdir,threshold = threshold)
  }
  
 
    
  if(Peptide_feature_summary){
      print("Peptide_feature_summary")
      if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
      Peptide_feature_summary_all_files_new(datafile,workdir,threshold = threshold)
      }

  if(plot_ion_image){
    print("plot_ion_image")
    if (dir.exists(paste(workdir,"/Ion images",sep=""))==FALSE){dir.create(paste(workdir,"/Ion images",sep=""))}
  setwd(paste(workdir,"/Ion images",sep=""))
  masslist<-read.csv(paste(workdir,"/Summary folder/Peptide_feature_summary_sl.csv",sep = ""))
  masslist<-unique(masslist[,c("Peptide","mz","adduct","moleculeNames")],by="moleculeNames")
  write.csv(masslist,paste(workdir,"/Ion images/Peptide_feature_summary_sl.csv",sep = ""))
  #simple_ion_image_cardinal(datafile,workdir=workdir,plotTolerance=ppm,interpolate =TRUE)
  simple_ion_image_cardinal(datafile=datafile,
                                      #workdir=WorkingDir(),
                                      Image_Type="",
                                      plotTolerance=ppm ,
                                      creat_new_file=TRUE,
                                      color = "black",
                                      smooth.image = "adaptive",
                                      contrast.enhance = "histogram")
  }

  if(Region_feature_summary){
    Spectrum_summary<-NULL
    for (i in 1:length(datafile)){
      datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
      currentdir<-paste0(datafile[i] ," ID")
      setwd(currentdir)
      name <-gsub(base::dirname(datafile[i]),"",datafile[i])
      message(paste("Region_feature_summary",datafile[i]))
      
      if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
      #Feature_summary_all_files(datafile,workdir,threshold = threshold)
      match_pattern <- "Peptide_region_file.csv"
      
      for (spectrum_file in dir()[str_detect(dir(), match_pattern)]){
        spectrum_file_table=fread(spectrum_file) 
        
        spectrum_file_table$ID=paste(spectrum_file_table$moleculeNames,spectrum_file_table$adduct,spectrum_file_table$mz,spectrum_file_table$Region,sep = "@")
        spectrum_file_table=spectrum_file_table[,c( "ID",	"Intensity")]
        colnames(spectrum_file_table)=c("ID",paste(name))
        if (is.null(Spectrum_summary)){
          Spectrum_summary=spectrum_file_table
        }else{
          Spectrum_summary=base::merge(Spectrum_summary,spectrum_file_table,by="ID",all=TRUE)
        }
      }
    }
    
    name <-gsub(base::dirname(datafile[1]),"",datafile)
    Spectrum_summary[is.na(Spectrum_summary)] <- 0
    colnames(Spectrum_summary)=c("ID",name)
    colnames(Spectrum_summary)=gsub("/","",colnames(Spectrum_summary))
    write.csv(Spectrum_summary,file = paste(workdir,"/Summary folder/Region_feature_summary.csv",sep=""),row.names = F)
    
    
  }
  
  if(Spectrum_validate){
    install.packages("OrgMassSpecR")
    require(OrgMassSpecR)
    Protein_feature_summary<-read.csv(paste(workdir,"/Summary folder",sep=""))
    
    SpectrumSimilarity
  }
  if(plot_cluster_image_grid){
    Protein_feature_list=fread(paste(workdir,"/Summary folder/Protein_feature_summary_sl.csv",sep=""))
    if (!is.null(rotateimg)){rotateimg=read.csv(rotateimg,stringsAsFactors = F)}
    imdata=list()
    combinedimdata=NULL
    #register(SerialParam())      
    if (!exists("mzrange")){
        mzrange=NULL
        testrange=c(0,0)
        for (i in 1:length(datafile)){
          
          rotate=rotateimg[rotateimg$filenames==datafile[i],"rotation"]
          rotate=as.numeric(rotate)
          
          if (length(rotate)==0){rotate=0}
          imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm*2,rotate = rotate,attach.only=F,as="MSImagingExperiment",mzrange=mzrange)
          #imdata[[i]]@elementMetadata@coord=imdata[[i]]@elementMetadata@coord[,c("x","y")]
          
          
          if (i==1) {
            testrange=c(min(imdata[[i]]@featureData@mz),max(imdata[[i]]@featureData@mz))
          }else{
            
            if (min(imdata[[i]]@featureData@mz)>testrange[1]) testrange[1]<-min(imdata[[i]]@featureData@mz)
            if (max(imdata[[i]]@featureData@mz)>testrange[2]) testrange[2]<-max(imdata[[i]]@featureData@mz)
          }
          imdata[[i]]=NULL
        }
        mzrange=testrange
    } 
    
    for (i in 1:length(datafile)){

      rotate=rotateimg[rotateimg$filenames==datafile[i],"rotation"]
      rotate=as.numeric(rotate)
      
      if (length(rotate)==0){rotate=0}
      imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm*2,rotate = rotate,attach.only=F,as="MSImagingExperiment",mzrange=mzrange)
      #imdata[[i]]@elementMetadata@coord=imdata[[i]]@elementMetadata@coord[,c("x","y")]
      max(imdata[[i]]@featureData@mz)
      min(imdata[[i]]@featureData@mz)
      if (i==1) {
        combinedimdata=imdata[[i]]
      }else{
        combinedimdata=cbind(combinedimdata,imdata[[i]]) 
      }
      imdata[[i]]=NULL
      
    }
      
    
    combinedimdata@elementMetadata@coord@listData[["z"]]=NULL
    
    imdata=combinedimdata
    
    outputfolder=paste(workdir,"/Summary folder/cluster Ion images/",sep="")
    
    if (dir.exists(outputfolder)==FALSE){dir.create(outputfolder)}
    setwd(outputfolder)
    
    lapply(unique(Protein_feature_list$Protein),
           cluster_image_grid,
           imdata=imdata,
           SMPLIST=Protein_feature_list,
           ppm=ppm,ClusterID_colname="Protein",
           componentID_colname="Peptide",
           plot_layout="line",
           Component_plot_threshold=4,
           export_Header_table=F)
    
    
    lapply(unique(Protein_feature_list$Protein),
           cluster_image_grid,
           imdata=NULL,
           SMPLIST=Protein_feature_list,
           ppm=ppm,ClusterID_colname="Protein",
           componentID_colname="Peptide",
           plot_layout="line",
           Component_plot_threshold=4,
           export_Header_table=T)
    
    
    
    Pngclusterkmean=NULL
    Pngclustervseg=NULL
    for (i in 1:length(datafile)){
      #imdata=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm)
      currentdir<-paste0(datafile[i] ," ID")
      setwd(paste(currentdir))
      name=gsub(base::dirname(datafile[i]),"",datafile[i])
      name=gsub("/","",name)
      pngimagekmean=magick::image_read("spatialKMeans_image.png")
      pngimagevseg=magick::image_read(paste0("Virtual_segmentation ",name,".png"))
      if( is.null(Pngclusterkmean)){
        Pngclusterkmean=pngimagekmean
        Pngclustervseg=pngimagevseg
      }else{
        Pngclusterkmean=magick::image_append(c(Pngclusterkmean,pngimagekmean),stack = T)
        Pngclustervseg=magick::image_append(c(Pngclustervseg,pngimagevseg),stack = T)
      }
    } 
    image_write(Pngclusterkmean, paste0(outputfolder,"datafiles_kmean.png"))
    image_write(Pngclustervseg, paste0(outputfolder,"datafiles_vseg.png"))
  }
  
  
}




#' imaging_Spatial_Quant
#'
#' This is a spatial quantitation function for maldi imaging data set
#' this function will read the candidate list file and generate quantification result
#' @param datafile
#' @param threshold specify the intensities threshold (0 to 1 in percentage)to report a identified molecule 
#' @param ppm the mz tolerance (in ppm) for peak integration
#' @param Quant_list the quantifiaction candidate list, spatial quantification will go through every datafile and collect the ion intensities for each listed component
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param cal.mz If set with \code{"true"}, the function will recalculate the mz value according to the column named "formular" in the \code{Quant_list} and the specified adducts.
#' @param mzlist_bypass  Set \code{"true"} if you want to bypass the mzlist generating process
#' @param Protein_feature_summary  \code{"PMF_analysis"} follow-up process that will collect all the identified peptide information and associate them with possible proteins 
#' @param plot_cluster_image  \code{"Protein_feature_summary"} follow-up process that will plot the protein cluster image 
#' 
#' @param Peptide_feature_summarya \code{"PMF_analysis"} follow-up process that will summarize all datafiles identified peptides and generats a \code{"peptide shortlist"} in the result summary folder
#' @param plot_ion_image  \code{"Peptide_feature_summarya"} follow-up process that will plot every connponents in the \code{"peptide shortlist"}
#' @param parallel the number of threads will be used in the PMF search, this option now only works for windows OS
#' @param spectra_segments_per_file optimal number of distinctive regions in the imaging, a virtual segmentation will be applied to the image files with this value. To have a better PMF result you may set a value that in the sweet point of sensitivety and false discovery rate (FDR).
#' @param spatialKMeans set true to enable a \code{"spatialKMeans"}  method for the automatic virtual segmentation. If a region rank file was supplied, you can disable this to perform a mannual segmentation.
#' @param Smooth_range \code{"spatialKMeans"} pixel smooth range 
#' @param Virtual_segmentation set \code{"TRUE"} if you want to overide the automaitic segmentation
#' @param Virtual_segmentation_rankfile specify a region rank file contains region information for manualy region segmentation
#' @return None
#'
#' @examples
#' imaging_Spatial_Quant(threshold=0.05, ppm=5,Digestion_site="[G]",
#'                        missedCleavages=0:1,Fastadatabase="murine_matrisome.fasta",
#'                        adducts=c("M+H","M+NH4","M+Na"),PMF_analysis=TRUE,
#'                        Protein_feature_summary=TRUE,plot_cluster_image=TRUE,
#'                        Peptide_feature_summary=TRUE,plot_ion_image=FALSE,
#'                        parallel=3,spectra_segments_per_file=5,spatialKMeans=TRUE
#'                        )
#'
#' @export

imaging_Spatial_Quant<-function(
  #==============Choose the imzml raw data file(s) to process  make sure the fasta file in the same folder
  datafile=tk_choose.files(filter =  matrix(c( "imzml file", ".imzML",
                                                         "Text", ".txt", "All files", "*"),
                                                      3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis"),
  threshold=0.00, 
  ppm=5,
  #Fastadatabase="murine_matrisome.fasta",
  Quant_list="lipid candidates manual.csv",
  adducts=c("M-H","M+Cl"),
  cal.mz=T,
  mzlist_bypass=F,
  #==============TRUE if you want to plot protein PMF result
  PMF_analysis=TRUE,
  #==============TRUE if you want to generate protein summary in the Summary folder
  Protein_feature_summary=T,
  #==============TRUE if you want to generate protein cluster image in the Summary folder
  plot_cluster_image=T,
  plot_style="fleximaging",
  Peptide_feature_summary=T,
  #PLOT_PMF_Protein=FALSE,
  #==============TRUE if you want to plot peptide in the Ion images folder, make sure there's imzml file in the folder
  plot_ion_image=FALSE,
  #==============Set a number if you want a parallel processing
  parallel=detectCores()/2,
  #==============Set a number (1 to maximum pixels in the data file) if you want to dig more peaks in the raw data
  spectra_segments_per_file=5,
  spatialKMeans=F,
  Smooth_range=1,
  Virtual_segmentation=T,
  Virtual_segmentation_rankfile=tk_choose.files(default = "Z:/George skyline results/maldiimaging/Maldi_imaging - Copy/radius_rank.csv",caption  = "Choose Virtual segmentation rank info file"),
  Spectrum_feature_summary=T,
  Region_feature_summary=T,
  Region_feature_analysis=T,
  plot_each_metabolites=T,
  Cluster_level="High",
  Region_feature_analysis_bar_plot=T,
  norm_datafiles=T,
  norm_Type="Median",
  mzrang=NULL,
  BPPARAM=bpparam(),
  rotateimg=NULL,
  ...
){
  library(pacman)
  p_load(RColorBrewer,RCurl,bitops,magick,ggplot2,reticulate,dplyr,stringr,tcltk,data.table,doParallel,
         iterators,foreach,protViz,cleaver,MALDIquant,Biostrings,XVector,IRanges,Cardinal,Rdisop,
         ProtGenerics,S4Vectors,stats4,EBImage,
         BiocParallel,BiocGenerics,parallel,stats,graphics,grDevices,utils,datasets,methods)
  
  datafile<-gsub(".imzML", "", datafile)
  workdir<-base::dirname(datafile[1])
  setwd(workdir)
  closeAllConnections()
  #cl <- autoStopCluster(makeCluster(parallel))
  #setwd(workdir)
  cl <- makeCluster(parallel)
  message(paste(try(detectCores()), "Cores detected,",parallel, "threads will be used for computing"))
  
  message(paste(length(datafile), "files were selected and will be used for Searching"))
  
  Meta_feature_list<-Meta_feature_list_fun(workdir=workdir,
                                           database = Quant_list,
                                           adducts=adducts,
                                           cal.mz = cal.mz,
                                           bypass=mzlist_bypass)
  
  if(PMF_analysis){
    
    #Peptide_Summary_searchlist<-unique(Protein_feature_list[,c("Peptide","mz","adduct")])
    
    Peptide_Summary_file<-PMF_Cardinal_Datafilelist(datafile, 
                                                    Meta_feature_list,
                                                    SPECTRUM_for_average=spectra_segments_per_file,
                                                    threshold=threshold,
                                                    ppm=ppm,
                                                    spatialKMeans=spatialKMeans,
                                                    Virtual_segmentation=Virtual_segmentation,
                                                    Virtual_segmentation_rankfile=Virtual_segmentation_rankfile,
                                                    Smooth_range=Smooth_range,BPPARAM = BPPARAM,
                                                    PMFsearch = TRUE)
    
  }
  #Summarize the protein list across the datafiles
  if(Protein_feature_summary){
    for (i in 1:length(datafile)){
      datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
      currentdir<-paste0(datafile[i] ," ID")
      
      setwd(paste(currentdir,sep=""))
      
      standardcol=c("Peptide","Name","FA","Formula","mz","moleculeNames","mass","Intensity","adduct")
      Peptide_Summary_file<-fread("Peptide_Summary_file.csv")
      Peptide_Summary_file<-as.data.frame(Peptide_Summary_file)
      Peptide_Summary_file<-as.data.frame(Peptide_Summary_file[,c(as.character(intersect(colnames(Peptide_Summary_file),standardcol)))])
      Peptide_feature_list<-Peptide_Summary_file[Peptide_Summary_file$Intensity>=max(Peptide_Summary_file$Intensity)*threshold,]
      write.csv(Peptide_feature_list,"Peptide_feature_list.csv")
      uniques_intensity<-unique(Peptide_Summary_file[,c("mz","Intensity")])
      Meta_feature_list$Intensity<-0
      message("Iterating identification result")
      #Meta_feature_list$Intensity<-unlist(parLapply(cl=autoStopCluster(makeCluster(parallel)),Meta_feature_list$mz,intensity_sum_para,uniques_intensity))
      Meta_feature_list$Intensity<-unlist(bplapply(Meta_feature_list$mz,intensity_sum_para,uniques_intensity,BPPARAM=BPPARAM))
      write.csv(Meta_feature_list[Meta_feature_list$Intensity>0,],"Cluster.csv",row.names = F)
      
    }
    }
  #Summarize the peptide list across the datafiles
  if(Peptide_feature_summary){
    print("Peptide_feature_summary")
    if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
    Peptide_feature_summary_all_files_new(datafile,workdir,threshold = threshold)
    Protein_feature_summary_all_files_new(datafile,workdir,threshold = threshold)
    }
  #plot the peptide image across the datafiles
  if(plot_ion_image){
    print("plot_ion_image")
    if (dir.exists(paste(workdir,"/Ion images",sep=""))==FALSE){dir.create(paste(workdir,"/Ion images",sep=""))}
    setwd(paste(workdir,"/Ion images",sep=""))
    masslist<-read.csv(paste(workdir,"/Summary folder/Peptide_feature_summary_sl.csv",sep = ""))
    masslist<-unique(masslist[,c("Peptide","mz","adduct","moleculeNames")],by="moleculeNames")
    write.csv(masslist,paste(workdir,"/Ion images/Peptide_feature_summary_sl.csv",sep = ""))
    #simple_ion_image_cardinal(datafile,workdir=workdir,plotTolerance=ppm,interpolate =TRUE)
    simple_ion_image_cardinal(datafile=datafile,
                              #workdir=WorkingDir(),
                              Image_Type="",
                              plotTolerance=ppm ,
                              creat_new_file=TRUE,
                              color = "black",
                              smooth.image = "adaptive",
                              contrast.enhance = "histogram",
                              Neighbour=Smooth_range)
  }
  #Summarize the spectrum across the datafiles
  if(Spectrum_feature_summary){
    Spectrum_summary<-NULL
    for (i in 1:length(datafile)){
      datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
      currentdir<-paste0(datafile[i] ," ID")
      setwd(currentdir)
      name <-gsub(base::dirname(datafile[i]),"",datafile[i])
      message(paste("Spectrum_feature_summary",datafile[i]))
      
      if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
      #Feature_summary_all_files(datafile,workdir,threshold = threshold)
      match_pattern <- "Spectrum....csv"
      
      for (spectrum_file in dir()[str_detect(dir(), match_pattern)]){
        spectrum_file_table=fread(spectrum_file) 
        colnames(spectrum_file_table)=c("mz",paste(name,spectrum_file))
        if (is.null(Spectrum_summary)){
          Spectrum_summary=spectrum_file_table
        }else{
          Spectrum_summary=base::merge(Spectrum_summary,spectrum_file_table,by="mz",all=TRUE)
        }
      }
    }
    write.csv(Spectrum_summary,file = paste(workdir,"/Summary folder/Spectrum_summary.csv",sep=""),row.names = F)
  }
  #Summarize the features in regions across the datafiles
  if(Region_feature_summary){
    Spectrum_summary<-NULL
    for (i in 1:length(datafile)){
      datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".imzML", "", datafile[i]))
      currentdir<-paste0(datafile[i] ," ID")
      setwd(currentdir)
      name <-gsub(base::dirname(datafile[i]),"",datafile[i])
      message(paste("Region_feature_summary",datafile[i]))
      
      if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
      #Feature_summary_all_files(datafile,workdir,threshold = threshold)
      match_pattern <- "Peptide_region_file.csv"
      
      for (spectrum_file in dir()[str_detect(dir(), match_pattern)]){
        spectrum_file_table=fread(spectrum_file) 
        
        spectrum_file_table$ID=paste(spectrum_file_table$moleculeNames,spectrum_file_table$adduct,spectrum_file_table$mz,spectrum_file_table$Region,sep = "@")
        spectrum_file_table=spectrum_file_table[,c( "ID",	"Intensity")]
        colnames(spectrum_file_table)=c("ID",paste(name))
        if (is.null(Spectrum_summary)){
          Spectrum_summary=spectrum_file_table
        }else{
          Spectrum_summary=base::merge(Spectrum_summary,spectrum_file_table,by="ID",all=TRUE)
        }
      }
    }
    #Spectrum_summary=gsub("NA",0,Spectrum_summary)
    #for(col in names(Spectrum_summary)) set(Spectrum_summary, i=which(is.na(Spectrum_summary[[col]])), j=col, value=0)
    name <-gsub(base::dirname(datafile[1]),"",datafile)
    Spectrum_summary[is.na(Spectrum_summary)] <- 0
    colnames(Spectrum_summary)=c("ID",name)
    colnames(Spectrum_summary)=gsub("/","",colnames(Spectrum_summary))
    write.csv(Spectrum_summary,file = paste(workdir,"/Summary folder/Region_feature_summary.csv",sep=""),row.names = F)
  }
  #plot the features in regions across the datafiles
  if(Region_feature_analysis){
    
    setworkdir(paste(workdir,"/Summary folder/Region_feature_analysis/",sep=""))
    library(plotly)
    if (T){Spectrum_summary=read.csv(file = paste(workdir,"/Summary folder/Region_feature_summary.csv",sep=""),stringsAsFactors = F)}
    #clusterid=unlist(str_split(Spectrum_summary$ID,pattern = " "))
    #Spectrum_summary$ClusterID=str_split(Spectrum_summary$ID,pattern = " ")[[1]][1]
    radius_rank=read.csv(file =Virtual_segmentation_rankfile)
    if (norm_datafiles){
      
      #Spectrum_summary_norm=Spectrum_summary
      clusterID=as.character(Spectrum_summary$ID)
      clusterID=data.frame(str_split(clusterID," "),stringsAsFactors = F)
      clusterID=as.data.frame(t(clusterID))
      clusterID=as.character(clusterID[["V1"]])
      clusterID=data.frame(str_split(clusterID,"@"),stringsAsFactors = F)
      clusterID=as.data.frame(t(clusterID))
      clusterID=as.character(clusterID[["V1"]])
      
        par_norm_datacol_cluster<-function(col,Spectrum_summary,clusterID){
          Spectrum_summary_col_norm<-Spectrum_summary[,col]
          Spectrum_summary_col_norm_value<-Spectrum_summary[,col]
          Spectrum_summary_col<-Spectrum_summary[,col]
          for (uniqueclusterID in unique(clusterID)){
            
            Spectrum_summary_col_norm_value[clusterID==uniqueclusterID]= median.default(Spectrum_summary_col[clusterID==uniqueclusterID],na.rm = T)
          }
          Spectrum_summary_col_norm=Spectrum_summary_col/Spectrum_summary_col_norm_value
          #for (row in 1:length(Spectrum_summary_col)){
          #Spectrum_summary_col_norm[row]=(Spectrum_summary_col[row]/median.default(Spectrum_summary_col[clusterID==clusterID[row]],na.rm = T))
          #}  
          Spectrum_summary_col_norm
        }
        message("Region_feature_analysis")
        #Spectrum_summary_norm=unlist(parLapply(cl=autoStopCluster(makeCluster(detectCores())),which(colnames(Spectrum_summary)!="ID"),par_norm_datacol_cluster,Spectrum_summary,clusterID))
        Spectrum_summary_norm=unlist(bplapply(which(colnames(Spectrum_summary)!="ID"),par_norm_datacol_cluster,Spectrum_summary,clusterID,BPPARAM = BPPARAM))
        
        #a=as.list(Spectrum_summary[,which(colnames(Spectrum_summary)!="ID")])
        Spectrum_summary_norm=matrix(Spectrum_summary_norm,nrow = nrow(Spectrum_summary))
        Spectrum_summary[,which(colnames(Spectrum_summary)!="ID")]=Spectrum_summary_norm
      write.csv(Spectrum_summary,file = paste(workdir,"/Summary folder/Region_feature_summary_normalized.csv",sep=""),row.names = F)
      #Spectrum_summary=Spectrum_summary_norm
    
    }
    colnames_case=gsub("y . - root mean square","",colnames(Spectrum_summary))
    colnames_case=gsub("/","",colnames_case)
    colnames_case=gsub("X","",colnames_case)
    colnames_case=gsub("y.+","",colnames_case)
    colnames_case=gsub("ID","",colnames_case)
    
    #colnames_caset=gsub("y[:print:]+","",colnames_case)
    Spectrum_summary_tran=t(Spectrum_summary)
    
    Spectrum_summary_tran=data.frame(Spectrum_summary_tran)
    colnames(Spectrum_summary_tran)=Spectrum_summary$ID
    Spectrum_summary_tran[,"Class"]=colnames_case
    
    
    
    case_info=str_split(Spectrum_summary$ID,"@")[1:length(Spectrum_summary$ID)]
    case_info=as.data.frame(case_info)
    case_info=t(case_info)
    colnames(case_info)=c("ID","adducts","mz","Rank")
    case_info=merge(case_info,radius_rank[,c("Rank","Name")])
    
    colnames(case_info)=c("RegionRank","moleculeNames","adducts","mz","RegionName")
    case_info=merge(case_info,unique(Meta_feature_list[,c("moleculeNames","Name")]),all.x=T,all.y=F)
    case_info$ClusterID=case_info$Name
    if (Cluster_level=="High"){
    a=data.table(str_split(case_info$Name," ",simplify = T))
    case_info$ClusterID=a$V1
    }
    case_info=t(case_info)
    case_info=data.frame(case_info)
    colnames(case_info)=colnames(Spectrum_summary_tran)[1:(length(Spectrum_summary_tran)-1)]
    case_info$Class=""
    Spectrum_summary_tran=rbind(case_info,Spectrum_summary_tran)
    
    #case_info=t(case_info)
    #colnames(case_info)=c("ID","adducts","mz","Rank")
    #case_info<-data.frame(case_info)
    library(dplyr)
    library(plotly)
    library(stringr)
    moleculeNames=t(Spectrum_summary_tran["ID",])
    

    data=Spectrum_summary_tran[which(rownames(Spectrum_summary_tran)=="ID"):nrow(Spectrum_summary_tran),]
    data=Spectrum_summary_tran
    
    base::Sys.setenv("plotly_api_key"="FQ8kYs3JqmKLqKd0wGRv")
    base::Sys.setenv("plotly_username"="guoguodigital")
    base::Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1IjoiZ3VvZ3VvZGlnaXQiLCJhIjoiY2p1aHhheHM5MTBuYjQ0bnZzMzg0Mjd3aSJ9.XyqEayJi68xfGloNQQ28KA')
    
    plotly_for_region_with_ID<-function(i,data,moleculeNames){
      library(dplyr)
      library(plotly)
      library(stringr)
      if (!require("processx")) install.packages("processx")
      windows_filename<- function(stringX){
      stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
      stringX<-gsub("\"", "", stringX)
      return(stringX)}
      p <- plot_ly(data=data, x =data$Class[which(data$Class!="")] , y = data[which(data$Class!=""),moleculeNames[i]], name = moleculeNames[i], type = 'scatter', mode = 'markers') %>%
      add_trace(y =lowess(data$Class[which(data$Class!="")],as.numeric(as.character(data[[moleculeNames[i]]][which(data$Class!="")])),iter=3)$y, name = 'Moving average', mode = 'lines') %>%   
      layout(title = paste(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i]),
             xaxis = list(title = "Age"),
             yaxis = list (title = "Relative Conc."),
             showlegend = FALSE)
      
      
    #link=api_create(p, filename = paste(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"]))
    
    #plotly_IMAGE(p, format = "png", out_file = windows_filename(paste0(data["moleculeNames",i], data["adducts",i],"in",data["RegionName",i],".png")))
    ##png(windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".png")))
    #p
    #dev.off()
    htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".html")), selfcontained = F, libdir = "lib")
       
    #plotly::orca(p, windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".png")))
    windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".html"))
    }
    
    if(plot_each_metabolites){
      message("plot_each_metabolites")
    #p_names=parLapply(cl=autoStopCluster(makeCluster(4)),which(moleculeNames!=""),plotly_for_region_with_ID,data,moleculeNames)
    p_names=bplapply(which(moleculeNames!=""),plotly_for_region_with_ID,data,moleculeNames,BPPARAM=BPPARAM)
    zip("Result_lipid.zip", c(as.character(p_names), "lib"))}
    
    Spectrum_summary_tran_tran=as.data.frame(t(Spectrum_summary_tran),stringsAsFactors = FALSE)
    #Spectrum_summary_tran_tran=Spectrum_summary_tran_tran[which(rownames(Spectrum_summary_tran_tran)!="Class"),which(Spectrum_summary_tran_tran["Class",]!="")]
    Region_CLuster=paste(Spectrum_summary_tran_tran$ClusterID,Spectrum_summary_tran_tran$RegionName) 
    Spectrum_summary_tran_tran=Spectrum_summary_tran_tran[,which(Spectrum_summary_tran_tran["Class",]!="")]
      
    Spectrum_summary_tran_tran=Spectrum_summary_tran_tran[which(rownames(Spectrum_summary_tran_tran)!="Class"),]
    
    
    #Spectrum_summary_tran_tran=as.data.frame(as.numeric(as.character(Spectrum_summary_tran_tran)))

    Spectrum_summary_tran_tran <- mutate_all(Spectrum_summary_tran_tran, function(x) as.numeric(as.character(x)))
    #Spectrum_summary_tran_tran$ClusterID=Region_CLuster
    #Spectrum_summary_clustersum=aggregate(. ~ ClusterID, data=Spectrum_summary_tran_tran, sum)
    DT <- as.data.table(Spectrum_summary_tran_tran[])
    DT$ClusterID=Region_CLuster[which(Region_CLuster!=" ")]
    Spectrum_summary_tran_tran=DT[ , lapply(.SD, sum), by = "ClusterID"]
    rownames(Spectrum_summary_tran_tran)=Spectrum_summary_tran_tran$ClusterID
    data=as.data.frame(t(Spectrum_summary_tran_tran))
    colnames(data)=t(data["ClusterID",])
    data=data[which(rownames(data)!="ClusterID"),]
    data$Class=Spectrum_summary_tran$Class[which('&'(Spectrum_summary_tran$Class!="",Spectrum_summary_tran$Class!=" "))]
    
    Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
      
    plotly_for_region_with_ClusterID<-function(i,data,output_statice=F){
      library(dplyr)
      library(plotly)
      library(stringr)
      
      if (!require("processx")) install.packages("processx")
      windows_filename<- function(stringX){
        stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
        stringX<-gsub("\"", "", stringX)
        return(stringX)}
      moleculeNames=colnames(data)
      #fit3 <- lm(as.numeric(as.character( data[which(data$Class!=""),moleculeNames[i]]))~poly(as.numeric(as.character(data$Class[which(data$Class!="")])),3) )
      
      x <- as.factor(as.character(data$Class))
      y <- as.numeric(as.character(data[,moleculeNames[i]]))
      #plot(x,y)
      
      fit3 <- lm( y~poly(as.numeric(as.character(x)),3) )
      
      fitx=seq(min(as.numeric(data$Class)), max(as.numeric(data$Class)), length=1000)
      
      fitdata=as.data.frame(stats::predict(fit3, data.frame(x = fitx), interval="confidence"))
      newdata=data.frame(x=x,y=y)
      p <- plot_ly(data=newdata, x =data$Class , y = ~y, name = moleculeNames[i], type = 'scatter', mode = 'markers') %>%
        add_lines(x=unique(fitx),y=unique(fitdata$fit), name = 'Poly Fit', mode = 'lines',inherit=FALSE) %>%  
        add_lines(x=unique(fitx),y=unique(fitdata$lwr), name = 'Poly Fit lwr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>% 
        add_lines(x=unique(fitx),y=unique(fitdata$upr), name = 'Poly Fit upr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>% 
        #add_trace(lowess(x,y,iter=6)$y, name = 'Moving average', mode = 'lines') %>%   
        layout(title = moleculeNames[i],
               xaxis = list(title = "Age"),
               yaxis = list (title = "Conc."),
               showlegend = TRUE)
    
      #link=api_create(p, filename = paste(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"]))
      
      #plotly_IMAGE(p, format = "png", out_file = windows_filename(paste0(data["moleculeNames",i], data["adducts",i],"in",data["RegionName",i],".png")))
      ##png(windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".png")))
      #p
      #dev.off()
      if (output_statice){
       
      orca(p, file = windows_filename(paste0(moleculeNames[i],".png"))) 
        windows_filename(paste0(moleculeNames[i],".png"))
      }else{
      htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(moleculeNames[i],".html")), selfcontained = F, libdir = "lib")
      windows_filename(paste0(moleculeNames[i],".html"))  
        
      }
      
      
      #plotly::orca(p, windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".png")))
      #list(filename=windows_filename(paste0(moleculeNames[i],".html")), image= p %>% add_surface())
    }
    message("plotly_for_region_with_ClusterID")
    #p_names=parLapply(cl=autoStopCluster(makeCluster(1)),which(colnames(data)!="Class"),plotly_for_region_with_ClusterID,data,output_statice=T)
    p_names=bplapply(which(colnames(data)!="Class"),plotly_for_region_with_ClusterID,data,output_statice=T,BPPARAM = BPPARAM)
    zip("Result_clusterID.zip", c(as.character(p_names), "lib"))
    
    
    #Spectrum_summary_norm=Spectrum_summary
    #for (rownum in 1:nrow(Spectrum_summary)){
    # Spectrum_summary_norm[rownum,2:ncol(Spectrum_summary)]=Spectrum_summary[rownum,2:ncol(Spectrum_summary)]/max(Spectrum_summary[rownum,2:ncol(Spectrum_summary)])
    # #message(max(Spectrum_summary[rownum,2:ncol(Spectrum_summary)]))
    #}
    
    old_code<-function(){
    colnames_case=t(colnames_case)
    colnames(colnames_case)=colnames(Spectrum_summary_norm)
    Spectrum_summary_norm=rbind((colnames_case),Spectrum_summary_norm)
    org_row=nrow(Spectrum_summary_norm)
    org_col=ncol(Spectrum_summary_norm)
    Spectrum_summary_norm=as.data.frame(Spectrum_summary_norm)
    Spectrum_summary_norm_sum=NULL
    Spectrum_summary_norm_sum=data.frame(Spectrum_summary_norm$ID)
    colnames(Spectrum_summary_norm_sum)="ID"
    for( case in unique(paste(Spectrum_summary_norm[1,2:org_col]))){
      if(case!="ID"){
        #grep(case,colnames(Spectrum_summary_norm))
        Spectrum_summary_norm_sum[[case]]=0
        for (rownum in 2:org_row){
       Spectrum_summary_norm_sum[rownum,case]=mean(as.numeric(paste(Spectrum_summary_norm[rownum,grep(case,colnames(Spectrum_summary_norm))])),na.rm=TRUE)
        }
      }
      
    }
    
    rownames(Spectrum_summary_norm_sum)=Spectrum_summary_norm_sum[,1]
    Spectrum_summary_norm_sum=Spectrum_summary_norm_sum[1:nrow(Spectrum_summary_norm_sum),2:ncol(Spectrum_summary_norm_sum)]
    
    Spectrum_summary_norm_sum=t((Spectrum_summary_norm_sum))
    Spectrum_summary_norm_sum=as.data.frame(Spectrum_summary_norm_sum)
    Spectrum_summary_norm_sum$ID=as.numeric(rownames(Spectrum_summary_norm_sum))
    Spectrum_summary_norm_sum[,2:ncol(Spectrum_summary_norm_sum)]=round(Spectrum_summary_norm_sum[,2:ncol(Spectrum_summary_norm_sum)],digits = 3)
    #df=Spectrum_summary_norm_sum[,2:ncol(Spectrum_summary_norm_sum)]
    #class(df[2,])
    

     base::Sys.setenv("plotly_api_key"="FQ8kYs3JqmKLqKd0wGRv")
     base::Sys.setenv("plotly_username"="guoguodigital")
     base::Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1IjoiZ3VvZ3VvZGlnaXQiLCJhIjoiY2p1aHhheHM5MTBuYjQ0bnZzMzg0Mjd3aSJ9.XyqEayJi68xfGloNQQ28KA')
     library(plotly)
     case_info=str_split(colnames_case,"@")[2:length(colnames_case)]
     case_info=as.data.frame(case_info)
     case_info=t(case_info)
     colnames(case_info)=c("ID","adducts","mz","Rank")
     case_info<-data.frame(case_info)
     
     case_info=merge(case_info,radius_rank[,c("Rank","Name")])
     
      data <- data.frame(x, trace_0, trace_1, trace_2)
      colnames_case=colnames(Spectrum_summary_norm_sum)
      
      plotly_for_region_with_Clustered_ID<-function(i,Spectrum_summary_norm_sum,case_info){
        library(plotly)
        library(stringr)
        colnames_case=colnames(Spectrum_summary_norm_sum)
        windows_filename<- function(stringX){
          stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
          stringX<-gsub("\"", "", stringX)
          return(stringX)}
             p <- plot_ly(data=Spectrum_summary_norm_sum, x = ~ID, y = Spectrum_summary_norm_sum[[colnames_case[i]]], name = colnames_case[i], type = 'scatter', mode = 'markers') %>%
               add_trace(y =lowess(Spectrum_summary_norm_sum$ID,Spectrum_summary_norm_sum[[colnames_case[i]]],iter=6)$y, name = 'Moving average', mode = 'lines') %>%   
               layout(title = paste(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"]),
                   xaxis = list(title = "Age"),
                   yaxis = list (title = "Relative Conc."),
                   showlegend = FALSE)
       #link=api_create(p, filename = paste(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"]))
             
             plotly_IMAGE(p, format = "png", out_file = windows_filename(paste0(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"],".png")))
             htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"],".html")))
       #return(link[["api_urls"]][["plots"]])
      }
      listreturn=parallel::parSapply(cl=makeCluster(4),2:ncol(Spectrum_summary_norm_sum),plotly_for_region,Spectrum_summary_norm_sum,case_info)
       # Create a shareable link to your chart
       # Set up API credentials: https://plot.ly/r/getting-started
       chart_link = api_create(p, filename="line-mode1")
     chart_link}
    #Spectrum_summary=gsub("NA",0,Spectrum_summary)
    #for(col in names(Spectrum_summary)) set(Spectrum_summary, i=which(is.na(Spectrum_summary[[col]])), j=col, value=0)
    
    #Spectrum_summary[is.na(Spectrum_summary)] <- 0
    write.csv(Spectrum_summary,file = paste(workdir,"/Summary folder/Region_feature_summary_result.csv",sep=""),row.names = F)
  }
  #Plot bar chart of the features in regions across the datafiles
  if(Region_feature_analysis_bar_plot){
    options(scipen=999)
    setworkdir(paste(workdir,"/Summary folder/Region_feature_analysis/",sep=""))
    library(plotly)
    if (is.null(Spectrum_summary)){Spectrum_summary=read.csv(file = paste(workdir,"/Summary folder/Region_feature_summary.csv",sep=""))}
    radius_rank=read.csv(file =Virtual_segmentation_rankfile)
    if (norm_datafiles){

      Spectrum_summary_org=Spectrum_summary
      for (datacol in colnames(Spectrum_summary)[which(colnames(Spectrum_summary)!="ID")]){
        
        Spectrum_summary[[datacol]]=(Spectrum_summary[[datacol]]/max(Spectrum_summary[[datacol]]))
        
      }
      
      
    }
    colnames_case=gsub("y . - root mean square","",colnames(Spectrum_summary))
    colnames_case=gsub("RMS","",colnames(Spectrum_summary))
    colnames_case=gsub("istd","",colnames(Spectrum_summary))
    colnames_case=gsub("/","",colnames_case)
    colnames_case=gsub("ID","",colnames_case)
    Spectrum_summary_tran=t(Spectrum_summary)
    
    Spectrum_summary_tran=data.frame(Spectrum_summary_tran)
    colnames(Spectrum_summary_tran)=Spectrum_summary$ID
    Spectrum_summary_tran[,"Class"]=colnames_case
    
    
    
    case_info=str_split(Spectrum_summary$ID,"@")[1:length(Spectrum_summary$ID)]
    case_info=as.data.frame(case_info)
    
    case_info=as.data.frame(t(case_info))
    colnames(case_info)=c("ID","adducts","mz","Rank")
    case_info$orgID=Spectrum_summary$ID
    case_info=merge(case_info,radius_rank[,c("Rank","Name")])
    
    colnames(case_info)=c("RegionRank","moleculeNames","adducts","mz","orgID","RegionName")
    #case_info=merge(case_info,unique(Meta_feature_list[,c("moleculeNames","Name")]),all.x=T,all.y=F)
    case_info$Name=case_info$moleculeNames
    case_info$ClusterID=case_info$Name
    if (Cluster_level=="High"){
      a=data.table(str_split(case_info$Name," ",simplify = T))
      case_info$ClusterID=a$V1
    }
    rownames(case_info)=case_info$orgID
    case_info=t(case_info)
    case_info=as.data.frame(case_info)
    #colnames(case_info)=as.character(case_info["orgID",])
    case_info$Class=""
    Spectrum_summary_tran=rbind(case_info,Spectrum_summary_tran)
    Spectrum_summary_tran_plot=as.data.frame(t(Spectrum_summary_tran))
    Spectrum_summary_tran_plot$moleculeNames=paste(Spectrum_summary_tran_plot$moleculeNames,Spectrum_summary_tran_plot$adducts)
    Spectrum_summary_tran_plot["Class","moleculeNames"]=""
    plotlist=as.character(unique(Spectrum_summary_tran_plot$moleculeNames))[which(as.character(unique(Spectrum_summary_tran_plot$moleculeNames))!="")]
    p_names=sapply(plotlist,plotly_for_region_with_ClusterID_barchart,Spectrum_summary_tran_plot)
    zip("Result_clusterID.zip", c(as.character(p_names), "lib"))

    
  }
  #Plot cluster images across the datafiles
  if(plot_cluster_image){
    Protein_feature_list=fread(paste(workdir,"/Summary folder/Protein_feature_summary_sl.csv",sep=""))
    if (!is.null(rotateimg)){rotateimg=read.csv(rotateimg,stringsAsFactors = F)}
    for (i in 1:length(datafile)){
      rotate=rotateimg[rotateimg$filenames==datafile[i],"rotation"]
      rotate=as.numeric(rotate)
    imdata=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm,rotate = rotate,attach.only=F)
    
    
    currentdir<-paste0(datafile[i] ," ID")
    if (dir.exists(paste(currentdir,"/cluster Ion images",sep=""))==FALSE){dir.create(paste(currentdir,"/cluster Ion images",sep=""))}
    setwd(paste(currentdir,"/cluster Ion images",sep=""))
    
    lapply(unique(Protein_feature_list$Name),
           cluster_image_cardinal_allinone,
           imdata=imdata,
           SMPLIST=Protein_feature_list,
           ppm=ppm,ClusterID_colname="Name",
           componentID_colname="FA",
           plot_layout="line",
           Component_plot_threshold=1,
           plot_style=plot_style)
    }
      outputfolder=paste(workdir,"/Summary folder/cluster Ion images/",sep="")
      if (dir.exists(outputfolder)==FALSE){dir.create(outputfolder)}
      
      
      lapply(unique(Protein_feature_list$Name),
             cluster_image_cardinal_allinone,
             imdata=NULL,
             SMPLIST=Protein_feature_list,
             ppm=ppm,ClusterID_colname="Name",
             componentID_colname="FA",
             plot_layout="line",
             Component_plot_threshold=1,
             export_Header_table=T)
      
      currentdir<-paste0(datafile[6] ," ID")
      pngfiles=dir(paste(currentdir,"/cluster Ion images",sep=""),pattern = "_flex.png")
    for (pngfile in pngfiles){
      Pngcluster=NULL
    for (i in 1:length(datafile)){  
      currentdir<-paste0(datafile[i] ," ID")
      setwd(paste(currentdir,"/cluster Ion images",sep=""))
      pngimage=magick::image_read(pngfile)
      #pngimage<-magick::image_border(pngimage, "transparent", "30x0")
      #pngimage<-magick::image_annotate(pngimage,gsub(base::dirname(datafile[i]),"",datafile[i]),gravity = "northeast",color="white",size = 30,degrees = 90)
      #print(pngimage)
      #Pngcluster[[i]]=pngimage
      if( is.null(Pngcluster)){
       Pngcluster=pngimage
      }else{
        Pngcluster=magick::image_append(c(Pngcluster,pngimage),stack = T)
        
        }
      
    
    }
      #Pngclusterappend=magick::image_append(Pngcluster[[1:15]],stack = T)
     image_write(Pngcluster, paste0(outputfolder,"datafiles_",pngfile))
      files <- dir(tempfile())
    }
      Pngclusterkmean=NULL
      Pngclustervseg=NULL
      for (i in 1:length(datafile)){
       #imdata=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm)
        currentdir<-paste0(datafile[i] ," ID")
        setwd(paste(currentdir))
        name=gsub(base::dirname(datafile[i]),"",datafile[i])
        name=gsub("/","",name)
        pngimagekmean=magick::image_read("spatialKMeans_image.png")
        pngimagevseg=magick::image_read(paste0("Virtual_segmentation ",name,".png"))
        if( is.null(Pngclusterkmean)){
          Pngclusterkmean=pngimagekmean
          Pngclustervseg=pngimagevseg
        }else{
          Pngclusterkmean=magick::image_append(c(Pngclusterkmean,pngimagekmean),stack = T)
          Pngclustervseg=magick::image_append(c(Pngclustervseg,pngimagevseg),stack = T)
        }
      } 
      image_write(Pngclusterkmean, paste0(outputfolder,"datafiles_kmean.png"))
      image_write(Pngclustervseg, paste0(outputfolder,"datafiles_vseg.png"))
   }

  if(plot_cluster_image_grid){
    Protein_feature_list=fread(paste(workdir,"/Summary folder/Protein_feature_summary_sl.csv",sep=""))
    if (!is.null(rotateimg)){rotateimg=read.csv(rotateimg,stringsAsFactors = F)}
    imdata=list()
    combinedimdata=NULL
    register(SerialParam())
    for (i in 1:length(datafile)){
      
      rotate=rotateimg[rotateimg$filenames==datafile[i],"rotation"]
      rotate=as.numeric(rotate)
      imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm*2,rotate = rotate,attach.only=F,as="MSImagingExperiment",mzrange=mzrange)
      #imdata[[i]]@elementMetadata@coord=imdata[[i]]@elementMetadata@coord[,c("x","y")]
      if (i==1) {
        combinedimdata=imdata[[i]]
        }else{
        combinedimdata=cbind(combinedimdata,imdata[[i]]) 
        }
      imdata[[i]]=NULL
      
    }
    
    
    combinedimdata@elementMetadata@coord@listData[["z"]]=NULL
    
    imdata=combinedimdata
    
    outputfolder=paste(workdir,"/Summary folder/cluster Ion images/",sep="")
    
    if (dir.exists(outputfolder)==FALSE){dir.create(outputfolder)}
    setwd(outputfolder)
    
    lapply(unique(Protein_feature_list$Name),
           cluster_image_grid,
           imdata=NULL,
           SMPLIST=Protein_feature_list,
           ppm=ppm,ClusterID_colname="Name",
           componentID_colname="FA",
           plot_layout="line",
           Component_plot_threshold=1,
           export_Header_table=T)
    
    

    Pngclusterkmean=NULL
    Pngclustervseg=NULL
    for (i in 1:length(datafile)){
      #imdata=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm)
      currentdir<-paste0(datafile[i] ," ID")
      setwd(paste(currentdir))
      name=gsub(base::dirname(datafile[i]),"",datafile[i])
      name=gsub("/","",name)
      pngimagekmean=magick::image_read("spatialKMeans_image.png")
      pngimagevseg=magick::image_read(paste0("Virtual_segmentation ",name,".png"))
      if( is.null(Pngclusterkmean)){
        Pngclusterkmean=pngimagekmean
        Pngclustervseg=pngimagevseg
      }else{
        Pngclusterkmean=magick::image_append(c(Pngclusterkmean,pngimagekmean),stack = T)
        Pngclustervseg=magick::image_append(c(Pngclustervseg,pngimagevseg),stack = T)
      }
    } 
    image_write(Pngclusterkmean, paste0(outputfolder,"datafiles_kmean.png"))
    image_write(Pngclustervseg, paste0(outputfolder,"datafiles_vseg.png"))
  }
}











#I have wrapped cleave function without name each string vector with its AAsequence
Cleave_noname<-function(x, enzym = "trypsin", missedCleavages = 0,custom = NULL, unique = TRUE){
  returnlist_noname<-NULL
  returnlist<-cleave(x, enzym , missedCleavages,custom, unique )
  for (i in 1:length(x)){
    returnlist_noname[[i]]<- returnlist[[x[i]]]
  }
  return(returnlist_noname)  
}

#I have wrapped parentIonMass to process the list of peptides. we could simply use the default setting which will produce +1 (H+) precursor ion 
parentIonMasslist<-function(peplist,Index_of_protein_sequence){
AA<-rep(0,26)
  PIM<-NULL
  for (i in 1:length(peplist)){
    
    PIM[[Index_of_protein_sequence$desc[i]]] <- parentIonMass(peplist[[Index_of_protein_sequence$desc[i]]],fixmod=AA)}
  return(PIM)
}

Peptide_feature_summary_all_files<-function(listfile,workdir,threshold=0.1){
  Peptide_feature_summary<-NULL
  for (i in 1:length(listfile)){
    datafilename<-gsub(paste(workdir,"/",sep=""),"",gsub(".csv", "", listfile[i]))
    currentdir<-paste(workdir,"/ID/",datafilename,"/",sep = "")
    setwd(paste(currentdir,sep=""))
    
    Peptide_feature_list<-fread("Peptide_feature_list.csv",
                                select=c("Peptide","mz","Intensity","adduct","moleculeNames","freq"))
    Peptide_feature_list$sourcefile<-datafilename
    if (i==1){
      Peptide_feature_summary<-Peptide_feature_list
    }else{
      Peptide_feature_summary<-merge(Peptide_feature_summary,Peptide_feature_list,all = TRUE)
    }
    Peptide_feature_summary<-Peptide_feature_summary[,c("freq"):=NULL]
    setwd(paste(workdir,"/Summary folder/",sep = ""))
    write.csv(Peptide_feature_summary,"Peptide_feature_summary.csv")
    Peptide_feature_summary_sl<-Peptide_feature_summary[Peptide_feature_summary$Intensity>=max(Peptide_feature_summary$Intensity)*threshold,]
    write.csv(Peptide_feature_summary_sl,"Peptide_feature_summary_sl.csv")
    setwd(workdir)
  }

}

Peptide_feature_summary_all_files_new<-function(datafile,workdir,threshold=0.1){

  for (i in 1:length(datafile)){
    datafilename<-gsub(paste(base::dirname(datafile[i]),"/",sep=""),"",gsub(".imzML", "", datafile[i]))
    currentdir<-paste(datafile[i]," ID/",sep = "")
    setwd(paste(currentdir,sep=""))
    
    Peptide_feature_list<-fread("Peptide_feature_list.csv",
                                select=c("Peptide","mz","Intensity","adduct","moleculeNames"))
    Peptide_feature_list$sourcefile<-datafilename
    if (i==1){
      Peptide_feature_summary<-Peptide_feature_list
    }else{
      Peptide_feature_summary<-rbind.data.frame(Peptide_feature_summary,Peptide_feature_list)
    }}
    setwd(paste(base::dirname(datafile[i]),"/Summary folder/",sep = ""))
    write.csv(Peptide_feature_summary,"Peptide_feature_summary.csv",row.names = F)
    Peptide_feature_summary_sl<-Peptide_feature_summary[Peptide_feature_summary$Intensity>=max(Peptide_feature_summary$Intensity)*threshold,]
    write.csv(Peptide_feature_summary_sl,"Peptide_feature_summary_sl.csv",row.names = F)
    setwd(base::dirname(datafile[i]))
  
  
}

Protein_feature_summary_all_files_new<-function(datafile,workdir,threshold=0.1){
  
  for (i in 1:length(datafile)){
    datafilename<-gsub(paste(base::dirname(datafile[i]),"/",sep=""),"",gsub(".imzML", "", datafile[i]))
    currentdir<-paste(datafile[i]," ID/",sep = "")
    setwd(paste(currentdir,sep=""))
    
    Peptide_feature_list<-fread("Cluster.csv")
    Peptide_feature_list$sourcefile<-datafilename
    if (i==1){
      Peptide_feature_summary<-Peptide_feature_list
    }else{
      Peptide_feature_summary<-rbind.data.frame(Peptide_feature_summary,Peptide_feature_list)
    }}
  setwd(paste(base::dirname(datafile[i]),"/Summary folder/",sep = ""))
  write.csv(Peptide_feature_summary,"Protein_feature_summary.csv",row.names = F)
  Peptide_feature_summary_sl<-Peptide_feature_summary[Peptide_feature_summary$Intensity>=max(Peptide_feature_summary$Intensity)*threshold,]
  write.csv(Peptide_feature_summary_sl,"Protein_feature_summary_sl.csv",row.names = F)
  setwd(base::dirname(datafile[i]))
  
  
}

Feature_summary_all_files<-function(datafile,workdir,threshold=0.1){
  
  for (i in 1:length(datafile)){
    datafilename<-gsub(paste(base::dirname(datafile[i]),"/",sep=""),"",gsub(".imzML", "", datafile[i]))
    currentdir<-paste(datafile[i]," ID/",sep = "")
    setwd(paste(currentdir,sep=""))
    
    Peptide_feature_list<-fread("Peptide_region_file.csv")
    Peptide_feature_list$sourcefile<-datafilename
    if (i==1){
      Peptide_feature_summary<-Peptide_feature_list
    }else{
      Peptide_feature_summary<-rbind.data.frame(Peptide_feature_summary,Peptide_feature_list)
    }}
  setwd(paste(base::dirname(datafile[i]),"/Summary folder/",sep = ""))
  write.csv(Peptide_feature_summary,"Peptide_feature_summary.csv",row.names = F)
  Peptide_feature_summary_sl<-Peptide_feature_summary[Peptide_feature_summary$Intensity>=max(Peptide_feature_summary$Intensity)*threshold,]
  write.csv(Peptide_feature_summary_sl,"Peptide_feature_summary_sl.csv",row.names = F)
  setwd(base::dirname(datafile[i]))
  
  
}

Peptide_Feature_Summary<- function(peplist,pimlist,pimresultlist,Peptide_feature_list, mode="append"){
  tempdf<- NULL
  
for (Proteins in names(peplist)){
for (Peptides in  peplist[[Proteins]]){
  tempdf<- rbind(tempdf,data.frame("Protein"= Proteins,"Peptide"=Peptides,"mz" = pimlist[[Proteins]][peplist[[Proteins]]==Peptides], "Intensity" = pimresultlist[[Proteins]][peplist[[Proteins]]==Peptides]))
  }
}
  
if (mode=="append") {
  
  Peptide_feature_list<- rbind(Peptide_feature_list,tempdf)
  
}else{
  
  Peptide_feature_list<- tempdf
}
  return(Peptide_feature_list)
}

Peptide_Feature_Summary_fast<- function(peplist,pimlist,pimresultlist,Peptide_feature_list, mode="append"){
  print("Peptide Feature Summary")
  tempdf<- NULL
  
  for (Proteins in 1: length(names(peplist))){
    for (Peptides in  1: length(peplist[[Proteins]])){
      tempdf<- rbind(tempdf,data.frame("Protein"= names(peplist)[Proteins],"Peptide"=peplist[[Proteins]][Peptides],"mz" = pimlist[[Proteins]][Peptides], "Intensity" = pimresultlist[[Proteins]][Peptides]))
    }
    #print(Proteins)
  }
  
  if (mode=="append") {
    
    Peptide_feature_list<- rbind(Peptide_feature_list,tempdf)
    
  }else{
    
    Peptide_feature_list<- tempdf
  }
  return(Peptide_feature_list)
}

Peptide_Feature_Summary_data_table<- function(peplist,pimlist,pimresultlist,Peptide_feature_list, mode="append"){
  print("Protein Peptide Feature Summary")
  tempdf<- NULL
  
  for (Proteins in 1: length(names(peplist))){
    ptm <- proc.time()
    tempdf1<- NULL
    for (Peptides in  1: length(peplist[[Proteins]])){
      tempdf1<- rbind(tempdf1,c(names(peplist)[Proteins],peplist[[Proteins]][Peptides], pimlist[[Proteins]][Peptides], pimresultlist[[Proteins]][Peptides]))
    }
    tempdf<- rbind(tempdf,tempdf1)
    #print(paste(Proteins, "in", proc.time() - ptm))
  }
  
  if (mode=="append") {
    
    Peptide_feature_list<- rbind(Peptide_feature_list,tempdf)
    
  }else{
    
    Peptide_feature_list<- tempdf
  }
  Peptide_feature_list<-as.data.frame(Peptide_feature_list)
  colnames(Peptide_feature_list)<-c("Protein","Peptide",	"mz",	"Intensity")
  Peptide_feature_list$Intensity<- as.numeric(Peptide_feature_list$Intensity)
  Peptide_feature_list$mz<- as.numeric(Peptide_feature_list$mz)
  print("Protein Peptide Feature Summary finished")
  return(Peptide_feature_list)
  
}

Peptide_Feature_Summary_data_table_para<- function(Proteins,peplist,pimlist,pimresultlist){

    tempdf1<- NULL
    for (Peptides in  1: length(peplist[[Proteins]])){
      tempdf1<- rbind(tempdf1,c(names(peplist)[Proteins],peplist[[Proteins]][Peptides], pimlist[[Proteins]][Peptides], pimresultlist[[Proteins]][Peptides]))
    }

  tempdf1
}

Peptide_Summary_para<- function(Proteins,peplist){
  
  tempdf1<- NULL
  for (Peptides in  1: length(peplist[[Proteins]])){
    tempdf1<- rbind(tempdf1,c(names(peplist)[Proteins],peplist[[Proteins]][Peptides]))
  }

  tempdf1
}

Peptide_Feature_Summary_data_table_fast<- function(peplist,pimlist,pimresultlist,Peptide_feature_list, mode="append"){
  print("Peptide Feature Summary")
  tempdf<- NULL
  
  for (Proteins in 1: length(names(peplist))){
    ptm <- proc.time()
    tempdf1<- NULL
    for (Peptides in  1: length(peplist[[Proteins]])){
      tempdf1<- rbind(tempdf1,c(names(peplist)[Proteins],peplist[[Proteins]][Peptides], pimlist[[Proteins]][Peptides], pimresultlist[[Proteins]][Peptides]))
    }
    
    tempdf<- rbind(tempdf,tempdf1)
    #print(paste(Proteins, "in", proc.time() - ptm))
  }
  
  if (mode=="append") {
    
    Peptide_feature_list<- rbind(Peptide_feature_list,tempdf)
    
  }else{
    
    Peptide_feature_list<- tempdf
  }
  colnames(Peptide_feature_list)<-c("Protein","Peptide",	"mz",	"Intensity")
  return(Peptide_feature_list)
}

Peptide_Feature_Summary_array<- function(peplist,pimlist,pimresultlist,Peptide_feature_list=NULL, mode="append"){
  temparray<- NULL
  
  for (Proteins in 1: length(names(peplist))){
    for (Peptides in  1: length(peplist[[Proteins]])){
      temparray<- rbind(temparray,c(names(peplist)[Proteins],peplist[[Proteins]][Peptides],pimlist[[Proteins]][Peptides],  pimresultlist[[Proteins]][Peptides]))}
    #print(Proteins)
    }
  
  if (mode=="append") {
    
    Peptide_feature_list<- rbind(Peptide_feature_list,temparray)
    
  }else{
    
    Peptide_feature_list<- temparray
  }
  colnames(Peptide_feature_list)<-c("Protein","Peptide","mz","Intensity")
  return(Peptide_feature_list)
}

Peptide_Feature_UNIQUE<- function(Peptide_feature_list=NULL){
  temparray<- NULL
  
  for (Proteins in 1: length(names(peplist))){
    for (Peptides in  1: length(peplist[[Proteins]])){
      temparray<- rbind(temparray,c(names(peplist)[Proteins],peplist[[Proteins]][Peptides],pimlist[[Proteins]][Peptides],  pimresultlist[[Proteins]][Peptides]))}
    print(Proteins)}
  
  if (mode=="append") {
    
    Peptide_feature_list<- rbind(Peptide_feature_list,temparray)
    
  }else{
    
    Peptide_feature_list<- temparray
  }
  colnames(Peptide_feature_list)<-c("Protein","Peptide","mz","Intensity")
  return(Peptide_feature_list)
}

searchPMF<-function(pimlist,spectrumlist,ppm){
  pimresultlist<-pimlist
  print("Start PMF search")
  for (pim in 1:length(pimlist)){
    
    pimresultlist[[pim]]<-PMFsum(pimlist[[pim]],spectrumlist,ppm) 
  
  
  }
  
  return(pimresultlist)
}

searchPMF_para<-function(pimlist,spectrumlist,ppm,BPPARAM=bpparam()){
  pimresultlist<-pimlist
  message("Start PMF search")
    #pimresultlist<-parLapply(cl=cl,pimlist,PMFsum,spectrumlist,ppm)
    pimresultlist<-bplapply(pimlist,PMFsum,spectrumlist,ppm,BPPARAM = BPPARAM)
  
  return(pimresultlist)
}

searchPMF_data_frame<-function(pimlist,spectrumlist,ppm,BPPARAM = bpparam()){
  pimresultlist<-pimlist
  message("Start PMF search")
  #pimresultlist<-parLapply(cl=cl,pimlist,PMFsum,spectrumlist,ppm)
  pimresultlist<-bplapply(pimlist,PMFsum,spectrumlist,ppm,BPPARAM = BPPARAM)
  return(pimresultlist)
}

PlotPMFsig<-function(pimresultindex,spectrumlist,peplist,pimlist,pimresultlist, threshold=0.05){
  library(ggplot2)
  #library(ggplot)
  print("Ploting Sig PMF")
  pimresultsl<-pimresultindex[pimresultindex[,"Normalized Mean"]>threshold,]
  
  plotpeaklist<-spectrumlist
  plotpeaklist[,1]<-as.numeric(plotpeaklist[,1])
  plotpeaklist[,2]<-as.numeric(plotpeaklist[,2])
  plotpeaklist<-plotpeaklist[plotpeaklist[,"intensities"]>0,]
  
  sp<-ggplot(plotpeaklist, size=1 ,aes(x=m.z, y=intensities)) +  geom_line()
  for (rownamelist in rownames(pimresultsl)){
  df<- data.frame("mz" = pimlist[[rownamelist]],"mzend" = pimlist[[rownamelist]], "yintercept" = rep(0,length(pimresultlist[rownamelist])), "intensities" = pimresultlist[[rownamelist]],"peptide"= as.factor(peplist[[rownamelist]]))
 
  colnames(df)<-c("mz","mzend","yintercept","intensities","peptide")
  df<-df[df[,"intensities"]>0,]  
  png(paste(getwd(),"\\",windows_filename(str_split(rownamelist,"\\|")[[1]][2]),'.png',sep=""),width = 1980,height = 1080)
  tempsp<-sp+ggtitle(paste(rownamelist,"\nNormalized intensities:",pimresultsl[rownamelist,"Normalized Mean"])) +geom_segment(aes(x = as.numeric(mz), y = as.numeric(yintercept), xend = as.numeric(mzend), yend = as.numeric(intensities), colour =peptide),lineend = "round", data = df,size = 1)
  print(tempsp)
  dev.off()
  }
  }

Plot_PMF_all<-function(Protein_feature_list,peaklist,threshold=threshold,savename=""){
  message("Ploting PMF all in one spectrum")
  library(ggplot2)
  plotpeaklist<-as.data.frame(peaklist)
  plotpeaklist[,1]<-as.numeric(plotpeaklist[,1])
  plotpeaklist[,2]<-as.numeric(plotpeaklist[,2])
  plotpeaklist<-plotpeaklist[plotpeaklist[,2]>0,]
  Protein_feature_list<-Protein_feature_list[Protein_feature_list$Intensity>=max(Protein_feature_list$Intensity)*threshold,]
  Protein_feature_list<-unique(Protein_feature_list[,c( "mz" ,"Intensity","adduct","moleculeNames")])
  sp<-ggplot2::ggplot(plotpeaklist, size=1 ,aes(x=m.z, y=intensities,xend=m.z,yend=rep(0,length(plotpeaklist[,2])))) +  geom_segment()
    df<- data.frame("mz" = Protein_feature_list$mz,"mzend" = Protein_feature_list$mz, "yintercept" = rep(0,length(Protein_feature_list$mz)), "intensities" = Protein_feature_list$Intensity,"moleculeNames"=Protein_feature_list$moleculeNames)
    
    colnames(df)<-c("mz","mzend","yintercept","intensities","moleculeNames")
    df<-df[df[,"intensities"]>0,]  
    png(paste(getwd(),savename,"PMF spectrum match",'.png',sep=""),width = 1980,height = 1080)
    tempsp<-sp+ggtitle(paste("PMF spectrum","\nNormalized intensities:")) +geom_segment(aes(x = as.numeric(mz), y = as.numeric(yintercept), xend = as.numeric(mzend), yend = as.numeric(intensities), colour =moleculeNames),lineend = "round", data = df,size = 1,show.legend=FALSE)
    print(tempsp)
    dev.off()
  
}

PMFsum<-function(pimmzlist,spectrumlist,ppm){
  intensitysum<-0
  
  for (mz in 1:length(pimmzlist)){
    lowmz<-pimmzlist[mz]-pimmzlist[mz]*ppm/1000000
    highmz<-pimmzlist[mz]+pimmzlist[mz]*ppm/1000000
    intensitysum[mz]<-sum(spectrumlist[data.table::between(spectrumlist[,1], lowmz, highmz),2])
    }
  return(intensitysum)
}

PMFsum_para<-function(pimmzlist,spectrumlist,ppm){
 
  
    lowmz<-pimmzlist-pimmzlist*ppm/1000000
    highmz<-pimmzlist+pimmzlist*ppm/1000000
    #sum(spectrumlist[data.table::between(spectrumlist[,1], lowmz, highmz),2])
    return(sum(spectrumlist$intensities[data.table::between(spectrumlist[,"m.z"], lowmz, highmz)]))
}

PMFresultindex<-function(resultlist){
  pimresultindex<-resultlist
  maxmean<-0
  maxsum<-0
  results<-base::names(resultlist)
  
  for (result in results){
    pimresultindex[[result]]<-0
    resultlist[[result]]<-as.numeric(resultlist[[result]])
    pimresultindex[[result]]<-as.numeric(pimresultindex[[result]])
    
    pimresultindex[[result]]<-c(mean(resultlist[[result]]))
    pimresultindex[[result]]<-c(sum(resultlist[[result]]),pimresultindex[[result]])
  if (mean(resultlist[[result]])>maxmean){maxmean<-mean(resultlist[[result]])}
  if (sum(resultlist[[result]])>maxsum){maxsum<-sum(resultlist[[result]])}
    }
  
  for (result in results){
    pimresultindex[[result]]<-c(pimresultindex[[result]],sum(resultlist[[result]])/maxsum)
    pimresultindex[[result]]<-c(pimresultindex[[result]],mean(resultlist[[result]])/maxmean)
  }
  pimresultindex<-t(as.data.frame(pimresultindex))
  colnames(pimresultindex)<-c("Sum","mean","Normalized Sum","Normalized Mean")
  rownames(pimresultindex)<-base::names(resultlist)
  return(pimresultindex)
  }

peptideSearchX <- function (x, peptideSequence,pimIdx = parentIonMass(peptideSequence),peptideMassTolerancePPM = 5,framentIonMassToleranceDa = 0.01, FUN = .byIon){
     query.mass <- ((x$pepmass * x$charge)) - (1.007825 * (x$charge -
                                                              + 1))
     eps <- query.mass * peptideMassTolerancePPM * 1e-06
     lower <- findNN(query.mass - eps, pimIdx)
     upper <- findNN(query.mass + eps, pimIdx)
     rv <- lapply(peptideSequence[lower:upper], function(p) {
       psm(p, x, plot = FALSE, FUN = FUN)
       })
     rv.error <- sapply(rv, function(p) {
       sum(abs(p$mZ.Da.error) < framentIonMassToleranceDa)
       })
     idx.tophit <- which(rv.error == max(rv.error))[1]
     data.frame(mass_error = eps,
                  idxDiff = upper - lower,
                  charge = x$charge,
                  pepmass = query.mass,
                  peptideSequence = rv[[idx.tophit]]$sequence,
                  groundTrue.peptideSequence = x$peptideSequence,
                  ms2hit = (rv[[idx.tophit]]$sequence == x$peptideSequence), hit = (x$peptideSequence %in% peptideSequence[lower:upper]))
}

THERO_IONS_peptides<-function(peptides){
pim<-parentIonMass(peptides)
fi<-fragmentIon(peptides)
par(mfrow=c(3,1))
for (i in 1:length(peptides)){
   plot(0,0,
          xlab='m/Z',
         ylab='',
        xlim=range(c(fi[i][[1]]$b,fi[i][[1]]$y)),
        ylim=c(0,1),
        type='n',
        axes=FALSE,
        sub=paste( pim[i], "Da"));
  box()
  axis(1,fi[i][[1]]$b,round(fi[i][[1]]$b,2))
  pepSeq<-strsplit(peptides[i],"")
  axis(3,fi[i][[1]]$b,pepSeq[[1]])
  
     abline(v=fi[i][[1]]$b, col='red',lwd=2)
   abline(v=fi[i][[1]]$c, col='orange')
   abline(v=fi[i][[1]]$y, col='blue',lwd=2)
   abline(v=fi[i][[1]]$z, col='cyan')
   }

}

WorkingDir <-function(Mainfolder = tk_choose.dir(caption = "Select working directory")){
  Mainfolder <- Mainfolder
  setwd(Mainfolder)
  Mainfolder
}

Plot_Ion_image_Png<- function(WKdir, imagefile, png_filename, mz,adduct="M-H", Tolerance= 0.25, title="",Neighbour = 0,Creat_new_file=TRUE,color="black",interpolate =FALSE){
  #x11()

  
  #dev.copy(png,paste(dir,"\\",png_filename,'.png',sep=""))
  Tolerance<-round(Tolerance,digits = 5)
  pngfillewrite<-paste(WKdir,"\\",png_filename,'.png',sep="")
  #pngfillewrite<-nice_file_create(pngfillewrite)
  png(paste(WKdir,"\\","temp.png",sep=""), bg = "transparent")
  try(plotMsiSlice(imagefile ,mz , tolerance=Tolerance, legend=FALSE,colRamp=colorRamp(c("black", "blue", "green", "yellow", "red","#FF00FF","white")),interpolate =interpolate ),silent = TRUE)
  dev.off()
  try(pngfile<-image_read(paste(WKdir,"\\","temp.png",sep="")))
  res <- try(pngfile<-image_trim(pngfile),silent = TRUE)
  if (class(res) != "try-error"){

    kern=Kernel_convolution(Neighbour)
    if (title=="D-mode"){
      pngfile<-image_border(pngfile, "transparent", "0x30")
      pngfile<-image_annotate(pngfile,paste(png_filename,"D-mode","\n",adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 9)
      pngfile<-image_trim(pngfile)
      res2 <- try(pngfileoriginal<-image_read(pngfillewrite),silent = TRUE)
      if ('&'(class(res2) != "try-error",Creat_new_file==FALSE)){
        pngfile<-image_append(c(pngfileoriginal,pngfile))
        image_write(pngfile,pngfillewrite)}
      if ('&'(class(res2) != "try-error",Creat_new_file==TRUE)){    
        image_write(pngfile,pngfillewrite)}
      if (class(res2) == "try-error"){    
        image_write(pngfile,pngfillewrite)}    
    }else{
      pngfile<-image_convolve(pngfile, kern)
      pngfile<-image_border(pngfile, "transparent", "20x20")
      pngfile<-image_annotate(pngfile,paste(png_filename,"\n", adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 10,color = color)
      pngfile<-image_trim(pngfile)
      
      image_write(pngfile,pngfillewrite)
    }
  }
  file.remove(paste(WKdir,"\\","temp.png",sep=""))
}

Kernel_convolution<-function(Neighbour){
  Neighbour<-1+(Neighbour*2)
  if (Neighbour==1){kern = matrix(1, ncol = Neighbour, nrow = Neighbour)}else{
    kern = matrix(1/(Neighbour*Neighbour-1), ncol = Neighbour, nrow = Neighbour)  
    kern[(Neighbour+1)/2,(Neighbour+1)/2]=0}
  kern
}

simple_ion_image<-function(workdir=WorkingDir(),
                           Image_Type="",
                           plotTolerance=5 ,
                           creat_new_file=TRUE,
                           Denoising=0,
                           color = "black",
                           interpolate =FALSE,
                           ...){
  
  datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
  listfile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.csv$")
  masslist<-NULL
  for (z in 1:length(datafile)){
    
    datafile[z]<-gsub(".imzML", "", datafile[z])  
    options(warn=-1)
    MALDI_IMAGE <- quiet(importImzMl(paste(file.path(datafile[z]),".imzML",sep="")))
    options(warn=0)
    print(paste("Plot Ion Images from",datafile[z]))
    if (dir.exists(datafile[z])==FALSE){dir.create(datafile[z])}
    
    for (i in 1:length(listfile)){
      foldername<-str_remove(listfile[i],paste(workdir,"/",sep=""))
      foldername<-str_remove(foldername,".csv")
      outputDir<-paste(datafile[z],"/",foldername,sep="")
      if (dir.exists(outputDir)==FALSE){dir.create(outputDir)}
      masslist<-read.csv(listfile[i],header = TRUE)
      for (j in 1:length(masslist$moleculeNames)){
        if (is.na(masslist$mz[j])==FALSE){
          for (k in 1:length(str_split(masslist$moleculeNames[j],", ")[[1]])){
            Plot_Ion_image_Png(WKdir=outputDir, 
                               imagefile=MALDI_IMAGE, 
                               png_filename=windows_filename(str_split(masslist$moleculeNames[j],", ")[[1]][k]),
                               mz=masslist$mz[j],
                               adduct=masslist$adduct[j], 
                               Tolerance= plotTolerance*masslist$mz[j]/1000000, 
                               title="",Creat_new_file=TRUE,Neighbour = Denoising,
                               color = color,interpolate =interpolate ) 
          }}}
    }
  }
}



quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

cleave_para<-function(protein_sequence){
  peplistpara<-cleave(as.character(protein_sequence),custom=Digestion_site, missedCleavages=missedCleavages)
}

autoStopCluster <- function(cl) {
  stopifnot(inherits(cl, "cluster"))
  env <- new.env()
  env$cluster <- cl
  attr(cl, "gcMe") <- env
  reg.finalizer(env, function(e) {
    #message("Finalizing cluster ...")
    #message(capture.output(print(e$cluster)))
    try(parallel::stopCluster(e$cluster), silent = FALSE)
    #message("Finalizing cluster ... done")
  })
  cl
}




PMF_analysis_fun<-function(Peptide_Summary_searchlist,peaklist,ppm,BPPARAM =bpparam(),mzrange){
  if (missing(mzrange)){mzrange=c(min(peaklist$m.z),max(peaklist$m.z))}
  #Peptide_Summary_searchlist$Intensity<-parLapply(cl=cl,Peptide_Summary_searchlist$mz,PMFsum_para,peaklist,ppm)
  #Peptide_Summary_searchlist=Peptide_Summary_searchlist[`&`(Peptide_Summary_searchlist$mz>=mzrange[1],Peptide_Summary_searchlist$mz<=mzrange[2]),]
  Peptide_Summary_searchlist=as.data.frame(Peptide_Summary_searchlist)
  Peptide_Summary_searchlist[,"Intensity"]<-0
  message("PMF 1st search...")
  Peptide_Summary_searchlist$Intensity<-bplapply(Peptide_Summary_searchlist$mz,PMFsum_para,peaklist,ppm,BPPARAM = BPPARAM)
  #Peptide_Summary_searchlistIntensity<-lapply(Peptide_Summary_searchlistmz,PMFsum_para,peaklist,ppm)
  Peptide_Summary_searchlist <- as.data.frame(sapply(Peptide_Summary_searchlist,unlist),stringsAsFactors=FALSE)
  Peptide_Summary_searchlist$Intensity<-as.numeric(Peptide_Summary_searchlist$Intensity)
  Peptide_Summary_searchlist$mz<-as.numeric(Peptide_Summary_searchlist$mz)
  Peptide_Summary_searchlist
}

simple_ion_image_cardinal<-function(datafile=tk_choose.files(filter =  matrix(c( "imzml file", ".imzML",
                                                                                           "Text", ".txt", "All files", "*"),
                                                                                        3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis"),
                                    #workdir=WorkingDir(),
                                    Image_Type="",
                                    plotTolerance=5 ,
                                    creat_new_file=TRUE,
                                    color = "black",
                                    smooth.image = "gaussian",
                                    contrast.enhance = "none",
                                    ...){
  datafile<-gsub(".imzML", "", datafile)
  workdir=base::dirname(datafile[1])
  if (dir.exists(paste0(workdir ,"/Ion images/"))==FALSE){dir.create(paste0(workdir ,"/Ion images/"))}
  setwd(paste0(workdir ,"/Ion images"))
  listfile <- list.files(path=paste0(workdir ,"/Ion images"),full.names = TRUE,pattern = "\\.csv$")
  masslist<-NULL
  for (z in 1:length(datafile)){
    name <-gsub(paste0(base::dirname(datafile[z]),"/"),"",datafile[z])
    folder<-base::dirname(datafile[z])
    MALDI_IMAGE <- readImzML(name, folder, attach.only=F,as = "MSImageSet")
    print(paste("Plot Ion Images from",name))
    if (dir.exists(paste0(workdir ,"/Ion images/",name))==FALSE){dir.create(paste0(workdir ,"/Ion images/",name))}
    
    for (i in 1:length(listfile)){
      listfilename<-gsub(base::dirname(listfile[i]),"",listfile[i])
      outputDir<-paste0(workdir ,"/Ion images/",name,"/",listfilename)
      if (dir.exists(outputDir)==FALSE){dir.create(outputDir)}
      masslist<-read.csv(listfile[i],header = TRUE)
      masslist<-masslist[!is.na(masslist$mz),]
      for (j in 1:length(masslist$moleculeNames)){
        if (is.na(masslist$mz[j])==FALSE){
          for (k in 1:length(str_split(masslist$moleculeNames[j],", ")[[1]])){
            Plot_Ion_image_Png_Cardinal(WKdir=outputDir, 
                                        imagefile=MALDI_IMAGE, 
                                        png_filename=windows_filename(str_split(masslist$moleculeNames[j],", ")[[1]][k]),
                                        mz=masslist$mz[j],
                                        adduct=masslist$adduct[j], 
                                        Tolerance= plotTolerance*masslist$mz[j]/1000000, 
                                        title="",
                                        Creat_new_file=TRUE,
                                        color = color,
                                        smooth.image = smooth.image,
                                        contrast.enhance = contrast.enhance) 
          }}}
    }
  }
}

Plot_Ion_image_Png_Cardinal<- function(WKdir, imagefile, png_filename, mz,adduct="M-H", Tolerance= 5, title="",smooth.image = "adaptive",Creat_new_file=TRUE,color="black",contrast.enhance = "suppression"){
  #x11()
  #dev.copy(png,paste(dir,"\\",png_filename,'.png',sep=""))
  Tolerance<-round(Tolerance,digits = 5)
  pngfillewrite<-paste(WKdir,"\\",png_filename,'.png',sep="")
  #pngfillewrite<-nice_file_create(pngfillewrite)
  tempfilename=tempfile(pattern = "file", tmpdir = tempdir())
  png(paste(tempfilename,".png",sep=""), bg = "transparent")  
  par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(2, 0, 0, 0),
                          bty="n",pty="s",xaxt="n",
                          yaxt="n",
                          no.readonly = TRUE,ann=FALSE)
  try(Cardinal::image(imagefile ,mz=mz , 
                      plusminus=Tolerance,
                      contrast.enhance = contrast.enhance,
                      smooth.image = smooth.image ,
                      col.regions=intensity.colors_customize(),
                      asp = 1, 
                      add=F),
                      silent = TRUE
                      
                      )
  dev.off()
  try(pngfile<-image_read(paste(tempfilename,".png",sep="")))
  res <- try(pngfile<-image_trim(pngfile),silent = TRUE)
  if (class(res) != "try-error"){
    if (title=="D-mode"){
      pngfile<-image_border(pngfile, "transparent", "0x30")
      pngfile<-image_annotate(pngfile,paste(png_filename,"D-mode","\n",adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 9)
      pngfile<-image_trim(pngfile)
      res2 <- try(pngfileoriginal<-image_read(pngfillewrite),silent = TRUE)
      if ('&'(class(res2) != "try-error",Creat_new_file==FALSE)){
        pngfile<-image_append(c(pngfileoriginal,pngfile))
        image_write(pngfile,pngfillewrite)}
      if ('&'(class(res2) != "try-error",Creat_new_file==TRUE)){    
        image_write(pngfile,pngfillewrite)}
      if (class(res2) == "try-error"){    
        image_write(pngfile,pngfillewrite)}    
    }else{
      pngfile<-image_border(pngfile, "transparent", "30x30")
      pngfile<-image_annotate(pngfile,paste("\n",png_filename,"\n", adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 10,color = color)
      pngfile<-image_trim(pngfile)
      
      image_write(pngfile,pngfillewrite)
    }
  }
  file.remove(paste(tempfilename,".png",sep=""))
}

PMF_Cardinal_Datafilelist<-function(datafile,Peptide_Summary_searchlist, 
                                    SPECTRUM_for_average=5,threshold=0.1,
                                    ppm,spatialKMeans=FALSE,
                                    Smooth_range=1,
                                    colorstyle="Set1",
                                    Virtual_segmentation=FALSE,
                                    Virtual_segmentation_rankfile="Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy\\radius_rank.csv",
                                    PMFsearch=TRUE,
                                    rotate=NULL,
                                    matching_validation=T,
                                    BPPARAM=bpparam(),
                                    Bypass_generate_spectrum=F,
                                    scoring=T,
                                    score_method="balanced-SQRT",
                                    Rank=3,
                                    Protein_feature_list=NULL,
                                    Decoy_mode="",
                                    Decoy_search=T,
                                    ...){
  library(data.table)
  library(Cardinal)
  library(RColorBrewer)
  library(stringr)
  if (!is.null(rotate)){
    message("Found rotation info")
    #rotatedegrees=rotate[rotate$filenames==datafile,"rotation"]
    rotatedegrees=sapply(datafile,function(x,df){
      degree=df[df$filenames==x,"rotation"]
      if (length(degree)==0) {
        message("Missing rotation data please check the rotation configuration file: ",x)
        degree=0
      }
      degree
    },rotate)
    rotate=unlist(rotatedegrees)
  }else{rotate=rep(0,length(datafile))}
  #mycol <- color.map(map =c("black", "blue", "green", "yellow", "red","#FF00FF","white"), n = 100)
  #mycol <- colorRampPalette(c("black", "blue", "green", "yellow", "red","#FF00FF","white"))
  #mycol <- gradient.colors(100, start="white", end="blue")
  
  #imdata <- readImzML(name, folder, attach.only=FALSE,as="MSImageSet",resolution=10, units="ppm")
  Peptide_Summary_file<-Peptide_Summary_searchlist  
  for (z in 1:length(datafile)){

  Peptide_Summary_file$Intensity<-rep(0,nrow(Peptide_Summary_file))
  
  imdata <- Load_Cardinal_imaging(datafile[z],preprocessing = F,attach.only = T,resolution = 200,rotate = rotate[z],as="MSImageSet",BPPARAM = BPPARAM)
  name <-gsub(base::dirname(datafile[z]),"",datafile[z])
  folder<-base::dirname(datafile[z])
  coordata=imdata@pixelData@data
  if (dir.exists(paste0(datafile[z] ," ID"))==FALSE){dir.create(paste0(datafile[z] ," ID"))}
  setwd(paste0(datafile[z] ," ID")) 
    
    
  if (Bypass_generate_spectrum==F){
  message("Segmentation in progress...")
     #cl=autoStopCluster(makeCluster(6))

  
  if (spatialKMeans){
    set.seed(1)
    message(paste0("spatialKMeans computing for ",name))
  skm <- spatialKMeans(imdata, r=Smooth_range, k=SPECTRUM_for_average, method="adaptive")
  message(paste0("spatialKMeans finished for ",name))
  #skm@resultData[["new"]]=skm@resultData[["r = 1, k = 5"]]
  png(paste(getwd(),"\\","spatialKMeans_image_plot",'.png',sep=""),width = 1024,height = 720)
  #plot(skm, col=c("pink", "blue", "red","orange","navyblue"), type=c('p','h'), key=FALSE)
  #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1+ceiling(SPECTRUM_for_average/2), 2),
  par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
      bty="n",pty="s",xaxt="n",
      yaxt="n",
      no.readonly = TRUE,ann=FALSE)
  Cardinal::image(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), key=FALSE, ann=FALSE,axes=FALSE)
  legend("topright", legend=1:SPECTRUM_for_average, fill=brewer.pal(SPECTRUM_for_average,colorstyle), col=brewer.pal(SPECTRUM_for_average,"Paired"), bg="transparent",xpd=TRUE,cex = 1)
  
  #Cardinal::plot(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), type=c('p','h'), key=FALSE)
  #Cardinal::plot(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), type=c('p','h'), key=FALSE,mode="centers")
  #Cardinal::plot(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), type=c('p','h'), key=FALSE,mode="betweenss")
  Cardinal::plot(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), type=c('p','h'), key=FALSE,mode="withinss")
  legend("topright", legend=1:SPECTRUM_for_average, fill=brewer.pal(SPECTRUM_for_average,colorstyle), col=brewer.pal(SPECTRUM_for_average,"Paired"), bg="transparent",xpd=TRUE,cex = 1)
  dev.off()
  withinss=skm@resultData[[1]][["withinss"]]
  centers=skm@resultData[[1]][["centers"]]
  cluster=as.data.frame(skm@resultData[[1]][["cluster"]])
  cluster$Coor=rownames(cluster)
  write.csv(withinss,paste("spatialKMeans_RESULT","withinss.csv"),row.names = F)
  write.csv(centers,paste("spatialKMeans_RESULT","centers.csv"),row.names = F)
  write.csv(cluster,paste("spatialKMeans_RESULT","cluster.csv"),row.names = F)
  
  
 #cluster=skm@resultData[["r = 1, k = 5"]][["cluster"]]

  
  y=skm@resultData[[paste0("r = ",Smooth_range,", k = ", SPECTRUM_for_average)]][["cluster"]]
  y=as.data.frame(y)
  rownames(y)=1:nrow(y)
  y$pixel=1:nrow(y)
  regions=unique(y[,1])
  x=NULL
  for (i in 1:length(regions)){
    listname=as.character(regions[i])
  x[[listname]]<-y[y[, 1] == regions[i], "pixel"]
  }
  
  
  
  
  }else if(Virtual_segmentation){
  radius_rank=read.csv(file = Virtual_segmentation_rankfile)
  radius_rank=radius_rank[order(radius_rank$Rank),]
  coordist_para=function(i,coordata){
    
  coordist_para=NULL
  for( j in 1  :  nrow(coordata)){
  coordist_para[j]=sqrt((coordata$x[i]-coordata$x[j])^2+(coordata$y[i]-coordata$y[j])^2)
  }
  coordist_para
  }
  
  
  #coordistmatrix<-matrix(nrow=nrow(coordata),ncol=nrow(coordata))
  #coordistmatrix<-NULL
  #cl=autoStopCluster(cl)
  #coordistmatrix=parLapply(cl=cl,1:  nrow(coordata),coordist_para,coordata)
  coordistmatrix=bplapply(1:  nrow(coordata),coordist_para,coordata,BPPARAM = BPPARAM)
  coordistmatrix=matrix(unlist(coordistmatrix),nrow = nrow(coordata),ncol = nrow(coordata))
  coordistmatrix=as.data.table(coordistmatrix)
  coordistmatrix$sum=0
  #coordistmatrix$sum=base::unlist(parLapply(cl=cl,1:nrow(coordata),function(j,coordistmatrix,coordata){coordistmatrix$sum[j]=sum(coordistmatrix[j,1:nrow(coordata)])},coordistmatrix,coordata))
  coordistmatrix$sum=base::unlist(bplapply(1:nrow(coordata),function(j,coordistmatrix,coordata){coordistmatrix$sum[j]=sum(coordistmatrix[j,1:nrow(coordata)])},coordistmatrix,coordata,BPPARAM = BPPARAM))
  
  coorrange=max(coordistmatrix$sum)-min(coordistmatrix$sum)
  

  
  findedge<-function(coordata){
    #for (i in 1: nrow(coordata)){
      uniquex=unique(coordata$x)
      uniquey=unique(coordata$y)
      coordata$edge=FALSE
    for (x in uniquex){
      
      min=min(coordata[coordata$x==x,"y"])
      max=max(coordata[coordata$x==x,"y"])
      coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
      coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
    }
      
    for (y in uniquey){
      min=min(coordata[coordata$y==y,"x"])
      max=max(coordata[coordata$y==y,"x"])
      coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
      coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
    }
      coordata
  }
  
  coordata=findedge(coordata)
  
  
  
  #paste0("v",colnames(coordistmatrix[coordistmatrix$sum==min(coordistmatrix$sum),]))
  
  #center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),])
  
  
  #plot(rownames(coordistmatrix),coordistmatrix$sum)
  
  #write.csv(radius_rank,file = "radius_rank.csv", row.names = F)
  

  
  rank_pixel<-function(coordata,coordistmatrix){
    #coordata[coordata$edge==TRUE,]=coordata[coordata$edge==TRUE,]
    shape_center=coordata[coordistmatrix$sum==min(coordistmatrix$sum),]
    center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),1:nrow(coordata)])
    p_load(useful)
    From <- shape_center[rep(seq_len(nrow(shape_center)), each=nrow(coordata)),1:2]
    To <- coordata[,1:2]
    df=To-From
    center_edge_angle=cbind(coordata[,1:2],cart2pol(df$x, df$y, degrees = F),edge=coordata[,"edge"])
    center_edge_angle_sdge=center_edge_angle[center_edge_angle$edge==TRUE,]
     coordata$rank=0 
     coordata$pattern=""       

     for (i in 1: (nrow(coordata))){
       

       #From <- coordata[i,][rep(seq_len(nrow(coordata[i,])), each=nrow(coordata[coordata$edge==TRUE,])),1:2]
       #To <- coordata[coordata$edge==TRUE,][,1:2]

       if (coordata$edge[i]!=TRUE){      
         df=coordata[i,1:2]-shape_center[,1:2]
       point_center_angle=cbind(coordata[i,1:2],cart2pol(df$x, df$y, degrees = F))
       pointedge=center_edge_angle_sdge[which(abs(center_edge_angle_sdge$theta-point_center_angle$theta)==min(abs(center_edge_angle_sdge$theta-point_center_angle$theta))),]
       #message(pointedge)

         pointedge=pointedge[which.min(pointedge$r),]
         to_edge=coordistmatrix[[i]]['&'(coordata$x==pointedge$x,coordata$y==pointedge$y)]
      }else{to_edge=0}
       
       
       to_center=center_dist[i]
       total=to_edge+to_center
       max(radius_rank$Radius_U)
       norm_center_dist=to_center/total*max(radius_rank$Radius_U)
       coordata$rank[i]=as.character(radius_rank$Rank['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
       coordata$pattern[i]=as.character(radius_rank$Name['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
     }
     coordata
  }
  
  coordata=rank_pixel(coordata,coordistmatrix)
 
  
  x=NULL
  
  for (rank in coordata$rank){
    
    x[[rank]]=which(coordata$rank==rank)
    
  }
  write.csv(coordata,"coordata.csv",row.names = F)
  region_pattern <- factor(coordata$rank,levels=unique(coordata$rank), labels=unique(coordata$pattern))
  set.seed(1)
  skm <- spatialKMeans(imdata, r=Smooth_range, k=length(unique(coordata$rank)), method="adaptive")
  png(paste(getwd(),"\\","spatialKMeans_image",'.png',sep=""),width = 1024,height = 1024)
  #plot(skm, col=c("pink", "blue", "red","orange","navyblue"), type=c('p','h'), key=FALSE)
  #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1+ceiling(SPECTRUM_for_average/2), 2),
  par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 1),
      bty="n",pty="s",xaxt="n",
      yaxt="n",
      no.readonly = TRUE,ann=FALSE)
  Cardinal::image(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), key=FALSE, ann=FALSE,axes=FALSE)
  
  legend("topright", legend=1:SPECTRUM_for_average, fill=brewer.pal(SPECTRUM_for_average,colorstyle), col=brewer.pal(SPECTRUM_for_average,"Paired"), bg="transparent",xpd=TRUE,cex = 1)
  
  #Cardinal::plot(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), type=c('p','h'), key=FALSE)
  #Cardinal::plot(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), type=c('p','h'), key=FALSE,mode="centers")
  #Cardinal::plot(skm, col=brewer.pal(SPECTRUM_for_average,colorstyle), type=c('p','h'), key=FALSE,mode="betweenss")
  dev.off()
  
  for (i in 1:nrow(coordata)){
  skm@resultData[[paste0("r = ",Smooth_range, ", k = ",length(unique(coordata$rank)))]][["cluster"]][[rownames(coordata)[i]]]=factor(coordata[i,"rank"],levels=unique(coordata$rank))
    
  }
  
  #Cardinal::image(imdata, col=brewer.pal(SPECTRUM_for_average,colorstyle), key=FALSE, ann=FALSE,axes=FALSE,groups =pattern)
  
  #msset <- generateImage(region_pattern, coord=coordata[,1:2],
  #                        range=c(1000, 5000), centers=c(2000, 3000, 4000),
  #                        resolution=100, step=3.3, as="MSImageSet")
  #msset@pixelData@data[["sample"]]=region_pattern
  png(paste(getwd(),"\\","Virtual_segmentation",gsub("/"," ",name),'.png',sep=""),width = 1024,height = 1024)
  #plot(skm, col=c("pink", "blue", "red","orange","navyblue"), type=c('p','h'), key=FALSE)
  #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1+ceiling(SPECTRUM_for_average/2), 2),
  par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 1),
      bty="n",pty="s",xaxt="n",
      yaxt="n",
      no.readonly = TRUE,ann=FALSE)
  
  Cardinal::image(skm, col=brewer.pal(length(unique(coordata$rank)),colorstyle), key=FALSE, ann=FALSE,axes=FALSE)
  
  #Cardinal::image(msset,feature=1:3, col=brewer.pal(length(unique(coordata$pattern)),colorstyle),strip=T,superpose=T, key=FALSE, ann=FALSE,axes=FALSE,groups =region_pattern)
  legend("topright", legend=paste(radius_rank$Rank,radius_rank$Name), fill=brewer.pal(length(unique(coordata$pattern)),colorstyle), col=brewer.pal(length(unique(coordata$pattern)),"Paired"), bg="transparent",xpd=TRUE,cex = 1)

  dev.off()
  
  
  
  }else{

  x=1:length(pixels(imdata))
  x=split(x, sort(x%%1)) 
    
    
  } 
  
if(PMFsearch){   
    imdata <- Load_Cardinal_imaging(datafile[z],preprocessing = F,resolution = ppm,rotate = rotate[z],as="MSImageSet")
    if (dir.exists(paste0(datafile[z] ," ID"))==FALSE){dir.create(paste0(datafile[z] ," ID"))}
    setwd(paste0(datafile[z] ," ID"))
    #cl=makeCluster(8)
    Peptide_Summary_file_regions<-data.frame()
    message(paste("PMFsearch",name))
    message(paste( "region",names(x),sep=" ",collapse = "\n"))
    for (SPECTRUM_batch in names(x)){
      if (dir.exists(paste0(datafile[z] ," ID/",SPECTRUM_batch,"/"))==FALSE){dir.create(paste0(datafile[z] ," ID/",SPECTRUM_batch,"/"))}
    #PMFsearch_para<-function(SPECTRUM_batch,x,imdata,name,ppm,cl,Peptide_Summary_file,Peptide_Summary_file_regions){
    imdata_ed <- batchProcess(imdata, normalize=FALSE, smoothSignal=F, reduceBaseline=F,
                           peakPick=T, peakAlign=FALSE, pixel=pixels(imdata)[unlist(x[SPECTRUM_batch])])
    #selectedpixel<-get_coord_info(names(pixels(imdata)[unlist(x[SPECTRUM_batch])]))

    if(class(imdata_ed)=="MSImageSet"){
    spectrum_file_table=imdata_ed@featureData@data
    savename=paste(name,SPECTRUM_batch)
    message(paste("PMF_analysis",name,"region",SPECTRUM_batch))
    peaklist<-imdata_ed@featureData@data
    peaklist<-peaklist[peaklist$mean>0,]
    
    colnames(peaklist)<-c("m.z","intensities")
    
    deconv_peaklist<-isopattern_ppm_filter_peaklist(peaklist,ppm=ppm)
    
    message(paste(nrow(deconv_peaklist),"mz features found in the spectrum") )
    #MassSpecWavelet_fun(peaklist = peaklist)
    mz_feature_list<-Do_PMF_search(deconv_peaklist,Peptide_Summary_searchlist,BPPARAM=BPPARAM,ppm = ppm)
    
    mz_feature_list<-unique(mz_feature_list)
    
    message("Iterating peptide information...")
    
    Peptide_Summary_searchlist<-as.data.table(Peptide_Summary_searchlist)
    
    mz_feature_list<-as.data.table(mz_feature_list)
    #Peptide_Summary_searchlist<-Peptide_Summary_searchlist2
    Peptide_Summary_searchlist<-Peptide_Summary_searchlist[,c("mz","Peptide","adduct","formula","isdecoy")]
    
    mz_feature_list$mz<-as.character(mz_feature_list$mz)
    
    Peptide_Summary_searchlist$mz<-as.character(Peptide_Summary_searchlist$mz)
    Peptide_Summary_searchlist<-merge(Peptide_Summary_searchlist,mz_feature_list,by.x="mz",by.y="mz",all.x=T,sort=F)
    Peptide_Summary_searchlist$Intensity[is.na(Peptide_Summary_searchlist$Intensity)]<-0
    
    #bplapply( mz_feature_list$mz, function(x,))
   
    #Peptide_Summary_searchlist$Intensity<-
    #Peptide_Summary_searchlist$Intensity<-unlist(bplapply(Peptide_Summary_searchlist$mz,intensity_sum_para,mz_feature_list,BPPARAM = BPPARAM))
    
    #Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+Peptide_Summary_searchlist$Intensity
    #Peptide_plot_list<-Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity>0,]
    
    Peptide_plot_list<-Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity>=max(Peptide_Summary_searchlist$Intensity)*threshold,]
    #Peptide_plot_list<-unique(Peptide_plot_list)
    
    if(is.null(Peptide_plot_list$moleculeNames)){Peptide_plot_list$moleculeNames=Peptide_plot_list$Peptide}
    Peptide_plot_list$Region=SPECTRUM_batch
    Peptide_plot_list=Peptide_plot_list[(!is.na(Peptide_plot_list$Intensity)),]
    message(paste("1st run returns",nrow(Peptide_plot_list)))
    write.csv(Peptide_plot_list,paste0(datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_1st_ID.csv"),row.names = F)
    
    #Peptide_plot_list<-Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity>=(max(Peptide_Summary_searchlist$Intensity)*threshold),]
    
    
    
    #Peptide_plot_list$moleculeNames<-Peptide_plot_list$Peptide
    if (scoring){
      library(enviPat)
      library(ggplot2)
      data(isotopes)
      
      
      Peptide_plot_list<-read.csv(paste0(datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_1st_ID.csv"),stringsAsFactors = F)
      #message("Scoring peptide...")
      #message(head(Peptide_plot_list$formula))
      #message(head(peaklist))
      #message(head(isotopes))
      
      Peptide_plot_list_Score=bplapply(Peptide_plot_list$formula,SCORE_PMF,peaklist=peaklist,isotopes=isotopes,score_method=score_method,charge = 1,ppm=ppm,BPPARAM = BPPARAM)
      #testscore=lapply("C50H77N11O10S1Na1",SCORE_PMF,score_method=score_method,peaklist=peaklist,isotopes=isotopes,charge = 1,ppm=ppm)
      #Peptide_plot_list_Score=lapply(testdf$formula,SCORE_PMF,peaklist=peaklist,isotopes=isotopes,score_method=score_method,charge = 1,ppm=ppm)
      #Peptide_plot_list_Score_decoyisotope=bplapply(Peptide_plot_list$formula,SCORE_PMF,peaklist=peaklist,isotopes=decoy_isotopes,score_method=score_method,charge = 1,ppm=ppm,BPPARAM = BPPARAM)
      #Peptide_plot_list_Score_decoyisotope_N=bplapply(Peptide_plot_list$formula,SCORE_PMF,peaklist=peaklist,isotopes=decoy_isotopes_N,score_method=score_method,charge = 1,ppm=ppm,BPPARAM = BPPARAM)
      
      #Peptide_plot_list_Score=bplapply(Peptide_plot_list$formula,SCORE_PMF,peaklist,BPPARAM = BPPARAM)
      Peptide_plot_list_Score=unlist(Peptide_plot_list_Score)
      #Peptide_plot_list_Score_decoyisotope=unlist(Peptide_plot_list_Score_decoyisotope)
      #Peptide_plot_list_decoy=Peptide_plot_list
      #Peptide_plot_list_decoy$isdecoy=1
      #Peptide_plot_list_decoy$Score=Peptide_plot_list_Score_decoyisotope
      Peptide_plot_list$Score=Peptide_plot_list_Score
       #Peptide_plot_list_Score_frame=do.call(rbind,Peptide_plot_list_Score)
      #Peptide_plot_list<-rbind(Peptide_plot_list,Peptide_plot_list_decoy)
      if (Decoy_search && ("isotope" %in% Decoy_mode)){
      decoy_isotopes=isotopes
      decoy_isotopes[decoy_isotopes$isotope=="13C",]=data.frame(element="C",isotope="11C",mass=10.99664516,aboudance=0.0107,ratioC=0,stringsAsFactors = F)
      #decoy_isotopes_N=isotopes
      decoy_isotopes[decoy_isotopes$isotope=="15N",]=data.frame(element="N",isotope="13N",mass=13.00603905,aboudance=0.00364000,ratioC=4,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="2H",]=data.frame(element="H",isotope="0H",mass=0.001548286,aboudance=0.00011500,ratioC=6,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="17O",]=data.frame(element="O",isotope="15O",mass=14.99069774,aboudance=0.00038000,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="18O",]=data.frame(element="O",isotope="14O",mass=13.99066884,aboudance=0.00205000,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="33S",]=data.frame(element="S",isotope="31S",mass=30.97268292,aboudance=0.0075,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="34S",]=data.frame(element="S",isotope="30S",mass=29.9762745,aboudance=0.0425,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="35S",]=data.frame(element="S",isotope="29S",mass=28.94414146,aboudance=0,ratioC=3,stringsAsFactors = F)
      decoy_isotopes[decoy_isotopes$isotope=="36S",]=data.frame(element="S",isotope="28S",mass=27.97706058,aboudance=0.0001,ratioC=3,stringsAsFactors = F)
      Peptide_plot_list_Score_decoy=bplapply(Peptide_plot_list$formula,SCORE_PMF,peaklist=peaklist,isotopes=decoy_isotopes,score_method=score_method,charge = 1,ppm=ppm,BPPARAM = BPPARAM)
      Peptide_plot_list_Score_decoy=unlist(Peptide_plot_list_Score_decoy)
      Peptide_plot_list_decoy=Peptide_plot_list
      Peptide_plot_list_decoy$isdecoy=1
      Peptide_plot_list_decoy$Score=Peptide_plot_list_Score_decoy
      Peptide_plot_list<-rbind(Peptide_plot_list,Peptide_plot_list_decoy)
      Protein_feature_list_decoy<-Protein_feature_list
      Protein_feature_list_decoy$isdecoy=1
      Protein_feature_list<-rbind(Protein_feature_list,Protein_feature_list_decoy)
      }
      
      Peptide_plot_list_rank=rank_mz_feature(Peptide_plot_list,mz_feature=deconv_peaklist,BPPARAM = BPPARAM)
      #Peptide_plot_list_rank=Peptide_plot_list_rank[Peptide_plot_list_rank$Rank<=Rank,]
      Peptide_plot_list_rank<-Peptide_plot_list_rank[,c("mz","Peptide","adduct","formula","isdecoy","Intensity","moleculeNames","Region","Score",        
                                                        "mz_align","Rank")]
      
      Protein_feature_result<-protein_scoring(Protein_feature_list,Peptide_plot_list_rank)
      Score_cutoff_protein= FDR_cutoff_plot_protein(Protein_feature_result,FDR_cutoff=0.1,plot_fdr=T,outputdir=paste0(datafile[z] ," ID/",SPECTRUM_batch),adjust_score = F)
      Protein_feature_result_cutoff=Protein_feature_result[((Protein_feature_result$Score>=Score_cutoff_protein)&(!is.na(Protein_feature_result$Intensity))&(Protein_feature_result$isdecoy==0)),]
      
      
      write.csv(Protein_feature_result,paste0(datafile[z] ," ID/",SPECTRUM_batch,"/Protein_ID_score_rank_",score_method,".csv"),row.names = F)
      write.csv(Protein_feature_result_cutoff,paste("Protein_segment_PMF_RESULT",SPECTRUM_batch,".csv"),row.names = F)
      #group_by(c("Protein","isdecoy"))  %>%  summarize(sum("Intensity"))
      
      write.csv(Peptide_plot_list_rank,paste0(datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_1st_ID_score_rank",score_method,".csv"),row.names = F)
      Score_cutoff= FDR_cutoff_plot(Peptide_plot_list_rank,FDR_cutoff=0.1,plot_fdr=T,outputdir=paste0(datafile[z] ," ID/",SPECTRUM_batch),adjust_score = F)
      Peptide_plot_list_2nd=Peptide_plot_list_rank[((Peptide_plot_list_rank$Score>=Score_cutoff)&(!is.na(Peptide_plot_list_rank$Intensity))&(Peptide_plot_list_rank$isdecoy==0)),]

      write.csv(Peptide_plot_list_2nd,paste0(datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_2nd_ID_score_rank",score_method,"_Rank_above_",Rank,".csv"),row.names = F)
      plot_matching_score(Peptide_plot_list_2nd,peaklist,charge=1,ppm,outputdir=paste0(datafile[z] ," ID/",SPECTRUM_batch,"/ppm"))
      Peptide_plot_list_2nd$mz=as.numeric(as.character(Peptide_plot_list_2nd$mz))
    } 
    
    Plot_PMF_all(Peptide_plot_list_2nd,peaklist=peaklist,threshold=threshold,savename)
    #uniques_intensity<-unique(Peptide_plot_list_2nd[,c("mz","Intensity")])
    #uniques_intensity<-uniques_intensity[uniques_intensity$Intensity>0,]

    write.csv(Peptide_plot_list_2nd,paste("Peptide_segment_PMF_RESULT",SPECTRUM_batch,".csv"),row.names = F)
    write.csv(peaklist,paste("Spectrum",SPECTRUM_batch,".csv"),row.names = F)
    #Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+unlist(parLapply(cl=cl,Peptide_Summary_file$mz,intensity_sum_para,uniques_intensity))
    #Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+unlist(bplapply(Peptide_Summary_file$mz,intensity_sum_para,uniques_intensity,BPPARAM = BPPARAM))
   #Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+unlist(bplapply(Peptide_Summary_file$mz,intensity_sum_para,uniques_intensity,BPPARAM = BPPARAM))
    
    Peptide_Summary_file_regions<-rbind(Peptide_Summary_file_regions,Peptide_plot_list_2nd)
    #return(list(Peptide_Summary_file,Peptide_Summary_file_regions))
    }
    }
        
    #PMF_result_file=parallel::parLapply(cl,names(x),PMFsearch_para,x,imdata,name,ppm,cl,Peptide_Summary_file,Peptide_Summary_file_regions)
    
  if(is.null(Peptide_Summary_file$Peptide)){Peptide_Summary_file$Peptide=Peptide_Summary_file$moleculeNames}
  if(is.null(Peptide_Summary_file$moleculeNames)){Peptide_Summary_file$moleculeNames=Peptide_Summary_file$Peptide}
  
  Peptide_Summary_file<-Peptide_regions_Summary_fun(Peptide_Summary_file_regions)
  write.csv(Peptide_Summary_file,"Peptide_Summary_file.csv",row.names = F)
  write.csv(Peptide_Summary_file_regions,"Peptide_region_file.csv",row.names = F)
  }
  
  }
  else{
    message("Spectrum generation is bypassed, will perform PMF search via existing spectra in ID folder...")
    message(datafile[z])
    if(PMFsearch){   
      if (dir.exists(paste0(datafile[z] ," ID"))==FALSE){dir.create(paste0(datafile[z] ," ID"))}
      setwd(paste0(datafile[z] ," ID"))
      match_pattern <- "Spectrum.{2,}csv"
      name <-gsub(base::dirname(datafile[z]),"",datafile[z])
      folder<-base::dirname(datafile[z])
      message(paste("Found spectrum file:",dir()[str_detect(dir(), match_pattern)]),collapse="\n")
      
      for (SPECTRUM_batch in dir()[str_detect(dir(), match_pattern)]){
        peaklist=fread(SPECTRUM_batch)
        SPECTRUM_batch=gsub("^Spectrum.","",SPECTRUM_batch)
        SPECTRUM_batch=gsub(".csv$","",SPECTRUM_batch)
        SPECTRUM_batch=gsub(" ","",SPECTRUM_batch)
      message(paste("PMFsearch"))
      message(paste( "region",SPECTRUM_batch,sep=" ",collapse = "\n"))
      savename=paste(name,SPECTRUM_batch)
      if (dir.exists(paste0(datafile[z] ," ID/",SPECTRUM_batch,"/"))==FALSE){dir.create(paste0(datafile[z] ," ID/",SPECTRUM_batch,"/"))}
           #PMFsearch_para<-function(SPECTRUM_batch,x,imdata,name,ppm,cl,Peptide_Summary_file,Peptide_Summary_file_regions){
          message(paste("PMF_analysis","region",SPECTRUM_batch))
          peaklist<-peaklist[peaklist$intensities!=0,]
          colnames(peaklist)<-c("m.z","intensities")
          mz_feature_list<-Do_PMF_search(peaklist,Peptide_Summary_searchlist,BPPARAM=BPPARAM)
          message("Iterating peptide information...")
          Peptide_Summary_searchlist$Intensity<-unlist(bplapply(Peptide_Summary_searchlist$mz,intensity_sum_para,mz_feature_list,BPPARAM = BPPARAM))
          Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+Peptide_Summary_searchlist$Intensity
          Peptide_plot_list<-Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity>300,]
          #Peptide_plot_list<-Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity>=(max(Peptide_Summary_searchlist$Intensity)*threshold),]
          Peptide_plot_list$Region=SPECTRUM_batch
          write.csv(Peptide_plot_list,paste0(datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_1st_ID.csv"),row.names = F)
          if(is.null(Peptide_plot_list$moleculeNames)){Peptide_plot_list$moleculeNames=Peptide_plot_list$Peptide}
          #Peptide_plot_list$moleculeNames<-Peptide_plot_list$Peptide
          if (scoring){
            library(enviPat)
            library(ggplot2)
            data(isotopes)
            
            Peptide_plot_list_Score=bplapply(Peptide_plot_list$formula,SCORE_PMF,peaklist,isotopes)
            #Peptide_plot_list_Score=bplapply(Peptide_plot_list$formula,SCORE_PMF,peaklist,BPPARAM = BPPARAM)
            Peptide_plot_list_Score=unlist(Peptide_plot_list_Score)
            #Peptide_plot_list_Score_frame=do.call(rbind,Peptide_plot_list_Score)
            Peptide_plot_list$Score=Peptide_plot_list_Score
            FDR_cutoff=0.1
            Score_cutoff= FDR_cutoff_plot(Peptide_plot_list,FDR_cutoff=FDR_cutoff,FDR_strip=0.002,plot_fdr=T,outputdir=paste0(datafile[z] ," ID/",SPECTRUM_batch))
            
            Peptide_plot_list=Peptide_plot_list[Peptide_plot_list$Score>=Score_cutoff,]
            message(paste("found",nrow(Peptide_plot_list),"at FDR cut off:", FDR_cutoff, "and score cut off:", Score_cutoff) )
            
            write.csv(Peptide_plot_list,paste0(datafile[z] ," ID/",SPECTRUM_batch,"/Peptide_2nd_ID.csv"),row.names = F)
            plot_matching_score(Peptide_plot_list,peaklist,charge=1,ppm,outputdir=paste0(datafile[z] ," ID/",SPECTRUM_batch))
            
          } 
          Plot_PMF_all(Peptide_plot_list,peaklist=peaklist,threshold=threshold,savename)
         
          write.csv(Peptide_plot_list,paste("Peptide_segment_PMF_RESULT",SPECTRUM_batch,".csv"),row.names = F)
          #write.csv(peaklist,paste("Spectrum",SPECTRUM_batch,".csv"),row.names = F)
          #Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+unlist(parLapply(cl=cl,Peptide_Summary_file$mz,intensity_sum_para,uniques_intensity))
          #Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+unlist(bplapply(Peptide_Summary_file$mz,intensity_sum_para,uniques_intensity,BPPARAM = BPPARAM))
          #Peptide_Summary_file$Intensity<-Peptide_Summary_file$Intensity+unlist(bplapply(Peptide_Summary_file$mz,intensity_sum_para,uniques_intensity,BPPARAM = BPPARAM))
          Peptide_Summary_file_regions<-rbind(Peptide_Summary_file_regions,Peptide_plot_list)
          #return(list(Peptide_Summary_file,Peptide_Summary_file_regions))
    
      }
      
      #PMF_result_file=parallel::parLapply(cl,names(x),PMFsearch_para,x,imdata,name,ppm,cl,Peptide_Summary_file,Peptide_Summary_file_regions)
      
      if(is.null(Peptide_Summary_file$Peptide)){Peptide_Summary_file$Peptide=Peptide_Summary_file$moleculeNames}
      if(is.null(Peptide_Summary_file$moleculeNames)){Peptide_Summary_file$moleculeNames=Peptide_Summary_file$Peptide}
      
      Peptide_Summary_file<-Peptide_Summary_file[Peptide_Summary_file$Intensity>0,]
      
     
      
      write.csv(Peptide_Summary_file,"Peptide_Summary_file.csv",row.names = F)
      write.csv(Peptide_Summary_file_regions,"Peptide_region_file.csv",row.names = F)
    }  
  } 
  

  
}
  
}


intensity_sum_para<-function(mz,Peptide_plot_list){
  
  ifelse(length(Peptide_plot_list[Peptide_plot_list$mz==mz,"Intensity"])!=0,Peptide_plot_list[Peptide_plot_list$mz==mz,"Intensity"],0)
}

Cardinal<- function(){
  workdir <- WorkingDir()
  datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
  
  listfile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.csv$")
  datafile<-gsub(".imzML", "", datafile)
  name <-gsub(paste(workdir,"\\/",sep=""),"",datafile[3])
  
  #imdata <- readImzML(name, folder, as="MSImagingExperiment")
  imdata <- readImzML(name, folder, attach.only=TRUE, as="MSImagingExperiment")
  summarize(data, sum, .by="pixel")
  tmp <- imdata %>%
    smoothSignal() %>%
    reduceBaseline() %>%
    peakPick() %>%
    peakFilter() %>%
    select(x == 1, y == 1)
  process(plot=TRUE,
          par=list(layout=c(1,3)),
          BPPARAM=SerialParam())
  
  plot(data, pixel=1)
  plot(data, coord=list(x=2, y=2))
  
  
  pattern <- factor(c(0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0,
                      0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 2, 1, 1, 2,
                      2, 0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 0, 0, 0, 0, 1, 2, 2,
                      2, 2, 2, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 2,
                      2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0),
                    levels=c(0,1,2), labels=c("blue", "black", "red"))
  tmp <- imdata %>%
    smoothSignal() %>%
    reduceBaseline() %>%
    peakPick() %>%
    peakFilter() %>%
    select(x == 1, y == 1) %>%
    process(plot=TRUE,
            par=list(layout=c(1,3)),
            BPPARAM=SerialParam())
}

intensity.colors_customize <- function(n = 100, alpha = 1) {
  col2 <- rainbow(3*n, alpha=alpha)[(2*n):1]
  f <- colorRamp(colors=c("black", "blue", "green", "yellow", "red","#FF00FF","white"))
  alpha <- col2rgb(col2, alpha=TRUE)[[4]]
  col1 <- sapply(seq(from=0, to=1, length.out=n), function(i) do.call(rgb,
                                                                      c(as.list(f(i)), maxColorValue=255, alpha=alpha)))
  col1
}

Pathway_overview_graphite<-function(){
  p_load(graphite,graph )
  humanReactome <- graphite::pathways("hsapiens", "reactome")
  humanmeatbolome <- graphite::pathways("hsapiens", "smpdb")
  metab_url <-
    url("https://romualdi.bio.unipd.it/wp-uploads/2018/04/Terunuma_metabolite_expr.txt")
  mexpr <- read.table(metab_url, header = TRUE, sep = "\t", row.names = NULL,
                      stringsAsFactors = FALSE)
  
  p <- humanReactome[["Glycolysis"]]
  head(nodes(p))
  head(edges(p))
  head(nodes(p), which = "mixed")
  head(edges(p), which = "mixed")
  pathwayDatabases()
  g <- pathwayGraph(p)
  g <- pathwayGraph(p, which = "mixed")
  pSymbol <- convertIdentifiers(p, "SYMBOL")
  head(nodes(pSymbol))
  reactomeSymbol <- convertIdentifiers(humanReactome[1:5], "SYMBOL")
  cytoscapePlot(convertIdentifiers(reactome$`Unwinding of DNA`, "symbol"), which = "mixed")
  meta_g <- pathwayGraph(p, which = "metabolites")
  
}

SCORE_PMF<-function(formula,peaklist,isotopes=NULL,threshold=2.5,charge=1,ppm=5,print.graphic=F,output.list=F,outputfile=NULL,score_method="SQRT"){
  library(rcdk)
  library(rcdklibs)
  library(OrgMassSpecR)
  library(enviPat)
  library(data.table)
  library(rJava)
  library(grid)
  if (is.null(isotopes)){data("isotopes")}
  
  
  
  Spectrum_scoring <- function(spec.top, spec.bottom, t = 0.25, b = 0, top.label = NULL, 
                               bottom.label = NULL, xlim = c(50, 1200), x.threshold = 0, print.alignment = FALSE,
                               print.graphic = F, output.list = F,score_method="SQRT",Peak_intensity=0,outputfile=NULL) {
    
    
    ## format spectra and normalize intensitites
    
    top_tmp <- data.frame(mz = spec.top[,1], intensity = spec.top[,2])
    top_tmp$normalized <- round((top_tmp$intensity / max(top_tmp$intensity)) ,digits = 3)
    #top_tmp <- subset(top_tmp, top_tmp$mz >= xlim[1] & top_tmp$mz <= xlim[2])   
    top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)   # data frame for plotting spectrum
    top <- subset(top_plot, top_plot$intensity > b)   # data frame for similarity score calculation
    
    bottom_tmp <- data.frame(mz = spec.bottom[,1], intensity = spec.bottom[,2])
    bottom_tmp$normalized <- round((bottom_tmp$intensity / max(bottom_tmp$intensity)) ,digits = 3)
    #bottom_tmp <- subset(bottom_tmp, bottom_tmp$mz >= xlim[1] & bottom_tmp$mz <= xlim[2])   
    bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)   # data frame for plotting spectrum
    bottom <- subset(bottom_plot, bottom_plot$intensity > b)   # data frame for similarity score calculation
    #top <- spec.top
    #bottom <- spec.bottom
    
    ## align the m/z axis of the two spectra, the bottom spectrum is used as the reference
    
   # for(i in 1:nrow(bottom))
   #   top[,1][bottom[,1][i] >= top[,1] - t & bottom[,1][i] <= top[,1] + t] <- bottom[,1][i]
    alignment <- merge(top, bottom, by = 1, all = TRUE)
    #if(length(unique(alignment[,1])) != length(alignment[,1])) warning("the m/z tolerance is set too high")
    alignment[,c(2,3)][is.na(alignment[,c(2,3)])] <- 0   # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
    names(alignment) <- c("mz", "intensity.top", "intensity.bottom")
    
    if(print.alignment == TRUE) {
      print(alignment)
    }
    
    
    ## similarity score calculation
    
    #if(x.threshold < 0) stop("x.threshold argument must be zero or a positive number")
    
    #alignment <- alignment[alignment[,1] >= x.threshold, ]
    
    match_score <- sum(`&`(alignment$intensity.top>0,alignment$intensity.bottom>0))/sum(alignment$intensity.top>0)
    
   # u <- alignment[,2]; v <- alignment[,3]
    
    if (score_method=="SQRT"){
      u <- alignment[,2]; v <- alignment[,3]
      #u<-u/100; v<-v/100
      
      #sum( u*v)
      #similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))) 
      
      similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2)) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))))
      
      #similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2)) ))
      
    }else if (score_method=="SQRTP"){
      u <- alignment[,2]; v <- alignment[,3]
      #u<-u/100; v<-v/100
      
      #sum( u*v)
      #similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))) 
      
      similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2)) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))))
      match_score <- log(sum(`&`(alignment$intensity.top>0,alignment$intensity.bottom>0))/sum(alignment$intensity.top>0))
      #similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2)) ))
      similarity_score <- similarity_score + match_score
    }else if (score_method=="balanced-SQRT"){
      u <- alignment[,2]; v <- alignment[,3]
      
      #similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))) 
      if (Peak_intensity<10) (Peak_intensity=10)
      
      similarity_score <- -log(as.vector(sqrt(sum(((u-v)*u)^2)) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))) * log(Peak_intensity,10)
      
    }
    else if (score_method=="Equal-SQRT"){
      #score=log(v)/log(u)
      u <- alignment[,2]; v <- alignment[,3]
      score=-log(v/u)
      score=score[score!=Inf]
      similarity_score=-log(sum(abs(score))/length(score))
      
      #similarity_score <- as.vector((u %*% v) * length(v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))) )  
    }
    else if (score_method=="Equal-intensity-SQRT"){
      similarity_score <- as.vector((u %*% v) * length(v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))*log(Peak_intensity)) )  
    }
    else if (score_method=="Mix-SQRT"){

      alignment_0=ifelse(alignment==0,0,1)
      
      matching_score = sum(alignment_0[,2]==alignment_0[,3])/nrow(alignment_0)
      
      alignment1<-alignment[alignment_0[,2]==alignment_0[,3],]
      
      u <- alignment1[,2]; v <- alignment1[,3]
      
      similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))) 
    }
    
    ## generate plot
    
    if(print.graphic == TRUE) {
      library(ggplot2)
      #library(magick)
      #tempfilename=tempfile(pattern = "file", tmpdir = tempdir())
      tempfilename=outputfile
      png(tempfilename, bg = "white")
      plot.new()
      plot.window(xlim = xlim, ylim = c(-125, 125))
      ticks <- c(-100, -50, 0, 50, 100)
      for(i in 1:length(top_plot$mz)) lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i]*100), col = "blue")
      for(i in 1:length(bottom_plot$mz)) lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i]*100), col = "red")
      axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "Intensity")
      axis(1, pos = -125)
      lines(xlim, c(0,0))
      rect(xlim[1], -125, xlim[2], 125)
      mtext("m/z", side = 1, line = 2)
      mtext("intensity (%)", side = 2, line = 2)
      plot.window(xlim = c(0, 20), ylim = c(-10, 10))
      text(10, 9, top.label)
      text(10, -9, bottom.label)
      dev.off()
      #p<-magick::image_read(tempfilename)
    }
    
    if (abs(similarity_score)==Inf) similarity_score = 0
    if(output.list == TRUE) {
      
      # Grid graphics head to tail plot
      
      
    }
    
    if(output.list == FALSE) {
      
      return(similarity_score)
      
    }
    
    if(output.list == TRUE) {
      
      return(list(similarity.score = similarity_score, 
                  plot = p))
      
    }
    
    # simscore <- as.vector((u %*% v)^2 / (sum(u^2) * sum(v^2)))   # cos squared
    
  }
  isopattern_ppm_filter<-function(pattern,ppm){
    
    pattern_ppm=as.numeric(as.character(pattern[,1]))
    
    pattern_ppm_delta=numeric()
    
    
    filtered_pattern<-pattern[1,]
    filtered_pattern<-t(data.frame(filtered_pattern,stringsAsFactors = F))
    #dim(filtered_pattern)
    for (i in 1:(length(pattern_ppm)-1)){
      
      pattern_ppm_delta[i]=(pattern_ppm[i+1]-pattern_ppm[i])/pattern_ppm[i]
      
      if (pattern_ppm_delta[i]>(ppm/1000000)){
        filtered_pattern<-rbind(filtered_pattern,pattern[i+1,])
      } else {
        #previous_iso=filtered_pattern[nrow(filtered_pattern),]
        
        #newline=as.data.frame(list("m/z"=(filtered_pattern[nrow(filtered_pattern),1]*filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,1]*pattern[i+1,2])/(sum(filtered_pattern[nrow(filtered_pattern),2],pattern[i+1,2])),"abundance"=filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,2]))
        #names(newline)=c("m/z","abundance")
        newline=c((filtered_pattern[nrow(filtered_pattern),1]*filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,1]*pattern[i+1,2])/(sum(filtered_pattern[nrow(filtered_pattern),2],pattern[i+1,2])),
                  filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,2])
        filtered_pattern[nrow(filtered_pattern),]<-newline
      }
      
    }
    rownames(filtered_pattern)=1:nrow(filtered_pattern)
    
    return(filtered_pattern)
  }
  #getMolecule(formula)
  #isopattern(chemforms = formula)
  #data(isotopes)
  #message(formula)
  pattern<-isopattern(
    isotopes,
    formula,
    threshold=threshold,
    plotit=F,
    charge=charge,
    emass=0.00054858,
    algo=1,
    verbose = F
  )
  #message(pattern)
  pattern=pattern[[formula]]
  pattern=isopattern_ppm_filter(pattern = pattern[,1:2], ppm=ppm)
  
  #monomass=pattern[1,1]
  #m_1_pattern=data.frame("m/z"=pattern[1,1]-(1.003354840/ifelse(charge==0,1,abs(charge))),abundance=0)
  #colnames(m_1_pattern)=colnames(pattern[,1:2])
  #pattern=rbind(m_1_pattern,pattern[,1:2])
  #rownames(pattern)=1:nrow(pattern)
  #message(charge)
  spectrumintensity=unlist(lapply(1:nrow(pattern), function(x,pattern,peaklist,ppm){
    PMF_spectrum_intensity<-as.numeric(peaklist[between(peaklist$m.z,pattern[x,1]*(1-ppm/1000000),pattern[x,1]*(1+ppm/1000000)),2])
    if(sum(!is.na(PMF_spectrum_intensity),length(PMF_spectrum_intensity)!=0)>1) sum(PMF_spectrum_intensity,na.rm = T) else 0
  },pattern,peaklist,ppm))
  #message(spectrumintensity)
  spectrum=data.frame(mz=pattern[,1],Intensity=spectrumintensity)
  Peak_intensity<-max(spectrum$Intensity)
  #message(head(spectrum))
  #message(c(min(range(pattern[,1],na.rm = T))-1," ",max(range(pattern[,1],na.rm = T))+1))
  score=Spectrum_scoring(spec.top = pattern,spec.bottom = spectrum,b=0,t=mean(pattern[,1])*ppm/1000000,top.label = "Theoretical",score_method=score_method,bottom.label = "Observed spectrum",xlim = c(min(range(pattern[,1],na.rm = T))-1,max(range(pattern[,1],na.rm = T))+1),output.list=output.list,print.graphic = print.graphic,outputfile = outputfile,Peak_intensity=Peak_intensity)
  #score
  #message(score)
  
  return(score)
}



Do_PMF_search<-function(peaklist,Peptide_Summary_searchlist,BPPARAM=bpparam(),ppm=2){
peaklist<-peaklist[peaklist$intensities!=0,]
colnames(peaklist)<-c("m.z","intensities")
mzrange=c(min(peaklist$m.z),max(peaklist$m.z))
#peaklist=peaklist[peaklist$intensities>=(max(peaklist$intensities)/3*threshold),]
Peptide_Summary_searchlist_mz=NULL
uniquemz=as.numeric(as.character(unique(Peptide_Summary_searchlist$mz)))
Peptide_Summary_searchlist_mz$mz=uniquemz[`&`(uniquemz>mzrange[1],uniquemz<mzrange[2])]
Peptide_Summary_searchlist_mz$Intensity=rep(0,length(Peptide_Summary_searchlist_mz$mz))
Peptide_Summary_searchlist_mz=as.data.frame(Peptide_Summary_searchlist_mz)
#Peptide_Summary_searchlist<-PMF_analysis_fun(Peptide_Summary_searchlist=Peptide_Summary_searchlist,peaklist=peaklist,ppm=ppm,BPPARAM=BPPARAM)
Peptide_Summary_searchlist_mz<-PMF_analysis_fun(Peptide_Summary_searchlist=Peptide_Summary_searchlist_mz,peaklist=peaklist,ppm=ppm,BPPARAM=BPPARAM,mzrange=mzrange)
Peptide_Summary_searchlist_mz=as.data.frame(Peptide_Summary_searchlist_mz,stringsAsFactors = F)
#Peptide_feature_list<-Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity>0,]
#Peptide_plot_list<-Peptide_feature_list[Peptide_feature_list$Intensity>=max(Peptide_feature_list$Intensity)*threshold,]
mz_feature_list<-Peptide_Summary_searchlist_mz[Peptide_Summary_searchlist_mz$Intensity>0,]

return(mz_feature_list)
}

FDR_cutoff_plot_protein<-function(Protein_feature_result,FDR_cutoff=0.1,plot_fdr=T,outputdir=paste0(datafile[z] ," ID/",SPECTRUM_batch),adjust_score = F){
  
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(zoo)
  library(FTICRMS)
  Protein_feature_result=data_test_rename(c("isdecoy","Proscore"),Protein_feature_result)
  
  Protein_feature_result$isdecoy<-factor(Protein_feature_result$isdecoy)
  
  #unique_fomula<-unique(Protein_feature_result$formula)
  
  unique_fomula_ID<-unique(Protein_feature_result[,c("Protein","isdecoy","Intensity_norm","Proscore")])
  
  
  Protein_feature_result=unique_fomula_ID
  
  Protein_feature_result$Proscore=Protein_feature_result$Proscore
  

  
  if(adjust_score){
    
    Proscore_boundary<-sapply(1:max(Protein_feature_result$mz),function(x,Protein_feature_result){
      
      max(Protein_feature_result$Proscore[between(Protein_feature_result$mz,x-0.5,x+0.5)])
      
    },Protein_feature_result)
    
    names(Proscore_boundary)=as.character(1:max(Protein_feature_result$mz))
    
    Proscore_boundary<-Proscore_boundary[Proscore_boundary>0]
    
    Proscore_boundary<-Proscore_boundary[!is.na(Proscore_boundary)]
    
    Proscore_boundary<-data.frame(mz=names(Proscore_boundary),top=Proscore_boundary,stringsAsFactors = F)
    
    x <- as.numeric(Proscore_boundary$mz)
    y <- Proscore_boundary$top
    #plot(x,y)
    
    fit3 <- lm( y~poly(as.numeric(as.character(x)),3) )
    
    fitx=Protein_feature_result$mz
    
    fitdata=as.data.frame(stats::predict(fit3, data.frame(x = fitx), interval="confidence"))
    
    Protein_feature_result$Proscore_factor<-fitdata$upr
    
    Protein_feature_result$adjusted_Proscore<-Protein_feature_result$Proscore/Protein_feature_result$Proscore_factor
    Protein_feature_result$adjusted_Proscore<-Protein_feature_result$adjusted_Proscore/max(Protein_feature_result$adjusted_Proscore,na.rm = T)
    Protein_feature_result$original_Proscore<-Protein_feature_result$Proscore
    Protein_feature_result$Proscore<-Protein_feature_result$adjusted_Proscore
    if (!is.null(outputdir) && plot_fdr){
      png(paste0(outputdir,"/Matching_PROTEIN_Score_vs_mz_after_adjustment.png"))
      p<-ggplot(data=Protein_feature_result,aes(x=mz,y=Proscore,color=isdecoy)) + geom_point(size=1,alpha=1/10) + 
        ggtitle("Matching score vs mz after adjustment") +
        xlab("mz") + ylab("Matching score") + labs(fill = "isdecoy")
      print(p)
      dev.off() 
    }
    
  }
  
  target_Proscore<-Protein_feature_result$Proscore[Protein_feature_result$isdecoy!=1]
  decoy_Proscore<-Protein_feature_result$Proscore[Protein_feature_result$isdecoy==1]
  target_Proscore[is.na(target_Proscore)]=0
  decoy_Proscore[is.na(decoy_Proscore)]=0
  target_Proscore[target_Proscore==Inf]<-0
  decoy_Proscore[decoy_Proscore==Inf]<-0
  
  breaks = seq(min(Protein_feature_result$Proscore), max(Protein_feature_result$Proscore), by=max(Protein_feature_result$Proscore)/FDR_strip)
  target_Proscore.cut = cut(target_Proscore, breaks, right=T) 
  decoy_Proscore.cut = cut(decoy_Proscore, breaks,right=T) 
  target_Proscore.freq = table(target_Proscore.cut)
  decoy_Proscore.freq = table(decoy_Proscore.cut)
  target_Proscore.freq= target_Proscore.freq[rev(names(target_Proscore.freq))]
  target_Proscore.freq0 = c(0, cumsum(target_Proscore.freq))
  decoy_Proscore.freq= decoy_Proscore.freq[rev(names(decoy_Proscore.freq))]
  decoy_Proscore.freq0 = c(0, cumsum(decoy_Proscore.freq))
  
  
  
  FDR=decoy_Proscore.freq0/target_Proscore.freq0
  
  df=data.frame(breaks=rev(breaks),target_Proscore.freq0=target_Proscore.freq0 ,decoy_Proscore.freq0=decoy_Proscore.freq0,FDR=FDR)
  #df=df[order(breaks,decreasing = TRUE),]
  df$FDR_m.av<-rollmean(df$FDR, 3,fill = list(0, NaN, NaN))
  #df$logProscore=log(df$breaks)
  
  Proscore_cutoff<-min(df$breaks[(df$FDR_m.av<=FDR_cutoff)==T],na.rm = T)
  
  if (!is.null(outputdir) && plot_fdr){
    
    png(paste0(outputdir,"/protein_FDR.png"))
    plot.new()
    plot(df$breaks, df$FDR,            # plot the data 
         main="FDR plot",  # main title 
         xlab="Matching Proscore",        # x−axis label 
         ylab="FDR")   # y−axis label 
    lines(df$breaks, df$FDR_m.av)
    dev.off()  
    
    write.csv(df,paste0(outputdir,"/PROTEIN_FDR.CSV"),row.names = F)  
    
  }
  
  if (!is.null(outputdir) && plot_fdr){
    Protein_feature_result_plot<-Protein_feature_result[,c("isdecoy","Proscore")]
    target_decoy<-factor(ifelse(Protein_feature_result$isdecoy==1,"Decoy","Tagret"))
    Protein_feature_result_plot$target_decoy=target_decoy
    png(paste0(outputdir,"/PROTEIN_Score_histogram.png"))
    p<-ggplot(data=Protein_feature_result_plot,aes(x=Proscore,color=target_decoy, fill=target_decoy)) + geom_histogram( fill="white",alpha=0.5, bins = 500, position="identity")  +
      ggtitle("Protein score vs Counts") +
      xlab("Protein score") + ylab("Counts") + labs(fill = "Is_Decoy") + theme_classic() + facet_grid(target_decoy ~ .)
    print(p)
    dev.off() 
  }
  
  return(Proscore_cutoff)
  
  

}


FDR_cutoff_plot<-function(Peptide_plot_list,FDR_cutoff=0.1,FDR_strip=500,plot_fdr=F,outputdir=NULL,adjust_score=F){
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(zoo)
  library(FTICRMS)
  Peptide_plot_list=data_test_rename(c("isdecoy","Score"),Peptide_plot_list)
  
  Peptide_plot_list$isdecoy<-factor(Peptide_plot_list$isdecoy)
  
  #unique_fomula<-unique(Peptide_plot_list$formula)
  
  unique_fomula_ID<-unique(Peptide_plot_list[,c("mz","formula","isdecoy","Intensity","Score")])
  
  
  Peptide_plot_list=unique_fomula_ID
  
  if (!is.null(outputdir) && plot_fdr){
    png(paste0(outputdir,"/Matching_Score_vs_mz_before_adjustment.png"))
    p<-ggplot(data=Peptide_plot_list,aes(x=mz,y=Score,color=isdecoy)) + geom_point(size=1,alpha=1/round(nrow(Peptide_plot_list)/1000)) + 
      ggtitle("Matching score vs mz before adjustment") +
      xlab("mz") + ylab("Matching score") + labs(fill = "isdecoy")+ theme(axis.text.x = element_text(angle=45))
    print(p)
    dev.off() 
  }
  
  if(adjust_score){
    
    score_boundary<-sapply(1:max(Peptide_plot_list$mz),function(x,Peptide_plot_list){
      
      max(Peptide_plot_list$Score[between(Peptide_plot_list$mz,x-0.5,x+0.5)])
      
    },Peptide_plot_list)
    
    names(score_boundary)=as.character(1:max(Peptide_plot_list$mz))
    
    score_boundary<-score_boundary[score_boundary>0]
    
    score_boundary<-score_boundary[!is.na(score_boundary)]
    
    score_boundary<-data.frame(mz=names(score_boundary),top=score_boundary,stringsAsFactors = F)
    
    x <- as.numeric(score_boundary$mz)
    y <- score_boundary$top
    #plot(x,y)
    
    fit3 <- lm( y~poly(as.numeric(as.character(x)),3) )
    
    fitx=Peptide_plot_list$mz
    
    fitdata=as.data.frame(stats::predict(fit3, data.frame(x = fitx), interval="confidence"))
    
    Peptide_plot_list$score_factor<-fitdata$upr
    
    Peptide_plot_list$adjusted_score<-Peptide_plot_list$Score/Peptide_plot_list$score_factor
    Peptide_plot_list$adjusted_score<-Peptide_plot_list$adjusted_score/max(Peptide_plot_list$adjusted_score,na.rm = T)
    Peptide_plot_list$original_score<-Peptide_plot_list$Score
    Peptide_plot_list$Score<-Peptide_plot_list$adjusted_score
    if (!is.null(outputdir) && plot_fdr){
      png(paste0(outputdir,"/Matching_Score_vs_mz_after_adjustment.png"))
      p<-ggplot(data=Peptide_plot_list,aes(x=mz,y=Score,color=isdecoy)) + geom_point(size=1,alpha=1/10) + 
        ggtitle("Matching score vs mz after adjustment") +
        xlab("mz") + ylab("Matching score") + labs(fill = "isdecoy")
      print(p)
      dev.off() 
    }
    
  }
  
  target_score<-Peptide_plot_list$Score[Peptide_plot_list$isdecoy!=1]
  decoy_score<-Peptide_plot_list$Score[Peptide_plot_list$isdecoy==1]
  target_score[is.na(target_score)]=0
  decoy_score[is.na(decoy_score)]=0
  target_score[target_score==Inf]<-0
  decoy_score[decoy_score==Inf]<-0
  
  breaks = seq(0, max(Peptide_plot_list$Score), by=max(Peptide_plot_list$Score)/FDR_strip)
  target_score.cut = cut(target_score, breaks, right=T) 
  decoy_score.cut = cut(decoy_score, breaks,right=T) 
  target_score.freq = table(target_score.cut)
  decoy_score.freq = table(decoy_score.cut)
  target_score.freq= target_score.freq[rev(names(target_score.freq))]
  target_score.freq0 = c(0, cumsum(target_score.freq))
  decoy_score.freq= decoy_score.freq[rev(names(decoy_score.freq))]
  decoy_score.freq0 = c(0, cumsum(decoy_score.freq))
  

  
  FDR=decoy_score.freq0/target_score.freq0
  
  df=data.frame(breaks=rev(breaks),target_score.freq0=target_score.freq0 ,decoy_score.freq0=decoy_score.freq0,FDR=FDR)
  #df=df[order(breaks,decreasing = TRUE),]
  df$FDR_m.av<-rollmean(df$FDR, 3,fill = list(0, NaN, NaN))
  df$logscore=log(df$breaks)
  
  Score_cutoff<-min(df$breaks[(df$FDR_m.av<=FDR_cutoff)==T],na.rm = T)
  
  if (!is.null(outputdir) && plot_fdr){
  
  png(paste0(outputdir,"/FDR.png"))
  plot.new()
  plot(df$breaks, df$FDR,            # plot the data 
       main="FDR plot",  # main title 
       xlab="Matching score",        # x−axis label 
       ylab="FDR")   # y−axis label 
  lines(df$breaks, df$FDR_m.av)
  dev.off()  
  
  write.csv(df,paste0(outputdir,"/FDR.CSV"),row.names = F)  
    
  }

  
  return(Score_cutoff)
  
}

plot_matching_score<-function(Peptide_plot_list,peaklist,charge,ppm,outputdir=getwd()){
  library(enviPat)
  data("isotopes")
  if (dir.exists(outputdir)==FALSE){dir.create(outputdir)}
  Peptide_plot_list1=Peptide_plot_list
  Peptide_plot_list=Peptide_plot_list[!is.na(Peptide_plot_list$Intensity),]
  for (i in 1:nrow(Peptide_plot_list)){
  SCORE_PMF(Peptide_plot_list$formula[i],peaklist,isotopes=isotopes,threshold=0.01,charge=charge,ppm=ppm,print.graphic=T,output.list=F,outputfile=paste0(outputdir,"/",Peptide_plot_list$Peptide[i]," ",Peptide_plot_list$adduct[i],".png"))
  }
  
  
}

Mass_defect_plot<-function(Protein_feature_list,outputdir=NULL){
  
 library(enviPat)
 library(ggplot2)
 data("isotopes")
  isotopes$nominalmass=as.numeric(sapply(isotopes$isotope,function(x,pattern){
    library(stringr)
    paste0(unlist(str_extract_all(x,pattern)),collapse = "")
    },pattern="[0-9]"))
  nomi_isotopes=isotopes[isotopes$abundance>0.5,]
  Protein_feature_list1=data_test_rename(c("isdecoy","mz","adduct","formula"),Protein_feature_list) 
  nomial_mass=Protein_feature_list1$formula
  formulalist=bplapply(nomial_mass,get_atoms,BPPARAM = BPPARAM)
  
  nominalmass=bplapply(formulalist,function(x,isotopes){
    nominalmass=0
    for( atom in names(x)){
      nominalmass=nominalmass+(isotopes$nominalmass[isotopes$element==atom]*x[[atom]])
    }
    nominalmass
  },isotopes=nomi_isotopes,BPPARAM = BPPARAM)
  
  nominalmass=unlist(nominalmass)
  Protein_feature_list1$nominalmass=nominalmass
  Protein_feature_list1$massdefect=Protein_feature_list1$mz-Protein_feature_list1$nominalmass
  
  Protein_feature_list2=Protein_feature_list1[Protein_feature_list1$mz<=4000,]
  Protein_feature_list2$isdecoy=as.factor(Protein_feature_list2$isdecoy)
  p<-ggplot(data=Protein_feature_list2,aes(x=mz,y=massdefect,color=isdecoy)) + geom_point(size=0.1,alpha=1/100) 
  
  if (is.null(outputdir)){
    outputfile="mass_defect_plot.png"
  }else{
    outputfile=paste0(outputdir,"/","mass_defect_plot.png")
  }
  
  png(outputfile,width = 1080,height = 680) 
  p
  dev.off()
}

Peptide_regions_Summary_fun<-function(Peptide_Summary_file_regions){
  
  library(data.table)
  #rowsorg<-nrow(Peptide_Summary_file_regions)
  #Peptide_Summary_file_regions<-rbind(Peptide_Summary_file_regions,Peptide_Summary_file_regions)
  #Peptide_Summary_file_regions$Region<-c(rep(3,rowsorg),rep(2,rowsorg))
  
  peptideregion<-list()
  unique_region<-unique(Peptide_Summary_file_regions$Region)
  peptideregion<-lapply(unique_region, function(x,Peptide_Summary_file_regions){
    Peptide_Summary_file_regions[Peptide_Summary_file_regions$Region==x,]
  },Peptide_Summary_file_regions)
  
  names(peptideregion)<-unique_region
  
  Peptide_regions_Summary<-data.frame()
  
  #colnames(Peptide_regions_Summary)<-c("mz","Peptide","adduct","formula","isdecoy","moleculeNames","Score")
  #assign(paste0(x,y), data.frame(matrix(nrow = 0, ncol = length(col_names))))
  
  for (region_name in names(peptideregion)){
  
  peptideregion[[region_name]][,paste("Intensity",region_name)]=peptideregion[[region_name]]$Intensity
  
  peptideregion[[region_name]]$Intensity<-NULL
  peptideregion[[region_name]]$Score<-NULL  
  peptideregion[[region_name]]$Region<-NULL
  peptideregion[[region_name]]$Rank<-NULL
  if (nrow(Peptide_regions_Summary)==0){
    Peptide_regions_Summary=peptideregion[[region_name]]
  }else{
   Peptide_regions_Summary<-merge(Peptide_regions_Summary,peptideregion[[region_name]],by=c("mz","Peptide","adduct","formula","isdecoy","moleculeNames"),all=T) 
    
  }
    
  }
  
  a=colnames(Peptide_regions_Summary) %in% c(as.character(paste("Intensity",as.character(unique_region))))
  region_intensity=Peptide_regions_Summary[,a]
  region_intensity[is.na(region_intensity)]=0
  if (!is.null(dim(region_intensity))){
    
    Peptide_regions_Summary$Intensity=rowSums(region_intensity)}else{
      
      Peptide_regions_Summary$Intensity=region_intensity
    
  }
  #region_intensity[1,1]=NaN
  return(Peptide_regions_Summary)
}

MassSpecWavelet_fun<-function(peaklist,SNR.Th=3){
  library(MassSpecWavelet)
  #data("exampleMS")
  exampleMS=peaklist
  sortpeak=peaklist
  sortpeak=sortpeak[order(peaklist$intensities,decreasing = T),]
  sigthreshold=sortpeak$intensities[round(nrow(sortpeak)*0.1)]
  amp.Th=sigthreshold/max(sortpeak$intensities)
  
  rownames(exampleMS)=exampleMS$m.z
  exampleMS$m.z<-NULL
  exampleMS$intensities<-NULL
  exampleMS$V1=peaklist$intensities
  exampleMS<-as.matrix(exampleMS)
  SNR.Th <- 3
  peakInfo <- peakDetectionCWT(exampleMS, SNR.Th=SNR.Th,nearbyPeak =F, amp.Th = amp.Th,)
  majorPeakInfo = peakInfo$majorPeakInfo
  peakIndex <- majorPeakInfo$peakIndex
  plotPeak(exampleMS, peakIndex, main=paste('Identified peaks with SNR >', SNR.Th)) 
}

get_coord_info<-function(pixelinfo){
  library(stringr)
  #pixelinfo=names(pixels(imdata)[unlist(x[SPECTRUM_batch])])
  xinfo<-str_extract(pixelinfo,"x =.{0,4},")
  yinfo<-str_extract(pixelinfo,"y =.{0,4},")
  zinfo<-str_extract(pixelinfo,"z =.{0,4}")
  
  xinfo<-str_remove(xinfo,"x = ")
  xinfo<-str_remove(xinfo,",")
  xinfo<-as.numeric(str_remove(xinfo," "))
  
  yinfo<-str_remove(yinfo,"y = ")
  yinfo<-str_remove(yinfo,",")
  yinfo<-as.numeric(str_remove(yinfo," "))
  
  zinfo<-str_remove(zinfo,"z = ")
  zinfo<-str_remove(zinfo,",")
  zinfo<-as.numeric(str_remove(zinfo," "))
  
  return(list(x=xinfo,y=yinfo,z=zinfo))
  
}

spec_peakdetect<-function(x){
  library(MALDIrppa)
  data(spectra) # list of MassSpectra class objects
  data(type) # metadata
  # Summary of spectra features (results for 20 first mass spectra)
  summarySpectra(spectra[1:20])
  # Some pre-processing
  sc.results <- screenSpectra(spectra, meta = type)
  spectra <- sc.results$fspectra
  type <- sc.results$fmeta
  spectra <- transformIntensity(spectra, method = "sqrt")
  spectra <- wavSmoothing(spectra)
  spectra <- removeBaseline(spectra)
  names(spectra) <- type$SpectID # spectra IDs are lost with removeBaseline()
  # Summary of spectra features (results for positions 10 to 20)
  summarySpectra(spectra[10:20])
}

protein_scoring<-function(Protein_feature_list,Peptide_plot_list_rank){
  
  library(dplyr)
  
  
  
  Protein_feature_list_rank=merge(Protein_feature_list,Peptide_plot_list_rank,by=c("mz","Peptide","adduct","formula","isdecoy"))
  
  sum_pro_int<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Intensity=mean(Intensity)) 
  
  sum_pro_score<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Score=mean(Score)) 
  
  sum_pro_pep_count<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(peptide_count=length(unique(Peptide))) 
  
  #if (missing(sum_pro_pep_count_thero)){
  #  sum_pro_pep_count_thero<-Protein_feature_list %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(peptide_count_thero=length(unique(Peptide)))
  #}
  
  list_of_protein_sequence<-get("list_of_protein_sequence", envir = .GlobalEnv)
  Index_of_protein_sequence<-get("Index_of_protein_sequence", envir = .GlobalEnv)
  #Proteinlist=sum_pro_int$Protein
  sum_pro_coverage <- unlist(bplapply(sum_pro_int$Protein, function(x,Protein_feature_list_rank,list_of_protein_sequence,Index_of_protein_sequence,peptide_map_to_protein){
    #proteinid=Index_of_protein_sequence$desc[Index_of_protein_sequence$recno==x]
    #message(proteinid)
    protein_seq<-try(list_of_protein_sequence[x],silent = T)
    peptides<-unique(Protein_feature_list_rank$Peptide[Protein_feature_list_rank$Protein==x])
    if (class(protein_seq)!="try-error"){
     if (length(peptides)==1){
       width(peptides)/width(protein_seq)
     } else if (length(peptides)>=2){
       return(peptide_map_to_protein(peptides = peptides, protein_seq = protein_seq,map_type = "grep"))
     }
     
    }else{return(1)}
  },Protein_feature_list_rank=Protein_feature_list_rank,list_of_protein_sequence=list_of_protein_sequence,Index_of_protein_sequence=Index_of_protein_sequence,peptide_map_to_protein=peptide_map_to_protein,BPPARAM = BPPARAM))
  
  
  
  #sum_pro_pep_count<-merge(sum_pro_pep_count,sum_pro_pep_count_thero,by=c("Protein","isdecoy"),all.x=T)
  sum_pro_pep_count$peptide_coverage<-sum_pro_coverage
  
  #sum_pro_pep_count$peptide_coverage=sum_pro_pep_count$peptide_count/sum_pro_pep_count$peptide_count_thero
  
  Protein_feature_result<-merge(sum_pro_int,sum_pro_score,by=c("Protein","isdecoy"))
  
  Protein_feature_result<-merge(Protein_feature_result,sum_pro_pep_count,by=c("Protein","isdecoy"))
  
  Protein_feature_result$Intensity_norm<-Protein_feature_result$Intensity/median(Protein_feature_result$Intensity)
  
  Protein_feature_result$Proscore=(log(Protein_feature_result$Intensity_norm,10) + Protein_feature_result$Score) * Protein_feature_result$peptide_coverage
  
  Protein_feature_result=merge(Protein_feature_result,Index_of_protein_sequence[,c("recno","desc")],by.x="Protein",by.y="recno",all.x=T)
  
  return(Protein_feature_result)
}

peptide_map_to_protein<-function(peptides, protein_seq ,map_type="grep"){
  #peptides<-Protein_feature_list_rank[Protein_feature_list_rank$Protein=="sp|A2ASS6|TITIN_MOUSE Titin OS=Mus musculus OX=10090 GN=Ttn PE=1 SV=1",]
  #protein_seq=list_of_protein_sequence["sp|A2ASS6|TITIN_MOUSE Titin OS=Mus musculus OX=10090 GN=Ttn PE=1 SV=1"]
  #protein_seq=seq(protein_seq)
  #if(map_type c("grep","alignment")){map_type="grep"}
  suppressMessages(library(IRanges))
  if (map_type=="alignment"){
      suppressMessages(library(Biostrings))
  
  s1=AAStringSet(peptides)
  
  s2=protein_seq
  palign1 <- pairwiseAlignment(s1, s2,type="overlap")
  coverage_ranges<-palign1@subject@range
  cov <- coverage(coverage_ranges)
  
  #plotRanges(cov)
  cov <- as.vector(cov)
  cov_len<-length(which(cov==1))
  
  coverage_percentage=cov_len/width(protein_seq)
  }else if (map_type=="grep"){
    s1=(peptides)
    s2=as.character(protein_seq)
    palign2 <- sapply(s1,regexpr , s2)
    coverage_ranges<-IRanges(start = palign2,end = palign2+width(s1)-1,width = width(s1))
    coverage_ranges<-coverage_ranges[coverage_ranges@start>0,]
    cov <- coverage(coverage_ranges)
    
    #plotRanges(cov)
    cov <- as.vector(cov)
    cov_len<-length(which(cov!=0))
    
    coverage_percentage=cov_len/width(protein_seq)
    
  }

  return(coverage_percentage)
  
}

seq.cov <- function(x,plot_coverage=F){
  suppressMessages(library(IRanges))
  if(!class(x)=="PairwiseAlignmentsSingleSubject"){
    message('Only supports objects from class PairwiseAlignedFixedSubject')
    stop()
  }
  
  coverage_ranges<-x@subject@range
  #coverage_table<-NULL
  #coverage_table$start<-coverage_ranges@start
  #coverage_table$end<-coverage_ranges@start+coverage_ranges@width-1
  #coverage_table$width<-coverage_ranges@width
  #coverage_table<-data.frame(coverage_table,stringsAsFactors = F)
  cov <- coverage(coverage_ranges)
  
  #plotRanges(cov)
  cov <- as.vector(cov)
  cov_len<-length(which(cov==1))
  #mat <- cbind(seq_along(cov)-0.5, cov)
  #d <- diff(cov) != 0
  #mat <- rbind(cbind(mat[d,1]+1, mat[d,2]), mat)
  #mat <- mat[order(mat[,1]),]
  #lines(mat, col="red", lwd=4)
  #axis(2)
  
  
  return(cov_len)
}

isopattern_ppm_filter_peaklist<-function(pattern,ppm){
  
  pattern_ppm=as.numeric(as.character(pattern[,1]))
  
  pattern_ppm_delta=numeric()
  
  
  filtered_pattern<-pattern[1,]
  #filtered_pattern<-t(data.frame(filtered_pattern,stringsAsFactors = F))
  #dim(filtered_pattern)
  for (i in 1:(length(pattern_ppm)-1)){
    
    pattern_ppm_delta[i]=(pattern_ppm[i+1]-pattern_ppm[i])/pattern_ppm[i]
    
    if (pattern_ppm_delta[i]>(ppm/1000000)){
      filtered_pattern<-rbind(filtered_pattern,pattern[i+1,])
    } else {
      #previous_iso=filtered_pattern[nrow(filtered_pattern),]
      
      #newline=as.data.frame(list("m/z"=(filtered_pattern[nrow(filtered_pattern),1]*filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,1]*pattern[i+1,2])/(sum(filtered_pattern[nrow(filtered_pattern),2],pattern[i+1,2])),"abundance"=filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,2]))
      #names(newline)=c("m/z","abundance")
      newline=c((filtered_pattern[nrow(filtered_pattern),1]*filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,1]*pattern[i+1,2])/(sum(filtered_pattern[nrow(filtered_pattern),2],pattern[i+1,2])),
                filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,2])
      filtered_pattern[nrow(filtered_pattern),]<-newline
    }
    
  }
  rownames(filtered_pattern)=1:nrow(filtered_pattern)
  
  return(filtered_pattern)
}