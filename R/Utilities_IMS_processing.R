

#' imaging_Spatial_Quant
#'
#' This is a spatial quantitation function for maldi imaging data set
#' this function will read the candidate list file and generate quantification result
#' @param datafile specify the imzML data files
#' @param threshold specify the intensities threshold (0 to 1 in percentage)to report a identified molecule
#' @param ppm the mz tolerance (in ppm) for peak integration
#' @param Quant_list the quantifiaction candidate list, spatial quantification will go through every datafile and collect the ion intensities for each listed component
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param cal.mz If set with \code{"true"}, the function will recalculate the mz value according to the column named "formular" in the \code{Quant_list} and the specified adducts.
#' @param mzlist_bypass  Set \code{"true"} if you want to bypass the mzlist generating process
#' @param Protein_feature_summary  \code{"IMS_analysis"} follow-up process that will collect all the identified peptide information and associate them with possible proteins
#' @param plot_cluster_image  \code{"Protein_feature_summary"} follow-up process that will plot the protein cluster image
#'
#' @param Peptide_feature_summarya \code{"IMS_analysis"} follow-up process that will summarize all datafiles identified peptides and generats a \code{"peptide shortlist"} in the result summary folder
#' @param plot_ion_image  \code{"Peptide_feature_summarya"} follow-up process that will plot every connponents in the \code{"peptide shortlist"}
#' @param parallel the number of threads will be used in the PMF search, this option now only works for windows OS
#' @param spectra_segments_per_file optimal number of distinctive regions in the imaging, a virtual segmentation will be applied to the image files with this value. To have a better PMF result you may set a value that in the sweet point of sensitivety and false discovery rate (FDR).
#' @param Segmentation set as "spatialKMeans" to enable a \code{"spatialKMeans"} Segmentation; set as "spatialShrunkenCentroids" to enable a \code{"spatialShrunkenCentroids"} Segmentation; If a region rank file was supplied, you can set this as "Virtual_segmentation" to perform a manual segmentation; Set it as "none" to bypass the segmentation.
#' @param Smooth_range \code{"Segmentation"} pixel smooth range
#' @param Virtual_segmentation_rankfile specify a region rank file contains region information for manualy region segmentation
#' @return None
#'
#' @examples
#' imaging_Spatial_Quant(threshold=0.05, ppm=5,Digestion_site="[G]",
#'                        missedCleavages=0:1,Fastadatabase="murine_matrisome.fasta",
#'                        adducts=c("M+H","M+NH4","M+Na"),IMS_analysis=TRUE,
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
  ppm=2.5,
  Quant_list="Metabolites of Interest.csv",
  adducts=c("M-H","M+Cl"),
  cal.mz=F,
  mzlist_bypass=T,
  #==============TRUE if you want to plot protein PMF result
  IMS_analysis=TRUE,
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
  Smooth_range=1,
  Segmentation=c("spatialKMeans","spatialShrunkenCentroids","Virtual_segmentation","none"),
  Virtual_segmentation_rankfile=tk_choose.files(default = "Z:/George skyline results/maldiimaging/Maldi_imaging - Copy/radius_rank.csv",caption  = "Choose Virtual segmentation rank info file"),
  Spectrum_feature_summary=T,
  Region_feature_summary=T,
  Region_feature_analysis=T,
  plot_each_metabolites=T,
  Cluster_level="High",
  ClusterID_colname="Name",
  Region_feature_analysis_bar_plot=T,
  norm_datafiles=T,
  norm_Type="Median",
  mzrange="auto-detect",
  BPPARAM=bpparam(),
  Rotate_IMG=NULL,
  ...
){
  library(pacman)
  suppressMessages(suppressWarnings(p_load(RColorBrewer,RCurl,bitops,magick,ggplot2,reticulate,dplyr,stringr,tcltk,data.table,doParallel,
                                           iterators,foreach,protViz,cleaver,MALDIquant,Biostrings,XVector,IRanges,Cardinal,Rdisop,
                                           ProtGenerics,S4Vectors,stats4,EBImage,
                                           BiocParallel,BiocGenerics,parallel,stats,graphics,grDevices,utils,datasets,methods)))

  datafile<-gsub(".imzML", "", datafile)
  workdir<-base::dirname(datafile[1])
  setwd(workdir)
  closeAllConnections()
  parallel=try(detectCores()/2)
  if (parallel<1 | is.null(parallel)){parallel=1}
  BPPARAM=bpparam()
  BiocParallel::bpworkers(BPPARAM)=parallel
  bpprogressbar(BPPARAM)=TRUE

  if (is.null(Rotate_IMG)){
    Rotate_IMG$filenames=datafile;Rotate_IMG$rotation=rep(0,length(datafile))
    Rotate_IMG=as.data.frame(Rotate_IMG)
    write.csv(Rotate_IMG,"image_rotation.csv",row.names = F)
    Rotate_IMG=paste0(workdir,"/image_rotation.csv")
  }

  message(paste(try(detectCores()), "Cores detected,",parallel, "threads will be used for computing"))

  message(paste(length(datafile), "files were selected and will be used for Searching"))

  Meta_feature_list<-Meta_feature_list_fun(workdir=workdir,
                                           database = Quant_list,
                                           adducts=adducts,
                                           cal.mz = cal.mz,
                                           bypass=mzlist_bypass)

  if(IMS_analysis){

    #Peptide_Summary_searchlist<-unique(Protein_feature_list[,c("Peptide","mz","adduct")])

    Peptide_Summary_file<-IMS_data_process_quant(datafile,
                                                 Meta_feature_list,
                                                 segmentation_num=spectra_segments_per_file,
                                                 threshold=threshold,
                                                 ppm=ppm,
                                                 Segmentation=Segmentation,
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

      standardcol=c("Peptide","Name","FA","formula","mz","moleculeNames","mass","Intensity","adduct")
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

      Spectrum_summary_norm=Spectrum_summary
      clusterID=as.character(Spectrum_summary$ID)
      clusterID=data.frame(str_split(clusterID,"@"),stringsAsFactors = F)
      clusterID=as.data.frame(t(clusterID))
      clusterID=as.character(clusterID[["V1"]])


      par_norm_datacol_cluster<-function(col,Spectrum_summary,clusterID,norm_method=c("TIC","mTIC")){
        Spectrum_summary_col_norm<-Spectrum_summary[,col]
        Spectrum_summary_col_norm_value<-Spectrum_summary[,col]
        Spectrum_summary_col<-Spectrum_summary[,col]

        if( norm_method=="mTIC"){
          for (uniqueclusterID in unique(clusterID)){
            Spectrum_summary_col_norm_value[clusterID==uniqueclusterID]= median.default(Spectrum_summary_col[clusterID==uniqueclusterID],na.rm = T)
          }
        }else if (norm_method=="TIC"){
          Spectrum_summary_col_norm_value=median(Spectrum_summary_col[],na.rm = T)
        }

        Spectrum_summary_col_norm=Spectrum_summary_col/Spectrum_summary_col_norm_value

        Spectrum_summary_col_norm
      }

      message("Region_feature_analysis normalizaation")
      #Spectrum_summary_norm=unlist(parLapply(cl=autoStopCluster(makeCluster(detectCores())),which(colnames(Spectrum_summary)!="ID"),par_norm_datacol_cluster,Spectrum_summary,clusterID))
      Spectrum_summary_norm=unlist(bplapply(which(colnames(Spectrum_summary)!="ID"),par_norm_datacol_cluster,Spectrum_summary,clusterID,norm_method="TIC",BPPARAM = BPPARAM))
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

    colnames_case=gsub(".hr","",colnames_case)
    colnames_case=gsub("20.ctrl","-1",colnames_case)
    colnames_case=c("ID",as.numeric(colnames_case[2:length(colnames_case)]))
    #colnames_caset=gsub("y[:print:]+","",colnames_case)
    colnames(Spectrum_summary)=colnames_case
    colnames_case=c(colnames_case[1],sort(as.numeric(colnames_case[2:length(colnames_case)])))
    Spectrum_summary=Spectrum_summary[,colnames_case]
    Spectrum_summary_tran=t(Spectrum_summary)

    Spectrum_summary_tran=data.frame(Spectrum_summary_tran)
    colnames(Spectrum_summary_tran)=Spectrum_summary$ID

    colnames_case=as.numeric(colnames_case)
    colnames_case[1]=""
    Spectrum_summary_tran[,"Class"]=(colnames_case)



    case_info=str_split(Spectrum_summary$ID,"@")[1:length(Spectrum_summary$ID)]
    case_info=as.data.frame(case_info)
    case_info=t(case_info)
    colnames(case_info)=c("ID","adducts","mz","Rank")

    case_info=merge(case_info,radius_rank[,c("Rank","Name")],by="Rank",sort=F)

    colnames(case_info)=c("RegionRank","moleculeNames","adducts","mz","RegionName")
    case_info=merge(case_info,unique(Meta_feature_list[,c("moleculeNames",ClusterID_colname)]),all.x=T,all.y=F,by="moleculeNames")

    case_info$ClusterID=case_info[,`ClusterID_colname`]
    if (Cluster_level=="High"){
      a=data.table(str_split(case_info$ClusterID," ",simplify = T))
      case_info$ClusterID=a$V1
    }

    colnames_case_info=paste0(case_info$moleculeNames,"@",case_info$adducts,"@",case_info$mz,"@",case_info$RegionRank)

    case_info=t(case_info)
    case_info=data.frame(case_info)
    colnames(case_info)=colnames_case_info
    case_info$Class=""
    Spectrum_summary_tran=rbind(case_info,Spectrum_summary_tran)

    library(dplyr)
    library(plotly)
    library(stringr)

    moleculeNames=c(colnames_case_info)


    data=Spectrum_summary_tran[which(rownames(Spectrum_summary_tran)=="ID"):nrow(Spectrum_summary_tran),]
    data=Spectrum_summary_tran

    base::Sys.setenv("plotly_api_key"="FQ8kYs3JqmKLqKd0wGRv")
    base::Sys.setenv("plotly_username"="guoguodigital")
    base::Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1IjoiZ3VvZ3VvZGlnaXQiLCJhIjoiY2p1aHhheHM5MTBuYjQ0bnZzMzg0Mjd3aSJ9.XyqEayJi68xfGloNQQ28KA')


    plotly_for_region_with_ID<-function(i,data,moleculeNames,output_statice=T){
      library(dplyr)
      library(plotly)
      library(stringr)
      if (!require("processx")) install.packages("processx")
      windows_filename<- function(stringX){
        stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
        stringX<-gsub("\"", "", stringX)
        return(stringX)}


      x <- as.numeric(as.character(as.numeric(data$Class[which(data$Class!="")])))
      y <- as.numeric(as.character(data[which(data$Class!=""),moleculeNames[i]]))
      #plot(x,y)

      fit3 <- lm( y~poly(as.numeric(as.character(x)),3) )

      fitx=seq(min(x), max(x), length=1000)

      fitdata=as.data.frame(stats::predict(fit3, data.frame(x = fitx), interval="confidence"))
      newdata=data.frame(x=x,y=y)

      p <- plot_ly(data=newdata, x =~x , y = ~y, name = moleculeNames[i], type = 'scatter', mode = 'markers') %>%
        add_lines(x=unique(fitx),y=unique(fitdata$fit), name = 'Poly Fit', mode = 'lines',inherit=FALSE) %>%
        add_lines(x=unique(fitx),y=unique(fitdata$lwr), name = 'Poly Fit lwr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>%
        add_lines(x=unique(fitx),y=unique(fitdata$upr), name = 'Poly Fit upr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>%

        layout(title = moleculeNames[i],
               xaxis = list(title = "Age"),
               yaxis = list (title = "Conc."),
               showlegend = TRUE)

      p <- plot_ly(data=data, x =as.numeric(data$Class[which(data$Class!="")]) , y = data[which(data$Class!=""),moleculeNames[i]], name = data["moleculeNames",moleculeNames[i]], type = 'scatter', mode = 'markers') %>%
        add_trace(y =lowess(as.numeric(data$Class[which(data$Class!="")]),as.numeric(as.character(data[[moleculeNames[i]]][which(data$Class!="")])),iter=3)$y, name = 'Moving average', mode = 'lines') %>%

        layout(title = paste(data["moleculeNames",data["moleculeNames",moleculeNames[i]]], data["adducts",moleculeNames[i]]," in ",data["RegionName",moleculeNames[i]]),
               xaxis = list(title = "Cases"),
               yaxis = list (title = "Relative Conc."),
               showlegend = FALSE)


      if (output_statice){

        orca(p, file = windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".png")))
        windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".png"))
      }else{
        htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".html")), selfcontained = F, libdir = "lib")
        windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".html"))

      }
      windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".html"))
    }

    if(plot_each_metabolites){
      message("plot_each_metabolites")
      p_names=bplapply(which(moleculeNames!=""),plotly_for_region_with_ID,data,moleculeNames,BPPARAM=BPPARAM)

      zip("Result_lipid.zip", c(as.character(p_names), "lib"))
    }

    plotly_for_SUM_region_with_ID<-function(i,data,moleculeNames,output_statice=T){
      library(dplyr)
      library(plotly)
      library(stringr)
      if (!require("processx")) install.packages("processx")
      windows_filename<- function(stringX){
        stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
        stringX<-gsub("\"", "", stringX)
        return(stringX)}

      candidatelist<-data[,data["moleculeNames",]==moleculeNames[i]]
      region_list<-as.vector(t(candidatelist)[,"RegionName"])
      candidatelist<-cbind(candidatelist,Class=data[,"Class"])
      candidatelist_p=candidatelist[candidatelist$Class=="",]
      candidatelist_e=candidatelist[candidatelist$Class!="",]

      p <- plot_ly(data=candidatelist, x =~x , y = ~y, name = moleculeNames[i], type = 'scatter', mode = 'markers')

      x <- as.numeric(as.character(as.numeric(data$Class[which(data$Class!="")])))
      y <- as.numeric(as.character(data[which(data$Class!=""),moleculeNames[i]]))

      fit3 <- lm( y~poly(as.numeric(as.character(x)),3) )

      fitx=seq(min(x), max(x), length=1000)

      fitdata=as.data.frame(stats::predict(fit3, data.frame(x = fitx), interval="confidence"))
      newdata=data.frame(x=x,y=y)

      p <-  p %>%
        add_lines(x=unique(fitx),y=unique(fitdata$fit), name = 'Poly Fit', mode = 'lines',inherit=FALSE) %>%
        add_lines(x=unique(fitx),y=unique(fitdata$lwr), name = 'Poly Fit lwr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>%
        add_lines(x=unique(fitx),y=unique(fitdata$upr), name = 'Poly Fit upr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>%
        layout(title = moleculeNames[i],
               xaxis = list(title = "Age"),
               yaxis = list (title = "Conc."),
               showlegend = TRUE)

      p <- plot_ly(data=data, x =as.numeric(data$Class[which(data$Class!="")]) , y = data[which(data$Class!=""),moleculeNames[i]], name = data["moleculeNames",moleculeNames[i]], type = 'scatter', mode = 'markers') %>%
        add_trace(y =lowess(as.numeric(data$Class[which(data$Class!="")]),as.numeric(as.character(data[[moleculeNames[i]]][which(data$Class!="")])),iter=3)$y, name = 'Moving average', mode = 'lines') %>%

        layout(title = paste(data["moleculeNames",data["moleculeNames",moleculeNames[i]]], data["adducts",moleculeNames[i]]," in ",data["RegionName",moleculeNames[i]]),
               xaxis = list(title = "Cases"),
               yaxis = list (title = "Relative Conc."),
               showlegend = FALSE)


      if (output_statice){

        orca(p, file = windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".png")))
        windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".png"))
      }else{
        htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".html")), selfcontained = F, libdir = "lib")
        windows_filename(paste0(data["moleculeNames",moleculeNames[i]],".html"))

      }
      #plotly::orca(p, windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".png")))
      windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".html"))
    }

    moleculeNames=unique(as.vector(t(data)[,"moleculeNames"]))
    moleculeNames=moleculeNames[moleculeNames!=""]

    if(plot_each_metabolites){
      message("plot_each_metabolites")
      #p_names=parLapply(cl=autoStopCluster(makeCluster(4)),which(moleculeNames!=""),plotly_for_region_with_ID,data,moleculeNames)
      p_names=bplapply(which(moleculeNames!=""),plotly_for_SUM_region_with_ID,data,moleculeNames,BPPARAM=BPPARAM)
      zip("Result_lipid.zip", c(as.character(p_names), "lib"))
    }

    Spectrum_summary_tran_tran=as.data.frame(t(Spectrum_summary_tran),stringsAsFactors = FALSE)
    Region_CLuster=paste(Spectrum_summary_tran_tran$ClusterID,Spectrum_summary_tran_tran$RegionName)
    Spectrum_summary_tran_tran=Spectrum_summary_tran_tran[,which(Spectrum_summary_tran_tran["Class",]!="")]
    Spectrum_summary_tran_tran=Spectrum_summary_tran_tran[which(rownames(Spectrum_summary_tran_tran)!="Class"),]
    Spectrum_summary_tran_tran <- mutate_all(Spectrum_summary_tran_tran, function(x) as.numeric(as.character(x)))
    DT <- as.data.table(Spectrum_summary_tran_tran[])
    DT$ClusterID=Region_CLuster[which(Region_CLuster!=" ")]
    Spectrum_summary_tran_tran=DT[ , lapply(.SD, sum), by = "ClusterID"]
    rownames(Spectrum_summary_tran_tran)=Spectrum_summary_tran_tran$ClusterID
    data=as.data.frame(t(Spectrum_summary_tran_tran))
    colnames(data)=t(data["ClusterID",])
    data=data[which(rownames(data)!="ClusterID"),]
    data$Class=Spectrum_summary_tran$Class[which('&'(Spectrum_summary_tran$Class!="",Spectrum_summary_tran$Class!=" "))]


    Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/Anaconda/orca_app", sep = .Platform$path.sep))


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
      x <- as.factor(as.character(data$Class))
      y <- as.numeric(as.character(data[,moleculeNames[i]]))

      fit3 <- lm( y~poly(as.numeric(as.character(x)),3) )

      fitx=seq(min(as.numeric(data$Class)), max(as.numeric(data$Class)), length=1000)

      fitdata=as.data.frame(stats::predict(fit3, data.frame(x = fitx), interval="confidence"))
      newdata=data.frame(x=x,y=y)
      p <- plot_ly(data=newdata, x =data$Class , y = ~y, name = moleculeNames[i], type = 'scatter', mode = 'markers') %>%
        add_lines(x=unique(fitx),y=unique(fitdata$fit), name = 'Poly Fit', mode = 'lines',inherit=FALSE) %>%
        add_lines(x=unique(fitx),y=unique(fitdata$lwr), name = 'Poly Fit lwr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>%
        add_lines(x=unique(fitx),y=unique(fitdata$upr), name = 'Poly Fit upr', mode = 'lines',inherit=FALSE, line = list(color="grey", dash = 'dash')) %>%
        layout(title = moleculeNames[i],
               xaxis = list(title = "Age"),
               yaxis = list (title = "Conc."),
               showlegend = TRUE)

      if (output_statice){

        orca(p, file = windows_filename(paste0(moleculeNames[i],".png")))
        windows_filename(paste0(moleculeNames[i],".png"))
      }else{
        htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(moleculeNames[i],".html")), selfcontained = F, libdir = "lib")
        windows_filename(paste0(moleculeNames[i],".html"))

      }


    }
    message("plotly_for_region_with_ClusterID")
    p_names=bplapply(which(colnames(data)!="Class"),plotly_for_region_with_ClusterID,data,output_statice=T,BPPARAM = BPPARAM)
    zip("Result_clusterID.zip", c(as.character(p_names), "lib"))

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

        plotly_IMAGE(p, format = "png", out_file = windows_filename(paste0(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"],".png")))
        htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(case_info[i,"ID"], case_info[i,"adducts"],"in",case_info[i,"Name"],".html")))
      }
      listreturn=parallel::parSapply(cl=makeCluster(4),2:ncol(Spectrum_summary_norm_sum),plotly_for_region,Spectrum_summary_norm_sum,case_info)
      # Create a shareable link to your chart
      # Set up API credentials: https://plot.ly/r/getting-started
      chart_link = api_create(p, filename="line-mode1")
      chart_link}
    write.csv(Spectrum_summary,file = paste(workdir,"/Summary folder/Region_feature_summary_result.csv",sep=""),row.names = F)
  }
  if(Region_feature_analysis_bar_plot){
    options(scipen=999)
    setworkdir(paste(workdir,"/Summary folder/Region_feature_analysis/",sep=""))
    library(plotly)
    if (is.null(Spectrum_summary)){Spectrum_summary=read.csv(file = paste(workdir,"/Summary folder/Region_feature_summary.csv",sep=""))}
    radius_rank=read.csv(file =paste0(workdir[1],"/",Virtual_segmentation_rankfile))
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

    case_info$ClusterID=case_info$moleculeNames
    case_info$ClusterID=case_info$ClusterID
    if (Cluster_level=="High"){
      a=data.table(str_split(case_info$ClusterID," ",simplify = T))
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

    #Protein_feature_list=merge(Protein_feature_list,Meta_feature_list[,c("moleculeNames","Pathway")],by="moleculeNames",allow.cartesian=TRUE)
    #Protein_feature_list$Pathway=Protein_feature_list$Pathway.y
    if (!is.null(Rotate_IMG)){Rotate_IMG=read.csv(Rotate_IMG,stringsAsFactors = F)}
    for (i in 1:length(datafile)){
      rotate=Rotate_IMG[Rotate_IMG$filenames==datafile[i],"rotation"]
      if(is.null(rotate)) rotate=0
      rotate= as.numeric(rotate)
      imdata=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm,rotate = rotate,attach.only=F)


      currentdir<-paste0(datafile[i] ," ID")
      if (dir.exists(paste(currentdir,"/cluster Ion images",sep=""))==FALSE){dir.create(paste(currentdir,"/cluster Ion images",sep=""))}
      setwd(paste(currentdir,"/cluster Ion images",sep=""))

      lapply(unique(Protein_feature_list$Name),
             cluster_image_cardinal_allinone,
             imdata=imdata,
             SMPLIST=Protein_feature_list,

             ppm=ppm,ClusterID_colname=ClusterID_colname,
             componentID_colname="moleculeNames",
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

           ppm=ppm,ClusterID_colname=ClusterID_colname,
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

    Protein_feature_list$Pathway=NULL
    Protein_feature_list=merge(Protein_feature_list,Meta_feature_list[,c("moleculeNames","Pathway")],by="moleculeNames",allow.cartesian=TRUE)
    if (!is.null(Rotate_IMG)){Rotate_IMG=read.csv(Rotate_IMG,stringsAsFactors = F)}
    imdata=list()
    combinedimdata=NULL
    #register(SerialParam())
    if (is.null(mzrange)){
      message("Detecting mz range...")
      for (i in 1:length(datafile)){
        rotate=Rotate_IMG[Rotate_IMG$filenames==datafile_imzML[i],"rotation"]
        if(is.null(rotate)) rotate=0
        rotate= as.numeric(rotate)
        imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm*2,rotate = rotate,attach.only=F,as="MSImagingExperiment")
        mzrangetemp=range(imdata[[1]]@featureData@mz)
        if (is.null(mzrange)) {mzrange=mzrangetemp}else{
          if (mzrange[1]<mzrangetemp[1]){mzrange[1]<-mzrangetemp[1]}
          if (mzrange[2]<mzrangetemp[2]){mzrange[2]<-mzrangetemp[2]}
        }

      }
    }
    for (i in 1:length(datafile)){

      rotate=Rotate_IMG[Rotate_IMG$filenames==datafile_imzML[i],"rotation"]
      if(is.null(rotate)) rotate=0
      rotate=as.numeric(rotate)
      imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm*2,rotate = rotate,attach.only=F,as="MSImagingExperiment",mzrange=mzrange)
      #imdata[[i]]@elementMetadata@coord=imdata[[i]]@elementMetadata@coord[,c("x","y")]
      if (i==1) {
        combinedimdata<-imdata[[i]]
      }else{
        combinedimdata<-cbind(combinedimdata,imdata[[i]])
      }
      imdata[[i]]=NULL

    }


    combinedimdata@elementMetadata@coord@listData[["z"]]=NULL

    imdata=combinedimdata

    outputfolder=paste(workdir,"/Summary folder/cluster Ion images/",sep="")

    if (dir.exists(outputfolder)==FALSE){dir.create(outputfolder)}
    setwd(outputfolder)

    Protein_feature_list=as.data.frame(Protein_feature_list)
    if(is.null(Protein_feature_list$formula)){Protein_feature_list$formula=Protein_feature_list$formular}

    clusterID=unique(Protein_feature_list[,ClusterID_colname])

    lapply(clusterID,
           cluster_image_grid,
           imdata=imdata,
           SMPLIST=Protein_feature_list,
           ppm=ppm,
           ClusterID_colname=ClusterID_colname,
           Component_plot_threshold=2,
           componentID_colname="moleculeNames",
           plot_layout="line",
           export_Header_table=F,
           plot_style="fleximaging")

    lapply(clusterID,
           cluster_image_grid,
           imdata=NULL,
           SMPLIST=Protein_feature_list,
           ppm=ppm,
           ClusterID_colname=ClusterID_colname,
           Component_plot_threshold=2,
           componentID_colname="moleculeNames",
           plot_layout="line",
           export_Header_table=T,
           plot_style="fleximaging")




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




searchPMF<-function(pimlist,spectrumlist,ppm){
  pimresultlist<-pimlist
  print("Start PMF search")
  for (pim in 1:length(pimlist)){

    pimresultlist[[pim]]<-PMFsum(pimlist[[pim]],spectrumlist,ppm)


  }

  return(pimresultlist)
}



IMS_analysis_fun<-function(Peptide_Summary_searchlist,peaklist,ppm,BPPARAM =bpparam(),mzrange){
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

IMS_analysis_fun_2<-function(Peptide_Summary_searchlist,peaklist,ppm,BPPARAM =bpparam(),mzrange){
  if (missing(mzrange)){mzrange=c(min(peaklist$m.z),max(peaklist$m.z))}
  Peptide_Summary_searchlist$Intensity<-rep(0,nrow(Peptide_Summary_searchlist))
  message("PMF 1st search...")
  match_ress<-lapply(1:nrow(peaklist), function(x,peaklist,ppm,Peptide_Summary_searchlist){
    lowmz<-peaklist$m.z[x]-peaklist$m.z[x]*ppm/1000000
    highmz<-peaklist$m.z[x]+peaklist$m.z[x]*ppm/1000000
    #Peptide_Summary_searchlist$Intensity[which(data.table::between(Peptide_Summary_searchlist$mz, lowmz, highmz))]<<-peaklist$intensities[x]
    return(list(mz_which=which(data.table::between(Peptide_Summary_searchlist$mz, lowmz, highmz)),
                mz_Intensity=peaklist$intensities[x]))
    },peaklist,ppm,Peptide_Summary_searchlist)

  for (match_res in match_ress){
    Peptide_Summary_searchlist$Intensity[match_res$mz_which]<-rep(match_res$mz_Intensity,length(match_res$mz_which))
  }
  return(Peptide_Summary_searchlist)
}

virtual_segmentation<-function(imdata,Virtual_segmentation_rankfile="~/GitHub/HiTMaP/inst/data/radius_rank_bovin.csv"){
  library(reshape2)
  library(BiocParallel)
  coordatadf=imdata@elementMetadata@coord
  coordatadf$run<-imdata@elementMetadata@run
  coordatadf$indx<-as.numeric(rownames(coordatadf))
  coordatadf<-as.data.frame(coordatadf)
  coordatalist=list()
  for (i in unique(coordatadf$run)){

    coordatalist[[i]]<-coordatadf[coordatadf$run==i,]

  }


  radius_rank=read.csv(file = Virtual_segmentation_rankfile)
  radius_rank=radius_rank[order(radius_rank$Rank),]


  coordata_rank_list<- function(coordata,radius_rank,BPPARAM=bpparam()){
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

    rank_pixel<-function(coordata,coordistmatrix){
      shape_center=coordata[coordistmatrix$sum==min(coordistmatrix$sum),]
      center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),1:nrow(coordata)])
      library(useful)
      From <- shape_center[rep(seq_len(nrow(shape_center)), each=nrow(coordata)),1:2]
      To <- coordata[,1:2]
      df=To-From
      center_edge_angle<-cbind(coordata[,1:2],cart2pol(df$x, df$y, degrees = F),edge=coordata[,"edge"])
      center_edge_angle_sdge=center_edge_angle[center_edge_angle$edge==TRUE,]
      coordata$rank=0
      coordata$pattern=""

      for (i in 1: (nrow(coordata))){

        if (coordata$edge[i]!=TRUE){
          df=coordata[i,1:2]-shape_center[,1:2]
          point_center_angle<-cbind(coordata[i,1:2],cart2pol(df$x, df$y, degrees = F))
          pointedge=center_edge_angle_sdge[which(abs(center_edge_angle_sdge$theta-point_center_angle$theta)==min(abs(center_edge_angle_sdge$theta-point_center_angle$theta))),]

          pointedge=pointedge[which.min(pointedge$r),]
          to_edge=coordistmatrix[[i]]['&'(coordata$x==pointedge$x,coordata$y==pointedge$y)]
        }else{to_edge=0}


        to_center=center_dist[i]
        total=to_edge+to_center
        max(radius_rank$Radius_U)
        norm_center_dist=to_center/total*max(radius_rank$Radius_U)
        coordata$rank[i]=as.character(radius_rank$Rank['&'(radius_rank$Radius_L<norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
        coordata$pattern[i]=as.character(radius_rank$Name['&'(radius_rank$Radius_L<norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
      }
      coordata
    }

    coordata=rank_pixel(coordata,coordistmatrix)

    coordata

  }


  coordatalist_result<-lapply(coordatalist,coordata_rank_list,radius_rank)

  coordatalist_result_merge<-do.call(rbind,coordatalist_result)
  coordatadf_merge<-merge(coordatadf,coordatalist_result_merge[,c("x","y","run","rank","pattern")],by=c("x","y","run"))
  rownames(coordatadf_merge)<-coordatadf_merge$indx
  return(coordatadf_merge)

}




affine_grid <- function(x, translate=c(0,0), rotate=0,
                        angle=c("degrees", "radians"), grid=TRUE)
{
  x=x[,c("x","y")]
  angle <- match.arg(angle)
  theta <- -rotate
  if ( angle == "degrees" ) theta <- theta * pi / 180
  new.x=spdep::Rotation(x,angle = theta)
  #
  # translate center of mass to be near origin
  #tt <- sapply(x, function(xs) mean(xs))
  #new.x <- t(as.matrix(x)) - tt
  # rotate around origin
  #A <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2)
  #new.x <- A %*% new.x
  # translate back and do requested translation
  #new.x <- t(new.x + tt + translate)
  # remove negative coordinates and round to integers
  if ( grid ) {
    new.x <- round(new.x)
  }
  # return data.frame of new coordinates
  new.x <- as.data.frame(new.x)
  names(new.x) <- names(x)
  new.x
}

rotate_image<-function(imzdata,degree=0){

  imaxmldata<- readLines(paste0(imzdata,".imzML"))
  filename=paste0(imzdata,".imzML")

  doc = xmlRoot(xmlTreeParse(filename, trim = FALSE, ignoreBlanks = FALSE))
  print(doc, indent = FALSE, tagSeparator = "")


  rotatetmp<-imdata@pixelData@data

  rotatenew<-affine_grid(rotatetmp[,c("x","y")])
  rownames(rotatenew)<-paste0("x = ",rotatenew$x,", ","y = ",rotatenew$y,", ","z = ",rotatenew$z)
  rotatenew$z=rotatetmp$z
  rotatenew$sample=rotatetmp$sample
  timdatapositionArray<-data.frame(imdata@imageData@positionArray,stringsAsFactors = F)
  class(timdatapositionArray)
  new_timdatapositionArray=data.frame(row.names = 1:max(rotatenew$y))

  new_timdatapositionArray=1
  imdata@pixelData@data<-rotatenew
  imdata@imageData@positionArray<-ttimdatapositionArray
  imdata@imageData@coord<-rotatenew
  rownames(imdata@imageData@coord)<-paste0("x = ",rotatenew$x,", ","y = ",rotatenew$y,", ","z = ",rotatenew$z)
  imdata@imageData@dimnames[[2]]=paste0("x = ",rotatenew$x,", ","y = ",rotatenew$y,", ","z = ",rotatenew$z)

  image_rotate(image, degrees)
  plot(1:10, rnorm(10))

  suppressMessages(suppressWarnings(require(gridGraphics)))

  grab_grob <- function(){
    grid.echo()
    grid.grab()
  }

  g <- grab_grob()
  grid.newpage()
  pushViewport(viewport(width=0.7,angle=30))
  grid.draw(g)
}

rank_mz_feature<-function(Peptide_plot_list,mz_feature,BPPARAM=bpparam()){

  Peptide_plot_list<-as.data.frame(Peptide_plot_list)

  suppressMessages(suppressWarnings(require(data.table)))

  mz=unique(Peptide_plot_list$mz)

  Peptide_plot_list$mz_align<-NULL

  mz_align<-sapply(mz,function(x,mz_feature){

    round(mz_feature$m.z[which.min(abs(mz_feature$m.z-x))],digits = 4)

  },mz_feature)

  mz_feature_cluster<-data.table(mz=mz,mz_align=mz_align)

  Peptide_plot_list<-merge(Peptide_plot_list,mz_feature_cluster,by="mz",all=T)

  mz=unique(Peptide_plot_list$mz_align)

  message(paste("Ranking mz feature:",length(unique(Peptide_plot_list$mz)), "unique candidates mz,",length(unique(Peptide_plot_list$mz_align)),"aligned mz feature"))

  if(F){Peptide_plot_list_rank<-bplapply(mz,function(x,Peptide_plot_list){

    randklist <- Peptide_plot_list[Peptide_plot_list$mz_align==x,]

    randklist$Rank<- as.numeric(factor(randklist$Score,levels=sort(unique(randklist$Score),decreasing = T)))

    return(randklist)

  },Peptide_plot_list,BPPARAM=BPPARAM)}



  Peptide_plot_list_rank<-lapply(X=mz,FUN = function(x,Peptide_plot_list){
    #message(x)
    randklist <- Peptide_plot_list[Peptide_plot_list$mz_align==x,]
    ranking=as.numeric(factor(randklist$Score,levels=sort(unique(randklist$Score),decreasing = T)))
    randklist$Rank<- ranking
    #message(ranking)
    return(randklist)

  },Peptide_plot_list=Peptide_plot_list[,c("mz_align","Score")])

  if(F){Peptide_plot_list_rank_temp<-NULL
  randklist<-list()
  for (x in mz[1:10]){
    randklist[[x]] <- Peptide_plot_list[Peptide_plot_list$mz_align==x,]
    ranking=as.numeric(factor(randklist[[x]]$Score,levels=sort(unique(randklist[[x]]$Score),decreasing = T)))
    randklist[[x]]$Rank<- ranking
  }
  Peptide_plot_list_rank_temp<-do.call(rbind,randklist)
  #Peptide_plot_list
  identical(Peptide_plot_list_rank_temp,Peptide_plot_list_rank)}

  Peptide_plot_list_rank<-do.call(rbind,Peptide_plot_list_rank)

  Peptide_plot_list_rank<-merge(Peptide_plot_list_rank,Peptide_plot_list,by=c("mz_align","Score"))

  return(Peptide_plot_list_rank)

}

rotateMSI<-function(imdata,rotation_degree=0){

  rotatetmp<-as.data.frame(imdata@elementMetadata@coord@listData)

  rotatenew<-affine_grid(rotatetmp[,c("x","y")],rotate=rotation_degree,grid = T)

  #sum(duplicated(rotatenew))

  rotatenew$x<-rotatenew$x-min(rotatenew$x)+1

  rotatenew$y<-rotatenew$y-min(rotatenew$y)+1

  rownames(rotatenew)<-paste0("x = ",rotatenew$x,", ","y = ",rotatenew$y,", ","z = ",rotatenew$z)
  rotatenew$z=rotatetmp$z
  imdata@elementMetadata@coord@listData<-as.list(rotatenew)
  return(imdata)
}

PCA_ncomp_selection<-function(imdata,variance_coverage=0.80,outputdir=NULL){
  suppressMessages(suppressWarnings(library(Cardinal)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  suppressMessages(suppressWarnings(library(gtable)))
  suppressMessages(suppressWarnings(library(egg)))
  percent<-function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }

  PCA_imdata<-Cardinal::PCA(imdata,ncomp=12)
  # if (!is.null(outputdir)){
  #   saveRDS(PCA_imdata,paste0(outputdir,"PCA_imdata.rds"))
  #   saveRDS(as.data.frame(summary(PCA_imdata)),paste0(outputdir,"PCA_imdata_df.rds"))
  # }
  PCA_imdata_df<-data.frame(Component=1:length(PCA_imdata@resultData@listData[[1]][["sdev"]]) , Standard.deviation=PCA_imdata@resultData@listData[[1]][["sdev"]])
  PCA_imdata_df$Standard.deviation<-PCA_imdata_df$Standard.deviation/sum(PCA_imdata_df$Standard.deviation)
  PCA_imdata_df$Component<-as.factor(PCA_imdata_df$Component)
  PCA_imdata_df$Percentage<-percent(PCA_imdata_df$Standard.deviation,digits =1)
  PCA_imdata_df$Percent<-PCA_imdata_df$Standard.deviation/sum(PCA_imdata_df$Standard.deviation)
  PCA_imdata_df <- within(PCA_imdata_df, cumulate <- cumsum(Percent))
  for (PC_cum in PCA_imdata_df$cumulate){
    if (PC_cum<variance_coverage){
      PCA_imdata_df$Include[which(PCA_imdata_df$cumulate==PC_cum)]<-"Yes"
    }else if(`&`(PC_cum>variance_coverage,which(PCA_imdata_df$cumulate==PC_cum)==1)){
      PCA_imdata_df$Include[which(PCA_imdata_df$cumulate==PC_cum)]<-"Yes"
    }else if(`&`(PC_cum>variance_coverage,PCA_imdata_df$cumulate[which(PCA_imdata_df$cumulate==PC_cum)-1]<variance_coverage)){
      PCA_imdata_df$Include[which(PCA_imdata_df$cumulate==PC_cum)]<-"Yes"
    }else{
      PCA_imdata_df$Include[which(PCA_imdata_df$cumulate==PC_cum)]<-"No"
    }
  }

  ncomp<-max(which(PCA_imdata_df$Include=="Yes"))
  PCA_imdata_df$Include<-factor(PCA_imdata_df$Include,levels = c("Yes","No"))
  actual_coverage<-percent(PCA_imdata_df$cumulate[ncomp])
  if (!is.null(outputdir)){
    png(paste(outputdir,"/","PCA_image.png",sep=""),width = 1024,height = 720)

    print(image(PCA_imdata, values="scores", superpose=FALSE, layout=c(3,4),normalize.image = c("linear"),contrast.enhance = c("histogram")))

    dev.off()
    png(paste(outputdir,"/","PCA_plot.png",sep=""),width = 1024,height = 720 ,res=150)
    print(ggplot(PCA_imdata_df, aes(x = Component,y = Percent,label=Percentage,fill=Include)) +
            geom_histogram(aes(y = Standard.deviation),stat="identity",color="black",show.legend=T)+
            labs(title ="",x = "Component",y = "Variation coverage")+
            geom_label(
              colour = "Black",
              position = position_stack(vjust = 0.5),
              angle = 90,show.legend = FALSE)+
            scale_y_continuous(labels = scales::percent)+
            theme_article()+ theme(legend.position = c(0.85, 0.8))+ ggtitle(paste("Principle components variance Actual coverage:",actual_coverage))
    )
    dev.off()


  }

  return(ncomp)
}




Preprocessing_segmentation<-function(datafile,
                                     workdir=NULL,
                                     segmentation_num=5,
                                     ppm=5,import_ppm=0,Bypass_Segmentation=F,
                                     mzrange="auto-detect",
                                     Segmentation=c("spatialKMeans","spatialShrunkenCentroids","Virtual_segmentation","none","def_file"),
                                     Segmentation_def="segmentation_def.csv",
                                     Segmentation_ncomp="auto-detect",
                                     Segmentation_variance_coverage=0.8,
                                     Smooth_range=1,
                                     colorstyle="Set1",
                                     Virtual_segmentation_rankfile=NULL,
                                     rotate=NULL,rotate_img=F,
                                     BPPARAM=bpparam(),
                                     preprocess=list(force_preprocess=FALSE,use_preprocessRDS=TRUE,mz_bin_list=NULL,smoothSignal=list(method="gaussian"),
                                                     reduceBaseline=list(method="locmin"),
                                                     peakPick=list(method="Default",SNR=6,window=5),
                                                     peakAlign=list(tolerance=ppm/2, units="ppm", level=c("global","local")),
                                                     peakFilter=list(freq.min=0.05),
                                                     normalize=list(method=c("rms","tic","reference")[1],mz=1))){
  
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(Cardinal)))
  suppressMessages(suppressWarnings(require(RColorBrewer)))
  suppressMessages(suppressWarnings(require(stringr)))
  setCardinalBPPARAM(BPPARAM)
  
  datafile <- paste0(workdir,"/",datafile)
  workdir <- dirname(datafile)
  datafile <- basename(datafile)
  datafile <-str_remove(datafile,regex(".imzML$"))
  rotate = Parse_rotation(datafile,rotate)
  datafile_imzML<-paste0(datafile,".imzML")
  
  for (z in 1:length(datafile)){
    Segmentation_ncomp_running<-Segmentation_ncomp
    name <-basename(datafile[z])
    name <-gsub(".imzML$","",name)
    name <-gsub("/$","",name)
    folder<-base::dirname(datafile[z])
    
    if (!str_detect(datafile[z],".imzML$")){
      datafile_imzML[z]<-paste0(datafile[z],".imzML")
    }
    setwd(workdir[z])
    
    if(import_ppm==0)
      if (ppm>=25) {
        instrument_ppm=50
        import_ppm = floor(instrument_ppm*2/5)[1]
      }else{
        instrument_ppm=10
        import_ppm = instrument_ppm
      }
    
    #setup import ppm which ensure pickpicking has correct number of data points (halfwindow>=2) per peak to work with
    if (import_ppm > ppm ) import_ppm = ppm
    
    imdata_org<-NULL
    imdata<-NULL
    imdata_ed<-NULL
    
    if (dir.exists(paste0(gsub(".imzML$","",datafile[z]) ," ID"))==FALSE){
      dir.create(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
    }
    
    message("Loading raw image data for statistical analysis: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
    if(mzrange[1]=="auto-detect"){
      imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=F,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam())
      imdata@centroided<-F
    }else {
      imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=F,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam(),mass.range =mzrange)
      imdata <-imdata[between(mz(imdata),mzrange[1],mzrange[2]),]
      imdata@centroided<-F
    }
    
    if(rotate_img){
      
      if(!is.na(rotate[datafile_imzML[z]])){
        imdata_r <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
      }else if(!is.na(rotate[datafile[z]])){
        imdata_r <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
      }
      
      run(imdata)<-factor(rep("Before",length(run(imdata))))
      Cardinal::cbind(imdata,imdata_r)->imdata_rimdata_new
      imdata_rimdata_new@elementMetadata@coord@listData[["z"]]<-NULL
      png(paste0(gsub(".imzML$","",datafile[z])  ," ID/","IMG_roation.png"),width = 1600,height = 800, res=200)
      print(image(imdata_rimdata_new,factor(run(imdata_rimdata_new)) ~ x * y ,superpose=F, key=T,
                  layout=c(1,2)))
      dev.off()
      imdata_r->imdata
      rm(imdata_r)
    }
    
    
    
    if ('|'(!file.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS")),
            if(!is.null(preprocess[["force_preprocess"]])){
              preprocess$force_preprocess
            }else{T})) {
      
      message("Preparing image data for statistical analysis: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
      if (!is.null(preprocess)){
        imdata_ed<-imdata
        rm(imdata)

        setCardinalBPPARAM(BPPARAM)
        if  ( ppm<25){

          
          if (!is.null(preprocess$normalize)){
            if (preprocess$normalize$method=="Disable") {
            } else if (preprocess$normalize$method %in% c("rms","tic")){
              imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
            } else if ('&'(preprocess$normalize$method == "reference", !is.null(preprocess$normalize$mz))){
              norm_feature<-which(dplyr::between(imdata_ed@featureData@mz,
                                                 preprocess$normalize$mz*(1-ppm/1000000),
                                                 preprocess$normalize$mz*(1+ppm/1000000)))
              if (length(norm_feature)>=1){
                imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method, feature = norm_feature) %>% process(BPPARAM=SerialParam())
              }
            } else {
              imdata_ed<- imdata_ed %>% normalize(method="rms") %>% process()
            }
          }

          
          #peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
          #write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
          

          #smoothSignal(method="gaussian") %>%
          #reduceBaseline(method="locmin") %>%
          if (!is.null(preprocess$mz_bin_list)){
            preprocess$mz_bin_list->mz_bin_list_s
            sort(as.numeric(unlist(mz_bin_list_s)))->ref_mz
            peaklist_deco<-data.frame(mz=ref_mz,intensities=rep(1,length(ref_mz)))
            peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
            unlist(peaklist_deco[,1])->ref_mz
            ref_mz<-ref_mz[!is.na(ref_mz)]
            range(mz(imdata_ed))->ref_mz_range
            ref_mz<-ref_mz[ref_mz*(1-ppm/1000000)>=ref_mz_range[1]]
            ref_mz<-ref_mz[ref_mz*(1+ppm/1000000)<=ref_mz_range[2]]
            ref_mz<-unique(ref_mz)
            ref_mz<-ref_mz[!is.na(ref_mz)]
            imdata_ed<-imdata_ed %>% mzBin(ref=ref_mz, tolerance=ppm, units="ppm") %>% process()
            write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
            
          }else{
            
          if (preprocess$smoothSignal$method=="Disable") {
          }else if (!is.null(preprocess$smoothSignal$method)){
            imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
          }else{
            imdata_ed<- imdata_ed %>% smoothSignal(method="gaussian")
          }
          
          if (preprocess$reduceBaseline$method=="Disable") {
          }else if (!is.null(preprocess$reduceBaseline$method)){
            imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
          }else{
            imdata_ed<- imdata_ed %>% reduceBaseline(method="locmin")
          }
          
          if (preprocess$peakPick$method=="Disable") {
          }else if (preprocess$peakPick$method %in% c("adaptive","mad","simple")){
            if (is.null(preprocess$peakPick$SNR)) preprocess$peakPick$SNR=6
            if (is.null(preprocess$peakPick$window)) preprocess$peakPick$window=5
            imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method, SNR=preprocess$peakPick$SNR, window=preprocess$peakPick$window) %>% process()
          }else if (preprocess$peakPick$method == "Default|default"){
            #add an peak picking function other than the Cardinal options.
            peaklist<-summarizeFeatures(imdata_ed,FUN = "sum", as="DataFrame")
            peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
            peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
            peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
            write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
            imdata_ed<-imdata_ed %>% peakBin(peaklist_deco$mz, tolerance=ppm, units="ppm") %>% process()
          }
            
          }
          
          saveRDS(imdata_ed,paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_peakpicked_imdata.RDS"))
          
          if (is.null(preprocess$peakAlign$level)) preprocess$peakAlign$level<-"local"
          if (preprocess$peakAlign$level=="global"){
            if (preprocess$peakAlign$tolerance==0 ) {
              message("preprocess$peakAlign$tolerance set as zero, step bypassed")
            }else if ('&'(!is.null(preprocess$peakAlign$tolerance),!is.null(preprocess$peakAlign$units))){
              message("preprocess$peakAlign$tolerance set as ", preprocess$peakAlign$tolerance)
              imdata_ed<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units)
              peaklist<-summarizeFeatures(imdata_ed,FUN = "sum", as="DataFrame")
              peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
              peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
              write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec.csv"),row.names = F)
              peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
              write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
            }else {
              message("preprocess$peakAlign$tolerance missing, use default tolerance in ppm ", ppm/2)
              imdata_ed<- imdata_ed %>% peakAlign(tolerance=ppm/2, units="ppm")
              peaklist<-summarizeFeatures(imdata_ed,FUN = "sum", as="DataFrame")
              peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
              peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
              write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec.csv"),row.names = F)
              peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
              write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
            }
          }
          
          imdata_ed<- imdata_ed %>% process()
          #imdata_ed_proc<-as(imdata_ed,"MSProcessedImagingExperiment")
          gc()
          

          
          #peaklist<-summarizeFeatures(imdata_ed,"sum", as="DataFrame")
          #peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
          #peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
          
          #peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=instrument_ppm,threshold=0)
          #imdata_ed<-imdata_ed %>% peakBin(peaklist_deco$mz, resolution=instrument_ppm, units="ppm") %>% process()
          
        } 
        else if(ppm>=25){
          #smoothSignal(method="gaussian") %>%
          #reduceBaseline(method="locmin") %>%
          if (!is.null(preprocess$normalize)){
            if (preprocess$normalize$method=="Disable") {
            } else if (preprocess$normalize$method %in% c("rms","tic")){
              imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
            } else if ('&'(preprocess$normalize$method == "reference", !is.null(preprocess$normalize$mz))){
              norm_feature<-which(dplyr::between(imdata_ed@featureData@mz,
                                                 preprocess$normalize$mz*(1-ppm/1000000),
                                                 preprocess$normalize$mz*(1+ppm/1000000)))
              if (length(norm_feature)>=1){
                imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method, feature = norm_feature) %>% process(BPPARAM=SerialParam())
              }
            } else {
              imdata_ed<- imdata_ed %>% normalize(method="rms") %>% process()
            }
          }
          
          if (preprocess$smoothSignal$method=="Disable") {
          }else if (!is.null(preprocess$smoothSignal$method)){
            imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
          }else{
            imdata_ed<- imdata_ed %>% smoothSignal(method="gaussian")
          }
          
          if (preprocess$reduceBaseline$method=="Disable") {
          }else if (!is.null(preprocess$reduceBaseline$method)){
            imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
          }else{
            imdata_ed<- imdata_ed %>% reduceBaseline(method="locmin")
          }
          
          if (preprocess$peakPick$method=="Disable") {
          }else if (preprocess$peakPick$method %in% c("adaptive","mad","simple")){
            imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method, window=4) %>% process()
          }else if (preprocess$peakPick$method == "Default|default"){
            #add an peak picking function other than the Cardinal options.
            imdata_ed<-imdata_ed %>% peakBin(peaklist_deco$mz, tolerance=ppm, units="ppm") %>% process()
          }
          
          if(is.null(preprocess$peakAlign$level)) preprocess$peakAlign$level<-"local"
          if (preprocess$peakAlign$level=="global"){
            if (preprocess$peakAlign$tolerance==0 ) {
              message("preprocess$peakAlign$tolerance set as zero, step bypassed")
            }else if ('&'(!is.null(preprocess$peakAlign$tolerance),!is.null(preprocess$peakAlign$units))){
              message("preprocess$peakAlign$tolerance set as ", preprocess$peakAlign$tolerance)
              imdata_ed<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units)
            }else {
              message("preprocess$peakAlign$tolerance missing, use default tolerance in ppm ", ppm/2)
              imdata_ed<- imdata_ed %>% peakAlign(tolerance=ppm/2, units="ppm")
            }
          }
          
          imdata_ed<- imdata_ed %>% process()
          
          
          
        }
        
        saveRDS(imdata_ed,paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
      }else{
        if (!is.null(preprocess$mz_bin_list)){
          imdata_ed<-imdata
        }
        imdata_ed<-imdata
        saveRDS(imdata_ed,paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
      }
    }else if (file.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))){
      imdata_ed<-readRDS(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
    }
    
    if ('&'(preprocess$use_preprocessRDS,file.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS")))){
      message("Using image data: ",paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
      
      imdata_ed<-readRDS(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
      #imdata_ed<-imdata
      if(mzrange[1]=="auto-detect"){
        imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam())
      }else {
        imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam(),mass.range =mzrange)
      }
    }else{
      message("Using image data: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
      if(mzrange[1]=="auto-detect"){
        imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam())
      }else {
        imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam(),mass.range =mzrange)
      }
      if(!is.na(rotate[datafile_imzML[z]])){
        imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
      }else if(!is.na(rotate[datafile[z]])){
        imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
      }
      if (!is.null(preprocess)){
        #imdata_ed<-imdata
        if ('|'(imdata@metadata[["ibd binary type"]]!="processed",preprocess$force_preprocess)){
          imdata_ed<-readRDS(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
        }else{
          imdata_ed<-imdata
        }
        
      }else{
        imdata_ed<-imdata
      }
    }
    
    
    
    #imdata_org<-imdata
    imdata<-imdata_ed
    gc()
    coordata=as.data.frame(imdata@elementMetadata@coord)
    setCardinalBPPARAM(BPPARAM)
    setwd(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
    
    if (Bypass_Segmentation!=T){
      message("Segmentation in progress...")
      #cl=autoStopCluster(makeCluster(6))
      
      #Prepare imdata for segmentation
      if (Segmentation[1] %in% c("PCA","spatialKMeans","spatialShrunkenCentroids")){
        if (exists("preprocess$peakAlign$tolerance")){
          if (preprocess$peakAlign$tolerance==0 ) {
            message("preprocess$peakAlign$tolerance set as zero, step bypassed")
          }else if ('&'(!is.null(preprocess$peakAlign$tolerance),!is.null(preprocess$peakAlign$units))){
            message("preprocess$peakAlign$tolerance set as ", preprocess$peakAlign$tolerance)
            imdata_stats<- imdata %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units)
          }
        else{
          message("preprocess$peakAlign$tolerance missing, use default tolerance in ppm ", ppm/2)
          imdata_stats <- imdata %>% peakAlign(tolerance=ppm/2, units="ppm")
        }
        }
        if (is.null(preprocess$peakFilter$freq.min)){
          imdata_stats<-imdata_stats %>% peakFilter(freq.min=0.05) %>% process()
        }else{
          imdata_stats<-imdata_stats %>% peakFilter(freq.min=preprocess$peakFilter$freq.min) %>% process()
        }
      }
      
      if (Segmentation[1]=="PCA") {
        if ('&'(Segmentation_ncomp=="auto-detect",Segmentation_variance_coverage>0)){
          Segmentation_ncomp_running<-PCA_ncomp_selection(imdata_stats,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
        }else{Segmentation_ncomp_running<-Segmentation_ncomp}
      }
      else if (Segmentation[1]=="spatialKMeans" && segmentation_num!=1) {
        if ('&'(Segmentation_ncomp=="auto-detect",Segmentation_variance_coverage>0)){
          Segmentation_ncomp_running<-PCA_ncomp_selection(imdata_stats,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
        }else{Segmentation_ncomp_running<-Segmentation_ncomp}
        set.seed(1)
        
        skm <-  suppressMessages(suppressWarnings(spatialKMeans(imdata_stats, r=Smooth_range, k=segmentation_num, method="adaptive",ncomp=Segmentation_ncomp_running,BPPARAM =BPPARAM )))
        message(paste0(Segmentation[1], " finished: ",name))
        png(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=F)
        imagefile<-Cardinal::image(skm, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], key=T, ann=FALSE,axes=FALSE)
        print(imagefile)
        dev.off()
        suppressMessages(suppressWarnings(require(magick)))
        skmimg<-image_read(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""))
        
        
        png(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""),width = 1024,height = 480*((segmentation_num)))
        centers=skm@resultData[[1]][["centers"]]
        
        centers_mz<-skm@featureData@mz
        suppressMessages(suppressWarnings(require(ggplot2)))
        sp<-NULL
        sp_plot<-NULL
        for (region in colnames(centers)){
          sp_plot[[region]]<-data.frame(mz=centers_mz,intensity=centers[,region])
          sp[[region]]<-ggplot2::ggplot(data=sp_plot[[region]], size=1 ,aes(x=mz, y=intensity,xend=mz,yend=rep(0,length( sp_plot[[region]]$intensity)),colour =brewer.pal_n(segmentation_num,colorstyle)[as.numeric(region)])) +
            geom_segment(show.legend=F,colour =brewer.pal_n(segmentation_num,colorstyle)[as.numeric(region)]) +
            theme_classic() +
            ggtitle(paste("Mean spectrum"," Segmentation:",region),)
        }
        suppressMessages(suppressWarnings(require(gridExtra)))
        grid.arrange( grobs = sp,ncol=1, nrow = ceiling(segmentation_num) )
        dev.off()
        
        skmimg_spec<-image_read(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""))
        skmimg<-image_append(c(skmimg,skmimg_spec),stack = T)
        image_write(skmimg,paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs_append.png",sep=""))
        correlation=as.data.frame(skm@resultData@listData[[1]][["correlation"]])
        correlation[,"mz"]<-as.numeric(gsub("m/z = ","",rownames(correlation)))
        centers=as.data.frame(skm@resultData[[1]][["centers"]])
        centers[,"mz"]<-as.numeric(gsub("m/z = ","",rownames(centers)))
        cluster=as.data.frame(skm@resultData[[1]][["cluster"]])
        cluster$Coor=rownames(cluster)
        cluster_df<-data.frame(Coor=rownames(cluster),class=skm@resultData[[1]][["cluster"]])
        write.csv(cluster_df,paste(Segmentation[1],"_RESULT","cluster",segmentation_num,"segs.csv"),row.names = F)
        y=skm@resultData[[1]][["cluster"]]
        y=as.data.frame(y)
        rownames(y)=1:nrow(y)
        y$pixel=1:nrow(y)
        regions=unique(y[,1])
        x=NULL
        for (i in 1:length(regions)){
          listname=as.character(regions[i])
          x[[listname]]<-y[y[, 1] == regions[i], "pixel"]
        }
      }
      else if (Segmentation[1]=="spatialShrunkenCentroids" && segmentation_num!=1) {
        set.seed(1)
        if ('&'(Segmentation_ncomp=="auto-detect",Segmentation_variance_coverage>0)){
          Segmentation_ncomp_running<-PCA_ncomp_selection(imdata_stats,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
        }else{Segmentation_ncomp_running<-Segmentation_ncomp}
        skm <-  suppressMessages(suppressWarnings(spatialShrunkenCentroids(imdata_stats, r=Smooth_range, k=segmentation_num, method="adaptive",s=3,BPPARAM =BPPARAM)))
        message(paste0(Segmentation[1], " finished: ",name))
        png(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=F)
        imagefile<-Cardinal::image(skm, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], key=T, ann=FALSE,axes=FALSE, model=list(s=3), values="class")
        print(imagefile)
        dev.off()
        suppressMessages(suppressWarnings(require(magick)))
        skmimg<-image_read(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""))
        
        
        png(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""),width = 1024,height = 480*((segmentation_num)))
        centers=skm@resultData[[1]][["centers"]]
        
        centers_mz<-skm@featureData@mz
        suppressMessages(suppressWarnings(require(ggplot2)))
        sp<-NULL
        sp_plot<-NULL
        for (region in colnames(centers)){
          sp_plot[[region]]<-data.frame(mz=centers_mz,intensity=centers[,region])
          sp[[region]]<-ggplot2::ggplot(data=sp_plot[[region]], size=1 ,aes(x=mz, y=intensity,xend=mz,yend=rep(0,length( sp_plot[[region]]$intensity)),colour =brewer.pal_n(segmentation_num,colorstyle)[as.numeric(region)])) +
            geom_segment(show.legend=F,colour =brewer.pal_n(segmentation_num,colorstyle)[as.numeric(region)]) +
            theme_classic() +
            ggtitle(paste("Mean spectrum"," Segmentation:",region),)
        }
        suppressMessages(suppressWarnings(require(gridExtra)))
        grid.arrange( grobs = sp,ncol=1, nrow = ceiling(segmentation_num) )
        dev.off()
        
        skmimg_spec<-image_read(paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""))
        skmimg<-image_append(c(skmimg,skmimg_spec),stack = T)
        image_write(skmimg,paste(getwd(),"/",Segmentation[1],"_image_plot_",segmentation_num,"_segs_append.png",sep=""))
        correlation=as.data.frame(skm@resultData@listData[[1]][["correlation"]])
        correlation[,"mz"]<-as.numeric(gsub("m/z = ","",rownames(correlation)))
        centers=as.data.frame(skm@resultData[[1]][["centers"]])
        centers[,"mz"]<-as.numeric(gsub("m/z = ","",rownames(centers)))
        cluster=as.data.frame(skm@resultData[[1]][["class"]])
        cluster$Coor=rownames(cluster)
        cluster_df<-data.frame(Coor=rownames(cluster),class=skm@resultData[[1]][["class"]])
        #write.csv(correlation,paste(Segmentation[1],"_RESULT","correlation",segmentation_num,"segs.csv"),row.names = F)
        write.csv(cluster_df,paste(Segmentation[1],"_RESULT","centers",segmentation_num,"segs.csv"),row.names = F)
        #write.csv(cluster,paste(Segmentation[1],"_RESULT","cluster",segmentation_num,"segs.csv"),row.names = F)
        
        
        
        y=skm@resultData@listData[[1]][["class"]]
        y=as.data.frame(y)
        rownames(y)=1:nrow(y)
        y$pixel=1:nrow(y)
        regions=unique(y[,1])
        x=NULL
        for (i in 1:length(regions)){
          listname=as.character(regions[i])
          x[[listname]]<-y[y[, 1] == regions[i], "pixel"]
        }
        
        
        
        
        
        
      }
      else if (Segmentation[1]=="Virtual_segmentation"){
        radius_rank=read.csv(file = paste0(workdir[1],"/",Virtual_segmentation_rankfile))
        radius_rank=radius_rank[order(radius_rank$Rank),]
        if (is.null(radius_rank$Core)) radius_rank$Core="central"
        
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
        coordistmatrix$sum=base::unlist(bplapply(1:nrow(coordata),function(j,coordistmatrix,coordata){
          coordistmatrix$sum[j]=sum(coordistmatrix[j,1:nrow(coordata)])
        },coordistmatrix,coordata,BPPARAM = BPPARAM))
        coorrange=max(coordistmatrix$sum)-min(coordistmatrix$sum)
        
        
        
        findedge<-function(coordata,center_type=c("central","southeast","northeast","northwest","southwest")){
          if (is.null(center_type)) center_type="central"
          
          if (center_type=="central"){
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
            
          }else {
            if(center_type=="southeast"){
              uniquex=unique(coordata$x)
              uniquey=unique(coordata$y)
              coordata$edge=FALSE
              for (x in uniquex){
                
                min=min(coordata[coordata$x==x,"y"])
                max=max(coordata[coordata$x==x,"y"])
                #coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
                coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
              }
              
              for (y in uniquey){
                min=min(coordata[coordata$y==y,"x"])
                max=max(coordata[coordata$y==y,"x"])
                #coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
                coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
              }
            } else if(center_type=="northeast"){
              uniquex=unique(coordata$x)
              uniquey=unique(coordata$y)
              coordata$edge=FALSE
              for (x in uniquex){
                
                min=min(coordata[coordata$x==x,"y"])
                max=max(coordata[coordata$x==x,"y"])
                coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
                #coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
              }
              
              for (y in uniquey){
                min=min(coordata[coordata$y==y,"x"])
                max=max(coordata[coordata$y==y,"x"])
                #coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
                coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
              }
            }else if(center_type=="northwest"){
              uniquex=unique(coordata$x)
              uniquey=unique(coordata$y)
              coordata$edge=FALSE
              for (x in uniquex){
                
                min=min(coordata[coordata$x==x,"y"])
                max=max(coordata[coordata$x==x,"y"])
                coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
                #coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
              }
              
              for (y in uniquey){
                min=min(coordata[coordata$y==y,"x"])
                max=max(coordata[coordata$y==y,"x"])
                coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
                #coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
              }
            }else if(center_type=="southwest"){
              uniquex=unique(coordata$x)
              uniquey=unique(coordata$y)
              coordata$edge=FALSE
              for (x in uniquex){
                
                min=min(coordata[coordata$x==x,"y"])
                max=max(coordata[coordata$x==x,"y"])
                #coordata['&'(coordata$y==max,coordata$x==x),"edge"]=TRUE
                coordata['&'(coordata$y==min,coordata$x==x),"edge"]=TRUE
              }
              
              for (y in uniquey){
                min=min(coordata[coordata$y==y,"x"])
                max=max(coordata[coordata$y==y,"x"])
                coordata['&'(coordata$x==max,coordata$y==y),"edge"]=TRUE
                #coordata['&'(coordata$x==min,coordata$y==y),"edge"]=TRUE
              }
            }}
          
          return(coordata)
          
        }
        
        coordata=findedge(coordata,center_type=unique(radius_rank$Core))
        
        rank_pixel<-function(coordata,coordistmatrix,radius_rank){
          if (unique(radius_rank$Core)=="central"){
            shape_center=coordata[coordistmatrix$sum==min(coordistmatrix$sum),]
            center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),1:nrow(coordata)])
            library(useful)
            From <- shape_center[rep(seq_len(nrow(shape_center)), each=nrow(coordata)),1:2]
            To <- coordata[,1:2]
            df=To-From
            center_edge_angle<-cbind(coordata[,1:2],cart2pol(df$x, df$y, degrees = F),edge=coordata[,"edge"])
            center_edge_angle_sdge=center_edge_angle[center_edge_angle$edge==TRUE,]
            coordata$rank=0
            coordata$pattern=""
            radius_rank$center<-(radius_rank$Radius_L+radius_rank$Radius_U)/2
            radius_rank<-radius_rank[order(radius_rank$center),]
            rank_interval<-c(unlist(radius_rank$Radius_L),radius_rank$Radius_U[nrow(radius_rank)])
            
            for (i in 1: (nrow(coordata))){
              
              if (coordata$edge[i]!=TRUE){
                df=coordata[i,1:2]-shape_center[,1:2]
                point_center_angle<-cbind(coordata[i,1:2],cart2pol(df$x, df$y, degrees = F))
                pointedge=center_edge_angle_sdge[which(abs(center_edge_angle_sdge$theta-point_center_angle$theta)==min(abs(center_edge_angle_sdge$theta-point_center_angle$theta))),]
                pointedge=pointedge[which.min(pointedge$r),]
                to_edge=coordistmatrix[[i]]['&'(coordata$x==pointedge$x,coordata$y==pointedge$y)]
              }else{to_edge=0}
              
              
              to_center=center_dist[i]
              total=to_edge+to_center
              
              norm_center_dist=to_center/total*max(radius_rank$Radius_U)
              coordata$rank[i]<-radius_rank$Rank[findInterval(norm_center_dist,rank_interval,all.inside=T)]
              coordata$pattern[i]<-radius_rank$Rank[findInterval(norm_center_dist,rank_interval,all.inside=T)]  
              
            }
            coordata$rank<-factor(coordata$rank)
            coordata
          }else {
            if(unique(radius_rank$Core)=="southeast"){
              shape_center=coordata[1,]
              shape_center$x=max(coordata$x)
              shape_center$y=max(coordata$y)
              coordata_new<-rbind(coordata,shape_center)
              center_dist=bplapply(nrow(coordata_new),coordist_para,coordata_new,BPPARAM = BPPARAM)[[1]]
            } else if(unique(radius_rank$Core)=="northeast"){
              shape_center=coordata[1,]
              shape_center$x=max(coordata$x)
              shape_center$y=min(coordata$y)
              
            }else if(unique(radius_rank$Core)=="northwest"){
              shape_center=coordata[1,]
              shape_center$x=min(coordata$x)
              shape_center$y=min(coordata$y)
              
            }else if(unique(radius_rank$Core)=="southwest"){
              shape_center=coordata[1,]
              shape_center$x=min(coordata$x)
              shape_center$y=max(coordata$y)
              
            }
            
            library(useful)
            coordata$rank=0
            coordata$pattern=""
            From <- shape_center[rep(seq_len(nrow(shape_center)), each=nrow(coordata)),1:2]
            To <- coordata[,1:2]
            df=To-From
            center_edge_angle<-cbind(coordata[,1:2],cart2pol(df$x, df$y, degrees = F),edge=coordata[,"edge"])
            center_edge_angle_sdge=center_edge_angle[center_edge_angle$edge==TRUE,]
            coordata$rank=0
            coordata$pattern=""
            radius_rank$center<-(radius_rank$Radius_L+radius_rank$Radius_U)/2
            radius_rank<-radius_rank[order(radius_rank$center),]
            rank_interval<-c(unlist(radius_rank$Radius_L),radius_rank$Radius_U[nrow(radius_rank)])
            for (i in 1: (nrow(coordata))){
              
              df=coordata[i,1:2]-shape_center[,1:2]
              point_center_angle<-cbind(coordata[i,1:2],cart2pol(df$x, df$y, degrees = F))
              pointedge=center_edge_angle_sdge[which(abs(center_edge_angle_sdge$theta-point_center_angle$theta)==min(abs(center_edge_angle_sdge$theta-point_center_angle$theta))),]
              pointedge=pointedge[which.min(pointedge$r),]
              to_edge=coordistmatrix[[i]]['&'(coordata$x==pointedge$x,coordata$y==pointedge$y)]
              to_center=center_dist[i]
              total=to_edge+to_center
              
              norm_center_dist=to_center/total*max(radius_rank$Radius_U)
              
              coordata$rank[i]<-radius_rank$Rank[findInterval(norm_center_dist,rank_interval,all.inside=T)]
              coordata$pattern[i]<-radius_rank$Rank[findInterval(norm_center_dist,rank_interval,all.inside=T)]  
            }
            coordata$rank<-factor(coordata$rank)
            coordata
          }
        }
        coordata=rank_pixel(coordata,coordistmatrix,radius_rank)
        x=NULL
        
        for (rank in unique(coordata$rank)){
          
          x[[rank]]=which(coordata$rank==rank)
          
        }
        write.csv(coordata,"coordata.csv",row.names = F)
        
        region_pattern <- factor(coordata$pattern,levels=unique(coordata$pattern), labels=unique(coordata$pattern))
        
        png(paste(getwd(),"/","Virtual_segmentation",gsub("/"," ",name),'.png',sep=""),width = 1024,height = 1024)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 1),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)
        print(image(imdata, region_pattern ~ x * y,col=brewer.pal_n(length(unique(coordata$rank)),colorstyle), key=TRUE))
        
        dev.off()
        
        
        
      }
      
      else if (Segmentation[1]=="def_file"){
        
        Segmentation_def_tbl<-read.csv(paste0(workdir[z],"/",paste0(gsub(".imzML$","",datafile[z])  ," ID/"), Segmentation_def))
        
        if(c("datafile") %in% colnames(Segmentation_def_tbl) ){
          Segmentation_def_tbl<-Segmentation_def_tbl[Segmentation_def_tbl$datafile==datafile[z],]
        }
        
        if(c("pixel") %in% colnames(Segmentation_def_tbl) ){
          Segmentation_def_tbl<-Segmentation_def_tbl[order(Segmentation_def_tbl$pixel),]
        }else if (c("Coor") %in% colnames(Segmentation_def_tbl)){
          Segmentation_def_tbl<-Segmentation_def_tbl[order(Segmentation_def_tbl$Coor),]
        }
        if(c("label") %in% colnames(Segmentation_def_tbl) ){
          Segmentation_def_tbl$label<-as.factor(Segmentation_def_tbl$label)
        }else if (c("class") %in% colnames(Segmentation_def_tbl)){
          Segmentation_def_tbl$label<-as.factor(Segmentation_def_tbl$class)
        }
        x=NULL
        
        for (label in levels(Segmentation_def_tbl$label)){
          
          x[[label]]<-Segmentation_def_tbl$pixel[Segmentation_def_tbl$label==label]
          
        }
        
        message("Segmentation_def_tbl loaded, labelled regions: ", paste(names(x),collapse = ", "))
        png(paste(getwd(),"/","Segmentation_def_file","_image_plot_",length(levels(Segmentation_def_tbl$label)),"_segs.png",sep=""),width = 1024,height = 720)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=F)
        imagefile<-Cardinal::image(imdata, factor(Segmentation_def_tbl$label) ~ x * y, key=T, ann=FALSE,axes=FALSE)
        print(imagefile)
        dev.off()
        
      }
      else {
        
        x=1:length(pixels(imdata))
        x=split(x, sort(x%%1))
        names(x)="1"
        png(paste(getwd(),"/","Segmentation_none","_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=F)
        imagefile<-Cardinal::image(imdata, factor(rep(1,length(pixels(imdata)))) ~ x * y, key=T, ann=FALSE,axes=FALSE)
        print(imagefile)
        dev.off()
      }
    } else if (Segmentation[1]=="none"){
      
      x=1:length(pixels(imdata))
      x=split(x, sort(x%%1))
      names(x)="1"
      png(paste(getwd(),"/","Segmentation_none","_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
      
      par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
          bty="n",pty="s",xaxt="n",
          yaxt="n",
          no.readonly = TRUE,ann=F)
      imagefile<-Cardinal::image(imdata, factor(rep(1,length(pixels(imdata)))) ~ x * y, key=T, ann=FALSE,axes=FALSE)
      print(imagefile)
      dev.off()
    }else {
      
      x=1:length(pixels(imdata))
      x=split(x, sort(x%%1))
      names(x)="1"
      png(paste(getwd(),"/","Segmentation_none","_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
      
      par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
          bty="n",pty="s",xaxt="n",
          yaxt="n",
          no.readonly = TRUE,ann=F)
      imagefile<-Cardinal::image(imdata, factor(rep(1,length(pixels(imdata)))) ~ x * y, key=T, ann=FALSE,axes=FALSE)
      print(imagefile)
      dev.off()
    }
  }
  message("workflow successfully completed")
  return(list(segmentation_label=x,imdata_org=NULL,imdata=imdata,imdata_ed=imdata_ed))
}



Parse_rotation<-function(datafile,rotate,wd=getwd()){

  if (!is.null(rotate)){

    message("Found rotation info")

    if (typeof(rotate)=="character"){
      if (file.exists(rotate)){
        rotate=read.csv(paste0(rotate),stringsAsFactors = F)
      }else{
        rotate=read.csv(paste0(wd,"/",rotate),stringsAsFactors = F)
      }

    }else if ('|'(typeof(rotate)=="double",typeof(rotate)=="integer")){
    return(rotate)
    }
    if (typeof(rotate)=="list"){
    rotatedegrees=sapply(datafile,function(x,df){
      library(stringr)
      df$filenames<-str_remove(df$filenames,"\\.imzML$")

      degree=df$rotation[df$filenames==(x)]
      if (length(degree)==0) {
        message("Missing rotation data please check the rotation configuration file: ",x)
        degree=0
      }
      
      unname(degree)
    },rotate)
    rotate=unlist(rotatedegrees)
    }

  }else {rotate=rep(0,length(datafile))}
  return(rotate)
}

Load_IMS_combine<-function(datafile,rotate=NULL,ppm=5,...){
  if ('&'(!is.null(Rotate_IMG),typeof(Rotate_IMG)=="character")){
    Rotate_IMG=read.csv(Rotate_IMG,stringsAsFactors = F)
  }
  imdata=list()
  combinedimdata=NULL
  #register(SerialParam())
  Rotate_IMG$filenames<-gsub(".imzML$","",Rotate_IMG$filenames)
  if (!exists("mzrange")){
    mzrange=NULL
    testrange=c(0,0)
    for (i in 1:length(datafile)){
      rotate=Rotate_IMG[Rotate_IMG$filenames==basename(datafile[i]),"rotation"]
      rotate=as.numeric(rotate)
      if (length(rotate)==0){rotate=0}
      imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm,rotate = rotate,attach.only=F,as=as,mzrange=mzrange)
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
    rotate=Rotate_IMG[Rotate_IMG$filenames==basename(datafile[i]),"rotation"]
    rotate=as.numeric(rotate)
    if (length(rotate)==0){rotate=0}
    #if (length(datafile)==1){
    #imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm*2,rotate = rotate,attach.only=F,as="MSImageSet",mzrange=mzrange)
    #}else{
    imdata[[i]]=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm*2,rotate = rotate,attach.only=F,as=as,mzrange=mzrange)
    #}
    #imdata[[i]]@elementMetadata@coord=imdata[[i]]@elementMetadata@coord[,c("x","y")]
    #max(imdata[[i]]@featureData@mz)
    #min(imdata[[i]]@featureData@mz)
    if (i==1) {
      combinedimdata=imdata[[i]]
    }else{
      combinedimdata<-cbind(combinedimdata,imdata[[i]])
    }
    imdata[[i]]=NULL
  }

  combinedimdata@elementMetadata@coord@listData[["z"]]=NULL

  imdata=combinedimdata

  return(imdata)
}

Load_IMS_decov_combine<-function(datafile,workdir,ppm=5,import_ppm=ppm/2,SPECTRUM_batch="overall",mass_correction_tol_ppm=12,mzAlign_runs=c("TopNfeature_mean","Combined_features"),
                                 threshold=0.0000001,rotate=NULL,mzrange="auto-detect", ppm_aligment=ppm, 
                                 deconv_peaklist=c("Load_exist","New","Target"),preprocessRDS_rotated=T,use_rawdata=F,target_mzlist=NULL){
  suppressMessages(suppressWarnings(library(matter)))
  suppressMessages(suppressWarnings(library(stringr)))
  suppressMessages(suppressWarnings(library(HiTMaP)))
  suppressMessages(suppressWarnings(library(Cardinal)))
  
  datafile_base<-basename(datafile)
  datafile <- str_remove(datafile_base,"\\.imzML$")
  if(length(workdir)!=length(datafile)){
    workdir=rep(workdir[1],length(datafile))
  }
  
  datafile_imzML=paste0(datafile,".imzML")
  rotate=HiTMaP:::Parse_rotation(datafile,rotate)
  if (sum(deconv_peaklist[1]=="New",!file.exists(paste0(workdir[1],"/","ClusterIMS_deconv_Spectrum.csv")))>=1){
    for (z in 1:length(datafile)){
      name <-basename(datafile[z])
      name <-gsub(".imzML$","",name)
      name <-gsub("/$","",name)
      setwd(workdir[z])
      folder<-base::dirname(datafile[z])
      #imdata <- Cardinal::readImzML(datafile[z],preprocessing = F,attach.only = T,resolution = 200,rotate = rotate[z],as="MSImageSet",BPPARAM = BPPARAM)
      if (!str_detect(datafile[z],".imzML$")){
        datafile_imzML[z]<-paste0(datafile[z],".imzML")
      }
      
      if (dir.exists(paste0(gsub(".imzML$","",datafile[z]) ," ID"))==FALSE){
        dir.create(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
      }
      
      
      if ('&'(file.exists(paste0(workdir[z],"/",gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS")),!use_rawdata)){
        imdata<-readRDS(paste0(workdir[z],"/",datafile[z]," ID/","preprocessed_imdata.RDS"))
        if (!preprocessRDS_rotated){
          if(!is.na(rotate[datafile_imzML[z]])){
            imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
          }else if(!is.na(rotate[datafile[z]])){
            imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
          }
        }
      }else if (use_rawdata){
        message("Loading raw image data for statistical analysis: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
        if(mzrange[1]=="auto-detect"){
          imdata <- Cardinal::readMSIData(paste0(workdir[z],"/",datafile_imzML[z]),  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam())
        }else {
          imdata <- Cardinal::readMSIData(paste0(workdir[z],"/",datafile_imzML[z]),  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam(),mass.range =mzrange)
        }
        if(!is.na(rotate[datafile_imzML[z]])){
          imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
        }else if(!is.na(rotate[datafile[z]])){
          imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
        }
      }
      
      file_deconv_Spectrum<-paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"_deconv_Spectrum.csv")
      message(file_deconv_Spectrum)
      
      if (`|`(deconv_peaklist[1]=="New",!file.exists(file_deconv_Spectrum[1]))){
      spectrum_file_table<- summarizeFeatures(imdata, FUN = "mean")
      spectrum_file_table<-data.frame(mz=spectrum_file_table@featureData@mz,mean=spectrum_file_table@featureData@listData[["mean"]])
      peaklist<-spectrum_file_table
      colnames(peaklist)<-c("m.z","intensities")
      savename=paste(name,SPECTRUM_batch)
      message(paste("Mean spectrum generated",name,"region",SPECTRUM_batch))
      write.csv(peaklist,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"_Spectrum.csv"),row.names = F)
      peaklist<-peaklist[peaklist$intensities>0,]
      deconv_peaklist_df<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist,ppm=ppm,threshold=threshold)
      write.csv(deconv_peaklist_df,paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"_deconv_Spectrum.csv"),row.names = F)
    
      }
    }
    deconv_peaklist_list<-NULL
    
    for (z in 1:length(datafile)){
      
      deconv_peaklist_list[[datafile[z]]] <- read.csv(paste0(workdir[z],"/",datafile[z] ," ID/",SPECTRUM_batch,"_deconv_Spectrum.csv"))
      
    }
    
    if (mzAlign_runs=="TopNfeature_mean"){
      message("Performing mz alignment across all runs..." )
      
      mz.ref.list.top.quantile<-lapply(deconv_peaklist_list,function(x,quantile=0.15){
        maxs<-matter::locmax(x$intensities)
        cutoff <- quantile(x$intensities[maxs], 1 - quantile)
        maxs <- maxs[x$intensities[maxs] >= cutoff]
        mz.ref <- x$m.z[maxs]
        return(mz.ref)
      })
      
      
      
      mz.ref.list.top.quantile.spec<-lapply(names(mz.ref.list.top.quantile),function(x,mz.ref.list.top.quantile,deconv_peaklist_list){
        
        deconv_peaklist_list[[x]][deconv_peaklist_list[[x]]$m.z %in% mz.ref.list.top.quantile[[x]],] 
        
      },mz.ref.list.top.quantile,deconv_peaklist_list)
      names(mz.ref.list.top.quantile.spec)<-names(mz.ref.list.top.quantile)
      mz.ref.list.top.quantile.spec.bind<- do.call(rbind,mz.ref.list.top.quantile.spec)
      mz.ref.list.top.quantile.spec.bind$intensities<-1
      mz.ref.list.top.quantile.spec.bind<-mz.ref.list.top.quantile.spec.bind[order(mz.ref.list.top.quantile.spec.bind$m.z),]
      rownames(mz.ref.list.top.quantile.spec.bind)<-1:nrow(mz.ref.list.top.quantile.spec.bind)
      mz.ref.list.top.quantile.bin<-HiTMaP:::isopattern_ppm_filter_peaklist(mz.ref.list.top.quantile.spec.bind,ppm=ppm_aligment,threshold=0.00)
      mz.ref.list.top.quantile.final<-mz.ref.list.top.quantile.bin$m.z[mz.ref.list.top.quantile.bin$intensities>(0.55*max(mz.ref.list.top.quantile.bin$intensities))]
      #match features
      
      deconv_peaklist_ref_match <- NULL
      deconv_peaklist_ref_match_locmax <- NULL
      for (z in 1:length(datafile)){
        IMS_datafile_aligment<-function(mz,mz.ref,mz.test,control = loess.control(),span = 0.75,tol_ppm=5,tol_ppm_ext=tol_ppm * 1.8,tol.ref="key"){
          suppressMessages(suppressWarnings(library(matter)))
            tol = tol_ppm * 1e-6
          i=bsearch(mz.ref,mz.test, tol=tol, tol.ref="key")
        found <- !is.na(i)
        if ( sum(found) < 1 ) {
          warning("no matching peaks found; try a larger tolerance")
          return(x)
        }
        if ( sum(found) < 0.5*length(mz.ref) ) {
          message("Less than 50% ref peaks found in: ",z," ",datafile[z],"; retrying with a larger tolerance")
          tol = tol_ppm_ext * 1e-6
          i=bsearch(mz.ref,mz.test, tol=tol, tol.ref="key")
          found <- !is.na(i)
          if ( sum(found) < 0.5*length(mz.ref) ) {
            message("Less than 50% ref peaks found in: ",z," ",datafile[z],"; consider using a larger tolerance")
            
          }
        }
        mz.ref <- mz.ref[found]
        i <- i[found]
        mz.test <- mz.test[i]
        diff <- mz.ref - mz.test
        diff_org <- diff
        mz.test <- c(mz[1], mz.test, mz[length(mz)])
        diff <- c(diff[1], diff, diff[length(diff)])
        shift <- suppressWarnings(loess(diff ~ mz.test, span=span, control=control))
        dmz <- predict(shift, mz)
        #warp <- splinefun(mz + dmz, x)
        return(list(final.mz = mz + dmz, dmz=dmz, mz.ref=mz.ref, diff = diff_org, shift = shift))
        #return(list(final.mz = mz + dmz, dmz=dmz))
        }
        # deconv_peaklist_ref_match[[datafile[z]]] <- IMS_datafile_aligment(mz = deconv_peaklist_list[[datafile[z]]]$m.z,
        #                                                                   mz.ref = mz.ref.list.top.quantile.final,
        #                                                                   mz.test = deconv_peaklist_list[[datafile[z]]]$m.z,
        #                                                                   control = loess.control(),
        #                                                                   tol_ppm=5,
        #                                                                   span = 0.75)
         
        deconv_peaklist_ref_match_locmax[[datafile[z]]] <- try(IMS_datafile_aligment(mz = mz.ref.list.top.quantile.spec[[datafile[z]]]$m.z,
                                                                          mz.ref = mz.ref.list.top.quantile.final,
                                                                          mz.test = mz.ref.list.top.quantile.spec[[datafile[z]]]$m.z,
                                                                          control = loess.control(surface = "direct"),
                                                                          tol_ppm=ppm_aligment,
                                                                          span = 0.75))
      }
      
     
      suppressMessages(suppressWarnings(library(ggplot2)))
      #deconv_peaklist_ref_match_df<-dplyr::bind_rows(lapply(deconv_peaklist_ref_match, function(x) list(mz=x$final.mz, dmz=x$dmz)), .id = 'file')
      #deconv_peaklist_ref_match_df_ref<-dplyr::bind_rows(lapply(deconv_peaklist_ref_match, function(x) list(mz=x$mz.ref, dmz=x$diff)), .id = 'file')
      
      
      
      deconv_peaklist_ref_match_df<-dplyr::bind_rows(lapply(deconv_peaklist_ref_match_locmax, function(x) list(mz=x$final.mz, dmz=x$dmz)), .id = 'file')
      deconv_peaklist_ref_match_df_ref<-dplyr::bind_rows(lapply(deconv_peaklist_ref_match_locmax, function(x) list(mz=x$mz.ref, dmz=x$diff)), .id = 'file')
      
      
      deconv_peaklist_ref_match_df$type<-"org"
      deconv_peaklist_ref_match_df_ref$type<-"ref"
      
      deconv_peaklist_ref_match_df$Batch<-deconv_peaklist_ref_match_df$file
      deconv_peaklist_ref_match_df_ref$Batch<-deconv_peaklist_ref_match_df_ref$file
      deconv_peaklist_ref_match_df$Batch[str_detect(deconv_peaklist_ref_match_df$Batch,"24Hr|48Hr|72Hr")]<-"Batch 2"
      deconv_peaklist_ref_match_df$Batch[str_detect(deconv_peaklist_ref_match_df$Batch,"aah_|ff_")]<-"Batch 3"
      deconv_peaklist_ref_match_df$Batch[str_detect(deconv_peaklist_ref_match_df$Batch,"SILGlucose|20201116|20201109")]<-"Batch 1"
      deconv_peaklist_ref_match_df_ref$Batch[str_detect(deconv_peaklist_ref_match_df_ref$Batch,"24Hr|48Hr|72Hr")]<-"Batch 2"
      deconv_peaklist_ref_match_df_ref$Batch[str_detect(deconv_peaklist_ref_match_df_ref$Batch,"aah_|ff_")]<-"Batch 3"
      deconv_peaklist_ref_match_df_ref$Batch[str_detect(deconv_peaklist_ref_match_df_ref$Batch,"SILGlucose|20201116|20201109")]<-"Batch 1"
      #deconv_peaklist_ref_match_df<-rbind(deconv_peaklist_ref_match_df,deconv_peaklist_ref_match_df_ref)
      
      mz.ref.list.top.quantile.spec.corrected<-mz.ref.list.top.quantile.spec
      for (z in 1:length(datafile)){
        mz.ref.list.top.quantile.spec.corrected[[datafile[z]]]$m.z <- predict(deconv_peaklist_ref_match_locmax[[datafile[z]]]$shift,mz.ref.list.top.quantile.spec.corrected[[datafile[z]]]$m.z) + mz.ref.list.top.quantile.spec.corrected[[datafile[z]]]$m.z 
      }
      mz.ref.list.top.quantile.spec.corrected_df<-dplyr::bind_rows(lapply(mz.ref.list.top.quantile.spec.corrected, function(x) list(mz=x$m.z)), .id = 'file')
      mz.ref.list.top.quantile.spec.org_df<-dplyr::bind_rows(lapply(mz.ref.list.top.quantile.spec, function(x) list(mz=x$m.z)), .id = 'file')
      mz.ref.list.top.quantile.spec.crossvalid<-mz.ref.list.top.quantile.spec.corrected_df
      mz.ref.list.top.quantile.spec.crossvalid$Corrected<-mz.ref.list.top.quantile.spec.corrected_df$mz
      mz.ref.list.top.quantile.spec.crossvalid$Original<-mz.ref.list.top.quantile.spec.org_df$mz
      
      png(paste0(workdir[1],"/","ClusterIMS_mz_correction.png"),width = 2160, height = 1080, res = 300)
      suppressMessages(suppressWarnings(library(egg)))
      getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Spectral"))
      g<-ggplot(deconv_peaklist_ref_match_df,aes(x=mz,y=dmz/mz* 1e6,group=file,colour=Batch)) + 
        geom_line() + 
        geom_point(data=deconv_peaklist_ref_match_df_ref,mapping=aes(x=mz,y=dmz/mz * 1e6,group=file,colour=Batch),size=.5) +
        scale_fill_manual(values = getPalette(length(unique(deconv_peaklist_ref_match_df$Batch))))+
        egg::theme_article()+
        labs(title ="",x = "m/z",y = "mass error (in ppm)") 
      g2<-ggplot(mz.ref.list.top.quantile.spec.crossvalid,aes(x=Original,y=Corrected,group=file,colour=file)) + 
        geom_line() + 
        scale_fill_manual(values = getPalette(length(unique(deconv_peaklist_ref_match_df$file))))+
        egg::theme_article()+
        labs(title ="",x = "m/z",y = "mass error (in ppm)") +
        guides(colour = "none") 
      print(g)
      dev.off()
      
      mz.ref.list.top.quantile.spec.corrected_df<-dplyr::bind_rows(lapply(mz.ref.list.top.quantile.spec.corrected, function(x) list(mz=x$m.z,
                                                                                                                                    intensity=x$intensities)), .id = 'file')
      mz.ref.list.top.quantile.spec.corrected_df$intensity_log<-log(mz.ref.list.top.quantile.spec.corrected_df$intensity)
      
      mz.ref.list.top.quantile.spec.corrected_df
      
      saveRDS(deconv_peaklist_ref_match_locmax,paste0(workdir[1],"/deconv_peaklist_ref_match_locmax.rds"))
      
      for (z in 1:length(datafile)){
        
        deconv_peaklist_list[[datafile[z]]]$m.z <- deconv_peaklist_list[[datafile[z]]]$m.z + 
          predict(deconv_peaklist_ref_match_locmax[[datafile[z]]]$shift, deconv_peaklist_list[[datafile[z]]]$m.z)
        
      }
    deconv_peaklist_bind<-do.call(rbind,deconv_peaklist_list)
    
    deconv_peaklist_bind<-deconv_peaklist_bind[order(deconv_peaklist_bind$m.z),]
    
    rownames(deconv_peaklist_bind)<-1:nrow(deconv_peaklist_bind)
    
    deconv_peaklist_decov<-HiTMaP:::isopattern_ppm_filter_peaklist(deconv_peaklist_bind,ppm=ppm,threshold=threshold)
    
    write.csv(deconv_peaklist_decov,paste0(workdir[z],"/","ClusterIMS_deconv_Spectrum.csv"),row.names = F)
    }else if (mzAlign_runs=="Combined_features"){
      rm(deconv_peaklist_ref_match_locmax)
      
      deconv_peaklist_bind<-do.call(rbind,deconv_peaklist_list)
      
      deconv_peaklist_bind<-deconv_peaklist_bind[order(deconv_peaklist_bind$m.z),]
      
      rownames(deconv_peaklist_bind)<-1:nrow(deconv_peaklist_bind)
      
      deconv_peaklist_decov<-HiTMaP:::isopattern_ppm_filter_peaklist(deconv_peaklist_bind,ppm=ppm,threshold=threshold)
      
      write.csv(deconv_peaklist_decov,paste0(workdir[z],"/","ClusterIMS_deconv_Spectrum.csv"),row.names = F)
      
      }
    
    
    
    
    
  }else if (sum(deconv_peaklist[1]=="Load_exist",file.exists(paste0(workdir[1],"/","ClusterIMS_deconv_Spectrum.csv")))==2){
    
    if (file.exists(paste0(workdir[1],"/deconv_peaklist_ref_match_locmax.rds"))) deconv_peaklist_ref_match_locmax<-readRDS(paste0(workdir[1],"/deconv_peaklist_ref_match_locmax.rds"))
    
    deconv_peaklist_decov<-read.csv(paste0(workdir[1],"/","ClusterIMS_deconv_Spectrum.csv"))
    
    deconv_peaklist_decov<-HiTMaP:::isopattern_ppm_filter_peaklist(deconv_peaklist_decov,ppm=ppm,threshold=threshold)
    
  }else if (sum(deconv_peaklist[1]=="Target")==1){
    if (file.exists(paste0(workdir[1],"/","ClusterIMS_deconv_Spectrum_target.csv"))){
          deconv_peaklist_decov<-read.csv(paste0(workdir[1],"/","ClusterIMS_deconv_Spectrum_target.csv"))
          deconv_peaklist_decov<-HiTMaP:::isopattern_ppm_filter_peaklist(deconv_peaklist_decov,ppm=ppm,threshold=threshold)
    }else{
     if (!is.null(target_mzlist)){
       deconv_peaklist_decov=data.frame(m.z=target_mzlist,intensities=1)
     }else{
      message("No target mz found, Please check the protein result and try again.")
    }

    
  }
  
  }
  
  
  #imdata_list<-list()
  imdata<<-NULL
  imdata<-NULL
  gc()
  combinedimdata_list<-NULL
  
  for (z in 1:length(datafile)){
    deconv_peaklist_decov_plot<-deconv_peaklist_decov
    if (!("m.z" %in% colnames(deconv_peaklist_decov_plot))) deconv_peaklist_decov_plot$m.z<-deconv_peaklist_decov_plot[,1]
    setwd(workdir[z])
    
    if ('&'(file.exists(paste0(workdir[z],"/",gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS")),!use_rawdata)){
      imdata<-readRDS(paste0(workdir[z],"/",datafile[z]," ID/","preprocessed_imdata.RDS"))
      if (!preprocessRDS_rotated){
        if(!is.na(rotate[datafile_imzML[z]])){
          imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
        }else if(!is.na(rotate[datafile[z]])){
          imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
        }
      }
          }else if (use_rawdata){
      message("Loading raw image data for statistical analysis: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
      if(mzrange[1]=="auto-detect"){
        imdata <- Cardinal::readMSIData(paste0(workdir[z],"/",datafile_imzML[z]),  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam())
      }else {
        imdata <- Cardinal::readMSIData(paste0(workdir[z],"/",datafile_imzML[z]),  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam(),mass.range =mzrange)
      }
      if(!is.na(rotate[datafile_imzML[z]])){
        imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
      }else if(!is.na(rotate[datafile[z]])){
        imdata <-HiTMaP:::rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
      }
          }
    message("mzbin for ",datafile[z])
    setCardinalBPPARAM(SerialParam())
    
    if ( exists("deconv_peaklist_ref_match_locmax")){
     message("Performing m/z correction...")
      if (class(imdata)=="MSProcessedImagingExperiment"){
        mz(imdata)<-predict(deconv_peaklist_ref_match_locmax[[datafile[z]]]$shift,mz(imdata)) + mz(imdata)
      }else{
        imdata@featureData@mz<-predict(deconv_peaklist_ref_match_locmax[[datafile[z]]]$shift,imdata@featureData@mz) + imdata@featureData@mz
      }
     
    }else(message("No batch level m/z drifting model found, m/z correction bypassed."))
    
    New_fdata_LB<-imdata[1,]
    New_fdata_LB@featureData@mz[1]<-min(deconv_peaklist_decov_plot$m.z, na.rm = T)-2*ppm/1000000*min(deconv_peaklist_decov_plot$m.z, na.rm = T)
    New_fdata_UB<-New_fdata_LB
    New_fdata_UB@featureData@mz[1]<-max(deconv_peaklist_decov_plot$m.z, na.rm = T)+2*ppm/1000000*min(deconv_peaklist_decov_plot$m.z, na.rm = T)
    
    if(min(mz(imdata))>min(New_fdata_LB@featureData@mz[1])){
      message("removing ",sum(deconv_peaklist_decov_plot$m.z<min(mz(imdata)))," m/z features from target list, out of LB of ", datafile[z])
      #imdata<-Cardinal::rbind(New_fdata_LB,imdata)
      deconv_peaklist_decov_plot<-deconv_peaklist_decov_plot[deconv_peaklist_decov_plot$m.z>=min(mz(imdata)),]
      deconv_peaklist_decov_plot->deconv_peaklist_decov
    }
    if(max(mz(imdata))<max(New_fdata_UB@featureData@mz[1])){
      message("removing ",sum(deconv_peaklist_decov_plot$m.z>max(mz(imdata)))," m/z features from target list, out of UB of ", datafile[z])
      #imdata<-Cardinal::rbind(imdata,New_fdata_UB)
      deconv_peaklist_decov_plot<-deconv_peaklist_decov_plot[deconv_peaklist_decov_plot$m.z<=max(mz(imdata)),]
      deconv_peaklist_decov_plot->deconv_peaklist_decov
    }

    
    imdata <- imdata %>%
    mzBin(deconv_peaklist_decov_plot$m.z, resolution=ppm, units="ppm")%>%
      process()
    imdata@elementMetadata@coord@listData[["z"]]<-NULL
    imdata@elementMetadata@resolution=c(x=1,y=1)

      

      
      imdata <- imdata %>%
        mzBin(deconv_peaklist_decov_plot$m.z, resolution=ppm, units="ppm")%>%
        process(BPPARAM=SnowParam(workers = 3))
      imdata@elementMetadata@coord@listData[["z"]]<-NULL
      imdata@elementMetadata@resolution=c(x=1,y=1)
    

    combinedimdata_list[[z]]<-imdata
    rm(imdata)

  }
      
  gc()
  
  for (z in 1:length(datafile)){
   
    combinedimdata_list[[z]]<-combinedimdata_list[[z]][mz(combinedimdata_list[[z]]) %in% deconv_peaklist_decov_plot$m.z,]
     
  }
  
  gc()
  
  saveRDS(combinedimdata_list,paste0(workdir[1],"/",mzAlign_runs,"_",deconv_peaklist,"combinedimdata_list.rds"),compress = T)
  
  for (z in 1:length(datafile)){
    
    combinedimdata_list[[z]]@centroided<-TRUE
    
  if (z==1) {
    
    combinedimdata<-combinedimdata_list[[z]]

  }else{

    combinedimdata<-Cardinal::cbind(combinedimdata,combinedimdata_list[[z]])

  }
    
    combinedimdata_list[[z]]<-1
    gc()
  }
  
  
  saveRDS(combinedimdata,paste0(workdir[1],"/",mzAlign_runs,"_",deconv_peaklist,"combinedimdata.rds"),compress = T)
  
  return(paste0(workdir[1],"/",mzAlign_runs,"_",deconv_peaklist,"combinedimdata.rds"))
}

sort_run_msi<-function(combinedimdata,datafiles,norm_coord=T){
  library(Cardinal)
  combinedimdata_sort<-NULL
  imdata<-NULL
  for (z in 1:length(datafile)){

    imdata <- combinedimdata[,combinedimdata@elementMetadata@run==tolower(datafile[z])]
    if(norm_coord){
      imdata <- rotateMSI(imdata,rotation_degree = 0)
    }

    if (z==1) {
      combinedimdata_sort<-imdata
    }else{
      combinedimdata_sort<-cbind(combinedimdata_sort,imdata)
    }

  }
  return(combinedimdata_sort)
}

load_pixel_label<-function(combinedimdata,datafile,workdir,coordata_file="coordata.csv",pixel_idx_col=base::row.names,label_col="pattern",...){
  library(Cardinal)
  library(stringr)
  library(HiTMaP)
  datafile_base<-basename(datafile)
  datafile <- str_remove(datafile_base,"\\.imzML$")
  if(length(workdir)!=length(datafile)){
    workdir=rep(workdir[1],length(datafile))
  }

  datafile_imzML=datafile
  coordata_file_tb<-NULL
  for (z in 1:length(datafile)){
    name <-basename(datafile[z])
    name <-gsub(".imzML$","",name)
    name <-gsub("/$","",name)
    setwd(workdir[z])
    folder<-base::dirname(datafile[z])
    #imdata <- Cardinal::readImzML(datafile[z],preprocessing = F,attach.only = T,resolution = 200,rotate = rotate[z],as="MSImageSet",BPPARAM = BPPARAM)
    if (!str_detect(datafile[z],".imzML$")){
      datafile_imzML[z]<-paste0(datafile[z],".imzML")
    }

    coordata_file_tb[[datafile[z]]]<-read.csv(paste0(workdir[z],"/",datafile[z]," ID/",coordata_file))
    coordata_file_tb[[datafile[z]]]$run<-datafile[z]
    coordata_file_tb[[datafile[z]]]$pixel_idx<-pixel_idx_col(coordata_file_tb[[datafile[z]]])
  }

  coordata_file_bind<-do.call(rbind,coordata_file_tb)
  coordata_file_bind$run<-as.factor(tolower(coordata_file_bind$run))
  coordata_file_bind$pixel_idx<-as.numeric(coordata_file_bind$pixel_idx)
  Pixel_run <-run(combinedimdata)
  Pixel_run <- data.frame(run=Pixel_run)
  Pixel_run <- Pixel_run %>% group_by(run) %>% summarise(pixel_idx=1:length(run))

  label_run <- merge(Pixel_run,coordata_file_bind,by=c("run","pixel_idx"),all.x = T,sort=F)

  Pixel_label<-label_run[[label_col]]
  return(Pixel_label)
}



intensity_sum_para<-function(mz,Peptide_plot_list){

  ifelse(length(Peptide_plot_list[Peptide_plot_list$mz==mz,"Intensity"])!=0,Peptide_plot_list[Peptide_plot_list$mz==mz,"Intensity"],0)
}

Mass_defect_plot<-function(Protein_feature_list,outputdir=NULL){

  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(ggplot2)))
  data("isotopes")
  isotopes$nominalmass=as.numeric(sapply(isotopes$isotope,function(x,pattern){
    suppressMessages(suppressWarnings(require(stringr)))
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

MassSpecWavelet_fun<-function(peaklist,SNR.Th=3){
  suppressMessages(suppressWarnings(require(MassSpecWavelet)))
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
  suppressMessages(suppressWarnings(require(stringr)))
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

# spec_peakdetect<-function(x){
#   suppressMessages(suppressWarnings(require(MALDIrppa)))
#   data(spectra) # list of MassSpectra class objects
#   data(type) # metadata
#   # Summary of spectra features (results for 20 first mass spectra)
#   summarySpectra(spectra[1:20])
#   # Some pre-processing
#   sc.results <- screenSpectra(spectra, meta = type)
#   spectra <- sc.results$fspectra
#   type <- sc.results$fmeta
#   spectra <- transformIntensity(spectra, method = "sqrt")
#   spectra <- wavSmoothing(spectra)
#   spectra <- removeBaseline(spectra)
#   names(spectra) <- type$SpectID # spectra IDs are lost with removeBaseline()
#   # Summary of spectra features (results for positions 10 to 20)
#   summarySpectra(spectra[10:20])
# }


pick.peaks <- function(peaklist, ppm) {
  span.width <- ppm * 2 + 1
  loc.max <- span.width + 1 -
    apply(embed(peaklist, span.width), 1, which.max)
  loc.max[loc.max == 1 | loc.max == span.width] <- NA

  pks <- loc.max + 0:(length(loc.max)-1)
  unique(pks[!is.na(pks)])
}

isopattern_ppm_filter_peaklist_dep<-function(pattern,ppm,threshold=0.001,verbose=T){

  org_feature=nrow(pattern)

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

  filtered_pattern<-filtered_pattern[filtered_pattern$intensities>=max(filtered_pattern$intensities)*threshold,]

  rownames(filtered_pattern)=1:nrow(filtered_pattern)
  if(verbose==T){
    message(paste("Origional Features:",org_feature,"Filtered Features:",nrow(filtered_pattern)))
  }


  return(filtered_pattern)
}

isopattern_ppm_filter_peaklist<-function(pattern,ppm,threshold=0.001,verbose=F){
  pattern<-pattern[pattern[,2]!=Inf,]
  org_feature=nrow(pattern)

  pattern<-as.data.frame(pattern)

  pattern_ppm=as.numeric(as.character(pattern[,1]))

  pattern_ppm_delta=numeric()


  filtered_pattern<-pattern[1,]

  mzrange<-range(pattern[,1])

  mzbinsecs<-ceiling((mzrange[2]-mzrange[1])/100)

  mzbins<-c()

  # mzbin<-mzrange[1]
  #
  # for (mzbinsec in mzbinsecs){
  #
  #   mzbin_u<-mzbin+100
  #
  # }
  #
  # mzbins<-c()
  #
  # mzbin<-mzrange[1]
  #
  # ppm_step<-ppm/1000000
  #
  # system.time(repeat{
  #
  #   if(mzbin>=mzrange[2]){
  #     break
  #   }
  #
  #   mzbins<-c(mzbins,mzbin)
  #
  #   mzbin=mzbin+(mzbin*ppm_step)
  #
  # })
  #filtered_pattern<-t(data.frame(filtered_pattern,stringsAsFactors = F))
  #dim(filtered_pattern)
  for (i in 1:(length(pattern_ppm)-1)){

    pattern_ppm_delta[i]=(pattern_ppm[i+1]-filtered_pattern[nrow(filtered_pattern),1])/filtered_pattern[nrow(filtered_pattern),1]

    if (pattern_ppm_delta[i]>(ppm/1000000)){
      #filtered_pattern<-rbind(filtered_pattern,pattern[i+1,])
      filtered_pattern[nrow(filtered_pattern)+1,]<-pattern[i+1,]
    } else {
      #previous_iso=filtered_pattern[nrow(filtered_pattern),]

      #newline=as.data.frame(list("m/z"=(filtered_pattern[nrow(filtered_pattern),1]*filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,1]*pattern[i+1,2])/(sum(filtered_pattern[nrow(filtered_pattern),2],pattern[i+1,2])),"abundance"=filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,2]))
      #names(newline)=c("m/z","abundance")
      newline=c((filtered_pattern[nrow(filtered_pattern),1]*filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,1]*pattern[i+1,2])/(sum(filtered_pattern[nrow(filtered_pattern),2],pattern[i+1,2])),
                filtered_pattern[nrow(filtered_pattern),2]+pattern[i+1,2])
      filtered_pattern[nrow(filtered_pattern),]<-newline
    }

  }

  filtered_pattern<-filtered_pattern[filtered_pattern[,2]>=max(filtered_pattern[,2])*threshold,]

  rownames(filtered_pattern)=1:nrow(filtered_pattern)
  if(verbose==T){
    message(paste("Origional Features:",org_feature,"Filtered Features:",nrow(filtered_pattern)))
  }


  return(filtered_pattern)
}

isopattern_ppm_filter_peaklist_par<-function(pattern,ppm,threshold=0.001,verbose=F,BPPARAM=bpparam()){
  library(BiocParallel)
  if (!require(OneR)) install.packages("OneR")
  library(OneR)
  org_feature=nrow(pattern)

  pattern<-as.data.frame(pattern)

  pattern_ppm=as.numeric(as.character(pattern[,1]))

  pattern_ppm_delta=numeric()


  filtered_pattern<-pattern[1,]

  mzrange<-range(pattern[,1])

  mzbinsecs<-ceiling((mzrange[2]-mzrange[1])/100)

  patternlist<-list()

  n_worker<-bpworkers(BPPARAM)

  pattern$mzbin<-(bin(pattern$mz, nbins = n_worker, method = "content"))

  for (mzbin in unique(pattern$mzbin)){

    patternlist[[mzbin]]<-pattern[pattern$mzbin==mzbin,]

    patternlist[[mzbin]]$mzbin<-NULL

  }

  filtered_pattern<-bplapply(patternlist,isopattern_ppm_filter_peaklist,ppm=ppm,BPPARAM=BPPARAM,verbose=F)

  filtered_pattern<-do.call(rbind,filtered_pattern)

  if(verbose==T){
    message(paste("Origional Features:",org_feature,"Filtered Features:",nrow(filtered_pattern),"m/z tolerance:"),ppm)
  }
  #filtered_pattern$mzbin<-(bin(filtered_pattern$mz, nbins = n_worker/2, method = "content"))

  return(filtered_pattern)
}

IMS_data_process_quant<-function (datafile, Peptide_Summary_searchlist, segmentation_num = 5,
                                  threshold = 0.1, ppm,  Segmentation=c("spatialKMeans","spatialShrunkenCentroids","Virtual_segmentation","none") , Smooth_range = 1,
                                  colorstyle = "Set1",
                                  Virtual_segmentation_rankfile = "Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy\\radius_rank.csv",
                                  PMFsearch = TRUE, rotate = NULL, BPPARAM = bpparam())
{
  library(data.table)
  library(Cardinal)
  library(RColorBrewer)
  if (!is.null(rotate)) {
    message("Found rotation info")
    rotatedegrees = sapply(datafile, function(x, df) {
      degree = df[df$filenames == x, "rotation"]
      if (length(degree) == 0) {
        message("Missing rotation data please check the rotation configuration file: ",
                x)
        degree = 0
      }
      degree
    }, rotate)
    rotate = unlist(rotatedegrees)
  }
  else {
    rotate = rep(0, length(datafile))
  }

  datafile_imzML<-datafile
  for (z in 1:length(datafile)){
    name <-gsub(base::dirname(datafile[z]),"",datafile[z])
    name <-gsub(".imzML$","",name)
    name <-gsub("/$","",name)
    folder<-base::dirname(datafile[z])
    #imdata <- Cardinal::readImzML(datafile[z],preprocessing = F,attach.only = T,resolution = 200,rotate = rotate[z],as="MSImageSet",BPPARAM = BPPARAM)
    if (!str_detect(datafile[z],".imzML$")){
      datafile_imzML[z]<-paste0(datafile[z],".imzML")
    }
    imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImageSet",resolution=200, units="ppm")


    #imdata <- Load_Cardinal_imaging(datafile[z], preprocessing = F,
    #                               attach.only = T, resolution = 200, rotate = rotate[z],
    #                                as = "MSImageSet", BPPARAM = BPPARAM)
    name <- gsub(base::dirname(datafile[z]), "", datafile[z])
    folder <- base::dirname(datafile[z])
    coordata = imdata@pixelData@data
    Peptide_Summary_file <- Peptide_Summary_searchlist
    Peptide_Summary_file$Intensity <- 0
    Peptide_Summary_file_regions <- NULL
    if (dir.exists(paste0(datafile[z], " ID")) == FALSE) {
      dir.create(paste0(datafile[z], " ID"))
    }
    setwd(paste0(datafile[z], " ID"))
    if (Segmentation=="spatialKMeans") {
      set.seed(1)
      message(paste0("spatialKMeans Clustering for ",
                     name))
      skm <- spatialKMeans(imdata, r = Smooth_range, k = segmentation_num,
                           method = "adaptive")
      message(paste0("spatialKMeans finished for ",
                     name))
      png(paste(getwd(), "/", "spatialKMeans_image_plot",
                ".png", sep = ""), width = 1024,
          height = 720)
      par(oma = c(0, 0, 0, 0), tcl = NA, mar = c(0, 0,
                                                 1, 1), mfrow = c(1, 2), bty = "n", pty = "s",
          xaxt = "n", yaxt = "n", no.readonly = TRUE,
          ann = FALSE)
      Cardinal::image(skm, col = brewer.pal_n(segmentation_num,
                                              colorstyle), key = FALSE, ann = FALSE, axes = FALSE)
      legend("topright", legend = 1:segmentation_num,
             fill = brewer.pal_n(segmentation_num, colorstyle),
             col = brewer.pal_n(segmentation_num, "Paired"),
             bg = "transparent", xpd = TRUE, cex = 1)
      Cardinal::plot(skm, col = brewer.pal_n(segmentation_num,
                                             colorstyle), type = c("p", "h"),
                     key = FALSE, mode = "withinss")
      legend("topright", legend = 1:segmentation_num,
             fill = brewer.pal_n(segmentation_num, colorstyle),
             col = brewer.pal_n(segmentation_num, "Paired"),
             bg = "transparent", xpd = TRUE, cex = 1)
      dev.off()
      withinss = skm@resultData[[1]][["withinss"]]
      centers = skm@resultData[[1]][["centers"]]
      cluster = as.data.frame(skm@resultData[[1]][["cluster"]])
      cluster$Coor = rownames(cluster)
      write.csv(withinss, paste("spatialKMeans_RESULT",
                                "withinss.csv"), row.names = F)
      write.csv(centers, paste("spatialKMeans_RESULT",
                               "centers.csv"), row.names = F)
      write.csv(cluster, paste("spatialKMeans_RESULT",
                               "cluster.csv"), row.names = F)
      y = skm@resultData[[paste0("r = ", Smooth_range,
                                 ", k = ", segmentation_num)]][["cluster"]]
      y = as.data.frame(y)
      rownames(y) = 1:nrow(y)
      y$pixel = 1:nrow(y)
      regions = unique(y[, 1])
      x = NULL
      for (i in 1:length(regions)) {
        listname = as.character(regions[i])
        x[[listname]] <- y[y[, 1] == regions[i], "pixel"]
      }
    }
    else if (Segmentation=="spatialShrunkenCentroids"){

      set.seed(1)
      message(paste0("spatialShrunkenCentroids Clustering for ",
                     name))
      skm <- spatialShrunkenCentroids(imdata, r = Smooth_range, k = segmentation_num,
                                      method = "adaptive",...)
      message(paste0("spatialShrunkenCentroids finished for ",
                     name))
      png(paste(getwd(), "/", "spatialShrunkenCentroids_image_plot",
                ".png", sep = ""), width = 1024,
          height = 720)
      par(oma = c(0, 0, 0, 0), tcl = NA, mar = c(0, 0,
                                                 1, 1), mfrow = c(1, 2), bty = "n", pty = "s",
          xaxt = "n", yaxt = "n", no.readonly = TRUE,
          ann = FALSE)
      Cardinal::image(skm, col = brewer.pal_n(segmentation_num,
                                              colorstyle), key = FALSE, ann = FALSE, axes = FALSE)
      legend("topright", legend = 1:segmentation_num,
             fill = brewer.pal_n(segmentation_num, colorstyle),
             col = brewer.pal_n(segmentation_num, "Paired"),
             bg = "transparent", xpd = TRUE, cex = 1)
      Cardinal::plot(skm, col = brewer.pal_n(segmentation_num,
                                             colorstyle), type = c("p", "h"),
                     key = FALSE, mode = "withinss")
      legend("topright", legend = 1:segmentation_num,
             fill = brewer.pal_n(segmentation_num, colorstyle),
             col = brewer.pal_n(segmentation_num, "Paired"),
             bg = "transparent", xpd = TRUE, cex = 1)
      dev.off()
      withinss = skm@resultData[[1]][["withinss"]]
      centers = skm@resultData[[1]][["centers"]]
      cluster = as.data.frame(skm@resultData[[1]][["cluster"]])
      cluster$Coor = rownames(cluster)
      write.csv(withinss, paste("spatialShrunkenCentroids_RESULT",
                                "withinss.csv"), row.names = F)
      write.csv(centers, paste("spatialShrunkenCentroids_RESULT",
                               "centers.csv"), row.names = F)
      write.csv(cluster, paste("spatialShrunkenCentroids_RESULT",
                               "cluster.csv"), row.names = F)
      y = skm@resultData[[paste0("r = ", Smooth_range,
                                 ", k = ", segmentation_num)]][["cluster"]]
      y = as.data.frame(y)
      rownames(y) = 1:nrow(y)
      y$pixel = 1:nrow(y)
      regions = unique(y[, 1])
      x = NULL
      for (i in 1:length(regions)) {
        listname = as.character(regions[i])
        x[[listname]] <- y[y[, 1] == regions[i], "pixel"]
      }
    }
    else if (Segmentation=="Virtual_segmentation") {
      radius_rank = read.csv(file = Virtual_segmentation_rankfile)
      radius_rank = radius_rank[order(radius_rank$Rank),]
      coordist_para = function(i, coordata) {
        coordist_para = NULL
        for (j in 1:nrow(coordata)) {
          coordist_para[j] = sqrt((coordata$x[i] - coordata$x[j])^2 +
                                    (coordata$y[i] - coordata$y[j])^2)
        }
        coordist_para
      }
      coordistmatrix = bplapply(1:nrow(coordata), coordist_para,
                                coordata, BPPARAM = BPPARAM)
      coordistmatrix = matrix(unlist(coordistmatrix), nrow = nrow(coordata),
                              ncol = nrow(coordata))
      coordistmatrix = as.data.table(coordistmatrix)
      coordistmatrix$sum = 0
      coordistmatrix$sum = base::unlist(bplapply(1:nrow(coordata),
                                                 function(j, coordistmatrix, coordata) {
                                                   coordistmatrix$sum[j] = sum(coordistmatrix[j,
                                                                                              1:nrow(coordata)])
                                                 }, coordistmatrix, coordata, BPPARAM = BPPARAM))
      coorrange = max(coordistmatrix$sum) - min(coordistmatrix$sum)
      findedge <- function(coordata) {
        uniquex = unique(coordata$x)
        uniquey = unique(coordata$y)
        coordata$edge = FALSE
        for (x in uniquex) {
          min = min(coordata[coordata$x == x, "y"])
          max = max(coordata[coordata$x == x, "y"])
          coordata[coordata$y == max & coordata$x ==
                     x, "edge"] = TRUE
          coordata[coordata$y == min & coordata$x ==
                     x, "edge"] = TRUE
        }
        for (y in uniquey) {
          min = min(coordata[coordata$y == y, "x"])
          max = max(coordata[coordata$y == y, "x"])
          coordata[coordata$x == max & coordata$y ==
                     y, "edge"] = TRUE
          coordata[coordata$x == min & coordata$y ==
                     y, "edge"] = TRUE
        }
        coordata
      }
      coordata = findedge(coordata)
      rank_pixel <- function(coordata, coordistmatrix) {
        shape_center = coordata[coordistmatrix$sum ==
                                  min(coordistmatrix$sum), ]
        center_dist = t(coordistmatrix[which.min(coordistmatrix$sum),
                                       1:nrow(coordata)])
        p_load(useful)
        From <- shape_center[rep(seq_len(nrow(shape_center)),
                                 each = nrow(coordata)), 1:2]
        To <- coordata[, 1:2]
        df = To - From
        center_edge_angle = cbind(coordata[, 1:2], cart2pol(df$x,
                                                            df$y, degrees = F), edge = coordata[, "edge"])
        center_edge_angle_sdge = center_edge_angle[center_edge_angle$edge ==
                                                     TRUE, ]
        coordata$rank = 0
        coordata$pattern = ""
        for (i in 1:(nrow(coordata))) {
          if (coordata$edge[i] != TRUE) {
            df = coordata[i, 1:2] - shape_center[, 1:2]
            point_center_angle = cbind(coordata[i, 1:2],
                                       cart2pol(df$x, df$y, degrees = F))
            pointedge = center_edge_angle_sdge[which(abs(center_edge_angle_sdge$theta -
                                                           point_center_angle$theta) == min(abs(center_edge_angle_sdge$theta -
                                                                                                  point_center_angle$theta))), ]
            pointedge = pointedge[which.min(pointedge$r),
            ]
            to_edge = coordistmatrix[[i]][coordata$x ==
                                            pointedge$x & coordata$y == pointedge$y]
          }
          else {
            to_edge = 0
          }
          to_center = center_dist[i]
          total = to_edge + to_center
          max(radius_rank$Radius_U)
          norm_center_dist = to_center/total * max(radius_rank$Radius_U)
          coordata$rank[i] = as.character(radius_rank$Rank[radius_rank$Radius_L <=
                                                             norm_center_dist & radius_rank$Radius_U >=
                                                             norm_center_dist])
          coordata$pattern[i] = as.character(radius_rank$Name[radius_rank$Radius_L <=
                                                                norm_center_dist & radius_rank$Radius_U >=
                                                                norm_center_dist])
        }
        coordata
      }
      coordata = rank_pixel(coordata, coordistmatrix)
      x = NULL
      for (rank in coordata$rank) {
        x[[rank]] = which(coordata$rank == rank)
      }
      write.csv(coordata, "coordata.csv", row.names = F)
      region_pattern <- factor(coordata$rank, levels = unique(coordata$rank),
                               labels = unique(coordata$pattern))
      set.seed(1)
      skm <- spatialKMeans(imdata, r = Smooth_range, k = length(unique(coordata$rank)),
                           method = "adaptive")
      png(paste(getwd(), "/", "spatialKMeans_image",
                ".png", sep = ""), width = 1024,
          height = 1024)
      par(oma = c(0, 0, 0, 0), tcl = NA, mar = c(0, 0,
                                                 1, 1), mfrow = c(1, 1), bty = "n", pty = "s",
          xaxt = "n", yaxt = "n", no.readonly = TRUE,
          ann = FALSE)
      Cardinal::image(skm, col = brewer.pal_n(segmentation_num,
                                              colorstyle), key = FALSE, ann = FALSE, axes = FALSE)
      legend("topright", legend = 1:segmentation_num,
             fill = brewer.pal_n(segmentation_num, colorstyle),
             col = brewer.pal_n(segmentation_num, "Paired"),
             bg = "transparent", xpd = TRUE, cex = 1)
      dev.off()
      for (i in 1:nrow(coordata)) {
        skm@resultData[[paste0("r = ", Smooth_range,
                               ", k = ", length(unique(coordata$rank)))]][["cluster"]][[rownames(coordata)[i]]] = factor(coordata[i,
                                                                                                                                  "rank"], levels = unique(coordata$rank))
      }
      png(paste(getwd(), "/", "Virtual_segmentation",
                gsub("/", " ", name), ".png",
                sep = ""), width = 1024, height = 1024)
      par(oma = c(0, 0, 0, 0), tcl = NA, mar = c(0, 0,
                                                 1, 1), mfrow = c(1, 1), bty = "n", pty = "s",
          xaxt = "n", yaxt = "n", no.readonly = TRUE,
          ann = FALSE)
      Cardinal::image(skm, col = brewer.pal_n(length(unique(coordata$rank)),
                                              colorstyle), key = FALSE, ann = FALSE, axes = FALSE)
      legend("topright", legend = paste(radius_rank$Rank,
                                        radius_rank$Name), fill = brewer.pal_n(length(unique(coordata$pattern)),
                                                                               colorstyle), col = brewer.pal_n(length(unique(coordata$pattern)),
                                                                                                               "Paired"), bg = "transparent", xpd = TRUE,
             cex = 1)
      dev.off()
    }
    else {
      x = 1:length(pixels(imdata))
      x = split(x, sort(x%%segmentation_num))
    }
    if (PMFsearch) {
      imdata <- Load_Cardinal_imaging(datafile[z], preprocessing = F,
                                      resolution = ppm, rotate = rotate[z], as = "MSImageSet")
      if (dir.exists(paste0(datafile[z], " ID")) ==
          FALSE) {
        dir.create(paste0(datafile[z], " ID"))
      }
      setwd(paste0(datafile[z], " ID"))
      message(paste("PMFsearch"))
      message(paste("region", names(x), sep = " ",
                    collapse = "\n"))
      for (SPECTRUM_batch in names(x)) {
        imdata_ed <- batchProcess(imdata, normalize = FALSE,
                                  smoothSignal = FALSE, reduceBaseline = FALSE,
                                  peakPick = list(SNR = 10), peakAlign = FALSE,
                                  pixel = pixels(imdata)[unlist(x[SPECTRUM_batch])])
        if (class(imdata_ed) == "MSImageSet") {
          savename = paste(name, SPECTRUM_batch)
          message(paste("IMS_analysis", name, "region",
                        SPECTRUM_batch))
          peaklist <- imdata_ed@featureData@data
          colnames(peaklist) <- c("m.z", "intensities")
          Peptide_Summary_searchlist <- IMS_analysis_fun(Peptide_Summary_searchlist = Peptide_Summary_searchlist,
                                                         peaklist = peaklist, ppm = ppm, BPPARAM = BPPARAM)
          Peptide_feature_list <- Peptide_Summary_searchlist[Peptide_Summary_searchlist$Intensity >
                                                               0, ]
          Peptide_plot_list <- Peptide_feature_list[Peptide_feature_list$Intensity >=
                                                      max(Peptide_feature_list$Intensity) * threshold,
          ]
          if (is.null(Peptide_plot_list$moleculeNames)) {
            Peptide_plot_list$moleculeNames = Peptide_plot_list$Peptide
          }
          Plot_PMF_all(Peptide_plot_list, peaklist, threshold = threshold,
                       savename)
          uniques_intensity <- unique(Peptide_plot_list[,
                                                        c("mz", "Intensity")])
          uniques_intensity <- uniques_intensity[uniques_intensity$Intensity >
                                                   0, ]
          Peptide_plot_list$Region = SPECTRUM_batch
          write.csv(Peptide_plot_list, paste("Peptide_segment_PMF_RESULT",
                                             SPECTRUM_batch, ".csv"), row.names = F)
          write.csv(peaklist, paste("Spectrum",
                                    SPECTRUM_batch, ".csv"), row.names = F)
          Peptide_Summary_file$Intensity <- Peptide_Summary_file$Intensity +
            unlist(bplapply(Peptide_Summary_file$mz,
                            intensity_sum_para, uniques_intensity,
                            BPPARAM = BPPARAM))
          Peptide_Summary_file_regions <- rbind(Peptide_Summary_file_regions,
                                                Peptide_plot_list)
        }
      }
      if (is.null(Peptide_Summary_file$Peptide)) {
        Peptide_Summary_file$Peptide = Peptide_Summary_file$moleculeNames
      }
      if (is.null(Peptide_Summary_file$moleculeNames)) {
        Peptide_Summary_file$moleculeNames = Peptide_Summary_file$Peptide
      }
      Peptide_Summary_file <- Peptide_Summary_file[Peptide_Summary_file$Intensity >
                                                     0, ]
      write.csv(Peptide_Summary_file, "Peptide_Summary_file.csv",
                row.names = F)
      write.csv(Peptide_Summary_file_regions, "Peptide_region_file.csv",
                row.names = F)
    }
  }
}

radius_segmentation<-function(i,z_dist,radius_rank){
  radius_rank['&'(z_dist[i] >= radius_rank$Radius_L, z_dist[i] < radius_rank$Radius_U),"Name"]
}

radius_segmentation_dist<-function(i,z_dist,radius_rank){
  radius_rank['&'(z_dist[i] >= radius_rank$Radius_L, z_dist[i] < radius_rank$Radius_U),"Name"]
}

Center_of_gravity_and_contour<-function(imdata,BPPARAM=NULL){
  library(BiocParallel)
  library(plotly)
  library(SDMTools)
  if(missing(BPPARAM)){
    BPPARAM=bpparam()
  }
  #radius_rank=data.frame(matrix(c(1:4),c("Core","Inner","Barrier","Outer"),c(0,4.5,6,8),ncol=3,nrow=4))
  COG=SDMTools::COGravity(coord(imdata)$x,coord(imdata)$y,coord(imdata)$z,wt=coord(imdata)$z)

  z=sqrt((coord(imdata)$y-COG[3])^2*COG[2] +(coord(imdata)$x-COG[1])^2*COG[4])
  z_dist=0.15*sqrt(((coord(imdata)$y-COG[3]))^2 *COG[2]/COG[4]+((coord(imdata)$x-COG[1]))^2)
  radius_rank=data.frame("Rank"=1:4,"Name"=c("Core","Inner","Barrier","Outer"),"Radius_L"=c(0,2.25,3,4),"Radius_U"=c(2.25,3,4,max(z_dist)))

  #z_radius_segmentation=as.data.frame(unlist(parallel::parLapply(cl=autoStopCluster(makeCluster(detectCores())),1:length(z_dist),fun = radius_segmentation,z_dist,radius_rank)))
  z_radius_segmentation=as.data.frame(unlist(bplapply(1:length(z_dist),fun = radius_segmentation,z_dist,radius_rank,BPPARAM = BPPARAM)))

  #z_radius_segmentation[i]=radius_rank['&'(z_dist[i] >= radius_rank$Radius_L, z_dist[i] < radius_rank$Radius_U),"Name"]

  p <- plot_ly(
    x = coord(imdata)$x,
    y = coord(imdata)$y,
    z = z,
    type = "contour"
  )

  p <- plot_ly(
    x=coord(imdata)$x,
    y=coord(imdata)$y,
    type = 'contour',
    z = z,
    colorscale = 'Jet',
    autocontour = F,
    contours = list(
      start = 0,
      end = max(z),
      size = 2
    )
  )
  p

}
