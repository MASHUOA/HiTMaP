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
  
  
  
  #paste0("v",colnames(coordistmatrix[coordistmatrix$sum==min(coordistmatrix$sum),]))
  
  #center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),])
  
  
  #plot(rownames(coordistmatrix),coordistmatrix$sum)
  
  #write.csv(radius_rank,file = "radius_rank.csv", row.names = F)
  
  
  
  rank_pixel<-function(coordata,coordistmatrix){
    #coordata[coordata$edge==TRUE,]=coordata[coordata$edge==TRUE,]
    shape_center=coordata[coordistmatrix$sum==min(coordistmatrix$sum),]
    center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),1:nrow(coordata)])
    library(useful)
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
  
  coordata
  
  }
  
  
  coordatalist_result<-lapply(coordatalist,coordata_rank_list,radius_rank)
  
  coordatalist_result_merge<-do.call(rbind,coordatalist_result)
  coordatadf_merge<-merge(coordatadf,coordatalist_result_merge[,c("x","y","run","rank","pattern")],by=c("x","y","run"))
  rownames(coordatadf_merge)<-coordatadf_merge$indx
  return(coordatadf_merge)
  
  }



Load_Cardinal_imaging<-function(datafile=tk_choose.files(filter = Filters,
                                                         caption  = "Choose single or multiple file(s) for analysis",
                                                         multi = F),
                                BPPARAM = bpparam(),
                                resolution=2.5,
                                attach.only=TRUE,
                                preprocessing=F,
                                rotate=0,as="MSImagingExperiment",
                                mzrange=NULL,is_centroided=F){
  #imdata_meta <- importImzMl(datafile, coordinates = matrix(c(2, 4),nrow=1, ncol=2),removeEmptySpectra = F, centroided = T)
  suppressMessages(suppressWarnings(library(stringr)))
  suppressMessages(suppressWarnings(library(BiocParallel)))
  datafile<-gsub("\\\\", "/", datafile)
  datafiles<-str_remove(datafile,regex("\\.imzML$", ignore_case = T))
  workdir<-base::dirname(datafiles) 
  name <-basename(datafiles)
  folder<-base::dirname(datafiles)
  if (rotate==0){
    imdata <-  suppressMessages(suppressWarnings(Cardinal::readImzML(name, folder, attach.only=attach.only,as=as,resolution=resolution, units="ppm",BPPARAM=BPPARAM,mass.range=mzrange,is_centroided=is_centroided)))
    
  }else  {
    imdata <-  suppressMessages(suppressWarnings(Cardinal::readImzML(name, folder, attach.only=attach.only,as=as,resolution=resolution, units="ppm",rotate = rotate,BPPARAM=BPPARAM,mass.range=mzrange,is_centroided=is_centroided)))
    imdata <-  rotateMSI(imdata=imdata,rotation_degree=rotate)
    }
  
  if (preprocessing){
    imdata <- try(batchProcess(imdata, normalize=FALSE, smoothSignal=TRUE, reduceBaseline=list(method = "median",blocks=500, fun=min, spar=1),
                               peakPick=list(SNR=12), peakAlign=TRUE,BPPARAM=BPPARAM))
  }
  imdata@centroided=is_centroided
  imdata
  
}

windows_filename<- function(stringX){
  stringX<-str_remove_all(stringX,"[><*?:\\/\\\\|]")
  stringX<-gsub("\"", "", stringX)
  return(stringX)
  
}
windows_filename2<- function(stringX){
  stringX<-str_replace_all(stringX, "[^[:alnum:]]", "")
  return(stringX)
  
}

Advanced_building<-function(){
  install.packages("ggplot2")
  install.packages("pryr")
  install.packages("devtools")
  devtools::install_github("hadley/lineprof")
  library(pryr)
  object_size(1:10)
  object_size(mean)
  sizes <- sapply(0:50, function(n) object_size(seq_len(n)))
  plot(0:50, sizes, xlab = "Length", ylab = "Size (bytes)", 
       type = "s")
  object_size(numeric())
  plot(0:50, sizes - 40, xlab = "Length", 
       ylab = "Bytes excluding overhead", type = "n")
  abline(h = 0, col = "grey80")
  abline(h = c(8, 16, 32, 48, 64, 128), col = "grey80")
  abline(a = 0, b = 4, col = "grey90", lwd = 4)
  lines(sizes - 40, type = "s")
  mem_used()
  
  read_delim <- function(file, header = TRUE, sep = ",") {
    # Determine number of fields by reading first line
    first <- scan(file, what = character(1), nlines = 1,
                  sep = sep, quiet = TRUE)
    p <- length(first)
    
    # Load all fields as character vectors
    all <- scan(file, what = as.list(rep("character", p)),
                sep = sep, skip = if (header) 1 else 0, quiet = TRUE)
    
    # Convert from strings to appropriate types (never to factors)
    all[] <- lapply(all, type.convert, as.is = TRUE)
    
    # Set column names
    if (header) {
      names(all) <- first
    } else {
      names(all) <- paste0("V", seq_along(all))
    }
    
    # Convert list into data frame
    as.data.frame(all)
  }
  
  library(ggplot2)
  write.csv(diamonds, "diamonds.csv", row.names = FALSE)
  
  library(lineprof)
  
  source("Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy/R/read-delim.R")
  prof <- lineprof(read_delim("diamonds.csv"))
  shine(prof)
}
string_slice <- function(string, size) {
  pat <- paste0('(?<=.{',size,'})')
  strsplit(string, pat, perl=TRUE)
}

initialization<-function(){
  #.libPaths("~/R/win-library/peptidegidest353")
  library("pacman")
  #p_unlock()
  p_unload("all",negate = T)
  library("pacman")
  library("devtools")
  
  #devtools::install_github("tidyverse")
  
 if (!requireNamespace("BiocManager")) install.packages("BiocManager")
  BiocManager::install(c("purrr","RColorBrewer","RCurl","bitops","magick","ggplot2",
                         "reticulate","dplyr","stringr","data.table","iterators","foreach",
                         "protViz","cleaver","MALDIquant","Biostrings","XVector","IRanges","Cardinal","ProtGenerics",
                         "S4Vectors","EBImage","BiocParallel","BiocGenerics",
                         "Rdisop","Rcpp"))
  p_load(purrr,fs,processx,RColorBrewer,RCurl,bitops,magick,ggplot2,Rdisop,Rcpp)
  p_load(reticulate,dplyr,stringr,tcltk,data.table,doParallel,iterators,foreach)
  p_load(protViz,cleaver,MALDIquant,Biostrings,XVector,IRanges,Cardinal,ProtGenerics)
  p_load(S4Vectors,stats4,EBImage,BiocParallel,BiocGenerics,parallel,stats,graphics)
  p_load(grDevices,utils,datasets,methods)
}

getMonomass <- function(formular){
  Rdisop::getMolecule(formular)[["exactmass"]][1]}
getMonomass_para <- function(i,formularlist){
  return(Rdisop::getMolecule(formularlist[i])$exactmass)}
setworkdir<-function(workdir){
  if (dir.exists(workdir)==FALSE){dir.create(workdir)}
  setwd(workdir)
}

cluster_image_cardinal_allinone<-function(clusterID,
                                          SMPLIST,
                                          ppm=20,
                                          imdata=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm),
                                          ClusterID_colname="Protein",
                                          componentID_colname="Peptide",
                                          Component_plot_threshold=2,
                                          smooth.image="gaussian",
                                          contrast.enhance = "suppression",
                                          colorpallet="Set1",
                                          plot_layout="grid",
                                          image_rotation_degree=NULL,
                                          export_Header_table=F,
                                          rotate_image=F,
                                          plot_style=c("fleximaging","ClusterOnly","rainbow")){
  #complementary(color="red", plot = TRUE, bg = "white", labcol = NULL, cex = 0.8, title = TRUE)
  windows_filename<- function(stringX){
    stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
    stringX<-gsub("\"", "", stringX)
    stringX<-gsub("\\|", " ", stringX)
    return(stringX)
    
  }
  Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
  library(grid)
  library(plotly)
  #rotate the image
  #imdata@pixelData@data<-rotatetmp

  
  SMPLIST=as.data.frame(SMPLIST)
  outputpng=paste(getwd(),"\\",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot",'.png',sep="")  
  #message(outputpng)
  candidate=SMPLIST[SMPLIST[[ClusterID_colname]]==clusterID,]
  #candidate=candidate[order(as.character())]
  candidateunique=as.numeric(as.character(unique(candidate[,"mz"])))
  if (length(candidateunique)>9){
    candidate.dt <- data.table(candidate)
    candidatet=candidate.dt[,list(Intensity=sum(Intensity)), by='mz']
    candidatet=candidatet[order(-candidatet$Intensity)]
    selections=as.numeric(t(candidatet[1:9,"mz"]))
    candidate=candidate[candidate$mz %in% selections,]
    candidateunique=as.numeric(unique(candidate[,"mz"]))
    candidateunique=candidateunique[order(as.character(candidateunique))]
    mycol <- factor(RColorBrewer::brewer.pal(length(candidateunique),colorpallet))
    mycol <- factor(mycol,levels(mycol)) 
  }else if (length(candidateunique)<3){
  candidateunique=candidateunique[order(as.character(candidateunique))]
  mycol <- factor(RColorBrewer::brewer.pal(3,colorpallet))
  mycol <- factor(mycol,levels(mycol)) 
  
  }else{
    candidateunique=candidateunique[order(as.character(candidateunique))]
    mycol <- factor(RColorBrewer::brewer.pal(length(candidateunique),colorpallet))
    mycol <- factor(mycol,levels(mycol)) 
  }
  
  
  
  if (length(candidateunique)>=Component_plot_threshold){
    
  if (is.null(imdata)){
    message("No imaging data")
    
  }else{
    
    
    
    library(RColorBrewer)
    library(Cardinal)
    library(EBImage)
    #library(colortools)
    
    if (plot_style=="ClusterOnly"){
    png(outputpng,width = 10,height = 10, bg = "black",units = "in",res = 300)
    par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,1),
        bty="n",pty="s",xaxt="n",
        yaxt="n",
        no.readonly = TRUE,ann=FALSE)
    image(imdata, mz=candidateunique, 
          col=mycol,
          contrast.enhance = contrast.enhance,
          smooth.image = smooth.image ,
          superpose=TRUE,normalize.image="linear",
          plusminus=median(ppm*candidateunique/1000000))
    
    
    l<-function(x,y,z,t=x+y,f=t+z){
      paste(t,f)
    }
    l(1,2,3)
    
    dev.off()
    pngfile<-image_read(outputpng)
    pngfile<-image_border(pngfile, "black", "30x30")
    pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
    pngfile<-image_trim(pngfile)
    image_write(pngfile,outputpng)
    
    }else if (plot_style=="rainbow"){
    
    
    
    if (plot_layout=="line"){
      png(outputpngsum,width = 5*((length(candidateunique)+1)),height = 5, bg = "black",units = "in",res = 75)
      par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
          bty="n",pty="s",xaxt="n",
          yaxt="n",
          no.readonly = TRUE,ann=FALSE)   
      
      
    }else{
    png(outputpngsum,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 75)
    par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
        bty="n",pty="s",xaxt="n",
        yaxt="n",
        no.readonly = TRUE,ann=FALSE)  
    }

    
    
    
    image(imdata, mz=candidateunique, 
          col=levels(mycol),
          contrast.enhance = contrast.enhance,
          smooth.image = smooth.image ,
          superpose=TRUE,normalize.image="linear",
          plusminus=median(ppm*candidateunique/1000000))
    
    for (i in 1:length(candidateunique)){
      #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
      col.regions <- gradient.colors(100, start="black", end=levels(mycol)[i])
      image(imdata, mz=candidateunique[i], 
            contrast.enhance=contrast.enhance,
            smooth.image = smooth.image ,
            col.regions=col.regions,
            
            normalize.image="none",
            plusminus=ppm*candidateunique[i]/1000000)
      componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
      for (component in componentname){
        text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
      }
    }
    dev.off()
    #pngfile<-image_read(outputpngsum)
    #pngfile<-image_border(pngfile, "black", "30x30")
    #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
    #pngfile<-image_trim(pngfile)
    #image_write(pngfile,outputpngsum)
    
    
    }else if(plot_style=="fleximaging"){
    
    ##################################
    
    
    
    if (plot_layout=="line"){
      png(outputpngsum,width = 5*((length(candidateunique)+1)),height = 5, bg = "black",units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
          bty="n",pty="s",xaxt="n",
          yaxt="n",
          no.readonly = TRUE,ann=FALSE)   
      
      
    }else{
      png(outputpngsum,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
          bty="n",pty="s",xaxt="n",
          yaxt="n",
          no.readonly = TRUE,ann=FALSE)  
    }
    
    
    
    
    image(imdata, mz=candidateunique, 
          col=levels(mycol),
          contrast.enhance = contrast.enhance,
          smooth.image = smooth.image ,
          superpose=TRUE,normalize.image="linear",
          plusminus=median(ppm*candidateunique/1000000))
    

    for (i in 1:length(candidateunique)){
      #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
      col.regions <- gradient.colors(100, start="black", end=levels(mycol)[i])
      image(imdata, mz=candidateunique[i], 
            contrast.enhance=contrast.enhance,
            smooth.image = smooth.image,#smooth.image ,
            col.regions=intensity.colors_customize1(),
            normalize.image="none",
            plusminus=ppm*candidateunique[i]/1000000,
            key=F,
            xlab=NULL,
            ylab=NULL,
            )
      componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
      for (component in componentname){
        text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
      }
    }
    dev.off()
    pngfile<-image_read(outputpngsum)
    pngfile<-image_border(pngfile, "black", "30x30")
    pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
    pngfile<-image_trim(pngfile)
    image_write(pngfile,outputpngsum)
    }
  }
    
    
    
    
  }
  
  if(export_Header_table){
    library(gridExtra)
    library(grid)
    candidate_unique_table=unique(candidate[,c("mz",ClusterID_colname,"formula","moleculeNames" , "adduct")])
    Header_table<-NULL
    Header_table$mz=candidateunique
    Header_table<-data.frame(Header_table)
    Header_table=base::merge(Header_table,candidate_unique_table,all.x=T,by="mz",all.y=F)
    #Header_table$ID=candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique)]
    #componentnames=unique(Header_table[[componentID_colname]][Header_table$mz==as.numeric(candidateunique[i])])
    p <- plot_ly(
      type = 'table',
      columnwidth = 20,
      header = list(
        values = c(paste0("<b>","Cluster","</b>"),Header_table$FA),
        colspan = I(30),
        align = c('center'),
        line = list(width = 1, color = 'black'),
        fill = list(color = "grey"),
        font = list(family = "Arial", size = 15, color = "white")
      ),
      cells = list(
        values = cbind(
          #rbind(paste0("<b>","Cluster","</b>"), as.matrix(Header_table$FA)),
          rbind("",as.matrix(Header_table$adduct)), 
          rbind("",as.matrix(round(Header_table$mz,digits = 3)))
        ),
        align = c('left', rep('center', ncol(mtcars))),
        line = list(color = "black", width = 1),
        fill = list(color = "white"),
        font = list(family = "Arial", size = 15, color = c("black"))
      
      ),
      width=120*(nrow(Header_table)+1),
      height=400,
      
      
      ) %>% layout(
        title = paste0("<b>",clusterID,"</b>"),
        autosize = T,
        margin=0,
        font = list(family = "Arial", size = 20, color = "black",align = "bottom")
        )
      
      
      
      
      
    orca(p, file = windows_filename(paste0(clusterID,"header.png")),width=120*(nrow(Header_table)+1),height=540) 
    
  }
  
  
  
  
}

#' cluster_image_grid
#'
#' This function renders the clustered images for maldi imaging data set
#'
#' @export

cluster_image_grid<-function(clusterID,
                             SMPLIST,
                             imdata,
                             ClusterID_colname="Protein",
                             componentID_colname="Peptide",
                             Component_plot_threshold=2,
                             Component_plot_coloure=c("mono","as.cluster"),
                             smooth.image="gaussian",
                             contrast.enhance = "suppression",
                             colorpallet="Set1",
                             plot_layout=c("line","grid"),
                             export_Header_table=F,
                             export_footer_table=F,
                             plot_style=c("fleximaging","ClusterOnly","rainbow"),
                             protein_coverage=F,
                             footer_style="Length",
                             output_png_width_limit=1980,
                             attach_summary_cluster=T,
                             cluster_color_scale=c("blackwhite","fleximaging"),
                             remove_cluster_from_grid=T,
                             img_brightness= 100,ppm=20,
                             list_of_protein_sequence,
                             workdir=getwd(),
                             pixel_size_um=50,
                             Score_thres=NULL){
  #complementary(color="red", plot = TRUE, bg = "white", labcol = NULL, cex = 0.8, title = TRUE)
  windows_filename<- function(stringX){
    stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
    stringX<-gsub("\"", "", stringX)
    stringX<-gsub("\\|", " ", stringX)
    return(stringX)
    
  }
   suppressMessages(suppressWarnings(require(magick)))
   suppressMessages(suppressWarnings(require(stringr)))
   suppressMessages(suppressWarnings(require(data.table)))
   #Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
   suppressMessages(suppressWarnings(require(grid)))
   suppressMessages(suppressWarnings(require(plotly)))
   suppressMessages(suppressWarnings(require(dplyr)))
   suppressMessages(suppressWarnings(require(colortools)))
   suppressMessages(suppressWarnings(require(data.table)))
   suppressMessages(suppressWarnings(require(Cardinal)))
  #rotate the image
  #imdata@pixelData@data<-rotatetmp
  outputpngsum=paste(workdir,"/",windows_filename(substr(clusterID, 1, 15)),"_cluster_imaging",'.png',sep="")
  
  #message(outputpngsum)
  SMPLIST=as.data.frame(SMPLIST)
  outputpng=paste(workdir,"/",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot",'.png',sep="")  
  #message(outputpng)
  candidate=SMPLIST[SMPLIST[[ClusterID_colname]]==clusterID,]
  
  if(!is.null(Score_thres)){candidate=candidate[candidate$Score>=Score_thres,]}
  
  #candidate=candidate[order(as.character())]
  if (componentID_colname %in% colnames(candidate)){
  candidate_u<- candidate %>% dplyr::group_by(mz) %>% dplyr::summarise(Peptide=Peptide[1], .groups = 'drop')

  candidate<-merge(data.frame(candidate_u),candidate,by=c("mz",componentID_colname),sort=F)
  }
  candidateunique=as.numeric(as.character(unique(candidate[,"mz"])))
  candidateunique=candidateunique[order(as.character(candidateunique))]
  
  #check if candidates' mz within the dataset range
  candidateunique<-candidateunique[which(between(candidateunique,range(mz(imdata))[1],range(mz(imdata))[2]))]
  candidate<-candidate[candidate$mz %in% candidateunique,]
  
  if(is.null(candidate$desc)){
    candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
    cluster_desc<-"-"
  }else{
    candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
    cluster_desc<-unique(candidate$desc)[1]
  } 
  if (length(candidateunique)>4){
    
    mycol=setColors("red", length(candidateunique))#wheel("red", num = length(candidateunique),bg = "white")
  } else if (length(candidateunique)==4){
    tmp_cols = setColors("red", 12)
    tetrad_colors <- tmp_cols[c(1, 3, 7, 9)]
    mycol=tetrad_colors
    #mycol=c("#FF0000", "#00FF00", "#0000FF")
  } else if (length(candidateunique)==3){
    #mycol=splitComp("red")
    mycol=c("#FF0000", "#00FF00", "#0000FF")
  } else if (length(candidateunique)==2){
    #mycol=complementary("red")
    mycol=c("#FF0000", "#00FFFF")
  } else if (length(candidateunique)==1){
    #mycol=complementary("red")
    mycol=c("#FF0000")
  }
    
  mycol <- as.factor(as.character(mycol))  
  mycol=mycol[order(mycol)]
  #mycolrgb<-hex2RGB(mycol)
  #meancoldf<-colSums(mycolrgb@coords)/max(colSums(mycolrgb@coords))
  #meancol<-RGB(R=meancoldf[1],G=meancoldf[2],B=meancoldf[3])
  #colorspace::hex(meancol)
  #print(length(candidateunique))
  #print(Component_plot_threshold)
  #print(length(candidateunique)>=Component_plot_threshold)
  if (length(candidateunique)>=Component_plot_threshold){
    
    if (is.null(imdata)){
      message("No imaging data")
      
    }else
      {
       suppressMessages(suppressWarnings(require(RColorBrewer)))
       suppressMessages(suppressWarnings(require(Cardinal)))
       suppressMessages(suppressWarnings(require(EBImage)))
      #library(colortools)
      
      if (plot_style=="ClusterOnly"){
        png(outputpng,width = 10,height = 10, bg = "black",units = "in",res = 300)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,1),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)
        image(imdata, mz=candidateunique, 
              col=mycol,
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              plusminus=median(ppm*candidateunique/1000000))
        
        
        l<-function(x,y,z,t=x+y,f=t+z){
          paste(t,f)
        }
        l(1,2,3)
        
        dev.off()
        pngfile<-image_read(outputpng)
        pngfile<-image_border(pngfile, "black", "30x30")
        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        pngfile<-image_trim(pngfile)
        image_write(pngfile,outputpng)
        
      }
      else if (plot_style=="rainbow"){
        
        
        if (plot_layout=="line"){
          png(outputpngsum,width = 5*((length(candidateunique)+1)),height = 5, bg = "black",units = "in",res = 75)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)   
          
          
        }else{
          png(outputpngsum,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 75)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)  
        }
        
        
        
        
        image(imdata, mz=candidateunique, 
              col=levels(mycol),
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              plusminus=round(median(ppm*candidateunique/1000000),digits = 4),key=F)
        
        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- gradient.colors(100, start="black", end=levels(mycol)[i])
          image(imdata, mz=candidateunique[i], 
                contrast.enhance=contrast.enhance,
                smooth.image = smooth.image ,
                col.regions=col.regions,
                
                normalize.image="none",
                plusminus=ppm*candidateunique[i]/1000000)
          componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          for (component in componentname){
            text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          }
        }
        dev.off()
        #pngfile<-image_read(outputpngsum)
        #pngfile<-image_border(pngfile, "black", "30x30")
        #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        #pngfile<-image_trim(pngfile)
        #image_write(pngfile,outputpngsum)
        
        
      }
      else if(plot_style=="fleximaging"){
        
        ##################################
        lightmode()
        tmp_dir <- tempdir()
        temp_cluster_png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        temp_component_png=list()
        temp_component_png.mono=list()
        if (plot_layout=="line"){
          #png(temp_cluster_png,width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = "black",units = "in",res = 300)
          temp_cluster_png<-image_graph(width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = "black",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)   
          
          
        }else{
          #png(temp_cluster_png,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 300)
          temp_cluster_png<-image_graph(width = 750*2,height = 750 * (ceiling((length(candidateunique)+1)/2)), bg = "black",res = 150)
          
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)  
        }
        
        
        col.regions <- gradient.colors(100, start="black", end=levels(mycol)[1])
        

        
        dev.off()
        #if (length(candidateunique)==2){
        #  bg = "grey9"
        #}else if (length(candidateunique)==3){
        #  bg = "grey8"
        #}else if (length(candidateunique)==3){
        #  
        #}
        bg = paste0("grey",ifelse((15-length(candidateunique))>1,(15-length(candidateunique)),2))
        #bg = "transparent"
        componentimg=list()
        componentimg.mono=list()
        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- gradient.colors(100, start="black", end=levels(mycol)[i])
          if  (Component_plot_coloure=="mono"){
            bg = paste0("grey",29)
            col.regions.mono=intensity.colors_customize1(colset = 2)
            col.regions=gradient.colors(100, start="black", end=levels(mycol)[i])
          }else if(Component_plot_coloure=="as.cluster"){
            col.regions=gradient.colors(100, start="black", end=levels(mycol)[i])
            #col.regions=gradient.colors(100, start="#030303", end=levels(mycol)[i])
            col=levels(mycol)[i]
          }
          
          if(cluster_color_scale=="blackwhite"){
            col.regions.mono=gradient.colors(100, start="black", end="white")
            col.regions=gradient.colors(100, start="black", end="white")
          }
          if (Component_plot_coloure=="mono"){
            componentimg.mono[[i]]=image(imdata, mz=candidateunique[i], 
                                    #contrast.enhance=contrast.enhance,
                                    contrast.enhance = "none",
                                    smooth.image = smooth.image,
                                    colorscale=col.regions.mono,
                                    col=col,
                                    normalize.image="linear",
                                    plusminus=round(ppm*candidateunique[i]/1000000*2,digits = 5),
                                    key=F,
                                    xlab=NULL,
                                    layout=c( length(levels(Cardinal::run(imdata))),1),
                                    bg = bg)
            #temp_component_png[[i]]=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
            componentimg.mono[[i]][["par"]][["ann"]]=F
            componentimg.mono[[i]][["par"]][["bty"]]="n"
            componentimg.mono[[i]][["par"]][["pty"]]="s"
            componentimg.mono[[i]][["par"]][["xaxt"]]="n"
            componentimg.mono[[i]][["par"]][["yaxt"]]="n"
            componentimg.mono[[i]][["par"]][["fg"]]="white"
            componentimg.mono[[i]][["par"]][["oma"]]=c(0, 0, 0, 0)
            #componentimg.mono[[i]][["par"]][["mar"]]=c(0, 0, 0, 1)
            attr(componentimg.mono[[i]][["facets"]][[1]],"strip")$strip=F
            #attr(componentimg.mono[[i]][["facets"]][[1]],"colorkey")$colorkey=c("0%","100%")
            attr(componentimg.mono[[i]][["facets"]][[1]],"colorkey")$colorkey=T
            
            # componentimg.mono[[i]][["facets"]][[1]][[1]][["dpage"]]="a"
            # componentimg.mono[[i]][["facets"]][[1]][[1]][["facet"]][[".feature.groups"]]="a"
            # componentimg.mono[[i]][["fids"]][[".feature.groups"]]="a"
            # componentimg.mono[[i]][["dpages"]]="a"
            #png(temp_component_png.mono[[i]],width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = bg,units = "in",res = 300)
            temp_component_png[[i]]<-image_graph(width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = bg,res = 150)
            par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
                #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
                bty="n",pty="s",xaxt="n",
                yaxt="n",
                no.readonly = TRUE,ann=FALSE)  
            
            try(tryCatch(print(componentimg.mono[[i]])),silent = T)
            dev.off()
          }else{
        componentimg[[i]]=image(imdata, mz=candidateunique[i], 
                                      #contrast.enhance=contrast.enhance,
                                      contrast.enhance = "none",
                                      smooth.image = smooth.image,
                                      colorscale=col.regions,
                                      col=col,
                                      normalize.image="linear",
                                      plusminus=round(ppm*candidateunique[i]/1000000*2,digits = 5),
                                      key=F,
                                      xlab=NULL,
                                      layout=c( length(levels(Cardinal::run(imdata))),1)
                                      )
        
        temp_component_png[[i]]=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

        #png(temp_component_png[[i]],width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = bg,units = "in",res = 300)
        temp_component_png[[i]]<-image_graph(width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = bg,res = 150)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
            #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)  

        try(tryCatch(print(componentimg[[i]])),silent = T)
        
        dev.off()
          }
          

        

          #componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          #for (component in componentname){
          #  text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          #}
        }
        
        temp_component_cover=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        
        #png(temp_component_cover,width = 5,height = 0.5, bg = bg,units = "in",res = 300)
        temp_component_cover<-image_graph(width = 750,height = 750, bg = bg,res = 150)
        plot.new()
        text("")
        #temp_component_cover<-image_draw(temp_component_cover)
        dev.off()
        temp_component_cover<-image_crop(temp_component_cover, "750x80+0")
        #temp_component_cover<-image_background(temp_component_cover, "hotpink")
        makeacover<-(temp_component_cover)
        #print(makeacover)
        #pngfile<-(temp_cluster_png)
        #pngfile<-image_border(pngfile, bg, "30x30")
        #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "north",size = 50,color = "white")
        
        pngcompfile_org=list()
        pngcompfile=list()
        for (i in 1:length(candidateunique)){
          pngcompfile_org[[i]]<-((temp_component_png[[i]]))
          pngcompfile_org[[i]]<-image_flatten(c(pngcompfile_org[[i]],makeacover))
          #if (Component_plot_coloure=="mono"){
          #pngcompfile[[i]]<-(temp_component_png.mono[[i]] ) 
          #pngcompfile[[i]]<-image_flatten(c(pngcompfile[[i]],makeacover))
          #}else{
          pngcompfile[[i]]<-pngcompfile_org[[i]] 
          #}
          
          pngcompfile[[i]]<-image_border(pngcompfile[[i]], bg, "30x30")
          pngcompfile[[i]]<-image_annotate(pngcompfile[[i]],paste(unique(candidate[candidate$mz==candidateunique[i],"moleculeNames"]),candidateunique[i]),gravity = "north",size = 30,color = "white")
          
        }
        

        
        
        pngcompfile_output<-pngcompfile[[1]]
        if (length(pngcompfile)>=2){
        for (i in 2: length(pngcompfile)){
          pngcompfile_output<-c(pngcompfile_output,unlist(pngcompfile[[i]]))
        } }
        
        #img_com<-(temp_component_png[[1]])
        #img_com<-image_flatten(c(img_com,makeacover))
        #pngcompfile_org[[1]]<-NULL
        img_com<-pngcompfile_org[[1]]
        if (length(pngcompfile)>=2){
        for (i in 2: length(pngcompfile_org)){
          img_com<-c(img_com,unlist(pngcompfile_org[[i]]))
        }}
        #img_com<-rbind(unlist(pngcompfile_org))
        
        
        #lapply(channel_types(),function(x,pngfile){
        #  pngfile %>% image_threshold(type = "white", threshold = "50%",channel = x) %>% image_scale( "x600") %>% image_write(path = paste0(x,".png"))
        #},pngfile)
        cluster_desc<-gsub(stringr::str_extract(cluster_desc,"OS=.{1,}"),"",cluster_desc)
        
       if ((Component_plot_coloure=="mono")) {
         
         pngfile_big<-image_average(img_com)
         pngfile_big_info=magick::image_info(pngfile_big)
         
         pngfile_big_ratio<-c(diff(range(imdata@elementMetadata@coord@listData[["x"]])),diff(range(imdata@elementMetadata@coord@listData[["y"]])))/c(pngfile_big_info$width-1*150,pngfile_big_info$height-1*150)
         

         
         pngfile_big <- image_draw(pngfile_big)
         #rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
         pixel_size_um_bar<-200
         pixel_size<-pixel_size_um_bar*pixel_size_um/(1/max(pngfile_big_ratio))/1000
         if(!is.integer(pixel_size)){
           pixel_size<-ceiling(pixel_size)
           pixel_size_um_bar<-pixel_size*1000/pixel_size_um*(1/max(pngfile_big_ratio))
         }
         
         rect(690-pixel_size_um_bar-60, 680,690-pixel_size_um_bar-60+pixel_size_um_bar,  680+13, border  = 'white',col  = 'white', lwd = 0)
         text(690, 690, paste0(pixel_size ," mm"),col  = 'white', cex = 3, srt = 0)
         #palette(rainbow(11, end = 0.9))
         #symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
         #        bg = 1:11, inches = FALSE, add = TRUE)
         dev.off()
         
         
         pngfile<-pngfile_big
         pngfile_big<-image_border(pngfile_big, bg, "30x30")
         #pngfile_big<-image_modulate(pngfile_big, brightness = 100 + 25 * length(candidateunique), saturation = 100)
         #pngfile_big<-image_threshold(pngfile_big,type = "white", threshold = "50%",channel = "All")
         #pngfile_big<-image_annotate(pngfile_big,paste(cluster_desc),gravity = "north",size = 50,color = "white")
         pngfile_big<-image_annotate(pngfile_big,cluster_desc,gravity = "north",size = 30*40/nchar(cluster_desc),color = "white")
         
         pngfile_big<-image_annotate(pngfile_big,"0%        Relative intensity        100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
         #pngfile_big<-image_annotate(pngfile_big,"0                                               100",location = "+52+165",gravity = "northeast",size = 30,color = "white",degree=270)
         #pngfile_big<-image_annotate(pngfile_big,"                 Relative intensity                  ",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
  
         #pngfile<-image_average(img_com)
         #pngfile<-image_threshold(pngfile, type = "black",  threshold = "50%")
         pngfile<-image_border(pngfile, bg, "30x30")
         #pngfile<-image_modulate(pngfile, brightness = 150)
         
         
         }
        else if ((Component_plot_coloure=="as.cluster")){
        pngfile_big<-image_average(img_com)
        pngfile_big_info=magick::image_info(pngfile_big)
        
        pngfile_big_ratio<-c(diff(range(imdata@elementMetadata@coord@listData[["x"]])),diff(range(imdata@elementMetadata@coord@listData[["y"]])))/c(pngfile_big_info$width-1*150,pngfile_big_info$height-1*150)
        
        
        
        pngfile_big <- image_draw(pngfile_big)
        #rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
        pixel_size_um_bar<-200
        pixel_size<-pixel_size_um_bar*pixel_size_um/(1/max(pngfile_big_ratio))/1000
        if(!is.integer(pixel_size)){
          pixel_size<-ceiling(pixel_size)
          pixel_size_um_bar<-pixel_size*1000/pixel_size_um*(1/max(pngfile_big_ratio))
        }
        
        rect(690-pixel_size_um_bar-60, 680,690-pixel_size_um_bar-60+pixel_size_um_bar,  680+13, border  = 'white',col  = 'white', lwd = 0)
        text(690, 690, paste0(pixel_size ," mm"),col  = 'white', cex = 3, srt = 0)
        #palette(rainbow(11, end = 0.9))
        #symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
        #        bg = 1:11, inches = FALSE, add = TRUE)
        dev.off()
        pngfile<-pngfile_big
        pngfile_big<-image_border(pngfile_big, bg, "30x30")
        pngfile_big<-image_modulate(pngfile_big, brightness = img_brightness + 25 * length(candidateunique), saturation = 100)
        pngfile_big<-image_threshold(pngfile_big,type = "white", threshold = "50%",channel = "All")
        #pngfile_big<-image_annotate(pngfile_big,paste(cluster_desc),gravity = "north",size = 50,color = "white")
        pngfile_big<-image_annotate(pngfile_big,cluster_desc,gravity = "north",size = 30*40/nchar(cluster_desc),color = "white")
        pngfile_big<-image_annotate(pngfile_big,"0%         Relative intensity        100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
        
        #pngfile_big<-image_annotate(pngfile_big,"0%                                               100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
        #pngfile_big<-image_annotate(pngfile_big,"                Relative intensity                 ",location = "+52+160",gravity = "northeast",size = 30,color = "white",degree=270)
        
        #pngfile<-image_average(img_com)
        #pngfile<-image_threshold(pngfile, type = "black",  threshold = "50%")
        pngfile<-image_border(pngfile, bg, "30x30")
        pngfile<-image_modulate(pngfile, brightness = img_brightness + 50)
      }
        
        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "north",size = 30,color = "white")
        if (remove_cluster_from_grid){
          pngfile=image_append(c((pngcompfile_output)))
        }else{
          pngfile=image_append(c(pngfile,(pngcompfile_output)))
        }
        
        bg = paste0("grey29")
        #property_png<-image_attributes(pngfile)
        #width_height<-as.numeric(unlist(stringr::str_split(property_png[property_png$property=="png:IHDR.width,height",2],", ")))
        #pngfile<-c(image_blank(width_height[1], width_height[2], color = bg, pseudo_image = ""),pngfile)
        image_write(pngfile,outputpngsum,flatten =T)
        
        
      
    }
    
    
    
    
  }
  
  if(export_Header_table){
    
     suppressMessages(suppressWarnings(require(gridExtra)))
     suppressMessages(suppressWarnings(require(grid)))
    
    if(is.null(candidate$desc)){
      candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
    }else{
      candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
      cluster_desc<-unique(candidate$desc)[1]
    }
     candidate_unique_table_Score<-candidate %>% group_by(mz) %>% summarise(Score=max(Score))
     candidate_unique_table_Score$Score<-round(candidate_unique_table_Score$Score,digits = 2)
     candidate_unique_table<-merge(candidate_unique_table,candidate_unique_table_Score)
     suppressMessages(suppressWarnings(require(reshape2)))
    
    Header_table<-NULL
    Header_table$mz=candidateunique
    Header_table<-data.frame(Header_table)
    Header_table=base::merge(Header_table,candidate_unique_table,all.x=T,by="mz",all.y=F)
    Header_table=Header_table[order(as.character(Header_table$mz)),]
    Header_table<-as.data.frame(Header_table[,c(componentID_colname,"mz","formula","adduct","Score")])
    
    #Header_table$ID=candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique)]
    #componentnames=unique(Header_table[[componentID_colname]][Header_table$mz==as.numeric(candidateunique[i])])
    t_Header_table<-as.data.frame(t(Header_table))
    t_Header_table<-sapply(t_Header_table,as.character)
    t_Header_table<-as.data.frame(t_Header_table)
    names(t_Header_table)<-as.character(t_Header_table[1,])
    tt3 <- ttheme_minimal()
    header_table_op<-tableGrob(t_Header_table,theme = tt3,cols = NULL,rows=NULL)
    
    find_cell <- function(table, row, col, name="core-fg"){
      l <- table$layout
      which(l$t==(row) & l$l==(col) & l$name==name)
    }
    if (remove_cluster_from_grid){
      t_Header_table<-cbind(t_Header_table)
      
    }else{ 
      t_Header_table<-cbind(c("ClusterID",clusterID,str_sub(cluster_desc,end = regexpr(" ",cluster_desc)),rep("",nrow(t_Header_table)-3)),t_Header_table)
    }    
    header_table_op<-tableGrob(t_Header_table,cols = NULL,rows=NULL)
    

    if (remove_cluster_from_grid){
    for (mzfeatures in 1:(length(candidateunique))){
      for (rows in 3:3){
        ind <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
        #ind2 <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
        header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures], col = levels(mycol)[mzfeatures], lwd=5)
        #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
        #header_table_op$grobs[ind2][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
        
        #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
      }}
    }else{
    for (mzfeatures in 2:(length(candidateunique)+1)){
    for (rows in 3:3){
    ind <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
    #ind2 <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
    header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = levels(mycol)[mzfeatures-1], lwd=5)
    #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
    #header_table_op$grobs[ind2][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
    
    #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
    }
    }
    }
    
    suppressMessages(suppressWarnings(require(gtable)))
    
    header_table_op <- gtable_add_grob(header_table_op,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 2, b = nrow(header_table_op), l = 1, r = ncol(header_table_op))
    header_table_op <- gtable_add_grob(header_table_op,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(header_table_op))
    header_table_op$widths <- rep(max(header_table_op$widths), length(header_table_op$widths))
    
    header_file_png = windows_filename(paste0(clusterID,"_header.png"))
    if (remove_cluster_from_grid){
      png(header_file_png,width = 5*length(candidateunique),height = 5,units = "in",res = 75)
      
    }else{
      png(header_file_png,width = 5*length(candidateunique+1),height = 5,units = "in",res = 75)
      
    }
    
    grid.arrange(header_table_op,nrow=1)
    
    dev.off()
    
   if(F){ p <- plot_ly(
      type = 'table',
      columnwidth = 20,
      header = list(
        values = c(paste0("<b>","Cluster","</b>"),Header_table[,componentID_colname]),
        colspan = I(30),
        align = c('center'),
        line = list(width = 1, color = 'black'),
        fill = list(color = "grey"),
        font = list(family = "Arial", size = c(15,9+6/(width(Header_table[,componentID_colname])/7)), color = "white")
      ),
      cells = list(
        values = cbind(
          #rbind(paste0("<b>","Cluster","</b>"), as.matrix(Header_table$FA)),
          rbind("",as.matrix(Header_table$adduct)), 
          rbind("",as.matrix(round(Header_table$mz,digits = 3)))
        ),
        align = c('left', rep('center', ncol(mtcars))),
        line = list(color = "black", width = 1),
        fill = list(color = "white"),
        font = list(family = "Arial", size = 15, color = c("black"))
        
      ),
      width=120*(nrow(Header_table)+10),
      height=400,
      
      
    ) %>% layout(
      title = paste0("<b>",clusterID,"</b>"),
      autosize = F,
      margin=0,
      font = list(family = "Arial",  color = "black",align = "bottom")
    )
    
    
    
    
    
    orca(p, file = windows_filename(paste0(clusterID,"_header.png")),width=120*(nrow(Header_table))+200,height=540) }
    
  if(file.exists(outputpngsum)){
      
      clusterpng<-image_read(outputpngsum)
      if (image_info(clusterpng)[2]>output_png_width_limit){
        clusterpng<-image_resize(clusterpng,paste0(output_png_width_limit,"x"))
      }
      header_file_pngfile<-header_file_png
      header_file_png<-image_read(header_file_png)
      header_file_png<-image_trim(header_file_png)
      header_file_png<-image_border(header_file_png, "white", "00x70")
      header_file_png<-image_resize(header_file_png,paste0(image_info(clusterpng)[2],"x"))
      clusterpng<-image_append(c(header_file_png,clusterpng),stack = T)
      
      image_write(clusterpng,outputpngsum)
      try(file.remove(header_file_pngfile))
    }
  }
  
  if(export_footer_table){
     suppressMessages(suppressWarnings(require(colorspace)))
     suppressMessages(suppressWarnings(require(stringr)))
     suppressMessages(suppressWarnings(require(ggplot2)))
    
    if(is.null(candidate$desc)){
      candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
      cluster_desc=ClusterID
    }else{
      candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
      cluster_desc<-unique(candidate$desc)[1]
    }
    
      prosequence<-list_of_protein_sequence[clusterID]
      candidate_unique_table=unique(candidate[,c(ClusterID_colname,componentID_colname,"Intensity","mz")])
      component_int<-candidate_unique_table %>% group_by_at((componentID_colname)) %>% dplyr::summarise(int=sum(Intensity))
      component_int$int<-component_int$int/max(component_int$int)
      
      s1=as.character(component_int$Peptide)
      s2=as.character(prosequence)
      
      palign2 <- sapply(s1,regexpr , s2)
      width_com<-str_length(s1)
      component_int$start=palign2
      component_int$end=component_int$start+width_com-1
      component_int<-merge(component_int,unique(candidate_unique_table[,c(componentID_colname,"mz")]),by=componentID_colname)
      component_int<-component_int[order(as.character(component_int$mz)),]
      
      component_int$mycol<-levels(mycol)
      
      
      
      transcolor<-rgb(0, 0, 0, max = 255, alpha = 0)
      pro_length<-unname(width(s2))
      pro_int<-rep(0,pro_length)
      pro_col<-rep(transcolor,pro_length)
      pro_col<-rep("grey93",pro_length)
      for(y in 1:nrow(component_int)){
        for( t in component_int$start[y]:component_int$end[y]){
          #pro_col[t]<-mixcolor(component_int$int[y]/(pro_int[t]+component_int$int[y]), col2RGB(pro_col[t]), col2RGB(component_int$mycol[y]))
          mixedcolor<-colorRamp(colors=c(pro_col[t],component_int$mycol[y]),  space ="rgb",
                    interpolate = "linear")((pro_int[t]+component_int$int[y])/(pro_int[t]+component_int$int[y]))
          pro_col[t]<-rgb(mixedcolor[,1],mixedcolor[,2],mixedcolor[,3], maxColorValue = 255)
          pro_int[t]<-pro_int[t]+component_int$int[y]
        }
      }
      
      ncharrow<-ceiling(width(s2)/length(candidateunique)/5)
      ncharw<-floor(width(s2)/ncharrow)
      component_int_plot<-data.frame(site=1:pro_length,int=pro_int,col=pro_col)
      component_int_plot$x<-(component_int_plot$site-1) %% ncharw
      component_int_plot$y<--((component_int_plot$site-1) %/% ncharw)
      component_int_plot$char<-str_sub(s2,start = component_int_plot$site,end = component_int_plot$site)
      
      
      if(F){p<-ggplot(component_int_plot, aes(x=site, y=int,fill=col)) + geom_area() + theme(legend.position="none") + ggtitle("Plot of length \n by dose") +
        theme(
          plot.title = element_text(color="red", size=14, face="bold.italic"),
          axis.title.x = element_text(color="blue", size=14, face="bold"),
          axis.title.y = element_text(color="#993333", size=14, face="bold")
        )}
      
      footerpng<-paste(workdir,"/",windows_filename(substr(clusterID, 1, 10)),"_footer.png",sep="")
      
      if (footer_style=="Protein"){
      png(footerpng,width = 5*length(candidateunique+1),height = 5*ceiling(ncharrow/10),units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),mar=c(1, 0, 0, 0))
      p <- ggplot(component_int_plot, aes(x, y, label = char)) + 
        geom_label(fill=component_int_plot$col,family = "mono",size=20) + 
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
              xlab(cluster_desc) + ylim(min(component_int_plot$y)-0.4,max(component_int_plot$y)+0.4)
      
      print(p)
      dev.off()
      } else if (footer_style=="Length"){
              
      png(footerpng,width = 5*length(candidateunique+1),height = 1,units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),mar=c(1, 0, 0, 0))
      component_int_plot$x=as.factor(1)
      component_int_plot$y=1
      p<- ggplot(data=component_int_plot, aes(x=site, y=1,group=site,fill=component_int_plot$col,col=component_int_plot$col)) +
        geom_bar(stat="identity",fill=component_int_plot$col,col=component_int_plot$col)+ 
        theme(axis.line=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
         
        ylab(cluster_desc)
      print(p)
      dev.off()
      
      }


      
      if(file.exists(outputpngsum)){
        
        clusterpng<-image_read(outputpngsum)
        footerpngfile<-footerpng
        footerpng<-image_read(footerpng)
        footerpng<-image_border(footerpng, "white", "00x70")
        footerpng<-image_resize(footerpng,paste0(image_info(clusterpng)[2],"x"))
        clusterpng<-image_append(c(clusterpng,footerpng),stack = T)
        
        image_write(clusterpng,outputpngsum)
        try(file.remove(footerpngfile))
      }
      
      if(F){
              wrap_strings <- function(vector_of_strings,width){as.character(sapply(vector_of_strings,FUN=function(x){paste(strwrap(x,width=width), collapse="\n")}))}
      
      plot(component_int_plot$site,component_int_plot$int,col=component_int_plot$col)
      title(bquote(wrap_strings(s2,50)),col.main="black",cex.main=0.25,adj  = 0,line = -1) 
      
      
      for (components in 1:nrow(component_int)){
        pep_start<-component_int$start[components]
        pep_end<-component_int$end[components]
        pep_body<-str_sub(s2,pep_start,pep_end)
        if(pep_start==1){
          title(bquote(.(pep_body)*phantom(.(str_sub(s2,pep_end+1,pro_length)))),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }else if(pep_end==pro_length){
          title(bquote(phantom(.(str_sub(s2,1,pep_start-1)))*.(pep_body)),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }else{
          title(bquote(phantom(.(str_sub(s2,1,pep_start-1)))*.(pep_body)*phantom(.(str_sub(s2,pep_end+1,pro_length)))),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }
      }
      
      dev.off()
      }

    }
  
  if(attach_summary_cluster){
    pngfile<-image_read(outputpngsum)
    #temp_component_png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
    
    #png(temp_component_png,width = 5,height = 5 , bg = bg,units = "in",res = 300)
    #print(pngfile_big)
    #dev.off()
    #pngfile_foot<-image_read(temp_component_png)image_scale(img, "100")
    
    #bg = paste0("grey15")
    #property_png<-image_attributes(pngfile_big)
    #width_height<-as.numeric(unlist(stringr::str_split(property_png[property_png$property=="png:IHDR.width,height",2],", ")))
    #pngfile_big<-image_flatten(c(image_blank(width_height[1], width_height[2], color = bg, pseudo_image = ""),pngfile_big))
    
    #image_write(pngfile,outputpngsum,flatten =T)
    
    pngfiletry=tryCatch(image_append(c(pngfile,pngfile_big),stack = T))
   # 
    if (class(pngfiletry)=="magick-image") image_write(pngfiletry,outputpngsum)
    #///rm(temp_component_png)
    
  }
    
    removetempfile(tmp_dir,matchword=c(".png$","^magick-"))
    
    rm(pngcompfile)
    
    message("Cluster image rendering done:", clusterID ,cluster_desc)
    
  }else{
    
  message("Cluster image rendering Skipped:", clusterID ,cluster_desc)
    
  }
  
  
  }

orca_initial<-function(){
  

Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/Anaconda/orca_app", sep = .Platform$path.sep))

  
}

affine <- function(x, translate=c(0,0), rotate=0,
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
  
  rotatenew<-affine(rotatetmp[,c("x","y")])
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

Ctypeof <- function(type) {
  vapply(type, function(t) {
    switch(t,
           `16-bit integer` = "short",
           `32-bit integer` = "int",
           `64-bit integer` = "long",
           `32-bit float` = "float",
           `64-bit float` = "double",
           stop("unrecognized binary type"))
  }, character(1))
}

seq.ppm <- function(from, to, ppm) {
  length.out <- (log(to) - log(from)) / log((1 + 1e-6 * ppm) / (1 - 1e-6 *ppm))
  length.out <- floor(1 + length.out)
  i <- seq_len(length.out)
  from * ((1 + 1e-6 * ppm) / (1 - 1e-6 * ppm))^(i-1)
}

intensity.colors_customize1 <- function(n = 100, alpha = 1,colset=1) {
  col2 <- rainbow(3*n, alpha=alpha)[(2*n):1]
  if (colset==1){
  f <- colorRamp(colors=c("darkorchid4", "dodgerblue","green", "greenyellow","yellow"))  
  } else if(colset==2){
  f <- colorRamp(colors=c("black", "dodgerblue","green", "greenyellow","yellow"))
  }
  else {
  f <- colorRamp(colors=c("darkorchid4", "blue", "darkseagreen1", "yellow"))  
  }
  
  alpha <- col2rgb(col2, alpha=TRUE)[[4]]
  col1 <- sapply(seq(from=0, to=1, length.out=n), function(i) do.call(rgb,
                                                                      c(as.list(f(i)), maxColorValue=255, alpha=alpha)))
  col1
}

Scilslab_mass_interval<-function(masslist,ppm){
  
  writelines=list()
  writelines[[1]]="<ImagingResults flexImagingVersion=\"3.0\" last_modified=\"\" created=\"SCiLS export\">"
  line=2
  for (mz in rownames(masslist)){
    
    writelines[[line]]=paste0("<Result Type=\"PkFilter\" Name=\"",masslist[mz,"moleculeNames"],
                              "\" Color=\"#5555ff\" Show=\"0\" MinIntensity=\"0\" IntensityThreshold=\"100\" AbsIntens=\"0\" LogScale=\"0\" MinMass=\"",
                              masslist[mz,"mz"]-(masslist[mz,"mz"]*ppm/1000000),
                              "\" MaxMass=\"",
                              masslist[mz,"mz"]+(masslist[mz,"mz"]*ppm/1000000),
                              "\" Integrate=\"0\" FindMass=\"0\" RelMass=\"1\"></Result>")
    line=line+1
  }
  writelines[[line]]="</ImagingResults>"
    
  writeLines(unlist(writelines),"massinterval.MIR")  
    
    
    
  
  
}

removetempfile<-function(temp_dir=tempdir(),matchword=c(".png$","^magick-"), intersect=F){
  filelist=list.files(temp_dir)
  if (intersect){
    test=grep(paste0(matchword,collapse = "+"),filelist,ignore.case = T)
  }else{
    test=grep(paste0(matchword,collapse = "|"),filelist,ignore.case = T)
  }
  file.remove(paste0(temp_dir,"/",filelist[test]))
  
  
}

testcolume<-function(df,testcolnames){
  
   suppressMessages(suppressWarnings(require(stringr)))
  
  testcolnames=unique(testcolnames)
  
  test=sapply(testcolnames,FUN = grepl,names(df))
  
  testresult=sapply(colnames(test),FUN = function(x,df){sum(df[,x])},test) == 0
  
  case_test=sapply(testcolnames,FUN = grepl,names(df),ignore.case = T)
  
  case_testresult=sapply(testcolnames,FUN = function(x,df){sum(df[,x])},case_test) > 1
  
  duplicatecol=names(case_testresult)[case_testresult==T]
  
  if (length(duplicatecol)==0){duplicatecol=NULL}
  
  testfail=names(testresult)[testresult==T]
  
  passcol=names(testresult)[testresult==F]
  
  renamecol=NULL
  
  if (length(passcol)<length(testcolnames)){
    
  case_test=sapply(testfail,FUN = grepl,names(df),ignore.case = T)
  
  case_testresult=sapply(colnames(case_test),FUN = function(x,df){sum(df[,x])},case_test) == 1
  
  renamecol=names(case_testresult)[case_testresult==T]
  
  if (length(renamecol)>0){
    
  failcol=testfail[!grep(renamecol,testfail)]  
    
  }else{
    
    failcol=NULL
    
  }
  
  
    
  } else{
    
  failcol=NULL
    
  }

  
  
  return(list(passcol=passcol,renamecol=renamecol,failcol=failcol,duplicatecol=duplicatecol))
}

data_test_rename<-function(required_col,df){
  
testcolumeresult = testcolume(df,required_col)

if (length(testcolumeresult$failcol)>0){
  
 stop(paste(testcolumeresult$failcol,"column is missing in the input file, please check your datafile"))  
  
}

if (length(testcolumeresult$renamecol)>0){
  
lapply(testcolumeresult$renamecol,gsub,testcolumeresult$renamecol,ignore.case = T,x=names(df))
  
for (i in length(testcolumeresult$renamecol)){
  
  names(df)=gsub(testcolumeresult$renamecol,testcolumeresult$renamecol,ignore.case = T,x=names(df))
  
}

}

if (length(testcolumeresult$duplicatecol)>0){
  
  message(paste("Found ambigous or duplicate column: ",testcolumeresult$duplicatecol))
  
}

return(df)

}

Build_adduct_list<-function(){
  
  Name=as.character(c("M+H","M+NH4","M+Na","M+K","M+","M-H","M-2H","M-3H","M+FA-H","M+Hac-H",
                      "M-","M+3H","M+2H+Na","M+H+2Na","M+3Na","M+2H","M+H+NH4","M+H+Na","M+H+K",
                      "M+ACN+2H","M+2Na","M+2ACN+2H","M+3ACN+2H","M+CH3OH+H","M+ACN+H","M+2Na-H",
                      "M+IsoProp+H","M+ACN+Na","M+2K-H","M+DMSO+H","M+2ACN+H","M+IsoProp+Na+H","2M+H",
                      "2M+NH4","2M+Na","2M+3H2O+2H","2M+K","2M+ACN+H","2M+ACN+Na","M-H2O-H","M+Na-2H",
                      "M+Cl","M+K-2H","M+Br","M+TFA-H","2M-H","2M+FA-H","2M+Hac-H","3M-H","M+He",
                      "M+Ne","M+Ar","M+Kr","M+Xe","M+Rn","M+Cu","M+Co","M+Ag"
  ))
  calc=c("M+1.007276","M+18.033823","M+22.989218","M+38.963158","M-0.00054858","M-1.007276","M/2-1.007276",
         "M/3-1.007276","M+44.998201","M+59.013851","M+0.00054858","M/3+1.007276","M/3+8.334590","M/3+15.7661904",
         "M/3+22.989218","M/2+1.007276","M/2+9.520550","M/2+11.998247","M/2+19.985217","M/2+21.520550","M/2+22.989218",
         "M/2+42.033823","M/2+62.547097","M+33.033489","M+42.033823","M+44.971160","M+61.06534","M+64.015765",
         "M+76.919040","M+79.02122","M+83.060370","M+84.05511","2M+1.007276","2M+18.033823","2M+22.989218",
         "M+28.02312","2M+38.963158","2M+42.033823","2M+64.015765","M-19.01839","M+20.974666","M+34.969402",
         "M+36.948606","M+78.918885","M+112.985586","2M-1.007276","2M+44.998201","2M+59.013851","3M-1.007276",
         "M+4.002606","M+19.992439","M+39.962383","M+83.911507","M+131.90416","M+222.017563","M+62.9285022",
         "M+58.9321032","M+106.9034432"
  )
  Charge=c(" 1"," 1"," 1"," 1"," 1","-1","-2","-3","-1","-1","-1"," 3"," 3"," 3"," 3"," 2"," 2"," 2",
           " 2"," 2"," 2"," 2"," 2"," 1"," 1"," 1"," 1"," 1"," 1"," 1"," 1"," 1"," 1"," 1"," 1"," 1"," 1",
           " 1"," 1","-1","-1","-1","-1","-1","-1","-1","-1","-1","-1","0","0","0","0","0","0","2","2","2"
  )
  Mult=c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1",
         "1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","1","1","1","1","1","1","2","2","2","3",
         "1","1","1","1","1","1","1","1","1")
  Mass=as.numeric(as.character(c("  1.00727600"," 18.03382300"," 22.98921800"," 38.96315800"," -0.00054858",
                                 " -1.00727600"," -1.00727600"," -1.00727600"," 44.99820100"," 59.01385100",
                                 "  0.00054858","  1.00727600","  8.33459000"," 15.76619000"," 22.98921800",
                                 "  1.00727600","  9.52055000"," 11.99824700"," 19.98521700"," 21.52055000",
                                 " 22.98921800"," 42.03382300"," 62.54709700"," 33.03348900"," 42.03382300",
                                 " 44.97116000"," 61.06534000"," 64.01576500"," 76.91904000"," 79.02122000",
                                 " 83.06037000"," 84.05511000","  1.00727600"," 18.03382300"," 22.98921800",
                                 " 28.02312000"," 38.96315800"," 42.03382300"," 64.01576500","-19.01839000",
                                 " 20.97466600"," 34.96940200"," 36.94860600"," 78.91888500","112.98558600",
                                 " -1.00727600"," 44.99820100"," 59.01385100","  1.00727600","4.002606","19.992439",
                                 "39.962383", "83.911507","131.90416","222.017563","62.9285022","58.9321032",
                                 "106.9034432"
                                 )))
  Ion_mode=c("positive","positive","positive","positive","positive","negative","negative","negative"
             ,"negative","negative","negative","positive","positive","positive","positive","positive"
             ,"positive","positive","positive","positive","positive","positive","positive","positive"
             ,"positive","positive","positive","positive","positive","positive","positive","positive"
             ,"positive","positive","positive","positive","positive","positive","positive","negative"
             ,"negative","negative","negative","negative","negative","negative","negative","negative"
             ,"negative","positive","positive","positive","positive","positive","positive","positive"
             ,"positive","positive")
  formula_add=c("H1","N1H4","Na1","K1","FALSE","FALSE","FALSE","FALSE","C1O2H2","C2O2H4","FALSE","H3",
                "H2Na1","H1Na2","Na3","H2","H1N1H4","H1Na1","H1K1","C2H5N1","Na2","C4H8N2","C6H11N3",
                "C1H5O1","C2H4N1","Na2","C3H9O1","C2H3N1Na1","K2","C2H7S1O1","C4H7N2","C3H9O1Na1","H1",
                "N1H4","Na1","H8O6","K1","C2H4N1","C2H3N1Na1","FALSE","Na1","Cl1","K1","Br1","C2F3O2H1",
                "FALSE","C1O2H2","C2O2H4","FALSE","He","Ne","Ar","Kr","Xe","Rn","Cu","Co","Ag"
  )
  formula_ded=c("FALSE","FALSE","FALSE","FALSE","FALSE","H1","H2","H3","H1","H1","FALSE","FALSE","FALSE",
                "FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE",
                "FALSE","H1","FALSE","FALSE","H1","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE",
                "FALSE","FALSE","FALSE","H3O1","H2","FALSE","H2","FALSE","H1","H1","H1","H1","H1","FALSE",
                "FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE")
  Multi=c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"
          ,"1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","1","1","1","1","1","1","2","2","2",
          "3","1","1","1","1","1","1","1","1","1")
  adductslist<-cbind.data.frame(Name,calc,Charge,Mult,Mass,Ion_mode,formula_add,formula_ded,Multi)
  adductslist
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


plotRanges <- function(ranged,labels=NULL,do.labs=T,skip.plot.new=F,lty="solid", alt.y=NULL,
                       v.lines=FALSE,ylim=NULL,xlim=NULL,scl=c("b","Kb","Mb","Gb"),
                       col=NULL,srt=0,pos=4,pch=1,lwd=1,cex=1,...) {
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) { 
    warning("ranged needs to be a RangedData or GRanges object, plot likely to fail") ; return(NULL) }
  chk <- chrNums(ranged)
  typ <- is(ranged)[1]
  if(!is.null(alt.y)) {
    if(is.numeric(alt.y)) {
      if(length(alt.y)==1 | length(alt.y)==length(ranged)) {
        yy <- alt.y
      } else {
        warning("alt.y ignored, must be same length as ranged, or else length 1"); alt.y <- NULL
      }
    } else {
      if(is.character(alt.y)) {
        if(typ=="GRanges") { 
          cn <- colnames(mcols(ranged)); df <- mcols(ranged)
        } else {
          cn <- colnames(ranged); df <- ranged
        }
        if(!alt.y %in% cn) { stop("alternative y.axis column name ",alt.y," not found in 'ranged'") }
        yy <- df[,alt.y]; rm(df)
      } else { 
        warning("invalid value for alt.y, ignoring"); alt.y <- NULL
      }
    }
  }
  if(!is.null(labels)) {
    labels <- paste(labels)
    if(is.character(labels)) {
      if(length(labels)==1 | length(labels)==length(ranged)) {
        if(length(labels)==1) {
          if(typ=="GRanges") { 
            cn <- colnames(mcols(ranged)); df <- mcols(ranged)
          } else {
            cn <- colnames(ranged); df <- ranged
          }
          if(!labels %in% cn) { stop("labels column name ",labels," not found in 'ranged'") }
          lab <- df[,labels]; rm(df)
        } else {
          lab <- labels
        }
      } else {
        warning("labels ignored, must be same length as ranged, or else length 1"); labels <- NULL
      }
    } else {
      warning("invalid value for labels, ignoring"); labels <- NULL
    }
  } else {
    lab <- rownames(ranged) 
  } 
  if(length(chk)>1) { 
    warning(length(chk)," chromosomes in 'ranged', only using the first, chr",chk[1]) 
    ranged <- chrSel(ranged,1) 
  }
  if(all(width(ranged)<=1)) { theyAreSnps <- TRUE } else { theyAreSnps <- FALSE }
  scl <- make.divisor(scl)
  xl <- range(c(start(ranged),end(ranged)),na.rm=T)
  xl <- xl + ((diff(xl)*0.1)*c(-1,1))
  xl <- xl/scl
  nr <- nrow(ranged); if(is.null(nr)) { nr <- length(ranged) }
  if(is.null(alt.y)) {
    yl <- c(0,(nr+2))
  } else {
    yl <- range(yy,na.rm=T)
  }
  if(is.numeric(ylim) & length(ylim)==2) {
    ylim <- range(ylim,na.rm=T)
    ydif <- diff(ylim)
    yl <- ylim
  }
  if(is.numeric(xlim) & length(xlim)==2) {
    xlim <- range(xlim,na.rm=T)
    xdif <- diff(xlim)
    xl <- xlim
  }
  if(is.null(alt.y)) {
    YY <- seq(from=yl[1],to=yl[2],length.out=nr+2)[-1]
  } else {
    if(length(yy)==1) { YY <- rep(yy,length(nr)) } else { YY <- yy }
  }
  #print(YY)
  if(!is.null(col)) {
    if(length(col)==1) {
      col <- rep(col,times=nr) 
    } else {
      if(length(col)!=nr) { warning("col was not the same length as ranged, using first only"); col <- rep(col[1],nr) }
    }
  }
  if(is.null(col)) {
    if(nr>22) { colz <- rep("black",nr) } else { colz <- get.distinct.cols(nr) }
  } else { colz <- col[1:nr] }
  if(is.null(lab) & do.labs) { lab <- paste(1:nr) } # last resort
  if(!skip.plot.new) {
    position <- c(start(ranged[1,]),end(ranged[1,]))/scl
    Y <- YY[c(1,1)]
    #prv(position,Y)
    TY <- if(theyAreSnps) { "p" } else { "l" }
    if(v.lines) {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col="white", lty=lty, ...)
      abline(v=position,col=colz[1])
    } else {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col=colz[1], lty=lty, lwd=lwd, cex=cex, ...)
    }
    st <- 2
  } else {
    st <- 1
  }
  if(nr>1 | st==1) {
    for (cc in st:nr) {
      if(v.lines) {
        abline(v=c(start(ranged[cc,]),end(ranged[cc,]))/scl,col=colz[cc],lty=lty)
      } else {
        if(theyAreSnps) { 
          points(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], pch=pch, cex=cex)
        } else { 
          lines(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], lty=lty, lwd=lwd)
        }
      }
    }
  }
  if(do.labs) {
    for (cc in 1:nr) {
      if(v.lines) { YY <- rep(tail(YY,1),length(YY)) }
      V.scale <- (diff(head(YY,2))*0.5)
      if(length(V.scale)<1 | srt!=90) { V.scale <- 0 }
      text(x=start(ranged[cc,])/scl,y=YY[cc]+V.scale,labels=lab[cc],cex=0.6,pos=pos,offset=0,srt=srt)
    }
  }
}

recheck_peptide_score<-function(formula="AGLQFPVGR",peaklist=read.csv(paste0(getwd(),"/Spectrum 2 .csv"))){
 peaklist 
}

IAnnotatedDataFrame_H<-function (data, varMetadata, dimLabels = c("pixelNames", 
                                           "pixelColumns"), ...) 
{
  reqLabelTypes <- c("dim", "sample", "pheno")
  if (missing(data)) 
    data <- data.frame(sample = factor())
  if (missing(varMetadata)) 
    varMetadata <- data.frame(labelType = factor(rep(NA, ncol(data)), levels = reqLabelTypes), row.names = names(data))
  Cardinal:::IAnnotatedDataFrame(data = data, varMetadata = varMetadata, 
                       dimLabels = dimLabels)
}

col2RGB<-function(x){
   suppressMessages(suppressWarnings(require(colorspace)))
  x_RGB<-t(as.numeric(col2rgb(x)))
  return(RGB(x_RGB))
}
rgb2hex <- function(r,g,b) sprintf('#%s',paste(as.hexmode(c(r,g,b)),collapse = ''))

get_protein_from_interest_desc<-function(Protein_desc_of_interest){
  Index_of_protein_sequence<-get("Index_of_protein_sequence", envir = .GlobalEnv)
  if (Protein_desc_of_interest!="."){
    Protein_feature_list_interest<-NULL
    for (interest_desc in Protein_desc_of_interest){
      Protein_feature_list_interest<-rbind(Protein_feature_list_interest,Index_of_protein_sequence[grepl(paste0(" ",interest_desc),Index_of_protein_sequence$desc,ignore.case = T),])
      Protein_feature_list_interest<-rbind(Protein_feature_list_interest,Index_of_protein_sequence[grepl(paste0("-",interest_desc),Index_of_protein_sequence$desc,ignore.case = T),])
    }
    #Protein_feature_list_crystallin$Protein=as.character(Protein_feature_list_crystallin$desc)
    return(unique(Protein_feature_list_interest$recno))
  }else{
    return(unique(Index_of_protein_sequence$recno))
    }
  
}
percent <- function(x, digits = 2, format = "f", ...) {
  paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
}

#' Cleavage_rules_fun
#'
#' This is a function will return a list of pre-defined enzyme digestion specificity. 
#' 
#' @return a table of Cleavage rules 
#'
#' @examples
#' Cleavage_rules_fun()
#'
#' @export

Cleavage_rules_fun <- function(){
  c(
  ## Arg-C proteinase
  "arg-c proteinase"="R(?=\\w)",
  ## Asp-N endopeptidase
  "asp-n endopeptidase"="\\w(?=D)",
  ## BNPS-Skatole
  "bnps-skatole"="W(?=\\w)",
  ## Caspase 1
  "caspase1"="(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])",
  ## Caspase 2
  "caspase2"="(?<=DVA)D(?=[^PEDQKR])",
  ## Caspase 3
  "caspase3"="(?<=DMQ)D(?=[^PEDQKR])",
  ## Caspase 4
  "caspase4"="(?<=LEV)D(?=[^PEDQKR])",
  ## Caspase 5
  "caspase5"="(?<=[LW]EH)D(?=\\w)",
  ## Caspase 6
  "caspase6"="(?<=VE[HI])D(?=[^PEDQKR])",
  ## Caspase 7
  "caspase7"="(?<=DEV)D(?=[^PEDQKR])",
  ## Caspase 8
  "caspase8"="(?<=[IL]ET)D(?=[^PEDQKR])",
  ## Caspase 9
  "caspase9"="(?<=LEH)D(?=\\w)",
  ## Caspase 10
  "caspase10"="(?<=IEA)D(?=\\w)",
  ## Chymotrypsin - high specifity
  "chymotrypsin-high"="([FY](?=[^P]))|(W(?=[^MP]))",
  ## Chymotrypsin - low specifity
  "chymotrypsin-low"="([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))",
  ## Clostripain
  "clostripain"="R(?=\\w)",
  ## CNBr
  "cnbr"="M(?=\\w)",
  ## Enterokinase
  "enterokinase"="(?<=[DE][DE][DE])K(?=\\w)",
  ## Factor Xa
  "factor xa"="(?<=[AFGILTVM][DE]G)R(?=\\w)",
  ## Formic acid
  "formic acid"="D(?=\\w)",
  ## Glutamyl endopeptidase
  "glutamyl endopeptidase"="E(?=\\w)",
  ## Granzyme B
  "granzyme-b"="(?<=IEP)D(?=\\w)",
  ## Hydroxylamine
  "hydroxylamine"="N(?=G)",
  ## Iodosobenzoic acid
  "iodosobenzoic acid"="W(?=\\w)",
  ## LysC
  "lysc"="K(?=\\w)",
  ## LysN
  "lysn"="\\w(?=K)",
  ## Neutrophil elastase
  "neutrophil elastase"="[AV](?=\\w)",
  ## NTCB (2-nitro-5-thiocyanobenzoic acid)
  "ntcb"="\\w(?=C)",
  ## Pepsin (pH 1.3)
  "pepsin1.3"="((?<=([^HKR][^P])|(^[^P]))[^R](?=[FL][^P]))|((?<=([^HKR][^P])|(^[^P]))[FL](?=\\w[^P]))",
  ## Pepsin (pH > 2.0)
  "pepsin"="((?<=([^HKR][^P])|(^[^P]))[^R](?=[FLWY][^P]))|((?<=([^HKR][^P])|(^[^P]))[FLWY](?=\\w[^P]))",
  ## Proline endopeptidase
  "proline endopeptidase"="(?<=[HKR])P(?=[^P])",
  ## Proteinase K
  "proteinase k"="[AEFILTVWY](?=\\w)",
  ## Staphylococcal Peptidase I
  "staphylococcal peptidase i"="(?<=[^E])E(?=\\w)",
  ## Thermolysin
  "thermolysin"="[^DE](?=[AFILMV])",
  ## Thrombin
  "thrombin"="((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVW]P)R(?=[^DE][^DE]))",
  ## Trypsin
  "trypsin"="([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))")
}

#' grid.ftable
#'
#' This is a function will plot table. 
#' 
#' @return a table image 
#'
#' @examples
#' grid.ftable()
#'
#' @export
grid.ftable <- function(d, padding = unit(4, "mm"), ...) {
  library(gridExtra)
  library(grid)
  nc <- ncol(d)
  nr <- nrow(d)
  
  ## character table with added row and column names
  extended_matrix <- cbind(c("", rownames(d)),
                           rbind(colnames(d),
                                 as.matrix(d)))
  
  ## string width and height
  w <- apply(extended_matrix, 2, strwidth, "inch")
  h <- apply(extended_matrix, 2, strheight, "inch")
  
  widths <- apply(w, 2, max)
  heights <- apply(h, 1, max)
  
  padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)
  
  x <- cumsum(widths + padding) - 0.5 * padding
  y <- cumsum(heights + padding) - padding
  
  rg <- rectGrob(x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 width = unit(widths + padding, "in"),
                 height = unit(heights + padding, "in"))
  
  tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 just = "center")
  
  g <- gTree(children = gList(rg, tg), ...,
             x = x, y = y, widths = widths, heights = heights)
  
  grid.draw(g)
  invisible(g)
}

brewer.pal_n<-function(n,colorstyle){
  library(RColorBrewer)
  if (n>9){
    getPalette = colorRampPalette(brewer.pal(8, colorstyle))
    return(getPalette(n))
  }else{
    getPalette = brewer.pal(n, colorstyle)
    return(getPalette)
    
  }
  
  
}

Parallel.OS<-function(Thread=1,bpprogressbar_t=TRUE,override_type=NULL){
  library(BiocParallel)
  if(Thread!=1){
    if(!is.null(override_type)){
      BPPARAM=override_type
    }else{
          switch(Sys.info()[['sysname']],
         Windows= {BPPARAM=BiocParallel::SnowParam()},
         Linux  = {BPPARAM=BiocParallel::MulticoreParam()},
         Darwin = {BPPARAM=BiocParallel::MulticoreParam()}) 
    }

    BiocParallel::bpworkers(BPPARAM)=Thread
  }else{
    BPPARAM=BiocParallel::SerialParam() 
    } 
  
  bpprogressbar(BPPARAM)=bpprogressbar_t
  return(BPPARAM)
}


rotateMSI<-function(imdata,rotation_degree=0){
  
  rotatetmp<-as.data.frame(imdata@elementMetadata@coord@listData)
  
  rotatenew<-affine(rotatetmp[,c("x","y")],rotate=rotation_degree,grid = T)
  
  #sum(duplicated(rotatenew))
  
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
  png(paste(outputdir,"\\","PCA_image.png",sep=""),width = 1024,height = 720)
    
  print(image(PCA_imdata, values="scores", superpose=FALSE, layout=c(3,4),normalize.image = c("linear"),contrast.enhance = c("histogram")))
  
  dev.off()
  png(paste(outputdir,"\\","PCA_plot.png",sep=""),width = 1024,height = 720 ,res=150) 
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
                                     ppm=3,import_ppm=3,Bypass_Segmentation=F,
                                     mzrange="auto-detect",
                                     Segmentation=c("spatialKMeans","spatialShrunkenCentroids","Virtual_segmentation","none","def_file"),
                                     Segmentation_def="segmentation_def.csv",
                                     Segmentation_ncomp="auto-detect",
                                     Segmentation_variance_coverage=0.8,
                                     Smooth_range=1,
                                     colorstyle="Set1",
                                     Virtual_segmentation_rankfile=NULL,
                                     rotate=NULL,
                                     BPPARAM=bpparam(),
                                     preprocess=list(force_preprocess=FALSE,use_preprocessRDS=TRUE,smoothSignal=list(method="gaussian"),
                                                     reduceBaseline=list(method="locmin"),
                                                     peakPick=list(method="adaptive"),
                                                     peakAlign=list(tolerance=ppm, units="ppm"),
                                                     normalize=list(method=c("rms","tic","reference")[1],mz=1)),
                                     ...){
  
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(Cardinal)))
  suppressMessages(suppressWarnings(require(RColorBrewer)))
  suppressMessages(suppressWarnings(require(stringr)))
  #getPalette = colorRampPalette(brewer.pal_n(9, colorstyle))
  setCardinalBPPARAM(BPPARAM)
  
  # if (is.null(workdir)){
  #   workdir<-base::dirname(datafile[1])
  #   datafile<-basename(datafile)
  # }else{ workdir<-workdir }
  datafile <- paste0(workdir,"/",datafile)
  workdir <- dirname(datafile)
  datafile <- basename(datafile)
  if (!is.null(rotate)){
    message("Found rotation info")
    if (typeof(rotate)=="character"){rotate=read.csv(paste0(workdir[1],"/",rotate),stringsAsFactors = F)}
    
    #rotatedegrees=rotate[rotate$filenames==datafile,"rotation"]
    rotatedegrees=sapply(datafile,function(x,df){
      library(stringr)
      if (!str_detect(x,".imzML$")){
        x<-paste0(x,".imzML")
      }
      degree=df[df$filenames==(x),"rotation"]
      if (length(degree)==0) {
        message("Missing rotation data please check the rotation configuration file: ",x)
        degree=0
      }
      degree
    },rotate)
    rotate=unlist(rotatedegrees)
  }else{rotate=rep(0,length(datafile))}

  

  datafile_imzML<-datafile
  for (z in 1:length(datafile)){
    Segmentation_ncomp_running<-Segmentation_ncomp
    name <-basename(datafile[z])
    name <-gsub(".imzML$","",name)
    name <-gsub("/$","",name)
    folder<-base::dirname(datafile[z])
    #imdata <- Cardinal::readImzML(datafile[z],preprocessing = F,attach.only = T,resolution = 200,rotate = rotate[z],as="MSImageSet",BPPARAM = BPPARAM)
    if (!str_detect(datafile[z],".imzML$")){
      datafile_imzML[z]<-paste0(datafile[z],".imzML")
    }
    setwd(workdir[z])
    #message("Porject dir",workdir)
    # imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=200, units="ppm",BPPARAM=BPPARAM,mass.range =mzrange)
    # if(!is.na(rotate[datafile_imzML[z]])){
    #   imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
    # }else if(!is.na(rotate[datafile[z]])){
    #   imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
    # }
    # 
    if (ppm>=25) {
      instrument_ppm=50
    }else{
      instrument_ppm=8
    }      
    
    
    if (dir.exists(paste0(gsub(".imzML$","",datafile[z]) ," ID"))==FALSE){
      dir.create(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
    }
    
    message("Preparing image data for statistical analysis: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
    if(mzrange=="auto-detect"){
      imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam())
    }else {
      imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam(),mass.range =mzrange)
    }
    if(!is.na(rotate[datafile_imzML[z]])){
      imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
    }else if(!is.na(rotate[datafile[z]])){
      imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
    }
    
    
    
    if ('|'(!file.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS")),
            preprocess$force_preprocess)){
      message("Preparing image data for statistical analysis: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
      if (!is.null(preprocess)){
        if  ( ppm<25){
          imdata_ed<-imdata 
            #smoothSignal(method="gaussian") %>% 
            #reduceBaseline(method="locmin") %>%
            if (preprocess$smoothSignal$method=="Disable") {
            }else if (!is.null(preprocess$smoothSignal$method)){
              imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
            }else{
              imdata_ed<- imdata_ed %>% smoothSignal(method="gaussian")
            }
          
            if (!is.null(preprocess$reduceBaseline$method)){
              imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
            }else{
              
            }
            
            if (preprocess$peakPick$method=="Disable") {
            }else if (!is.null(preprocess$peakPick$method)){
              imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method)
            }else{
              imdata_ed<- imdata_ed %>% peakPick(method="adaptive")
            }
          
            if (preprocess$peakAlign$tolerance==0) {
            }else if ('&'(!is.null(preprocess$peakAlign$tolerance),!is.null(preprocess$peakAlign$tolerance))){
              imdata_ed<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units)
            }else {
              imdata_ed<- imdata_ed %>% peakAlign(tolerance=ppm, units="ppm")
            }
          
            imdata_ed<- imdata_ed %>% process()
          
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
            
        } else if(ppm>=25){
          imdata_ed<-imdata 
          #smoothSignal(method="gaussian") %>% 
          #reduceBaseline(method="locmin") %>%
          if (!is.null(preprocess$smoothSignal$method)){
            imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
          }else{
            imdata_ed<- imdata_ed %>% smoothSignal(method="gaussian")
          }
          
          if (!is.null(preprocess$reduceBaseline$method)){
            imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
          }else{
            imdata_ed<- imdata_ed %>% reduceBaseline(method="locmin")
          }
          
          if (!is.null(preprocess$peakPick$method)){
            imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method)
          }else{
            imdata_ed<- imdata_ed %>% peakPick(method="adaptive")
          }
          
          if ('&'(!is.null(preprocess$peakAlign$tolerance),!is.null(preprocess$peakAlign$tolerance))){
            imdata_ed<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units)
          }else{
            imdata_ed<- imdata_ed %>% peakAlign(tolerance=ppm, units="ppm")
          }
          
          imdata_ed<- imdata_ed %>% process()
          
          if (!is.null(preprocess$normalize)){
            if (preprocess$normalize$method %in% c("rms","tic")){
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
        }
        saveRDS(imdata_ed,paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"),compress = F)  
      }}else{
        imdata_ed<-readRDS(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
      }
    
    if ('&'(preprocess$use_preprocessRDS,file.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS")))){
      message("Using image data: ",paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
      
      imdata_ed<-readRDS(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
      #imdata_ed<-imdata
      if(mzrange=="auto-detect"){
        imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam())
      }else {
        imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,as="MSImagingExperiment",resolution=import_ppm, units="ppm",BPPARAM=SerialParam(),mass.range =mzrange)
      }
    }else{
      message("Using image data: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
      if(mzrange=="auto-detect"){
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
    imdata_org<-imdata  
    imdata<-imdata_ed
    
    coordata=as.data.frame(imdata@elementMetadata@coord)
    
    setwd(paste0(gsub(".imzML$","",datafile[z])  ," ID")) 
    
    
    if (Bypass_Segmentation!=T){
      message("Segmentation in progress...")
      #cl=autoStopCluster(makeCluster(6))
      if (Segmentation[1]=="PCA") {
        if ('&'(Segmentation_ncomp=="auto-detect",Segmentation_variance_coverage>0)){

          Segmentation_ncomp_running<-PCA_ncomp_selection(imdata=imdata,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
        }else{Segmentation_ncomp_running<-Segmentation_ncomp}
      }
      else if (Segmentation[1]=="spatialKMeans" && segmentation_num!=1) {
        if ('&'(Segmentation_ncomp=="auto-detect",Segmentation_variance_coverage>0)){
          Segmentation_ncomp_running<-PCA_ncomp_selection(imdata=imdata,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
        }else{Segmentation_ncomp_running<-Segmentation_ncomp}
        set.seed(1)
        
        skm <-  suppressMessages(suppressWarnings(spatialKMeans(imdata, r=Smooth_range, k=segmentation_num, method="adaptive",ncomp=Segmentation_ncomp_running,BPPARAM =BPPARAM )))
        message(paste0(Segmentation[1], " finished: ",name))
        png(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=F)
        imagefile<-Cardinal::image(skm, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], key=T, ann=FALSE,axes=FALSE)
        print(imagefile)
        dev.off()
        suppressMessages(suppressWarnings(require(magick)))
        skmimg<-image_read(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""))
        
        
        png(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""),width = 1024,height = 480*((segmentation_num)))
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
        
        skmimg_spec<-image_read(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""))
        skmimg<-image_append(c(skmimg,skmimg_spec),stack = T)
        image_write(skmimg,paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs_append.png",sep=""))
        correlation=as.data.frame(skm@resultData@listData[[1]][["correlation"]])
        correlation[,"mz"]<-as.numeric(gsub("m/z = ","",rownames(correlation)))
        centers=as.data.frame(skm@resultData[[1]][["centers"]])
        centers[,"mz"]<-as.numeric(gsub("m/z = ","",rownames(centers)))
        cluster=as.data.frame(skm@resultData[[1]][["cluster"]])
        cluster$Coor=rownames(cluster)
        cluster_df<-data.frame(Coor=rownames(cluster),class=skm@resultData[[1]][["cluster"]])
        ##write.csv(correlation,paste(Segmentation[1],"_RESULT","correlation",segmentation_num,"segs.csv"),row.names = F)
        #write.csv(centers,paste(Segmentation[1],"_RESULT","centers",segmentation_num,"segs.csv"),row.names = F)
        write.csv(cluster_df,paste(Segmentation[1],"_RESULT","cluster",segmentation_num,"segs.csv"),row.names = F)
        
        
        #cluster=skm@resultData[["r = 1, k = 5"]][["cluster"]]
        
        
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
          Segmentation_ncomp_running<-PCA_ncomp_selection(imdata=imdata,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
        }else{Segmentation_ncomp_running<-Segmentation_ncomp}
        #message(paste0("spatialShrunkenCentroids computing for ",name))
        skm <-  suppressMessages(suppressWarnings(spatialShrunkenCentroids(imdata, r=Smooth_range, k=segmentation_num, method="adaptive",s=3,BPPARAM =BPPARAM)))
        message(paste0(Segmentation[1], " finished: ",name))
        png(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=F)
        imagefile<-Cardinal::image(skm, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], key=T, ann=FALSE,axes=FALSE, model=list(s=3), values="class")
        print(imagefile)
        dev.off()
        suppressMessages(suppressWarnings(require(magick)))
        skmimg<-image_read(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs.png",sep=""))
        
        
        png(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""),width = 1024,height = 480*((segmentation_num)))
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
        
        skmimg_spec<-image_read(paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs_spec.png",sep=""))
        skmimg<-image_append(c(skmimg,skmimg_spec),stack = T)
        image_write(skmimg,paste(getwd(),"\\",Segmentation[1],"_image_plot_",segmentation_num,"_segs_append.png",sep=""))
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
        coordistmatrix$sum=base::unlist(bplapply(1:nrow(coordata),function(j,coordistmatrix,coordata){coordistmatrix$sum[j]=sum(coordistmatrix[j,1:nrow(coordata)])},coordistmatrix,coordata,BPPARAM = BPPARAM))
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
        
        
        
        #paste0("v",colnames(coordistmatrix[coordistmatrix$sum==min(coordistmatrix$sum),]))
        
        #center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),])
        
        
        #plot(rownames(coordistmatrix),coordistmatrix$sum)
        
        #write.csv(radius_rank,file = "radius_rank.csv", row.names = F)
        
        
        
        rank_pixel<-function(coordata,coordistmatrix,radius_rank){
          #coordata[coordata$edge==TRUE,]=coordata[coordata$edge==TRUE,]
          
          if (unique(radius_rank$Core)=="central"){
            shape_center=coordata[coordistmatrix$sum==min(coordistmatrix$sum),]
            center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),1:nrow(coordata)])
            library(useful)
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
              
              norm_center_dist=to_center/total*max(radius_rank$Radius_U)
              if(length(radius_rank$Rank['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])==1){
                coordata$rank[i]=as.character(radius_rank$Rank['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
                coordata$pattern[i]=as.character(radius_rank$Name['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
              }else{
                coordata$rank[i]="Undefined"
                coordata$pattern[i]="Undefined"
              }
            }
            coordata$rank<-factor(coordata$rank)
            coordata
          }else { 
            if(unique(radius_rank$Core)=="southeast"){
              #shape_center=coordata['&'(coordata$x>=max(coordata$x)-1,coordata$y>=max(coordata$y)-1),]
              shape_center=coordata[1,]
              shape_center$x=max(coordata$x)
              shape_center$y=max(coordata$y)
              coordata_new<-rbind(coordata,shape_center)
              center_dist=bplapply(nrow(coordata_new),coordist_para,coordata_new,BPPARAM = BPPARAM)[[1]]
              #center_dist=t(coordistmatrix[which.min(coordistmatrix$sum),1:nrow(coordata)])
            } else if(unique(radius_rank$Core)=="northeast"){
              #shape_center=coordata['&'(coordata$x>=max(coordata$x)-1,coordata$y>=max(coordata$y)-1),]
              shape_center=coordata[1,]
              shape_center$x=max(coordata$x)
              shape_center$y=min(coordata$y)
              
            }else if(unique(radius_rank$Core)=="northwest"){
              #shape_center=coordata['&'(coordata$x>=max(coordata$x)-1,coordata$y>=max(coordata$y)-1),]
              shape_center=coordata[1,]
              shape_center$x=min(coordata$x)
              shape_center$y=min(coordata$y)
              
            }else if(unique(radius_rank$Core)=="southwest"){
              #shape_center=coordata['&'(coordata$x>=max(coordata$x)-1,coordata$y>=max(coordata$y)-1),]
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
            center_edge_angle=cbind(coordata[,1:2],cart2pol(df$x, df$y, degrees = F),edge=coordata[,"edge"])
            center_edge_angle_sdge=center_edge_angle[center_edge_angle$edge==TRUE,]
            coordata$rank=0 
            coordata$pattern=""      
            for (i in 1: (nrow(coordata))){
              
              
              #From <- coordata[i,][rep(seq_len(nrow(coordata[i,])), each=nrow(coordata[coordata$edge==TRUE,])),1:2]
              #To <- coordata[coordata$edge==TRUE,][,1:2]
              
              df=coordata[i,1:2]-shape_center[,1:2]
              point_center_angle=cbind(coordata[i,1:2],cart2pol(df$x, df$y, degrees = F))
              pointedge=center_edge_angle_sdge[which(abs(center_edge_angle_sdge$theta-point_center_angle$theta)==min(abs(center_edge_angle_sdge$theta-point_center_angle$theta))),]
              #message(pointedge)
              
              pointedge=pointedge[which.min(pointedge$r),]
              to_edge=coordistmatrix[[i]]['&'(coordata$x==pointedge$x,coordata$y==pointedge$y)]
              
              
              to_center=center_dist[i]
              total=to_edge+to_center
              
              norm_center_dist=to_center/total*max(radius_rank$Radius_U)
              if(length(radius_rank$Rank['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])==1){
                coordata$rank[i]=as.character(radius_rank$Rank['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
                coordata$pattern[i]=as.character(radius_rank$Name['&'(radius_rank$Radius_L<=norm_center_dist,radius_rank$Radius_U>=norm_center_dist)])
              }else{
                coordata$rank[i]="Undefined"
                coordata$pattern[i]="Undefined"
              }
            }
            coordata$rank<-factor(coordata$rank)
            coordata
            
          }
          
          
          
        }
        
        coordata=rank_pixel(coordata,coordistmatrix,radius_rank)
        
        
        x=NULL
        
        for (rank in coordata$rank){
          
          x[[rank]]=which(coordata$rank==rank)
          
        }
        write.csv(coordata,"coordata.csv",row.names = F)
        
        region_pattern <- factor(coordata$pattern,levels=unique(coordata$pattern), labels=unique(coordata$pattern))
        #image(imdata, region_pattern ~ x * y, key=T)
        
        #set.seed(1)
        #skm <- spatialKMeans(imdata, r=Smooth_range, k=length(unique(coordata$rank)), method="adaptive")
        #png(paste(getwd(),"\\","spatialKMeans_image",'.png',sep=""),width = 1024,height = 1024)
        #plot(skm, col=c("pink", "blue", "red","orange","navyblue"), type=c('p','h'), key=FALSE)
        #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1+ceiling(segmentation_num/2), 2),
        #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 1),
        #    bty="n",pty="s",xaxt="n",
        #    yaxt="n",
        #    no.readonly = TRUE,ann=FALSE)
        #print(Cardinal::image(imdata, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], key=FALSE, ann=FALSE,axes=FALSE))
        
        #legend("topright", legend=1:segmentation_num, fill=brewer.pal_n(segmentation_num,colorstyle), col=brewer.pal_n(segmentation_num,"Paired"), bg="transparent",xpd=TRUE,cex = 1)
        
        #Cardinal::plot(skm, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], type=c('p','h'), key=FALSE)
        #Cardinal::plot(skm, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], type=c('p','h'), key=FALSE,mode="centers")
        #Cardinal::plot(skm, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], type=c('p','h'), key=FALSE,mode="betweenss")
        #dev.off()
        
        #for (i in 1:nrow(coordata)){
        #skm@resultData[[paste0("r = ",Smooth_range, ", k = ",length(unique(coordata$rank)))]][["cluster"]][[rownames(coordata)[i]]]=factor(coordata[i,"rank"],levels=unique(coordata$rank))
        
        #}
        
        #Cardinal::image(imdata, col=brewer.pal_n(segmentation_num,colorstyle)[1:segmentation_num], key=FALSE, ann=FALSE,axes=FALSE,groups =pattern)
        
        #msset <- generateImage(region_pattern, coord=coordata[,1:2],
        #                        range=c(1000, 5000), centers=c(2000, 3000, 4000),
        #                        resolution=import_ppm00, step=3.3, as="MSImageSet")
        #msset@pixelData@data[["sample"]]=region_pattern
        png(paste(getwd(),"\\","Virtual_segmentation",gsub("/"," ",name),'.png',sep=""),width = 1024,height = 1024)
        #plot(skm, col=c("pink", "blue", "red","orange","navyblue"), type=c('p','h'), key=FALSE)
        #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1+ceiling(segmentation_num/2), 2),
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 1),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)
        print(image(imdata, region_pattern ~ x * y,col=brewer.pal_n(length(unique(coordata$rank)),colorstyle), key=TRUE))
        
        dev.off()
        
        
        
      }
      else if (Segmentation=="def_file"){
        
        Segmentation_def_tbl<-read.csv(paste0(workdir[z],"/", Segmentation_def))
        
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
        png(paste(getwd(),"\\","Segmentation_def_file","_image_plot_",length(levels(Segmentation_def_tbl$label)),"_segs.png",sep=""),width = 1024,height = 720)
        
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
        png(paste(getwd(),"\\","Segmentation_none","_image_plot_",segmentation_num,"_segs.png",sep=""),width = 1024,height = 720)
        
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 2),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=F)
        imagefile<-Cardinal::image(imdata, factor(rep(1,length(pixels(imdata)))) ~ x * y, key=T, ann=FALSE,axes=FALSE)
        print(imagefile)
        dev.off()
      } 
    }
  }
  return("workflow successfully completed")
}


remove_pep_score_outlier<-function(SMPLIST,IQR_LB=0.75,outputdir=getwd(),abs_cutoff=-2){
  #if (!require(OneR)) install.packages("OneR")
  suppressMessages(suppressWarnings(library(OneR)))
  suppressMessages(suppressWarnings(library(dplyr)))
  suppressMessages(suppressWarnings(library(ggplot2)))
  nbins = floor(length(SMPLIST$mz)/500)
  if(nbins>=2){
    SMPLIST$mzbin<-(bin(SMPLIST$mz, nbins = floor(length(SMPLIST$mz)/500),method = "content"))
  }else{
    message("Insufficient mz features to find the outlier")
    SMPLIST<-SMPLIST[SMPLIST$Score>=abs_cutoff,]
    return(SMPLIST)
  }
  

  remove_outliers <- function(x, na.rm = TRUE,IQR_LB=0.75,...) {
    qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
    H <- IQR_LB * IQR(x, na.rm = na.rm)
    y <- x
    y[x < (qnt[1] - H)] <- NA
    #y[x > (qnt[2] + H)] <- NA
    y
  }


  
  
  SMPLIST_bin_quantile <-SMPLIST %>% group_by(mzbin) %>% summarise(mz=mz,Score=remove_outliers(Score),mzbin=unique(mzbin)) %>% ungroup
  #SMPLIST_bin_1SD <-SMPLIST %>% group_by(mzbin) %>% summarise(mz=min(mz),Score_mean=mean(Score),ScoreSD=sd(Score),ScoreLB= mean(Score)-1*sd(Score),mzbin=unique(mzbin))
  SMPLIST_bin_quantile<-SMPLIST_bin_quantile[!is.na(SMPLIST_bin_quantile$Score),]
  SMPLIST_abs_cutoff<-SMPLIST_bin_quantile[SMPLIST_bin_quantile$Score<abs_cutoff,]
  SMPLIST_bin_quantile<-SMPLIST_bin_quantile[SMPLIST_bin_quantile$Score>=abs_cutoff,]
  
  p<-ggplot(data = SMPLIST, mapping = aes(x=mz,y=Score,color=mzbin,group=mzbin)) + 
    #geom_jitter(aes(color='blue'),alpha=0.2) +
    #geom_boxplot(fill="bisque",color="black",alpha=0.3,coef=0.75) + 
    #stat_boxplot(geom ='errorbar')+
    #geom_point(data=SMPLIST_bin_1.5SD,mapping = aes(x=mz,y=ScoreLB),color="blue",group="Cut-off") +
    geom_point(data=SMPLIST,mapping = aes(x=mz,y=Score),
               color="darkblue"
               ,group="Cut-off",col="Removed Features") +
    geom_point(data=SMPLIST_bin_quantile,mapping = aes(x=mz,y=Score),
               color="red",
               group="Cut-off",col="Kept Features") +
    # geom_point(data=SMPLIST_abs_cutoff,mapping = aes(x=mz,y=Score),
    #            color="green",
    #            group="Cut-off",col="Kept Features") +
    #stat_smooth(data=SMPLIST_bin_1SD,mapping = aes(x=mz,y=ScoreLB),method = "lm", formula = y ~ poly(x, 4), se = FALSE,fill="black") +
    guides(color=guide_legend("m/z bin size")) +
    labs(x='m/z') +
    
    theme_classic() 
  
  
  if(!is.null(outputdir)){
    png(paste(outputdir,"/Summary folder/Peptide_score_outlier.png",sep=""),width = 1200,height = 800,res = 150)
    print(p)
    dev.off()
    
  }
  return(SMPLIST[SMPLIST$mz %in% SMPLIST_bin_quantile$mz,])
}
