
Load_Cardinal_imaging<-function(datafile=tk_choose.files(filter = Filters,
                                                         caption  = "Choose single or multiple file(s) for analysis",
                                                         multi = F),
                                BPPARAM = bpparam(),
                                resolution=2.5,
                                attach.only=TRUE,
                                preprocessing=F,
                                rotate=0,as="MSImagingExperiment",
                                mzrange=NULL){
  #imdata_meta <- importImzMl(datafile, coordinates = matrix(c(2, 4),nrow=1, ncol=2),removeEmptySpectra = F, centroided = T)
  
  datafiles<-gsub(".imzML", "", datafile)
  workdir<-base::dirname(datafiles) 
  name <-gsub(base::dirname(datafiles),"",datafiles)
  folder<-base::dirname(datafiles)
  imdata <- readImzML(name, folder, attach.only=attach.only,as=as,resolution=resolution, units="ppm",rotate = rotate,BPPARAM=BPPARAM,mass.range=mzrange)
  if (preprocessing){
    imdata <- try(batchProcess(imdata, normalize=FALSE, smoothSignal=TRUE, reduceBaseline=list(method = "median",blocks=500, fun=min, spar=1),
                               peakPick=list(SNR=12), peakAlign=TRUE,BPPARAM=BPPARAM))
  }
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
    
    outputpngsum=paste(getwd(),"\\",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot_sum",'.png',sep="")
    
    
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
    
    outputpngsum=paste(getwd(),"\\",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot_sum_flex",'.png',sep="")
    
    
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
    candidate_unique_table=unique(candidate[,c("mz","FA","Formula","moleculeNames" , "adduct")])
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

cluster_image_grid<-function(clusterID,
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
                                          export_Header_table=F,
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
        
        outputpngsum=paste(getwd(),"\\",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot_sum",'.png',sep="")
        
        
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
        darkmode()
        tmp_dir <- tempdir()
        outputpngsum=paste(getwd(),"\\",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot_sum_flex",'.png',sep="")
        temp_cluster_png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        temp_component_png=list()
        if (plot_layout=="line"){
          png(temp_cluster_png,width = 5,height = 5 * length( levels(run(imdata))), bg = "black",units = "in",res = 300)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(run(imdata))),1),
          #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)   
          
          
        }else{
          png(temp_cluster_png,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 300)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)  
        }
        
        
        
        
        clusterimg=image(imdata, mz=candidateunique, 
              col=levels(mycol),
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              plusminus=median(ppm*candidateunique/1000000),
              layout=c( length(levels(run(imdata))),1))
        print(clusterimg)
        
        dev.off()
        
        
        
        componentimg=list()
        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- gradient.colors(100, start="black", end=levels(mycol)[i])
        componentimg[[i]]=image(imdata, mz=candidateunique[i], 
                contrast.enhance=contrast.enhance,
                smooth.image = smooth.image,#smooth.image ,
                col.regions=intensity.colors_customize1(),
                normalize.image="none",
                plusminus=ppm*candidateunique[i]/1000000,
                key=F,
                xlab=NULL,
                ylab=NULL,
                layout=c( length(levels(run(imdata))),1)
                #xlim=c(0,50),
                #ylim=c(0,40)
          )
        temp_component_png[[i]]=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        png(temp_component_png[[i]],width = 5,height = 5 * length( levels(run(imdata))), bg = "black",units = "in",res = 300)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(run(imdata))),1),
            #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)  
        print(componentimg[[i]])
        dev.off()
          #componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          #for (component in componentname){
          #  text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          #}
        }

        pngfile<-image_read(temp_cluster_png)
        pngfile<-image_border(pngfile, "black", "30x30")
        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "north",size = 50,color = "white")
        
        
        pngcompfile=list()
        for (i in 1:length(candidateunique)){
          pngcompfile[[i]]<-image_read(temp_component_png[[i]])
          pngcompfile[[i]]<-image_border(pngcompfile[[i]], "black", "30x30")
          pngcompfile[[i]]<-image_annotate(pngcompfile[[i]],paste(unique(candidate[candidate$mz==candidateunique[i],"moleculeNames"]),candidateunique[i]),gravity = "north",size = 50,color = "white")
         
          }
         pngfile=image_append(c(pngfile,unlist(pngcompfile)))
        image_write(pngfile,outputpngsum)
        
        removetempfile(tmp_dir,matchword=c(".png$","^magick-"))
        
        rm(pngcompfile)
      
    }
    
    
    
    
  }
  
  if(export_Header_table){
    candidate_unique_table=unique(candidate[,c("mz","FA","Formula","moleculeNames" , "adduct")])
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
      width=120*(nrow(Header_table)+10),
      height=400,
      
      
    ) %>% layout(
      title = paste0("<b>",clusterID,"</b>"),
      autosize = F,
      margin=0,
      font = list(family = "Arial",  color = "black",align = "bottom")
    )
    
    
    
    
    
    orca(p, file = windows_filename(paste0(clusterID,"header.png")),width=120*(nrow(Header_table))+200,height=540) 
    
  }
  
  
  
  
}}

orca_initial<-function(){
  
Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
  
  
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
  
  library(gridGraphics)
  
  grab_grob <- function(){
    grid.echo()
    grid.grab()
  }
  
  g <- grab_grob()
  grid.newpage()
  pushViewport(viewport(width=0.7,angle=30))
  grid.draw(g)
}

#### Read imzML files ####
## ----------------------

readImzML <- function(name, folder = getwd(), attach.only = TRUE,
                      mass.range = NULL, resolution = 200, units = c("ppm", "mz"),
                      as = c("MSImagingExperiment", "MSImageSet"), parse.only=FALSE,
                      BPPARAM = bpparam(),rotate=0, ...)
{
  library(Cardinal)
  # check input
  dots <- list(...)
  if ( "mass.accuracy" %in% names(dots) ) {
    .stop("'mass.accuracy' is defunct.\n",
          "Use 'resolution' instead.")
  }
  if ( "units.accuracy" %in% names(dots) ) {
    .stop("'units.accuracy' is defunct.\n",
          "Use 'units' instead.")
  }
  # get output format
  outclass <- match.arg(as)
  # check for files
  xmlpath <- normalizePath(file.path(folder, paste(name, ".imzML", sep="")),
                           mustWork=FALSE)
  if ( !file.exists(xmlpath) ) .stop("expected file ", xmlpath, " does not exist")
  ibdpath <- normalizePath(file.path(folder, paste(name, ".ibd", sep="")),
                           mustWork=FALSE)
  if ( !file.exists(ibdpath) ) .stop("expected file ", ibdpath, " does not exist")
  # read imzML file
  message("reading imzML file: '", xmlpath, "'")
  info <- .readImzML(xmlpath)
  coord=data.frame(x=info@scanList@listData[["position x"]],y=info@scanList@listData[["position y"]])
  newcoord=affine(coord,rotate = rotate)
  
  newcoord$x=newcoord$x-min(newcoord$x)+1
  newcoord$y=newcoord$y-min(newcoord$y)+1
  info@scanList@listData[["position x"]]=as.integer(as.character(newcoord$x))
  info@scanList@listData[["position y"]]=as.integer(as.character(newcoord$y))
  info@metadata[["max count of pixels x"]]=max(newcoord$x)
  info@metadata[["max count of pixels y"]]=max(newcoord$y)
  
  
  if ( parse.only )    return(info)
  # read ibd file
  info@metadata[["files"]] <- c(xmlpath, ibdpath)
  info@metadata[["name"]] <- name
  units <- match.arg(units)
  message("reading ibd file: '", ibdpath, "'")
  object <- .readIbd(ibdpath, info, outclass=outclass, attach.only=attach.only,
                     mass.range=mass.range, resolution=resolution, units=units, BPPARAM=BPPARAM)
  #.log.collapse("loaded dataset:", capture.output(print(object)))
  if ( validObject(object) ) {
    message("done.")
    object
  }
}

.readImzML <- function(file) {
  parse <- .Call("C_readImzML", normalizePath(file), PACKAGE="Cardinal")
  len <- sapply(parse$experimentMetadata, nchar, type="bytes")
  experimentMetadata <- parse$experimentMetadata[len > 0]
  .MSImagingInfo(
    scanList=as(parse$scanList, "DataFrame"),
    mzArrayList=as(parse$mzArrayList, "DataFrame"),
    intensityArrayList=as(parse$intensityArrayList, "DataFrame"),
    metadata=experimentMetadata)
}
intensityData<-function(x){
  x@intensityData
}

metadata<-function(x){
  x@metadata
}
.readIbd <- function(file, info, outclass, attach.only,
                     mass.range, resolution, units,...)
{
  library(matter)
  file <- normalizePath(file)
  ibdtype <- metadata(info)[["ibd binary type"]]
  mz.ibdtype <- mzData(info)[["binary data type"]]
  intensity.ibdtype <- imageData(info)[["binary data type"]]
  # read binary data
  if ( ibdtype == "continuous" ) {
    mz <- matter::matter_vec(paths=file,
                     datamode=Ctypeof(mz.ibdtype[1]),
                     offset=mzData(info)[["external offset"]][1],
                     extent=mzData(info)[["external array length"]][1])
    intensity <- matter::matter_mat(paths=file,
                            datamode=Ctypeof(intensity.ibdtype[1]),
                            offset=imageData(info)[["external offset"]],
                            extent=imageData(info)[["external array length"]])
    if ( attach.only ) {
      spectra <- intensity
    } else {
      spectra <- intensity[]
    }
    mz <- mz[]
  } else if ( ibdtype == "processed" ) {
    mz <- matter::matter_list(paths=file,
                      datamode=Ctypeof(mz.ibdtype),
                      offset=mzData(info)[["external offset"]],
                      extent=mzData(info)[["external array length"]])
    intensity <- matter::matter_list(paths=file,
                             datamode=Ctypeof(intensity.ibdtype),
                             offset=imageData(info)[["external offset"]],
                             extent=imageData(info)[["external array length"]])
    if ( is.null(mass.range) ) {
      mzvec <- as(mz, "matter_vec")
      chunksize(mzvec) <- 1e8L # read chunks of 800 MB
      mz.range <- range(mzvec)
    } else {
      mz.range <- mass.range
    }
    mz.min <- mz.range[1]
    mz.max <- mz.range[2]
    if ( units == "ppm" ) {
      if ( floor(mz.min) <= 0 )
        .stop("readImzML: m/z values must be positive for units='ppm'")
      mzout <- seq.ppm(
        from=floor(mz.min),
        to=ceiling(mz.max),
        ppm=resolution) # ppm == half-bin-widths
      error <- resolution * 1e-6 * mzout
      tol <- c(relative = resolution * 1e-6)
    } else {
      mzout <- seq(
        from=floor(mz.min),
        to=ceiling(mz.max),
        by=resolution * 2)  # by == full-bin-widths
      error <- rep(resolution, length(mzout))
      tol <- c(absolute = resolution)
    }
    mz.bins <- c(mzout[1] - error[1], mzout + error)
    if ( attach.only ) {
      data <- list(keys=mz, values=intensity)
      mz <- mzout
      spectra <- sparse_mat(data, keys=mz,
                            nrow=length(mz), ncol=length(intensity),
                            tolerance=tol, combiner="sum")
    } else {
      if ( outclass == "MSImageSet" ) {
        data <- list(keys=list(), values=list())
        for ( i in seq_along(mz) ) {
          mzi <- mz[[i]]
          wh <- findInterval(mzi, mz.bins)
          s <- as.vector(tapply(intensity[[i]], wh, sum))
          data$keys[[i]] <- mzout[unique(wh)]
          data$values[[i]] <- s
        }
      } else if ( outclass == "MSImagingExperiment") {
        data <- list(keys=mz[], values=intensity[])
      }
      mz <- mzout
      spectra <- sparse_mat(data, keys=mz,
                            nrow=length(mz), ncol=length(intensity),
                            tolerance=tol, combiner="sum")
    }
  }
  # set up coordinates
  x <- scans(info)[["position x"]]
  y <- scans(info)[["position y"]]
  z <- scans(info)[["position z"]]
  x3d <- scans(info)[["3DPositionX"]]
  y3d <- scans(info)[["3DPositionY"]]
  z3d <- scans(info)[["3DPositionZ"]]
  if ( all(is.na(z)) && all(is.na(z3d)) ) {
    coord <- data.frame(x=x, y=y)
  } else if ( all(is.na(z3d)) ) {
    coord <- data.frame(x=x, y=y, z=z)
  } else {
    z <- as.integer(as.factor(z3d))
    coord <- data.frame(x=x, y=y, z=z)
  }
  if ( outclass == "MSImageSet" ) {
    experimentData <- new("MIAPE-Imaging")
    processingData <- new("MSImageProcess", files=metadata(info)[["files"]])
    object <- MSImageSet(spectra=spectra, mz=mz, coord=coord,
                         processingData=processingData,
                         experimentData=experimentData)
    sampleNames(object) <- metadata(info)[["name"]]
  } else if ( outclass == "MSImagingExperiment" ) {
    centroided <- isTRUE(spectrumRepresentation(info) == "centroid spectrum")
    object <- MSImagingExperiment(spectra,
                                  featureData=MassDataFrame(mz=mz),
                                  pixelData=PositionDataFrame(coord=coord,
                                                              run=metadata(info)[["name"]]),
                                  metadata=metadata(info),
                                  centroided=centroided)
  } else {
    stop("unrecognized outclass")
  }
  object
}

.MSImagingInfo <- setClass("MSImagingInfo",
                           contains = "Vector",
                           slots = c(
                             scanList = "DataTable",
                             mzArrayList = "DataTable",
                             intensityArrayList = "DataTable"),
                           prototype = prototype(
                             scanList = DataFrame(),
                             mzArrayList = DataFrame(),
                             intensityArrayList = DataFrame()))
.message <- function(..., progress=c("none", "start", "stop", "increment"), min=0, max=1) {
  progress <- match.arg(progress)
  if ( progress == "none" ) {
    .log(...)
    for ( f in .Cardinal$message ) {
      f(...)
    }
  } else if ( progress == "start" ) {
    if ( length(list(...)) > 1 )
      .log(...)
    .Cardinal$progress$i <- min
    .Cardinal$progress$min <- min
    .Cardinal$progress$max <- max
    for ( f in .Cardinal$progress$start ) {
      f(..., min=min, max=max)
    }
  } else if ( progress == "increment" ) {
    .Cardinal$progress$i <- .Cardinal$progress$i + 1
    for ( f in .Cardinal$progress$increment ) {
      f()
    }
  } else if ( progress == "stop" ) {
    if ( length(list(...)) > 1 )
      .log(...)
    for ( f in .Cardinal$progress$stop ) {
      f(...)
    }
  }
}
.log <- function(...) {
  msg <- paste(date(), paste0(..., collapse="\n  "))
  .Cardinal$log <- append(.Cardinal$log, msg)
  elapsed <- proc.time()[3] - .Cardinal$time$flush
  if ( elapsed > getOption("Cardinal.flush") )
    .log.flush()
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
  } else {
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
  file.remove(paste0(temp_dir,"\\",filelist[test]))
  
  
}

testcolume<-function(df,testcolnames){
  
  library(stringr)
  
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
  Formula_add=c("H1","N1H4","Na1","K1","FALSE","FALSE","FALSE","FALSE","C1O2H2","C2O2H4","FALSE","H3",
                "H2Na1","H1Na2","Na3","H2","H1N1H4","H1Na1","H1K1","C2H5N1","Na2","C4H8N2","C6H11N3",
                "C1H5O1","C2H4N1","Na2","C3H9O1","C2H3N1Na1","K2","C2H7S1O1","C4H7N2","C3H9O1Na1","H1",
                "N1H4","Na1","H8O6","K1","C2H4N1","C2H3N1Na1","FALSE","Na1","Cl1","K1","Br1","C2F3O2H1",
                "FALSE","C1O2H2","C2O2H4","FALSE","He","Ne","Ar","Kr","Xe","Rn","Cu","Co","Ag"
  )
  Formula_ded=c("FALSE","FALSE","FALSE","FALSE","FALSE","H1","H2","H3","H1","H1","FALSE","FALSE","FALSE",
                "FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE",
                "FALSE","H1","FALSE","FALSE","H1","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE",
                "FALSE","FALSE","FALSE","H3O1","H2","FALSE","H2","FALSE","H1","H1","H1","H1","H1","FALSE",
                "FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE","FALSE")
  Multi=c("1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1","1"
          ,"1","1","1","1","1","1","1","1","2","2","2","2","2","2","2","1","1","1","1","1","1","2","2","2",
          "3","1","1","1","1","1","1","1","1","1")
  adductslist<-cbind.data.frame(Name,calc,Charge,Mult,Mass,Ion_mode,Formula_add,Formula_ded,Multi)
  adductslist
}

rank_mz_feature<-function(Peptide_plot_list,mz_feature,BPPARAM=bpparam()){
  
  library(data.table)
  
  mz=unique(Peptide_plot_list$mz)
  
  Peptide_plot_list$mz_align<-NULL
  
  mz_align<-sapply(mz,function(x,mz_feature){
    
    mz_feature$m.z[which.min(abs(mz_feature$m.z-x))]
    
  },mz_feature)
  
  mz_feature_cluster<-data.table(mz=mz,mz_align=mz_align)
  
  Peptide_plot_list<-merge(Peptide_plot_list,mz_feature_cluster,by="mz",all=T)
  
  mz=unique(Peptide_plot_list$mz_align)
  
  Peptide_plot_list_rank<-bplapply(mz,function(x,Peptide_plot_list){
    
    randklist <- Peptide_plot_list[Peptide_plot_list$mz_align==x,]
    
    randklist$Rank<- as.numeric(factor(randklist$Score,levels=sort(unique(randklist$Score),decreasing = T)))
    
    return(randklist)
    
  },Peptide_plot_list,BPPARAM=BPPARAM)
  
  Peptide_plot_list_rank<-do.call(rbind,Peptide_plot_list_rank)
  
  return(Peptide_plot_list_rank)
}
