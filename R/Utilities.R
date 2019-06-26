Load_Cardinal_imaging<-function(datafile=tk_choose.files(filter = Filters,
                                                         caption  = "Choose single or multiple file(s) for analysis",
                                                         multi = F),
                                resolution=2.5,
                                attach.only=TRUE,
                                preprocessing=F){
  #imdata_meta <- importImzMl(datafile, coordinates = matrix(c(2, 4),nrow=1, ncol=2),removeEmptySpectra = F, centroided = T)
  datafiles<-gsub(".imzML", "", datafile)
  workdir<-base::dirname(datafiles) 
  name <-gsub(base::dirname(datafiles),"",datafiles)
  folder<-base::dirname(datafiles)
  imdata <- readImzML(name, folder, attach.only=attach.only,as="MSImageSet",resolution=resolution, units="ppm")
  if (preprocessing){
    imdata <- try(batchProcess(imdata, normalize=FALSE, smoothSignal=TRUE, reduceBaseline=list(method = "median",blocks=500, fun=min, spar=1),
                               peakPick=list(SNR=12), peakAlign=TRUE))
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
  if(rotate_image){
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
          plusminus=ppm)
    
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
          plusminus=ppm)
    
    for (i in 1:length(candidateunique)){
      #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
      col.regions <- gradient.colors(100, start="black", end=levels(mycol)[i])
      image(imdata, mz=candidateunique[i], 
            contrast.enhance=contrast.enhance,
            smooth.image = smooth.image ,
            col.regions=col.regions,
            normalize.image="linear",
            plusminus=ppm)
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
    
    image_rotate<-function(...){
      library(gridGraphics)
      library(magick)
      temp.png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
      png(temp.png,width = 5,height = 5, bg = "black",units = "in",res = 300)
      image(...)
      dev.off()
      image_read(temp.png)
      
      
      grab_grob <- function(...){
        grid.echo()
        grid.grab()
      }
      
      g <- grab_grob()
      grid.newpage()
      pushViewport(viewport(width=0.7,angle=30))
      
      grid.draw(g)
      
    }
    
    
    image(imdata, mz=candidateunique, 
          col=levels(mycol),
          contrast.enhance = contrast.enhance,
          smooth.image = smooth.image ,
          superpose=TRUE,normalize.image="linear",
          plusminus=ppm)
    

    for (i in 1:length(candidateunique)){
      #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
      col.regions <- gradient.colors(100, start="black", end=levels(mycol)[i])
      image(imdata, mz=candidateunique[i], 
            contrast.enhance=contrast.enhance,
            smooth.image = smooth.image,#smooth.image ,
            col.regions=intensity.colors_customize1(),
            normalize.image="none",
            plusminus=ppm,
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

orca_initial<-function(){
  
Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
  
  
}

affine <- function(x, translate=c(0,0), rotate=180,
                   angle=c("degrees", "radians"), grid=TRUE)
{
  x=x[,c("x","y")]
  angle <- match.arg(angle)
  theta <- -rotate
  if ( angle == "degrees" ) theta <- theta * pi / 180
  # translate center of mass to be near origin
  tt <- sapply(x, function(xs) mean(xs))
  new.x <- t(as.matrix(x)) - tt
  # rotate around origin
  A <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), nrow=2)
  new.x <- A %*% new.x
  # translate back and do requested translation
  new.x <- t(new.x + tt + translate)
  # remove negative coordinates and round to integers
  if ( grid ) {
    new.x <- round(new.x)
    new.x[new.x < 1] <- 1
  }
  # return data.frame of new coordinates
  new.x <- as.data.frame(new.x)
  names(new.x) <- names(x)
  new.x
}
