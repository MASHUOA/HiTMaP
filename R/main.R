
#' \tabular{ll}{
#' Package: \tab Metwork \cr
#' License: \tab GPL (>= 2)\cr
#' URL: \tab 
#' }
#'
#' @docType package
#' @name Metwork
#' @author George GUO \email{George.GUO@@auckland.ac.nz}
#' @references TBD
#' @keywords package
#' @import tcltk
#' @import Rcpp
#' @import tidyr
#' @import randomForest
#' @import rfviz
#' @import extrafont
#' @import mzR
#' @import xcms
#' @import KEGGREST
#' @import RColorBrewer
#' @import rgl
#' @import dplyr
#' @import tkrplot
#' @import multtest
#' @import XML
#' @import CAMERA
#' @import doParallel
#' @import BiocManager
#' @import pacman
#' @import BiocParallel
#' @import spdep
#' @import ggplot2
#' @import rcdk
#' @import magick
#' @import enviPat
#' @import grid
#' @import ggplot2
#' @import zoo
#' @import FTICRMS
#' @import OrgMassSpecR
#' @importFrom stats na.omit runif
#' @importFrom utils download.file modifyList packageVersion read.table tail
#' @import S4Vectors



Generatespectrum_workflow_maldiquant<-function(){
  table(sapply(MALDI_IMAGE, length))
  plot(MALDI_IMAGE[[1]])
  
  workdir <- WorkingDir()
  datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
  datafile<-gsub(".imzML", "", datafile)
  MALDI_IMAGE <- quiet(importImzMl(paste(file.path( datafile[i]),".imzml",sep="")))
  #use the square root transformation to simplify graphical visualization and to overcome the potential dependency of the variance from the mean.
  spectra <- transformIntensity(MALDI_IMAGE, method="sqrt")
  #use a 21 point Savitzky-Golay-Filter (Savitzky and Golay, 1964) to smooth the spectra.
  spectra <- smoothIntensity(spectra, method="SavitzkyGolay", halfWindowSize=10)
  #Before we correct the baseline we visualize it. Here we use the SNIP algorithm (Ryan et al., 1988).
  baseline <- estimateBaseline(spectra[[16]], method="SNIP", iterations=100)
  plot(spectra[[16]])
  lines(baseline, col="red", lwd=2)
  #If we are satisfied with our estimated baseline we remove it.
  spectra <- removeBaseline(spectra, method="SNIP",iterations=100)
  plot(spectra[[1]])
  #spectra <- calibrateIntensity(spectra, method="TIC")
  #do the alignment if there's a huge mass shift between each spectrum
  spectra <- alignSpectra(spectra,halfWindowSize=50,SNR=3,tolerance=0.05,warpingMethod="lowess")
  
  samples <- factor(sapply(spectra, function(x)metaData(x)$sampleName))
  
  avgSpectra <- averageMassSpectra(spectra, labels=samples, method="mean")
  
  noise <- estimateNoise(avgSpectra[[1]])
  plot(avgSpectra[[1]], xlim=c(500, 4000), ylim=c(0, 0.002))
  lines(noise, col="red")
  lines(noise[,1], noise[, 2]*2, col="blue")
  peaks <- detectPeaks(avgSpectra, method="MAD",halfWindowSize=20, SNR=2)
  plot(avgSpectra[[1]], xlim=c(4000, 5000), ylim=c(0, 0.002))
  points(peaks[[1]], col="red", pch=4)
  peaks <- binPeaks(peaks, tolerance=0.002)
  peaks <- filterPeaks(peaks, minFrequency=0.25)
  
}


Cardinal_utilities<- function(){
  folder <- WorkingDir()
  datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
  
  listfile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.csv$")
  datafile<-gsub(".imzML", "", datafile)
  name <-gsub(paste(folder,"\\/",sep=""),"",datafile[1])
  
  #imdata <- readImzML(name, folder, as="MSImagingExperiment")
  imdata <- readImzML(name, folder, attach.only=TRUE, as="MSImagingExperiment")
  summarize(data, sum, .by="pixel")
  tmp <- data %>%
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

####load pathway database
#Map_SMP_MP <- read.csv("DB//smpdb_pathways.csv",as.is = TRUE)
####load raw data annotation

####load adducts data
#adducts<-read.csv("DB//adducts.csv",as.is = TRUE)
####load SMPDB data
#SMPDB <- read.csv(file="DB/SMPdb.csv",stringAsFactors=FALSE,as.is = TRUE)
#SMPDBFREQ <- read.csv(file="DB/SMPdb_ID_FREQ.csv",col.names= c("","SMPDB.ID","Freq"),as.is = TRUE)
##work dir

##load raw data file
#MALDI_IMAGE <- importImzMl(file.path( "NEDC.imzML"))
######Example steps

metwork_maldiquant<- function(workdir= tk_choose.dir(caption = "Select working directory")){  
  
datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
datafile<-gsub(".imzML", "", datafile)
setwd("Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy")
HMDBids <- read.csv("combine.csv", fill = TRUE,skip=2,as.is = TRUE)
SMPDB <- read.csv(file="DB/SMPdb.csv",as.is = TRUE)

for (i in 1: length(datafile)){

  if (dir.exists(datafile[i])==FALSE){dir.create(datafile[i])}
  HMDBDATA <- Generate_HMDBIDS(HMDBids)
  SMPLIST <- META_pathway_search(HMDBDATA, SMPDB, datafile[i], PWtype="Metabolic")
  SMPRESULT <-read.csv(file=paste(datafile[i],sep="","\\","SMPRESULT.csv"),as.is = TRUE)
  options(warn=-1)
  MALDI_IMAGE <- quiet(importImzMl(paste(file.path( datafile[i]),".imzml",sep="")))
  
  options(warn=0)
  META_pathway_ion_image(SMPLIST, SMPRESULT, datafile[i] ,MALDI_IMAGE, Image_Type="",PWpScore=0.26,plotTolerance=0.125 )
  META_pathway_Whole_Picture(SMPLIST, SMPRESULT, datafile[i], PWpScore=0.26,pathway_mode="")
  META_pathway_Whole_Picture_discovery_mode(SMPLIST, SMPRESULT, datafile[i], SMPDB, adductslist=c("M-H","M+Cl"), PWpScore=0.26)
}
}

metwork_cardinal<- function(workdir= tk_choose.dir(caption = "Select working directory")){  
  
  datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
  datafile<-gsub(".imzML", "", datafile)
  setwd("Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy")
  HMDBids <- read.csv("combine.csv", fill = TRUE,skip=2,as.is = TRUE)
  SMPDB <- read.csv(file="DB/SMPdb.csv",as.is = TRUE)
  
  for (i in 1: length(datafile)){
    
    if (dir.exists(datafile[i])==FALSE){dir.create(datafile[i])}
    HMDBDATA <- Generate_HMDBIDS(HMDBids)
    SMPLIST <- META_pathway_search(HMDBDATA, SMPDB, datafile[i], PWtype="Metabolic")
    SMPRESULT <-read.csv(file=paste(datafile[i],sep="","\\","SMPRESULT.csv"),as.is = TRUE)
    options(warn=-1)
    imdata <- quiet(Load_Cardinal_imaging(paste(file.path( datafile[i]),sep=""),resolution = 20))
    
    options(warn=0)
    META_pathway_ion_image_cardinal(SMPLIST, SMPRESULT, datafile[i] ,imdata, Image_Type="",PWpScore=0.26,plotTolerance=20 ,)
    filename=gsub(dirname(datafile[i]),"",datafile[i])
    for (SMPDB.ID in SMPRESULT$SMPDB.ID[SMPRESULT$pScore>=PWpScore]){
      
      cluster_image_cardinal(SMPDB.ID,imdata,SMPLIST,SMPRESULT$Name[SMPRESULT$SMPDB.ID==SMPDB.ID],ppm=ppm)
      #cluster_image_cardinal(SMPDB.ID,imdata_istd,SMPLIST,SMPRESULT$Name[SMPRESULT$SMPDB.ID==SMPDB.ID])
    }
    #lapply(SMPRESULT$SMPDB.ID[SMPRESULT$pScore>=PWpScore],cluster_image_cardinal,imdata,SMPLIST)
    
    META_pathway_Whole_Picture(SMPLIST, SMPRESULT, datafile[i], PWpScore=0.26,pathway_mode="")
    META_pathway_Whole_Picture_discovery_mode(SMPLIST, SMPRESULT, datafile[i], SMPDB, adductslist=c("M-H","M+Cl"), PWpScore=0.26)
  }
}
##get HMDB IDs list for pathway search
intensity.colors_customize_mono <- function(n = 100, alpha = 1,colorstart="transparent",colorend="red") {
  
  col1 <- colorRamp(colors=c(colorstart, colorend),interpolate = "linear")
 
  col1
}

cluster_image_cardinal<-function(clusterID,imdata,SMPLIST,clustername,ppm=20){
  library(RColorBrewer)
  outputpng=paste(getwd(),"\\",clustername,"_cluster_plot",'.png',sep="")  
  
  candidate=SMPLIST[SMPLIST$SMPDB.ID==clusterID,]
  candidateunique=unique(candidate[,"mz"])
  mycol <- RColorBrewer::brewer.pal(length(candidateunique),"Set1")
  png(outputpng,width = 720,height = 720, bg = "transparent")
  
  image(imdata, mz=candidateunique, 
        col=mycol,
        #contrast.enhance = contrast.enhance,
        smooth.image = smooth.image ,
        superpose=TRUE,normalize.image="linear")
  
  
  dev.off()
  pngfile<-image_read(outputpng)
  pngfile<-image_border(pngfile, "transparent", "30x30")
  pngfile<-image_annotate(pngfile,paste(clustername),gravity = "south",size = 40,color = "white")
  pngfile<-image_trim(pngfile)
  image_write(pngfile,outputpng)
  
  
for (i in 1:length(candidateunique)){
    outputpngsum=paste(getwd(),"\\",clustername,"_cluster_plot_",candidateunique[i],'.png',sep="")
  png(outputpngsum,width = 720,height = 720, bg = "transparent")
  par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1, 1),
      bty="n",pty="s",xaxt="n",
      yaxt="n",
      no.readonly = TRUE,ann=FALSE)
  #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
    col.regions <- gradient.colors(100, start="transparent", end=mycol[i])
    image(imdata, mz=candidateunique[i], 
          contrast.enhance=contrast.enhance,
          smooth.image = smooth.image ,
          col.regions=col.regions,normalize.image="linear")
    dev.off()
    }
  
}

#plot cluster image all in one file



Generate_HMDBIDS <- function(HMDBids){
  HMDBDATA<-NULL
  HMDBSearchlist<-str_split(paste(HMDBids$moleculeIds,collapse=", "),pattern=", ")[[1]]
  Formular<-str_split(paste(HMDBids$moleculeNames,collapse=", "),pattern=", ")[[1]]
  HMDBDATA$HMDB.ID<-HMDBSearchlist
  HMDBDATA$moleculeNames<-Formular
  MZLIST <- NULL
  adduct <- NULL
  options(warn=-1)
  for (i in 1: length(HMDBSearchlist)){
    MZLIST[i]<-as.numeric(HMDBids[grep(HMDBSearchlist[i],HMDBids[,"moleculeIds"]),"mz"])
    adduct[i]<-as.character(HMDBids[grep(HMDBSearchlist[i],HMDBids[,"moleculeIds"]),"adduct"])
  }
  HMDBDATA$mz<-MZLIST
  HMDBDATA$adduct<-adduct
  options(warn=0)
  HMDBDATA<-as.data.frame(HMDBDATA)
  HMDBDATA<-unique(HMDBDATA)
  return(HMDBDATA)
  
}

##simple_ion_image: define the work dir read all csv and enerate the ion image respectively
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
  MALDI_IMAGE <- importImzMl(paste(file.path(datafile[z]),".imzML",sep=""))
  options(warn=0)

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
Plot_Ion_image_Png3(outputDir, MALDI_IMAGE, windows_filename(str_split(masslist$moleculeNames[j],", ")[[1]][k]), masslist$mz[j],masslist$adduct[j], Tolerance= plotTolerance*masslist$mz[j]/1000000, title="",Creat_new_file=TRUE,Neighbour = Denoising,color = color,interpolate =interpolate) 
}}}
}
}
}

WorkingDir <-function(Mainfolder = tk_choose.dir(caption = "Select working directory")){
  Mainfolder <- Mainfolder
  setwd(Mainfolder)
  Mainfolder
}

##Pathway search return search result
META_pathway_search <- function(HMDBDATA, SMPDB, WKdir, PWtype="Metabolic", ...) {
  SMPLIST<- NULL
  SMPMZlist<-NULL
  for (i in 1: length(HMDBDATA[,"HMDB.ID"])){
    SMPquery<-NULL
    
    SMPquery<-SMPDB[`&`(SMPDB[,"Pathway.Subject"]==PWtype,SMPDB[,"HMDB.ID"] == HMDBDATA[,"HMDB.ID"][i]),"SMPDB.ID"]
    SMPLIST$SMPDB.ID<-c(SMPLIST$SMPDB.ID,SMPquery)
    SMPLIST$HMDB.ID<-c(SMPLIST$HMDB.ID, rep(as.character(HMDBDATA[,"HMDB.ID"][i]),length(SMPquery)))
    SMPLIST$mz<-c(SMPLIST$mz, rep(HMDBDATA[,"mz"][i],length(SMPquery)))
    SMPLIST$moleculeNames<-c(SMPLIST$moleculeNames, rep(as.character(HMDBDATA[,"moleculeNames"][i]),length(SMPquery)))
    SMPLIST$KEGG.ID<-c(SMPLIST$KEGG.ID, rep(SMPDB[SMPDB[,"HMDB.ID"]==HMDBDATA[,"HMDB.ID"][i],"KEGG.ID"][1],length(SMPquery)))
    SMPLIST$Metabolite.ID<-c(SMPLIST$Metabolite.ID, rep(SMPDB[SMPDB[,"HMDB.ID"]==HMDBDATA[,"HMDB.ID"][i],"Metabolite.ID"][1],length(SMPquery)))
    SMPLIST$adduct<-c(SMPLIST$adduct, rep(as.character(HMDBDATA[HMDBDATA[,"HMDB.ID"]==HMDBDATA[,"HMDB.ID"][i],"adduct"][1]),length(SMPquery)))
  }
  if (dir.exists(paste(WKdir,sep="","\\"))==FALSE){dir.create(paste(WKdir,sep="","\\"))}
  
  SMPLIST<-as.data.frame(SMPLIST)
  write.csv(SMPLIST,paste(WKdir,sep="","\\","SMPqueryRESULT.csv"),row.names=FALSE)
  SMPRESULT<-as.data.frame(table(SMPLIST$SMPDB.ID))
  write.csv(SMPRESULT,paste(WKdir,sep="","\\","SMPRESULT.csv"),row.names=FALSE)
  SMPRESULT<-read.csv(file=paste(WKdir,sep="","\\","SMPRESULT.csv"),col.names= c("SMPDB.ID","Freq"))
  SMPRESULTfreq<-NULL
  
  for (i in 1: length(SMPRESULT$SMPDB.ID)){
    SMPRESULTfreq<-c(SMPRESULTfreq,SMPDBFREQ[SMPDBFREQ[,"SMPDB.ID"]==as.character(SMPRESULT[i,"SMPDB.ID"]),"Freq"])
    SMPRESULT$Name[i]<-as.character(Map_SMP_MP[Map_SMP_MP[,"SMPDB.ID"] %in% SMPRESULT[,"SMPDB.ID"][i],"Name"])
  }
  SMPRESULT$DBfreq<-SMPRESULTfreq
  SMPRESULT$pScore<-SMPRESULT$Freq/SMPRESULT$DBfreq
  write.csv(SMPRESULT,paste(WKdir,sep="","\\","SMPRESULT.csv"),row.names=FALSE)
  return(SMPLIST)
}


###Ion Image manipulation
META_pathway_ion_image<- function(SMPLIST, SMPRESULT, datafile,MALDI_IMAGE, Image_Type="",PWpScore=0.26,plotTolerance=0.125 ,creat_new_file=TRUE,...) {
  if (dir.exists(datafile)==FALSE){dir.create(datafile)}
  for (j in 1:length(SMPRESULT$SMPDB.ID)){
    if (SMPRESULT$pScore[j]>=PWpScore){
      WKdir<-paste(datafile,"\\",SMPRESULT$SMPDB.ID[j],sep="")
            PWSUBLIST<-SMPLIST[SMPLIST$SMPDB.ID==SMPRESULT$SMPDB.ID[j],]
      if (dir.exists(WKdir)==FALSE){dir.create(WKdir)}
      for (i in 1:length(PWSUBLIST$mz)){
        Plot_Ion_image_Png(WKdir,MALDI_IMAGE,PWSUBLIST$moleculeNames[i],PWSUBLIST$mz[i],PWSUBLIST$adduct[i],Tolerance= plotTolerance, title=Image_Type,Creat_new_file=creat_new_file)
      }}
  }}

META_pathway_ion_image_cardinal<- function(SMPLIST, SMPRESULT, datafile,MALDI_IMAGE, Image_Type="",PWpScore=0.26,plotTolerance=0.125 ,creat_new_file=TRUE,...) {
  if (dir.exists(datafile)==FALSE){dir.create(datafile)}
  for (j in 1:length(SMPRESULT$SMPDB.ID)){
    if (SMPRESULT$pScore[j]>=PWpScore){
      WKdir<-paste(datafile,"\\",SMPRESULT$SMPDB.ID[j],sep="")
      PWSUBLIST<-SMPLIST[SMPLIST$SMPDB.ID==SMPRESULT$SMPDB.ID[j],]
      if (dir.exists(WKdir)==FALSE){dir.create(WKdir)}
      for (i in 1:length(PWSUBLIST$mz)){
        Plot_Ion_image_Png_Cardinal(WKdir,MALDI_IMAGE,PWSUBLIST$moleculeNames[i],PWSUBLIST$mz[i],PWSUBLIST$adduct[i],Tolerance= plotTolerance, title=Image_Type,Creat_new_file=creat_new_file,contrast.enhance = "suppression")
      }}
  }}
###DRAW META_pathway_Whole_Picture
META_pathway_Whole_Picture<- function(SMPLIST, SMPRESULT, datafile, PWpScore=0.26,pathway_mode="", ...) {
  MPID <-NULL
  if (dir.exists(datafile)==FALSE){dir.create(datafile)}
  for (j in 1:length(SMPRESULT$SMPDB.ID)){
    if (SMPRESULT$pScore[j]>=PWpScore){
      WKdir<-paste(datafile,"\\",SMPRESULT$SMPDB.ID[j],sep="")
      if (exists("Map_SMP_MP")==FALSE){Map_SMP_MP <- read.csv("DB/smpdb_pathways.csv")}
      MPID <- sprintf("PW%06d", match(SMPRESULT$SMPDB.ID[j],Map_SMP_MP[,1]))
      #PWsvg <- htmlParse(paste("DB/smpdb_svg/",MPID ,".svg",sep=""))
      PWSUBLIST<-SMPLIST[SMPLIST$SMPDB.ID==SMPRESULT$SMPDB.ID[j],c("HMDB.ID","mz","moleculeNames","Metabolite.ID")]
      PW_svg_line<- readLines(paste("DB/smpdb_svg/",MPID ,".svg",sep=""))
      for (i in 1:length(PWSUBLIST$mz)){
        pngfile<- NULL
        res <- try(pngfile<-image_read(paste(WKdir,"\\",PWSUBLIST$moleculeNames[i],'.png',sep="")),silent = TRUE)
        if (class(res) != "try-error"){
          Ion_image_svg_base64 <- base64Encode(readBin(paste(WKdir,"\\",PWSUBLIST$moleculeNames[i],'.png',sep=""), "raw", file.info(paste(WKdir,"\\",PWSUBLIST$moleculeNames[i],'.png',sep=""))[1, "size"]), "txt")
          res<-try(Start_moldb_cache<-grep(paste('id=','"',PWSUBLIST$Metabolite.ID[i],'_moldb_cache',sep=""),as.character(PW_svg_line),value=FALSE),silent = TRUE)
          res2<-try(End_moldb_cache<-grep("pointer-events=",as.character(PW_svg_line[Start_moldb_cache:length(PW_svg_line)]),value=FALSE)[1],silent = TRUE)
          if ("&"(class(res) != "try-error",class(res2) != "try-error")){
            
            PW_svg_line<-c(PW_svg_line[1:Start_moldb_cache], paste("<image xlink:href=\"data:image/png;base64,",Ion_image_svg_base64,aep=""),PW_svg_line[(Start_moldb_cache+End_moldb_cache-1):length(PW_svg_line)])
            if (pathway_mode=="D-mode"){write(PW_svg_line,file=paste(datafile,"\\",SMPRESULT$SMPDB.ID[j]," D-mode.svg",sep=""))}else
            {write(PW_svg_line,file=paste(datafile,"\\",SMPRESULT$SMPDB.ID[j],".svg",sep=""))}}
        } 
          
          
          
        }
      }
    }
}


#Metapathway Image discovery_mode_list
META_pathway_Whole_Picture_discovery_mode<- function(SMPLIST, SMPRESULT, datafile, SMPDB, adductslist="M-H", PWpScore=0.26, ...) {
  
  MPID <-NULL
  SMPLIST_D_MODE<-SMPDB[`&`(SMPDB[,"SMPDB.ID"] %in% SMPRESULT[SMPRESULT[,"pScore"]>0.26,"SMPDB.ID"],!(SMPDB[,"HMDB.ID"] %in% unique(SMPLIST[,"HMDB.ID"]))),]
  
  SMPLIST_D_MODE$mass<- t(as.data.frame(lapply(SMPLIST_D_MODE$Formula, getMonomass)))
  SMPLIST_D_MODE$moleculeNames<-SMPLIST_D_MODE$Metabolite.Name
  for (adduction in 1:length(adductslist)){
    SMPLIST_D_MODE$adduct<-adductslist[adduction]
    SMPLIST_D_MODE$mz <- SMPLIST_D_MODE$mass+ adducts[adducts[,"Name"]==adductslist[adduction],"Mass"]
    
    if (adduction==1){
      META_pathway_ion_image(SMPLIST_D_MODE, SMPRESULT, datafile,MALDI_IMAGE, Image_Type="D-mode", PWpScore=0.26,plotTolerance=0.125,creat_new_file=TRUE)
    }else{
      META_pathway_ion_image(SMPLIST_D_MODE, SMPRESULT, datafile,MALDI_IMAGE, Image_Type="D-mode", PWpScore=0.26,plotTolerance=0.125,creat_new_file=FALSE)
      } 
  }
  options(warn=-1)
  SMPLIST_D_MODE_combine<-bind_rows(SMPLIST_D_MODE[,colnames(SMPLIST)],SMPLIST)
  options(warn=0)
  META_pathway_Whole_Picture(SMPLIST_D_MODE_combine, SMPRESULT, datafile, PWpScore=0.26,pathway_mode="D-mode")
}

#Cluster_mode_list
Whole_Picture_Cluster_mode<- function(workdir=WorkingDir(), 
                                      cluster_list_file="Cluster.csv",
                                      clustercol="Protein", 
                                      Sample_list_file="SMAPLE.csv", 
                                      adductslist="M-H", 
                                      Threshold=0.15,...
                                      ) {

  cluster_list<-read.csv(cluster_list_file)
  Sample_list<-read.csv(Sample_list)
  for (i in 1:length(Sample_list$Sample)){
    MALDI_IMAGE_list[as.numeric(Sample_list$Order[i])] <- importImzMl(paste(file.path(Sample_list$Sample[i]),".imzml",sep=""))
  }
  
  
  MPID <-NULL
  SMPLIST_D_MODE<-SMPDB[`&`(SMPDB[,"SMPDB.ID"] %in% SMPRESULT[SMPRESULT[,"pScore"]>0.26,"SMPDB.ID"],!(SMPDB[,"HMDB.ID"] %in% unique(SMPLIST[,"HMDB.ID"]))),]
  
  SMPLIST_D_MODE$mass<- t(as.data.frame(lapply(SMPLIST_D_MODE$Formula, getMonomass)))
  SMPLIST_D_MODE$moleculeNames<-SMPLIST_D_MODE$Metabolite.Name
  for (adduction in 1:length(adductslist)){
    SMPLIST_D_MODE$adduct<-adductslist[adduction]
    SMPLIST_D_MODE$mz <- SMPLIST_D_MODE$mass+ adducts[adducts[,"Name"]==adductslist[adduction],"Mass"]
    
    if (adduction==1){
      META_pathway_ion_image(SMPLIST_D_MODE, SMPRESULT, datafile,MALDI_IMAGE, Image_Type="D-mode", PWpScore=0.26,plotTolerance=0.125,creat_new_file=TRUE)
    }else{
      META_pathway_ion_image(SMPLIST_D_MODE, SMPRESULT, datafile,MALDI_IMAGE, Image_Type="D-mode", PWpScore=0.26,plotTolerance=0.125,creat_new_file=FALSE)
    } 
  }
  options(warn=-1)
  SMPLIST_D_MODE_combine<-bind_rows(SMPLIST_D_MODE[,colnames(SMPLIST)],SMPLIST)
  options(warn=0)
  META_pathway_Whole_Picture(SMPLIST_D_MODE_combine, SMPRESULT, datafile, PWpScore=0.26,pathway_mode="D-mode")
}

#for (j in 1:length(SMPRESULT$SMPDB.ID)){

#if (SMPRESULT$pScore[j]>=PWpScore){
#dir<-paste(datafile,"\\",SMPRESULT$SMPDB.ID[j],sep="")
#discovery_moldb_list<-NULL
#if (exists("Map_SMP_MP")==FALSE){Map_SMP_MP <- read.csv("DB/smpdb_pathways.csv")}

#MPID <- sprintf("PW%06d", match(SMPRESULT$SMPDB.ID[j],Map_SMP_MP[,1]))
#PWsvg <- htmlParse(paste("DB/smpdb_svg/",MPID ,".svg",sep=""))
#PWSUBLIST<-SMPLIST[SMPLIST$SMPDB.ID==SMPRESULT$SMPDB.ID[j],c("HMDB.ID","mz","moleculeNames","Metabolite.ID")]
#PW_svg_line<- readLines(paste("DB/smpdb_svg/",MPID ,".svg",sep=""))

#discovery_moldb_cache<-grep('<g id="location" transform=',as.character(PW_svg_line),value=FALSE)
#discovery_moldb_cache<-str_locate(PW_svg_line[discovery_moldb_cache],"data-element-id=")

#for (i in 1 : length(discovery_moldb_cache)){
#temp<-substr(PW_svg_line[discovery_moldb_cache[i]],discovery_moldb_startpos[i,"end"]+2,discovery_moldb_startpos[i,"end"]+2+9)
#discovery_moldb_cache<-grep('data-element-type="compound" data-element-id=',as.character(PW_svg_line),value=FALSE)
#discovery_moldb_startpos<-str_locate(PW_svg_line[discovery_moldb_cache],"data-element-id=")
#discovery_moldb_list<-substr(PW_svg_line[discovery_moldb_cache],discovery_moldb_startpos[,"end"]+2,discovery_moldb_startpos[,"end"]+2+9)
#discovery_moldb_list <- union(discovery_moldb_list,discovery_moldb_list)
#discovery_moldb_list<- fetch_ID_from_DB(discovery_moldb_list,SMPDB,SMPRESULT$SMPDB.ID[j])



##fetch_ID_from SMPDB
fetch_ID_from_DB <- function(Fetch_id_list, SMPDB, SMPDB.ID){
  if (SMPDB.ID != "") {for (z in 1: length(Fetch_id_list)){temp3[z,]<-subset(SMPDB , Metabolite.ID == as.character(Fetch_id_list)[z] & SMPDB.ID ==SMPRESULT$SMPDB.ID[j])}}else{
    for (z in 1: length(Fetch_id_list)){temp3[z,]<-subset(SMPDB , Metabolite.ID == as.character(Fetch_id_list)[z] )}}
  return(temp3)
}



nice_file_create <- function(filename){
  ok <- file.create(filename)
    if(!ok)
    {return(windows_filename2(filename))
    }else{return(filename)}
  
}

Kernel_convolution<-function(Neighbour,resolution=1){
  Neighbour<-1+(Neighbour*2*resolution)
  if (Neighbour==1){kern = matrix(1, ncol = Neighbour, nrow = Neighbour)}else{
  kern = matrix(1/(Neighbour*Neighbour-1), ncol = Neighbour, nrow = Neighbour)  
  kern[(Neighbour+1)/2,(Neighbour+1)/2]=0}
  kern
}
  
  
  
  
## Basic Function PNG ion image plot
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
    
    #resolution<-imagefile[[1]]@metaData[["imaging"]][["Size"]]["x"]
    kern=Kernel_convolution(Neighbour,resolution)
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

Plot_Ion_image_Png3<- function(WKdir, imagefile, png_filename, mz,adduct="M-H", Tolerance= 0.25, title="",Neighbour = 0,Creat_new_file=TRUE,color="black",interpolate =FALSE){
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
    a<-image_attributes(pngfile)
    pixelx<-as.numeric(str_split(a[a$property=="png:IHDR.width,height","value"],", ")[[1]][1])
    pixely<-as.numeric(str_split(a[a$property=="png:IHDR.width,height","value"],", ")[[1]][2])
    try(resolutionx<-pixelx/imagefile[[1]]@metaData[["imaging"]][["size"]]["x"])
    try(resolutiony<-pixely/imagefile[[1]]@metaData[["imaging"]][["size"]]["y"])
    try(resolution<-sqrt(as.numeric((resolutionx + resolutiony) / 2)))
    try(kern<-Kernel_convolution(Neighbour,resolution))
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

Plot_Ion_image_Png2<- function(WKdir, imagefile, png_filename, mz,adduct="M-H", Tolerance= 0.25, title="",Creat_new_file=TRUE){
  #x11()
  pngfillewrite<-paste(WKdir,"\\",png_filename,'.png',sep="")
  #pngfillewrite<-nice_file_create(pngfillewrite)
  #dev.copy(png,paste(dir,"\\",png_filename,'.png',sep=""))
  png(paste(WKdir,"\\","temp.png",sep=""))
  try(plotMsiSlice(imagefile ,mz , tolerance=Tolerance,colRamp=colorRamp(c("black", "blue", "green", "yellow", "red")), legend=FALSE),silent = TRUE)
  dev.off()
  try(pngfile<-image_read(paste(WKdir,"\\","temp.png",sep="")))
  res <- try(pngfile<-image_trim(pngfile),silent = TRUE)
  if (class(res) != "try-error"){
    if (title=="D-mode"){
      pngfile<-image_annotate(pngfile,paste(png_filename,"D-mode","\n",adduct,sprintf("%g",mz),"±",Tolerance),size = 9)
      res2 <- try(pngfileoriginal<-image_read(pngfillewrite),silent = TRUE)
      if ('&'(class(res2) != "try-error",Creat_new_file==FALSE)){
        pngfile<-image_append(c(pngfileoriginal,pngfile))
        image_write(pngfile,pngfillewrite)}
      if ('&'(class(res2) != "try-error",Creat_new_file==TRUE)){    
        image_write(pngfile,pngfillewrite)}
      if (class(res2) == "try-error"){    
        image_write(pngfile,pngfillewrite)}    
    }else{
      pngfile<-image_annotate(pngfile,paste(png_filename,"\n", adduct,sprintf("%g",mz),"±",Tolerance),size = 9)
      image_write(pngfile,pngfillewrite)
    }
  }
  file.remove(paste(WKdir,"\\","temp.png",sep=""))
}

#wrap the getMolecule function


##prepare SMPDB

#prepare_SMPDB<- function(DBpath="DB"){
#SMPDB <- list.files(path=paste(DBpath,"/metapath",sep=""),full.names = TRUE) %>%  lapply(read.csv) %>%  bind_rows 
#write.csv(SMPDB, file=paste(DBpath,"/SMPdb.csv",sep=","))
#if (exists("SMPDB")==FALSE){SMPDB <- read.csv(file="DB/SMPdb.csv", as.is = TRUE)}
#SMPDBFREQ<-as.data.frame(table(SMPDB[,"SMPDB.ID"]))
#write.csv(SMPDBFREQ, file=(paste(DBpath,"/SMPdb_ID_FREQ.csv",sep="")))
#DBpath="DB"
#SMPDB <- read.csv(file=paste(DBpath,"/SMPdb.csv",sep=""), as.is = TRUE)
#SMPDBFREQ <- read.csv(file=paste(DBpath,"/SMPdb_ID_FREQ.csv",sep="") ,col.names= c("","SMPDB.ID","Freq"),as.is = TRUE)
#pathwayDB <- read.csv(file=paste(DBpath,"/smpdb_pathways.csv",sep=""), as.is = TRUE)
#R Python interface
#install.packages("reticulate")


#install.packages("R.matlab")
#}
quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

Cardinal_old<- function(){
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
  p_load(org.Hs.eg.db)
  p_load(Rgraphviz)
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
  plot(g)
}


Cardinal_cluster_classification<-function(){
  set.seed(1)
  skm <- spatialKMeans(imdata, r=2, k=4, method="adaptive")
  plot(skm, col=c("pink", "blue", "red","orange","navyblue"), type=c('p','h'), key=FALSE)
  image(skm, col=c("pink", "blue", "red","orange","navyblue"), key=FALSE)
  set.seed(1)
  ssc <- spatialShrunkenCentroids(imdata, r=2, k=5, s=3, method="gaussian")
  plot(ssc, col=c("pink", "blue", "red","orange","navyblue"), type=c('p','h'), key=FALSE)
  image(ssc, col=c("pink", "blue", "red","orange","navyblue"), key=FALSE)
  
  
  
  
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
