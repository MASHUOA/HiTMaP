library(magrittr)
library(shiny)
library(shinyFiles)
library(future)
plan(multisession)
options(shiny.maxRequestSize = 3000*1024^2)

list.dirs.depth.n <- function(p, n) {
    res <- list.dirs(p, recursive = FALSE)
    if (n > 1) {
        add <- list.dirs.depth.n(res, n-1)
        c(res, add)
    } else {
        res
    }
}

current_states_filesum<-function(wd=getwd()){
    library(stringr)
    folder_file_DT=dir(wd,recursive = F)
    folder_file_DT_folder<-list.dirs.depth.n(wd,1)
    folder_file_DT_folder<-gsub(wd,"",folder_file_DT_folder)
    folder_file_DT_folder<-gsub("^/","",folder_file_DT_folder)
    folder_file_DT_folder_ID<-folder_file_DT_folder[stringr::str_detect(folder_file_DT_folder,regex(' ID$', ignore_case = T))]
    folder_file_DT_folder_Summary<-folder_file_DT_folder[stringr::str_detect(folder_file_DT_folder,regex('Summary folder$', ignore_case = T))]
    folder_file_DT_folder_Source<-folder_file_DT_folder[stringr::str_detect(folder_file_DT_folder,regex('Source$', ignore_case = T))]
    
    folder_file_DT<-folder_file_DT[!str_detect(folder_file_DT,"source/")]
    #paste0(unlist(input$Projectdir[[2]]),paste0(input$Projectdir[[1]],collapse = '/',sep=""),collapse = '/',sep="")
    folder_file_DT_size<-file.size(folder_file_DT)
    folder_file_DT_imzml<-folder_file_DT[stringr::str_detect(folder_file_DT,regex('.imzml$', ignore_case = T))]
    folder_file_DT_ibd<-folder_file_DT[stringr::str_detect(folder_file_DT,regex('.ibd$', ignore_case = T))]
    
    
    folder_file_DT_csv<-folder_file_DT[stringr::str_detect(folder_file_DT,regex('.csv$|.txt$|.xls|.xlsx', ignore_case = T))]
    folder_file_DT_rda<-folder_file_DT[stringr::str_detect(folder_file_DT,regex('.RDA$|.RDS$|.Rdata$', ignore_case = T))]
    folder_file_DT_fasta<-folder_file_DT[stringr::str_detect(folder_file_DT,regex('.fasta$', ignore_case = T))]
    
    folder_file_DT_imzml_valid<-intersect(tools::file_path_sans_ext(folder_file_DT_imzml),tools::file_path_sans_ext(folder_file_DT_ibd))
    folder_file_DT_imzml<-folder_file_DT_imzml[tools::file_path_sans_ext(folder_file_DT_imzml) %in% folder_file_DT_imzml_valid]
    folder_file_DT_ibd<-folder_file_DT_ibd[tools::file_path_sans_ext(folder_file_DT_ibd) %in% folder_file_DT_imzml_valid]
    
    folder_file_DT_size_imzml<-file.size((paste0(wd,"/",folder_file_DT_imzml)))
    folder_file_DT_size_ibd<-file.size((paste0(wd,"/",folder_file_DT_ibd)))
    folder_file_DT_size_fasta<-file.size((paste0(wd,"/",folder_file_DT_fasta)))
    folder_file_DT_size_csv<-file.size((paste0(wd,"/",folder_file_DT_csv)))
    folder_file_DT_size_rda<-file.size((paste0(wd,"/",folder_file_DT_rda)))
    
    
    resultdf<- data.frame(filename=c(folder_file_DT_imzml_valid,
                                     folder_file_DT_csv,
                                     folder_file_DT_fasta,
                                     folder_file_DT_rda,
                                     folder_file_DT_folder_ID,
                                     folder_file_DT_folder_Summary,
                                     folder_file_DT_folder_Source),
                          filetype=c(rep("IMS data",length(folder_file_DT_imzml_valid)),
                                     rep("Config file",length(folder_file_DT_csv)),
                                     rep("Database",length(folder_file_DT_fasta)),
                                     rep("RDA",length(folder_file_DT_rda)),
                                     rep("ID folder",length(folder_file_DT_folder_ID)),
                                     rep("Summary",length(folder_file_DT_folder_Summary)),
                                     rep("Summary",length(folder_file_DT_folder_Source))
                                     ))                                     
    
    
    if (nrow(resultdf)!=0){
        resultdf$size=c(if(length(folder_file_DT_imzml_valid)) (utils:::format.object_size(folder_file_DT_size_imzml+folder_file_DT_size_ibd, standard = "SI",units = "MB")),
                        if(length(folder_file_DT_csv)) utils:::format.object_size(folder_file_DT_size_csv, standard = "SI",units = "MB"),
                        if(length(folder_file_DT_fasta)) utils:::format.object_size(folder_file_DT_size_fasta, standard = "SI",units = "MB"),
                        if(length(folder_file_DT_rda)) utils:::format.object_size(folder_file_DT_size_rda, standard = "SI",units = "MB"),
                        if(length(folder_file_DT_folder_ID)) rep("",length(folder_file_DT_folder_ID)),
                        if(length(folder_file_DT_folder_Summary)) rep("",length(folder_file_DT_folder_Summary)),
                        if(length(folder_file_DT_folder_Source)) rep("",length(folder_file_DT_folder_Source))
                        )
    }
    
    return(resultdf)
}
Preprocessing_segmentation<-function(datafile,
                                     workdir=NULL,
                                     segmentation_num=5,threshold=0.001,
                                     ppm,
                                     mzrange=c(500,4000),
                                     Segmentation=c("spatialKMeans","spatialShrunkenCentroids","Virtual_segmentation","none","def_file"),
                                     Segmentation_def="segmentation_def.csv",
                                     Segmentation_ncomp="auto-detect",
                                     Segmentation_variance_coverage=0.8,
                                     Smooth_range=1,
                                     colorstyle="Set1",
                                     Virtual_segmentation_rankfile=NULL,
                                     rotate_file=NULL,
                                     BPPARAM=bpparam(),
                                     preprocess=list(force_preprocess=FALSE,use_preprocessRDS=TRUE,smoothSignal=list(method="gaussian"),
                                                     reduceBaseline=list(method="locmin"),
                                                     peakPick=list(method="adaptive"),
                                                     peakAlign=list(tolerance=ppm, units="ppm")),
                                     ...){
    
    suppressMessages(suppressWarnings(require(data.table)))
    suppressMessages(suppressWarnings(require(Cardinal)))
    suppressMessages(suppressWarnings(require(RColorBrewer)))
    suppressMessages(suppressWarnings(require(stringr)))
    setCardinalBPPARAM(BPPARAM)
    if ((!is.null(rotate_file))){
        message("Found rotation info")
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
    if (is.null(workdir)){
        workdir<-base::dirname(datafile[1])
        datafile<-basename(datafile)
    }else{ workdir<-workdir }
    
    
    datafile_imzML<-datafile
    for (z in 1:length(datafile)){
        Peptide_Summary_file$Intensity<-rep(0,nrow(Peptide_Summary_file))
        name <-basename(datafile[z])
        name <-gsub(".imzML$","",name)
        name <-gsub("/$","",name)
        folder<-base::dirname(datafile[z])
        #imdata <- Cardinal::readImzML(datafile[z],preprocessing = F,attach.only = T,resolution = 200,rotate = rotate[z],as="MSImageSet",BPPARAM = BPPARAM)
        if (!str_detect(datafile[z],".imzML$")){
            datafile_imzML[z]<-paste0(datafile[z],".imzML")
        }
        setwd(workdir)
        message("Porject dir",workdir)
        # imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,resolution=200, units="ppm",BPPARAM=BPPARAM,mass.range =mzrange)
        # if(!is.na(rotate[datafile_imzML[z]])){
        #   imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
        # }else if(!is.na(rotate[datafile[z]])){
        #   imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
        # }
        # 
        if (ppm>=25) {
            instrument_ppm=50
        }else{
            instrument_ppm=10
        }      
        
        
        if (dir.exists(paste0(gsub(".imzML$","",datafile[z]) ," ID"))==FALSE){
            dir.create(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
        }else{
            tryCatch(unlink(paste0(gsub(".imzML$","",datafile[z]) ," ID"),recursive = T))
            dir.create(paste0(gsub(".imzML$","",datafile[z])  ," ID"))
        }
        
        if (!file.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))){
            message("Preparing image data for statistical analysis: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
            imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,resolution=ppm, units="ppm",BPPARAM=BPPARAM,mass.range =mzrange)
            if(!is.na(rotate[datafile_imzML[z]])){
                imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile_imzML[z]])
            }else if(!is.na(rotate[datafile[z]])){
                imdata <-rotateMSI(imdata=imdata,rotation_degree=rotate[datafile[z]])
            }
            if (!is.null(preprocess)){
                if  ( ppm<25){
                    imdata_ed<-imdata %>% 
                        #smoothSignal(method="gaussian") %>% 
                        #reduceBaseline(method="locmin") %>%
                        peakPick(method=peakPick$method) %>%
                        peakAlign(tolerance=ppm, units="ppm") %>%
                        #peakFilter(mse_pre, freq.min=0.00) %>%
                        process()
                } else if(ppm>=25){
                    imdata_ed<-imdata %>% 
                        smoothSignal(method=smoothSignal$method) %>% 
                        reduceBaseline(method=reduceBaseline$method) %>%
                        peakPick(method=peakPick$method) %>%
                        peakAlign(tolerance=ppm, units="ppm") %>%
                        #peakFilter(mse_pre, freq.min=0.00) %>%
                        process()
                }
                saveRDS(imdata_ed,paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"),compress = F)  
            }}else{
                imdata_ed<-readRDS(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
            }
        
        if ('&'(preprocess$use_preprocessRDS,file.exists(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS")))){
            message("Using image data: ",paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
            
            imdata_ed<-readRDS(paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_imdata.RDS"))
            #imdata_ed<-imdata
            imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,resolution=ppm, units="ppm",BPPARAM=BPPARAM,mass.range =mzrange)
            
        }else{
            message("Using image data: ",paste0(gsub(".imzML$","",datafile[z]), ".imzML"))
            imdata <- Cardinal::readMSIData(datafile_imzML[z],  attach.only=T,resolution=ppm, units="ppm",BPPARAM=BPPARAM,mass.range =mzrange)
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
        
        
        if (Segmentation==F){
            message("Segmentation in progress...")
            #cl=autoStopCluster(makeCluster(6))
            
            
            if (Segmentation[1]=="spatialKMeans" && segmentation_num!=1) {
                if ('&'(Segmentation_ncomp=="auto-detect",Segmentation_variance_coverage>0)){
                    Segmentation_ncomp<-PCA_ncomp_selection(imdata=imdata,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
                }
                set.seed(1)
                
                skm <-  suppressMessages(suppressWarnings(spatialKMeans(imdata, r=Smooth_range, k=segmentation_num, method="adaptive",ncomp=Segmentation_ncomp,BPPARAM =BPPARAM )))
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
                write.csv(correlation,paste(Segmentation[1],"_RESULT","correlation",segmentation_num,"segs.csv"),row.names = F)
                write.csv(centers,paste(Segmentation[1],"_RESULT","centers",segmentation_num,"segs.csv"),row.names = F)
                write.csv(cluster,paste(Segmentation[1],"_RESULT","cluster",segmentation_num,"segs.csv"),row.names = F)
                
                
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
                    Segmentation_ncomp<-PCA_ncomp_selection(imdata=imdata,variance_coverage=Segmentation_variance_coverage,outputdir=paste0(getwd(),"/"))
                }
                #message(paste0("spatialShrunkenCentroids computing for ",name))
                skm <-  suppressMessages(suppressWarnings(spatialShrunkenCentroids(imdata, r=Smooth_range, k=segmentation_num, method="adaptive",ncomp=Segmentation_ncomp,BPPARAM =BPPARAM)))
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
                write.csv(correlation,paste(Segmentation[1],"_RESULT","correlation",segmentation_num,"segs.csv"),row.names = F)
                write.csv(centers,paste(Segmentation[1],"_RESULT","centers",segmentation_num,"segs.csv"),row.names = F)
                write.csv(cluster,paste(Segmentation[1],"_RESULT","cluster",segmentation_num,"segs.csv"),row.names = F)
                
                
                
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
                
                
                
                
                
                
            }else if
            (Segmentation[1]=="Virtual_segmentation"){
                radius_rank=read.csv(file = Virtual_segmentation_rankfile)
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
                #                        resolution=100, step=3.3, as="MSImageSet")
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
                
                
                
            }else if(Segmentation=="def_file"){
                
                Segmentation_def_tbl<-read.csv(paste0(workdir,"/", Segmentation_def))
                
                Segmentation_def_tbl<-Segmentation_def_tbl[Segmentation_def_tbl$datafile==datafile[z],]
                
                Segmentation_def_tbl<-Segmentation_def_tbl[order(Segmentation_def_tbl$pixel),]
                
                Segmentation_def_tbl$label<-as.factor(Segmentation_def_tbl$label)
                
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
                
            }else{
                
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
}

dput_df <- function(df,
                        name=as.character(substitute(df)),
                        sep='\t',
                        header=TRUE,
                        stringsAsFactors = FALSE,
                        n= nrow(df),
                        random=FALSE,
                        seed = 1){
    name
    if(random) {
        set.seed(seed)
        df <- df[sample(1:nrow(df),n),]
    } else {
        df <- head(df,n)
    }
    cat(sep='',name,' <- read.table(sep="',sub('\t','\\\\t',sep),'", text="\n  ',
        paste(colnames(df),collapse=sep))
    df <- head(df,n)
    apply(df,1,function(x) cat(sep='','\n  ',paste(x,collapse=sep)))
    cat(sep='','", header=',header,', stringsAsFactors=',stringsAsFactors,')')
    
    sapply(names(df), function(x){
        if(is.character(df[[x]]) & suppressWarnings(identical(as.character(as.numeric(df[[x]])),df[[x]]))){ # if it's a character column containing numbers
            cat(sep='','\n',name,'$',x,' <- as.character(', name,'$',x,')')
        } else if(is.factor(df[[x]]) & !stringsAsFactors) { # if it's a factor and conversion is not automated
            cat(sep='','\n',name,'$',x,' <- factor(', name,'$',x,')')
        } else if(inherits(df[[x]], "POSIXct")){
            cat(sep='','\n',name,'$',x,' <- as.POSIXct(', name,'$',x,')')
        } else if(inherits(df[[x]], "Date")){
            cat(sep='','\n',name,'$',x,' <- as.Date(', name,'$',x,')')
        }})
    invisible(NULL)
}
   
dput_list<- function(x,
                       name=as.character(substitute(x)),
                       multiline = TRUE,
                       n=if ('list' %in% class(x)) length(x) else nrow(x),
                       random=FALSE,
                       seed = 1){
    name
    if('tbl_df' %in% class(x)) create_fun <- "tibble::tibble" else
        if('list' %in% class(x)) create_fun <- "list" else
            if('data.table' %in% class(x)) create_fun <- "data.table::data.table" else
                create_fun <- "data.frame"
            
            if(random) {
                set.seed(seed)
                if(create_fun == "list") x <- x[sample(1:length(x),n)] else 
                    x <- x[sample(1:nrow(x),n),]
            } else {
                x <- head(x,n)
            }
            
            line_sep <- if (multiline) "\n    " else ""
            paste(sep='',name," <- ",create_fun,"(\n  ",
                paste0(unlist(
                    Map(function(item,nm) paste0(nm,if(nm=="") "" else " = ",paste(capture.output(dput(item)),collapse=line_sep)),
                        x,if(is.null(names(x))) rep("",length(x)) else names(x))),
                    collapse=",\n  "),
                if(create_fun == "data.frame") ",\n  stringsAsFactors = FALSE)" else "\n)")
} 

dput_tibble <- function(df,
                        name=as.character(substitute(df)),
                        n= nrow(df),
                        random=FALSE,
                        seed = 1){
    name
    if(random) {
        set.seed(seed)
        df <- df[sample(1:nrow(df),n),]
    } else {
        df <- head(df,n)
    }
    df1 <- lapply(df,function(col) if(is.factor(col)) as.character(col) else col)
    dputs   <- sapply(df1,function(col){
        col_dputs <- sapply(col,function(elt) paste(capture.output(dput(elt)),collapse=""))
        max_char <- max(nchar(unlist(col_dputs)))
        sapply(col_dputs,function(elt) paste(c(rep(" ",max_char-nchar(elt)),elt),collapse=""))
    })
    lines   <- paste(apply(dputs,1,paste,collapse=", "),collapse=",\n  ")
    output  <- paste0(name," <- tibble::tribble(\n  ",
                      paste0("~",names(df),collapse=", "),
                      ",\n  ",lines,"\n)")
    cat(output)
    sapply(names(df), function(x) if(is.factor(df[[x]])) cat(sep='','\n',name,'$',x,' <- factor(', name,'$',x,')'))
    invisible(NULL)
}
