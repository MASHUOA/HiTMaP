memory_profile<-function(){
  #install.packages("ggplot2")
  #install.packages("pryr")
  #install.packages("devtools")
  #devtools::install_github("hadley/lineprof")
  #install.packages("tictoc")
  library(Cardinal)
  library(HiTMaP)
  library(pryr)
  library(stringr)
  library(tictoc)
  
  tic.clearlog() 
  for(x in 1:10) { 
    tic(x) 
    Sys.sleep(.001) 
    toc(log=TRUE,quiet=TRUE) 
    } 
  log.txt<-tic.log(format=TRUE) 
  log.lst<-tic.log(format=FALSE) 
  tic.clearlog()
  
  
  
  wd="~/expdata/"
  wd->workdir
  z=1
  ppm=10
  instrument_ppm=10
  parallel=4
 
  
  BPPARAM=HiTMaP:::Parallel.OS(parallel,bpexportglobals=F,bpforceGC=T, bpprogressbar_t=TRUE)
  
  setCardinalBPPARAM(BPPARAM = BPPARAM)
  
  set.seed(2020)
 # mse <- simulateImage(preset=1, npeaks=10, nruns=2, baseline=1)
  #mse
  
  
  
  
  import_ppm = floor(instrument_ppm*2/5)[1]
  
  datafile=c("MouseBrain_Trypsin_FT/Mouse_brain.imzML")
  preprocess=list(force_preprocess=T,
                  use_preprocessRDS=TRUE,
                  smoothSignal=list(method="Disable"),
                  reduceBaseline=list(method="Disable"),
                  peakPick=list(method="adaptive"),
                  peakAlign=list(tolerance=5, units="ppm"),
                  normalize=list(method=c("rms","tic","reference")[1],mz=1))
  
  
  
  
  datafile <- paste0(workdir,"/",datafile)
  workdir <- dirname(datafile)
  datafile <- basename(datafile)
  datafile <-str_remove(datafile,regex(".imzML$"))
  datafile_imzML<-paste0(datafile,".imzML")
  setwd(workdir[z])
  
  
  
  log.txt<-NULL
  
  
  imdata <- Cardinal::readMSIData(datafile_imzML[z],   as="MSImagingExperiment",resolution=10, units="ppm",BPPARAM=SerialParam())
  imdata@centroided<-F
  imdata->imdata_org
  imdata[,2800:3000]->imdata

  tic.clearlog()
  for (resolution in c(10, 200,400,800)){
  resolution(imdata) <- c(ppm=resolution)
  tic(resolution)
  image(imdata)
  toc(log=TRUE,quiet=TRUE) 
  }
  log.txt<-c(log.txt,tic.log(format=TRUE)) 
  
  
  
  preprocess=list(force_preprocess=T,
                  use_preprocessRDS=TRUE,
                  smoothSignal=list(method="gaussian"),
                  reduceBaseline=list(method="median"),
                  peakPick=list(method="adaptive"),
                  peakAlign=list(tolerance=5, units="ppm"),
                  normalize=list(method=c("rms","tic","reference")[1],mz=1))
  imdata_ed<-imdata
  log.txt<-NULL
  imdatalist<-NULL
  tic.clearlog()
  for (resolution in c(800,100,10)){
    imdata_ed<-imdata
    resolution(imdata_ed) <- c(ppm=resolution)
    imdatalist[[paste(resolution,"ORG")]]<-imdata_ed
    tic(paste(resolution,"NR"))
    imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
    imdatalist[[paste(resolution,"NR")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP"))
    #imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
    #imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
    imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_PA"))
    imdata_ed_pa<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_PA")]]<-imdata_ed_pa
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_NR"))
    imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_NR")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    # tic(paste(resolution,"SS_RB_PP_PA_NR"))
    # #imdata_ed_pa<- imdata_ed_pa %>% normalize(method=preprocess$normalize$method) %>% process()
    # imdatalist[[paste(resolution,"SS_RB_PP_PA_NR")]]<-imdata_ed_pa
    # toc(log=TRUE,quiet=TRUE) 
  }
  log.txt<-c(log.txt,tic.log(format=TRUE)) 
  
  imdatalist->imdatalist_param_native
  
  lapply(imdatalist_param_native,object_size)
  
  lapply(imdatalist_param_native,mz)
  
  preprocess=list(force_preprocess=T,
                  use_preprocessRDS=TRUE,
                  smoothSignal=list(method="Disable"),
                  reduceBaseline=list(method="Disable"),
                  peakPick=list(method="adaptive"),
                  peakAlign=list(tolerance=5, units="ppm"),
                  normalize=list(method=c("rms","tic","reference")[1],mz=1))

  imdatalist<-NULL
  tic.clearlog()
  for (resolution in c(800,100,10)){  
    imdata_ed<-imdata
    resolution(imdata_ed) <- c(ppm=resolution)
    imdatalist[[paste(resolution,"ORG")]]<-imdata_ed
    tic(paste(resolution,"NR"))
    imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
    imdatalist[[paste(resolution,"NR")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    
    tic(paste(resolution,"SS_RB_PP"))
    #imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
    #imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
    imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_PA"))
    imdata_ed_pa<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_PA")]]<-imdata_ed_pa
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_NR"))
    imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_NR")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    # tic(paste(resolution,"SS_RB_PP_PA_NR"))
    # #imdata_ed_pa<- imdata_ed_pa %>% normalize(method=preprocess$normalize$method) %>% process()
    # imdatalist[[paste(resolution,"SS_RB_PP_PA_NR")]]<-imdata_ed_pa
    # toc(log=TRUE,quiet=TRUE) 
  }
  log.txt<-c(log.txt,tic.log(format=TRUE)) 
  
  imdatalist->imdatalist_param_default
  
  size1<-lapply(imdatalist_param_native,object_size)
  size2<-lapply(imdatalist_param_default,object_size)
  mz1 <- lapply(imdatalist_param_native,function(x) length(mz(x))) 
  mz2 <- lapply(imdatalist_param_default, function(x) length(mz(x)))
  size<-c(size1,size2)
  
  sizedf<-data.frame(name=names(size1),native_size=unlist(size1),default_size=unlist(size2),native_mz=unlist(mz1),default_mz=unlist(mz2))
  sizedf$resolution=str_split_fixed(sizedf$name," ",2)[,1]
  sizedf$type=str_split_fixed(sizedf$name," ",2)[,2]
  write.csv(sizedf,"sizedf.csv")
  writeLines(unlist(log.txt),"log_subset.txt")
  
  
  
  
  preprocess=list(force_preprocess=T,
                  use_preprocessRDS=TRUE,
                  smoothSignal=list(method="gaussian"),
                  reduceBaseline=list(method="locmin"),
                  peakPick=list(method="adaptive"),
                  peakAlign=list(tolerance=5, units="ppm"),
                  normalize=list(method=c("rms","tic","reference")[1],mz=1))
  imdata_ed<-imdata_org
  log.txt<-NULL
  imdatalist<-NULL
  tic.clearlog()
  for (resolution in c(800,100,10,4)){
    resolution(imdata_ed) <- c(ppm=resolution)
    imdatalist[[paste(resolution,"ORG")]]<-imdata_ed
    tic(paste(resolution,"SS_RB_PP"))
    imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
    imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
    imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method, SNR=preprocess$peakPick$SNR, window=preprocess$peakPick$window) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_PA"))
    imdata_ed_pa<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_PA")]]<-imdata_ed_pa
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_NR"))
    imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_NR")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    # tic(paste(resolution,"SS_RB_PP_PA_NR"))
    # #imdata_ed_pa<- imdata_ed_pa %>% normalize(method=preprocess$normalize$method) %>% process()
    # imdatalist[[paste(resolution,"SS_RB_PP_PA_NR")]]<-imdata_ed_pa
    # toc(log=TRUE,quiet=TRUE) 
  }
  log.txt<-c(log.txt,tic.log(format=TRUE)) 
  
  imdatalist->imdatalist_param_native
  
  lapply(imdatalist_param_native,object_size)
  
  
  
  preprocess=list(force_preprocess=T,
                  use_preprocessRDS=TRUE,
                  smoothSignal=list(method="Disable"),
                  reduceBaseline=list(method="Disable"),
                  peakPick=list(method="adaptive"),
                  peakAlign=list(tolerance=5, units="ppm"),
                  normalize=list(method=c("rms","tic","reference")[1],mz=1))
  imdata_ed<-imdata_org
  imdatalist<-NULL
  tic.clearlog()
  for (resolution in c(800,100,10,4)){
    resolution(imdata_ed) <- c(ppm=resolution)
    imdatalist[[paste(resolution,"ORG")]]<-imdata_ed
    tic(paste(resolution,"SS_RB_PP"))
    #imdata_ed<- imdata_ed %>% smoothSignal(method=preprocess$smoothSignal$method)
    #imdata_ed<- imdata_ed %>% reduceBaseline(method=preprocess$reduceBaseline$method)
    imdata_ed<- imdata_ed %>% peakPick(method=preprocess$peakPick$method) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_PA"))
    imdata_ed_pa<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_PA")]]<-imdata_ed_pa
    toc(log=TRUE,quiet=TRUE) 
    tic(paste(resolution,"SS_RB_PP_NR"))
    imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
    imdatalist[[paste(resolution,"SS_RB_PP_NR")]]<-imdata_ed
    toc(log=TRUE,quiet=TRUE) 
    # tic(paste(resolution,"SS_RB_PP_PA_NR"))
    # #imdata_ed_pa<- imdata_ed_pa %>% normalize(method=preprocess$normalize$method) %>% process()
    # imdatalist[[paste(resolution,"SS_RB_PP_PA_NR")]]<-imdata_ed_pa
    # toc(log=TRUE,quiet=TRUE) 
  }
  log.txt<-c(log.txt,tic.log(format=TRUE)) 
  
  imdatalist->imdatalist_param_default
  
  size1<-lapply(imdatalist_param_native,object_size)
  
  size2<-lapply(imdatalist_param_default,object_size)
  
  size<-c(size1,size2)
  
  sizedf<-data.frame(name=names(size1),native=unlist(size1),default=unlist(size2))
  sizedf$resolution=str_split_fixed(sizedf$name," ",2)[,1]
  sizedf$type=str_split_fixed(sizedf$name," ",2)[,2]
  
  writeLines(unlist(log.txt),"log_all.txt")
  
  
  object_size(imdata_ed)
  
    setCardinalBPPARAM(BPPARAM)
    
    peaklist<-summarizeFeatures(imdata_ed,"sum")
    peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
    peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
    write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec.csv"),row.names = F)
    
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
        peaklist<-summarizeFeatures(imdata_ed,"sum")
        peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
        peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
        peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
        write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
        imdata_ed<-imdata_ed %>% peakBin(peaklist_deco$mz, tolerance=ppm, units="ppm") %>% process()
      }
      
    }
    
    saveRDS(imdata_ed,paste0(gsub(".imzML$","",datafile[z])  ," ID/preprocessed_peakpicked_imdata.RDS"))
    
    if (is.null(preprocess$peakAlign$level)) preprocess$peakAlign$level<-"global"
    if (preprocess$peakAlign$level=="global"){
      if (preprocess$peakAlign$tolerance==0 ) {
        message("preprocess$peakAlign$tolerance set as zero, step bypassed")
      }else if ('&'(!is.null(preprocess$peakAlign$tolerance),!is.null(preprocess$peakAlign$tolerance))){
        message("preprocess$peakAlign$tolerance set as ", preprocess$peakAlign$tolerance)
        imdata_ed<- imdata_ed %>% peakAlign(tolerance=preprocess$peakAlign$tolerance, units=preprocess$peakAlign$units)
        peaklist<-summarizeFeatures(imdata_ed,"sum")
        peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
        peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
        write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec.csv"),row.names = F)
        peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
        write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
      }else {
        message("preprocess$peakAlign$tolerance missing, use default tolerance in ppm ", ppm/2)
        imdata_ed<- imdata_ed %>% peakAlign(tolerance=ppm/2, units="ppm")
        peaklist<-summarizeFeatures(imdata_ed,"sum")
        peaklist_deco<-data.frame(mz=peaklist@mz,intensities=peaklist$sum)
        peaklist_deco<-peaklist_deco[peaklist_deco$intensities>0,]
        write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec.csv"),row.names = F)
        peaklist_deco<-HiTMaP:::isopattern_ppm_filter_peaklist(peaklist_deco,ppm=ppm,threshold=0)
        write.csv(peaklist_deco,paste0(gsub(".imzML$","",datafile[z])  ," ID/Sum_spec_decov.csv"),row.names = F)
      }
    }
    
    imdata_ed<- imdata_ed %>% process()
    
    gc()
    
    if (!is.null(preprocess$normalize)){
      if (preprocess$normalize$method=="Disable") {
      } else if (preprocess$normalize$method %in% c("rms","tic")){
        imdata_ed<- imdata_ed %>% normalize(method=preprocess$normalize$method) %>% process()
      } else if ('&'(preprocess$normalize$method == "reference", !is.null(preprocess$normalize$mz))){
        norm_feature<-which(dplyr::between(imdata_ed@featureData@listData[["mz"]],
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