Load_Cardinal_imaging<-function(datafile=tk_choose.files(filter = Filters,
                                                         caption  = "Choose single or multiple file(s) for analysis",
                                                         multi = F),
                                BPPARAM = bpparam(),
                                resolution=2.5,
                                memory=TRUE,
                                preprocessing=F,
                                rotate=0,
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
    imdata <-  suppressMessages(suppressWarnings(Cardinal::readImzML(name, folder,resolution=resolution, units="ppm",memory=memory,BPPARAM=BPPARAM,mass.range=mzrange,as="MSImagingExperiment")))
    
  }else  {
    imdata <-  suppressMessages(suppressWarnings(Cardinal::readImzML(name, folder,resolution=resolution, units="ppm",memory=memory,BPPARAM=BPPARAM,mass.range=mzrange,as="MSImagingExperiment")))
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




setworkdir<-function(workdir){
  if (dir.exists(workdir)==FALSE){dir.create(workdir)}
  setwd(workdir)
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

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

seq.ppm <- function(from, to, ppm) {
  length.out <- (log(to) - log(from)) / log((1 + 1e-6 * ppm) / (1 - 1e-6 *ppm))
  length.out <- floor(1 + length.out)
  i <- seq_len(length.out)
  from * ((1 + 1e-6 * ppm) / (1 - 1e-6 * ppm))^(i-1)
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


Parallel.OS<-function(Thread=1,bpprogressbar_t=TRUE,override_type=NULL,bpexportglobals=F,bpforceGC=F){
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
    bpexportglobals(BPPARAM)=bpexportglobals
    BPPARAM$force.GC<-bpforceGC
  }else{
    BPPARAM=BiocParallel::SerialParam() 
  } 
  
  bpprogressbar(BPPARAM)=bpprogressbar_t
  return(BPPARAM)
}

WorkingDir <-function(Mainfolder = tk_choose.dir(caption = "Select working directory")){
  Mainfolder <- Mainfolder
  setwd(Mainfolder)
  Mainfolder
}


Kernel_convolution<-function(Neighbour){
  Neighbour<-1+(Neighbour*2)
  if (Neighbour==1){kern = matrix(1, ncol = Neighbour, nrow = Neighbour)}else{
    kern = matrix(1/(Neighbour*Neighbour-1), ncol = Neighbour, nrow = Neighbour)
    kern[(Neighbour+1)/2,(Neighbour+1)/2]=0}
  kern
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

#initiallization

Filters <- matrix(c( "imzml file", ".imzML",
                     "Text", ".txt", "All files", "*"),
                  3, 2, byrow = TRUE)



resolve_multi_modes<-function(mm,adjust=0.25){
  
  
  library(dplyr)
  library(tidyr)
  get.modes2 <- function(x,adjust,signifi,from,to) {
    den <- density(x, kernel=c("gaussian"),adjust=adjust,from=from,to=to)
    den.s <- smooth.spline(den$x, den$y, all.knots=TRUE, spar=0.1)
    s.1 <- predict(den.s, den.s$x, deriv=1)
    s.0 <- predict(den.s, den.s$x, deriv=0)
    den.sign <- sign(s.1$y)
    a<-c(1,1+which(diff(den.sign)!=0))
    b<-rle(den.sign)$values
    df<-data.frame(a,b)
    df = df[which(df$b %in% -1),]
    modes<-s.1$x[df$a]
    density<-s.0$y[df$a]
    df2<-data.frame(modes,density)
    df2$sig<-signif(df2$density,signifi)
    df2<-df2[with(df2, order(-sig)), ]
    #print(df2)
    df<-as.data.frame(df2 %>%
                        mutate(m = min_rank(desc(sig)) ) %>% #, count = sum(n)) %>%
                        group_by(m) %>%
                        summarize(a = paste(format(round(modes,2),nsmall=2), collapse = ',')) %>%
                        spread(m, a, sep = ''))
    colnames(df)<-paste0("m",1:length(colnames(df)))
    return(df)
  }
  mmdf<-data.frame(mm=mm)
  library(ggplot2)
  #0.25 defines the number of peaks.
  p<-ggplot(mmdf,aes(mm)) + geom_density(adjust=adjust) + xlim((min(mm)-1),(max(mm)+1) )
  #2 defines ties
  modes<-get.modes2(mm,adjust=adjust,2,min(mm)-1,max(mm)+1)
  return(list(p=p,modes=modes))
}

topN_feature<-function (vec, n, dec = T) 
{
  inx <- order(vec, decreasing = dec)[1:n]
  vec <- rep(F, length = length(vec))
  vec[inx] <- T
  return(vec)
}

Imzml_temp_fix<-function(folder=getwd(),name,
                         representation = c("profile spectrum","centroid spectrum")[1],
                         ibdbinarytype= c("processed","continuous")[1]){

  xmlpath <- normalizePath(file.path(folder, paste(name, ".imzML", 
                                                   sep = "")), mustWork = FALSE)
  
  if (!file.exists(xmlpath)) .stop("expected file ", xmlpath, " does not exist")
  ibdpath <- normalizePath(file.path(folder, paste(name, ".ibd", 
                                                   sep = "")), mustWork = FALSE)
  
  if (!file.exists(ibdpath)) .stop("expected file ", ibdpath, " does not exist")
  
  message("reading imzML file: '", xmlpath, "'")
  message("setting imzML file: '", xmlpath, "'"," to the format as ", ibdbinarytype ," ",representation)
  
  parse <- CardinalIO::parseImzML(xmlpath, ibd = TRUE)
  fileContent <- parse[["fileDescription"]][["fileContent"]]
  if ("IMS:1000030" %in% names(fileContent)) {
    ibdbinarytype <- "continuous"
  }
  else if ("IMS:1000031" %in% names(fileContent)) {
    ibdbinarytype <- "processed"
  }
  else {
    ibdbinarytype <- NULL
  }
  if (!is.null(ibdbinarytype)) 
    .message("detected ibd binary type: ", sQuote(ibdbinarytype))
  if ("MS:1000127" %in% names(fileContent)) {
    representation <- "centroid spectrum"
  }
  else if ("MS:1000128" %in% names(fileContent)) {
    representation <- "profile spectrum"
  }
  else {
    representation <- NULL
  }
  if (!is.null(representation)) 
    .message("detected representation: ", sQuote(representation))
  if (parse.only) 
    return(parse)
  units <- match.arg(units)
  .message("reading ibd file: '", ibdpath, "'")
  object <- .attachIbd(parse, name, ibdbinarytype, representation, 
                       attach.only, mass.range, resolution, units, guess.max)
  .log.collapse("loaded dataset:", capture.output(print(object)))
}
  percent<-function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  