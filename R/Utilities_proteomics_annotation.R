PlotPMFsig<-function(pimresultindex,spectrumlist,peplist,pimlist,pimresultlist, threshold=0.05){
  suppressMessages(suppressWarnings(require(ggplot2)))
  #library(ggplot)
  print("Plotting Sig PMF")
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
    png(paste(getwd(),"/",windows_filename(str_split(rownamelist,"\\|")[[1]][2]),'.png',sep=""),width = 1980,height = 1080)
    tempsp<-sp+ggtitle(paste(rownamelist,"\nNormalized intensities:",pimresultsl[rownamelist,"Normalized Mean"])) +geom_segment(aes(x = as.numeric(mz), y = as.numeric(yintercept), xend = as.numeric(mzend), yend = as.numeric(intensities), colour =peptide),lineend = "round", data = df,size = 1)
    print(tempsp)
    dev.off()
  }
}

Plot_PMF_all<-function(Protein_feature_list,peaklist,threshold=threshold,savename=""){
  message("Plotting PMF all in one spectrum")
  suppressMessages(suppressWarnings(require(ggplot2)))
  plotpeaklist<-as.data.frame(peaklist)
  plotpeaklist[,1]<-as.numeric(plotpeaklist[,1])
  plotpeaklist[,2]<-as.numeric(plotpeaklist[,2])
  plotpeaklist<-plotpeaklist[plotpeaklist[,2]>0,]
  Protein_feature_list<-Protein_feature_list[Protein_feature_list$Intensity>=max(Protein_feature_list$Intensity)*threshold,]
  Protein_feature_list<-unique(Protein_feature_list[,c( "mz" ,"Intensity","adduct","moleculeNames")])
  sp<-ggplot2::ggplot(plotpeaklist, size=1 ,aes(x=m.z, y=intensities,xend=m.z,yend=rep(0,length(plotpeaklist[,2])))) +  geom_segment() +theme_classic()
  df<- data.frame("mz" = Protein_feature_list$mz,"mzend" = Protein_feature_list$mz, "yintercept" = rep(0,length(Protein_feature_list$mz)), "intensities" = Protein_feature_list$Intensity,"moleculeNames"=as.factor(Protein_feature_list$moleculeNames))
  colnames(df)<-c("mz","mzend","yintercept","intensities","moleculeNames")
  df<-df[df[,"intensities"]>0,]
  png(paste(getwd(),"/",savename," PMF spectrum match",'.png',sep=""),width = 1980,height = 1080)
  tempsp<-sp+ggtitle(paste("PMF spectrum","\nNormalized intensities:")) +geom_segment(aes(x = as.numeric(mz), y = as.numeric(yintercept), xend = as.numeric(mzend), yend = as.numeric(intensities), colour =moleculeNames),lineend = "round", data = df,size = 1,show.legend=F)
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


searchPMF_para<-function(pimlist,spectrumlist,ppm,BPPARAM=bpparam()){
  pimresultlist<-pimlist
  message("Start PMF search")
  pimresultlist<-bplapply(pimlist,PMFsum,spectrumlist,ppm,BPPARAM = BPPARAM)

  return(pimresultlist)
}

searchPMF_data_frame<-function(pimlist,spectrumlist,ppm,BPPARAM = bpparam()){
  pimresultlist<-pimlist
  message("Start PMF search")
  pimresultlist<-bplapply(pimlist,PMFsum,spectrumlist,ppm,BPPARAM = BPPARAM)
  return(pimresultlist)
}

recheck_peptide_score<-function(formula="AGLQFPVGR",peaklist=read.csv(paste0(getwd(),"/Spectrum 2 .csv"))){
  peaklist
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

SCORE_PMF<-function(formula,peaklist,isotopes=NULL,threshold=1,charge=1,ppm=5,print.graphic=F,output.list=F,outputfile=NULL,score_method="SQRTP",similarity_plot_res=72,anno_info="",isolabel="Normal",output_monomz=F){
  suppressMessages(suppressWarnings(require(rcdk)))
  suppressMessages(suppressWarnings(require(rcdklibs)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(rJava)))
  suppressMessages(suppressWarnings(require(grid)))
  if (is.null(isotopes)){data("isotopes")}


  formula<-as.character(formula)

  # define the matching score calculation function
  Spectrum_scoring <- function(spec.top, spec.bottom, t = 0.25, b = 0, top.label = NULL,
                               bottom.label = NULL, xlim = c(50, 1200), x.threshold = 0, print.alignment = FALSE,
                               print.graphic = F, output.list = F,score_method="SQRT",Peak_intensity=0,outputfile=NULL,formula=NULL,pattern_ppm=NULL,ppm_error=0,res=72,anno_info="") {


    top_tmp <- data.frame(mz = spec.top[,1], intensity = spec.top[,2])
    top_tmp$normalized <- round((top_tmp$intensity / max(top_tmp$intensity)) ,digits = 3)
    top_plot <- data.frame(mz = top_tmp$mz, intensity = top_tmp$normalized)   # data frame for plotting spectrum
    top <- subset(top_plot, top_plot$intensity > b)   # data frame for similarity score calculation

    bottom_tmp <- data.frame(mz = spec.bottom[,1], intensity = spec.bottom[,2])
    bottom_tmp$normalized <- round((bottom_tmp$intensity / max(bottom_tmp$intensity)) ,digits = 3)
    bottom_plot <- data.frame(mz = bottom_tmp$mz, intensity = bottom_tmp$normalized)   # data frame for plotting spectrum
    bottom <- subset(bottom_plot, bottom_plot$intensity > b)   # data frame for similarity score calculation

    alignment <- merge(top, bottom, by = 1, all = TRUE)
    alignment[,c(2,3)][is.na(alignment[,c(2,3)])] <- 0   # convert NAs to zero (R-Help, Sept. 15, 2004, John Fox)
    names(alignment) <- c("mz", "intensity.top", "intensity.bottom")

    if(print.alignment == TRUE) {
      print(alignment)
    }

    ## similarity score calculation

    match_score <- sum(`&`(alignment$intensity.top>0,alignment$intensity.bottom>0))/sum(alignment$intensity.top>0)

    if (score_method=="SQRT"){
      u <- alignment[,2]; v <- alignment[,3]
      similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2)) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))))

    }else if (score_method=="SQRTP"){
      u <- alignment[,2]; v <- alignment[,3]
      similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2) / (sum(u^2) * sum(v^2)))))
      match_score <- log(sum(`&`(alignment$intensity.top>0,alignment$intensity.bottom>0))/sum(alignment$intensity.top>0))
      similarity_score <- similarity_score + match_score
    }else if (score_method=="balanced-SQRT"){
      u <- alignment[,2]; v <- alignment[,3]
      if (Peak_intensity<10) (Peak_intensity=10)

      similarity_score <- -log(as.vector(sqrt(sum(((u-v)*u)^2)) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))) * log(Peak_intensity,10)

    }else if (score_method=="Equal-SQRT"){
      u <- alignment[,2]; v <- alignment[,3]
      score=-log(v/u)
      score=score[score!=Inf]
      similarity_score=-log(sum(abs(score))/length(score))

    }    else if (score_method=="Equal-intensity-SQRT"){
      u <- alignment[,2]; v <- alignment[,3]
      similarity_score <- as.vector((u %*% v) * length(v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))*log(Peak_intensity)) )
    }    else if (score_method=="Mix-SQRT"){

      alignment_0=ifelse(alignment==0,0,1)

      matching_score = sum(alignment_0[,2]==alignment_0[,3])/nrow(alignment_0)

      alignment1<-alignment[alignment_0[,2]==alignment_0[,3],]

      u <- alignment1[,2]; v <- alignment1[,3]

      similarity_score <- as.vector((u %*% v) / (sqrt(sum(u^2)) * sqrt(sum(v^2))))
    }

    ## generate plot
    similarity_score=round((similarity_score-ppm_error),digits = 7)
    if(print.graphic == TRUE) {
      suppressMessages(suppressWarnings(require(ggplot2)))
      tempfilename=outputfile
      if (!is.null(dev.list())) dev.off()
      png(tempfilename, bg = "white", width = 6.67, height = 6.67, res = res, units = "in")

      plot.new()
      plot.window(xlim = xlim, ylim = c(-125, 125))
      ticks <- c(-100, -50, 0, 50, 100)
      for(i in 1:length(top_plot$mz)) {
        lines(rep(top_plot$mz[i], 2), c(0, top_plot$intensity[i]*100), col = "blue")
        points(top_plot$mz[i],top_plot$intensity[i]*100, col = "darkblue")
      }
      for(i in 1:length(bottom_plot$mz)) {

        lines(rep(bottom_plot$mz[i], 2), c(0, -bottom_plot$intensity[i]*100), col = "red")
        if (bottom_plot$intensity[i]!=0){
          points(bottom_plot$mz[i],-bottom_plot$intensity[i]*100, col = "darkred")
          if(length(pattern_ppm$plotppm[pattern_ppm$mz==bottom_plot$mz[i]])==1){
            text(bottom_plot$mz[i],-bottom_plot$intensity[i]*100-7, pattern_ppm$plotppm[pattern_ppm$mz==bottom_plot$mz[i]],cex=0.6,col="gray21")
          }}
      }
      axis(2, at = ticks, labels = abs(ticks), pos = xlim[1], ylab = "Intensity")
      axis(1, pos = -125)
      lines(xlim, c(0,0))
      rect(xlim[1], -125, xlim[2], 125)
      mtext("m/z", side = 1, line = 2)
      mtext("intensity (%)", side = 2, line = 2)
      plot.window(xlim = c(0, 20), ylim = c(-10, 10))
      text(10, 9, top.label)
      text(10, -9, bottom.label)
      if(F){
        atome_replace<-isotopes$element[isotopes$isotope==as.character(isolabel)]
        formula_label=stringr::str_replace(formula,atome_replace[1],atome_replace[2])
      }else{
        formula_label=formula
      }
      
      if (anno_info==""){
        mtext(paste("Formula:",formula_label,"Score:",round(similarity_score,digits = 2),"Method:",score_method),cex=0.9)
      }else{
        mtext(paste(anno_info,"\nFormula:",formula_label,"Score:",round(similarity_score,digits = 2),"Method:",score_method),cex=0.9)
        
      }
      if(!is.null(pattern_ppm$norm_ppm)){
        text(pattern_ppm[,1], rep(-8.5,nrow(pattern_ppm)), pattern_ppm$norm_ppm)
      }
      dev.off()
    }
    
    if (sum(abs(similarity_score)==Inf,is.nan(similarity_score),is.na(similarity_score),na.rm = T)) similarity_score = 0

    if(output.list == FALSE) {

      return(similarity_score)

    }

    if(output.list == TRUE) {

      return(list(similarity.score = similarity_score,
                  plot = p))

    }

  }

  # define the function to merge the un-resolvable mz of a given fomula by instrument resolution setting
  isopattern_ppm_filter_peaklist<-function(pattern,ppm,threshold=0.001,verbose=F){

    org_feature=nrow(pattern)

    pattern<-as.data.frame(pattern)

    pattern_ppm=as.numeric(as.character(pattern[,1]))

    pattern_ppm_delta=numeric()

    filtered_pattern<-pattern[1,]

    for (i in 1:(length(pattern_ppm)-1)){

      pattern_ppm_delta[i]=(pattern_ppm[i+1]-pattern_ppm[i])/pattern_ppm[i]

      if (pattern_ppm_delta[i]>(ppm/1000000)){
        filtered_pattern<-rbind(filtered_pattern,pattern[i+1,])
      } else {

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

  # define the function to merge the un-resolvable mz of a given fomula by instrument resolution setting
  isopattern_ppm_filter<-function(pattern,ppm){

    pattern_ppm=as.numeric(as.character(pattern[,1]))

    pattern_ppm_delta=numeric()


    filtered_pattern<-pattern[1,]
    filtered_pattern<-t(data.frame(filtered_pattern,stringsAsFactors = F))

    for (i in 1:(length(pattern_ppm)-1)){

      pattern_ppm_delta[i]=(pattern_ppm[i+1]-pattern_ppm[i])/pattern_ppm[i]

      if (length(pattern_ppm_delta[i] > (ppm/1e+06))==1 && pattern_ppm_delta[i]>(ppm/1000000)){
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
  
  #generate the isotopic pattern
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

  if (ppm>=25) {
    instrument_ppm=50
  }else{
    instrument_ppm=max(c(8,ppm),na.rm = T)
  }

  #Filter and merge the isotopic pattern using instrument resolution
  pattern=pattern[[formula]]
  #pattern=isopattern_ppm_filter(pattern = pattern[,1:2], ppm=instrument_ppm)
  #pattern=pattern[pattern[,2]>=1,]
  pattern=isopattern_ppm_filter(pattern = matrix(pattern[,1:2],ncol=2), ppm=instrument_ppm)
  pattern= matrix(pattern[topN_feature(pattern[,2],10),],ncol=2)
  pattern<-matrix(pattern[,1:2],ncol=2)
  #peptide similarity score calculation
  spectrumintensity=unlist(lapply(1:nrow(pattern), function(x,pattern,peaklist,ppm){
    PMF_spectrum_intensity<-as.numeric(peaklist[between(peaklist$m.z,pattern[x,1]*(1-ppm/1000000),pattern[x,1]*(1+ppm/1000000)),2])
    if(sum(!is.na(PMF_spectrum_intensity),length(PMF_spectrum_intensity)!=0)>1) sum(PMF_spectrum_intensity,na.rm = T) else 0
  },pattern,peaklist,instrument_ppm))

  spectrum=data.frame(mz=pattern[,1],Intensity=spectrumintensity)
  Peak_intensity<-max(spectrum$Intensity)

  spectrum_pk<-lapply(1:nrow(pattern), function(x,pattern,peaklist,instrument_ppm){
    peaklist[between(peaklist$m.z,pattern[x,1]*(1-(instrument_ppm)/1000000),pattern[x,1]*(1+(instrument_ppm)/1000000)),]
  },pattern,peaklist,instrument_ppm)
  spectrum_pk<-do.call(rbind,spectrum_pk)
  spectrum_pk<-unique(spectrum_pk)
  if (nrow(spectrum_pk)>1){
    spectrum_pk=isopattern_ppm_filter_peaklist(pattern = spectrum_pk, ppm=ppm)
  }

  #peptide ppm-error score calculation
  pattern_ppm<-do.call(rbind,(lapply(1:nrow(pattern), function(x,pattern,spectrum_pk,instrument_ppm){
    PMF_spectrum<-spectrum_pk[between(spectrum_pk$m.z,pattern[x,1]*(1-instrument_ppm/1000000),pattern[x,1]*(1+instrument_ppm/1000000)),]
    if(nrow(PMF_spectrum)>=1){
      spectrummz<-sum(PMF_spectrum[,1]*PMF_spectrum[,2])/sum(PMF_spectrum[,2])
      return(data.frame(mz=pattern[x,1],delta_ppm=(spectrummz-pattern[x,1])/pattern[x,1]*1000000))
    }
  },pattern,spectrum_pk,instrument_ppm)))

  if (nrow(spectrum_pk)>1 && nrow(pattern_ppm)>1 && !is.null(pattern_ppm)){
    pattern_ppm$norm_ppm<-pattern_ppm$delta_ppm-mean(pattern_ppm$delta_ppm)
    meanppm=mean(abs(pattern_ppm$norm_ppm))
    pattern_ppm$plotppm<-ifelse(abs(pattern_ppm$norm_ppm)>ppm,paste0(">",ppm),round((pattern_ppm$norm_ppm),digits = 1))
  }else if (nrow(spectrum_pk)==1 && nrow(pattern_ppm)==1 && !is.null(pattern_ppm)){
    pattern_ppm$norm_ppm<-pattern_ppm$delta_ppm
    meanppm=mean(abs(pattern_ppm$norm_ppm))
    pattern_ppm$plotppm<-ifelse(abs(pattern_ppm$norm_ppm)>ppm,paste0(">",ppm),round((pattern_ppm$norm_ppm),digits = 1))
  }else{
    meanppm=instrument_ppm
  }

  ppm_error=abs(pnorm(meanppm/ppm)-0.5)

  #peptide final score calculation

  finalscore=Spectrum_scoring(spec.top = pattern,spec.bottom = spectrum,b=0,t=mean(pattern[,1])*ppm/1000000,top.label = "Theoretical",
                              score_method=score_method,bottom.label = "Observed spectrum",
                              xlim = c(min(range(pattern[,1],na.rm = T))-1,max(range(pattern[,1],na.rm = T))+1),
                              output.list=output.list,print.graphic = print.graphic,
                              outputfile = outputfile,
                              Peak_intensity=Peak_intensity,res = similarity_plot_res,
                              formula=formula,pattern_ppm=pattern_ppm,ppm_error=ppm_error,
                              anno_info = anno_info)

  Peak_intensity<-sum(spectrum$Intensity)
  return(as.numeric(data.frame(finalscore,meanppm,Peak_intensity)))
}

Do_PMF_search<-function(peaklist,Peptide_Summary_searchlist,BPPARAM=bpparam(),ppm=2){

  #Summarize the intensity for a given candidates between m/z tolerance window

  peaklist<-peaklist[peaklist$intensities!=0,]
  colnames(peaklist)<-c("m.z","intensities")
  mzrange=c(min(peaklist$m.z),max(peaklist$m.z))
  Peptide_Summary_searchlist_mz=NULL
  uniquemz=as.numeric(as.character(unique(Peptide_Summary_searchlist$mz)))
  Peptide_Summary_searchlist_mz$mz=uniquemz[`&`(uniquemz>mzrange[1],uniquemz<mzrange[2])]
  Peptide_Summary_searchlist_mz$Intensity=rep(0,length(Peptide_Summary_searchlist_mz$mz))
  Peptide_Summary_searchlist_mz=as.data.frame(Peptide_Summary_searchlist_mz)
  Peptide_Summary_searchlist_mz<-IMS_analysis_fun_2(Peptide_Summary_searchlist=Peptide_Summary_searchlist_mz,peaklist=peaklist,ppm=ppm,BPPARAM=BPPARAM,mzrange=mzrange)
  Peptide_Summary_searchlist_mz=as.data.frame(Peptide_Summary_searchlist_mz,stringsAsFactors = F)
  mz_feature_list<-Peptide_Summary_searchlist_mz[Peptide_Summary_searchlist_mz$Intensity>0,]

  return(mz_feature_list)
}

Peptide_regions_Summary_fun<-function(Peptide_Summary_file_regions){

  library(data.table)
  peptideregion<-list()
  unique_region<-unique(Peptide_Summary_file_regions$Region)
  peptideregion<-lapply(unique_region, function(x,Peptide_Summary_file_regions){
    Peptide_Summary_file_regions[Peptide_Summary_file_regions$Region==x,]
  },Peptide_Summary_file_regions)

  names(peptideregion)<-unique_region

  Peptide_regions_Summary<-data.frame()

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


  return(Peptide_regions_Summary)
}

FDR_cutoff_plot_protein<-function(Protein_feature_result,FDR_cutoff=0.1,plot_fdr=T,outputdir=paste0(datafile[z] ," ID/",SPECTRUM_batch),FDR_strip=500,adjust_score = F){

  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(zoo)))
  Protein_feature_result=data_test_rename(c("isdecoy","Proscore"),Protein_feature_result)


  #formating inout df and vars
  Protein_feature_result$isdecoy<-factor(Protein_feature_result$isdecoy)

  Protein_feature_result$Proscore=Protein_feature_result$Proscore

  #Disabled for later revision
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

  #isolate target score and decoy score
  target_Proscore<-Protein_feature_result$Proscore[Protein_feature_result$isdecoy!=1]
  decoy_Proscore<-Protein_feature_result$Proscore[Protein_feature_result$isdecoy==1]
  target_Proscore[is.na(target_Proscore)]=0
  decoy_Proscore[is.na(decoy_Proscore)]=0
  target_Proscore[target_Proscore==Inf]<-0
  decoy_Proscore[decoy_Proscore==Inf]<-0
  if (length(unique(Protein_feature_result$Score))==1 && !(Protein_feature_result$isdecoy %in% 1)) {
    return(unique(Protein_feature_result$Score))
  }else if (length(unique(Protein_feature_result$Score))==1 && (Protein_feature_result$isdecoy %in% 1)){
    return(unique(Protein_feature_result$Score)+100)
  }

  #define FDR test score bins
  breaks = seq(min(Protein_feature_result$Proscore,na.rm = T), max(Protein_feature_result$Proscore,na.rm = T), by=abs(max(Protein_feature_result$Proscore,na.rm = T)/FDR_strip))
  target_Proscore.cut = cut(target_Proscore, breaks, right=T)
  decoy_Proscore.cut = cut(decoy_Proscore, breaks,right=T)
  target_Proscore.freq = table(target_Proscore.cut)
  decoy_Proscore.freq = table(decoy_Proscore.cut)
  target_Proscore.freq= target_Proscore.freq[rev(names(target_Proscore.freq))]
  target_Proscore.freq0 = c(0, cumsum(target_Proscore.freq))
  decoy_Proscore.freq= decoy_Proscore.freq[rev(names(decoy_Proscore.freq))]
  decoy_Proscore.freq0 = c(0, cumsum(decoy_Proscore.freq))


  #define FDR test score bins
  FDR=decoy_Proscore.freq0/target_Proscore.freq0
  df=data.frame(breaks=rev(breaks),target_Proscore.freq0=target_Proscore.freq0 ,decoy_Proscore.freq0=decoy_Proscore.freq0,FDR=FDR)
  df$FDR_m.av<-rollmean(df$FDR, 3,fill = list(0, NaN, NaN))

  #Find the lowest score cut-off that make the final FDR above a given FDR threshold
  Proscore_cutoff<-min(c(df$breaks[(df$FDR_m.av<=FDR_cutoff)==T],df$breaks[(df$FDR<=FDR_cutoff)==T]),na.rm = T)

  #Output protein scoring result
  message(paste("Protein score cutoff:",round(Proscore_cutoff,digits = 4)))
  if (!is.null(outputdir) && plot_fdr){

    png(paste0(outputdir,"/Protein_FDR.png"))
    plot.new()
    try(plot(df$breaks, df$FDR_m.av,            # plot the data
             main="FDR plot",  # main title
             xlab="Matching Proscore",        # x−axis label
             ylab="FDR"))   # y−axis label
    try(lines(df$breaks, df$FDR_m.av))
    try(abline(v=Proscore_cutoff))
    dev.off()

    write.csv(df,paste0(outputdir,"/PROTEIN_FDR.CSV"),row.names = F)

  }

  if (!is.null(outputdir) && plot_fdr){
    Protein_feature_result_plot<-Protein_feature_result[,c("isdecoy","Proscore","Score")]
    target_decoy<-factor(ifelse(Protein_feature_result$isdecoy==0,"Target","Decoy"),levels = c("Target","Decoy"))
    Protein_feature_result_plot$target_decoy=target_decoy
    Protein_feature_result_plot$color=ifelse(Protein_feature_result$isdecoy==0,"red","blue")

    png(paste0(outputdir,"/PROTEIN_Score_histogram.png"))
    p<-ggplot(data=Protein_feature_result_plot,aes(x=Proscore,color=target_decoy, fill=target_decoy)) + geom_histogram( fill="white",alpha=0.5, bins = 200)  +
      ggtitle("Protein score vs Counts") + xlim(-0.05,max(quantile(Protein_feature_result_plot$Proscore, c(0.99), na.rm=T),Proscore_cutoff+0.05)) +
      xlab("Protein score") + ylab("Counts") + labs(fill = "Is_Decoy") + theme_classic()+ geom_vline(xintercept = Proscore_cutoff)  #+ facet_grid(target_decoy ~ .)
    print(p)
    dev.off()
    if(F){
      pep_cali<-Protein_feature_result_plot$Proscore[grep("Pep_cali",Protein_feature_result_plot$desc,ignore.case = Proscore_cutoff)]
      png(paste0("PROTEIN_Score_histogram_full_scale.png"))
      gg_color_hue <- function(n) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100)[1:n]
      }
      p<-ggplot(data=Protein_feature_result_plot,aes(x=Proscore,color=target_decoy, fill=target_decoy)) + geom_histogram( fill="white",alpha=0.5, bins = 200)  +
        ggtitle("Protein score vs Counts")  +
        xlab("Protein score") + ylab("Counts") + labs(fill = "Is_Decoy") + theme_classic()+
        scale_color_manual(values=c("grey50",gg_color_hue(2)[2:1])) +
        geom_vline(xintercept = Proscore_cutoff,linetype="dashed") + #geom_segment(aes(x=Proscore,xend = Proscore,y=rep(0.3,length(Proscore)),yend=rep(0.1,length(Proscore))))
        geom_segment(data=data.frame(pep_cali=pep_cali),aes(x=pep_cali,color=" Peptide calibrant",fill="black",xend = pep_cali,y=rep(300,length(pep_cali)),yend=rep(0.1,length(pep_cali))),arrow = arrow(angle = 15, ends = "last", type = "closed",length = unit(0.15, "inches")))
      #geom_vline(,xend = ,show.legend = T)# Protein_feature_result_plot$desc[grep("Pep_cali",Protein_feature_result_plot$desc,ignore.case = T)]) #+ facet_grid(target_decoy ~ .)
      print(p)
      dev.off()
    }

  }

  return(Proscore_cutoff)



}

FDR_cutoff_plot<-function(Peptide_plot_list,FDR_cutoff=0.1,FDR_strip=500,plot_fdr=F,outputdir=NULL,adjust_score=F){
  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(zoo)))
  #suppressMessages(suppressWarnings(require(FTICRMS)))
  Peptide_plot_list=HiTMaP:::data_test_rename(c("isdecoy","Score"),Peptide_plot_list)

  Peptide_plot_list$isdecoy<-factor(Peptide_plot_list$isdecoy)

  #unique_fomula<-unique(Peptide_plot_list$formula)
  Peptide_plot_list_output<-Peptide_plot_list

  #unique_fomula_ID<-unique(Peptide_plot_list[,c("mz","formula","isdecoy","Intensity","Score")])

  plot_fdr_histogram<-function(Peptide_plot_list,plot_name="target-decoy",outputdir=outputdir){
    suppressMessages(suppressWarnings(require(dplyr)))
    if (nrow(Peptide_plot_list)!=0) {
      Peptide_plot_list_plot<-Peptide_plot_list[,c("isdecoy","Score")]
      target_decoy<-factor(ifelse(Peptide_plot_list_plot$isdecoy==0,"Target","Decoy"),levels = c("Target","Decoy"))
      Peptide_plot_list_plot$target_decoy=target_decoy
      Peptide_plot_list$target_decoy=target_decoy
      if(F){ png(paste0(outputdir,"/Peptide_Score_histogram_",plot_name,".png"))
        p<-ggplot(data=Peptide_plot_list_plot,aes(x=Peptide_plot_list_plot$Score,color=target_decoy, fill=target_decoy)) +
          geom_histogram( fill="white",alpha=0.5, bins = 50, position="Dodge")  +
          ggtitle("Peptide score vs Counts") +
          xlab("Score") + ylab("Counts") + labs(fill = "Is_Decoy") + theme_classic() #+ facet_grid(target_decoy ~ .)
        print(p)
        dev.off() }

      mu <- Peptide_plot_list_plot %>% group_by(target_decoy) %>% summarize(mean=mean(Score))

      png(paste0(outputdir,"/Peptide_Score_histogram_",plot_name,".png"))
      p<-ggplot(Peptide_plot_list_plot, aes(x=Peptide_plot_list_plot$Score, color=target_decoy, fill=target_decoy)) +
        geom_histogram(aes(y=after_stat(density)), position="Dodge", alpha=0.5, bins = 50)+
        geom_density(alpha=0.6)+ xlim(quantile(Peptide_plot_list_plot$Score, c(0.005,0.995)))+
        geom_vline(data=mu, aes(xintercept=mu$mean, color=mu$target_decoy), linetype="dashed")+
        ggtitle("Peptide score vs Counts") +
        xlab("Score") + ylab("Counts") +
        #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_classic()
      print(p)
      dev.off()


      if (log(nrow(Peptide_plot_list_plot))>0){
        png(paste0(outputdir,"/Matching_Score_vs_mz_",plot_name,".png"))
        alpha=1/20/(round(log(nrow(Peptide_plot_list_plot),base = 2))/40)
        if (alpha>=1 ) alpha=1
        p<-ggplot(data=Peptide_plot_list,aes(x=mz,y=Score,colour=target_decoy)) + geom_point(size=1,alpha=alpha,aes(colour = target_decoy)) +
          ggtitle("Matching score vs mz") +
          xlab("mz") + ylab("Matching score") + labs(fill = "isdecoy")+ theme(axis.text.x = element_text(angle=45))
        print(p)
        dev.off()

      }

    }
  }
  #Peptide_plot_list=unique_fomula_ID

  if (!is.null(outputdir) && plot_fdr){

    plot_fdr_histogram(Peptide_plot_list,outputdir=outputdir)

  }


  if(F){

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

  if(adjust_score){

    Peptide_plot_list_sep<-Peptide_plot_list_output[Peptide_plot_list_output$isdecoy==0,]
    Peptide_plot_list_sep_d<-Peptide_plot_list_output[Peptide_plot_list_output$isdecoy==1,]
    Peptide_plot_list_mer<-merge(Peptide_plot_list_sep,Peptide_plot_list_sep_d,by=names(Peptide_plot_list_sep)[!names(Peptide_plot_list_sep) %in% c("isdecoy","Score" ,"Rank")])
    Peptide_plot_list_mer$Score=Peptide_plot_list_mer$Score.x-Peptide_plot_list_mer$Score.y
    Peptide_plot_list_mer$isdecoy=0
    Peptide_plot_list_adj<-Peptide_plot_list_mer
    Peptide_plot_list_mer$Score=0
    Peptide_plot_list_mer$isdecoy=1
    Peptide_plot_list_adj<-rbind(Peptide_plot_list_adj,Peptide_plot_list_mer)
    Peptide_plot_list_adj$isdecoy<-as.factor(Peptide_plot_list_adj$isdecoy)
    plot_fdr_histogram(Peptide_plot_list_adj,plot_name = "adj",outputdir=outputdir)
    Peptide_plot_list<-Peptide_plot_list_adj


    plot(density(Peptide_plot_list_adj$Score[Peptide_plot_list_adj$isdecoy==0]))
    den<-density(Peptide_plot_list_adj$Score[Peptide_plot_list_adj$isdecoy==0])
    find_modes<- function(x) {
      modes <- NULL
      for ( i in 2:(length(x)-1) ){
        if ( (x[i] > x[i-1]) & (x[i] > x[i+1]) ) {
          modes <- c(modes,i)
        }
      }
      if ( length(modes) == 0 ) {
        modes = 'This is a monotonic distribution'
      }
      return(modes)
    }
    mymodes_indices <- find_modes(den$y)
    den$x[find_modes(den$y)]

    mm<-resolve_multi_modes(Peptide_plot_list_adj$Score[Peptide_plot_list_adj$isdecoy==0],adjust = 2)

    unique_fomula_ID<-unique(Peptide_plot_list[,c("mz","formula","isdecoy","Intensity","Score")])
    target_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy!=1]
    decoy_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy==1]
    target_score[is.na(target_score)]=0
    decoy_score[is.na(decoy_score)]=0
    target_score[target_score==Inf]<-0
    decoy_score[decoy_score==Inf]<-0

    breaks = seq( min(Peptide_plot_list$Score), max(Peptide_plot_list$Score) + (max(Peptide_plot_list$Score)-min(Peptide_plot_list$Score))/FDR_strip , by=(max(Peptide_plot_list$Score)-min(Peptide_plot_list$Score))/FDR_strip)
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
    #df$logscore=log(df$breaks)

    Score_cutoff<-max(as.numeric(mm$modes))-0.75

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



  }else{

    unique_fomula_ID<-unique(Peptide_plot_list[,c("mz","formula","isdecoy","Intensity","Score")])
    target_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy!=1]
    decoy_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy==1]
    target_score[is.na(target_score)]=0
    decoy_score[is.na(decoy_score)]=0
    target_score[target_score==Inf]<-0
    decoy_score[decoy_score==Inf]<-0
    if (length(unique(Peptide_plot_list$Score))==1 && !(Peptide_plot_list$isdecoy %in% 1)) {
      return(unique(Peptide_plot_list$Score))
    }else if (length(unique(Peptide_plot_list$Score))==1 && (Peptide_plot_list$isdecoy %in% 1)){
      return(unique(Peptide_plot_list$Score)+100)
    }
    if (length(unique(Peptide_plot_list$Score))==1) return(unique(Peptide_plot_list$Score))
    breaks = seq( min(Peptide_plot_list$Score), max(Peptide_plot_list$Score) + (max(Peptide_plot_list$Score)-min(Peptide_plot_list$Score))/FDR_strip , by=abs(max(Peptide_plot_list$Score)-min(Peptide_plot_list$Score))/FDR_strip)
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
    #df$logscore=log(df$breaks)

    Score_cutoff<-min(c(df$breaks[(df$FDR_m.av<=FDR_cutoff)==T],df$breaks[(df$FDR<=FDR_cutoff)==T]),na.rm = T)

    if (!is.null(outputdir) && plot_fdr){
      if ((sum(is.na(df$FDR)) + sum(is.infinite(df$FDR)))<length(df$FDR)){
        png(paste0(outputdir,"/FDR.png"))
        plot.new()
        plot(df$breaks, df$FDR,            # plot the data
             main="FDR plot",  # main title
             xlab="Matching score",        # x−axis label
             ylab="FDR")   # y−axis label
        lines(df$breaks, df$FDR_m.av)
        dev.off()

      }


      write.csv(df,paste0(outputdir,"/FDR.CSV"),row.names = F)

    }
  }

  if (adjust_score){
    return(list(Score_cutoff,Peptide_plot_list))
  }else{
    return(Score_cutoff)
  }




}

plot_cpd_matching_score<-function(Peptide_plot_list,peaklist,charge,ppm,outputdir=getwd(),filename_col=c("formula","adduct"),anno_col=filename_col,isotopes="default",similarity_plot_res = 72){
  if (ppm>=25) {
    instrument_ppm=50
  }else{
    instrument_ppm=8
  }
  message("plot matching isotopic pattern")
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(dplyr)))
  Peptide_plot_list=data_test_rename(unique(c("formula","isdecoy","adduct",filename_col,anno_col)),Peptide_plot_list)
  if(unlist(isotopes)[1]=="default") data("isotopes")
  decoy_isotopes=isotopes
  decoy_isotopes[decoy_isotopes$isotope=="13C",]=data.frame(element="C",isotope="11C",mass=10.99664516,abundance=0.0107,ratioC=0,stringsAsFactors = F)
  #decoy_isotopes_N=isotopes
  decoy_isotopes[decoy_isotopes$isotope=="15N",]=data.frame(element="N",isotope="13N",mass=13.00603905,abundance=0.00364000,ratioC=4,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="2H",]=data.frame(element="H",isotope="0H",mass=0.001548286,abundance=0.00011500,ratioC=6,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="17O",]=data.frame(element="O",isotope="15O",mass=14.99069774,abundance=0.00038000,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="18O",]=data.frame(element="O",isotope="14O",mass=13.99066884,abundance=0.00205000,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="33S",]=data.frame(element="S",isotope="31S",mass=30.97268292,abundance=0.0075,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="34S",]=data.frame(element="S",isotope="30S",mass=29.9762745,abundance=0.0425,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="35S",]=data.frame(element="S",isotope="29S",mass=28.94414146,abundance=0,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="36S",]=data.frame(element="S",isotope="28S",mass=27.97706058,abundance=0.0001,ratioC=3,stringsAsFactors = F)
  if (dir.exists(outputdir)==FALSE){dir.create(outputdir)}
  #Peptide_plot_list1=Peptide_plot_list
  Peptide_plot_list=Peptide_plot_list[!is.na(Peptide_plot_list$Intensity),]
  Peptide_plot_list=Peptide_plot_list %>% group_by(.dots = unique(c("formula","isdecoy","adduct",filename_col,anno_col))) %>% summarise(Metabolite.Name=paste(Metabolite.Name,collapse = "_"))
  
  if(nrow(Peptide_plot_list)!=0){
    Peptide_plot_list_score<-NULL
    for (i in 1:nrow(Peptide_plot_list)){
      if ("Isotype" %in% anno_col){
      anno_info=paste(Peptide_plot_list[i,anno_col[anno_col!="Isotype"]],collapse = "\n")
      isolabel=Peptide_plot_list[i,"Isotype"]
      }else{
      anno_info=paste(Peptide_plot_list[i,anno_col],collapse = "\n") 
      isolabel="Normal"
      }
      
      #anno_info=stringr::str_replace_all(anno_info,";","\n")
      if(Peptide_plot_list$isdecoy[i]==0){
        Peptide_plot_list_score[i]<-try(SCORE_PMF(Peptide_plot_list$formula[i],
                                                  peaklist,isotopes=isotopes,
                                                  charge=charge,ppm=ppm,isolabel=isolabel,
                                                  print.graphic=T,output.list=F,
                                                  outputfile=paste0(outputdir,"/",paste(Peptide_plot_list[i,filename_col],collapse = "_"),"_target",".png"),
                                                  similarity_plot_res = similarity_plot_res,
                                                  anno_info=anno_info))
      
      
      }else{
        Peptide_plot_list_score[i]<-try(SCORE_PMF(Peptide_plot_list$formula[i],
                                                  peaklist,isotopes=decoy_isotopes,
                                                  charge=charge,ppm=ppm,isolabel=isolabel,
                                                  print.graphic=T,output.list=F,
                                                  outputfile=paste0(outputdir,"/",
                                                                    paste(Peptide_plot_list[i,filename_col],collapse = "_"),"_decoy",".png"),
                                                  similarity_plot_res = similarity_plot_res,
                                                  anno_info=anno_info))
      }
    }
    Peptide_plot_list$Score<-Peptide_plot_list_score
    return(Peptide_plot_list)
  }
  
  
  
}

plot_matching_score<-function(Peptide_plot_list,peaklist,charge,ppm,outputdir=getwd(),filename_col=c("formula","adduct"),anno_col=filename_col,isotopes="default",similarity_plot_res = 72){
  if (ppm>=25) {
    instrument_ppm=50
  }else{
    instrument_ppm=8
  }
  message("plot matching isotopic pattern")
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(dplyr)))
  Peptide_plot_list=data_test_rename(unique(c("formula","isdecoy","adduct",filename_col,anno_col)),Peptide_plot_list)
  if(unlist(isotopes)[1]=="default") data("isotopes")
  decoy_isotopes=isotopes
  decoy_isotopes[decoy_isotopes$isotope=="13C",]=data.frame(element="C",isotope="11C",mass=10.99664516,abundance=0.0107,ratioC=0,stringsAsFactors = F)
  #decoy_isotopes_N=isotopes
  decoy_isotopes[decoy_isotopes$isotope=="15N",]=data.frame(element="N",isotope="13N",mass=13.00603905,abundance=0.00364000,ratioC=4,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="2H",]=data.frame(element="H",isotope="0H",mass=0.001548286,abundance=0.00011500,ratioC=6,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="17O",]=data.frame(element="O",isotope="15O",mass=14.99069774,abundance=0.00038000,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="18O",]=data.frame(element="O",isotope="14O",mass=13.99066884,abundance=0.00205000,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="33S",]=data.frame(element="S",isotope="31S",mass=30.97268292,abundance=0.0075,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="34S",]=data.frame(element="S",isotope="30S",mass=29.9762745,abundance=0.0425,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="35S",]=data.frame(element="S",isotope="29S",mass=28.94414146,abundance=0,ratioC=3,stringsAsFactors = F)
  decoy_isotopes[decoy_isotopes$isotope=="36S",]=data.frame(element="S",isotope="28S",mass=27.97706058,abundance=0.0001,ratioC=3,stringsAsFactors = F)
  if (dir.exists(outputdir)==FALSE){dir.create(outputdir)}
  #Peptide_plot_list1=Peptide_plot_list
  Peptide_plot_list=Peptide_plot_list[!is.na(Peptide_plot_list$Intensity),]
  Peptide_plot_list=Peptide_plot_list %>% group_by(.dots = unique(c("formula","isdecoy","adduct",filename_col,anno_col))) %>% summarise(Peptide=paste(Peptide,collapse = "_"))
  
  if(nrow(Peptide_plot_list)!=0){
    Peptide_plot_list_score<-NULL
    for (i in 1:nrow(Peptide_plot_list)){
      if ("Isotype" %in% anno_col){
        anno_info=paste(Peptide_plot_list[i,anno_col[anno_col!="Isotype"]],collapse = "\n")
        isolabel=Peptide_plot_list[i,"Isotype"]
      }else{
        anno_info=paste(Peptide_plot_list[i,anno_col],collapse = "\n") 
        isolabel="Normal"
      }
      
      #anno_info=stringr::str_replace_all(anno_info,";","\n")
      if(Peptide_plot_list$isdecoy[i]==0){
        Peptide_plot_list_score[i]<-try(SCORE_PMF(Peptide_plot_list$formula[i],
                                                  peaklist,isotopes=isotopes,
                                                  charge=charge,ppm=ppm,isolabel=isolabel,
                                                  print.graphic=T,output.list=F,
                                                  outputfile=paste0(outputdir,"/",paste(Peptide_plot_list[i,filename_col],collapse = "_"),"_target",".png"),
                                                  similarity_plot_res = similarity_plot_res,
                                                  anno_info=anno_info))
        
        
      }else{
        Peptide_plot_list_score[i]<-try(SCORE_PMF(Peptide_plot_list$formula[i],
                                                  peaklist,isotopes=decoy_isotopes,
                                                  charge=charge,ppm=ppm,isolabel=isolabel,
                                                  print.graphic=T,output.list=F,
                                                  outputfile=paste0(outputdir,"/",
                                                                    paste(Peptide_plot_list[i,filename_col],collapse = "_"),"_decoy",".png"),
                                                  similarity_plot_res = similarity_plot_res,
                                                  anno_info=anno_info))
      }
    }
    Peptide_plot_list$Score<-Peptide_plot_list_score
    return(Peptide_plot_list)
  }
  
  
  
}


protein_scoring<-function(Protein_feature_list,
                          Peptide_plot_list_rank,
                          scoretype=c("median","sum","mean"),
                          BPPARAM = bpparam(),
                          protein_nr_grouping=T,
                          prioritize_protein=T,
                          compete_decoy=T,
                          peptide_ID_filter=2,
                          use_top_rank=NULL){

  message("Start protein scoring...")
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(IRanges)))

  protein_coverage<-function(query_protein_list,Protein_feature_list_rank){
    suppressMessages(suppressWarnings(require(IRanges)))
    peptides_entry<-Protein_feature_list_rank[`&`(Protein_feature_list_rank$Protein==query_protein_list$Protein,Protein_feature_list_rank$isdecoy==query_protein_list$isdecoy),c("start","end","pro_end")]
    coverage_ranges<-IRanges(start = peptides_entry$start,end = peptides_entry$end,width = 1+peptides_entry$end-peptides_entry$start)
    coverage_ranges<-coverage_ranges[coverage_ranges@start>0,]
    cov <- coverage(coverage_ranges)
    cov <- as.vector(cov)
    cov_len<-length(which(cov!=0))
    coverage_percentage=cov_len/peptides_entry$pro_end[1]
  }

  sum_pro_pep_count_fun<-function(Protein_feature_list_rank){
    suppressMessages(suppressWarnings(require(data.table)))
    suppressMessages(suppressWarnings(require(dplyr)))
    sum_pro_pep_count<- Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarise(peptide_count=length(unique(Peptide)))
    return(sum_pro_pep_count)
  }

  Sum_protein_info<-function(Protein_feature_list_rank,scoretype="mean"){
    suppressMessages(suppressWarnings(require(data.table)))
    suppressMessages(suppressWarnings(require(dplyr)))
    e <- environment()
    p <- parent.env(e)

    sum_pro_pep_count<- Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarise(peptide_count=length(unique(Peptide)))



    if (scoretype=="sum_wi_int_norm"){

      sum_pro_int<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Intensity=mean(Intensity))

      sum_pro_score<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Score=sum(Score*(Intensity))/sum((Intensity)))

    }else if(scoretype=="mean_wi_int_norm"){

      sum_pro_int<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Intensity=mean(Intensity))

      sum_pro_score<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Score=mean(Score*(Intensity),na.rm=T)/mean((Intensity),na.rm=T))

    }else if(scoretype=="mean"){

      sum_pro_int<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Intensity=mean(Intensity))

      sum_pro_score<-Protein_feature_list_rank %>% group_by(.dots=c("Protein","isdecoy")) %>% summarize(Score=mean(Score*rank(Intensity),na.rm=T)/mean(rank(Intensity)))

    }
    assign((("sum_pro_int")), sum_pro_int, env = p)
    assign((("sum_pro_score")), sum_pro_score, env = p)
    assign((("sum_pro_pep_count")), sum_pro_pep_count, env = p)
  }


  # formating the input dataframe
  Protein_feature_list<-as.data.frame(Protein_feature_list)
  Peptide_plot_list_rank<-as.data.frame(Peptide_plot_list_rank)
  Peptide_plot_list_rank$Modification[is.na(Peptide_plot_list_rank$Modification)]<-""
  Protein_feature_list$Modification[is.na(Protein_feature_list$Modification)]<-""
  Peptide_plot_list_rank<-Peptide_plot_list_rank[Peptide_plot_list_rank$Intensity>0,]
  if (!is.null(use_top_rank)){
    Peptide_plot_list_rank<-Peptide_plot_list_rank[,Peptide_plot_list_rank$Rank<=use_top_rank]
  }

  # link the peptide result to the protein grouping information
  Protein_feature_list_rank=merge(Protein_feature_list,Peptide_plot_list_rank,by=intersect(colnames(Protein_feature_list),colnames(Peptide_plot_list_rank)))
  Protein_feature_list_rank$Peptide<-as.character(Protein_feature_list_rank$Peptide)

  # Define the protein score factor
  # Define protein coverage score
  Sum_protein_info(Protein_feature_list_rank,scoretype=scoretype)

  #message(length(unique(Protein_feature_list_rank$Protein)))
  # Filter protein candidates
  Protein_feature_list_rank<-Protein_feature_list_rank[Protein_feature_list_rank$Protein %in% sum_pro_pep_count$Protein[sum_pro_pep_count$peptide_count>=peptide_ID_filter],]
  #message(length(unique(Protein_feature_list_rank$Protein)))
  Protein_feature_list_rank_filtered<-Protein_feature_list_rank
  if(nrow(Protein_feature_list_rank)==0){

    df <- data.frame(matrix(ncol = 12, nrow = 0))
    x <- c("Protein","isdecoy","Intensity","Score","peptide_count","Protein_coverage","Protein_coverage.Protein","Protein_coverage.isdecoy",
           "Protein_coverage.Protein_coverage","Intensity_norm","Proscore","desc")
    colnames(df) <- x

    #return(list(df,Protein_feature_list_rank))
    return(list(Protein_feature_result,
                Protein_feature_list_rank,
                Protein_feature_list_rank,
                Protein_feature_list_rank,
                Protein_feature_list_rank))
  }

  # get protein candidates index from global enviornment
  Index_of_protein_sequence<-get("Index_of_protein_sequence", envir = .GlobalEnv)
  message("perform Peptide feature grouping...")


  # prioritize protein ID with higher coverage if same annotated peptide set found in multiple proteins
  if(prioritize_protein){

    Sum_protein_info(Protein_feature_list_rank,scoretype=scoretype)

    query_protein = sum_pro_int[,c("Protein","isdecoy")]

    query_protein_list = split(query_protein, seq(nrow(query_protein)))
    #sum_pro_coverage <-unlist(bplapply(query_protein_list,protein_coverage,Protein_feature_list_rank[,c("isdecoy","Protein","start","end","pro_end")],BPPARAM = BPPARAM))
    sum_pro_coverage <-unlist(lapply(query_protein_list,protein_coverage,Protein_feature_list_rank[,c("isdecoy","Protein","start","end","pro_end")]))
    sum_pro_coverage=data.frame(Protein=(sum_pro_int$Protein),isdecoy=sum_pro_int$isdecoy,Protein_coverage=sum_pro_coverage)
    sum_pro_pep_count<-merge(sum_pro_pep_count,sum_pro_coverage,by=c("Protein","isdecoy"),sort=F)
    Protein_feature_list_rank$peptide_count<-NULL
    Protein_feature_list_rank$Protein_coverage<-NULL
    Protein_feature_list_rank<-merge(Protein_feature_list_rank,sum_pro_pep_count,by=c("Protein","isdecoy"))
    mz_max_peptide<-Protein_feature_list_rank %>% group_by(mz) %>% summarize(Protein_coverage=max(Protein_coverage))
    Protein_feature_list_rank_trim<-merge(mz_max_peptide,Protein_feature_list_rank,by=c("mz","Protein_coverage"))
    Protein_feature_list_rank<-Protein_feature_list_rank_trim
  }


  # generate protein grouping info from matched protein ID list
  Protein_feature_list_rank$peptide_count<-NULL
  Protein_feature_list_rank<-as.data.table(Protein_feature_list_rank)

  Sum_protein_info(Protein_feature_list_rank,scoretype=scoretype)
  Protein_feature_list_rank<-as.data.frame(Protein_feature_list_rank)
  Protein_feature_list_rank<-merge(Protein_feature_list_rank,unique(sum_pro_pep_count[,c("Protein","peptide_count","isdecoy")]),by=c("Protein","isdecoy"))
  Protein_feature_list_rank<-Protein_feature_list_rank[Protein_feature_list_rank$Protein %in% sum_pro_pep_count$Protein[sum_pro_pep_count$peptide_count>=peptide_ID_filter],]
  Protein_feature_list_rank_filtered_grouped<-Protein_feature_list_rank

  #message(length(unique(Protein_feature_list_rank$Protein)))
  # Get the non-redundant protein grouping info
  if (protein_nr_grouping){
    protein_nr<-function(Protein_feature_list_rank){
      #suppressMessages(suppressWarnings(require(igraph)))
      suppressMessages(suppressWarnings(require(Biostrings)))
      suppressMessages(suppressWarnings(require(dplyr)))

      protein_ids=unique(Protein_feature_list_rank$Protein)

      peptide_links<-Protein_feature_list_rank  %>% group_by(.dots=c("Peptide")) %>% summarize(Protein=paste(unique(Protein),collapse = ","))

      nr_pro_pep<-unique(Protein_feature_list_rank[,c("Protein","Peptide")])

      protein_links<-nr_pro_pep %>% group_by(.dots=c("Protein")) %>% summarize(Peptide=list(Peptide))

      matched_proteins<-lapply(protein_ids,function(x,Protein_feature_list_rank,peptide_links){

        testvectors=unique(Protein_feature_list_rank$Peptide[Protein_feature_list_rank$Protein==x])

        matched_proteins<-peptide_links$Protein[peptide_links$Peptide %in% testvectors]

        matched_proteins<-strsplit(matched_proteins,",")

        if (length(matched_proteins)>1){

          matched_proteins_ol<-Reduce(intersect,matched_proteins)

        }else if (length(matched_proteins)==1){matched_proteins_ol<-matched_proteins[[1]]}

        matched_proteins_ol

      },Protein_feature_list_rank,peptide_links)

      names(matched_proteins)<-protein_ids

      protein_links$test_list_length<-unlist(lapply(protein_links$Peptide,length))

      matched_proteins_dup<-matched_proteins[duplicated(matched_proteins)]

      matched_proteins_final<-matched_proteins[!duplicated(matched_proteins)]

      message("Iterating Protein information...")

      matched_proteins_final_hP<-lapply(names(matched_proteins_final),function(x,protein_links){

        test_length<- protein_links$test_list_length[protein_links$Protein==as.numeric(x)]+1

        test_peptides<-protein_links$Peptide[protein_links$Protein==as.numeric(x)][[1]]

        protein_t_within<-protein_links[protein_links$test_list_length>=test_length,]
        length_test_peptides=length(test_peptides)
        findhigher=NA

        if (nrow(protein_t_within)>0){
          for (test_peptide_list in 1:length(protein_t_within$Peptide)){

            if (length(intersect(test_peptides,protein_t_within$Peptide[test_peptide_list][[1]]))==length_test_peptides){
              findhigher=protein_t_within$Protein[test_peptide_list]
              break
            }

          }

        }

        return(findhigher)

      },protein_links)

      matched_proteins_final_hP<-unlist(matched_proteins_final_hP)

      unique_protein<-names(matched_proteins_final[is.na(matched_proteins_final_hP)])

      return(unique_protein)

    }
    Protein_feature_list_rank_ID<-protein_nr(Protein_feature_list_rank)
    Protein_feature_list_rank<-Protein_feature_list_rank[Protein_feature_list_rank$Protein %in% Protein_feature_list_rank_ID,]
    Sum_protein_info(Protein_feature_list_rank,scoretype=scoretype)
  }


  # calculate protein score from final grouping info
  Sum_protein_info(Protein_feature_list_rank,scoretype=scoretype)

  #message(length(unique(Protein_feature_list_rank$Protein)))

  Protein_feature_list_rank<-Protein_feature_list_rank[Protein_feature_list_rank$Protein %in% sum_pro_pep_count$Protein[sum_pro_pep_count$peptide_count>=peptide_ID_filter],]
  Protein_feature_list_rank_filtered_grouped_final<-Protein_feature_list_rank

  query_protein=sum_pro_int[,c("Protein","isdecoy")]
  query_protein_list=  split(query_protein, seq(nrow(query_protein)))
  sum_pro_coverage <-unlist(lapply(query_protein_list,protein_coverage,Protein_feature_list_rank))
  sum_pro_coverage_df=data.frame(Protein=(sum_pro_int$Protein),isdecoy=sum_pro_int$isdecoy,Protein_coverage=sum_pro_coverage)
  sum_pro_pep_count<-merge(sum_pro_pep_count,sum_pro_coverage_df,by=c("Protein","isdecoy"),sort=F)
  Protein_feature_result<-merge(sum_pro_int,sum_pro_score,by=c("Protein","isdecoy"))
  Protein_feature_result<-merge(Protein_feature_result,sum_pro_pep_count,by=c("Protein","isdecoy"))
  Protein_feature_result$Intensity_norm<-log(Protein_feature_result$Intensity)/mean(log(Protein_feature_result$Intensity))
  Protein_feature_result$Intensity_norm<-(Protein_feature_result$Intensity_norm-min(Protein_feature_result$Intensity_norm))/max(Protein_feature_result$Intensity_norm)
  Protein_feature_result$Proscore=( Protein_feature_result$Score) * Protein_feature_result$Protein_coverage * Protein_feature_result$Intensity_norm
  Protein_feature_result$Protein=as.numeric(Protein_feature_result$Protein)
  Protein_feature_list_rank$Protein=as.numeric(Protein_feature_list_rank$Protein)
  Protein_feature_result=merge(Protein_feature_result,Index_of_protein_sequence[,c("recno","desc")],by.x="Protein",by.y="recno",all.x=T)
  Protein_feature_list_rank=merge(Protein_feature_list_rank,Index_of_protein_sequence[,c("recno","desc")],by.x="Protein",by.y="recno",all.x=T)
    # mask the mz featuers while decoy ID get higher scores, not recommended while no On-tissue MS/MS fragmentation is acquired.
  if(compete_decoy==T){
    protein_decoy_select<-Protein_feature_result %>% group_by(.dots=c("Protein")) %>% summarize(Proscore=max(Proscore))
    Protein_feature_result_decoy_compete<-merge(protein_decoy_select,Protein_feature_result,by=c("Protein","Proscore"))
    Protein_feature_result<-Protein_feature_result_decoy_compete
  }

  Protein_feature_result<-Protein_feature_result[!is.na(Protein_feature_result$Proscore),]
  #message(length(unique(Protein_feature_list_rank$Protein)))
  return(list(Protein_feature_result,
              Protein_feature_list_rank,
              Protein_feature_list_rank_filtered,
              Protein_feature_list_rank_filtered_grouped,
              Protein_feature_list_rank_filtered_grouped_final))
}


seq.cov <- function(x,plot_coverage=F){
  suppressMessages(suppressWarnings(require(IRanges)))
  if(!class(x)=="PairwiseAlignmentsSingleSubject"){
    message('Only supports objects from class PairwiseAlignedFixedSubject')
    stop()
  }

  coverage_ranges<-x@subject@range
  cov <- coverage(coverage_ranges)
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


#' export_pixel_level_data
#'
#' This is a peptide mass fingerprint search function for maldi imaging data analysis
#' @param datafile the data files' path for the analysis, leave it as blank to enable a graphical user interface to select the data
#' @param projectfolder optional, if NULL script will extract the path from datafile(s), and use the first workdir as project folder
#' @param ppm the mz tolerance (in ppm) for peak integration
#' @param Rotate_IMG specify a configuration file to further change the rotation of the images
#' @param Protein_desc_of_interest Specify a list of protein descriptions for cluster image plotting. Default setting will plot all reported proteins.
#' @param Protein_desc_of_exclusion Specify a list of protein descriptions to be excluded from cluster image plotting.
#' 
#' @return None
#'
#' @examples
#'
#' @export
#'
projectfolder="G:\\Documents\\expdata\\MouseBrain_Trypsin_FT\\"
Protein_peptide_file="/Summary folder/Protein_peptide_Summary.csv"
datafile=c("Mouse_brain_trimmed")
ppm=10
Rotate_IMG=NULL
Protein_desc_of_interest="."
Protein_desc_of_exclusion=NULL
Thread=4
export_pixel_level_data<-function(projectfolder=NULL,Protein_peptide_file,
                                  datafile=c(),
                                  ppm=5,
                                  Rotate_IMG=NULL,
                                  Protein_desc_of_interest=".",
                                  Protein_desc_of_exclusion=NULL,
                                  Thread=4,
                                  ...
){
  suppressMessages(suppressWarnings(library("pacman")))
  suppressMessages(suppressWarnings(p_load(stringr,BiocParallel,data.table,Cardinal,parallel)))
  
  if (missing(datafile)) stop("Missing data file, Choose single or multiple imzml file(s) for analysis")
  
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
    BPPARAM=HiTMaP:::Parallel.OS(parallel)
    setCardinalBPPARAM(BPPARAM = BPPARAM)
  }else{
    parallel=Thread
    BPPARAM=HiTMaP:::Parallel.OS(parallel)
    setCardinalBPPARAM(BPPARAM = BPPARAM)
  }
  
  message(paste(try(detectCores()), "Cores detected,",parallel, "threads will be used for computing"))
  
  message(paste(length(datafile), "files were selected and will be used for data extraction"))
  
  if(is.null(Protein_peptide_file)) {
    message("Missing proteomics data file, will use ID file in each folder respectively")
    Protein_peptide_file=paste0("/",datafile," ID/Peptide_region_file.csv")
  }
  
  if(file.exists(paste0(workdir,Protein_peptide_file))){
    read.csv()
  }
  
  
}

