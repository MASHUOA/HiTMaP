MonoisotopicMass_sil<-function (formula = list(), isotopes = list(), charge = 0, Isotype="Normal",c_num=6) 
{
  defaultFormula <- list(C = 0, H = 0, N = 0, O = 0, S = 0, 
                         P = 0, Br = 0, Cl = 0, F = 0, Si = 0, x = 0,'[13]C'=0)
  defaultFormula[names(formula)] <- formula
  defaultIsotopes <- list(C = 12, H = 1.0078250321, N = 14.0030740052, 
                          O = 15.9949146221, S = 31.97207069, P = 30.97376151, 
                          Br = 78.9183376, Cl = 34.96885271, F = 18.9984032, Si = 27.9769265327, 
                          x = 0,'[13]C'=13.003355)
  
  
  defaultIsotopes_13c <- list(C = 13.003355, H = 1.0078250321, N = 14.0030740052, 
                          O = 15.9949146221, S = 31.97207069, P = 30.97376151, 
                          Br = 78.9183376, Cl = 34.96885271, F = 18.9984032, Si = 27.9769265327, 
                          x = 0,'[13]C'=13.003355)
  defaultIsotopes[names(isotopes)] <- isotopes
  if(Isotype=="Normal"){
      if (charge < 0 & abs(charge) > defaultFormula$H) 
    stop("the number of negative charges exceeds the number of hydrogens in the formula list")
  
  mass <- (defaultFormula$C * defaultIsotopes$C + defaultFormula$H * 
             defaultIsotopes$H + defaultFormula$N * defaultIsotopes$N + 
             defaultFormula$O * defaultIsotopes$O + defaultFormula$S * 
             defaultIsotopes$S + defaultFormula$P * defaultIsotopes$P + 
             defaultFormula$Br * defaultIsotopes$Br + defaultFormula$Cl * 
             defaultIsotopes$Cl + defaultFormula$F * defaultIsotopes$F + 
             defaultFormula$Si * defaultIsotopes$Si + defaultFormula$x * 
             defaultIsotopes$x + defaultIsotopes$'[13]C' * defaultFormula$'[13]C')
  if (charge != 0) 
    mass <- abs((mass + charge * 1.007276466)/charge)
  }else if(Isotype=="SIL" ){
    if (charge < 0 & abs(charge) > defaultFormula$H) 
      stop("the number of negative charges exceeds the number of hydrogens in the formula list")
    if (defaultFormula$C>=6){
          mass <- (6 * defaultIsotopes_13c$C + (defaultFormula$C-6) * defaultIsotopes$C + defaultFormula$H * 
               defaultIsotopes_13c$H + defaultFormula$N * defaultIsotopes_13c$N + 
               defaultFormula$O * defaultIsotopes_13c$O + defaultFormula$S * 
               defaultIsotopes_13c$S + defaultFormula$P * defaultIsotopes_13c$P + 
               defaultFormula$Br * defaultIsotopes_13c$Br + defaultFormula$Cl * 
               defaultIsotopes_13c$Cl + defaultFormula$F * defaultIsotopes_13c$F + 
               defaultFormula$Si * defaultIsotopes_13c$Si + defaultFormula$x * 
               defaultIsotopes_13c$x)
    }else{
      mass <- (defaultFormula$C * defaultIsotopes$C + defaultFormula$H * 
                 defaultIsotopes$H + defaultFormula$N * defaultIsotopes$N + 
                 defaultFormula$O * defaultIsotopes$O + defaultFormula$S * 
                 defaultIsotopes$S + defaultFormula$P * defaultIsotopes$P + 
                 defaultFormula$Br * defaultIsotopes$Br + defaultFormula$Cl * 
                 defaultIsotopes$Cl + defaultFormula$F * defaultIsotopes$F + 
                 defaultFormula$Si * defaultIsotopes$Si + defaultFormula$x * 
                 defaultIsotopes$x)
    }

    if (charge != 0) 
      mass <- abs((mass + charge * 1.007276466)/charge)
  }

  return(mass)
}


Cpd_spectrum_match_rescore<-function(cpd_list,peaklist,wd=getwd(),
                                     tolerance_ppm=8,
                                     SIL_cpd=T,
                                     C13_SIL_incorporation=0.995,    
                                     SPECTRUM_batch=1,
                                     plot_matching_score_t=T,    
                                     adducts=c("M-H","M+Cl"),
                                     BPPARAM=bpparam(),
                                     atom_isotope_sub=NULL,
                                     ppm=10,
                                     score_method="SQRTP",
                                     adjust_score=F,
                                     FDR_cutoff=0.05,
                                     Decoy_search=T,
                                     Decoy_mode="isotope"){
  

    
    workdir<-wd
    adductslist<-Build_adduct_list()
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(stringr)))
  suppressMessages(suppressWarnings(require("Rcpp")))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(Rdisop)))
  suppressMessages(suppressWarnings(require(Biostrings)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(BiocParallel)))
  suppressMessages(suppressWarnings(require(data.table)))
  data(isotopes)
  cpd_list$formula_mono<-cpd_list$formula
  cpd_list$Matched.Form->cpd_list$adduct
  peaklist$intensities<-as.numeric(peaklist$intensities)
  peaklist <- peaklist[peaklist$intensities>0,]
   
  if(is.null(cpd_list$Isotype)) cpd_list$Isotype<-"Normal"
  if (SIL_cpd){
    cpd_list$Isotype[str_detect(cpd_list$Matched.Form,"M_13C")]<-"SIL"
    cpd_list$adduct<-str_replace(cpd_list$adduct,".M_13C..","M")
    #cpd_list$formula_mono[cpd_list$Isotype=="SIL"]<-str_replace(cpd_list$formula_mono[cpd_list$Isotype=="SIL"],"C","\\[13\\]C")
   #if(is.na(C13_SIL_incorporation) ){ stop("Please specify the incorporation rate of C(13)")}
   
   if (length(cpd_list$Isotype=="SIL")==0){
     cpd_list_sil<-cpd_list
     cpd_list$Isotype="SIL"
     cpd_list<-rbind(cpd_list,cpd_list_sil)
   }
   isotopes$abundance['&'(isotopes$isotope=="13C",isotopes$element=="[13]C")]=C13_SIL_incorporation
   isotopes['&'(isotopes$isotope=="13C",isotopes$element=="[13]C"),]->addrow
   addrow$isotope[1]="12C"
   addrow$mass[1]=12
   addrow$abundance[1]=1-C13_SIL_incorporation
   isotopes<-rbind(isotopes,addrow)
   rownames(isotopes)<-1:nrow(isotopes)
  }
  
  cpd_list<-cpd_list[cpd_list$formula!="",]
  cpd_list<-cpd_list[!grepl(")n",cpd_list$formula),]
  cpd_list->database->candidates
  
  required_col=c("formula","Isotype")
  candidates=data_test_rename(required_col,candidates)
  candidates$formula_isotype<-paste0(candidates$formula_mono,"_",candidates$Isotype)
  unique_formula<-as.data.frame(unique(candidates[,c("formula_mono","Isotype")]))
  unique_formula_list<-lapply(unique_formula$formula,get_atoms)
  names(unique_formula_list)<-paste0(unique_formula$formula,"_",unique_formula$Isotype)
  


  masslist<-lapply(1:nrow(unique_formula),function(x,unique_formula,unique_formula_list){
    MonoisotopicMass_sil(unique_formula_list[[paste0(unique_formula$formula[x],"_",unique_formula$Isotype[x])]],Isotype=unique_formula$Isotype[x])->mass
    return(mass)
    },unique_formula,unique_formula_list)
  
  if (SIL_cpd){
  unique_formula_list<-lapply(names(unique_formula_list),function(x){
    unique_formula_list[[x]]->y
    if(str_detect(x,"_SIL$")){
      if(y$C<6){
        y<-NULL
      }else{
        y$C<-y$C-6
        y$'[13]C'=6
      }
    }  
    return(y)
  })
  names(unique_formula_list)<-paste0(unique_formula$formula,"_",unique_formula$Isotype)
  sil_valid_idx<-which(!unlist(lapply(unique_formula_list,is.null)))
  unique_formula_list<-unique_formula_list[sil_valid_idx]
  unique_formula<-unique_formula[sil_valid_idx,]
  masslist<-masslist[sil_valid_idx]
  }
  
  unique_formula_list_reduce<-atomes.to.formula(unique_formula_list)
  
  uniquemass<-unlist(masslist)
  mass_DF<-data.frame(formula=unique_formula$formula,Isotype=unique_formula$Isotype,mass=uniquemass,stringsAsFactors = F)
  candidates<-base::merge(candidates,mass_DF,by=c("formula","Isotype"))
  candidates$formula<-unique_formula_list_reduce[match(paste0(candidates$formula,"_",candidates$Isotype),names(unique_formula_list))]
  
  candidates<-candidates[duplicated(names(candidates))==FALSE]
  
  if (is.null(candidates$Metabolite.Name)){candidates$Metabolite.Name<-candidates$Matched.Compound}

  
  required_col=c("mz","Metabolite.Name","mass","formula","Isotype")
  candidates=data_test_rename(required_col,candidates) 
  candidates$formula<-as.character(candidates$formula)
  meta_symbol<-unique_formula_list
  
  symbol_adducts=bplapply(adducts,convert_peptide_adduct_list,meta_symbol,BPPARAM = BPPARAM,adductslist=adductslist)
  
  
  symbol_adducts_df=lapply(symbol_adducts,
                           function(x,names_x){
                             data.frame(formula_isotype=names_x,formula_adduct=x)
                           },names(meta_symbol))
  names(symbol_adducts_df)<-adducts
  
    Meta_Summary<-NULL
  for (i in 1:length(adducts)){
    candidates_adduct<-unique(candidates[,c("Metabolite.Name","mass","formula_mono","Isotype")])
    candidates_adduct<-candidates_adduct[!duplicated(candidates_adduct[,c("Metabolite.Name","formula_mono","Isotype")]),]
    adductmass<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Mass"]))
    adducts_Charge<-as.numeric(adductslist[adductslist$Name==adducts[i],"Charge"])
    candidates_adduct$mz<-(as.numeric(candidates_adduct$mass+adductmass))/abs(adducts_Charge)
    candidates_adduct$adduct<-adducts[i]
    candidates_adduct$charge<-adducts_Charge
    candidates_adduct$formula_isotype<-paste0(candidates_adduct$formula_mono,"_",candidates_adduct$Isotype)
    candidates_adduct<-merge(candidates_adduct,symbol_adducts_df[[adducts[i]]])
    candidates_adduct$formula<-candidates_adduct$formula_adduct
    candidates_adduct$formula_adduct<-NULL
    Meta_Summary<-rbind.data.frame(Meta_Summary,unique(candidates_adduct))
  }
  
  #Meta_Summary$Isotype<-"Normal"
  #Meta_Summary<-Meta_Summary[!stringr::str_detect(Meta_Summary$formula,"-"),]
  # cpd_list$adduct <- str_remove_all(cpd_list$adduct,"\\[-\\]|\\[+\\]")
  # cpd_list$formula <- NULL
  # cpd_list$Matched.Compound->cpd_list$Metabolite.Name
  # cpd_list<-merge(cpd_list,Meta_Summary,by=c("Metabolite.Name","adduct"),all.x=T)
  cpd_list<-Meta_Summary
  
  #Do first round of peptide search to get putative result
  mz_feature_list<-Do_PMF_search(peaklist,cpd_list,BPPARAM=BPPARAM,ppm = ppm)
  mz_feature_list<-unique(mz_feature_list)
  message("Summarizing cpd information...",length(unique(mz_feature_list$mz))," were found in the first search.")
  cpd_list<-as.data.table(cpd_list)
  mz_feature_list<-as.data.table(mz_feature_list)
  cpd_list$Intensity<-NULL
  mz_feature_list$mz<-as.character(mz_feature_list$mz)
  
  cpd_list$mz<-as.character(cpd_list$mz)
  cpd_list->cpd_list_full
  cpd_list<-merge(cpd_list,mz_feature_list,by.x="mz",by.y="mz",all.x=T,sort=F)
  cpd_list$Intensity[is.na(cpd_list$Intensity)]<-0
  cpd_list<-cpd_list[cpd_list$Intensity>0,]
  cpd_list$formula<-as.character(cpd_list$formula)
  
  unique_formula<-unique(cpd_list[,c("formula","Isotype")])
  
  unique_formula<-unique_formula[!is.na(unique_formula$formula),]
  unique_formula<<-unique_formula
  
  cpd_list_Score=(lapply(unique_formula$formula,try(SCORE_CPD),peaklist=peaklist,isotopes=isotopes,score_method=score_method,charge = -1,ppm=ppm,threshold=1))
   cpd_list_Score<<-cpd_list_Score
  cpd_list_Score_m=as.data.frame(do.call(rbind, cpd_list_Score))
  names(cpd_list_Score_m)<-c("Score", "delta_ppm","Intensity")
  cpd_list_Score_m$Isotype<-unique_formula$Isotype
  formula_score<-data.frame(formula=unique_formula$formula,Score=cpd_list_Score_m$Score,Delta_ppm=cpd_list_Score_m$delta_ppm,Intensity=cpd_list_Score_m$Intensity,Isotype=unique_formula$Isotype)
  # if (SIL_cpd){
  # 
  #   cpd_list_Score=(lapply(unique_formula$formula[unique_formula$Isotype=="SIL"],SCORE_CPD,peaklist=peaklist,isotopes=SIL_isotopes,score_method=score_method,charge = -1,ppm=10))
  #   cpd_list_Score_m=as.data.frame(do.call(rbind, cpd_list_Score))
  #   #cpd_list_Score_m_sil<<-cpd_list_Score_m
  #   names(cpd_list_Score_m)<-c("Score", "delta_ppm","Intensity")
  #   cpd_list_Score_m$Isotype<-"SIL"
  #   formula_score_sil<-data.frame(formula=unique_formula$formula[unique_formula$Isotype=="SIL"],Score=cpd_list_Score_m$Score,Delta_ppm=cpd_list_Score_m$delta_ppm,Intensity=cpd_list_Score_m$Intensity)
  #   formula_score_sil$Isotype<-"SIL"
  #   formula_score<-rbind(formula_score,formula_score_sil)
  #   message("got",nrow(formula_score_sil),"Sil CPD scored.")
  # }
  
  cpd_list$Score<-NULL
  cpd_list$Delta_ppm<-NULL
  cpd_list$Intensity<-NULL
  cpd_list<-merge(cpd_list,formula_score,by=c("formula","Isotype"),all.x=T)
  cpd_list$isdecoy=0
  #generate decoy candidate list if running in "isotope" decoy mode 
  if (Decoy_search && ("isotope" %in% Decoy_mode)){
    decoy_isotopes=isotopes
    C_rows<-data.frame(element=c("C","C","[13]C","[13]C"),isotope=c("11C","12C","13C","14C"),mass=c(10.99664516,12,13.003355,14),abundance=c(0.0107,1-0.0107,1-C13_SIL_incorporation,C13_SIL_incorporation),ratioC=0,stringsAsFactors = F)
    decoy_isotopes[!(decoy_isotopes$element %in% c("[13]C","C")),]->decoy_isotopes
    decoy_isotopes<-rbind(decoy_isotopes,C_rows)
    #decoy_isotopes[decoy_isotopes$isotope %in% "13C",]=data.frame(element=c("C","[13]C","[13]C"),isotope=c("11C","13C","14C"),mass=c(10.99664516,13.003355,14.0067),aboudance=c(0.0107,0.0107,1-0.0107),ratioC=0,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "15N",]=data.frame(element="N",isotope="13N",mass=13.00603905,aboudance=0.00364000,ratioC=4,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "2H",]=data.frame(element="H",isotope="0H",mass=0.001548286,aboudance=0.00011500,ratioC=6,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "17O",]=data.frame(element="O",isotope="15O",mass=14.99069774,aboudance=0.00038000,ratioC=3,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "18O",]=data.frame(element="O",isotope="14O",mass=13.99066884,aboudance=0.00205000,ratioC=3,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "33S",]=data.frame(element="S",isotope="31S",mass=30.97268292,aboudance=0.0075,ratioC=3,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "34S",]=data.frame(element="S",isotope="30S",mass=29.9762745,aboudance=0.0425,ratioC=3,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "35S",]=data.frame(element="S",isotope="29S",mass=28.94414146,aboudance=0,ratioC=3,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "36S",]=data.frame(element="S",isotope="28S",mass=27.97706058,aboudance=0.0001,ratioC=3,stringsAsFactors = F)
    decoy_isotopes[decoy_isotopes$isotope %in% "37Cl",]=data.frame(element="Cl",isotope="33Cl",mass=32.965903,aboudance=0.24240000,ratioC=3,stringsAsFactors = F)
    
    cpd_list_decoy<-cpd_list[cpd_list$isdecoy==0,]
    cpd_list_decoy$isdecoy=1
    unique_formula<-unique(cpd_list_decoy[,c("formula","Isotype")])
    
    unique_formula<-unique_formula[!is.na(unique_formula$formula),]
    unique_formula<<-unique_formula
    
    cpd_list_Score=(lapply(unique_formula$formula,SCORE_CPD,peaklist=peaklist,isotopes=decoy_isotopes,score_method=score_method,charge = -1,ppm=ppm,threshold=1))
    cpd_list_Score<<-cpd_list_Score
    cpd_list_Score_m=as.data.frame(do.call(rbind, cpd_list_Score))
    names(cpd_list_Score_m)<-c("Score", "delta_ppm","Intensity")
    cpd_list_Score_m$Isotype<-unique_formula$Isotype
    formula_score<-data.frame(formula=unique_formula$formula,Score=cpd_list_Score_m$Score,Delta_ppm=cpd_list_Score_m$delta_ppm,Intensity=cpd_list_Score_m$Intensity,Isotype=unique_formula$Isotype)
    # if (SIL_cpd){
    # 
    #   cpd_list_Score=(lapply(unique_formula$formula[unique_formula$Isotype=="SIL"],SCORE_CPD,peaklist=peaklist,isotopes=SIL_isotopes,score_method=score_method,charge = -1,ppm=10))
    #   cpd_list_Score_m=as.data.frame(do.call(rbind, cpd_list_Score))
    #   #cpd_list_Score_m_sil<<-cpd_list_Score_m
    #   names(cpd_list_Score_m)<-c("Score", "delta_ppm","Intensity")
    #   cpd_list_Score_m$Isotype<-"SIL"
    #   formula_score_sil<-data.frame(formula=unique_formula$formula[unique_formula$Isotype=="SIL"],Score=cpd_list_Score_m$Score,Delta_ppm=cpd_list_Score_m$delta_ppm,Intensity=cpd_list_Score_m$Intensity)
    #   formula_score_sil$Isotype<-"SIL"
    #   formula_score<-rbind(formula_score,formula_score_sil)
    #   message("got",nrow(formula_score_sil),"Sil CPD scored.")
    # }
    cpd_list_decoy$Score<-NULL
    cpd_list_decoy$Delta_ppm<-NULL
    cpd_list_decoy$Intensity<-NULL
    cpd_list_decoy<-merge(cpd_list_decoy,formula_score,by=c("formula","Isotype"),all.x=T)
    cpd_list<-rbind(cpd_list,cpd_list_decoy)
  }
  cpd_list$mz->cpd_list$Query.Mass
  cpd_list$mz<-as.numeric(as.character(cpd_list$mz))
  #cpd_list<<-cpd_list
  if(all(unique(cpd_list$Score) %in% c(0,NA))){
    message("No matched feature, proceed to next spectrum")
    return(1)
    }
  cpd_list_rank=rank_mz_feature(cpd_list,mz_feature=peaklist,BPPARAM = BPPARAM)
  cpd_list_rank$mz_align<-unlist(cpd_list_rank$mz_align)
  
  cpd_list_rank<-as.data.frame(unique(cpd_list_rank))
  cpd_list_rank<-cpd_list_rank[!is.na(cpd_list_rank$Score),]
  

  
  dir.create(paste0(wd,"/",SPECTRUM_batch), recursive = T)
  cpd_list_rank$Delta_ppm[is.na(cpd_list_rank$Delta_ppm)]<-ppm
  cpd_list_rank<-cpd_list_rank[cpd_list_rank$Intensity>=0,]
  cpd_list_rank<-cpd_list_rank[!is.na(cpd_list_rank$Intensity),]
  
  write.csv(cpd_list_rank,paste0(wd,"/",SPECTRUM_batch,"/","CPD_1st_ID_score_rank_",score_method,".csv"),row.names = F)
  
  cpd_list_rank<<-cpd_list_rank
  #Peptide score filtering (default=False) 
  if(nrow(cpd_list_rank)>=2){
  if (adjust_score==F){
    dir.create(paste0(wd,"/",SPECTRUM_batch))
    Score_cutoff= FDR_cutoff_plot(cpd_list_rank,FDR_cutoff=FDR_cutoff,plot_fdr=T,outputdir=paste0(wd,"/",SPECTRUM_batch),adjust_score = adjust_score)
    cpd_list_2nd=cpd_list_rank[((cpd_list_rank$Score>=Score_cutoff)&(!is.na(cpd_list_rank$Intensity))),]
    if(SIL_cpd){
      dir.create(paste0(wd,"/",SPECTRUM_batch,"/","SIL"))
      cpd_list_rank[cpd_list_rank$Isotype=="SIL",]->cpd_list_2nd_sil
      Score_cutoff_sil= FDR_cutoff_plot(cpd_list_2nd_sil,FDR_cutoff=FDR_cutoff,plot_fdr=T,outputdir=paste0(wd,"/",SPECTRUM_batch,"/","SIL"),adjust_score = adjust_score)
      cpd_list_2nd_sil=cpd_list_2nd_sil[((cpd_list_2nd_sil$Score>=Score_cutoff_sil)&(!is.na(cpd_list_2nd_sil$Intensity))),]
      if (nrow(cpd_list_2nd_sil)>=1){
        cpd_list_2nd<-rbind(cpd_list_2nd,cpd_list_2nd_sil)
        cpd_list_2nd<-unique(cpd_list_2nd)
      }
    }
    }else{
      dir.create(paste0(wd,"/",SPECTRUM_batch))
      cpd_list_rank$Score_org<-cpd_list_rank$Score
      mz.score.correction <-function(mz.ref,Score.ref,mz,Score,
                               control = loess.control(),
                               span = 0.75,fit_type="loess",force=F) {
        suppressMessages(suppressWarnings(library(matter)))
        
        if(fit_type=="loess"){
          shift <-
            suppressWarnings(loess(Score.ref ~ mz.ref, span = span, control = control))
          dScore <- predict(shift, mz)
        }
        
        if(fit_type=="lm"){
          shift <-
            suppressWarnings(lm(Score.ref ~ mz.ref))
          
          dScore <- predict.lm(shift, data.frame(mz.ref=mz))
        }
        
        if(fit_type=="poly"){
          shift <-
            suppressWarnings(lm(Score.ref ~ poly(mz.ref,2)))
          
          dScore <- predict.lm(shift, data.frame(mz.ref=mz))
        }
        #warp <- splinefun(mz + dScore, x)
        return(list(
          final.Score = Score + dScore,
          dScore = dScore,
          shift = shift
        ))
        #return(list(final.mz = mz + dmz, dmz=dmz))
      }
      unique(cpd_list_rank[,c("formula","mz","Score_org","Isotype","isdecoy")])->unique_formula_score
      unique_formula_score_decoy<-unique_formula_score[unique_formula_score$isdecoy==1,]
      cpd_list_rank$Score_adj<-mz.score.correction(mz.ref = unique_formula_score_decoy$mz,
                          Score.ref = unique_formula_score_decoy$Score_org,
                          mz = cpd_list_rank$mz,
                          Score = cpd_list_rank$Score_org,
                          fit_type = "lm")$final.Score
      cpd_list_rank$Score<-cpd_list_rank$Score_adj
      Score_cutoff= FDR_cutoff_plot(cpd_list_rank,FDR_cutoff=FDR_cutoff,plot_fdr=T,outputdir=paste0(wd,"/",SPECTRUM_batch),adjust_score = F)
      cpd_list_2nd=cpd_list_rank[((cpd_list_rank$Score>=Score_cutoff)&(!is.na(cpd_list_rank$Intensity))),]
      if(SIL_cpd){
        dir.create(paste0(wd,"/",SPECTRUM_batch,"/","SIL"))
        cpd_list_rank[cpd_list_rank$Isotype=="SIL",]->cpd_list_2nd_sil
        unique(cpd_list_2nd_sil[,c("formula","mz","Score_org","Isotype","isdecoy")])->unique_formula_score
        unique_formula_score_decoy<-unique_formula_score[unique_formula_score$isdecoy==1,]
        cpd_list_2nd_sil$Score_adj<-mz.score.correction(mz.ref = unique_formula_score_decoy$mz,
                                                     Score.ref = unique_formula_score_decoy$Score_org,
                                                     mz = cpd_list_2nd_sil$mz,
                                                     Score = cpd_list_2nd_sil$Score_org,
                                                     fit_type = "loess")$final.Score
        cpd_list_2nd_sil$Score<-cpd_list_2nd_sil$Score_adj
        #plot(cpd_list_2nd_sil$mz,cpd_list_2nd_sil$Score_adj-cpd_list_2nd_sil$Score_org)
        Score_cutoff= FDR_cutoff_plot(cpd_list_2nd_sil,FDR_cutoff=FDR_cutoff,plot_fdr=T,outputdir=paste0(wd,"/",SPECTRUM_batch,"/","SIL"),adjust_score = F)
        cpd_list_2nd_sil=cpd_list_2nd_sil[((cpd_list_2nd_sil$Score>=Score_cutoff)&(!is.na(cpd_list_2nd_sil$Intensity))),]
        if (nrow(cpd_list_2nd_sil)>=1){
          cpd_list_2nd<-rbind(cpd_list_2nd,cpd_list_2nd_sil)
          cpd_list_2nd$Score_org->cpd_list_2nd$Score
          cpd_list_2nd$Score_org<-NULL
          cpd_list_2nd$Score_adj<-NULL
          cpd_list_2nd<-unique(cpd_list_2nd)
        }
      }
     }
    }            
  write.csv(cpd_list_2nd_sil,paste0(wd,"/",SPECTRUM_batch,"/","CPD_2nd_ID_score_rank_SIL_",score_method,".csv"),row.names = F)
  cpd_list_rank$mz=as.numeric(as.character(cpd_list_rank$mz))
  cpd_list_2nd$mz=as.numeric(as.character(cpd_list_2nd$mz))
  cpd_list_2nd<-cpd_list_2nd[cpd_list_2nd$isdecoy==0,]
  write.csv(cpd_list_2nd,paste0(wd,"/",SPECTRUM_batch,"/","CPD_2nd_ID_score_rank_",score_method,".csv"),row.names = F)

  #Plot peptide matching spectrum and score 
  if (plot_matching_score_t && nrow(cpd_list_2nd)!=0){
    dir.create(paste0(wd,"/",SPECTRUM_batch,"/cpd_spectrum_match"))
    try(plot_cpd_matching_score(cpd_list_2nd[cpd_list_2nd$Isotype=="Normal",],peaklist,isotopes = isotopes,charge=-1,ppm,outputdir=paste0(wd,"/",SPECTRUM_batch,"/cpd_spectrum_match")))
    if (SIL_cpd){
      dir.create(paste0(wd,"/",SPECTRUM_batch,"/cpd_spectrum_match_sil"))
    try(plot_cpd_matching_score(cpd_list_2nd[cpd_list_2nd$Isotype=="SIL",],peaklist,isotopes = isotopes,charge=-1,ppm,outputdir=paste0(wd,"/",SPECTRUM_batch,"/cpd_spectrum_match_sil")))
    }
    
    }


}

Plot_cpd_spectrum_match_sig<-function(pimresultindex,spectrumlist,peplist,pimlist,pimresultlist, threshold=0.05){
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

Plot_cpd_spectrum_match<-function(pimresultindex,spectrumlist,peplist,pimlist,pimresultlist, threshold=0.05){
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

remove_cpd_score_outlier<-function(SMPLIST,IQR_LB=0.75,outputdir=getwd(),abs_cutoff=-2){
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
  

FDR_cutoff_plot_cpd<-function(cpd_plot_list,FDR_cutoff=0.1,FDR_strip=500,plot_fdr=F,outputdir=NULL,adjust_score=F){
  suppressMessages(suppressWarnings(require(ggplot2)))
  suppressMessages(suppressWarnings(require(data.table)))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(zoo)))
  #suppressMessages(suppressWarnings(require(FTICRMS)))
  cpd_plot_list=data_test_rename(c("isdecoy","Score"),cpd_plot_list)
  
  cpd_plot_list$isdecoy<-factor(cpd_plot_list$isdecoy)
  
  #unique_fomula<-unique(cpd_plot_list$formula)
  cpd_plot_list_output<-cpd_plot_list
  
  #unique_fomula_ID<-unique(cpd_plot_list[,c("mz","formula","isdecoy","Intensity","Score")])
  
  plot_fdr_histogram<-function(cpd_plot_list,plot_name="target-decoy",outputdir=outputdir){
    suppressMessages(suppressWarnings(require(dplyr)))
    if (nrow(cpd_plot_list)!=0) {
      cpd_plot_list_plot<-cpd_plot_list[,c("isdecoy","Score")]
      target_decoy<-factor(ifelse(cpd_plot_list_plot$isdecoy==0,"Target","Decoy"),levels = c("Target","Decoy"))
      cpd_plot_list_plot$target_decoy=target_decoy
      cpd_plot_list$target_decoy=target_decoy
      if(F){ png(paste0(outputdir,"/Peptide_Score_histogram_",plot_name,".png"))
        p<-ggplot(data=cpd_plot_list_plot,aes(x=cpd_plot_list_plot$Score,color=target_decoy, fill=target_decoy)) +
          geom_histogram( fill="white",alpha=0.5, bins = 50, position="Dodge")  +
          ggtitle("Peptide score vs Counts") +
          xlab("Score") + ylab("Counts") + labs(fill = "Is_Decoy") + theme_classic() #+ facet_grid(target_decoy ~ .)
        print(p)
        dev.off() }
      
      mu <- cpd_plot_list_plot %>% group_by(target_decoy) %>% dplyr::summarize(mean=mean(Score))
      
      png(paste0(outputdir,"/Peptide_Score_histogram_",plot_name,".png"))
      p<-ggplot(cpd_plot_list_plot, aes(x=cpd_plot_list_plot$Score, color=target_decoy, fill=target_decoy)) +
        geom_histogram(aes(y=..density..), position="Dodge", alpha=0.5, bins = 50)+
        geom_density(alpha=0.6)+ xlim(quantile(cpd_plot_list_plot$Score, c(0.005,0.995)))+
        geom_vline(data=mu, aes(xintercept=mu$mean, color=mu$target_decoy), linetype="dashed")+
        ggtitle("Peptide score vs Counts") +
        xlab("Score") + ylab("Counts") +
        #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
        theme_classic()
      print(p)
      dev.off()
      
      
      if (log(nrow(cpd_plot_list_plot))>0){
        png(paste0(outputdir,"/Matching_Score_vs_mz_",plot_name,".png"))
        alpha=1/20/(round(log(nrow(cpd_plot_list_plot),base = 2))/40)
        if (alpha>=1 ) alpha=1
        p<-ggplot(data=cpd_plot_list,aes(x=mz,y=Score,colour=target_decoy)) + geom_point(size=1,alpha=alpha,aes(colour = target_decoy)) +
          ggtitle("Matching score vs mz") +
          xlab("mz") + ylab("Matching score") + labs(fill = "isdecoy")+ theme(axis.text.x = element_text(angle=45))
        print(p)
        dev.off()
        
      }
      
    }
  }
  #cpd_plot_list=unique_fomula_ID
  
  if (!is.null(outputdir) && plot_fdr){
    
    plot_fdr_histogram(cpd_plot_list,outputdir=outputdir)
    
  }
  
  
  if(F){
    
    score_boundary<-sapply(1:max(cpd_plot_list$mz),function(x,cpd_plot_list){
      
      max(cpd_plot_list$Score[between(cpd_plot_list$mz,x-0.5,x+0.5)])
      
    },cpd_plot_list)
    
    names(score_boundary)=as.character(1:max(cpd_plot_list$mz))
    
    score_boundary<-score_boundary[score_boundary>0]
    
    score_boundary<-score_boundary[!is.na(score_boundary)]
    
    score_boundary<-data.frame(mz=names(score_boundary),top=score_boundary,stringsAsFactors = F)
    
    x <- as.numeric(score_boundary$mz)
    y <- score_boundary$top
    #plot(x,y)
    
    fit3 <- lm( y~poly(as.numeric(as.character(x)),3) )
    
    fitx=cpd_plot_list$mz
    
    fitdata=as.data.frame(stats::predict(fit3, data.frame(x = fitx), interval="confidence"))
    
    cpd_plot_list$score_factor<-fitdata$upr
    
    cpd_plot_list$adjusted_score<-cpd_plot_list$Score/cpd_plot_list$score_factor
    cpd_plot_list$adjusted_score<-cpd_plot_list$adjusted_score/max(cpd_plot_list$adjusted_score,na.rm = T)
    cpd_plot_list$original_score<-cpd_plot_list$Score
    cpd_plot_list$Score<-cpd_plot_list$adjusted_score
    if (!is.null(outputdir) && plot_fdr){
      png(paste0(outputdir,"/Matching_Score_vs_mz_after_adjustment.png"))
      p<-ggplot(data=cpd_plot_list,aes(x=mz,y=Score,color=isdecoy)) + geom_point(size=1,alpha=1/10) +
        ggtitle("Matching score vs mz after adjustment") +
        xlab("mz") + ylab("Matching score") + labs(fill = "isdecoy")
      print(p)
      dev.off()
    }
    
  }
  
  if(adjust_score){
    
    cpd_plot_list_sep<-cpd_plot_list_output[cpd_plot_list_output$isdecoy==0,]
    cpd_plot_list_sep_d<-cpd_plot_list_output[cpd_plot_list_output$isdecoy==1,]
    cpd_plot_list_mer<-merge(cpd_plot_list_sep,cpd_plot_list_sep_d,by=names(cpd_plot_list_sep)[!names(cpd_plot_list_sep) %in% c("isdecoy","Score" ,"Rank")])
    cpd_plot_list_mer$Score=cpd_plot_list_mer$Score.x-cpd_plot_list_mer$Score.y
    cpd_plot_list_mer$isdecoy=0
    cpd_plot_list_adj<-cpd_plot_list_mer
    cpd_plot_list_mer$Score=0
    cpd_plot_list_mer$isdecoy=1
    cpd_plot_list_adj<-rbind(cpd_plot_list_adj,cpd_plot_list_mer)
    cpd_plot_list_adj$isdecoy<-as.factor(cpd_plot_list_adj$isdecoy)
    plot_fdr_histogram(cpd_plot_list_adj,plot_name = "adj",outputdir=outputdir)
    cpd_plot_list<-cpd_plot_list_adj
    
    
    plot(density(cpd_plot_list_adj$Score[cpd_plot_list_adj$isdecoy==0]))
    den<-density(cpd_plot_list_adj$Score[cpd_plot_list_adj$isdecoy==0])
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
    
    mm<-resolve_multi_modes(cpd_plot_list_adj$Score[cpd_plot_list_adj$isdecoy==0],adjust = 2)
    
    unique_fomula_ID<-unique(cpd_plot_list[,c("mz","formula","isdecoy","Intensity","Score")])
    target_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy!=1]
    decoy_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy==1]
    target_score[is.na(target_score)]=0
    decoy_score[is.na(decoy_score)]=0
    target_score[target_score==Inf]<-0
    decoy_score[decoy_score==Inf]<-0
    
    breaks = seq( min(cpd_plot_list$Score), max(cpd_plot_list$Score) + (max(cpd_plot_list$Score)-min(cpd_plot_list$Score))/FDR_strip , by=(max(cpd_plot_list$Score)-min(cpd_plot_list$Score))/FDR_strip)
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
    
    unique_fomula_ID<-unique(cpd_plot_list[,c("mz","formula","isdecoy","Intensity","Score")])
    target_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy!=1]
    decoy_score<-unique_fomula_ID$Score[unique_fomula_ID$isdecoy==1]
    target_score[is.na(target_score)]=0
    decoy_score[is.na(decoy_score)]=0
    target_score[target_score==Inf]<-0
    decoy_score[decoy_score==Inf]<-0
    if (length(unique(cpd_plot_list$Score))==1 && !(cpd_plot_list$isdecoy %in% 1)) {
      return(unique(cpd_plot_list$Score))
    }else if (length(unique(cpd_plot_list$Score))==1 && (cpd_plot_list$isdecoy %in% 1)){
      return(unique(cpd_plot_list$Score)+100)
    }
    if (length(unique(cpd_plot_list$Score))==1) return(unique(cpd_plot_list$Score))
    breaks = seq( min(cpd_plot_list$Score), max(cpd_plot_list$Score) + (max(cpd_plot_list$Score)-min(cpd_plot_list$Score))/FDR_strip , by=abs(max(cpd_plot_list$Score)-min(cpd_plot_list$Score))/FDR_strip)
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
    return(list(Score_cutoff,cpd_plot_list))
  }else{
    return(Score_cutoff)
  }
  
  
  
  
}

plot_matching_score_cpd<-function(cpd_plot_list,peaklist,charge,ppm,outputdir=getwd(),filename_col=c("formula","adduct"),anno_col=filename_col,isotopes="default",similarity_plot_res = 72){
  if (ppm>=25) {
    instrument_ppm=50
  }else{
    instrument_ppm=8
  }
  message("plot matching isotopic pattern")
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(dplyr)))
  cpd_plot_list=data_test_rename(unique(c("formula","isdecoy","adduct",filename_col,anno_col)),cpd_plot_list)
  if(isotopes=="default") data("isotopes")
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
  #cpd_plot_list1=cpd_plot_list
  cpd_plot_list=cpd_plot_list[!is.na(cpd_plot_list$Intensity),]
  cpd_plot_list=cpd_plot_list %>% group_by(.dots = unique(c("formula","isdecoy","adduct",filename_col,anno_col))) %>% summarise(Peptide=paste(Peptide,collapse = "_"))
  
  if(nrow(cpd_plot_list)!=0){
    cpd_plot_list_score<-NULL
    for (i in 1:nrow(cpd_plot_list)){
      if ("Isotype" %in% anno_col){
        anno_info=paste(cpd_plot_list[i,anno_col[anno_col!="Isotype"]],collapse = "\n")
        isolabel=cpd_plot_list[i,"Isotype"]
      }else{
        anno_info=paste(cpd_plot_list[i,anno_col],collapse = "\n") 
        isolabel="Normal"
      }
      
      #anno_info=stringr::str_replace_all(anno_info,";","\n")
      if(cpd_plot_list$isdecoy[i]==0){
        cpd_plot_list_score[i]<-try(SCORE_CPD(cpd_plot_list$formula[i],
                                                  peaklist,isotopes=isotopes,
                                                  charge=charge,ppm=ppm,isolabel=isolabel,
                                                  print.graphic=T,output.list=F,
                                                  outputfile=paste0(outputdir,"/",paste(cpd_plot_list[i,filename_col],collapse = "_"),"_target",".png"),
                                                  similarity_plot_res = similarity_plot_res,
                                                  anno_info=anno_info))
        
        
      }else{
        cpd_plot_list_score[i]<-try(SCORE_CPD(cpd_plot_list$formula[i],
                                                  peaklist,isotopes=decoy_isotopes,
                                                  charge=charge,ppm=ppm,isolabel=isolabel,
                                                  print.graphic=T,output.list=F,
                                                  outputfile=paste0(outputdir,"/",
                                                                    paste(cpd_plot_list[i,filename_col],collapse = "_"),"_decoy",".png"),
                                                  similarity_plot_res = similarity_plot_res,
                                                  anno_info=anno_info))
      }
    }
    cpd_plot_list$Score<-cpd_plot_list_score
    return(cpd_plot_list)
  }
  
  
  
}

annotate_PSEA<-function(mSet,mummichog.lib.new){
  lapply(mSet[["path.hits"]],function(x){
    data.frame(cpd_id=paste0(x,collapse = ";"),cpd_names=paste0(mummichog.lib.new$cpd.lib$name[match(x,mummichog.lib.new$cpd.lib$id)],collapse = ";"))
  })->res
  res<-do.call(rbind,res)
  cbind(mSet[["mummi.resmat"]],res)->resdf
  return(resdf)
}

annotate_cpd<-function(mummichog.lib.new){
  cpd_list<-read.csv("mummichog_matched_compound_all.csv")
  cpd_list$Name=mummichog.lib.new$cpd.lib$name[match(cpd_list$Matched.Compound,mummichog.lib.new$cpd.lib$id)]
  return(cpd_list)
}

mummichog.lib.mod<-function(lib="bta_kegg",lib.new="bta_kegg_13C",C13_number=c(3,6),C13_deltamass=1.003355,wd=getwd(),force_update_lib=F,adducts_list="all",method=c("new_cpd","new_adduct")){
  
  filenm <- paste(wd,"/",lib, ".qs", sep="")
  mum.url <- paste("https://www.metaboanalyst.ca/resources/libs/mummichog/", paste(lib, ".qs", sep=""), sep="")
  if(`|`(force_update_lib,!file.exists(filenm))){
    download.file(mum.url, destfile = filenm, method="libcurl", mode = "wb")
  }
  
  mummichog.lib <- qs::qread(filenm)
  cmpd.map <- MetaboAnalystR:::.get.my.lib("compound_db.qs")
  hit.inx <- match(tolower(mummichog.lib$cpd.lib$id), tolower(cmpd.map$kegg))
  nms <- cmpd.map[hit.inx, "smiles"]
  c_count <-stringr::str_count(nms,"C|c")
  nms_pub <- cmpd.map[hit.inx, "pubchem_id"]
  
  df<-data.frame(smile=nms,c_count,mw_c_ratio=mummichog.lib$cpd.lib$mw/c_count,id=mummichog.lib$cpd.lib$id,mw=mummichog.lib$cpd.lib$mw,name=mummichog.lib$cpd.lib$name)
  
  mummichog.lib.mod<-NULL
  mummichog.lib.new<-mummichog.lib
  
  if (method[1]=="new_cpd"){
    for (cnum in C13_number){
      (df$c_count>=cnum)->mod.item
      (df$mw>=26.09864*as.numeric(cnum))->mod.item2
      mod.item[is.na(mod.item)]<-F
      mod.item2[is.na(mod.item2)]<-F
      mod.item_final<-as.numeric(ifelse(mod.item,mod.item,mod.item2))
      mummichog.lib.mod[[cnum]]<-mummichog.lib
      lapply(mummichog.lib.mod[[cnum]][["pathways"]][["cpds"]],function(x,cnum){
        paste0(x,"_13C",cnum)
      },cnum)->mummichog.lib.mod[[cnum]][["pathways"]][["cpds"]]
      mummichog.lib.mod[[cnum]][["cpd.lib"]][["id"]]<-paste0(mummichog.lib.mod[[cnum]][["cpd.lib"]][["id"]],"_13C",cnum)
      mummichog.lib.mod[[cnum]][["cpd.lib"]][["name"]]<-paste0(mummichog.lib.mod[[cnum]][["cpd.lib"]][["name"]],"_13C",cnum)
      mummichog.lib.mod[[cnum]][["cpd.lib"]][["mw"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["mw"]] + C13_deltamass*as.numeric(cnum)*mod.item_final
      mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["dpj_positive"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["dpj_positive"]]+ C13_deltamass*as.numeric(cnum)*mod.item_final
      mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["positive"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["positive"]]+ C13_deltamass*as.numeric(cnum)*mod.item_final
      mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["negative"]]<-mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]][["negative"]]+ C13_deltamass*as.numeric(cnum)*mod.item_final
      
    }
    for (cnum in C13_number){
      
      lapply(1:length(mummichog.lib.new$pathways$cpds),function(x){
        c(mummichog.lib.new$pathways$cpds[[x]],mummichog.lib.mod[[cnum]]$pathways$cpds[[x]])
      })->mummichog.lib.new$pathways$cpds
      
      mummichog.lib.new$cpd.lib$id<-c(mummichog.lib.new$cpd.lib$id,mummichog.lib.mod[[cnum]]$cpd.lib$id)
      mummichog.lib.new$cpd.lib$name<-c(mummichog.lib.new$cpd.lib$name,mummichog.lib.mod[[cnum]]$cpd.lib$name)
      mummichog.lib.new$cpd.lib$mw<-c(mummichog.lib.new$cpd.lib$mw,mummichog.lib.mod[[cnum]]$cpd.lib$mw)
      mummichog.lib.new[["cpd.lib"]][["adducts"]]$dpj_positive<-rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$dpj_positive,mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$dpj_positive)
      mummichog.lib.new[["cpd.lib"]][["adducts"]]$positive<-rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$positive,mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$positive)
      mummichog.lib.new[["cpd.lib"]][["adducts"]]$negative<-rbind(mummichog.lib.new[["cpd.lib"]][["adducts"]]$negative,mummichog.lib.mod[[cnum]][["cpd.lib"]][["adducts"]]$negative)
      
    }
  }
  
  
  
  
  
  #MetaboAnalystR:::CreateLibFromKEGG(mummichog.lib.new$cpd.lib, mummichog.lib.new$pathways, lib.new)
  cpd.lib <- mummichog.lib.new$cpd.lib;
  ms_modes <- c('dpj_positive', 'positive', 'negative');
  adducts <- list();
  for (ms_mode in ms_modes){
    if(adducts_list[1]=="all"){
      adducts[[ms_mode]] <- MetaboAnalystR:::Compound_function_mzlist(ms_mode, cpd.lib$mw)
    }else{
      resdf <- MetaboAnalystR:::Compound_function_mzlist(ms_mode, cpd.lib$mw)
      if(sum(adducts_list %in% colnames(resdf))!=0){
        adducts[[ms_mode]]<-resdf[,adducts_list[adducts_list %in% colnames(resdf)]]
      }else{
        adducts[[ms_mode]]<-resdf
      }
    }
  }
  
  if (method[1]=="new_adduct"){
    
    ms_modes <- c('dpj_positive', 'positive', 'negative')
    for (ms_mode in ms_modes){
      new_adduct<-NULL
      for (cnum in (C13_number)){
        new_adduct[[cnum]]<-adducts[[ms_mode]]+C13_deltamass*as.numeric(cnum)
        colnames(new_adduct[[cnum]])<-stringr::str_replace(colnames(new_adduct[[cnum]]),"^M",paste0("[M_13C",cnum,"]"))
      }
      new_adduct<-do.call(cbind,new_adduct)
      adducts[[ms_mode]]<-cbind(adducts[[ms_mode]],new_adduct)
    }
  }
  cpd.lib$adducts <- adducts;
  
  # create a dictionary for look up in the range of 50-2000
  # now need to create ladder (tree) for each new mz
  # key is the mass 50 to 2000, values are the compounds (if any of their modified mw gives the value)
  # now create cpd tree for each mass pos
  # note, this can be slow, but this can be created before hand
  # for each species and for each mode
  # note l2 only stores the index of the cpd.lib
  
  cpd.tree <- list();
  for (ms_mode in ms_modes){
    l2 <- list();
    l2[[49]] <- "";
    l2[[2001]] <- "";
    mz.mat <- cpd.lib$adducts[[ms_mode]];
    floor.mzs <- floor(mz.mat);
    for(i in 1:nrow(floor.mzs)){
      neighbourhood <- floor.mzs[i,];
      for(n in neighbourhood){
        if((n>50) & (n<2000)){
          l2[[n]] <- append(l2[[n]], i);
        }
      }
    }
    cpd.tree[[ms_mode]] <- lapply(l2, unique);
  }
  
  # set up the variables
  mummichog.lib.new <- list(
    pathways = mummichog.lib.new$pathways,
    cpd.tree = cpd.tree,
    cpd.lib = cpd.lib
  )
  qs::qsave(mummichog.lib.new,paste0(wd,"/",lib.new,".qs"))
  return(mummichog.lib.new)
}

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
    
  }else if (score_method=="SQRTC"){
    u <- alignment[,2]; v <- alignment[,3]
    similarity_score <- -log(as.vector(sqrt(sum(((u-v))^2) / (sum(u^2) * sum(v^2)))))
    match_score <- (sum(`&`(alignment$intensity.top>0,alignment$intensity.bottom>0))/sum(alignment$intensity.top>0))
    similarity_score <- similarity_score * match_score
    
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
  if (length(pattern_ppm)>=2){
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
  }
  matrix(filtered_pattern[,1:2],ncol=2)->filtered_pattern
  rownames(filtered_pattern)=1:nrow(filtered_pattern)
  
  return(filtered_pattern)
}

SCORE_CPD<-function(formula,peaklist,isotopes=NULL,threshold=0.001,charge=1,ppm=5,print.graphic=F,output.list=F,outputfile=NULL,score_method="SQRTP",similarity_plot_res=72,anno_info="",isolabel="Normal",output_monomz=F){
  suppressMessages(suppressWarnings(require(rcdk)))
  suppressMessages(suppressWarnings(require(rcdklibs)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(data.table)))
 require_rjava("Mummichog SIL compound annotation")
  suppressMessages(suppressWarnings(require(grid)))
  if (is.null(isotopes)){data("isotopes")}
  
  
  formula<-as.character(formula)
  
  # define the matching score calculation function

  
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

atomes.to.formula<-function(unique_formula_list){
  lapply(names(unique_formula_list),function(x){
    unique_formula_list[[x]]->temp
    temp[temp==0]<-NULL
    temp_for<-paste0(names(temp),unlist(temp),collapse = "")
    return(temp_for)
  })->unique_formula_list_for
  return(unlist(unique_formula_list_for))
}
