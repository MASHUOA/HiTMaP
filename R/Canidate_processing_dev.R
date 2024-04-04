
Protein_feature_list_table_import<-function(workdir=getwd(),
                                            ptable,
                                            ptable_pipeline=c("fragpipe","UserTable","maxQuant","ProteomeDiscover","ProteinPilot"),
                                            target_GroupID=NULL,
                                            database=NULL,
                                            missedCleavages=0:1,
                                            adducts=c("M+H","M+NH4","M+Na"),
                                            BPPARAM=bpparam(),
                                            Decoy_search=T,
                                            Decoy_mode=c("isotope","adducts","elements","sequence"),
                                            Decoy_adducts=c("M+He","M+Ne","M+Ar","M+Kr","M+Xe","M+Rn"),
                                            Substitute_AA=list(AA=c(NULL),AA_new_formula=c(NULL),Formula_with_water=c(NULL)),
                                            mzrange=c(500,4000),
                                            output_candidatelist=T,
                                            mod_include_conversion=c('Methyl'="Methylation",
                                                                     'Thiazolidine'="formaldehyde adduct",
                                                                     'hydroxymethyl'="hydroxymethyl/formaldehyde induced modifications",
                                                                     'Deamidation'="Deamidation",
                                                                     'Formyl'="Formylation",
                                                                     'Hydroxylation'="Oxidation or Hydroxylation",
                                                                     'Pyro-glu'="Pyro-glu from Q/Loss of ammonia",
                                                                     'Acetyl'="N-term(42.0106)"
                                            ),
                                            use_previous_candidates=F,
                                            Protein_desc_of_exclusion=NULL,
                                            Database_stats=F
){
  
  library(curl)
  library(reshape2)
  library(dplyr)
  library(xml2)
  library(jsonlite)
  library(IRanges)
  library(readr)
  library(Biostrings)
  library(readxl)
  library(stringr)
  library(magick)
  library(scales)
  if (!dir.exists(workdir)) dir.create(workdir,recursive = T)
  setwd(workdir)
  final_res<-NULL
  ptable->ptable_path
  if ('|'(missing(target_GroupID),is.null(target_GroupID))) message("No input target found. Will use all protein groups in the file.")
  if ('|'(missing(ptable),is.null(ptable))) stop("No identification data found...")
  
  if (is.character(target_GroupID)){
    if( '&'(length(target_GroupID)==1,sum(file.exists(target_GroupID))==1) ){
      read_lines(target_GroupID)->target_GroupID
    }else{
      target_GroupID=target_GroupID
    }
  } 
  
  required_col<-c("Protein","Peptide","start","end",
                 "Modification","pro_end","pepmz",
                 "formula","mz","adduct","isdecoy","charge","desc")
  
  
  if (is.character(ptable)){
    if( '&'(length(ptable)==1,sum(file.exists(ptable))==1) ){
      if (ptable_pipeline[1]=="ProteinPilot") {
        ptable=read_excel(ptable, sheet = "Processed Peptide Summary")
        required_col<-c(Pro_col<-"Accessions",
                        Pep_col<-"Sequence" ,
                        Pep_mod_col<-"Modifications",
                        RT_col<-"Acq.Time",
                        Cleavage_col<-"Cleavages",
                        Chanrge_col<-"Theor.z",
                        Score_col<-"Conf",
                        Intensities_col<-"Intensity.Peptides")
      }
      if (ptable_pipeline[1]=="fragpipe") {
        ptable=read_delim(ptable,delim = "\t", escape_double = FALSE, trim_ws = TRUE)
        required_col<-c(Pro_col<-"Protein",
                        Pep_col<-"Peptide" ,
                        #RT_col<-"RT",
                        pro_start<-"Protein Start",
                        pro_end<-"Protein End",
                        Pep_mod_col<-"Assigned Modifications",
                        Pep_mod_col_obs<-"Observed Modifications",
                        Cleavage_col<-"Number of Missed Cleavages",
                        Pro_desc="Protein Description")
        if (sum(!(required_col %in% colnames(ptable)))>=1){
          stop("Required column(s) is missing in identification data:", paste0(required_col[!(required_col %in% colnames(ptable))],collapse = ","))
        }
        
        ptable=ptable[,required_col]
        ptable=unique(ptable)
        ptable[is.na(unlist(ptable[,Pep_mod_col])),Pep_mod_col]<-""
        ptable[,Pep_mod_col]->mod_excl_id
        lapply(unlist(mod_excl_id),str_split_1,pattern=', ')->mod_split
        lapply(mod_split,function(x){
          x<-unlist(x)
          names(mod_include_conversion[mod_include_conversion %in% x])->x
          paste0(x,collapse = ";")
          })->mod_split_final
        ptable[,Pep_mod_col]<-unlist(mod_split_final)
        
        unname(unlist(ptable[,Pep_mod_col_obs]))->mod_obs_id
        lapply((mod_obs_id),str_split,pattern='Mod.: ')->mod_split_obs
        lapply(mod_split_obs,function(x){
          unlist(x)->x
          if (!is.na(x[1])) {
            str_remove_all(x," \\(PeakApex.{1,1000}\\)|\\(Theoretical.{1,1000}\\)|\\;|,")->x
            str_remove_all(x," $")->x
            x<-x[!str_detect(x,"First isotopic peak")]
            x<-x[!str_detect(x,"Unannotated mass-shift")]
            (x[x!=""])->x
            names(mod_include_conversion[mod_include_conversion %in% x])->x
            paste0(x,collapse = ";")
          }else{""}
        })->mod_split_obs
        (unlist(str_split(mod_split_obs,";")))->all_ob_mod
        table(all_ob_mod)->all_ob_mod_tb
        png("mod.stats.png",width = 8000,height = 5000, res=300,units = "px")
        barplot(sort(all_ob_mod_tb,decreasing = T))
        dev.off()
        ptable[,Pep_mod_col]<-paste0(unlist(ptable[,Pep_mod_col]),";",unlist(mod_split_obs))
        ptable[,Pep_mod_col]<-str_remove(unlist(ptable[,Pep_mod_col]),";$|^;")
        ptable[,!(colnames(ptable) %in% Pep_mod_col_obs)]->ptable
        ptable=unique(ptable)
        colnames(ptable)<-c("Protein","Peptide","start","end","Modification","Cleavages","desc")
        ptable->ptable_candidate
        
      }
      if (ptable_pipeline[1]=="maxQuant") {
        ptable=read_delim(ptable,delim = "\t", escape_double = FALSE, trim_ws = TRUE)
        required_col<-c(Pro_col<-"Proteins",
                        Pep_col<-"Sequence" ,
                        #RT_col<-"RT",
                        pro_start<-"Protein Start",
                        pro_end<-"Protein End",
                        Pep_mod_col<-"Modifications",
                        Cleavage_col<-"Missed cleavages",
                        Pro_desc="Protein Names")
        
        required_col_org<-required_col
        misscol<-NULL
        if (sum(!(required_col %in% colnames(ptable)))>=1){
          message("Required column(s) is missing in identification data:", paste0(required_col[!(required_col %in% colnames(ptable))],collapse = ","),". Will replac these columns by default")
          which(!(required_col %in% colnames(ptable))) -> misscol    
          required_col<-required_col[-misscol]
          misscol<-required_col_org[misscol]
        }
        

        ptable=ptable[,required_col]
        ptable=unique(ptable)
        ptable[is.na(unlist(ptable[,Pep_mod_col])),Pep_mod_col]<-""
        ptable[,Pep_mod_col]->mod_excl_id
        lapply(unlist(mod_excl_id),str_split_1,pattern=', ')->mod_split
        lapply(mod_split,function(x){
          x<-unlist(x)
          names(mod_include_conversion[mod_include_conversion %in% x])->x
          paste0(x,collapse = ";")
        })->mod_split_final
        ptable[,Pep_mod_col]<-unlist(mod_split_final)
        
        ptable[,Pep_mod_col]<-str_remove(unlist(ptable[,Pep_mod_col]),";$|^;")
        ptable=unique(ptable)
        
        if (!is.null(misscol)){
          for (eachcol in misscol){
            ptable[,eachcol]<-1
          }
        }
        ptable<-ptable[,required_col_org]
        
        colnames(ptable)<-c("Protein","Peptide","start","end","Modification","Cleavages","desc")
        ptable->ptable_candidate
        
      }
      if (ptable_pipeline[1] %in% c("UserTable")) {
        required_col<-c(Pro_col<-"Protein",
                        Pep_col<-"Sequence" ,
                        Pep_mod_col<-"Modification",
                        RT_col<-"RT",
                        Cleavage_col<-"Cleavages",
                        Chanrge_col<-"Charge",
                        Score_col<-"Score",
                        Intensities_col<-"Area")
        ptable=read.csv(ptable)
      }
      
    } else {
      stop("Identification data is invalid.")
    }
    
  } else if (is.data.frame(ptable)){
    target_GroupID=target_GroupID
  } else {
    stop("Identification data is invalid.")
  }
  
  
  
  
  if (mzrange[1] == "auto-detect") {
    mzrange_proteomics_candidates_prefiltering<-c(500,4000)
    mzrange=mzrange_proteomics_candidates_prefiltering
  }
  
  suppressMessages(suppressWarnings(require(Biostrings)))
  suppressMessages(suppressWarnings(require(cleaver)))
  suppressMessages(suppressWarnings(require(protViz)))
  suppressMessages(suppressWarnings(require(rcdk)))
  suppressMessages(suppressWarnings(require(BiocParallel)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(rJava)))
  suppressMessages(suppressWarnings(require(rcdklibs)))
  suppressMessages(suppressWarnings(require(grid)))
  suppressMessages(suppressWarnings(require(stringr)))
  setwd(workdir)
  suppressMessages(suppressWarnings(require(enviPat)))
  suppressMessages(suppressWarnings(require(rcdklibs)))
  data("isotopes")
  # apply missedCleavages filter
  missedCleavages<<-missedCleavages
  ptable_candidate[ptable_candidate$Cleavages %in% missedCleavages,]->ptable_candidate
  
  
  Decoy_adducts=Decoy_adducts[!(Decoy_adducts %in% adducts)]
  Decoy_adducts=Decoy_adducts[1:length(adducts)]
  
  
  
  if (!is.null(database)){
  
    if (file.exists(database)){
    Index_of_protein_sequence <<- fasta.index(database,
                                            nrec=-1L, 
                                            skip=0L)  
  
  
  list_of_protein_sequence <- readAAStringSet(database,
                                              format="fasta",
                                              nrec=-1L, 
                                              skip=0L, 
                                              seek.first.rec=FALSE
  )}else{
    message("Protein DB is missing, Bypass Protein end modification.")
  }
  }
  
  if(use_previous_candidates){
    #if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
    if (sum(c("candidatelist.csv","protein_index.csv") %in% dir(paste(workdir,"/Summary folder",sep="")))==2){
      Protein_Summary<-read.csv(paste(workdir,"/Summary folder/candidatelist.csv",sep=""),)
      Protein_Summary$Modification[is.na(Protein_Summary$Modification)]<-""
      Index_of_protein_sequence<<-read.csv(paste(workdir,"/Summary folder/protein_index.csv",sep=""))
      message("Candidate list has been loaded.") 
    }else{
      message("Can not find the previously established candidate list.")
      use_previous_candidates=F
    }
    
  }
  
  if(use_previous_candidates!=T){
    
    
    AA<-c(71.037114, 0.000000, 103.009185, 115.026943, 129.042593, 147.068414, 
          57.021464, 137.058912, 113.084064, 0.000000, 128.094963, 113.084064, 
          131.040485, 114.042927, 0.000000, 97.052764, 128.058578, 156.101111, 
          87.032028, 101.047679, 0.000000, 99.068414, 186.079313, 0.000000, 
          163.063329, 100.994269)
    
    names(AA)<-LETTERS
    message(paste("Generated",length(unique(ptable_candidate$Protein)),"Proteins in total. Computing exact masses..."))

    ptable_candidate->tempdf
    if (length(grep("X",tempdf$Peptide))!=0 && !("X" %in% Substitute_AA$AA))  tempdf<-tempdf[-grep("X",tempdf$Peptide),]
    if (length(grep("U",tempdf$Peptide))!=0 && !("U" %in% Substitute_AA$AA))  tempdf<-tempdf[-grep("U",tempdf$Peptide),]
    
    mod.df_fix<-Peptide_modification(retrive_ID = names(mod_include_conversion))
    mod.df<-rbind(mod.df_fix)
    mod.df_unique<-unique(mod.df[,c("record_id","code_name","full_name","ex_code_name","composition", "mono_mass" )])
    
    
    
    if (is.null(mod.df)){
      min_mod_massdiff<-100
      max_mod_massdiff<-500
    }else{
      min_mod_massdiff<-min(as.numeric(mod.df$mono_mass))
      max_mod_massdiff<-max(as.numeric(mod.df$mono_mass))
      
      if(min_mod_massdiff>0){min_mod_massdiff=0}
      if(max_mod_massdiff<0){max_mod_massdiff=0}
    }
    
    #parentIonMass gives the M+H
    tempdf$pepmz <- as.numeric(parentIonMass(tempdf$Peptide,fixmod=AA)- 1.007276 )
    tempdf<-tempdf['&'(tempdf$pepmz>=mzrange[1]-max_mod_massdiff,tempdf$pepmz<=mzrange[2]+min_mod_massdiff),]
    Protein_Summary<-NULL
    adductslist<-Build_adduct_list()
    
    message(paste("Generating peptide formula..."))
    uniquepep=(tempdf$Peptide)
    uniqueAA=unique(strsplit(paste0(uniquepep,collapse = ""), "")[[1]])
    
    
    
    Element_tbl<-BuildElement(Substitute_AA=Substitute_AA,uniqueAA=uniqueAA)
    peptide_symbol=lapply(uniquepep,ConvertPeptide,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl)
    names(peptide_symbol)<-paste0(uniquepep,"_",tempdf$Modification)
    
    peptide_symbol_mod=convert_peptide_modsymbol(mod.df_unique,peptide_symbol,peptide_info=tempdf,BPPARAM = BPPARAM)
    tempdf$pepmz_mod<-round(unlist(lapply(peptide_symbol_mod,MonoisotopicMass)),4)
    tempdf$diff<-round(as.numeric(tempdf$pepmz_mod)-as.numeric(tempdf$pepmz),3)
    check_mz_df<-unique(tempdf[,c("Modification","diff")])
    message(paste("Generating peptide formula with adducts:",paste(adducts,collapse = " ")))
    peptides_symbol_mod_adducts=bplapply(adducts,convert_peptide_adduct_list,peptide_symbol_mod,BPPARAM = BPPARAM,adductslist=adductslist)
    
    for (i in 1:length(adducts)){
      adductmass <- as.numeric(as.character(adductslist[adductslist$Name == adducts[i], "Mass"]))
      charge=as.numeric(as.character(adductslist$Charge[adductslist$Name == adducts[i]]))
      tempdf$formula<-peptides_symbol_mod_adducts[[i]]
      message(paste("Calculating peptide mz with adducts:",adducts[i]))
      if (charge==0){actingcharge=1} else {actingcharge=abs(charge)}
      tempdf$mz<-(tempdf$pepmz+adductmass)/actingcharge
      tempdf$adduct<-adducts[i]
      tempdf$isdecoy<-rep(0,nrow(tempdf))
      tempdf$charge<-charge
      Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
    }
    
    if (length(Decoy_adducts)>0 && Decoy_search && ("adducts" %in% Decoy_mode)){
      message(paste("Generating peptide formula with Decoy adducts:",paste(Decoy_adducts,collapse = " ")))
      peptides_symbol_adducts=bplapply(Decoy_adducts,convert_peptide_adduct_list,peptide_symbol,BPPARAM = BPPARAM,adductslist=adductslist)
      for (i in 1:length(Decoy_adducts)){
        adductmass <- as.numeric(as.character(adductslist[adductslist$Name == Decoy_adducts[i], "Mass"]))
        charge=as.numeric(as.character(adductslist$Charge[adductslist$Name==Decoy_adducts[i]]))
        tempdf$formula<-peptides_symbol_adducts[[i]]
        message(paste("Calculating peptide mz with Decoy_adducts:",Decoy_adducts[i]))
        if (charge==0){actingcharge=1} else {actingcharge=abs(charge)}
        tempdf$mz<-(tempdf$pepmz+adductmass)/actingcharge
        tempdf$adduct<-Decoy_adducts[i]
        tempdf$isdecoy<-rep(1,nrow(tempdf))
        tempdf$charge<-charge
        Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
      } 
      
    }
    
    
    Protein_Summary$Modification[is.na(Protein_Summary$Modification)]<-""
    Protein_Summary$mz<-round(Protein_Summary$mz,digits = 4)
    Protein_Summary<-Protein_Summary[`&`(Protein_Summary$mz>=mzrange[1],Protein_Summary$mz<=mzrange[2]),]
    if (!exists(Index_of_protein_sequence)){
    Index_of_protein_sequence<-data.frame(unique(Protein_Summary[,c("Protein","desc")]))
    Index_of_protein_sequence$recno<-Index_of_protein_sequence$Protein
    Index_of_protein_sequence$fileno<-1
    Index_of_protein_sequence$offset<-0
    Index_of_protein_sequence$seqlength<-0
    Index_of_protein_sequence$filepath<-ptable_path
    Protein_Summary$pro_end<-0
    }
    
    Protein_Summary$desc<-NULL
    
    if (output_candidatelist){
      if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
      write.csv(Protein_Summary,paste(workdir,"/Summary folder/candidatelist.csv",sep=""),row.names = F)
      write.csv(Index_of_protein_sequence,paste(workdir,"/Summary folder/protein_index.csv",sep=""),row.names = F)
      message("Candidate list has been exported.")
    }
  }
  
  if (!is.null(Protein_desc_of_exclusion)){
    Protein_Summary<-merge(Protein_Summary,Index_of_protein_sequence[,c("recno","desc")],by.x="Protein",by.y="recno",all.x=T)
    Protein_Summary_exclusion<-NULL
    num_of_interest<-numeric(0)
    for (interest_desc in Protein_desc_of_exclusion){
      num_before=nrow(Protein_Summary_exclusion)
      if(is.null(num_before)) num_before=0
      Protein_Summary_exclusion<-rbind(Protein_Summary_exclusion,Protein_Summary[grepl(paste0(" ",interest_desc),Protein_Summary$desc,ignore.case = T),])
      Protein_Summary_exclusion<-rbind(Protein_Summary_exclusion,Protein_Summary[grepl(paste0("-",interest_desc),Protein_Summary$desc,ignore.case = T),])
      num_after=nrow(Protein_Summary_exclusion)
      if(is.null(num_after)) num_after=0
      num_of_interest[interest_desc]<-nrow(unique(Protein_Summary[grepl(paste0(" ",interest_desc),Protein_Summary$desc,ignore.case = T),]))+nrow(unique(Protein_Summary[grepl(paste0("-",interest_desc),Protein_Summary$desc,ignore.case = T),]))
    }
    #Protein_feature_list_crystallin$Protein=as.character(Protein_feature_list_crystallin$desc)
    Protein_Summary=Protein_Summary[!(Protein_Summary$Protein %in% unique(Protein_Summary_exclusion$Protein)),]
    Protein_Summary$desc<-NULL
    message(paste(num_of_interest,"Protein(s) found with annotations of exclusion:",Protein_desc_of_exclusion,collapse = "\n"))     
  }
  
  if(Database_stats){
    
    suppressMessages(suppressWarnings(library(dplyr)))
    suppressMessages(suppressWarnings(library(egg)))
    suppressMessages(suppressWarnings(library(RColorBrewer)))
    mz_vs_formula<-Protein_Summary %>% group_by(mz) %>% summarize(Intnsity=length(unique(formula)))
    
    bin_mz_ppm<-function(mz_vs_peptide,ppm_test_list=c(1,2,5)){
      mz_vs_peptide_filtered<-list() 
      ppm_test_list<-sort(ppm_test_list,decreasing = T)
      for (ppm_test in ppm_test_list){
        
        mz_vs_peptide_filtered[[as.character(paste(ppm_test,"ppm"))]]<-isopattern_ppm_filter_peaklist_par(mz_vs_peptide,ppm_test,threshold=0.00)
        
      }
      
      mz_vs_peptide_filtered<-lapply(mz_vs_peptide_filtered,`colnames<-`,c("mz","unique_formula"))
      
      mycol=RColorBrewer::brewer.pal(name = "Set1",n=length(names(mz_vs_peptide_filtered)))
      names(mycol)=names(mz_vs_peptide_filtered)
      mycol_name<-mycol
      names(mycol_name)<-mycol
      sp<-ggplot2::ggplot() 
      sp <-sp + geom_segment(data=mz_vs_peptide_filtered[[1]],mapping = aes(x=mz, y=unique_formula,xend=mz,yend=rep(0,length(unique_formula)),col=names(mycol)[1]),size=0.1)
      sp <-sp + geom_segment(data=mz_vs_peptide_filtered[[2]],mapping = aes(x=mz, y=unique_formula,xend=mz,yend=rep(0,length(unique_formula)),col=names(mycol)[2]),size=0.1)
      sp <-sp + geom_segment(data=mz_vs_peptide_filtered[[3]],mapping = aes(x=mz, y=unique_formula,xend=mz,yend=rep(0,length(unique_formula)),col=names(mycol)[3]),size=0.1)

      sp <- sp + theme_article() + 
        
        theme(legend.position = "top",axis.text=element_text(size=12),
              axis.title=element_text(size=14,face="bold"),
              legend.text =element_text(size=14,face="bold"),
              legend.title =element_text(size=14,face="bold"),
              legend.key = element_rect(fill = "white"))+
        guides(color=guide_legend("m/z bin size"),override.aes = list(linetype = 0, size = 5)) +
        labs(title = "",y = "Unique formula",x = "m/z")
      
      return(sp)
      
    }
    
    bin_mz_ppm_pic<-bin_mz_ppm(mz_vs_formula)
    
    png(paste(workdir,"/Summary folder/DB_stats_bin_mz_ppm.png",sep=""),width = 1200,height = 800,res = 150)
    
    print(bin_mz_ppm_pic)
    
    dev.off()
    
    mz_unqiue<-unique(Protein_Summary$mz)
    
    mz_unqiue_formula<-Protein_Summary %>% group_by(mz) %>% summarize(formula=paste(unique(formula),sep="; "))
    
    mz_unqiue<-sort(mz_unqiue)
    
    mz_unqiue_diff<-diff(mz_unqiue)
    
    mz_unqiue_center<-zoo::rollapply(mz_unqiue, 2, mean, by = 1, align = "left", partial = FALSE)
    
    #dif_kmeans=kmeans(1/(mz_vs_resolution$mz_unqiue_diff),centers = 200,iter.max = 500)
    
    min_formula_L<-mz_unqiue_formula[sort(which(mz_unqiue_diff==min(mz_unqiue_diff))),]
    min_formula_R<-mz_unqiue_formula[sort(which(mz_unqiue_diff==min(mz_unqiue_diff))+1),]
    
    png(paste(workdir,"/Summary folder/DB_stats_mz_diff_resolution.png",sep=""),width = 1200,height = 800,res = 150)
    
    mz_vs_resolution<-data.frame(mz=mz_unqiue_center,Resolution=mz_unqiue_center/mz_unqiue_diff,mz_unqiue_diff=mz_unqiue_diff)
    
    sp<-ggplot2::ggplot() 
    sp <-sp + geom_point(data=mz_vs_resolution,mapping = aes(x=mz, y=Resolution),size=0.3,alpha=0.4)+
      theme_article() + 
      
      theme(legend.position = "top",axis.text=element_text(size=12),
            axis.title=element_text(size=14,face="bold"),
            legend.text =element_text(size=14,face="bold"),
            legend.title =element_text(size=14,face="bold"),
            legend.key = element_rect(fill = "white"))+
      guides(color=guide_legend("m/z bin size"),override.aes = list(linetype = 0, size = 5)) +
      labs(title = "",y = "Required resolution",x = "m/z")
    print(sp)
    
    dev.off()
    
  }
  

  
  Protein_feature_list<<-Protein_Summary 
  
  return(Protein_Summary)
}

convert_peptide_modsymbol<-function(mod.df,peptide_symbol,peptide_info,BPPARAM=BPPARAM){
  suppressMessages(suppressWarnings(require(rcdk)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(stringr)))
  suppressMessages(suppressWarnings(require(Biostrings)))
  suppressMessages(suppressWarnings(require(stats)))
  pep_sequence=peptide_info$Peptide
  peptide_seqinfo=peptide_info[,c("start","end")]
  
  merge_atoms<-function(atoms,addelements,check_merge=T,mode=c("add","ded"),multiplier=c(1,1)){
    
    atomsorg=atoms
    
    if (missing(mode)) mode="add"
    
    if (mode=="ded"){
      for (x in names(addelements)){
        addelements[[x]]=-addelements[[x]]
      }
    }
    
    for (x in names(addelements)){
      if (is.null(atoms[[x]])){
        atoms[[x]]=addelements[[x]]
      }else {
        atoms[[x]]=(atoms[[x]]*multiplier[1]) + (addelements[[x]]*multiplier[2])
      } 
    }
    
    if (check_merge==F){
      for (x in names(atoms)){
        if (atoms[[x]]<0){
          stop("add atoms failed due to incorrect number of",x)
        } 
      }
      
    }
    return(as.list(atoms))
    
  }
  
  get_atoms_mod<-function(Symbol){
    #form = "C5H11BrO" 
    if (Symbol!=""){
      ups = c(gregexpr("[[:upper:]]", Symbol)[[1]], nchar(Symbol) + 1) 
      seperated = sapply(1:(length(ups)-1), function(x) substr(Symbol, ups[x], ups[x+1] - 1)) 
      elements =  gsub("[[:digit:]]", "", seperated) 
      elements =  gsub(" ", "", elements)
      elements =  gsub("\\)", "", elements)
      elements =  gsub("\\(", "", elements)
      elements =  gsub("-", "", elements)
      nums = gsub("[[:alpha:]]", "", seperated) 
      nums = gsub(" ", "", nums) 
      nums = gsub("\\(", "", nums) 
      nums = gsub("\\)", "", nums) 
      ans = data.frame(elements = as.character(elements), num = as.numeric(ifelse(nums == "", 1, nums)), stringsAsFactors = FALSE)
      list<-as.list(ans$num)
      names(list)=ans$elements
      return(list)
    }else{return(NULL)}
  }
  
  mod.df.list <-unique(mod.df[,c("record_id","code_name","full_name","ex_code_name","composition", "mono_mass")])
  formula_mod<-lapply(mod.df.list$composition,get_atoms_mod)
  names(formula_mod)<-mod.df.list$code_name
  pep_mod_list<-str_split(peptide_info$Modification,";")
  names(pep_mod_list)<-names(peptide_symbol)
  
  lapply(names(pep_mod_list),function(x){
    new_symbol<-peptide_symbol[[x]]
    for (i in pep_mod_list[[x]]){
      new_symbol<-merge_atoms(new_symbol,formula_mod[[i]],check_merge=T,mode="add",multiplier=c(1,1))
    }
    new_symbol_v<-unlist(new_symbol)
    names(new_symbol_v)<-names(new_symbol)
    return(new_symbol)
  })->peptide_symbol_new
  names(peptide_symbol_new)<-names(peptide_symbol)
  
  
  
  return(peptide_symbol_new)
  
}
