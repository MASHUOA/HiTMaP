
#' Meta_feature_list_fun
#'
#' This is a function that prepare the candiate list for maldi imaging data qualitative or quantitative analysis.
#' this function will read the candidate list file and generate mz for the adducts list defined in \code{"adducts"}. 
#' @param database the file name of candidate list
#' @param workdir the folder that contains candidate list
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param cal.mz If set with \code{"TRUE"}, the function will recalculate the mz value according to the column named "Formula" in the \code{database} and the specified adducts.
#' @param mzlist_bypass  Set \code{"TRUE"} if you want to bypass the mzlist generating process, the function will keep the mz and adduct as it is for the furture analysis. Be sure that the table contains "mz", "adduct" and "moleculeNames" as they are essential for later steps.
#' @param BPPARAM parallel processing parameter for BiocParallel
#' @return a table of candiate list
#'
#' @examples
#' Meta_feature_list_fun(database="lipid candidates.csv",adducts=c("M-H","M+Cl"))
#'
#' @export


Meta_feature_list_fun<-function(database,
                                workdir=getwd(),
                                adducts=adducts,
                                cal.mz=TRUE,
                                bypass=FALSE,
                                BPPARAM=bpparam()){
  library("Rcpp")
  library(dplyr)
  library(Rdisop)
  library(Biostrings)
  
  adductslist<-Build_adduct_list()
  candidates<-read.csv(paste0(workdir,"/",database),as.is = TRUE)
  
  
  if (bypass){
    required_col=c("mz","adduct","moleculeNames")
    candidates=data_test_rename(required_col,candidates)
    Meta_Summary=candidates[is.na(candidates$mz)!=T,]
    
  }else{
    required_col=c("moleculeNames")
    candidates=data_test_rename(required_col,candidates)
    
       candidates$mass=0
    if (cal.mz==F){
    required_col=c("mz")
    candidates=data_test_rename(required_col,candidates)
      if (sum(candidates$adducts!=0)>0){
        candidates$adducts_mass=0
        candidates$adducts_mass=unlist(bplapply(1:nrow(candidates),function(i,adducts,adductslist,adducts_mass){
          adducts_mass[i]<-adductslist[adductslist$Name==adducts[i],"Mass"]*abs(as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Charge"])))
        },candidates$adducts,adductslist,candidates$adducts_mass,BPPARAM=BPPARAM))
        candidates$mass<-candidates$mz-candidates$adducts_mass
      }else{
        candidates$mass<-candidates$mz
        }
      }else{
        
        required_col=c("Formula")
        candidates=data_test_rename(required_col,candidates)
        
      candidates_mass<-candidates$Formula %>% lapply(getMonomass) 
      candidates$mass<-as.numeric(unlist(candidates_mass))
    }
    candidates<-candidates[duplicated(names(candidates))==FALSE]
    
    Meta_Summary<-NULL
    
    for (i in 1:length(adducts)){
      adductmass<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Mass"]))
      candidates$mz<-as.numeric(candidates$mass+adductmass)
      candidates$adduct<-adducts[i]
      Meta_Summary<-rbind.data.frame(Meta_Summary,candidates)
    }
    } 
    return(Meta_Summary)                
  
} 


#' Protein_feature_list_fun
#'
#' This is a function that prepare the candiate list for maldi imaging data qualitative or quantitative analysis.
#' this function will read the fasta file and generate mz for the adducts list defined in \code{"adducts"}. 
#' @param workdir the folder that contains fasta file
#' @param database the file name of fasta file
#' @param Digestion_site Digestion pattern, this is the 
#' @param missedCleavages Define a number range for misscleaved peptides to be considered
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param BPPARAM parallel processing parameter for BiocParallel
#' @return a table of peptide candiate list
#'
#' @examples
#' Protein_feature_list_fun(database="lipid candidates.csv",adducts=c("M+H","M+NH4","M+Na"))
#'
#' @export

Protein_feature_list_fun<-function(workdir=getwd(),
                                   database,
                                   Digestion_site="[G]",
                                   missedCleavages=0:1,
                                   adducts=c("M+H","M+NH4","M+Na"),
                                   BPPARAM=bpparam(),
                                   Decoy_adducts=c("M+He","M+Ne","M+Ar","M+Kr","M+Xe","M+Rn"),
                                   Decoy_mode=c("adducts","elements","isotope"),
                                   Decoy_search=T,
                                   mzrange=c(500,4000),
                                   output_candidatelist=T,
                                   use_previous_candidates=F){
  library(Biostrings)
  library(cleaver)
  library(protViz)
  library(rcdk)
  library(BiocParallel)
  library(OrgMassSpecR)
  library(rJava)
  library(rcdklibs)
  library(grid)
  library(stringr)
  setwd(workdir)
  
  if(use_previous_candidates!=T){
  Decoy_adducts=Decoy_adducts[!(Decoy_adducts %in% adducts)]
  Decoy_adducts=Decoy_adducts[1:length(adducts)]
  list_of_protein_sequence<<-readAAStringSet(database,
                                            format="fasta",
                                            nrec=-1L, 
                                            skip=0L, 
                                            seek.first.rec=FALSE,
                                            use.names=TRUE, 
                                            with.qualities=FALSE)    
 
  #if (length(list_of_protein_sequence)<2000){bpworkers(BPPARAM)=3}
  
  #if (length(list_of_protein_sequence)<500){bpworkers(BPPARAM)=1}
  
  
  
  
  Index_of_protein_sequence<<-fasta.index(database,
                                         nrec=-1L, 
                                         skip=0L)   
  names_pro<-data.frame(desc=names(list_of_protein_sequence),stringsAsFactors = F)
  
  names_pro<-merge(names_pro,Index_of_protein_sequence,by="desc")
  
  names(list_of_protein_sequence)<-names_pro$recno
  
  Index_of_protein_sequence$Degestion=""
  
  peplist<-list()
  
  Index_of_protein_sequence_list<-data.frame() 
  
  parentIonMasslist<-function(peplist,Index_of_protein_sequence){
    AA<-rep(0,26)
    PIM<-NULL
    for (i in 1:length(peplist)){
      PIM[[Index_of_protein_sequence$recno[i]]] <- parentIonMass(peplist[[Index_of_protein_sequence$recno[i]]],fixmod=AA)}
    return(PIM)
  }
  
  for (i in 1:length(Digestion_site)){
    
    if (Digestion_site[i]==""){Digestion_site[i]="J"}
    
    peplist_option<-cleave(as.character(list_of_protein_sequence),custom=Digestion_site[i], missedCleavages=missedCleavages)
    
    Index_of_protein_sequence_option<-Index_of_protein_sequence
    
    #Index_of_protein_sequence_option$desc<-paste(names(peplist_option),Digestion_site[i])
    
    #names(peplist_option)=paste(names(peplist_option),Digestion_site[i])
    
    
    
    pimlist<-bplapply(peplist_option,function(x){
    AA<-rep(0,26)
    protViz::parentIonMass(x,fixmod=AA)},BPPARAM = BPPARAM)
    
    
    
    peplist_option<-peplist_option[duplicated(names(peplist_option))==FALSE]
    
    message(paste("Testing fasta sequances..."))
    
    peplist_option<-peplist_option[names(peplist_option) %in% names(pimlist) ==TRUE]
    
    Index_of_protein_sequence_option<-Index_of_protein_sequence_option[Index_of_protein_sequence_option$recno %in% names(pimlist),]
    
    Index_of_protein_sequence_option$Degestion=Digestion_site[i]
    
    Index_of_protein_sequence_option
    
    peplist<-c(peplist,peplist_option)
    
    Index_of_protein_sequence_list<-rbind(Index_of_protein_sequence_list,Index_of_protein_sequence_option)
    
  }
  
  #message(paste("Peptide list generated",length(peplist),"entries in total."))
  
  
  
  
  
  #pimlist<-parentIonMasslist(peplist,Index_of_protein_sequence_list)
  
  
  AA<-c(71.037114, 114.534940, 103.009185, 115.026943, 129.042593, 147.068414, 
        57.021464, 137.058912, 113.084064, 0.000000, 128.094963, 113.084064, 
        131.040485, 114.042927, 0.000000, 97.052764, 128.058578, 156.101111, 
        87.032028, 101.047679, 150.953630, 99.068414, 186.079313, 111.000000, 
        163.063329, 100.994269)
  #tempdf<-parLapply(cl=cl,  1: length(names(peplist)), Peptide_Summary_para,peplist)
  #bplapply()
  message(paste("Generated",length(peplist),"Proteins in total. Computing exact masses..."))
  tempdf<-bplapply( 1: length(names(peplist)), Peptide_Summary_para,peplist,BPPARAM = BPPARAM)
  tempdf <- do.call("rbind", tempdf)
  colnames(tempdf)<-c("Protein","Peptide")
  
  tempdf<-as.data.frame(tempdf)
  
  tempdf$Protein<-as.character(tempdf$Protein)
  tempdf$Peptide<-as.character(tempdf$Peptide)

  tempdf<-tempdf[-grep("X",tempdf$Peptide),]
  tempdf<-tempdf[-grep("U",tempdf$Peptide),]
  
  tempdf$pepmz <- as.numeric(parentIonMass(tempdf$Peptide,fixmod=AA)- 1.007276 )
  tempdf<-tempdf[`&`(tempdf$pepmz>=mzrange[1],tempdf$pepmz<=mzrange[2]),]
  #tempdf1<-tempdf
  Protein_Summary<-NULL
  adductslist<-Build_adduct_list()
  
  message(paste("Generating peptide formula..."))
  peptide_symbol=bplapply(tempdf$Peptide,ConvertPeptide,BPPARAM = BPPARAM)
  message(paste("Generating peptide formula with adducts:",paste(adducts,collapse = " ")))
  peptides_symbol_adducts=bplapply(adducts,convert_peptide_adduct_list,peptide_symbol,BPPARAM = BPPARAM,adductslist=adductslist)
  
  for (i in 1:length(adducts)){
    adductmass <- as.numeric(as.character(adductslist[adductslist$Name == adducts[i], "Mass"]))
    charge=as.numeric(as.character(adductslist$Charge[adductslist$Name==adducts[i]]))
    #tempdf$formula[1:1000]<-as.character(peptides_symbol_adducts[[i]])
    tempdf$formula<-peptides_symbol_adducts[[i]]
    message(paste("Calculating peptide mz with adducts:",adducts[i]))
    #3templist=bplapply(tempdf$formula,function(x,charge){
    #  rcdk::get.formula(x,charge = charge)@mass
    #  },charge,BPPARAM = BPPARAM)
    if (charge==0){actingcharge=1} else {actingcharge=abs(charge)}
    #tempdf$mz<-as.numeric(unlist(templist))
    tempdf$mz<-(tempdf$pepmz+adductmass)/actingcharge
    tempdf$adduct<-adducts[i]
    tempdf$isdecoy<-rep(0,nrow(tempdf))
    tempdf$charge<-charge
    #tempdf$isdecoy<-rep(0,nrow(tempdf))
    #convert_peptide_adduct(tempdf$Peptide[2],adductsname = adducts[i],multiplier = c(multiplier,1),adductslist = adductslist)
    #tempdf$formula<-bplapply(tempdf$Peptide,convert_peptide_adduct,adductsformula = adducts[i],multiplier = c(multiplier,1),BPPARAM = BPPARAM)
    Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
  }
  
  if (length(Decoy_adducts)>0 && Decoy_search && ("adducts" %in% Decoy_mode)){
    message(paste("Generating peptide formula with Decoy adducts:",paste(Decoy_adducts,collapse = " ")))
    peptides_symbol_adducts=bplapply(Decoy_adducts,convert_peptide_adduct_list,peptide_symbol,BPPARAM = BPPARAM,adductslist=adductslist)
  for (i in 1:length(Decoy_adducts)){
    adductmass <- as.numeric(as.character(adductslist[adductslist$Name == Decoy_adducts[i], "Mass"]))
    charge=as.numeric(as.character(adductslist$Charge[adductslist$Name==Decoy_adducts[i]]))
    #tempdf$formula[1:1000]<-as.character(peptides_symbol_adducts[[i]])
    tempdf$formula<-peptides_symbol_adducts[[i]]
    message(paste("Calculating peptide mz with Decoy_adducts:",Decoy_adducts[i]))
    #3templist=bplapply(tempdf$formula,function(x,charge){
    #  rcdk::get.formula(x,charge = charge)@mass
    #  },charge,BPPARAM = BPPARAM)
    if (charge==0){actingcharge=1} else {actingcharge=abs(charge)}
    #tempdf$mz<-as.numeric(unlist(templist))
    tempdf$mz<-(tempdf$pepmz+adductmass)/actingcharge
    tempdf$adduct<-Decoy_adducts[i]
    tempdf$isdecoy<-rep(1,nrow(tempdf))
    tempdf$charge<-charge
    #tempdf$isdecoy<-rep(0,nrow(tempdf))
    #convert_peptide_adduct(tempdf$Peptide[2],adductsname = adducts[i],multiplier = c(multiplier,1),adductslist = adductslist)
    #tempdf$formula<-bplapply(tempdf$Peptide,convert_peptide_adduct,adductsformula = adducts[i],multiplier = c(multiplier,1),BPPARAM = BPPARAM)
    Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
  } 
  }
  if (Decoy_search && ("elements" %in% Decoy_mode)){
    
    message(paste("Generating peptide formula with Decoy elemental composition."))
    library(enviPat)
    library(rcdk)
    total_target_mols<-unique(Protein_Summary$formula)

    target_mz=unique(round(as.numeric(Protein_Summary$mz),digits = 3))
    select_mol_num=ceiling(length(total_target_mols)/length(target_mz))
    
    decoymol<-lapply(target_mz,function(x,ppm,select_mol_num){
      require(rcdk,quietly = T)
     mit <- generate.formula.iter(x, window = x*ppm/(1000000),
                          elements = list(
                            C=c(0,100),
                            H=c(0,100),
                            N=c(0,100),
                            O=c(0,100),
                            S=c(0,100)),
                          validation = FALSE,
                          charge = 0.0,
                          as.string=TRUE)

    i=1
    moleresult<-as.character()
    hit <- itertools::ihasNext(mit)
    while (`&`(itertools::hasNext(hit),i<=select_mol_num))  {
      moleresult<-c(moleresult,iterators::nextElem(hit))
      
      i=i+1
    }  
    #message(x)
    return(moleresult)
    
    },ppm,select_mol_num)
    
    names(decoymol)=as.character(target_mz)
    
    #saveRDS(decoymol, file = "decoymol.rds")
    
    #decoymol<-readRDS(file = "decoymol.rds")
    
    decoymolvec<-do.call(c,decoymol)
    
    length_decoydf<-lapply(decoymol,length)
    
    decoy_mz<-names(length_decoydf)
    
    decoymol<-decoymol[!(decoymol %in% Protein_Summary$formula)]
    
    decoy_tempdf=data.frame(stringsAsFactors = F,
                            Protein=paste("Decoy_protein",sep=""),
                            Peptide=paste("Decoy_",rep(decoy_mz,length_decoydf),sep=""),
                            pepmz=rep(decoy_mz,length_decoydf),
                            formula=decoymolvec,mz=rep(decoy_mz,length_decoydf),adduct="M",
                            isdecoy=1,charge=1)
    
    Protein_Summary<-rbind(Protein_Summary,decoy_tempdf)
    
  }
  
  #Protein_Summary$Protein=Index_of_protein_sequence[Index_of_protein_sequence$desc==Protein_Summary$Protein]
  #temp_index=Index_of_protein_sequence
  #temp_index$Protein=temp_index$desc
  #temp_index=temp_index[,c("Protein","recno")]
  #Protein_Summary=merge(Protein_Summary,temp_index,by="Protein",all.x=T)
  #Protein_Summary$Protein=Protein_Summary$recno
  if(output_candidatelist){
    if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
    write.csv(Protein_Summary,paste(workdir,"/Summary folder/candidatelist.csv",sep=""),row.names = F)
    write.csv(Index_of_protein_sequence,paste(workdir,"/Summary folder/protein_index.csv",sep=""),row.names = F)
    message("Candidate list has been exported.")
  }
  }

  if(use_previous_candidates){
    #if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
    if (sum(c("candidatelist.csv","protein_index.csv") %in% dir(paste(workdir,"/Summary folder",sep="")))==2){
     Protein_Summary<-read.csv(paste(workdir,"/Summary folder/candidatelist.csv",sep=""),row.names = F)
    Index_of_protein_sequence<-read.csv(paste(workdir,"/Summary folder/protein_index.csv",sep=""),row.names = F)
    message("Candidate list has been loaded.") 
    }else{stop("Can not find the previously established candidate list.")}
    
  }
  return(Protein_Summary)
}


vendiagram<-function(){
  library(stringr)
  library(tcltk)
  folder1<-tkchooseDirectory()
  folder2<-tkchooseDirectory()
  Protein_feature_summary_uniport <- read.csv(paste0(as.character(paste0(folder1,collapse = " ")),"/Summary folder/Protein_feature_summary.csv"))
  Protein_feature_summary_marine <- read.csv(paste0(as.character(paste0(folder2,collapse = " ")),"/Summary folder/Protein_feature_summary.csv"))
  #Protein_feature_summary_nonmarine <- read.csv("D:/Tumour test/test9 new matrisome/Summary folder/Protein_feature_summary.csv")
  
  #Protein_feature_summary_uniport <- read.csv("D:/Tumour test/test12 new non matrisome gnil/Summary folder/Protein_feature_summary.csv")
  #Protein_feature_summary_marine <- read.csv("D:/Tumour test/test11 new matrisome gnil/Summary folder/Protein_feature_summary.csv")
  #Protein_feature_summary_nonmarine <- read.csv("D:/Tumour test/test9 new matrisome/Summary folder/Protein_feature_summary.csv")
  
  result<-list()
  maxinte=min(max(Protein_feature_summary_marine$Intensity),max(Protein_feature_summary_uniport$Intensity))
  uniport_peptides0.05<-unique(Protein_feature_summary_uniport$Peptide[Protein_feature_summary_uniport$Intensity>=maxinte*0.05])
  matrisome_peptides0.05<-unique(Protein_feature_summary_marine$Peptide[Protein_feature_summary_marine$Intensity>=maxinte*0.05])
  result[["peptides0.05"]][1]<-length(matrisome_peptides0.05)
  result[["peptides0.05"]][2]<-length(uniport_peptides0.05)
  result[["peptides0.05"]][3]<-sum(matrisome_peptides0.05 %in% uniport_peptides0.05)
  
  uniport_peptides0.005<-unique(Protein_feature_summary_uniport$Peptide[Protein_feature_summary_uniport$Intensity>=maxinte*0.005])
  matrisome_peptides0.005<-unique(Protein_feature_summary_marine$Peptide[Protein_feature_summary_marine$Intensity>=maxinte*0.005])
  result[["peptides0.005"]][1]<-length(matrisome_peptides0.005)
  result[["peptides0.005"]][2]<-length(uniport_peptides0.005)
  result[["peptides0.005"]][3]<-sum(matrisome_peptides0.005 %in% uniport_peptides0.005)
  
  uniport_mz0.05<-unique(Protein_feature_summary_uniport$mz[Protein_feature_summary_uniport$Intensity>=maxinte*0.05])
  matrisome_mz0.05<-unique(Protein_feature_summary_marine$mz[Protein_feature_summary_marine$Intensity>=maxinte*0.05])
  result[["mz0.05"]][1]<-length(matrisome_mz0.05)
  result[["mz0.05"]][2]<-length(uniport_mz0.05)
  result[["mz0.05"]][3]<-sum(uniport_mz0.05 %in% matrisome_mz0.05)
  
  uniport_mz0.005<-unique(Protein_feature_summary_uniport$mz[Protein_feature_summary_uniport$Intensity>=maxinte*0.005])
  matrisome_mz0.005<-unique(Protein_feature_summary_marine$mz[Protein_feature_summary_marine$Intensity>=maxinte*0.005])
  result[["mz0.005"]][1]<-length(matrisome_mz0.005)
  result[["mz0.005"]][2]<-length(uniport_mz0.005)
  result[["mz0.005"]][3]<-sum(uniport_mz0.005 %in% matrisome_mz0.005)
  
  uniport_prot<-unique(Protein_feature_summary_uniport$Protein[Protein_feature_summary_uniport$Intensity>=maxinte*0.05])
  matrisome_prot<-unique(Protein_feature_summary_marine$Protein[Protein_feature_summary_marine$Intensity>=maxinte*0.05])
  
  
  uniport_prot_a=str_extract(uniport_prot,"\\|.{2,}\\|")
  matrisome_prot_a=str_extract(matrisome_prot,"\\|.{2,}\\|")
  result[["prot0.05"]][2]<-length(matrisome_prot)
  result[["prot0.05"]][1]<-length(uniport_prot)
  result[["prot0.05"]][3]<-sum(uniport_prot %in% matrisome_prot)
  result<-data.frame(result,stringsAsFactors = F)
  return(result)
  #protein_unique=unique(Protein_feature_list$mz)
}

fastafile_utils<-function(){
  
  library(Biostrings)
  library(cleaver)
  library(protViz)
  library(rcdk)
  
  list_of_protein_sequence<-readAAStringSet(database,
                                            format="fasta",
                                            nrec=-1L, 
                                            skip=0L, 
                                            seek.first.rec=FALSE,
                                            use.names=TRUE, 
                                            with.qualities=FALSE) 
  
  Index_of_protein_sequence_uniport<-fasta.index("D:/Tumour test/test7/uniprot.fasta",
                                         nrec=-1L, 
                                         skip=0L) 
  Index_of_protein_sequence_matrisome<-fasta.index("D:/Tumour test/test8_matrisome/murine_matrisome.fasta",
                                                 nrec=-1L, 
                                                 skip=0L)
  
  Index_of_protein_sequence_uniport$ID<-str_extract(Index_of_protein_sequence_uniport$desc,"\\|.{2,}\\|")
  
  Index_of_protein_sequence_matrisome$ID<-str_extract(Index_of_protein_sequence_matrisome$desc,"\\|.{2,}\\|")
  
  overlapID<-intersect(Index_of_protein_sequence_uniport$ID,Index_of_protein_sequence_matrisome$ID)

  uniquematrisomeID<-Index_of_protein_sequence_matrisome[!(Index_of_protein_sequence_matrisome$ID %in% overlapID),]
  
  matrisome_def<-fread(file="matrisome_mm_masterlist.csv",header = T)
  
  matrisome_def_id<-paste(matrisome_def$UniProt_IDs, collapse = ":")
  
  matrisome_def_id<-(str_split(matrisome_def_id,":"))[[1]]
  
  matrisome_def_id<-paste0("|",matrisome_def_id,"|")
  
  sum(matrisome_def_id %in% Index_of_protein_sequence_uniport$ID)
  
  sum(!(matrisome_def_id %in% Index_of_protein_sequence_uniport$ID))
  
  matrisome_def_unique_id<-matrisome_def_id[!(matrisome_def_id %in% Index_of_protein_sequence_uniport$ID)]
  
  matrisome_def_overlap_id<-matrisome_def_id[(matrisome_def_id %in% Index_of_protein_sequence_uniport$ID)]
  
  Index_of_protein_sequence_matrisome_subset<-Index_of_protein_sequence_uniport[Index_of_protein_sequence_uniport$ID %in% matrisome_def_overlap_id,]
  
  Index_of_protein_sequence_nonmatrisome_subset<-Index_of_protein_sequence_uniport[!(Index_of_protein_sequence_uniport$ID %in% matrisome_def_overlap_id),]
  
  list_of_matrisome<-readAAStringSet(Index_of_protein_sequence_matrisome_subset) 
  
  list_of_nonmatrisome<-readAAStringSet(Index_of_protein_sequence_nonmatrisome_subset)
  
  
  
  writeXStringSet(list_of_matrisome,"matrisome.fasta")
  
  writeXStringSet(list_of_nonmatrisome,"non-matrisome.fasta")
  
  sum(!(matrisome_def_unique_id %in% uniquematrisomeID$ID))
  
  sum(!(uniquematrisomeID$ID %in% matrisome_def_unique_id))
  
}

convert_peptide_adduct<-function(peptide,adductsname,multiplier=c(1,1),adductslist=Build_adduct_list()){
  library(rcdk)
  library(OrgMassSpecR)
  ConvertPeptide<-function (sequence, output = "elements") {
    peptideVector <- strsplit(sequence, split = "")[[1]]
    if (output == "elements") {
      FindElement <- function(residue) {
        element<-c(C = 0, H = 0, N = 0, O = 0, S = 0)
        switch(residue,
               A={element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)},
               R={element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)},
               N={element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)},
               D={element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)},
               E={element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)},
               Q={element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)},
               G={element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)},
               H={element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)},
               I={element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)},
               L={element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)},
               K={element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)},
               M={element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)},
               F={element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)},
               P={element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)},
               S={element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)},
               T={element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)},
               W={element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)},
               Y={element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)},
               V={element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)},
               C={element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)})
        
        return(element)
      }
      resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
      for (i in 1:length(peptideVector)) {
        resultsVector <- FindElement(peptideVector[i]) + 
          resultsVector
      }
      resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, 
                                         O = 1, S = 0)
      return(as.list(resultsVector))
    }
    if (output == "3letter") {
      FindCode <- function(residue) {
        if (residue == "A") 
          let <- "Ala"
        if (residue == "R") 
          let <- "Arg"
        if (residue == "N") 
          let <- "Asn"
        if (residue == "D") 
          let <- "Asp"
        if (residue == "C") 
          let <- "Cys"
        if (residue == "E") 
          let <- "Glu"
        if (residue == "Q") 
          let <- "Gln"
        if (residue == "G") 
          let <- "Gly"
        if (residue == "H") 
          let <- "His"
        if (residue == "I") 
          let <- "Ile"
        if (residue == "L") 
          let <- "Leu"
        if (residue == "K") 
          let <- "Lys"
        if (residue == "M") 
          let <- "Met"
        if (residue == "F") 
          let <- "Phe"
        if (residue == "P") 
          let <- "Pro"
        if (residue == "S") 
          let <- "Ser"
        if (residue == "T") 
          let <- "Thr"
        if (residue == "W") 
          let <- "Trp"
        if (residue == "Y") 
          let <- "Tyr"
        if (residue == "V") 
          let <- "Val"
        return(let)
      }
      codes <- sapply(peptideVector, FindCode)
      return(paste(codes, collapse = ""))
    }
  }
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
  get_atoms<-function(Symbol){
    #form = "C5H11BrO" 
    if (Symbol!=""){
      ups = c(gregexpr("[[:upper:]]", Symbol)[[1]], nchar(Symbol) + 1) 
      seperated = sapply(1:(length(ups)-1), function(x) substr(Symbol, ups[x], ups[x+1] - 1)) 
      elements =  gsub("[[:digit:]]", "", seperated) 
      nums = gsub("[[:alpha:]]", "", seperated) 
      ans = data.frame(elements = as.character(elements), num = as.numeric(ifelse(nums == "", 1, nums)), stringsAsFactors = FALSE)
      list<-as.list(ans$num)
      names(list)=ans$elements
      return(list)
    }else{return(NULL)}
  }
  
  
  adductsformula_add= as.character(adductslist[adductslist$Name==adductsname,"Formula_add"])
  adductsformula_ded= as.character(adductslist[adductslist$Name==adductsname,"Formula_ded"])
  
  null_if_false<-function(x){if (x==FALSE){""} else{get.formula(x)@string}}
  
  adductsformula_add<-null_if_false(adductsformula_add)
  adductsformula_ded<-null_if_false(adductsformula_ded)
  
  adductsformula_add<-get_atoms(adductsformula_add)
  adductsformula_ded<-get_atoms(adductsformula_ded)
  
  adductsformula<-merge_atoms(adductsformula_add,adductsformula_ded,check_merge = F,mode = "ded")
  
  formula<-ConvertPeptide(peptide)
  
  formula_with_adducts<-merge_atoms(formula,adductsformula,check_merge = T,mode = "add", multiplier = multiplier)
  
  for (name in names(formula_with_adducts)){
    if(formula_with_adducts[[name]]==0){ formula_with_adducts[[name]]=NULL}
  }

  return(paste0(names(formula_with_adducts),formula_with_adducts,collapse = ""))
  
}

convert_peptide_adduct_list<-function(adductsname,peptide_symbol,multiplier=c(1,1),adductslist=Build_adduct_list()){
  library(rcdk)
  library(OrgMassSpecR)
  ConvertPeptide<-function (sequence, output = "elements") {
    peptideVector <- strsplit(sequence, split = "")[[1]]
    if (output == "elements") {
      FindElement <- function(residue) {
        element<-c(C = 0, H = 0, N = 0, O = 0, S = 0)
        switch(residue,
               A={element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)},
               R={element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)},
               N={element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)},
               D={element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)},
               E={element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)},
               Q={element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)},
               G={element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)},
               H={element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)},
               I={element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)},
               L={element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)},
               K={element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)},
               M={element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)},
               F={element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)},
               P={element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)},
               S={element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)},
               T={element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)},
               W={element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)},
               Y={element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)},
               V={element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)},
               C={element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)},
               Z={element <- c(C = 3, H = 3, N = 1, O = 1, S = 1)})
        
        return(element)
      }
      resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
      for (i in 1:length(peptideVector)) {
        resultsVector <- FindElement(peptideVector[i]) + 
          resultsVector
      }
      resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, 
                                         O = 1, S = 0)
      return(as.list(resultsVector))
    }
    if (output == "3letter") {
      FindCode <- function(residue) {
        if (residue == "A") 
          let <- "Ala"
        if (residue == "R") 
          let <- "Arg"
        if (residue == "N") 
          let <- "Asn"
        if (residue == "D") 
          let <- "Asp"
        if (residue == "C") 
          let <- "Cys"
        if (residue == "E") 
          let <- "Glu"
        if (residue == "Q") 
          let <- "Gln"
        if (residue == "G") 
          let <- "Gly"
        if (residue == "H") 
          let <- "His"
        if (residue == "I") 
          let <- "Ile"
        if (residue == "L") 
          let <- "Leu"
        if (residue == "K") 
          let <- "Lys"
        if (residue == "M") 
          let <- "Met"
        if (residue == "F") 
          let <- "Phe"
        if (residue == "P") 
          let <- "Pro"
        if (residue == "S") 
          let <- "Ser"
        if (residue == "T") 
          let <- "Thr"
        if (residue == "W") 
          let <- "Trp"
        if (residue == "Y") 
          let <- "Tyr"
        if (residue == "V") 
          let <- "Val"
        return(let)
      }
      codes <- sapply(peptideVector, FindCode)
      return(paste(codes, collapse = ""))
    }
  }
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
  
  get_atoms<-function(Symbol){
    #form = "C5H11BrO" 
    if (Symbol!=""){
      ups = c(gregexpr("[[:upper:]]", Symbol)[[1]], nchar(Symbol) + 1) 
      seperated = sapply(1:(length(ups)-1), function(x) substr(Symbol, ups[x], ups[x+1] - 1)) 
      elements =  gsub("[[:digit:]]", "", seperated) 
      nums = gsub("[[:alpha:]]", "", seperated) 
      ans = data.frame(elements = as.character(elements), num = as.numeric(ifelse(nums == "", 1, nums)), stringsAsFactors = FALSE)
      list<-as.list(ans$num)
      names(list)=ans$elements
      return(list)
    }else{return(NULL)}
  }
   if (missing(multiplier)){
     multiplier[1]= as.numeric(as.character(adductslist[adductslist$Name==adductsname,"Mult"]))
     multiplier[2]=1
   }
  adductsformula_add= as.character(adductslist[adductslist$Name==adductsname,"Formula_add"])
  adductsformula_ded= as.character(adductslist[adductslist$Name==adductsname,"Formula_ded"])
  
  null_if_false<-function(x){if (x==FALSE){""} else{get.formula(x)@string}}
  
  adductsformula_add<-null_if_false(adductsformula_add)
  adductsformula_ded<-null_if_false(adductsformula_ded)
  
  adductsformula_add<-get_atoms(adductsformula_add)
  adductsformula_ded<-get_atoms(adductsformula_ded)
  
  adductsformula<-merge_atoms(adductsformula_add,adductsformula_ded,check_merge = F,mode = "ded")
  peptide_formula<-as.character()
  #formula<-ConvertPeptide(peptide)
  for (formula in 1:length(peptide_symbol)){
   formula_with_adducts<-merge_atoms(peptide_symbol[[formula]],adductsformula,check_merge = T,mode = "add", multiplier = multiplier) 
   for (name in names(formula_with_adducts)){
    if(formula_with_adducts[[name]]==0){ formula_with_adducts[[name]]=NULL}
   }
   peptide_formula[[formula]]=(paste0(names(formula_with_adducts),formula_with_adducts,collapse = ""))
  }
  
  return(peptide_formula)
  
}
#' get_atoms
#'
#' This is a function that prepare the candiate list for maldi imaging data qualitative or quantitative analysis.
#' this function will read the candidate list file and generate mz for the adducts list defined in \code{"adducts"}. 
#' @param database the file name of candidate list
#' @param workdir the folder that contains candidate list
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param cal.mz If set with \code{"TRUE"}, the function will recalculate the mz value according to the column named "Formula" in the \code{database} and the specified adducts.
#' @param mzlist_bypass  Set \code{"TRUE"} if you want to bypass the mzlist generating process, the function will keep the mz and adduct as it is for the furture analysis. Be sure that the table contains "mz", "adduct" and "moleculeNames" as they are essential for later steps.
#' @param BPPARAM parallel processing parameter for BiocParallel
#' @return a table of candiate list
#'
#' @examples
#' Meta_feature_list_fun(database="lipid candidates.csv",adducts=c("M-H","M+Cl"))
#'
#' @export
#' 
get_atoms<-function(Symbol){
  #form = "C5H11BrO" 
  if (Symbol!=""){
      ups = c(gregexpr("[[:upper:]]", Symbol)[[1]], nchar(Symbol) + 1) 
  seperated = sapply(1:(length(ups)-1), function(x) substr(Symbol, ups[x], ups[x+1] - 1)) 
  elements =  gsub("[[:digit:]]", "", seperated) 
  nums = gsub("[[:alpha:]]", "", seperated) 
  ans = data.frame(element = as.character(elements), num = as.numeric(ifelse(nums == "", 1, nums)), stringsAsFactors = FALSE)
  list<-as.list(ans$num)
  names(list)=ans$element
  return(list)
  }else{return(NULL)}
}


#' merge_atoms
#'
#' This is a function that prepare the candiate list for maldi imaging data qualitative or quantitative analysis.
#' this function will read the candidate list file and generate mz for the adducts list defined in \code{"adducts"}. 
#' @param database the file name of candidate list
#' @param workdir the folder that contains candidate list
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param cal.mz If set with \code{"TRUE"}, the function will recalculate the mz value according to the column named "Formula" in the \code{database} and the specified adducts.
#' @param mzlist_bypass  Set \code{"TRUE"} if you want to bypass the mzlist generating process, the function will keep the mz and adduct as it is for the furture analysis. Be sure that the table contains "mz", "adduct" and "moleculeNames" as they are essential for later steps.
#' @param BPPARAM parallel processing parameter for BiocParallel
#' @return a table of candiate list
#'
#' @examples
#' Meta_feature_list_fun(database="lipid candidates.csv",adducts=c("M-H","M+Cl"))
#'
#' @export
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



ConvertPeptide<-function (sequence, output = "elements") 
{
  peptideVector <- strsplit(sequence, split = "")[[1]]
  if (output == "elements") {
    FindElement <- function(residue) {
      element<-c(C = 0, H = 0, N = 0, O = 0, S = 0)
      switch(residue,
             A={element <- c(C = 3, H = 5, N = 1, O = 1, S = 0)},
             R={element <- c(C = 6, H = 12, N = 4, O = 1, S = 0)},
             N={element <- c(C = 4, H = 6, N = 2, O = 2, S = 0)},
             D={element <- c(C = 4, H = 5, N = 1, O = 3, S = 0)},
             E={element <- c(C = 5, H = 7, N = 1, O = 3, S = 0)},
             Q={element <- c(C = 5, H = 8, N = 2, O = 2, S = 0)},
             G={element <- c(C = 2, H = 3, N = 1, O = 1, S = 0)},
             H={element <- c(C = 6, H = 7, N = 3, O = 1, S = 0)},
             I={element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)},
             L={element <- c(C = 6, H = 11, N = 1, O = 1, S = 0)},
             K={element <- c(C = 6, H = 12, N = 2, O = 1, S = 0)},
             M={element <- c(C = 5, H = 9, N = 1, O = 1, S = 1)},
             F={element <- c(C = 9, H = 9, N = 1, O = 1, S = 0)},
             P={element <- c(C = 5, H = 7, N = 1, O = 1, S = 0)},
             S={element <- c(C = 3, H = 5, N = 1, O = 2, S = 0)},
             T={element <- c(C = 4, H = 7, N = 1, O = 2, S = 0)},
             W={element <- c(C = 11, H = 10, N = 2, O = 1, S = 0)},
             Y={element <- c(C = 9, H = 9, N = 1, O = 2, S = 0)},
             V={element <- c(C = 5, H = 9, N = 1, O = 1, S = 0)},
             C={element <- c(C = 3, H = 5, N = 1, O = 1, S = 1)},
             Z={element <- c(C = 3, H = 3, N = 1, O = 1, S = 1)})

      return(element)
    }
    resultsVector <- c(C = 0, H = 0, N = 0, O = 0, S = 0)
    for (i in 1:length(peptideVector)) {
      resultsVector <- FindElement(peptideVector[i]) + 
        resultsVector
    }
    resultsVector <- resultsVector + c(C = 0, H = 2, N = 0, 
                                       O = 1, S = 0)
    return(as.list(resultsVector))
  }
  if (output == "3letter") {
    FindCode <- function(residue) {
      if (residue == "A") 
        let <- "Ala"
      if (residue == "R") 
        let <- "Arg"
      if (residue == "N") 
        let <- "Asn"
      if (residue == "D") 
        let <- "Asp"
      if (residue == "C") 
        let <- "Cys"
      if (residue == "E") 
        let <- "Glu"
      if (residue == "Q") 
        let <- "Gln"
      if (residue == "G") 
        let <- "Gly"
      if (residue == "H") 
        let <- "His"
      if (residue == "I") 
        let <- "Ile"
      if (residue == "L") 
        let <- "Leu"
      if (residue == "K") 
        let <- "Lys"
      if (residue == "M") 
        let <- "Met"
      if (residue == "F") 
        let <- "Phe"
      if (residue == "P") 
        let <- "Pro"
      if (residue == "S") 
        let <- "Ser"
      if (residue == "T") 
        let <- "Thr"
      if (residue == "W") 
        let <- "Trp"
      if (residue == "Y") 
        let <- "Tyr"
      if (residue == "V") 
        let <- "Val"
      return(let)
    }
    codes <- sapply(peptideVector, FindCode)
    return(paste(codes, collapse = ""))
  }
}
