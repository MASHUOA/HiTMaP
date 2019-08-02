
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
                                   Decoy_adducts=c("M+He","M+Li","M+Cu","M+Co","M+Ag","M+Ar")){
  library(Biostrings)
  library(cleaver)
  library(protViz)
  library(rcdk)
  library(BiocParallel)
  Decoy_adducts=Decoy_adducts[1:length(adducts)]
  list_of_protein_sequence<-readAAStringSet(database,
                                            format="fasta",
                                            nrec=-1L, 
                                            skip=0L, 
                                            seek.first.rec=FALSE,
                                            use.names=TRUE, 
                                            with.qualities=FALSE)    
 
  
  
  
  Index_of_protein_sequence<-fasta.index(database,
                                         nrec=-1L, 
                                         skip=0L)   
  
  peplist<-cleave(as.character(list_of_protein_sequence),custom=Digestion_site, missedCleavages=missedCleavages)
  
  
  peplist<-peplist[duplicated(names(peplist))==FALSE]
  
  
  
  pimlist<-parentIonMasslist(peplist,Index_of_protein_sequence)
  
  peplist<-peplist[names(peplist) %in% names(pimlist) ==TRUE] 
  
  #tempdf<-parLapply(cl=cl,  1: length(names(peplist)), Peptide_Summary_para,peplist)
  #bplapply()
  tempdf<-bplapply( 1: length(names(peplist)), Peptide_Summary_para,peplist,BPPARAM = BPPARAM)
  tempdf <- do.call("rbind", tempdf)
  colnames(tempdf)<-c("Protein","Peptide")
  tempdf<-as.data.frame(tempdf)
  tempdf$Protein<-as.character(tempdf$Protein)
  tempdf$Peptide<-as.character(tempdf$Peptide)
  tempdf1<-tempdf
  Protein_Summary<-NULL
  adductslist<-Build_adduct_list()
  for (i in 1:length(adducts)){
    adductmass<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Mass"]))
    multiplier<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Multi"]))
    adductsformula<-c(as.character(adductslist[adductslist$Name==adducts[i],"Formula_add"]),
    as.character(adductslist[adductslist$Name==adducts[i],"Formula_ded"]))
    adductsformula=paste0(ifelse(adductsformula=="FALSE","",adductsformula),collapse = "")
    
    tempdf$mz<-as.numeric(parentIonMass(tempdf$Peptide)-1.007825+adductmass)
    tempdf$adduct<-adducts[i]
    #tempdf$isdecoy<-rep(0,nrow(tempdf))
    #convert_peptide_adduct(tempdf$Peptide[2],adductsname = adducts[i],multiplier = c(multiplier,1),adductslist = adductslist)
    #tempdf$formula<-bplapply(tempdf$Peptide,convert_peptide_adduct,adductsformula = adducts[i],multiplier = c(multiplier,1),BPPARAM = BPPARAM)
    Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
  }
  
  Protein_Summary
}

Protein_feature_list_fun1<-function(workdir=getwd(),
                                   database,
                                   Digestion_site="[G]",
                                   missedCleavages=0:1,
                                   adducts=c("M+H","M+NH4","M+Na"),
                                   BPPARAM=bpparam(),
                                   Decoy_adducts=c("M+He","M+Li","M+Cu","M+Co","M+Ag","M+Ar")){
  library(Biostrings)
  library(cleaver)
  library(protViz)
  library(rcdk)
  library(BiocParallel)
  Decoy_adducts=Decoy_adducts[1:length(adducts)]
  list_of_protein_sequence<-readAAStringSet(database,
                                            format="fasta",
                                            nrec=-1L, 
                                            skip=0L, 
                                            seek.first.rec=FALSE,
                                            use.names=TRUE, 
                                            with.qualities=FALSE)    
  
  
  
  
  Index_of_protein_sequence<-fasta.index(database,
                                         nrec=-1L, 
                                         skip=0L)   
  
  peplist<-cleave(as.character(list_of_protein_sequence),custom=Digestion_site, missedCleavages=missedCleavages)
  
  
  peplist<-peplist[duplicated(names(peplist))==FALSE]
  
  
  
  pimlist<-parentIonMasslist(peplist,Index_of_protein_sequence)
  
  peplist<-peplist[names(peplist) %in% names(pimlist) ==TRUE] 
  
  #tempdf<-parLapply(cl=cl,  1: length(names(peplist)), Peptide_Summary_para,peplist)
  #bplapply()
  tempdf<-bplapply( 1: length(names(peplist)), Peptide_Summary_para,peplist,BPPARAM = BPPARAM)
  tempdf <- do.call("rbind", tempdf)
  colnames(tempdf)<-c("Protein","Peptide")
  tempdf<-as.data.frame(tempdf)
  tempdf$Protein<-as.character(tempdf$Protein)
  tempdf$Peptide<-as.character(tempdf$Peptide)
  tempdf1<-tempdf
  Protein_Summary<-NULL
  adductslist<-Build_adduct_list()
  for (i in 1:length(adducts)){
    adductmass<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Mass"]))
    multiplier<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Multi"]))
    
    adductsformula<-c(as.character(adductslist[adductslist$Name==adducts[i],"Formula_add"]),
                      as.character(adductslist[adductslist$Name==adducts[i],"Formula_ded"]))
    
    adductsformula=paste0(ifelse(adductsformula=="FALSE","",adductsformula),collapse = "")
    
    tempdf$mz<-as.numeric(parentIonMass(tempdf$Peptide)-1.007825+adductmass)
    tempdf$adduct<-adducts[i]
    tempdf$isdecoy<-rep(0,nrow(tempdf))
    #convert_peptide_adduct(tempdf$Peptide[2],adductsname = adducts[i],multiplier = c(multiplier,1),adductslist = adductslist)
    message(paste("Building candidates list for adduct",adducts[i])) 
    tempdf$formula<-bplapply(tempdf$Peptide,convert_peptide_adduct,adductsname = adducts[i],adductslist=adductslist,multiplier = c(multiplier,1),BPPARAM = BPPARAM)
    Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
  }
  tempdf<-tempdf1
  for (i in 1:length(Decoy_adducts)){
    tempdf<-tempdf1
    adductmass<-as.numeric(as.character(adductslist[adductslist$Name==Decoy_adducts[i],"Mass"]))
    tempdf$mz<-as.numeric(parentIonMass(tempdf$Peptide)-1.007825+adductmass)
    tempdf$isdecoy<-rep(1,nrow(tempdf))
    tempdf$adduct<-Decoy_adducts[i]
    Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
  }
  Protein_Summary
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
  
  
  adductsformula_add= as.character(adductslist[adductslist$Name==adductsname,"Formula_add"])
  adductsformula_ded= as.character(adductslist[adductslist$Name==adductsname,"Formula_ded"])
  
  null_if_false<-function(x){if (x==FALSE){""} else{get.formula(x)@string}}
  
  adductsformula_add<-null_if_false(adductsformula_add)
  adductsformula_ded<-null_if_false(adductsformula_ded)
  
  adductsformula_add<-get_atoms(adductsformula_add)
  adductsformula_ded<-get_atoms(adductsformula_ded)
  
  adductsformula<-merge_atoms(adductsformula_add,adductsformula_ded,check_merge = F,mode = "ded")
  
  formula<-ConvertPeptide(peptide,IAA = F)
  
  formula_with_adducts<-merge_atoms(formula,adductsformula,check_merge = T,mode = "add", multiplier = multiplier)
  
  for (name in names(formula_with_adducts)){
    if(formula_with_adducts[[name]]==0){ formula_with_adducts[[name]]=NULL}
  }

  return(paste0(names(formula),formula,collapse = ""))
  
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
