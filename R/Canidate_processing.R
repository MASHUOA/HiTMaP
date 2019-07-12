
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
                                   BPPARAM=bpparam()){
  library(Biostrings)
  library(cleaver)
  library(protViz)
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
  
  Protein_Summary<-NULL
  adductslist<-Build_adduct_list()
  for (i in 1:length(adducts)){
    adductmass<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Mass"]))
    tempdf$mz<-as.numeric(parentIonMass(tempdf$Peptide)-1.00727600+adductmass)
    tempdf$adduct<-adducts[i]
    Protein_Summary<-rbind.data.frame(Protein_Summary,tempdf)
  }
  Protein_Summary
}