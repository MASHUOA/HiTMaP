#' imaging_Spatial_Quant
#'
#' This is a spatial quantitation function for maldi imaging data set
#' this function will read the candidate list file and generate quantification result
#' @param workdir the 
#' @param Quant_list the file path of candidate list, if 
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param cal.mz If set with \code{"true"}, the function will recalculate the mz value according to the column named "formular" in the \code{Quant_list} and the specified adducts.
#' @param mzlist_bypass  Set \code{"true"} if you want to bypass the mzlist generating process
#' @param Protein_feature_summary  \code{"PMF_analysis"} follow-up process that will collect all the identified peptide information and associate them with possible proteins 
#' @param plot_cluster_image  \code{"Protein_feature_summary"} follow-up process that will plot the protein cluster image 
#' 
#' @param Peptide_feature_summarya \code{"PMF_analysis"} follow-up process that will summarize all datafiles identified peptides and generats a \code{"peptide shortlist"} in the result summary folder
#' @param plot_ion_image  \code{"Peptide_feature_summarya"} follow-up process that will plot every connponents in the \code{"peptide shortlist"}
#' @param parallel the number of threads will be used in the PMF search, this option now only works for windows OS
#' @param spectra_segments_per_file optimal number of distinctive regions in the imaging, a virtual segmentation will be applied to the image files with this value. To have a better PMF result you may set a value that in the sweet point of sensitivety and false discovery rate (FDR).
#' @param spatialKMeans set true to enable a \code{"spatialKMeans"}  method for the automatic virtual segmentation. If a region rank file was supplied, you can disable this to perform a mannual segmentation.
#' @param Smooth_range \code{"spatialKMeans"} pixel smooth range 
#' @param Virtual_segmentation set \code{"TRUE"} if you want to overide the automaitic segmentation
#' @param Virtual_segmentation_rankfile specify a region rank file contains region information for manualy region segmentation
#' @return None
#'
#' @examples
#' imaging_Spatial_Quant(threshold=0.05, ppm=5,Digestion_site="[G]",
#'                        missedCleavages=0:1,Fastadatabase="murine_matrisome.fasta",
#'                        adducts=c("M+H","M+NH4","M+Na"),PMF_analysis=TRUE,
#'                        Protein_feature_summary=TRUE,plot_cluster_image=TRUE,
#'                        Peptide_feature_summary=TRUE,plot_ion_image=FALSE,
#'                        parallel=3,spectra_segments_per_file=5,spatialKMeans=TRUE
#'                        )
#'
#' @export

Meta_feature_list_fun<-function(database,
                                workdir=getwd(),
                                mode=c("molecule","peptide"),
                                adducts=adducts,
                                cl=autoStopCluster(makeCluster(parallel)),
                                cal.mz=F,
                                bypass=F){
  library("Rcpp")
  library(dplyr)
  library(Rdisop)
  library(Biostrings)
  
  adductslist<-Build_adduct_list()
  candidates<-read.csv(paste0(workdir,"/",Quant_list),as.is = TRUE)
  
  if (bypass){Meta_Summary=candidates[is.na(candidates$mz)!=T,]}else{
    
    if (mode=="molecule"){
       candidates$mass=0
    if (cal.mz==F){
      if (sum(candidates$adducts!=0)>0){
        candidates$adducts_mass=0
        candidates$adducts_mass=unlist(parallel::parLapply(cl,1:nrow(candidates),function(i,adducts,adductslist,adducts_mass){
          adducts_mass[i]<-adductslist[adductslist$Name==adducts[i],"Mass"]*abs(as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Charge"])))
        },candidates$adducts,adductslist,candidates$adducts_mass))
        candidates$mass<-candidates$mz-candidates$adducts_mass
      }else{
        candidates$mass<-candidates$mz
        }
      }else{
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
    } else if (mode=="peptide"){
      
      
      
      
    }
    
    
    
  }
  
  

  return(Meta_Summary)                
  
} 

Protein_feature_list_fun<-function(workdir=WorkingDir(),
                                     Fastadatabase,
                                     format="fasta",
                                     nrec=-1L, 
                                     skip=0L, 
                                     seek.first.rec=FALSE,
                                     use.names=TRUE, 
                                     with.qualities=FALSE,
                                     Digestion_site="[G]",
                                     missedCleavages=0:1,
                                     adducts=c("M+H","M+NH4","M+Na"),
                                     cl){
    library(Biostrings)
    list_of_protein_sequence<-readAAStringSet(Fastadatabase,
                                              format,
                                              nrec, 
                                              skip, 
                                              seek.first.rec,
                                              use.names, 
                                              with.qualities
    )    
    
    Index_of_protein_sequence<-fasta.index(Fastadatabase,
                                           nrec, 
                                           skip)   
    
    peplist<-cleave(as.character(list_of_protein_sequence),custom=Digestion_site, missedCleavages=missedCleavages)
    
    
    peplist<-peplist[duplicated(names(peplist))==FALSE]
    
    
    
    pimlist<-parentIonMasslist(peplist,Index_of_protein_sequence)
    
    peplist<-peplist[names(peplist) %in% names(pimlist) ==TRUE] 
    
    tempdf<-parLapply(cl=cl,  1: length(names(peplist)), Peptide_Summary_para,peplist)
    
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

Protein_feature_list_fun<-function(workdir=WorkingDir(),
                                   Fastadatabase,
                                   format="fasta",
                                   nrec=-1L, 
                                   skip=0L, 
                                   seek.first.rec=FALSE,
                                   use.names=TRUE, 
                                   with.qualities=FALSE,
                                   Digestion_site="[G]",
                                   missedCleavages=0:1,
                                   adducts=c("M+H","M+NH4","M+Na"),
                                   cl){
  library(Biostrings)
  list_of_protein_sequence<-readAAStringSet(Fastadatabase,
                                            format,
                                            nrec, 
                                            skip, 
                                            seek.first.rec,
                                            use.names, 
                                            with.qualities
  )    
  
  Index_of_protein_sequence<-fasta.index(Fastadatabase,
                                         nrec, 
                                         skip)   
  
  peplist<-cleave(as.character(list_of_protein_sequence),custom=Digestion_site, missedCleavages=missedCleavages)
  
  
  peplist<-peplist[duplicated(names(peplist))==FALSE]
  
  
  
  pimlist<-parentIonMasslist(peplist,Index_of_protein_sequence)
  
  peplist<-peplist[names(peplist) %in% names(pimlist) ==TRUE] 
  
  tempdf<-parLapply(cl=cl,  1: length(names(peplist)), Peptide_Summary_para,peplist)
  
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