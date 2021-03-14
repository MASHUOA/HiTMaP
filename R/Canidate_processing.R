
#' Meta_feature_list_fun
#'
#' This is a function that prepare the candiate list for maldi imaging data qualitative or quantitative analysis.
#' this function will read the candidate list file and generate mz for the adducts list defined in \code{"adducts"}. 
#' @param database the file name of candidate list
#' @param workdir the folder that contains candidate list
#' @param adducts  the adducts list to be used for generating the PMF search candidates
#' @param cal.mz If set with \code{"TRUE"}, the function will recalculate the mz value according to the column named "formula" in the \code{database} and the specified adducts.
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
                                adducts=c("M-H","M+Cl"),
                                cal.mz=TRUE,
                                bypass=FALSE,
                                BPPARAM=bpparam()){
  suppressMessages(suppressWarnings(require("Rcpp")))
  suppressMessages(suppressWarnings(require(dplyr)))
  suppressMessages(suppressWarnings(require(Rdisop)))
  suppressMessages(suppressWarnings(require(Biostrings)))
  suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  suppressMessages(suppressWarnings(require(BiocParallel)))
  adductslist<-Build_adduct_list()
  candidates<-read.csv(paste0(workdir,"/",database),as.is = TRUE)
  
  
  if (bypass){
    required_col=c("mz","adduct","moleculeNames")
    candidates=data_test_rename(required_col,candidates)
    Meta_Summary=candidates[is.na(candidates$mz)!=T,]
    Meta_Summary$Metabolite.Name=Meta_Summary$moleculeNames
    
  }else{
    Meta_Summary=candidates
    required_col=c("moleculeNames")
    candidates=data_test_rename(required_col,candidates)
    Meta_Summary$Metabolite.Name=Meta_Summary$moleculeNames
    candidates$mass<-NULL
    if (cal.mz==F){
      required_col=c("mz")        
      candidates$adduct=candidates$adducts
      candidates$adducts<-NULL
      candidates=data_test_rename(required_col,candidates)
      if (sum(candidates$adduct!=0)>0){
        
        adduct_df<-do.call(rbind,(lapply(unique(candidates$adduct),function(x,adductslist){
          adducts_mass<-adductslist[adductslist$Name==x,"Mass"]
          adducts_Charge<-adductslist[adductslist$Name==x,"Charge"]
          return(data.frame(adduct=x,adducts_mass=adducts_mass,charge=adducts_Charge))
        },adductslist)))
        
        candidates$adducts_mass<-NULL
        candidates<-merge(candidates,adduct_df,by="adduct",all.x=T)
        candidates$mass<-candidates$mz*abs(as.numeric(candidates$charge))-0.0005485799*as.numeric(candidates$charge)-(candidates$adducts_mass*abs(as.numeric(candidates$charge)))
      }else{
        candidates$mass<-candidates$mz
      }
    }else{
      
      required_col=c("formula")
      candidates=data_test_rename(required_col,candidates)
      unique_formula<-as.character(unique(candidates$formula))
      unique_formula_list<-lapply(unique_formula,get_atoms)
      masslist<-lapply(unique_formula_list,MonoisotopicMass)
      uniquemass<-unlist(masslist)
      mass_DF<-data.frame(formula=unique_formula,mass=uniquemass,stringsAsFactors = F)
      candidates<-base::merge(candidates,mass_DF,by="formula")
      
      
      #candidates_mass<-candidates$formula %>% lapply(getMonomass) 
      #candidates$mass<-as.numeric(unlist(candidates_mass))
    }
    
    candidates<-candidates[duplicated(names(candidates))==FALSE]
    if (is.null(candidates$Metabolite.Name)){candidates$Metabolite.Name<-candidates$moleculeNames}
    Meta_Summary<-NULL
    required_col=c("mz","Metabolite.Name","mass","formula")
    candidates=data_test_rename(required_col,candidates) 
    candidates$formula<-as.character(candidates$formula)
    candidates<-candidates[candidates$formula!="",]
    candidates<-candidates[!grepl(")n",candidates$formula),]
    meta_symbol<-lapply(unique(candidates$formula),get_atoms)
    symbol_adducts=bplapply(adducts,convert_peptide_adduct_list,meta_symbol,BPPARAM = BPPARAM,adductslist=adductslist,ConvertPeptide=ConvertPeptide)

    
    symbol_adducts_df=lapply(symbol_adducts,
                             function(x,names_x){
                               data.frame(formula=names_x,formula_adduct=x)
                             },unique(candidates$formula))

    
    for (i in 1:length(adducts)){
      candidates_adduct<-unique(candidates[,c("Metabolite.Name","mass","formula")])
      candidates_adduct<-candidates_adduct[!duplicated(candidates_adduct[,c("Metabolite.Name","formula")]),]
      adductmass<-as.numeric(as.character(adductslist[adductslist$Name==adducts[i],"Mass"]))
      adducts_Charge<-as.numeric(adductslist[adductslist$Name==adducts[i],"Charge"])
      candidates_adduct$mz<-(as.numeric(candidates_adduct$mass+adductmass))/abs(adducts_Charge)
      candidates_adduct$adduct<-adducts[i]
      candidates_adduct$charge<-adducts_Charge
      candidates_adduct<-merge(candidates_adduct,symbol_adducts_df[i])
      candidates_adduct$formula<-candidates_adduct$formula_adduct
      candidates_adduct$formula_adduct<-NULL
      Meta_Summary<-rbind.data.frame(Meta_Summary,unique(candidates_adduct))
      
      
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
                                   Multiple_mode=c("sequential","parallel"),
                                   adducts=c("M+H","M+NH4","M+Na"),
                                   BPPARAM=bpparam(),
                                   Decoy_search=T,
                                   Decoy_mode=c("adducts","elements","isotope","sequence"),
                                   Decoy_adducts=c("M+He","M+Ne","M+Ar","M+Kr","M+Xe","M+Rn"),
                                   Substitute_AA=list(AA=c(NULL),AA_new_formula=c(NULL),Formula_with_water=c(NULL)),
                                   mzrange=c(500,4000),
                                   output_candidatelist=T,
                                   Modifications=list(fixed=NULL,fixmod_position=NULL,variable=NULL,varmod_position=NULL),
                                   use_previous_candidates=F,
                                   Protein_desc_of_exclusion=NULL,
                                   Database_stats=F
                                   ){
  if (mzrange=="auto-detect") mzrange=c(500,4000)
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
   
   parse_cleavage_rule<-function(Digestion_site){
    Cleavage_rules<-Cleavage_rules_fun()
    found_enzyme<-Digestion_site[Digestion_site %in% names(Cleavage_rules)]
    not_found_enzyme<-Digestion_site[!(Digestion_site %in% names(Cleavage_rules))]
    found_rule<-not_found_enzyme[not_found_enzyme %in% (Cleavage_rules)]
    not_found_rule<-not_found_enzyme[!(not_found_enzyme %in% (Cleavage_rules))]
    Digestion_site_rule<-Cleavage_rules[found_enzyme] 
    Digestion_site_final<- unique(c(Digestion_site_rule,found_rule,not_found_rule))
     message("Found enzyme: ",paste(found_enzyme,collapse  = " "))
     message("Found rule: \"",paste(found_rule,collapse = " "),"\"")
     message("Found customized rule: \"",not_found_rule,"\"")
    Digestion_site_final
   }
   Digestion_site<-Digestion_site[Digestion_site!=""]
   Digestion_site<-parse_cleavage_rule(Digestion_site)
   
   if ('&'(length(Digestion_site)>=2,Multiple_mode[1]=="sequential")){
     Digestion_site<-paste0(Digestion_site,collapse = "|")
   }
  
  missedCleavages<<-missedCleavages
  Decoy_adducts=Decoy_adducts[!(Decoy_adducts %in% adducts)]
  Decoy_adducts=Decoy_adducts[1:length(adducts)]
     
 
  #if (length(list_of_protein_sequence)<2000){bpworkers(BPPARAM)=3}
  
  #if (length(list_of_protein_sequence)<500){bpworkers(BPPARAM)=1}


  
  Index_of_protein_sequence <<- fasta.index(database,
                                         nrec=-1L, 
                                         skip=0L)  
  
  
  list_of_protein_sequence <- readAAStringSet(database,
                                            format="fasta",
                                            nrec=-1L, 
                                            skip=0L, 
                                            seek.first.rec=FALSE
                                            ) 
  
  if (Decoy_search && ("sequence" %in% Decoy_mode)){
    list_of_protein_sequence_rev<-Biostrings::reverse(list_of_protein_sequence)
    names(list_of_protein_sequence_rev)<-paste0("Decoy_",names(list_of_protein_sequence_rev))
    list_of_protein_sequence<-c(list_of_protein_sequence,list_of_protein_sequence_rev)
    Index_of_protein_sequence_rev<-Index_of_protein_sequence
    Index_of_protein_sequence_rev$recno<-Index_of_protein_sequence_rev$recno+nrow(Index_of_protein_sequence_rev)
    Index_of_protein_sequence_rev$desc<-paste0("Decoy_", Index_of_protein_sequence_rev$desc)
    Index_of_protein_sequence<-rbind(Index_of_protein_sequence,Index_of_protein_sequence_rev)
    }
  
  #assign("list_of_protein_sequence", list_of_protein_sequence, envir=.GlobalEnv) 
  
  
  names_pro<-merge(data.frame(desc=names(list_of_protein_sequence),stringsAsFactors = F),Index_of_protein_sequence,by="desc",sort=F)
  
  names(list_of_protein_sequence) <- names_pro$recno
  Index_of_protein_sequence<<-Index_of_protein_sequence
  list_of_protein_sequence<<-list_of_protein_sequence
  
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
  
  Index_of_protein_sequence$Degestion=""
  
  peplist<-list()
  peplist_range<-list()
  Index_of_protein_sequence_list<-data.frame() 
  
  parentIonMasslist<-function(peplist,Index_of_protein_sequence){
    AA<-rep(0,26)
    PIM<-NULL
    for (i in 1:length(peplist)){
      PIM[[Index_of_protein_sequence$recno[i]]] <- parentIonMass(peplist[[Index_of_protein_sequence$recno[i]]],fixmod=AA)}
    return(PIM)
  }
  
  for (i in 1:length(Digestion_site)){
    message(paste("Testing fasta sequances for degestion site:",Digestion_site[i]))
    if (Digestion_site[i]==""){Digestion_site[i]="J"}
    peplist_range_option<-cleavageRanges(as.character(list_of_protein_sequence),custom=Digestion_site[i], missedCleavages=missedCleavages)
    peplist_option<-cleave(as.character(list_of_protein_sequence),custom=Digestion_site[i], missedCleavages=missedCleavages,unique =FALSE)
    Index_of_protein_sequence_option<-Index_of_protein_sequence
    pimlist<-lapply(peplist_option,function(x){
    AA<-rep(0,26)
    protViz::parentIonMass(x,fixmod=AA)})
    peplist_option<-peplist_option[duplicated(names(peplist_option))==FALSE]
    peplist_option<-peplist_option[names(peplist_option) %in% names(pimlist) ==TRUE]
    
    peplist_range_option<-peplist_range_option[duplicated(names(peplist_range_option))==FALSE]
    peplist_range_option<-peplist_range_option[names(peplist_range_option) %in% names(pimlist) ==TRUE]
    
    Index_of_protein_sequence_option<-Index_of_protein_sequence_option[Index_of_protein_sequence_option$recno %in% names(pimlist),]
    
    Index_of_protein_sequence_option$Degestion=Digestion_site[i]
    
    peplist_range<-c(peplist_range,peplist_range_option)
    
    peplist<-c(peplist,peplist_option)
    
    Index_of_protein_sequence_list<-rbind(Index_of_protein_sequence_list,Index_of_protein_sequence_option)
    
  }
  
  #message(paste("Peptide list generated",length(peplist),"entries in total."))
  
  
  
  
  
  #pimlist<-parentIonMasslist(peplist,Index_of_protein_sequence_list)
  
  
  AA<-c(71.037114, 0.000000, 103.009185, 115.026943, 129.042593, 147.068414, 
        57.021464, 137.058912, 113.084064, 0.000000, 128.094963, 113.084064, 
        131.040485, 114.042927, 0.000000, 97.052764, 128.058578, 156.101111, 
        87.032028, 101.047679, 0.000000, 99.068414, 186.079313, 0.000000, 
        163.063329, 100.994269)
  
  names(AA)<-LETTERS
  

  #tempdf<-parLapply(cl=cl,  1: length(names(peplist)), Peptide_Summary_para,peplist)
  #bplapply()
  message(paste("Generated",length(peplist),"Proteins in total. Computing exact masses..."))
  #tempdf<-bplapply( 1: length((peplist)), Peptide_Summary_para,peplist,BPPARAM = BPPARAM)
  #tempdf<-lapply( 1: length((peplist)), Peptide_Summary_para,peplist)
  #tempdf <- do.call("rbind", tempdf)
  #colnames(tempdf)<-c("Protein","Peptide")
  #do.call(rbind,peplist_range)
  start_end<-do.call(rbind,peplist_range)
  
  
  Protein.df=rep(names(peplist),unname(unlist(lapply(peplist,length))))
  Peptide.df=unname(unlist(peplist))
  tempdf<-data.frame(Protein=Protein.df,Peptide=Peptide.df,start=start_end[,1],end=start_end[,2],stringsAsFactors = FALSE)
  
  #do.call()
  
  #
  #list_of_protein_sequence<-get("list_of_protein_sequence", envir = .GlobalEnv)
  
  
  tempdf$Protein<-as.character(tempdf$Protein)
  tempdf$Peptide<-as.character(tempdf$Peptide)
  tempdf<-as.data.frame(tempdf)
  tempdf$Modification=""
  pro_end<-sapply(list_of_protein_sequence,length)
  pro_end<-data.frame(Protein=names(pro_end),pro_end=pro_end)
  tempdf<-merge(tempdf,pro_end,by="Protein")
  if (length(grep("X",tempdf$Peptide))!=0 && !("X" %in% Substitute_AA$AA))  tempdf<-tempdf[-grep("X",tempdf$Peptide),]
  if (length(grep("U",tempdf$Peptide))!=0 && !("U" %in% Substitute_AA$AA))  tempdf<-tempdf[-grep("U",tempdf$Peptide),]
  
  mod.df_fix<-Peptide_modification(retrive_ID = Modifications$fixed,mod_position=Modifications$fixmod_position)
  mod.df_var<-Peptide_modification(retrive_ID = Modifications$variable,mod_position=Modifications$varmod_position)
  mod.df<-rbind(mod.df_fix,mod.df_var)
  if (is.null(mod.df)){
    min_mod_massdiff<-100
    max_mod_massdiff<-500
  }else{
  min_mod_massdiff<-min(as.numeric(mod.df$mono_mass))
  max_mod_massdiff<-max(as.numeric(mod.df$mono_mass))
  
  if(min_mod_massdiff>0){min_mod_massdiff=0}
  if(max_mod_massdiff<0){max_mod_massdiff=0}
  }
  
  if (!is.null(Substitute_AA$AA)) {
    Substitute_AA_df<-data.frame(Substitute_AA[c("AA","AA_new_formula","Formula_with_water")],stringsAsFactors = F)
    for (AA_row in 1:nrow(Substitute_AA_df)){
      if(is.null(Substitute_AA_df$Formula_with_water[AA_row])){
        AA[Substitute_AA_df$AA[AA_row]]=rcdk::get.formula((Substitute_AA_df$AA_new_formula[AA_row]),charge = 0)@mass
      }else{
       if(Substitute_AA_df$Formula_with_water[AA_row]){
         AA[Substitute_AA_df$AA[AA_row]]=rcdk::get.formula(get_formula_from_atom(merge_atoms(get_atoms(Substitute_AA_df$AA_new_formula[AA_row]),get_atoms("H2O"),mode = "ded")),charge = charge)@mass
       }else { AA[Substitute_AA_df$AA[AA_row]]=rcdk::get.formula((Substitute_AA_df$AA_new_formula[AA_row]),charge = 0)@mass} 
      }
     Substitute_AA$aavector[[AA_row]]<-unlist(get_atoms(Substitute_AA$AA_new_formula[AA_row])) 
    }
    
  }
  
  #parentIonMass gives the M+H
  tempdf$pepmz <- as.numeric(parentIonMass(tempdf$Peptide,fixmod=AA)- 1.007276 )
  tempdf<-tempdf['&'(tempdf$pepmz>=mzrange[1]-max_mod_massdiff,tempdf$pepmz<=mzrange[2]+min_mod_massdiff),]
  #tempdf1<-tempdf
  Protein_Summary<-NULL
  adductslist<-Build_adduct_list()
  tempdf$Modification<-""
  
  
  message(paste("Generating peptide formula..."))
  uniquepep=(tempdf$Peptide)
  uniqueAA=unique(strsplit(paste0(uniquepep,collapse = ""), "")[[1]])
  
  Element_tbl<-BuildElement(Substitute_AA=Substitute_AA,uniqueAA=uniqueAA)
  peptide_symbol=lapply(uniquepep,ConvertPeptide,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl)
  # peptide_symbol=bplapply(uniquepep,ConvertPeptide,BPPARAM = BPPARAM,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl)
  # 
  # Element_tbl_m<-BuildElement(Substitute_AA=Substitute_AA,uniqueAA=uniqueAA,element_matrix = T)
  # peptide_symbol=lapply(uniquepep,ConvertPeptide.2,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl_m)
  # peptide_symbol=bplapply(uniquepep,ConvertPeptide.2,BPPARAM = BPPARAM,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl_m)
  # 
  # 
  # system.time({peptide_symbol1=lapply(uniquepep[1:200],ConvertPeptide,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl)})
  # system.time({peptide_symbol2=lapply(uniquepep[1:200],ConvertPeptide.2,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl_m)})
  # 
  # system.time({peptide_symbol3=bplapply(uniquepep[1:200],ConvertPeptide,BPPARAM = SerialParam(),Substitute_AA=Substitute_AA,Element_tbl=Element_tbl)})
  # system.time({peptide_symbol4=bplapply(uniquepep[1:200],ConvertPeptide.2,BPPARAM = SerialParam(),Substitute_AA=Substitute_AA,Element_tbl=Element_tbl_m)})
  # system.time({peptide_symbol5=lapply(uniquepep,ConvertPeptide,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl)})
  # system.time({peptide_symbol6=lapply(uniquepep,ConvertPeptide.2,Substitute_AA=Substitute_AA,Element_tbl=Element_tbl_m)})
  
  if (!is.null(Modifications$fixed)){
    mod.df<-Peptide_modification(retrive_ID = Modifications$fixed,mod_position=Modifications$fixmod_position)
    if(length(unique(mod.df$full_name))==length(Modifications$fixed)){
    message(paste("Fixed modifications:",unique(mod.df$full_name),"found in unimod DB",sep=" ",collapse = ", "))
    #peptide_symbol=bplapply(mod.df,convert_peptide_fixmod,peptide_symbol,BPPARAM = BPPARAM,pep_sequence=tempdf$Peptide,ConvertPeptide=ConvertPeptide)

    
    peptide_symbol_var=convert_peptide_fixmod(mod.df,peptide_symbol,peptide_info=tempdf,BPPARAM = BPPARAM)
    message(paste("Merge modification formula done."))
    tempdf_var<-tempdf
    reserve_entry<-rep(FALSE,nrow(tempdf_var))
    for (fixmod in mod.df$record_id){
      mods<-ifelse(peptide_symbol_var$multiplier[[fixmod]]>=1,mod.df$code_name[mod.df$record_id==fixmod],"")
      tempdf_var$Modification<-paste(tempdf_var$Modification,mods)
      tempdf_var$pepmz <-tempdf_var$pepmz + peptide_symbol_var$multiplier[[fixmod]]*as.numeric(mod.df$mono_mass[mod.df$record_id==fixmod])
      reserve_entry<-ifelse('&'(peptide_symbol_var$multiplier[[fixmod]]>=1,reserve_entry==FALSE),TRUE,reserve_entry)
    }
    
    tempdf<-tempdf_var
    peptide_symbol<-peptide_symbol_var$peptide_symbol
    }else{
    message(paste("warning:",
                  Modifications$fixed['&'((Modifications$fixed %in% unique(mod.df$code_name))==F , (Modifications$fixed %in% unique(mod.df$record_id) )==F)],
                  "not found in unimod DB. Please check the unimod DB again. Any Charactor string that matches code name or number that matches mod ID will be OKAY.",sep=" ",collapse = ", "))
    }
    
  }
  
  if (!is.null(Modifications$variable)){
    mod.df<-Peptide_modification(retrive_ID = Modifications$variable,mod_position=Modifications$varmod_position)
    if(length(unique(mod.df$full_name))==length(Modifications$variable)){
      message(paste("variable modifications:",unique(mod.df$full_name),"found in unimod DB",sep=" ",collapse = ", "))
      #peptide_symbol=bplapply(mod.df,convert_peptide_fixmod,peptide_symbol,BPPARAM = BPPARAM,pep_sequence=tempdf$Peptide,ConvertPeptide=ConvertPeptide)
      peptide_symbol_var<-convert_peptide_fixmod(mod.df,peptide_symbol,peptide_info=tempdf,BPPARAM = BPPARAM)
      message(paste("Merge modification formula done."))
      tempdf_var<-tempdf
      reserve_entry<-rep(FALSE,nrow(tempdf_var))
      for (fixmod in mod.df$record_id){
        mods<-ifelse(peptide_symbol_var$multiplier[[fixmod]]>=1,mod.df$code_name[mod.df$record_id==fixmod],"")
        tempdf_var$Modification<-paste(tempdf_var$Modification,mods)
        tempdf_var$pepmz <-tempdf_var$pepmz + peptide_symbol_var$multiplier[[fixmod]]*as.numeric(mod.df$mono_mass[mod.df$record_id==fixmod])
        reserve_entry<-ifelse('&'(peptide_symbol_var$multiplier[[fixmod]]>=1,reserve_entry==FALSE),TRUE,reserve_entry)
        }
      
      peptide_symbol_var<-peptide_symbol_var$peptide_symbol
      s1<-peptide_symbol[reserve_entry]
      s2<-peptide_symbol_var[reserve_entry]
      identical(peptide_symbol_var,peptide_symbol)
      peptide_symbol<-c(peptide_symbol,peptide_symbol_var[reserve_entry])
      tempdf<-rbind(tempdf,tempdf_var[reserve_entry,])
    }else{
      message(paste("warning:",
                    Modifications$variable['&'((Modifications$variable %in% unique(mod.df$code_name))==F , (Modifications$variable %in% unique(mod.df$record_id) )==F)],
                    "not found in unimod DB. Please check the unimod DB again. Any Charactor string that matches code name or number that matches mod ID will be OKAY.",sep=" ",collapse = ", "))
    }
    
  }
  
  message(paste("Generating peptide formula with adducts:",paste(adducts,collapse = " ")))
  peptides_symbol_adducts=bplapply(adducts,convert_peptide_adduct_list,peptide_symbol,BPPARAM = BPPARAM,adductslist=adductslist)
  
  for (i in 1:length(adducts)){
    adductmass <- as.numeric(as.character(adductslist[adductslist$Name == adducts[i], "Mass"]))
    charge=as.numeric(as.character(adductslist$Charge[adductslist$Name == adducts[i]]))
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
     suppressMessages(suppressWarnings(require(enviPat)))
     suppressMessages(suppressWarnings(require(rcdk)))
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
  
  Protein_Summary$Modification[is.na(Protein_Summary$Modification)]<-""
  Protein_Summary$mz<-round(Protein_Summary$mz,digits = 4)
  Protein_Summary<-Protein_Summary[`&`(Protein_Summary$mz>=mzrange[1],Protein_Summary$mz<=mzrange[2]),]
  #Protein_Summary$Protein=Index_of_protein_sequence[Index_of_protein_sequence$desc==Protein_Summary$Protein]
  #temp_index=Index_of_protein_sequence
  #temp_index$Protein=temp_index$desc
  #temp_index=temp_index[,c("Protein","recno")]
  #Protein_Summary=merge(Protein_Summary,temp_index,by="Protein",all.x=T)
  #Protein_Summary$Protein=Protein_Summary$recno
  if (output_candidatelist){
    if (dir.exists(paste(workdir,"/Summary folder",sep=""))==FALSE){dir.create(paste(workdir,"/Summary folder",sep=""))}
    write.csv(Protein_Summary,paste(workdir,"/Summary folder/candidatelist.csv",sep=""),row.names = F)
    write.csv(Index_of_protein_sequence,paste(workdir,"/Summary folder/protein_index.csv",sep=""),row.names = F)
    message("Candidate list has been exported.")
  }
  
  
  
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

      #mz_vs_peptide_filtered_order <- mz_vs_peptide_filtered %>%  rlist::list.subset(paste(ppm_test_list,"ppm"))
      
      mz_vs_peptide_filtered<-lapply(mz_vs_peptide_filtered,`colnames<-`,c("mz","unique_formula"))
      

      mycol=RColorBrewer::brewer.pal(name = "Set1",n=length(names(mz_vs_peptide_filtered)))
      names(mycol)=names(mz_vs_peptide_filtered)
      mycol_name<-mycol
      names(mycol_name)<-mycol
      sp<-ggplot2::ggplot() 
      sp <-sp + geom_segment(data=mz_vs_peptide_filtered[[1]],mapping = aes(x=mz, y=unique_formula,xend=mz,yend=rep(0,length(unique_formula)),col=names(mycol)[1]),size=0.1)
      sp <-sp + geom_segment(data=mz_vs_peptide_filtered[[2]],mapping = aes(x=mz, y=unique_formula,xend=mz,yend=rep(0,length(unique_formula)),col=names(mycol)[2]),size=0.1)
      sp <-sp + geom_segment(data=mz_vs_peptide_filtered[[3]],mapping = aes(x=mz, y=unique_formula,xend=mz,yend=rep(0,length(unique_formula)),col=names(mycol)[3]),size=0.1)
      
      # for(i in 1:length(mz_vs_peptide_filtered)){
      #   sp <-sp + geom_segment(data=mz_vs_peptide_filtered[[i]][1:10000,],mapping = aes(x=mz, y=unique_formula,xend=mz,yend=rep(0,length(unique_formula)),color=mycol[i]),color=mycol[i],size=0.1)
      # }

      #sp <-sp + scale_color_manual(name="m/z bin size") 
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
    
    dif_kmeans=kmeans(1/(mz_vs_resolution$mz_unqiue_diff),centers = 200,iter.max = 500)
    
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
  
  if(F){if (Decoy_search && ("isotope" %in% Decoy_mode)){
    message("attaching decoy IDs in isotope mode...")
    if(is.null(Protein_Summary$isdecoy)){
      Protein_Summary$isdecoy=0
    }
    
    Protein_feature_list_decoy<-Protein_Summary[Protein_Summary$isdecoy==0,]
    #Protein_feature_list_decoy
    Protein_feature_list_decoy$isdecoy=1
    Protein_Summary<-rbind(Protein_Summary[Protein_Summary$isdecoy==0,],Protein_feature_list_decoy)
    Protein_feature_list<<-Protein_Summary
    message("attaching decoy IDs in isotope mode...Done")
  }}
  
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
  
  Protein_feature_list<<-Protein_Summary 
  

  
  return(Protein_Summary)
}


vendiagram<-function(){
   suppressMessages(suppressWarnings(require(stringr)))
   suppressMessages(suppressWarnings(require(tcltk)))
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
  
   suppressMessages(suppressWarnings(require(Biostrings)))
   suppressMessages(suppressWarnings(require(cleaver)))
   suppressMessages(suppressWarnings(require(protViz)))
   suppressMessages(suppressWarnings(require(rcdk)))
  
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

convert_peptide_adduct<-function(peptide_formula,adductsname,multiplier=c(1,1),adductslist=Build_adduct_list()){
   suppressMessages(suppressWarnings(require(rcdk)))
   suppressMessages(suppressWarnings(require(OrgMassSpecR)))
  
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
  
  
  adductsformula_add= as.character(adductslist[adductslist$Name==adductsname,"formula_add"])
  adductsformula_ded= as.character(adductslist[adductslist$Name==adductsname,"formula_ded"])
  
  null_if_false<-function(x){if (x==FALSE){""} else{get.formula(x)@string}}
  
  adductsformula_add<-null_if_false(adductsformula_add)
  adductsformula_ded<-null_if_false(adductsformula_ded)
  
  adductsformula_add<-get_atoms(adductsformula_add)
  adductsformula_ded<-get_atoms(adductsformula_ded)
  
  adductsformula<-merge_atoms(adductsformula_add,adductsformula_ded,check_merge = F,mode = "ded")
  
  formula_with_adducts<-merge_atoms(peptide_formula,adductsformula,check_merge = T,mode = "add", multiplier = multiplier)
  
  for (name in names(formula_with_adducts)){
    if(formula_with_adducts[[name]]==0){ formula_with_adducts[[name]]=NULL}
  }

  return(paste0(names(formula_with_adducts),formula_with_adducts,collapse = ""))
  
}


convert_peptide_fixmod<-function(mod.df,peptide_symbol,peptide_info,BPPARAM=BPPARAM){
  
  
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
  
  #if (missing(multiplier)){
  #  multiplier[1]= as.numeric(as.character(adductslist[adductslist$Name==adductsname,"Mult"]))
  #  multiplier[2]=1
  #}
  
  mod.df.list <- split(mod.df, seq(nrow(mod.df)))
  multiplier_for_mod<-function(x,pep_sequence,peptide_info,BPPARAM=bpparam()){
    if(x$position_key==2){
      #multiplier_pep<-lapply(pep_sequence,grepl,x$one_letter)
      multiplier_pep<-str_locate_all(pep_sequence,x$one_letter)
      multiplier_pep<-sapply(multiplier_pep,nrow)
      return(multiplier_pep)
    } else if(x$position_key %in% c(3,4)){
      multiplier_pep=rep(1,length(pep_sequence))
      return(multiplier_pep)
    } else if(x$position_key %in% c(5)){
      message("Protein N-term modification selected")
      #list_of_protein_sequence<-get("list_of_protein_sequence", envir = .GlobalEnv)
      
      #multiplier_pep=rep(0,length(pep_sequence))
      #map_res<-sapply(list_of_protein_sequence,str_locate_all,pep_sequence)
      #map_res_found<-bplapply(1:length(pep_sequence),function(x,map_res){
      #  map_resdf<-(do.call(rbind,map_res[x,]))
      #  if (1 %in% map_resdf[,"start"]){
      #    if (sum(map_resdf[,"start"]!=1)>0){
      #      return(c(1,1))
      #    }else{return(c(1,0))}
      #  }else{return(c(0,0))}
      #},map_res,BPPARAM=BPPARAM)
      #map_res_found<-do.call(rbind,map_res_found)
      #multiplier_pep<-map_res_found[,1]
      multiplier_pep<-ifelse(peptide_info[,"start"]==1,1,0)
      return(multiplier_pep)
    } else if(x$position_key %in% c(6)){
      message("Protein C-term modification selected")
      multiplier_pep=rep(0,length(pep_sequence))
      #map_res<-sapply(list_of_protein_sequence,str_locate_all,pep_sequence)
      #map_res<-bplapply(list_of_protein_sequence,str_locate_all,pep_sequence)
      
      multiplier_pep<-ifelse(peptide_info$end==peptide_info$pro_end,1,0)
      #map_res_found<-bplapply(1:length(pep_sequence),function(x,map_res,pro_end){
      #  map_resdf<-do.call(rbind,map_res[x,])
      #  pro_label<-rep(pro_end,sapply(map_res[x,],nrow))
      #  map_resdf[,"end"]==pro_label
      #  if (sum(map_resdf[,"end"]==pro_label)>=1){
      #    if (sum(map_resdf[,"end"]!=pro_label)>0){
      #      return(c(1,1))
      #    }else{return(c(1,0))}
      #  }else{return(c(0,0))}
      #},map_res,pro_end,BPPARAM = BPPARAM)
      #map_res_found<-do.call(rbind,map_res_found)
      #multiplier_pep<-map_res_found[,1]
      return(multiplier_pep)
    }
    
  }
  
  formula_mod<-lapply(mod.df$composition,get_atoms_mod)
  
  names(formula_mod)<-mod.df$record_id
  
  multiplier<-lapply(mod.df.list,multiplier_for_mod,pep_sequence=pep_sequence,peptide_info=peptide_info,BPPARAM=BPPARAM)
  
  names(multiplier)<-as.character(mod.df$record_id)
  #formula<-ConvertPeptide(peptide)
  for (fixmod in mod.df$record_id){
    #for (formula in 1:length(peptide_symbol)){
    #peptide_symbol[[formula]]<-merge_atoms(peptide_symbol[[formula]],formula_mod[[fixmod]],check_merge = F,mode = "add", multiplier = c(1,multiplier[[fixmod]][formula]))
    #}
    #message(formula_mod[[as.character(fixmod)]])
    #message(multiplier[[as.character(fixmod)]])
    peptide_symbol[which(multiplier[[as.character(fixmod)]]>=1)]
    peptide_symbol[which(multiplier[[as.character(fixmod)]]>=1)]<-bplapply(1:length(peptide_symbol[which(multiplier[[as.character(fixmod)]]>=1)]),function(x,symbol,addelements,merge_atoms,multiplier_list){
      if (multiplier_list[x]!=0){
        return(merge_atoms(atoms = symbol[[x]],addelements = addelements,check_merge=F,mode="add",multiplier=c(1,multiplier_list[x])))
      }else{
        return(symbol[[x]])
        }
      },symbol=peptide_symbol[which(multiplier[[as.character(fixmod)]]>=1)],merge_atoms=merge_atoms,addelements=formula_mod[[as.character(fixmod)]], multiplier_list = multiplier[[as.character(fixmod)]][which(multiplier[[as.character(fixmod)]]>=1)],BPPARAM = BPPARAM)
  }
  
  return(list(peptide_symbol=peptide_symbol,multiplier=multiplier))
  
}

convert_peptide_adduct_list<-function(adductsname,peptide_symbol,multiplier=c(1,1),adductslist=Build_adduct_list()){
   suppressMessages(suppressWarnings(require(rcdk)))
   suppressMessages(suppressWarnings(require(OrgMassSpecR)))

  merge_atoms<-function(atoms,addelements,check_merge=F,mode=c("add","ded"),multiplier=c(1,1)){
    
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
    
    if (check_merge==T){
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
     multiplier<-NULL
     multiplier[1]= as.numeric(as.character(adductslist[adductslist$Name==adductsname,"Mult"]))
     multiplier[2]=1
   }
  adductsformula_add= as.character(adductslist[adductslist$Name==adductsname,"formula_add"])
  adductsformula_ded= as.character(adductslist[adductslist$Name==adductsname,"formula_ded"])
  
  null_if_false<-function(x){if (x==FALSE){""} else{get.formula(x)@string}}
  
  adductsformula_add<-null_if_false(adductsformula_add)
  adductsformula_ded<-null_if_false(adductsformula_ded)
  
  adductsformula_add<-get_atoms(adductsformula_add)
  adductsformula_ded<-get_atoms(adductsformula_ded)
  
  adductsformula<-merge_atoms(adductsformula_add,adductsformula_ded,check_merge = F,mode = "ded")
  peptide_formula<-as.character()
  #formula<-ConvertPeptide(peptide)
  for (formula in 1:length(peptide_symbol)){
  #if(is.null(peptide_symbol[[formula]])) peptide_symbol[[formula]]=list(H=0)
   formula_with_adducts<-merge_atoms(peptide_symbol[[formula]],adductsformula,check_merge = F,mode = "add", multiplier = multiplier) 
   for (name in names(formula_with_adducts)){
    if(formula_with_adducts[[name]]==0){ formula_with_adducts[[name]]=NULL}
   }
   peptide_formula[[formula]]=(paste0(names(formula_with_adducts),formula_with_adducts,collapse = ""))
  }
  
  return(peptide_formula)
  
}
#' get_atoms
#'
#' This is a function that prepare the atoms list for candidate molecules.
#' 
#' @param Symbol the Symbol of the molecule
#' @return a atoms list
#'
#' @examples
#' get_atoms(Symbol="C7H13O8P")
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
#' This is a function that merge the atoms to a melecule atom list . 
#' @param atoms the file name of candidate list
#' @param addelements the folder that contains candidate list
#' @param check_merge  Check the final formula for any negative number of atom(s).
#' @param mode add the elements to the formula by using "add", and deduct the elements by using "ded".
#' @param multiplier  Two number to define the atom and adding elements' proportion respectively.
#' @return an atom list
#'
#' @examples
#' merge_atoms(atoms=get_atoms("C7H13O8P"),addelements=list(H=1),check_merge=T,mode="ded",multiplier=c(1,2))
#'
#' @export
merge_atoms<-function(atoms,addelements,check_merge=F,mode=c("add","ded"),multiplier=c(1,1)){
  
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
  
  if (check_merge==T){
    for (x in names(atoms)){
      if (atoms[[x]]<0){
      stop("add atoms failed due to incorrect number of",x)
      } 
    }
    
  }
  return(as.list(atoms))
  
}

get_formula_from_atom<-function(atoms){
  paste0(names(atoms),atoms,collapse = "")
}
  

ConvertPeptide_depr<-function (sequence, output = "elements",Substitute_AA=NULL) 
{
  peptideVector <- strsplit(sequence, split = "")[[1]]
  if (output == "elements") {
    FindElement <- function(residue,Substitute_AA) {
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
      
      if (!is.null(Substitute_AA)){
        if(residue %in% Substitute_AA$AA){
          element <- Substitute_AA$aavector[Substitute_AA$AA==residue][[1]]
        }
      }
      
      return(element)
    }
    aa_compo<-union(c("C","H","N","O","S"),unique(unlist(lapply(Substitute_AA$aavector,names))))
    resultsVector<-rep(0,length(aa_compo))
    names(resultsVector)=aa_compo
    resultsVector_bak<-resultsVector
    
    for (i in 1:length(peptideVector)) {
      to_add<-resultsVector_bak
      to_add<-FindElement(peptideVector[i],Substitute_AA = Substitute_AA) 
      to_add[aa_compo[(!(aa_compo %in% names(to_add)))]]=0
      to_add<-to_add[aa_compo]
      resultsVector <- to_add +  resultsVector
    }
    
    resultsVector["H"] <- resultsVector["H"] + 2
    resultsVector["O"] <- resultsVector["O"] + 1
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

BuildElement <- function(Substitute_AA=NULL,uniqueAA=c("A"),element_matrix=F) {
  if (is.null(Substitute_AA)){
    Element_tbl<-list(
      A=c(C = 3L, H = 5L, N = 1L, O = 1L, S = 0L),
      R=c(C = 6L, H = 12L, N = 4L, O = 1L, S = 0L),
      N=c(C = 4L, H = 6L, N = 2L, O = 2L, S = 0L),
      D=c(C = 4L, H = 5L, N = 1L, O = 3L, S = 0L),
      E=c(C = 5L, H = 7L, N = 1L, O = 3L, S = 0L),
      Q=c(C = 5L, H = 8L, N = 2L, O = 2L, S = 0L),
      G=c(C = 2L, H = 3L, N = 1L, O = 1L, S = 0L),
      H=c(C = 6L, H = 7L, N = 3L, O = 1L, S = 0L),
      I=c(C = 6L, H = 11L, N = 1L, O = 1L, S = 0L),
      L=c(C = 6L, H = 11L, N = 1L, O = 1L, S = 0L),
      K=c(C = 6L, H = 12L, N = 2L, O = 1L, S = 0L),
      M=c(C = 5L, H = 9L, N = 1L, O = 1L, S = 1L),
      F=c(C = 9L, H = 9L, N = 1L, O = 1L, S = 0L),
      P=c(C = 5L, H = 7L, N = 1L, O = 1L, S = 0L),
      S=c(C = 3L, H = 5L, N = 1L, O = 2L, S = 0L),
      T=c(C = 4L, H = 7L, N = 1L, O = 2L, S = 0L),
      W=c(C = 11L, H = 10L, N = 2L, O = 1L, S = 0L),
      Y=c(C = 9L, H = 9L, N = 1L, O = 2L, S = 0L),
      V=c(C = 5L, H = 9L, N = 1L, O = 1L, S = 0L),
      C=c(C = 3L, H = 5L, N = 1L, O = 1L, S = 1L),
      B=c(C = 4L, H = 5L, N = 1L, O = 3L, S = 0L),
      Z=c(C = 5L, H = 7L, N = 1L, O = 3L, S = 0L),
      J=c(C = 6L, H = 11L, N = 1L, O = 1L, S = 0L))
  }
  if (!is.null(Substitute_AA)){
    common_ele<-c("C","H","N","O","S")
    aa_compo<-union(common_ele,unique(unlist(lapply(Substitute_AA$aavector,names))))
    aa_sub<-aa_compo[!(aa_compo %in% common_ele)]
    Element_tbl<-list(
      A=c(C = 3L, H = 5L, N = 1L, O = 1L, S = 0L),
      R=c(C = 6L, H = 12L, N = 4L, O = 1L, S = 0L),
      N=c(C = 4L, H = 6L, N = 2L, O = 2L, S = 0L),
      D=c(C = 4L, H = 5L, N = 1L, O = 3L, S = 0L),
      E=c(C = 5L, H = 7L, N = 1L, O = 3L, S = 0L),
      Q=c(C = 5L, H = 8L, N = 2L, O = 2L, S = 0L),
      G=c(C = 2L, H = 3L, N = 1L, O = 1L, S = 0L),
      H=c(C = 6L, H = 7L, N = 3L, O = 1L, S = 0L),
      I=c(C = 6L, H = 11L, N = 1L, O = 1L, S = 0L),
      L=c(C = 6L, H = 11L, N = 1L, O = 1L, S = 0L),
      K=c(C = 6L, H = 12L, N = 2L, O = 1L, S = 0L),
      M=c(C = 5L, H = 9L, N = 1L, O = 1L, S = 1L),
      F=c(C = 9L, H = 9L, N = 1L, O = 1L, S = 0L),
      P=c(C = 5L, H = 7L, N = 1L, O = 1L, S = 0L),
      S=c(C = 3L, H = 5L, N = 1L, O = 2L, S = 0L),
      T=c(C = 4L, H = 7L, N = 1L, O = 2L, S = 0L),
      W=c(C = 11L, H = 10L, N = 2L, O = 1L, S = 0L),
      Y=c(C = 9L, H = 9L, N = 1L, O = 2L, S = 0L),
      V=c(C = 5L, H = 9L, N = 1L, O = 1L, S = 0L),
      C=c(C = 3L, H = 5L, N = 1L, O = 1L, S = 1L),
      B=c(C = 4L, H = 5L, N = 1L, O = 3L, S = 0L),
      Z=c(C = 5L, H = 7L, N = 1L, O = 3L, S = 0L),
      J=c(C = 6L, H = 11L, N = 1L, O = 1L, S = 0L))
    
    for (AA in Substitute_AA$AA){
      aasub_ele_list<-Substitute_AA$aavector[Substitute_AA$AA==AA][[1]]
      if(Substitute_AA$Formula_with_water[Substitute_AA$AA==AA][[1]]){
        aasub_ele_list[["H"]]<-aasub_ele_list[["H"]]-2
        aasub_ele_list[["O"]]<-aasub_ele_list[["O"]]-1
      }
      mode(aasub_ele_list)<- "integer"
      Element_tbl[[AA]]<-aasub_ele_list
    }
    
    for (AA in names(Element_tbl)){      
      to_add<-Element_tbl[[AA]]
      to_add[aa_compo[(!(aa_compo %in% names(to_add)))]]=0L
      to_add<-to_add[aa_compo]
      #print(to_add)
      Element_tbl[[AA]] <- to_add
      
    }  
  }
  
  if(length(uniqueAA[!(uniqueAA %in% names(Element_tbl))])>0){
    for (AA in uniqueAA[!(uniqueAA %in% names(Element_tbl))]){
      to_add<-Element_tbl[["A"]]
      to_add=rep(0L,length(to_add))
      names(to_add)<-names(Element_tbl[["A"]])
      #print(to_add)
      Element_tbl[[AA]] <- to_add
      message(paste(AA,"is found in database as an uncommon amino acid, use Substitute_AA to define a formula for",AA))
    }
  }
  
  if (element_matrix) Element_tbl<-do.call(rbind,Element_tbl)
  
  return(Element_tbl)
}

BuildElement_1 <- function(Substitute_AA=NULL) {
  if (is.null(Substitute_AA)){
    Element_tbl<-list(
      A=c(C = 3, H = 5, N = 1, O = 1, S = 0),
      R=c(C = 6, H = 12, N = 4, O = 1, S = 0),
      N=c(C = 4, H = 6, N = 2, O = 2, S = 0),
      D=c(C = 4, H = 5, N = 1, O = 3, S = 0),
      E=c(C = 5, H = 7, N = 1, O = 3, S = 0),
      Q=c(C = 5, H = 8, N = 2, O = 2, S = 0),
      G=c(C = 2, H = 3, N = 1, O = 1, S = 0),
      H=c(C = 6, H = 7, N = 3, O = 1, S = 0),
      I=c(C = 6, H = 11, N = 1, O = 1, S = 0),
      L=c(C = 6, H = 11, N = 1, O = 1, S = 0),
      K=c(C = 6, H = 12, N = 2, O = 1, S = 0),
      M=c(C = 5, H = 9, N = 1, O = 1, S = 1),
      F=c(C = 9, H = 9, N = 1, O = 1, S = 0),
      P=c(C = 5, H = 7, N = 1, O = 1, S = 0),
      S=c(C = 3, H = 5, N = 1, O = 2, S = 0),
      T=c(C = 4, H = 7, N = 1, O = 2, S = 0),
      W=c(C = 11, H = 10, N = 2, O = 1, S = 0),
      Y=c(C = 9, H = 9, N = 1, O = 2, S = 0),
      V=c(C = 5, H = 9, N = 1, O = 1, S = 0),
      C=c(C = 3, H = 5, N = 1, O = 1, S = 1),
      Z=c(C = 3, H = 3, N = 1, O = 1, S = 1))
  }
  if (!is.null(Substitute_AA)){
    common_ele<-c("C","H","N","O","S")
    aa_compo<-union(common_ele,unique(unlist(lapply(Substitute_AA$aavector,names))))
    aa_sub<-aa_compo[!(aa_compo %in% common_ele)]
    Element_tbl<-list(
      A=c(C = 3, H = 5, N = 1, O = 1, S = 0),
      R=c(C = 6, H = 12, N = 4, O = 1, S = 0),
      N=c(C = 4, H = 6, N = 2, O = 2, S = 0),
      D=c(C = 4, H = 5, N = 1, O = 3, S = 0),
      E=c(C = 5, H = 7, N = 1, O = 3, S = 0),
      Q=c(C = 5, H = 8, N = 2, O = 2, S = 0),
      G=c(C = 2, H = 3, N = 1, O = 1, S = 0),
      H=c(C = 6, H = 7, N = 3, O = 1, S = 0),
      I=c(C = 6, H = 11, N = 1, O = 1, S = 0),
      L=c(C = 6, H = 11, N = 1, O = 1, S = 0),
      K=c(C = 6, H = 12, N = 2, O = 1, S = 0),
      M=c(C = 5, H = 9, N = 1, O = 1, S = 1),
      F=c(C = 9, H = 9, N = 1, O = 1, S = 0),
      P=c(C = 5, H = 7, N = 1, O = 1, S = 0),
      S=c(C = 3, H = 5, N = 1, O = 2, S = 0),
      T=c(C = 4, H = 7, N = 1, O = 2, S = 0),
      W=c(C = 11, H = 10, N = 2, O = 1, S = 0),
      Y=c(C = 9, H = 9, N = 1, O = 2, S = 0),
      V=c(C = 5, H = 9, N = 1, O = 1, S = 0),
      C=c(C = 3, H = 5, N = 1, O = 1, S = 1),
      Z=c(C = 3, H = 3, N = 1, O = 1, S = 1)) 
    
    for (AA in Substitute_AA$AA){
      aasub_ele_list<-Substitute_AA$aavector[Substitute_AA$AA==AA][[1]]
      if(Substitute_AA$Formula_with_water[Substitute_AA$AA==AA][[1]]){
        aasub_ele_list[["H"]]<-aasub_ele_list[["H"]]-2
        aasub_ele_list[["O"]]<-aasub_ele_list[["O"]]-1
      }
      #mode(aasub_ele_list)<- "integer"
      Element_tbl[[AA]]<-aasub_ele_list
    }
    
    for (AA in names(Element_tbl)){
      to_add<-Element_tbl[[AA]]
      to_add[aa_compo[(!(aa_compo %in% names(to_add)))]]=0
      to_add<-to_add[aa_compo]
      #print(to_add)
      Element_tbl[[AA]] <- to_add
    }  
  }
  
  return(Element_tbl)
}

ConvertPeptide<-function(sequence,Substitute_AA=NULL, Element_tbl=BuildElement(Substitute_AA)){
  peptideVector<- strsplit(sequence, split = "")[[1]]
  #L<-lapply(peptideVector, function(x,Element_tbl){Element_tbl[[x]]},Element_tbl)
  
  #resultsVector<-tapply(unlist(L), names(unlist(L)), sum)
  
  
  for (i in 1:length(peptideVector)) {
    if(i==1){
      resultsVector<-Element_tbl[[peptideVector[1]]]
    }else{
    to_add<-Element_tbl[[peptideVector[i]]]
    resultsVector <- to_add +  resultsVector
    }
  }
  
  
  resultsVector["H"] <- resultsVector["H"] + 2L
  resultsVector["O"] <- resultsVector["O"] + 1L
  
  resultsVector<- resultsVector[which(resultsVector!=0)]
  
  return(resultsVector)
}

ConvertPeptide.2<-function(sequence,Substitute_AA=NULL, Element_tbl=BuildElement(Substitute_AA,element_matrix = T)){
  peptideVector<- strsplit(sequence, split = "")[[1]]
  
  if (length(peptideVector)==1){
    resultsVector <- (Element_tbl[peptideVector,])
  }else{
    resultsVector <- colSums(Element_tbl[peptideVector,])
  }
  
  
  resultsVector["H"] <- resultsVector["H"] + 2L
  resultsVector["O"] <- resultsVector["O"] + 1L
  
  resultsVector<- resultsVector[which(resultsVector!=0)]
  
  return(resultsVector)
}

Peptide_Summary_para<- function(Proteins,peplist){
  
  tempdf1<- NULL
  for (Peptides in  1: length(peplist[[Proteins]])){
    tempdf1<- rbind(tempdf1,c(names(peplist)[Proteins],peplist[[Proteins]][Peptides]))
  }
  
  tempdf1
}

#' Peptide_modification
#'
#' This is a function that prepare the modification list for maldi imaging data qualitative or quantitative analysis.
#' this function will load a pre-built Unimod database (unimod.df) into global environment. 
#' @param retrive_ID the file name of candidate list
#' @param mod_position the folder that contains candidate list
#' @param update_unimod  the adducts list to be used for generating the PMF search candidates
#' @return a table of modification list 
#'
#' @examples
#' Peptide_modification()
#'
#' @export
Peptide_modification<-function(retrive_ID=NULL,mod_position=NULL,update_unimod=F){
   suppressMessages(suppressWarnings(require(protViz)))
   suppressMessages(suppressWarnings(require(XML)))
  if(update_unimod){
  message("Updating unimod database...")
  unimodurl <- url("http://www.unimod.org/xml/unimod_tables.xml")
  unimod.list <- XML::xmlToList(
    XML::xmlParse(
      scan(unimodurl, what = character())))
    save(unimod.list,file =paste0(path.package("HiTMaP"), "/data/unimod.list.rda"))
  
  } 
   
  if(!exists("unimod.df",envir = globalenv())){
      #try(data("unimod.list",package = "HiTMaP"))
  
  
    
  data("unimod.list",package = "HiTMaP")
  
  unimod.df<-lapply(unimod.list, function(x){
    x<-unname(x)
    x<-data.frame(do.call('rbind', lapply(x, function(y) y[match(names(x[[1]]), names(y))])),stringsAsFactors = F)
    if (".attrs" %in% colnames(x)) {
      attrs_tb<-data.frame(do.call('rbind', lapply(x$.attrs, function(y) y[match(names(x$.attrs[[1]]), names(y))])),stringsAsFactors = F)
      x<-cbind(x,attrs_tb)}
    return(x)
    })
  
  unimod.df <<-unimod.df
  
  }else{
    
    unimod.df <- get("unimod.df",envir = globalenv())
    
  }
  
  if (!is.null(retrive_ID)){
    retrive_mod<-data.frame(retrive_ID,mod_position,stringsAsFactors = F)
    unimod.modification.df<-merge(unimod.df$modifications,unimod.df$specificity,by.x=c("record_id"),by.y=c("mod_key"),all.x=T)
    unimod.modification.df=unimod.modification.df[unimod.modification.df$hidden==0,]
    unimod.df$positions_new<-unimod.df$positions
    unimod.df$positions_new$position_ext<-str_replace(unimod.df$positions_new$position,"Anywhere|Any N-term|Any C-term","Peptide")
    unimod.df$positions_new$position_ext<-str_replace(unimod.df$positions_new$position_ext,"Protein N-term|Protein C-term","Protein")
    
    unimod.modification.df<-merge(unimod.modification.df,unimod.df$positions_new,by.x="position_key",by.y="record_id")
    return_mod.df<-suppressWarnings(lapply(1:nrow(retrive_mod),function(x,retrive_mod,unimod.modification.df){
    if (is.na(as.numeric(retrive_mod$retrive_ID[x]))){
      unimod.modification.df['&'(unimod.modification.df$code_name==as.character(retrive_mod$retrive_ID[x]),unimod.modification.df$position_key==retrive_mod$mod_position[x]),]
    
      }else{
      unimod.modification.df[`&`(unimod.modification.df$record_id==as.character(retrive_mod$retrive_ID[x]),unimod.modification.df$position_key==retrive_mod$mod_position[x]),]
    
      }
    
    },retrive_mod,unimod.modification.df))
    
    return_mod.df=do.call(rbind,return_mod.df)
    
    #return_mod.df=merge(return_mod.df,unimod.df$specificity,by.x="record_id",by.y="mod_key",all.x=T)
    #return_mod.df=return_mod.df[return_mod.df$hidden==0,]
    return(return_mod.df)
  }

}
