

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