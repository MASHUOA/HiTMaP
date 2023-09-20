## As you already know we ideally want to exclude peptides with variably modified residues such as NQ (deamidation) M (oxidation), 
## N-terminal Q or E (pyro-glu formation) etc etc, and equally we are looking to avoid inadvertently selecting peptides 
## which are present as a variety of cleaved forms  - where those alternative forms have an intensity that is a significant 
## proportion of the “perfectly formed” sequence. 
## Here I am feeling that a “significant proportion” might be ~30% of the intensity of the perfect form. 
## When you manually conduct this process it is easy to sort the list of candidate peptides by sequence in order to 
## spot those peptides that are also represented by longer alternative forms on the C-terminus – 
## but sadly not on the N-terminus. Of course when there are lots of identified peptides to choose from for a particular protein 
## you can usually find about 5 good ones to select, 
## but when the protein of interest has only a few options present you tend to have to be prepared to 
## bend the rules a little to at least deliver something to the researcher… 
## Perhaps what we need in those cases is a way to mark such choices as being suboptimal.


PRM_target_slt<-function(GroupID,ID_data,TopN_Feat=5,wd=getwd(),
                         Id_pipeline=c("ProteinPilot","Fragpipe", "UserTable"),
                         Consider_RT=F,Consider_intensities=T,ppm_cutoff=15,
                         ref_fasta=NULL,Uniprot_mapping=T,Peptideatlas_mapping=T, 
                         pep_length=c(7,40),score_cutoff=50,
                         recluster_by_mis=T){
  
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
  if (!dir.exists(wd)) dir.create(wd,recursive = T)
  setwd(wd)
  final_res<-NULL
  
  if ('|'(missing(GroupID),is.null(GroupID))) stop("No input target found...")
  if ('|'(missing(ID_data),is.null(ID_data))) stop("No identification data found...")
    
  if (is.character(GroupID)){
    if( '&'(length(GroupID)==1,sum(file.exists(GroupID))==1) ){
      read_lines(GroupID)->GroupID
    }else{
      GroupID=GroupID
    }
  } 
  
  if (is.character(ID_data)){
    if( '&'(length(ID_data)==1,sum(file.exists(ID_data))==1) ){
      if (Id_pipeline[1]=="ProteinPilot") {
        ID_data=read_excel(ID_data, sheet = "Processed Peptide Summary")
        required_col<-c(Pro_col<-"Accessions",
                            Pep_col<-"Sequence" ,
                            Pep_mod_col<-"Modifications",
                            RT_col<-"Acq.Time",
                            Cleavage_col<-"Cleavages",
                            Chanrge_col<-"Theor.z",
                            Score_col<-"Conf",
                            Intensities_col<-"Intensity.Peptides")
      }
      if (Id_pipeline[1]=="Fragpipe") {
        ID_data=read_delim(ID_data, 
                                                         delim = "\t", escape_double = FALSE, 
                                                         trim_ws = TRUE)
        required_col<-c(Pro_col<-"Protein.Ids",
                         Pep_col<-"Stripped.Sequence" ,
                         Pep_mod_col<-"Modified.Sequence",
                         RT_col<-"RT",
                         Cleavage_col<-"Proteotypic",
                         Chanrge_col<-"Precursor.Charge",
                         Score_col<-"Global.Q.Value",
                         Intensities_col<-"Ms1.Area")
      }
         
      if (Id_pipeline[1]=="UserTable") {
        required_col<-c(Pro_col<-"Protein",
                           Pep_col<-"Sequence" ,
                           Pep_mod_col<-"Modification",
                           RT_col<-"RT",
                           Cleavage_col<-"Cleavages",
                           Chanrge_col<-"Charge",
                           Score_col<-"Score",
                           Intensities_col<-"Area")
        ID_data=read.csv(ID_data)
      }
      
    } else {
      stop("Identification data is invalid.")
    }
    
  } else if (is.data.frame(ID_data)){
      GroupID=GroupID
  } else {
      stop("Identification data is invalid.")
  }
  
  
  
  if (sum(!(required_col %in% colnames(ID_data)))>=1){
    stop("Required column(s) is missing in identification data:", paste0(required_col[!(required_col %in% colnames(ID_data))],collapse = ","))
  }
  
  ID_data=ID_data[,required_col]
  colnames(ID_data)<-c("Protein","Peptide","Modification","RT","Cleavages","Charge","Score","Intensity")
  ID_data$Intensity<-as.numeric(ID_data$Intensity)
  ID_data$Cleavages[is.na(ID_data$Cleavages)] <- ""
  ID_data<-ID_data[ID_data$Score>=score_cutoff,]
  for(GroupID_slt in GroupID){
  
    ID_data_slt<-ID_data[ID_data[,"Protein"]==GroupID_slt,]
    if (length(unique(str_split_fixed(unlist(str_split(unique(ID_data_slt$Modification),"; ")),"@",2)[,1]))==0){
      next
    }
    ID_data_slt<-ID_data_slt[ between(width(ID_data_slt$Peptide),pep_length[1],pep_length[2]),]
    
    # Biostrings::BStringSet(unique(ID_data_slt$Peptide))->BStringSet_pep
    # msa(BStringSet_pep,type="protein",method="ClustalW")->BStringSet_pep
    # 
    # tmpFile <- tempfile(pattern="msa", tmpdir=".", fileext=".pdf")
    # msaPrettyPrint(BStringSet_pep, file=tmpFile, output="pdf",
    #                showNames="left", showNumbering="none", showLogo="top",
    #                showConsensus="bottom", logoColors="rasmol",
    #                verbose=FALSE, askForOverwrite=FALSE)
    # paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),"_aligment.pdf")
    
    if(Uniprot_mapping){
    #Uniprotxml<-curl_download(paste0("https://www.ebi.ac.uk/proteins/api/proteins/",get_pro_acc(GroupID_slt,acc_loc="mid"),"?format=xml"),paste0(which(GroupID==GroupID_slt),".xml"))
      if (!file.exists(paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),".json")))  {
        Uniprotjson<-curl_download(paste0("https://www.ebi.ac.uk/proteins/api/proteins/",get_pro_acc(GroupID_slt,acc_loc="mid")),paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),".json"))
      }else{
        Uniprotjson<-paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),".json")
      }
      
    jsonlite::read_json(Uniprotjson)->Uniprotjson
    Uniprotjson[["sequence"]][["sequence"]]->pro_map_seq
    ID_data_slt_map<-str_locate(pro_map_seq,ID_data_slt$Peptide)
    ID_data_slt<-cbind(ID_data_slt,ID_data_slt_map)
    }
    
    ID_data_slt->ID_data_slt_withmiss
    ID_data_slt_withmiss<-ID_data_slt_withmiss[str_detect(ID_data_slt_withmiss$Cleavages,regex("missed|cleaved",ignore_case = T)),]
    
    ID_data_slt<-ID_data_slt[!str_detect(ID_data_slt$Cleavages,regex("missed|cleaved",ignore_case = T)),]
    ir1<-IRanges(start=ID_data_slt$start,end=ID_data_slt$end)
    ir1_r<-reduce(ir1,min.gapwidth=0L)
    ol<-findOverlaps(ir1,ir1_r)
    ir1_acc<-unlist(lapply(1:length(ir1),function(x) {ir1[x,]@start:(ir1[x,]@start+ir1[x,]@width)}))
    ir1_acc<-table(ir1_acc)
    ir1_acc_all<-data.frame(ir1_acc=1:width(pro_map_seq))
    ir1_acc_all<-merge(ir1_acc_all,ir1_acc,all.x=T)
    ID_data_slt$Cluster=as.character(ol@to)
    ID_data_slt<-ID_data_slt[order(ID_data_slt$Cluster),]
    ID_data_slt$Modification[is.na(ID_data_slt$Modification)]<-""
    ID_data_slt$Cleavages[is.na(ID_data_slt$Cleavages)]<-""
    
    #plotRanges(ir1,main=get_pro_acc(GroupID_slt,acc_loc="mid"),width = 1800,height = 1200,res = 300)->img
    #image_write(img,path=paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),"_range.png"))
    
    #plotRanges(ir1_r,main=get_pro_acc(GroupID_slt,acc_loc="mid"),width = 1800,height = 1200,res = 300)->img
    #image_write(img,path=paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),"_reduced_range.png"))

    
    # plotRanges(ir1,main=get_pro_acc(GroupID_slt,acc_loc="mid"),width = 1800,height = 1200,res = 300)->img
    # image_draw(img)->img_ol
    # plot(ir1_acc_all$ir1_acc,ir1_acc_all$Freq,col="red")
    # lines(ir1_acc_all$ir1_acc,ir1_acc_all$Freq,col="red")
    # dev.off()
    # image_write(img_ol,path=paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),"_ol_range.png"))
    
    
    
    if (recluster_by_mis){
      ID_data_slt_new<-NULL
        for (clusterID in unique(ID_data_slt$Cluster)){
          
          ID_data_slt[ID_data_slt$Cluster==clusterID,]->ID_data_slt_temp
          ID_data_slt_withmiss['&'(ID_data_slt_withmiss$start<=min(ID_data_slt_temp$start),ID_data_slt_withmiss$end>=max(ID_data_slt_temp$end)),]->ID_data_slt_withmiss_temp 
          
          miss_row_num=sum(str_detect(ID_data_slt_temp$Cleavages,regex("missed|cleaved",ignore_case = T)),na.rm = T)
          if('&'(miss_row_num>=1,miss_row_num!=nrow(ID_data_slt_temp))){
            ID_data_slt_withmiss_temp$Cluster=clusterID
            ID_data_slt_temp<-rbind(ID_data_slt_temp,ID_data_slt_withmiss_temp)
            
            ID_data_slt_temp->ID_data_slt_temp_org
            ID_data_slt_temp<-ID_data_slt_temp[!str_detect(ID_data_slt_temp$Cleavages,"missed"),]
            ID_data_slt_temp_org_miss<-ID_data_slt_temp_org[str_detect(ID_data_slt_temp_org$Cleavages,"missed"),]
            ir1<-IRanges(start=ID_data_slt_temp$start,end=ID_data_slt_temp$end)
            ir1_r<-reduce(ir1,min.gapwidth=0L)
            ol<-findOverlaps(ir1,ir1_r)
            ir1_acc<-unlist(lapply(1:length(ir1),function(x) {ir1[x,]@start:(ir1[x,]@start+ir1[x,]@width)}))
            ir1_acc<-table(ir1_acc)
            ir1_acc_all<-data.frame(ir1_acc=1:width(pro_map_seq))
            ir1_acc_all<-merge(ir1_acc_all,ir1_acc,all.x=T)
            ID_data_slt_temp$Cluster=ol@to
            ID_data_slt_temp<-ID_data_slt_temp[order(ID_data_slt_temp$Cluster),]
            
            for (x in unique(ID_data_slt_temp$Cluster)){
              ID_data_slt_temp_org_miss_x<-ID_data_slt_temp_org_miss
              ID_data_slt_temp_org_miss_x$Cluster<-x
              ID_data_slt_temp<-rbind(ID_data_slt_temp,ID_data_slt_temp_org_miss_x)
            }
            ID_data_slt_temp$Cluster<-paste0(clusterID,"_",ID_data_slt_temp$Cluster)
            ID_data_slt_new[[clusterID]]<-ID_data_slt_temp
          }else{
            ID_data_slt_new[[clusterID]]<-ID_data_slt_temp
          }
        }
      ID_data_slt_new<-do.call(rbind,ID_data_slt_new)
      ID_data_slt_new->ID_data_slt
    }

    
    
    
    ID_data_slt %>% group_by(Cluster) %>% mutate(Cluster_max_intensity=max(Intensity)) -> ID_data_slt
    ID_data_slt$Cluster_intensity_percentage<-(ID_data_slt$Intensity/ID_data_slt$Cluster_max_intensity)
    ID_data_slt %>% group_by(Cluster,Peptide,Modification) %>% mutate(Other_greater_than30='&'(max(Cluster_intensity_percentage)>=0.3,max(Cluster_intensity_percentage)!=1)) -> ID_data_slt
    ID_data_slt %>% group_by(Cluster) %>% mutate(No_other_greater_than30=!any(Other_greater_than30)) -> ID_data_slt
    
    
    
    
    No_other_greater_than30 <- ID_data_slt$No_other_greater_than30
    modification_slt <- !str_detect(ID_data_slt$Modification,"Deamidated\\(N\\)|Oxidation\\(M\\)|Glu->pyro-Glu")
    modification_slt[is.na(modification_slt)]<-T
    cleavage_slt <- !str_detect(ID_data_slt$Cleavages,"missed")
    cleavage_slt[is.na(cleavage_slt)]<-T
    ID_data_slt<-cbind(ID_data_slt,cleavage_slt=cleavage_slt)
    ID_data_slt<-cbind(ID_data_slt,modification_slt=modification_slt)
    
    AND1 <- function (...)  Reduce("&", list(...))
    ID_data_slt$Decision="Exclude"
    ID_data_slt$Decision[AND1(No_other_greater_than30,modification_slt,cleavage_slt)]<-"select"
    
    ID_data_slt %>% group_by(Cluster,Decision) %>% mutate(Cluster_rank=order(Intensity,decreasing = T)) -> ID_data_slt
    
    write.csv(ID_data_slt,paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),"_decision.csv"),row.names = F)
    
    ID_data_slt_reduce<-ID_data_slt[ID_data_slt$Decision=="select",]
    ID_data_slt_reduce<-ID_data_slt_reduce[ID_data_slt_reduce$Cluster_rank==1,]
    
 
  if(Peptideatlas_mapping){
    pep_atlas_info<-NULL
    for (pep in unique(ID_data_slt_reduce$Peptide)){
    if (!file.exists(paste0(pep,".tsv")))  {
    Patlasurl = paste0('https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetPeptide?',
                       'atlas_build_id=' , 446 ,
                       '&searchWithinThis=Peptide+Sequence' ,
                       '&searchForThis=' , pep , '&action=QUERY')
    curl_download(Patlasurl,paste0(pep,".html"))
    
    readLines(paste0(pep,".html"))->Patlasurl
    
    Patlasurl[str_detect(Patlasurl,"Download as: .{10,10000}tsv")]->tsv_line
    str_extract(tsv_line,"href=.{1,10000} TITLE=")->tsv_line
    tsv_line<-str_remove_all(tsv_line,"href=.|. TITLE=")
    curl_download(paste0("https://db.systemsbiology.net/",tsv_line),paste0(pep,".tsv"))
    }
    
    peptide_data<-read_delim(paste0(pep,".tsv"), 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)
    table(peptide_data$Instr)->Instru_stats
    all=sum(Instru_stats)
    Instru_stats<-paste(names(Instru_stats),Instru_stats,sep = ":",collapse = "; ")
    Instru_stats<-paste0("All:",all,"; ",Instru_stats)
    pep_atlas_info[[pep]]<-Instru_stats
    }
    ID_data_slt_reduce$Peptide_atlas<-unlist(pep_atlas_info[match(ID_data_slt_reduce$Peptide,names(pep_atlas_info))])
  }
    
    ID_data_slt_reduce<-ID_data_slt_reduce[order(ID_data_slt_reduce$Intensity,decreasing = T),]
    ID_data_slt_reduce<-ID_data_slt_reduce[1:TopN_Feat,]
    
    write.csv(ID_data_slt_reduce,paste0(get_pro_acc(GroupID_slt,acc_loc="mid"),"_final.csv"),row.names = F)
    
    final_res[[GroupID_slt]]<-ID_data_slt_reduce
    
  }


  final_res_bind<-do.call(rbind,final_res)
  
  write.csv(final_res_bind,paste0("all_selected","_final.csv"),row.names = F)
  
}



get_pro_acc<-function(accession,acc_loc=c("mid","end")){
  library(stringr)
  accession_org<-accession
  if(acc_loc=="end"){
    accession<-str_extract(accession,"\\|.{3,10}_")
    accession<-stringr::str_remove_all(accession,"_")
    accession<-stringr::str_remove_all(accession,"\\|")
  } else if (acc_loc=="mid"){
    accession<-str_extract(accession,"\\|.{3,10}\\|")
    accession<-stringr::str_remove_all(accession,"_")
    accession<-stringr::str_remove_all(accession,"\\|")
  }
  accession[is.na(accession)]<-accession_org[is.na(accession)]
  
  return(accession)
}

plotRanges<-function(x,xlim=x,main=deparse(substitute(x)), col="black",sep=0.5,...) { 
  library(magick)
  magick::image_graph(...)->img
  plot.new()
  height<-1 
  if(is(xlim,"IntegerRanges"))
    xlim<-c(min(start(xlim)),max(end(xlim))) 
  bins<-disjointBins(IRanges(start(x),end(x)+1)) 
 
  plot.window(xlim,c(0,max(bins)*(height+sep))) 
  ybottom<-bins*(sep+height)-height 
  rect(start(x)-0.5,ybottom,end(x)+0.5,ybottom+height,col=col,...) 
  title(main)
  axis(1) 
  dev.off()
  return(img)
}

