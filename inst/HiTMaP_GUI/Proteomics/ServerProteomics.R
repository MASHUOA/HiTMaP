library(magrittr)
library(shiny)
library(shinyFiles)
library(future)
plan(multisession)
options(shiny.maxRequestSize = 3000*1024^2)

parse_mod<-function(mod){
    library(stringr)
    if(!is.null(mod)){
    return_mod<-list()
    return_mod[[1]]<-as.numeric(unlist(lapply(str_split(mod,": "),function(x){
        x[1]})))
    return_mod[[2]]<-unlist(lapply(str_split(mod,"@ "),function(x){
        str_split(x[2],":")[[1]][2]
        }))
    return_mod<-as.data.frame(return_mod)
    colnames(return_mod)=c(1,2)
    return(return_mod)
    }else{
    return(data.frame("1"=NULL,"2"=NULL))  
    }

    }
