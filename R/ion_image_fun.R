plotly_for_region_with_ClusterID_barchart<-function(moleculeNames,data){
  library(dplyr)
  library(plyr)
  library(plotly)
  library(stringr)
  library(ggplot2)
  if (!require("processx")) install.packages("processx")
  
  plot_data=data[data$moleculeNames==moleculeNames,data["Class",]!=""]
  colnames(plot_data)=as.character(t(data["Class",data["Class",]!=""]))
  
  rownames(plot_data)=data[data$moleculeNames==moleculeNames,"RegionName"]
  
  
  stderr <- function(x, na.rm=FALSE) {
    if (na.rm) x <- na.omit(x)
    sqrt(var(x)/length(x))
  }
  plot_data=mutate_all(plot_data, function(x) as.numeric(as.character(x)))
  for (row in rownames(plot_data)){
    plot_data[row,"mean"]=mean(as.numeric(plot_data[row,as.character(t(data["Class",data["Class",]!=""]))]))
    plot_data[row,"se"]=stderr(as.numeric(plot_data[row,as.character(t(data["Class",data["Class",]!=""]))]))
  }
  rownames(plot_data)=data[data$moleculeNames==moleculeNames,"RegionName"]
  plot_data$mz=data[data$moleculeNames==moleculeNames,"mz"]
  plot_data$RegionRank=as.numeric(as.character(data[data$moleculeNames==moleculeNames,"RegionRank"]))
  plot_data=plot_data[order(-plot_data$RegionRank),]
  plot_data$region= as.factor(as.character(rownames(plot_data)))
  if (sum(plot_data$region==c("OC" ,  "IC" ,  "Core"))==length(plot_data$region)){plot_data$region=factor(plot_data$region,levels = c("OC" ,  "IC" ,  "Core"))}
  #data[data$mz==plot_data$mz[1],"moleculeNames"]
  p <- plot_ly(data = plot_data, x = ~region, y = ~mean, type = 'bar', name = moleculeNames,
               error_y = ~list(array = se,
                               color = '#000000')) %>%
    layout(title = paste(plot_data$mz[1]),
           xaxis = list(title = "Region"),
           yaxis = list (title = "Relative Conc.",ticksuffix = "%"),
           showlegend = F)
  #%>%
  #  add_trace(data = data[which(data$supp == 'VC'),], name = 'VC')
  
  
  
  
  windows_filename<- function(stringX){
    stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
    stringX<-gsub("\"", "", stringX)
    return(stringX)}
  #moleculeNames=colnames(data)
  #fit3 <- lm(as.numeric(as.character( data[which(data$Class!=""),moleculeNames[i]]))~poly(as.numeric(as.character(data$Class[which(data$Class!="")])),3) )
  
  
  
  htmlwidgets::saveWidget(as_widget(p), windows_filename(paste0(moleculeNames,".html")), selfcontained = F, libdir = "lib")
  
  #plotly::orca(p, windows_filename(paste0(data["moleculeNames",i], data["adducts",i]," in ",data["RegionName",i],".png")))
  windows_filename(paste0(moleculeNames,".html"))
}