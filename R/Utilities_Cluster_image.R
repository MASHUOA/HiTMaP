



#' cluster_image_grid
#'
#' This function renders the clustered images for maldi imaging data set
#'
#' @export

cluster_image_grid<-function(clusterID,
                             SMPLIST,
                             imdata,
                             ClusterID_colname="Protein",
                             componentID_colname="Peptide",
                             Component_plot_threshold=2,
                             Component_plot_coloure=c("mono","as.cluster"),
                             smooth.image="gaussian",
                             contrast.enhance = "suppression",
                             colorpallet="Set1",
                             plot_layout=c("line","grid"),
                             export_Header_table=F,
                             export_footer_table=F,
                             plot_style=c("fleximaging","ClusterOnly","rainbow"),
                             protein_coverage=F,
                             footer_style="Length",
                             output_png_width_limit=198000,
                             attach_summary_cluster=T,
                             cluster_color_scale=c("blackwhite","fleximaging"),
                             remove_cluster_from_grid=T,
                             img_brightness= 100,ppm=20,
                             list_of_protein_sequence,
                             workdir=getwd(),
                             pixel_size_um=50,
                             Score_thres=NULL){
  #complementary(color="red", plot = TRUE, bg = "white", labcol = NULL, cex = 0.8, title = TRUE)
  windows_filename<- function(stringX){
    stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
    stringX<-gsub("\"", "", stringX)
    stringX<-gsub("\\|", " ", stringX)
    return(stringX)

  }
  
   suppressMessages(suppressWarnings(require(magick)))
   suppressMessages(suppressWarnings(require(stringr)))
   suppressMessages(suppressWarnings(require(data.table)))
   #Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
   suppressMessages(suppressWarnings(require(grid)))
   suppressMessages(suppressWarnings(require(plotly)))
   suppressMessages(suppressWarnings(require(dplyr)))
   suppressMessages(suppressWarnings(if(!require(colortools)) {
     library(devtools)
     install_github('gastonstat/colortools')
   }))
   suppressMessages(suppressWarnings(require(data.table)))
   suppressMessages(suppressWarnings(require(Cardinal)))
  #rotate the image
  #imdata@pixelData@data<-rotatetmp
  outputpngsum=paste(workdir,"/",windows_filename(substr(clusterID, 1, 15)),"_cluster_imaging",'.png',sep="")

  #message(outputpngsum)
  SMPLIST=as.data.frame(SMPLIST)
  outputpng=paste(workdir,"/",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot",'.png',sep="")
  #message(outputpng)
  candidate=SMPLIST[SMPLIST[[ClusterID_colname]]==clusterID,]

  if(!is.null(Score_thres)){candidate=candidate[candidate$Score>=Score_thres,]}

  #candidate=candidate[order(as.character())]
  if (componentID_colname %in% colnames(candidate)){
  candidate_u<- candidate %>% dplyr::group_by(mz) %>% dplyr::summarise(Peptide=Peptide[1], .groups = 'drop')

  candidate<-merge(data.frame(candidate_u),candidate,by=c("mz",componentID_colname),sort=F)
  }
  candidateunique=as.numeric(as.character(unique(candidate[,"mz"])))
  candidateunique=candidateunique[order(as.character(candidateunique))]

  #check if candidates' mz within the dataset range
  candidateunique<-candidateunique[which(between(candidateunique,range(mz(imdata))[1],range(mz(imdata))[2]))]
  candidate<-candidate[candidate$mz %in% candidateunique,]

  if(is.null(candidate$desc)){
    candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
    cluster_desc<-"-"
  }else{
    candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
    cluster_desc<-unique(candidate$desc)[1]
  }
  if (length(candidateunique)>4){
 
    mycol=setColors("red", length(candidateunique))#wheel("red", num = length(candidateunique),bg = "white")
  } else if (length(candidateunique)==4){
    #tmp_cols = setColors("red", 12)
    #tetrad_colors <- tmp_cols[c(1, 3, 7, 9)]
    mycol=c("#FF0000", "#FFFF00", "#00FFFF", "#0000FF")
    #mycol=c("#FF0000", "#00FF00", "#0000FF")
  } else if (length(candidateunique)==3){
    #mycol=splitComp("red")
    mycol=c("#FF0000", "#00FF00", "#0000FF")
  } else if (length(candidateunique)==2){
    #mycol=complementary("red")
    mycol=c("#FF0000", "#00FFFF")
  } else if (length(candidateunique)==1){
    #mycol=complementary("red")
    mycol=c("#FF0000")
  }

  mycol <- as.factor(as.character(mycol))
  mycol=mycol[order(mycol)]
  #mycolrgb<-hex2RGB(mycol)
  #meancoldf<-colSums(mycolrgb@coords)/max(colSums(mycolrgb@coords))
  #meancol<-RGB(R=meancoldf[1],G=meancoldf[2],B=meancoldf[3])
  #colorspace::hex(meancol)
  #print(length(candidateunique))
  #print(Component_plot_threshold)
  #print(length(candidateunique)>=Component_plot_threshold)
  if (length(candidateunique)>=Component_plot_threshold){

    if (is.null(imdata)){
      message("No imaging data")

    }else
      {
       suppressMessages(suppressWarnings(require(RColorBrewer)))
       suppressMessages(suppressWarnings(require(Cardinal)))
       suppressMessages(suppressWarnings(require(EBImage)))
      #library(colortools)

      if (plot_style=="ClusterOnly"){
        png(outputpng,width = 10,height = 10, bg = "black",units = "in",res = 300)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,1),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)
        Cardinal::image(imdata, mz=candidateunique,
              col=mycol,
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              tolerance=ppm, units="ppm")


        l<-function(x,y,z,t=x+y,f=t+z){
          paste(t,f)
        }
        l(1,2,3)

        dev.off()
        pngfile<-image_read(outputpng)
        pngfile<-image_border(pngfile, "black", "30x30")
        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        pngfile<-image_trim(pngfile)
        image_write(pngfile,outputpng)

      }
      else if (plot_style=="rainbow"){


        if (plot_layout=="line"){
          png(outputpngsum,width = 5*((length(candidateunique)+1)),height = 5, bg = "black",units = "in",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)


        }else{
          png(outputpngsum,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)
        }




        Cardinal::image(imdata, mz=candidateunique,
              col=levels(mycol),
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              tolerance=ppm, units="ppm",key=F)

        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          Cardinal::image(imdata, mz=candidateunique[i],
                contrast.enhance=contrast.enhance,
                smooth.image = smooth.image ,
                col.regions=col.regions,

                normalize.image="none",
                tolerance=ppm, units="ppm")
          componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          for (component in componentname){
            text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          }
        }
        dev.off()
        #pngfile<-image_read(outputpngsum)
        #pngfile<-image_border(pngfile, "black", "30x30")
        #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        #pngfile<-image_trim(pngfile)
        #image_write(pngfile,outputpngsum)

        

      }
      else if(plot_style=="fleximaging"){

        ##################################
        matter::vizi_style(style = "light", dpal = "Tableau 10", cpal = "Viridis")
        tmp_dir <- tempdir()
        temp_cluster_png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        temp_component_png=list()
        temp_component_png.mono=list()
        if (plot_layout=="line"){
          #png(temp_cluster_png,width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = "black",units = "in",res = 300)
          temp_cluster_png<-tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
          png(temp_cluster_png,width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = "black",units = "px",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)


        }else{
          #png(temp_cluster_png,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 300)
          temp_cluster_png<-tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
          png(temp_cluster_png,width = 750*2,height = 750 * (ceiling((length(candidateunique)+1)/2)), bg = "black",units = "px",res = 150)

          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)
        }


        col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[1])
        #col.regions <- cpal(palette = "Viridis")


        dev.off()
        #if (length(candidateunique)==2){
        #  bg = "grey9"
        #}else if (length(candidateunique)==3){
        #  bg = "grey8"
        #}else if (length(candidateunique)==3){
        #
        #}
        bg = paste0("grey",ifelse((15-length(candidateunique))>1,(15-length(candidateunique)),2))
        #bg = "transparent"
        componentimg=list()
        componentimg.mono=list()
        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          if  (Component_plot_coloure=="mono"){
            bg = paste0("grey",29)
            col.regions.mono=intensity.colors_customize1(colset = 2)
            col.regions=create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          }else if(Component_plot_coloure=="as.cluster"){
            col.regions=create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
            #col.regions=create_gradient_colors_robust(100, start="#030303", end=levels(mycol)[i])
            col=levels(mycol)[i]
          }

          if(cluster_color_scale=="blackwhite"){
            col.regions.mono=create_gradient_colors_robust(100, start="black", end="white")
            col.regions=create_gradient_colors_robust(100, start="black", end="white")
          }
          if (Component_plot_coloure=="mono"){
            
            componentimg.mono[[i]]=Cardinal::image(imdata, mz=candidateunique[i],
                                    #contrast.enhance=contrast.enhance,
                                    contrast.enhance = "none",
                                    smooth.image = smooth.image,
                                    colorscale=col.regions.mono,
                                    normalize.image="linear",
                                    tolerance=ppm, units="ppm",
                                    key=F,
                                    xlab=NULL,
                                    layout=c( length(levels(Cardinal::run(imdata))),1),
                                    bg = bg)
            
            
            temp_component_png[[i]]=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
            
            componentimg.mono[[i]][["par"]][["ann"]]=F
            componentimg.mono[[i]][["par"]][["bty"]]="n"
            componentimg.mono[[i]][["par"]][["pty"]]="s"
            componentimg.mono[[i]][["par"]][["xaxt"]]="n"
            componentimg.mono[[i]][["par"]][["yaxt"]]="n"
            componentimg.mono[[i]][["par"]][["fg"]]="white"
            componentimg.mono[[i]][["par"]][["oma"]]=c(0, 0, 0, 0)
            #componentimg.mono[[i]][["par"]][["mar"]]=c(0, 0, 0, 1)
            if (!is.null(componentimg.mono[[i]][["facets"]]) && 
                length(componentimg.mono[[i]][["facets"]]) > 0 &&
                !is.null(componentimg.mono[[i]][["facets"]][[1]])) {
              attr(componentimg.mono[[i]][["facets"]][[1]],"strip")$strip=F
              attr(componentimg.mono[[i]][["facets"]][[1]],"colorkey")$colorkey=T
            }

            # componentimg.mono[[i]][["facets"]][[1]][[1]][["dpage"]]="a"
            # componentimg.mono[[i]][["facets"]][[1]][[1]][["facet"]][[".feature.groups"]]="a"
            # componentimg.mono[[i]][["fids"]][[".feature.groups"]]="a"
            # componentimg.mono[[i]][["dpages"]]="a"
            #png(temp_component_png.mono[[i]],width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = bg,units = "in",res = 300)
            png(temp_component_png[[i]],width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = bg,units = "px",res = 150)
            par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
                #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
                bty="n",pty="s",xaxt="n",
                yaxt="n",
                no.readonly = TRUE,ann=FALSE)

            try(tryCatch(print(componentimg.mono[[i]])),silent = T)
            dev.off()
          }else{
        componentimg[[i]]=Cardinal::image(imdata, mz=candidateunique[i],
                                      #contrast.enhance=contrast.enhance,
                                      contrast.enhance = "none",
                                      smooth.image = smooth.image,
                                      colorscale=col.regions,
                                      normalize.image="linear",
                                      tolerance=ppm, units="ppm",
                                      key=F,
                                      xlab=NULL,
                                      layout=c( length(levels(Cardinal::run(imdata))),1)
                                      )


        temp_component_png[[i]]=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

        #png(temp_component_png[[i]],width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = bg,units = "in",res = 300)
        png(temp_component_png[[i]],width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = bg,units = "px",res = 150)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
            #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)

        try(tryCatch(print(componentimg[[i]])),silent = T)

        dev.off()
          }




          #componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          #for (component in componentname){
          #  text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          #}
        }

        temp_component_cover=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

        #png(temp_component_cover,width = 5,height = 0.5, bg = bg,units = "in",res = 300)
        temp_component_cover<-tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        png(temp_component_cover,width = 750,height = 750, bg = bg,units = "px",res = 150)
        plot.new()
        text("")
        #temp_component_cover<-image_draw(temp_component_cover)
        dev.off()
        temp_component_cover<-image_read(temp_component_cover)
        temp_component_cover<-image_crop(temp_component_cover, "750x80+0")
        #temp_component_cover<-image_background(temp_component_cover, "hotpink")
        makeacover<-(temp_component_cover)
        #print(makeacover)
        #pngfile<-(temp_cluster_png)
        #pngfile<-image_border(pngfile, bg, "30x30")
        #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "north",size = 50,color = "white")

        pngcompfile_org=list()
        pngcompfile=list()
        for (i in 1:length(candidateunique)){
          pngcompfile_org[[i]]<-image_read(temp_component_png[[i]])
          pngcompfile_org[[i]]<-image_flatten(c(pngcompfile_org[[i]],makeacover))
          #if (Component_plot_coloure=="mono"){
          #pngcompfile[[i]]<-(temp_component_png.mono[[i]] )
          #pngcompfile[[i]]<-image_flatten(c(pngcompfile[[i]],makeacover))
          #}else{
          pngcompfile[[i]]<-pngcompfile_org[[i]]
          #}

          pngcompfile[[i]]<-image_border(pngcompfile[[i]], bg, "30x30")
          #pngcompfile[[i]]<-image_annotate(pngcompfile[[i]],paste(unique(candidate[candidate$mz==candidateunique[i],"moleculeNames"]),candidateunique[i]),gravity = "north",size = 30,color = "white")

        }




        pngcompfile_output<-pngcompfile[[1]]
        if (length(pngcompfile)>=2){
        for (i in 2: length(pngcompfile)){
          pngcompfile_output<-c(pngcompfile_output,unlist(pngcompfile[[i]]))
        } }

        #img_com<-(temp_component_png[[1]])
        #img_com<-image_flatten(c(img_com,makeacover))
        #pngcompfile_org[[1]]<-NULL
        img_com<-pngcompfile_org[[1]]
        if (length(pngcompfile)>=2){
        for (i in 2: length(pngcompfile_org)){
          img_com<-c(img_com,unlist(pngcompfile_org[[i]]))
        }}
        #img_com<-rbind(unlist(pngcompfile_org))


        #lapply(channel_types(),function(x,pngfile){
        #  pngfile %>% image_threshold(type = "white", threshold = "50%",channel = x) %>% image_scale( "x600") %>% image_write(path = paste0(x,".png"))
        #},pngfile)
        cluster_desc<-gsub(stringr::str_extract(cluster_desc,"OS=.{1,}"),"",cluster_desc)

       if ((Component_plot_coloure=="mono")) {

         pngfile_big<-image_average(img_com)
         pngfile_big_info=magick::image_info(pngfile_big)

         pngfile_big_ratio<-c(diff(range(imdata@elementMetadata@listData[["x"]])),diff(range(imdata@elementMetadata@listData[["y"]])))/c(pngfile_big_info$width-1*150,pngfile_big_info$height-1*150)



         pngfile_big <- image_draw(pngfile_big)
         #rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
         pixel_size_um_bar<-200
         pixel_size<-pixel_size_um_bar*pixel_size_um/(1/max(pngfile_big_ratio))/1000
         if(!is.integer(pixel_size)){
           pixel_size<-ceiling(pixel_size)
           pixel_size_um_bar<-pixel_size*1000/pixel_size_um*(1/max(pngfile_big_ratio))
         }

         rect(690-pixel_size_um_bar-60, 680,690-pixel_size_um_bar-60+pixel_size_um_bar,  680+13, border  = 'white',col  = 'white', lwd = 0)
         text(690, 690, paste0(pixel_size ," mm"),col  = 'white', cex = 3, srt = 0)
         #palette(rainbow(11, end = 0.9))
         #symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
         #        bg = 1:11, inches = FALSE, add = TRUE)
         dev.off()


         pngfile<-pngfile_big
         pngfile_big<-image_border(pngfile_big, bg, "30x30")
         #pngfile_big<-image_modulate(pngfile_big, brightness = 100 + 25 * length(candidateunique), saturation = 100)
         #pngfile_big<-image_threshold(pngfile_big,type = "white", threshold = "50%",channel = "All")
         #pngfile_big<-image_annotate(pngfile_big,paste(cluster_desc),gravity = "north",size = 50,color = "white")
         pngfile_big<-image_annotate(pngfile_big,cluster_desc,gravity = "north",size = 30*40/nchar(cluster_desc),color = "white")

         pngfile_big<-image_annotate(pngfile_big,"0%        Relative intensity        100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
         #pngfile_big<-image_annotate(pngfile_big,"0                                               100",location = "+52+165",gravity = "northeast",size = 30,color = "white",degree=270)
         #pngfile_big<-image_annotate(pngfile_big,"                 Relative intensity                  ",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)

         #pngfile<-image_average(img_com)
         #pngfile<-image_threshold(pngfile, type = "black",  threshold = "50%")
         pngfile<-image_border(pngfile, bg, "30x30")
         #pngfile<-image_modulate(pngfile, brightness = 150)


         }
        else if ((Component_plot_coloure=="as.cluster")){
        pngfile_big<-image_average(img_com)
        pngfile_big_info=magick::image_info(pngfile_big)

        pngfile_big_ratio<-c(diff(range(imdata@elementMetadata@listData[["x"]])),diff(range(imdata@elementMetadata@listData[["y"]])))/c(pngfile_big_info$width-1*150,pngfile_big_info$height-1*150)



        pngfile_big <- image_draw(pngfile_big)
        #rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
        pixel_size_um_bar<-200
        pixel_size<-pixel_size_um_bar*pixel_size_um/(1/max(pngfile_big_ratio))/1000
        if(!is.integer(pixel_size)){
          pixel_size<-ceiling(pixel_size)
          pixel_size_um_bar<-pixel_size*1000/pixel_size_um*(1/max(pngfile_big_ratio))
        }

        rect(690-pixel_size_um_bar-60, 680,690-pixel_size_um_bar-60+pixel_size_um_bar,  680+13, border  = 'white',col  = 'white', lwd = 0)
        text(690, 690, paste0(pixel_size ," mm"),col  = 'white', cex = 3, srt = 0)
        #palette(rainbow(11, end = 0.9))
        #symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
        #        bg = 1:11, inches = FALSE, add = TRUE)
        dev.off()
        pngfile<-pngfile_big
        pngfile_big<-image_border(pngfile_big, bg, "30x30")
        pngfile_big<-image_modulate(pngfile_big, brightness = img_brightness + 25 * length(candidateunique), saturation = 100)
        pngfile_big<-image_threshold(pngfile_big,type = "white", threshold = "50%",channel = "All")
        #pngfile_big<-image_annotate(pngfile_big,paste(cluster_desc),gravity = "north",size = 50,color = "white")
        pngfile_big<-image_annotate(pngfile_big,cluster_desc,gravity = "north",size = 30*40/nchar(cluster_desc),color = "white")
        pngfile_big<-image_annotate(pngfile_big,"0%         Relative intensity        100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)

        #pngfile_big<-image_annotate(pngfile_big,"0%                                               100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
        #pngfile_big<-image_annotate(pngfile_big,"                Relative intensity                 ",location = "+52+160",gravity = "northeast",size = 30,color = "white",degree=270)

        #pngfile<-image_average(img_com)
        #pngfile<-image_threshold(pngfile, type = "black",  threshold = "50%")
        pngfile<-image_border(pngfile, bg, "30x30")
        pngfile<-image_modulate(pngfile, brightness = img_brightness + 50)
      }

        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "north",size = 30,color = "white")
        if (remove_cluster_from_grid){
          pngfile=image_append(c((pngcompfile_output)))
        }else{
          pngfile=image_append(c(pngfile,(pngcompfile_output)))
        }

        bg = paste0("grey29")
        #property_png<-image_attributes(pngfile)
        #width_height<-as.numeric(unlist(stringr::str_split(property_png[property_png$property=="png:IHDR.width,height",2],", ")))
        #pngfile<-c(image_blank(width_height[1], width_height[2], color = bg, pseudo_image = ""),pngfile)
        image_write(pngfile,outputpngsum,flatten =T)



    }




  }

  if(export_Header_table){

     suppressMessages(suppressWarnings(require(gridExtra)))
     suppressMessages(suppressWarnings(require(grid)))

    if(is.null(candidate$desc)){
      candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
    }else{
      candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
      cluster_desc<-unique(candidate$desc)[1]
    }
     candidate_unique_table_Score<-candidate %>% group_by(mz) %>% summarise(Score=max(Score))
     candidate_unique_table_Score$Score<-round(candidate_unique_table_Score$Score,digits = 2)
     candidate_unique_table<-merge(candidate_unique_table,candidate_unique_table_Score)
     suppressMessages(suppressWarnings(require(reshape2)))

    Header_table<-NULL
    Header_table$mz=candidateunique
    Header_table<-data.frame(Header_table)
    Header_table=base::merge(Header_table,candidate_unique_table,all.x=T,by="mz",all.y=F)
    Header_table=Header_table[order(as.character(Header_table$mz)),]
    Header_table<-as.data.frame(Header_table[,c(componentID_colname,"mz","formula","adduct","Score")])

    #Header_table$ID=candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique)]
    #componentnames=unique(Header_table[[componentID_colname]][Header_table$mz==as.numeric(candidateunique[i])])
    t_Header_table<-as.data.frame(t(Header_table))
    t_Header_table<-sapply(t_Header_table,as.character)
    t_Header_table<-as.data.frame(t_Header_table)
    names(t_Header_table)<-as.character(t_Header_table[1,])
    tt3 <- ttheme_minimal()
    header_table_op<-tableGrob(t_Header_table,theme = tt3,cols = NULL,rows=NULL)

    find_cell <- function(table, row, col, name="core-fg"){
      l <- table$layout
      which(l$t==(row) & l$l==(col) & l$name==name)
    }
    if (remove_cluster_from_grid){
      t_Header_table<-cbind(t_Header_table)

    }else{
      t_Header_table<-cbind(c("ClusterID",clusterID,str_sub(cluster_desc,end = regexpr(" ",cluster_desc)),rep("",nrow(t_Header_table)-3)),t_Header_table)
    }
    header_table_op<-tableGrob(t_Header_table,cols = NULL,rows=NULL)


    if (remove_cluster_from_grid){
    for (mzfeatures in 1:(length(candidateunique))){
      for (rows in 3:3){
        ind <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
        #ind2 <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
        header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures], col = levels(mycol)[mzfeatures], lwd=5)
        #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
        #header_table_op$grobs[ind2][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)

        #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
      }}
    }else{
    for (mzfeatures in 2:(length(candidateunique)+1)){
    for (rows in 3:3){
    ind <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
    #ind2 <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
    header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = levels(mycol)[mzfeatures-1], lwd=5)
    #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
    #header_table_op$grobs[ind2][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)

    #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
    }
    }
    }

    suppressMessages(suppressWarnings(require(gtable)))

    header_table_op <- gtable_add_grob(header_table_op,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 2, b = nrow(header_table_op), l = 1, r = ncol(header_table_op))
    header_table_op <- gtable_add_grob(header_table_op,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(header_table_op))
    header_table_op$widths <- rep(max(header_table_op$widths), length(header_table_op$widths))

    header_file_png = windows_filename(paste0(clusterID,"_header.png"))
    if (remove_cluster_from_grid){
      png(header_file_png,width = 5*length(candidateunique),height = 5,units = "in",res = 300)

    }else{
      png(header_file_png,width = 5*length(candidateunique+1),height = 5,units = "in",res = 300)

    }

    grid.arrange(header_table_op,nrow=1)

    dev.off()

   if(F){ p <- plot_ly(
      type = 'table',
      columnwidth = 20,
      header = list(
        values = c(paste0("<b>","Cluster","</b>"),Header_table[,componentID_colname]),
        colspan = I(30),
        align = c('center'),
        line = list(width = 1, color = 'black'),
        fill = list(color = "grey"),
        font = list(family = "Arial", size = c(15,9+6/(width(Header_table[,componentID_colname])/7)), color = "white")
      ),
      cells = list(
        values = cbind(
          #rbind(paste0("<b>","Cluster","</b>"), as.matrix(Header_table$FA)),
          rbind("",as.matrix(Header_table$adduct)),
          rbind("",as.matrix(round(Header_table$mz,digits = 3)))
        ),
        align = c('left', rep('center', ncol(mtcars))),
        line = list(color = "black", width = 1),
        fill = list(color = "white"),
        font = list(family = "Arial", size = 15, color = c("black"))

      ),
      width=120*(nrow(Header_table)+10),
      height=400,


    ) %>% layout(
      title = list(text = paste0("<b>",clusterID,"</b>")),
      margin = list(t=0,b=0,l=0,r=0),
      font = list(family = "Arial",  color = "black",align = "bottom")
    )





    orca(p, file = windows_filename(paste0(clusterID,"_header.png")),width=120*(nrow(Header_table))+200,height=540) }

  if(file.exists(outputpngsum)){

      clusterpng<-image_read(outputpngsum)
      if (image_info(clusterpng)[2]>output_png_width_limit){
        clusterpng<-image_resize(clusterpng,paste0(output_png_width_limit,"x"))
      }
      header_file_pngfile<-header_file_png
      header_file_png<-image_read(header_file_png)
      header_file_png<-image_trim(header_file_png)
      header_file_png<-image_border(header_file_png, "white", "00x70")
      header_file_png<-image_resize(header_file_png,paste0(image_info(clusterpng)[2],"x"))
      clusterpng<-image_append(c(header_file_png,clusterpng),stack = T)

      image_write(clusterpng,outputpngsum)
      try(file.remove(header_file_pngfile))
    }
  }

  if(export_footer_table){
     suppressMessages(suppressWarnings(require(colorspace)))
     suppressMessages(suppressWarnings(require(stringr)))
     suppressMessages(suppressWarnings(require(ggplot2)))

    if(is.null(candidate$desc)){
      candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
      cluster_desc=ClusterID
    }else{
      candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
      cluster_desc<-unique(candidate$desc)[1]
    }

      prosequence<-list_of_protein_sequence[clusterID]
      candidate_unique_table=unique(candidate[,c(ClusterID_colname,componentID_colname,"Intensity","mz")])
      component_int<-candidate_unique_table %>% group_by_at((componentID_colname)) %>% dplyr::summarise(int=sum(Intensity))
      component_int$int<-component_int$int/max(component_int$int)

      s1=as.character(component_int$Peptide)
      s2=as.character(prosequence)

      palign2 <- sapply(s1,regexpr , s2)
      width_com<-str_length(s1)
      component_int$start=palign2
      component_int$end=component_int$start+width_com-1
      component_int<-merge(component_int,unique(candidate_unique_table[,c(componentID_colname,"mz")]),by=componentID_colname)
      component_int<-component_int[order(as.character(component_int$mz)),]

      component_int$mycol<-levels(mycol)



      transcolor<-rgb(0, 0, 0, max = 255, alpha = 0)
      pro_length<-unname(width(s2))
      pro_int<-rep(0,pro_length)
      pro_col<-rep(transcolor,pro_length)
      pro_col<-rep("grey93",pro_length)
      for(y in 1:nrow(component_int)){
        for( t in component_int$start[y]:component_int$end[y]){
          #pro_col[t]<-mixcolor(component_int$int[y]/(pro_int[t]+component_int$int[y]), col2RGB(pro_col[t]), col2RGB(component_int$mycol[y]))
          mixedcolor<-colorRamp(colors=c(pro_col[t],component_int$mycol[y]),  space ="rgb",
                    interpolate = "linear")((pro_int[t]+component_int$int[y])/(pro_int[t]+component_int$int[y]))
          pro_col[t]<-rgb(mixedcolor[,1],mixedcolor[,2],mixedcolor[,3], maxColorValue = 255)
          pro_int[t]<-pro_int[t]+component_int$int[y]
        }
      }

      ncharrow<-ceiling(width(s2)/length(candidateunique)/5)
      ncharw<-floor(width(s2)/ncharrow)
      component_int_plot<-data.frame(site=1:pro_length,int=pro_int,col=pro_col)
      component_int_plot$x<-(component_int_plot$site-1) %% ncharw
      component_int_plot$y<--((component_int_plot$site-1) %/% ncharw)
      component_int_plot$char<-str_sub(s2,start = component_int_plot$site,end = component_int_plot$site)


      if(F){p<-ggplot(component_int_plot, aes(x=site, y=int,fill=col)) + geom_area() + theme(legend.position="none") + ggtitle("Plot of length \n by dose") +
        theme(
          plot.title = element_text(color="red", size=14, face="bold.italic"),
          axis.title.x = element_text(color="blue", size=14, face="bold"),
          axis.title.y = element_text(color="#993333", size=14, face="bold")
        )}

      footerpng<-paste(workdir,"/",windows_filename(substr(clusterID, 1, 10)),"_footer.png",sep="")

      if (footer_style=="Protein"){
        png(footerpng,width = 5*length(candidateunique+1),height = 5*ceiling(ncharrow/10),units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),mar=c(1, 0, 0, 0))
      p <- ggplot(component_int_plot, aes(x, y, label = char)) +
        geom_label(fill=component_int_plot$col,family = "mono",size=20) +
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
              xlab(cluster_desc) + ylim(min(component_int_plot$y)-0.4,max(component_int_plot$y)+0.4)

      print(p)
      dev.off()
      } else if (footer_style=="Length"){

        png(footerpng,width = 5*length(candidateunique+1),height = 1,units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),mar=c(1, 0, 0, 0))
      component_int_plot$x=as.factor(1)
      component_int_plot$y=1
      p<- ggplot(data=component_int_plot, aes(x=site, y=1,group=site,fill=component_int_plot$col,col=component_int_plot$col)) +
        geom_bar(stat="identity",fill=component_int_plot$col,col=component_int_plot$col)+
        theme(axis.line=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +

        ylab(cluster_desc)
      print(p)
      dev.off()

      }



      if(file.exists(outputpngsum)){

        clusterpng<-image_read(outputpngsum)
        footerpngfile<-footerpng
        footerpng<-image_read(footerpng)
        footerpng<-image_border(footerpng, "white", "00x70")
        footerpng<-image_resize(footerpng,paste0(image_info(clusterpng)[2],"x"))
        clusterpng<-image_append(c(clusterpng,footerpng),stack = T)

        image_write(clusterpng,outputpngsum)
        try(file.remove(footerpngfile))
      }

      if(F){
              wrap_strings <- function(vector_of_strings,width){as.character(sapply(vector_of_strings,FUN=function(x){paste(strwrap(x,width=width), collapse="\n")}))}

      plot(component_int_plot$site,component_int_plot$int,col=component_int_plot$col)
      title(bquote(wrap_strings(s2,50)),col.main="black",cex.main=0.25,adj  = 0,line = -1)


      for (components in 1:nrow(component_int)){
        pep_start<-component_int$start[components]
        pep_end<-component_int$end[components]
        pep_body<-str_sub(s2,pep_start,pep_end)
        if(pep_start==1){
          title(bquote(.(pep_body)*phantom(.(str_sub(s2,pep_end+1,pro_length)))),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }else if(pep_end==pro_length){
          title(bquote(phantom(.(str_sub(s2,1,pep_start-1)))*.(pep_body)),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }else{
          title(bquote(phantom(.(str_sub(s2,1,pep_start-1)))*.(pep_body)*phantom(.(str_sub(s2,pep_end+1,pro_length)))),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }
      }

      dev.off()
      }

    }

  if(attach_summary_cluster){
    pngfile<-image_read(outputpngsum)
    #temp_component_png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

    #png(temp_component_png,width = 5,height = 5 , bg = bg,units = "in",res = 300)
    #print(pngfile_big)
    #dev.off()
    #pngfile_foot<-image_read(temp_component_png)image_scale(img, "100")

    #bg = paste0("grey15")
    #property_png<-image_attributes(pngfile_big)
    #width_height<-as.numeric(unlist(stringr::str_split(property_png[property_png$property=="png:IHDR.width,height",2],", ")))
    #pngfile_big<-image_flatten(c(image_blank(width_height[1], width_height[2], color = bg, pseudo_image = ""),pngfile_big))

    #image_write(pngfile,outputpngsum,flatten =T)

    pngfiletry=tryCatch(image_append(c(pngfile,pngfile_big),stack = T))
   #
    if (class(pngfiletry)=="magick-image") image_write(pngfiletry,outputpngsum)
    #///rm(temp_component_png)

  }

    removetempfile(tmp_dir,matchword=c(".png$","^magick-"))

    rm(pngcompfile)

    # TEMPORARY: Save complete componentimg.mono list for debugging
    # if (exists("componentimg.mono")) {
    #   rds_filename <- paste0(workdir, "/componentimg_mono_", windows_filename(substr(clusterID, 1, 20)), ".rds")
    #   saveRDS(componentimg.mono, file = rds_filename)
    #   cat("TEMP DEBUG: Saved complete componentimg.mono list to ", rds_filename, "\n")
    # }
    
    message("Cluster image rendering done:", clusterID ,cluster_desc)

  }else{

  message("Cluster image rendering Skipped:", clusterID ,cluster_desc)

  }


}

cluster_score_grid<-function(clusterID,
                             SMPLIST,
                             imdata,
                             ClusterID_colname="Protein",
                             componentID_colname="Peptide",
                             Component_plot_threshold=2,
                             Component_plot_coloure=c("mono","as.cluster"),
                             smooth.image="gaussian",
                             contrast.enhance = "suppression",
                             colorpallet="Set1",
                             plot_layout=c("line","grid"),
                             export_Header_table=F,
                             export_footer_table=F,
                             plot_style=c("fleximaging","ClusterOnly","rainbow"),
                             protein_coverage=F,
                             footer_style="Length",
                             output_png_width_limit=198000,
                             attach_summary_cluster=T,
                             cluster_color_scale=c("blackwhite","fleximaging"),
                             remove_cluster_from_grid=T,
                             img_brightness= 100,ppm=20,
                             list_of_protein_sequence,
                             workdir=getwd(),
                             pixel_size_um=50,
                             Score_thres=NULL){
  #complementary(color="red", plot = TRUE, bg = "white", labcol = NULL, cex = 0.8, title = TRUE)
  windows_filename<- function(stringX){
    stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
    stringX<-gsub("\"", "", stringX)
    stringX<-gsub("\\|", " ", stringX)
    return(stringX)

  }
   suppressMessages(suppressWarnings(require(magick)))
   suppressMessages(suppressWarnings(require(stringr)))
   suppressMessages(suppressWarnings(require(data.table)))
   #Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
   suppressMessages(suppressWarnings(require(grid)))
   suppressMessages(suppressWarnings(require(plotly)))
   suppressMessages(suppressWarnings(require(dplyr)))
   suppressMessages(suppressWarnings(if(!require(colortools)) {
     library(devtools)
     install_github('gastonstat/colortools')
   }))
   suppressMessages(suppressWarnings(require(data.table)))
   suppressMessages(suppressWarnings(require(Cardinal)))
  #rotate the image
  #imdata@pixelData@data<-rotatetmp
  outputpngsum=paste(workdir,"/",windows_filename(substr(clusterID, 1, 15)),"_score_imaging",'.png',sep="")

  #message(outputpngsum)
  SMPLIST=as.data.frame(SMPLIST)
  outputpng=paste(workdir,"/",windows_filename(substr(clusterID, 1, 10)),"_score_plot",'.png',sep="")
  #message(outputpng)
  candidate=SMPLIST[SMPLIST[[ClusterID_colname]]==clusterID,]

  if(!is.null(Score_thres)){candidate=candidate[candidate$Score>=Score_thres,]}

  #candidate=candidate[order(as.character())]
  if (componentID_colname %in% colnames(candidate)){
  candidate_u<- candidate %>% dplyr::group_by(mz) %>% dplyr::summarise(Peptide=Peptide[1], .groups = 'drop')

  candidate<-merge(data.frame(candidate_u),candidate,by=c("mz",componentID_colname),sort=F)
  }
  candidateunique=as.numeric(as.character(unique(candidate[,"mz"])))
  candidateunique=candidateunique[order(as.character(candidateunique))]

  #check if candidates' mz within the dataset range
  candidateunique<-candidateunique[which(between(candidateunique,range(mz(imdata))[1],range(mz(imdata))[2]))]
  candidate<-candidate[candidate$mz %in% candidateunique,]

  if(is.null(candidate$desc)){
    candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
    cluster_desc<-"-"
  }else{
    candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
    cluster_desc<-unique(candidate$desc)[1]
  }
  if (length(candidateunique)>4){
 
    mycol=setColors("red", length(candidateunique))#wheel("red", num = length(candidateunique),bg = "white")
  } else if (length(candidateunique)==4){
    #tmp_cols = setColors("red", 12)
    #tetrad_colors <- tmp_cols[c(1, 3, 7, 9)]
    mycol=c("#FF0000", "#FFFF00", "#00FFFF", "#0000FF")
    #mycol=c("#FF0000", "#00FF00", "#0000FF")
  } else if (length(candidateunique)==3){
    #mycol=splitComp("red")
    mycol=c("#FF0000", "#00FF00", "#0000FF")
  } else if (length(candidateunique)==2){
    #mycol=complementary("red")
    mycol=c("#FF0000", "#00FFFF")
  } else if (length(candidateunique)==1){
    #mycol=complementary("red")
    mycol=c("#FF0000")
  }

  mycol <- as.factor(as.character(mycol))
  mycol=mycol[order(mycol)]
  #mycolrgb<-hex2RGB(mycol)
  #meancoldf<-colSums(mycolrgb@coords)/max(colSums(mycolrgb@coords))
  #meancol<-RGB(R=meancoldf[1],G=meancoldf[2],B=meancoldf[3])
  #colorspace::hex(meancol)
  #print(length(candidateunique))
  #print(Component_plot_threshold)
  #print(length(candidateunique)>=Component_plot_threshold)
  if (length(candidateunique)>=Component_plot_threshold){

    if (is.null(imdata)){
      message("No imaging data")

    }else
      {
       suppressMessages(suppressWarnings(require(RColorBrewer)))
       suppressMessages(suppressWarnings(require(Cardinal)))
       suppressMessages(suppressWarnings(require(EBImage)))
      #library(colortools)

      if (plot_style=="ClusterOnly"){
        png(outputpng,width = 10,height = 10, bg = "black",units = "in",res = 300)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,1),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)
        Cardinal::image(imdata, mz=candidateunique,
              col=mycol,
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              tolerance=ppm, units="ppm")


        l<-function(x,y,z,t=x+y,f=t+z){
          paste(t,f)
        }
        l(1,2,3)

        dev.off()
        pngfile<-image_read(outputpng)
        pngfile<-image_border(pngfile, "black", "30x30")
        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        pngfile<-image_trim(pngfile)
        image_write(pngfile,outputpng)

      }
      else if (plot_style=="rainbow"){


        if (plot_layout=="line"){
          png(outputpngsum,width = 5*((length(candidateunique)+1)),height = 5, bg = "black",units = "in",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)


        }else{
          png(outputpngsum,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)
        }




        Cardinal::image(imdata, mz=candidateunique,
              col=levels(mycol),
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              tolerance=ppm, units="ppm",key=F)

        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          Cardinal::image(imdata, mz=candidateunique[i],
                contrast.enhance=contrast.enhance,
                smooth.image = smooth.image ,
                col.regions=col.regions,

                normalize.image="none",
                tolerance=ppm, units="ppm")
          componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          for (component in componentname){
            text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          }
        }
        dev.off()
        #pngfile<-image_read(outputpngsum)
        #pngfile<-image_border(pngfile, "black", "30x30")
        #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        #pngfile<-image_trim(pngfile)
        #image_write(pngfile,outputpngsum)

        

      }
      else if(plot_style=="fleximaging"){

        ##################################
        lightmode()
        tmp_dir <- tempdir()
        temp_cluster_png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        temp_component_png=list()
        temp_component_png.mono=list()
        if (plot_layout=="line"){
          #png(temp_cluster_png,width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = "black",units = "in",res = 300)
          temp_cluster_png<-tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
          png(temp_cluster_png,width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = "black",units = "px",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)


        }else{
          #png(temp_cluster_png,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 300)
          temp_cluster_png<-tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
          png(temp_cluster_png,width = 750*2,height = 750 * (ceiling((length(candidateunique)+1)/2)), bg = "black",units = "px",res = 150)

          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)
        }


        col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[1])



        dev.off()
        #if (length(candidateunique)==2){
        #  bg = "grey9"
        #}else if (length(candidateunique)==3){
        #  bg = "grey8"
        #}else if (length(candidateunique)==3){
        #
        #}
        bg = paste0("grey",ifelse((15-length(candidateunique))>1,(15-length(candidateunique)),2))
        #bg = "transparent"
        componentimg=list()
        componentimg.mono=list()
        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          if  (Component_plot_coloure=="mono"){
            bg = paste0("grey",29)
            col.regions.mono=intensity.colors_customize1(colset = 2)
            col.regions=create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          }else if(Component_plot_coloure=="as.cluster"){
            col.regions=create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
            #col.regions=create_gradient_colors_robust(100, start="#030303", end=levels(mycol)[i])
            col=levels(mycol)[i]
          }

          if(cluster_color_scale=="blackwhite"){
            col.regions.mono=create_gradient_colors_robust(100, start="black", end="white")
            col.regions=create_gradient_colors_robust(100, start="black", end="white")
          }
          if (Component_plot_coloure=="mono"){
            
            componentimg.mono[[i]]=Cardinal::image(imdata, mz=candidateunique[i],
                                    #contrast.enhance=contrast.enhance,
                                    contrast.enhance = "none",
                                    smooth.image = smooth.image,
                                    colorscale=col.regions.mono,
                                    normalize.image="linear",
                                    tolerance=ppm, units="ppm",
                                    key=F,
                                    xlab=NULL,
                                    layout=c( length(levels(Cardinal::run(imdata))),1),
                                    bg = bg)
            
            
            temp_component_png[[i]]=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
            
            componentimg.mono[[i]][["par"]][["ann"]]=F
            componentimg.mono[[i]][["par"]][["bty"]]="n"
            componentimg.mono[[i]][["par"]][["pty"]]="s"
            componentimg.mono[[i]][["par"]][["xaxt"]]="n"
            componentimg.mono[[i]][["par"]][["yaxt"]]="n"
            componentimg.mono[[i]][["par"]][["fg"]]="white"
            componentimg.mono[[i]][["par"]][["oma"]]=c(0, 0, 0, 0)
            #componentimg.mono[[i]][["par"]][["mar"]]=c(0, 0, 0, 1)
            if (!is.null(componentimg.mono[[i]][["facets"]]) && 
                length(componentimg.mono[[i]][["facets"]]) > 0 &&
                !is.null(componentimg.mono[[i]][["facets"]][[1]])) {
              attr(componentimg.mono[[i]][["facets"]][[1]],"strip")$strip=F
              attr(componentimg.mono[[i]][["facets"]][[1]],"colorkey")$colorkey=T
            }

            # componentimg.mono[[i]][["facets"]][[1]][[1]][["dpage"]]="a"
            # componentimg.mono[[i]][["facets"]][[1]][[1]][["facet"]][[".feature.groups"]]="a"
            # componentimg.mono[[i]][["fids"]][[".feature.groups"]]="a"
            # componentimg.mono[[i]][["dpages"]]="a"
            #png(temp_component_png.mono[[i]],width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = bg,units = "in",res = 300)
            png(temp_component_png[[i]],width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = bg,units = "px",res = 150)
            par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
                #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
                bty="n",pty="s",xaxt="n",
                yaxt="n",
                no.readonly = TRUE,ann=FALSE)

            try(tryCatch(print(componentimg.mono[[i]])),silent = T)
            dev.off()
          }else{
        componentimg[[i]]=Cardinal::image(imdata, mz=candidateunique[i],
                                      #contrast.enhance=contrast.enhance,
                                      contrast.enhance = "none",
                                      smooth.image = smooth.image,
                                      colorscale=col.regions,
                                      normalize.image="linear",
                                      tolerance=ppm, units="ppm",
                                      key=F,
                                      xlab=NULL,
                                      layout=c( length(levels(Cardinal::run(imdata))),1)
                                      )


        temp_component_png[[i]]=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

        #png(temp_component_png[[i]],width = 5,height = 5 * length( levels(Cardinal::run(imdata))), bg = bg,units = "in",res = 300)
        png(temp_component_png[[i]],width = 750,height = 750 * length( levels(Cardinal::run(imdata))), bg = bg,units = "px",res = 150)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(length( levels(Cardinal::run(imdata))),1),
            #par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)

        try(tryCatch(print(componentimg[[i]])),silent = T)

        dev.off()
          }




          #componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          #for (component in componentname){
          #  text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          #}
        }

        temp_component_cover=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

        #png(temp_component_cover,width = 5,height = 0.5, bg = bg,units = "in",res = 300)
        temp_component_cover<-tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")
        png(temp_component_cover,width = 750,height = 750, bg = bg,units = "px",res = 150)
        plot.new()
        text("")
        #temp_component_cover<-image_draw(temp_component_cover)
        dev.off()
        temp_component_cover<-image_read(temp_component_cover)
        temp_component_cover<-image_crop(temp_component_cover, "750x80+0")
        #temp_component_cover<-image_background(temp_component_cover, "hotpink")
        makeacover<-(temp_component_cover)
        #print(makeacover)
        #pngfile<-(temp_cluster_png)
        #pngfile<-image_border(pngfile, bg, "30x30")
        #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "north",size = 50,color = "white")

        pngcompfile_org=list()
        pngcompfile=list()
        for (i in 1:length(candidateunique)){
          pngcompfile_org[[i]]<-image_read(temp_component_png[[i]])
          pngcompfile_org[[i]]<-image_flatten(c(pngcompfile_org[[i]],makeacover))
          #if (Component_plot_coloure=="mono"){
          #pngcompfile[[i]]<-(temp_component_png.mono[[i]] )
          #pngcompfile[[i]]<-image_flatten(c(pngcompfile[[i]],makeacover))
          #}else{
          pngcompfile[[i]]<-pngcompfile_org[[i]]
          #}

          pngcompfile[[i]]<-image_border(pngcompfile[[i]], bg, "30x30")
          #pngcompfile[[i]]<-image_annotate(pngcompfile[[i]],paste(unique(candidate[candidate$mz==candidateunique[i],"moleculeNames"]),candidateunique[i]),gravity = "north",size = 30,color = "white")

        }




        pngcompfile_output<-pngcompfile[[1]]
        if (length(pngcompfile)>=2){
        for (i in 2: length(pngcompfile)){
          pngcompfile_output<-c(pngcompfile_output,unlist(pngcompfile[[i]]))
        } }

        #img_com<-(temp_component_png[[1]])
        #img_com<-image_flatten(c(img_com,makeacover))
        #pngcompfile_org[[1]]<-NULL
        img_com<-pngcompfile_org[[1]]
        if (length(pngcompfile)>=2){
        for (i in 2: length(pngcompfile_org)){
          img_com<-c(img_com,unlist(pngcompfile_org[[i]]))
        }}
        #img_com<-rbind(unlist(pngcompfile_org))


        #lapply(channel_types(),function(x,pngfile){
        #  pngfile %>% image_threshold(type = "white", threshold = "50%",channel = x) %>% image_scale( "x600") %>% image_write(path = paste0(x,".png"))
        #},pngfile)
        cluster_desc<-gsub(stringr::str_extract(cluster_desc,"OS=.{1,}"),"",cluster_desc)

       if ((Component_plot_coloure=="mono")) {

         pngfile_big<-image_average(img_com)
         pngfile_big_info=magick::image_info(pngfile_big)

         pngfile_big_ratio<-c(diff(range(imdata@elementMetadata@listData[["x"]])),diff(range(imdata@elementMetadata@listData[["y"]])))/c(pngfile_big_info$width-1*150,pngfile_big_info$height-1*150)



         pngfile_big <- image_draw(pngfile_big)
         #rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
         pixel_size_um_bar<-200
         pixel_size<-pixel_size_um_bar*pixel_size_um/(1/max(pngfile_big_ratio))/1000
         if(!is.integer(pixel_size)){
           pixel_size<-ceiling(pixel_size)
           pixel_size_um_bar<-pixel_size*1000/pixel_size_um*(1/max(pngfile_big_ratio))
         }

         rect(690-pixel_size_um_bar-60, 680,690-pixel_size_um_bar-60+pixel_size_um_bar,  680+13, border  = 'white',col  = 'white', lwd = 0)
         text(690, 690, paste0(pixel_size ," mm"),col  = 'white', cex = 3, srt = 0)
         #palette(rainbow(11, end = 0.9))
         #symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
         #        bg = 1:11, inches = FALSE, add = TRUE)
         dev.off()


         pngfile<-pngfile_big
         pngfile_big<-image_border(pngfile_big, bg, "30x30")
         #pngfile_big<-image_modulate(pngfile_big, brightness = 100 + 25 * length(candidateunique), saturation = 100)
         #pngfile_big<-image_threshold(pngfile_big,type = "white", threshold = "50%",channel = "All")
         #pngfile_big<-image_annotate(pngfile_big,paste(cluster_desc),gravity = "north",size = 50,color = "white")
         pngfile_big<-image_annotate(pngfile_big,cluster_desc,gravity = "north",size = 30*40/nchar(cluster_desc),color = "white")

         pngfile_big<-image_annotate(pngfile_big,"0%        Relative intensity        100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
         #pngfile_big<-image_annotate(pngfile_big,"0                                               100",location = "+52+165",gravity = "northeast",size = 30,color = "white",degree=270)
         #pngfile_big<-image_annotate(pngfile_big,"                 Relative intensity                  ",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)

         #pngfile<-image_average(img_com)
         #pngfile<-image_threshold(pngfile, type = "black",  threshold = "50%")
         pngfile<-image_border(pngfile, bg, "30x30")
         #pngfile<-image_modulate(pngfile, brightness = 150)


         }
        else if ((Component_plot_coloure=="as.cluster")){
        pngfile_big<-image_average(img_com)
        pngfile_big_info=magick::image_info(pngfile_big)

        pngfile_big_ratio<-c(diff(range(imdata@elementMetadata@listData[["x"]])),diff(range(imdata@elementMetadata@listData[["y"]])))/c(pngfile_big_info$width-1*150,pngfile_big_info$height-1*150)



        pngfile_big <- image_draw(pngfile_big)
        #rect(20, 20, 200, 100, border = "red", lty = "dashed", lwd = 5)
        pixel_size_um_bar<-200
        pixel_size<-pixel_size_um_bar*pixel_size_um/(1/max(pngfile_big_ratio))/1000
        if(!is.integer(pixel_size)){
          pixel_size<-ceiling(pixel_size)
          pixel_size_um_bar<-pixel_size*1000/pixel_size_um*(1/max(pngfile_big_ratio))
        }

        rect(690-pixel_size_um_bar-60, 680,690-pixel_size_um_bar-60+pixel_size_um_bar,  680+13, border  = 'white',col  = 'white', lwd = 0)
        text(690, 690, paste0(pixel_size ," mm"),col  = 'white', cex = 3, srt = 0)
        #palette(rainbow(11, end = 0.9))
        #symbols(rep(200, 11), seq(0, 400, 40), circles = runif(11, 5, 35),
        #        bg = 1:11, inches = FALSE, add = TRUE)
        dev.off()
        pngfile<-pngfile_big
        pngfile_big<-image_border(pngfile_big, bg, "30x30")
        pngfile_big<-image_modulate(pngfile_big, brightness = img_brightness + 25 * length(candidateunique), saturation = 100)
        pngfile_big<-image_threshold(pngfile_big,type = "white", threshold = "50%",channel = "All")
        #pngfile_big<-image_annotate(pngfile_big,paste(cluster_desc),gravity = "north",size = 50,color = "white")
        pngfile_big<-image_annotate(pngfile_big,cluster_desc,gravity = "north",size = 30*40/nchar(cluster_desc),color = "white")
        pngfile_big<-image_annotate(pngfile_big,"0%         Relative intensity        100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)

        #pngfile_big<-image_annotate(pngfile_big,"0%                                               100%",location = "+55+160",gravity = "northeast",size = 30,color = "white",degree=270)
        #pngfile_big<-image_annotate(pngfile_big,"                Relative intensity                 ",location = "+52+160",gravity = "northeast",size = 30,color = "white",degree=270)

        #pngfile<-image_average(img_com)
        #pngfile<-image_threshold(pngfile, type = "black",  threshold = "50%")
        pngfile<-image_border(pngfile, bg, "30x30")
        pngfile<-image_modulate(pngfile, brightness = img_brightness + 50)
      }

        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "north",size = 30,color = "white")
        if (remove_cluster_from_grid){
          pngfile=image_append(c((pngcompfile_output)))
        }else{
          pngfile=image_append(c(pngfile,(pngcompfile_output)))
        }

        bg = paste0("grey29")
        #property_png<-image_attributes(pngfile)
        #width_height<-as.numeric(unlist(stringr::str_split(property_png[property_png$property=="png:IHDR.width,height",2],", ")))
        #pngfile<-c(image_blank(width_height[1], width_height[2], color = bg, pseudo_image = ""),pngfile)
        image_write(pngfile,outputpngsum,flatten =T)



    }




  }

  if(export_Header_table){

     suppressMessages(suppressWarnings(require(gridExtra)))
     suppressMessages(suppressWarnings(require(grid)))

    if(is.null(candidate$desc)){
      candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
    }else{
      candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
      cluster_desc<-unique(candidate$desc)[1]
    }
     candidate_unique_table_Score<-candidate %>% group_by(mz) %>% summarise(Score=max(Score))
     candidate_unique_table_Score$Score<-round(candidate_unique_table_Score$Score,digits = 2)
     candidate_unique_table<-merge(candidate_unique_table,candidate_unique_table_Score)
     suppressMessages(suppressWarnings(require(reshape2)))

    Header_table<-NULL
    Header_table$mz=candidateunique
    Header_table<-data.frame(Header_table)
    Header_table=base::merge(Header_table,candidate_unique_table,all.x=T,by="mz",all.y=F)
    Header_table=Header_table[order(as.character(Header_table$mz)),]
    Header_table<-as.data.frame(Header_table[,c(componentID_colname,"mz","formula","adduct","Score")])

    #Header_table$ID=candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique)]
    #componentnames=unique(Header_table[[componentID_colname]][Header_table$mz==as.numeric(candidateunique[i])])
    t_Header_table<-as.data.frame(t(Header_table))
    t_Header_table<-sapply(t_Header_table,as.character)
    t_Header_table<-as.data.frame(t_Header_table)
    names(t_Header_table)<-as.character(t_Header_table[1,])
    tt3 <- ttheme_minimal()
    header_table_op<-tableGrob(t_Header_table,theme = tt3,cols = NULL,rows=NULL)

    find_cell <- function(table, row, col, name="core-fg"){
      l <- table$layout
      which(l$t==(row) & l$l==(col) & l$name==name)
    }
    if (remove_cluster_from_grid){
      t_Header_table<-cbind(t_Header_table)

    }else{
      t_Header_table<-cbind(c("ClusterID",clusterID,str_sub(cluster_desc,end = regexpr(" ",cluster_desc)),rep("",nrow(t_Header_table)-3)),t_Header_table)
    }
    header_table_op<-tableGrob(t_Header_table,cols = NULL,rows=NULL)


    if (remove_cluster_from_grid){
    for (mzfeatures in 1:(length(candidateunique))){
      for (rows in 3:3){
        ind <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
        #ind2 <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
        header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures], col = levels(mycol)[mzfeatures], lwd=5)
        #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
        #header_table_op$grobs[ind2][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)

        #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
      }}
    }else{
    for (mzfeatures in 2:(length(candidateunique)+1)){
    for (rows in 3:3){
    ind <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
    #ind2 <- find_cell(header_table_op, rows, mzfeatures, "core-fg")
    header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = levels(mycol)[mzfeatures-1], lwd=5)
    #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)
    #header_table_op$grobs[ind2][[1]][["gp"]] <- gpar(fill=levels(mycol)[mzfeatures-1], col = "black", lwd=5)

    #header_table_op$grobs[ind][[1]][["gp"]] <- gpar(fill="darkolivegreen1", col = "darkolivegreen4", lwd=5)
    }
    }
    }

    suppressMessages(suppressWarnings(require(gtable)))

    header_table_op <- gtable_add_grob(header_table_op,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 2, b = nrow(header_table_op), l = 1, r = ncol(header_table_op))
    header_table_op <- gtable_add_grob(header_table_op,
                         grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                         t = 1, l = 1, r = ncol(header_table_op))
    header_table_op$widths <- rep(max(header_table_op$widths), length(header_table_op$widths))

    header_file_png = windows_filename(paste0(clusterID,"_header.png"))
    if (remove_cluster_from_grid){
      png(header_file_png,width = 5*length(candidateunique),height = 5,units = "in",res = 300)

    }else{
      png(header_file_png,width = 5*length(candidateunique+1),height = 5,units = "in",res = 300)

    }

    grid.arrange(header_table_op,nrow=1)

    dev.off()

   if(F){ p <- plot_ly(
      type = 'table',
      columnwidth = 20,
      header = list(
        values = c(paste0("<b>","Cluster","</b>"),Header_table[,componentID_colname]),
        colspan = I(30),
        align = c('center'),
        line = list(width = 1, color = 'black'),
        fill = list(color = "grey"),
        font = list(family = "Arial", size = c(15,9+6/(width(Header_table[,componentID_colname])/7)), color = "white")
      ),
      cells = list(
        values = cbind(
          #rbind(paste0("<b>","Cluster","</b>"), as.matrix(Header_table$FA)),
          rbind("",as.matrix(Header_table$adduct)),
          rbind("",as.matrix(round(Header_table$mz,digits = 3)))
        ),
        align = c('left', rep('center', ncol(mtcars))),
        line = list(color = "black", width = 1),
        fill = list(color = "white"),
        font = list(family = "Arial", size = 15, color = c("black"))

      ),
      width=120*(nrow(Header_table)+10),
      height=400,


    ) %>% layout(
      title = list(text = paste0("<b>",clusterID,"</b>")),
      margin = list(t=0,b=0,l=0,r=0),
      font = list(family = "Arial",  color = "black",align = "bottom")
    )





    orca(p, file = windows_filename(paste0(clusterID,"_header.png")),width=120*(nrow(Header_table))+200,height=540) }

  if(file.exists(outputpngsum)){

      clusterpng<-image_read(outputpngsum)
      if (image_info(clusterpng)[2]>output_png_width_limit){
        clusterpng<-image_resize(clusterpng,paste0(output_png_width_limit,"x"))
      }
      header_file_pngfile<-header_file_png
      header_file_png<-image_read(header_file_png)
      header_file_png<-image_trim(header_file_png)
      header_file_png<-image_border(header_file_png, "white", "00x70")
      header_file_png<-image_resize(header_file_png,paste0(image_info(clusterpng)[2],"x"))
      clusterpng<-image_append(c(header_file_png,clusterpng),stack = T)

      image_write(clusterpng,outputpngsum)
      try(file.remove(header_file_pngfile))
    }
  }

  if(export_footer_table){
     suppressMessages(suppressWarnings(require(colorspace)))
     suppressMessages(suppressWarnings(require(stringr)))
     suppressMessages(suppressWarnings(require(ggplot2)))

    if(is.null(candidate$desc)){
      candidate_unique_table=unique(candidate[,c("mz",componentID_colname,"formula","adduct")])
      cluster_desc=ClusterID
    }else{
      candidate_unique_table=unique(candidate[,c(componentID_colname,"mz","formula","adduct")])
      cluster_desc<-unique(candidate$desc)[1]
    }

      prosequence<-list_of_protein_sequence[clusterID]
      candidate_unique_table=unique(candidate[,c(ClusterID_colname,componentID_colname,"Intensity","mz")])
      component_int<-candidate_unique_table %>% group_by_at((componentID_colname)) %>% dplyr::summarise(int=sum(Intensity))
      component_int$int<-component_int$int/max(component_int$int)

      s1=as.character(component_int$Peptide)
      s2=as.character(prosequence)

      palign2 <- sapply(s1,regexpr , s2)
      width_com<-str_length(s1)
      component_int$start=palign2
      component_int$end=component_int$start+width_com-1
      component_int<-merge(component_int,unique(candidate_unique_table[,c(componentID_colname,"mz")]),by=componentID_colname)
      component_int<-component_int[order(as.character(component_int$mz)),]

      component_int$mycol<-levels(mycol)



      transcolor<-rgb(0, 0, 0, max = 255, alpha = 0)
      pro_length<-unname(width(s2))
      pro_int<-rep(0,pro_length)
      pro_col<-rep(transcolor,pro_length)
      pro_col<-rep("grey93",pro_length)
      for(y in 1:nrow(component_int)){
        for( t in component_int$start[y]:component_int$end[y]){
          #pro_col[t]<-mixcolor(component_int$int[y]/(pro_int[t]+component_int$int[y]), col2RGB(pro_col[t]), col2RGB(component_int$mycol[y]))
          mixedcolor<-colorRamp(colors=c(pro_col[t],component_int$mycol[y]),  space ="rgb",
                    interpolate = "linear")((pro_int[t]+component_int$int[y])/(pro_int[t]+component_int$int[y]))
          pro_col[t]<-rgb(mixedcolor[,1],mixedcolor[,2],mixedcolor[,3], maxColorValue = 255)
          pro_int[t]<-pro_int[t]+component_int$int[y]
        }
      }

      ncharrow<-ceiling(width(s2)/length(candidateunique)/5)
      ncharw<-floor(width(s2)/ncharrow)
      component_int_plot<-data.frame(site=1:pro_length,int=pro_int,col=pro_col)
      component_int_plot$x<-(component_int_plot$site-1) %% ncharw
      component_int_plot$y<--((component_int_plot$site-1) %/% ncharw)
      component_int_plot$char<-str_sub(s2,start = component_int_plot$site,end = component_int_plot$site)


      if(F){p<-ggplot(component_int_plot, aes(x=site, y=int,fill=col)) + geom_area() + theme(legend.position="none") + ggtitle("Plot of length \n by dose") +
        theme(
          plot.title = element_text(color="red", size=14, face="bold.italic"),
          axis.title.x = element_text(color="blue", size=14, face="bold"),
          axis.title.y = element_text(color="#993333", size=14, face="bold")
        )}

      footerpng<-paste(workdir,"/",windows_filename(substr(clusterID, 1, 10)),"_footer.png",sep="")

      if (footer_style=="Protein"){
        png(footerpng,width = 5*length(candidateunique+1),height = 5*ceiling(ncharrow/10),units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),mar=c(1, 0, 0, 0))
      p <- ggplot(component_int_plot, aes(x, y, label = char)) +
        geom_label(fill=component_int_plot$col,family = "mono",size=20) +
        theme(axis.line=element_blank(),axis.text.x=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +
              xlab(cluster_desc) + ylim(min(component_int_plot$y)-0.4,max(component_int_plot$y)+0.4)

      print(p)
      dev.off()
      } else if (footer_style=="Length"){

        png(footerpng,width = 5*length(candidateunique+1),height = 1,units = "in",res = 300)
      par(oma=c(0, 0, 0, 0),mar=c(1, 0, 0, 0))
      component_int_plot$x=as.factor(1)
      component_int_plot$y=1
      p<- ggplot(data=component_int_plot, aes(x=site, y=1,group=site,fill=component_int_plot$col,col=component_int_plot$col)) +
        geom_bar(stat="identity",fill=component_int_plot$col,col=component_int_plot$col)+
        theme(axis.line=element_blank(),
              axis.text.y=element_blank(),axis.ticks=element_blank(),
              axis.title.y=element_blank(),legend.position="none",
              panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),plot.background=element_blank()) +

        ylab(cluster_desc)
      print(p)
      dev.off()

      }



      if(file.exists(outputpngsum)){

        clusterpng<-image_read(outputpngsum)
        footerpngfile<-footerpng
        footerpng<-image_read(footerpng)
        footerpng<-image_border(footerpng, "white", "00x70")
        footerpng<-image_resize(footerpng,paste0(image_info(clusterpng)[2],"x"))
        clusterpng<-image_append(c(clusterpng,footerpng),stack = T)

        image_write(clusterpng,outputpngsum)
        try(file.remove(footerpngfile))
      }

      if(F){
              wrap_strings <- function(vector_of_strings,width){as.character(sapply(vector_of_strings,FUN=function(x){paste(strwrap(x,width=width), collapse="\n")}))}

      plot(component_int_plot$site,component_int_plot$int,col=component_int_plot$col)
      title(bquote(wrap_strings(s2,50)),col.main="black",cex.main=0.25,adj  = 0,line = -1)


      for (components in 1:nrow(component_int)){
        pep_start<-component_int$start[components]
        pep_end<-component_int$end[components]
        pep_body<-str_sub(s2,pep_start,pep_end)
        if(pep_start==1){
          title(bquote(.(pep_body)*phantom(.(str_sub(s2,pep_end+1,pro_length)))),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }else if(pep_end==pro_length){
          title(bquote(phantom(.(str_sub(s2,1,pep_start-1)))*.(pep_body)),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }else{
          title(bquote(phantom(.(str_sub(s2,1,pep_start-1)))*.(pep_body)*phantom(.(str_sub(s2,pep_end+1,pro_length)))),col.main=component_int$mycol[components],cex.main=0.25,adj  = 0,line = -1)
        }
      }

      dev.off()
      }

    }

  if(attach_summary_cluster){
    pngfile<-image_read(outputpngsum)
    #temp_component_png=tempfile(pattern = "file", tmpdir = tempdir(), fileext = ".png")

    #png(temp_component_png,width = 5,height = 5 , bg = bg,units = "in",res = 300)
    #print(pngfile_big)
    #dev.off()
    #pngfile_foot<-image_read(temp_component_png)image_scale(img, "100")

    #bg = paste0("grey15")
    #property_png<-image_attributes(pngfile_big)
    #width_height<-as.numeric(unlist(stringr::str_split(property_png[property_png$property=="png:IHDR.width,height",2],", ")))
    #pngfile_big<-image_flatten(c(image_blank(width_height[1], width_height[2], color = bg, pseudo_image = ""),pngfile_big))

    #image_write(pngfile,outputpngsum,flatten =T)

    pngfiletry=tryCatch(image_append(c(pngfile,pngfile_big),stack = T))
   #
    if (class(pngfiletry)=="magick-image") image_write(pngfiletry,outputpngsum)
    #///rm(temp_component_png)

  }

    removetempfile(tmp_dir,matchword=c(".png$","^magick-"))

    rm(pngcompfile)

    # TEMPORARY: Save complete componentimg.mono list for debugging
    # if (exists("componentimg.mono")) {
    #   rds_filename <- paste0(workdir, "/componentimg_mono_", windows_filename(substr(clusterID, 1, 20)), ".rds")
    #   saveRDS(componentimg.mono, file = rds_filename)
    #   cat("TEMP DEBUG: Saved complete componentimg.mono list to ", rds_filename, "\n")
    # }
    
    message("Cluster image rendering done:", clusterID ,cluster_desc)

  }else{

  message("Cluster image rendering Skipped:", clusterID ,cluster_desc)

  }

if (Export_protein_pixel_data){
  
}
  
}

cluster_image_cardinal_allinone<-function(clusterID,
                                          SMPLIST,
                                          ppm=20,
                                          imdata=Load_Cardinal_imaging(datafile[i],preprocessing = F,resolution = ppm),
                                          ClusterID_colname="Protein",
                                          componentID_colname="Peptide",
                                          Component_plot_threshold=2,
                                          smooth.image="gaussian",
                                          contrast.enhance = "suppression",
                                          colorpallet="Set1",
                                          plot_layout="grid",
                                          image_rotation_degree=NULL,
                                          export_Header_table=F,
                                          rotate_image=F,
                                          plot_style=c("fleximaging","ClusterOnly","rainbow")){
  #complementary(color="red", plot = TRUE, bg = "white", labcol = NULL, cex = 0.8, title = TRUE)
  windows_filename<- function(stringX){
    stringX<-stringr::str_remove_all(stringX,"[><*?:\\/\\\\]")
    stringX<-gsub("\"", "", stringX)
    stringX<-gsub("\\|", " ", stringX)
    return(stringX)

  }
  Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/ProgramData/Anaconda3/orca_app", sep = .Platform$path.sep))
  library(grid)
  library(plotly)
  #rotate the image
  #imdata@pixelData@data<-rotatetmp


  SMPLIST=as.data.frame(SMPLIST)
  outputpng=paste(getwd(),"/",windows_filename(substr(clusterID, 1, 10)),"_cluster_plot",'.png',sep="")
  #message(outputpng)
  candidate=SMPLIST[SMPLIST[[ClusterID_colname]]==clusterID,]
  #candidate=candidate[order(as.character())]
  candidateunique=as.numeric(as.character(unique(candidate[,"mz"])))
  if (length(candidateunique)>9){
    candidate.dt <- data.table(candidate)
    candidatet=candidate.dt[,list(Intensity=sum(Intensity)), by='mz']
    candidatet=candidatet[order(-candidatet$Intensity)]
    selections=as.numeric(t(candidatet[1:9,"mz"]))
    candidate=candidate[candidate$mz %in% selections,]
    candidateunique=as.numeric(unique(candidate[,"mz"]))
    candidateunique=candidateunique[order(as.character(candidateunique))]
    mycol <- factor(RColorBrewer::brewer.pal(length(candidateunique),colorpallet))
    mycol <- factor(mycol,levels(mycol))
  }else if (length(candidateunique)<3){
    candidateunique=candidateunique[order(as.character(candidateunique))]
    mycol <- factor(RColorBrewer::brewer.pal(3,colorpallet))
    mycol <- factor(mycol,levels(mycol))

  }else{
    candidateunique=candidateunique[order(as.character(candidateunique))]
    mycol <- factor(RColorBrewer::brewer.pal(length(candidateunique),colorpallet))
    mycol <- factor(mycol,levels(mycol))
  }



  if (length(candidateunique)>=Component_plot_threshold){

    if (is.null(imdata)){
      message("No imaging data")

    }else{



      library(RColorBrewer)
      library(Cardinal)
      library(EBImage)
      #library(colortools)

      if (plot_style=="ClusterOnly"){
        png(outputpng,width = 10,height = 10, bg = "black",units = "in",res = 300)
        par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,1),
            bty="n",pty="s",xaxt="n",
            yaxt="n",
            no.readonly = TRUE,ann=FALSE)
        Cardinal::image(imdata, mz=candidateunique,
              col=mycol,
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              tolerance=ppm, units="ppm")


        l<-function(x,y,z,t=x+y,f=t+z){
          paste(t,f)
        }
        l(1,2,3)

        dev.off()
        pngfile<-image_read(outputpng)
        pngfile<-image_border(pngfile, "black", "30x30")
        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        pngfile<-image_trim(pngfile)
        image_write(pngfile,outputpng)

      }else if (plot_style=="rainbow"){



        if (plot_layout=="line"){
          png(outputpngsum,width = 5*((length(candidateunique)+1)),height = 5, bg = "black",units = "in",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)


        }else{
          png(outputpngsum,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 150)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)
        }




        Cardinal::image(imdata, mz=candidateunique,
              col=levels(mycol),
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              tolerance=ppm, units="ppm")

        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          Cardinal::image(imdata, mz=candidateunique[i],
                contrast.enhance=contrast.enhance,
                smooth.image = smooth.image ,
                col.regions=col.regions,

                normalize.image="none",
                tolerance=ppm, units="ppm")
          componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          for (component in componentname){
            text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          }
        }
        dev.off()
        #pngfile<-image_read(outputpngsum)
        #pngfile<-image_border(pngfile, "black", "30x30")
        #pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        #pngfile<-image_trim(pngfile)
        #image_write(pngfile,outputpngsum)


      }else if(plot_style=="fleximaging"){

        ##################################



        if (plot_layout=="line"){
          png(outputpngsum,width = 5*((length(candidateunique)+1)),height = 5, bg = "black",units = "in",res = 300)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(1,(length(candidateunique)+1)),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)


        }else{
          png(outputpngsum,width = 5*2,height = 5*(ceiling((length(candidateunique)+1)/2)), bg = "black",units = "in",res = 300)
          par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(0, 0, 1, 1),mfrow = c(ceiling((length(candidateunique)+1)/2), 2),
              bty="n",pty="s",xaxt="n",
              yaxt="n",
              no.readonly = TRUE,ann=FALSE)
        }




        Cardinal::image(imdata, mz=candidateunique,
              col=levels(mycol),
              contrast.enhance = contrast.enhance,
              smooth.image = smooth.image ,
              superpose=TRUE,normalize.image="linear",
              tolerance=ppm, units="ppm")


        for (i in 1:length(candidateunique)){
          #image(imdata, mz=candidateunique[i], col=mycol[i], superpose=F,normalize.image="linear")
          col.regions <- create_gradient_colors_robust(100, start="black", end=levels(mycol)[i])
          Cardinal::image(imdata, mz=candidateunique[i],
                contrast.enhance=contrast.enhance,
                smooth.image = smooth.image,#smooth.image ,
                col.regions=intensity.colors_customize1(),
                normalize.image="none",
                tolerance=ppm, units="ppm",
                key=F,
                xlab=NULL,
                ylab=NULL,
          )
          componentname=unique(candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique[i])])
          for (component in componentname){
            text(cex=30/ifelse(nchar(component)>30,nchar(component),30),labels=paste0(component),x=1,y=1+(which(componentname==component)-1)*3,adj=0,pos=4,offset=1,col = "white")
          }
        }
        dev.off()
        pngfile<-image_read(outputpngsum)
        pngfile<-image_border(pngfile, "black", "30x30")
        pngfile<-image_annotate(pngfile,paste(clusterID),gravity = "south",size = 50,color = "white")
        pngfile<-image_trim(pngfile)
        image_write(pngfile,outputpngsum)
      }
    }




  }

  if(export_Header_table){
    library(gridExtra)
    library(grid)
    candidate_unique_table=unique(candidate[,c("mz",ClusterID_colname,"formula","moleculeNames" , "adduct")])
    Header_table<-NULL
    Header_table$mz=candidateunique
    Header_table<-data.frame(Header_table)
    Header_table=base::merge(Header_table,candidate_unique_table,all.x=T,by="mz",all.y=F)
    #Header_table$ID=candidate[[componentID_colname]][candidate$mz==as.numeric(candidateunique)]
    #componentnames=unique(Header_table[[componentID_colname]][Header_table$mz==as.numeric(candidateunique[i])])
    p <- plot_ly(
      type = 'table',
      columnwidth = 20,
      header = list(
        values = c(paste0("<b>","Cluster","</b>"),Header_table$FA),
        colspan = I(30),
        align = c('center'),
        line = list(width = 1, color = 'black'),
        fill = list(color = "grey"),
        font = list(family = "Arial", size = 15, color = "white")
      ),
      cells = list(
        values = cbind(
          #rbind(paste0("<b>","Cluster","</b>"), as.matrix(Header_table$FA)),
          rbind("",as.matrix(Header_table$adduct)),
          rbind("",as.matrix(round(Header_table$mz,digits = 3)))
        ),
        align = c('left', rep('center', ncol(mtcars))),
        line = list(color = "black", width = 1),
        fill = list(color = "white"),
        font = list(family = "Arial", size = 15, color = c("black"))

      ),
      width=120*(nrow(Header_table)+1),
      height=400,


    ) %>% layout(
      title = list(text = paste0("<b>",clusterID,"</b>")),
      margin = list(t=0,b=0,l=0,r=0),
      font = list(family = "Arial", size = 20, color = "black",align = "bottom")
    )





    orca(p, file = windows_filename(paste0(clusterID,"header.png")),width=120*(nrow(Header_table)+1),height=540)

  }




}

orca_initial<-function(){


Sys.setenv("PATH" = paste(paste(unique(str_split(Sys.getenv("PATH"),.Platform$path.sep)[[1]]), sep = .Platform$path.sep,collapse = .Platform$path.sep), "C:/Anaconda/orca_app", sep = .Platform$path.sep))


}



intensity.colors_customize1 <- function(n = 100, alpha = 1,colset=1) {
  col2 <- rainbow(3*n, alpha=alpha)[(2*n):1]
  if (colset==1){
  f <- colorRamp(colors=c("darkorchid4", "dodgerblue","green", "greenyellow","yellow"))
  } else if(colset==2){
  f <- colorRamp(colors=c("black", "dodgerblue","green", "greenyellow","yellow"))
  }
  else {
  f <- colorRamp(colors=c("darkorchid4", "blue", "darkseagreen1", "yellow"))
  }

  alpha <- col2rgb(col2, alpha=TRUE)[[4]]
  col1 <- sapply(seq(from=0, to=1, length.out=n), function(i) do.call(rgb,
                                                                      c(as.list(f(i)), maxColorValue=255, alpha=alpha)))
  col1
}

Scilslab_mass_interval<-function(masslist,ppm){

  writelines=list()
  writelines[[1]]="<ImagingResults flexImagingVersion=\"3.0\" last_modified=\"\" created=\"SCiLS export\">"
  line=2
  for (mz in rownames(masslist)){

    writelines[[line]]=paste0("<Result Type=\"PkFilter\" Name=\"",masslist[mz,"moleculeNames"],
                              "\" Color=\"#5555ff\" Show=\"0\" MinIntensity=\"0\" IntensityThreshold=\"100\" AbsIntens=\"0\" LogScale=\"0\" MinMass=\"",
                              masslist[mz,"mz"]-(masslist[mz,"mz"]*ppm/1000000),
                              "\" MaxMass=\"",
                              masslist[mz,"mz"]+(masslist[mz,"mz"]*ppm/1000000),
                              "\" Integrate=\"0\" FindMass=\"0\" RelMass=\"1\"></Result>")
    line=line+1
  }
  writelines[[line]]="</ImagingResults>"

  writeLines(unlist(writelines),"massinterval.MIR")

}

removetempfile<-function(temp_dir=tempdir(),matchword=c(".png$","^magick-"), intersect=F){
  filelist=list.files(temp_dir)
  if (intersect){
    test=grep(paste0(matchword,collapse = "+"),filelist,ignore.case = T)
  }else{
    test=grep(paste0(matchword,collapse = "|"),filelist,ignore.case = T)
  }
  file.remove(paste0(temp_dir,"/",filelist[test]))


}






plotRanges <- function(ranged,labels=NULL,do.labs=T,skip.plot.new=F,lty="solid", alt.y=NULL,
                       v.lines=FALSE,ylim=NULL,xlim=NULL,scl=c("b","Kb","Mb","Gb"),
                       col=NULL,srt=0,pos=4,pch=1,lwd=1,cex=1,...) {
  if(!is(ranged)[1] %in% c("RangedData","GRanges")) {
    warning("ranged needs to be a RangedData or GRanges object, plot likely to fail") ; return(NULL) }
  chk <- chrNums(ranged)
  typ <- is(ranged)[1]
  if(!is.null(alt.y)) {
    if(is.numeric(alt.y)) {
      if(length(alt.y)==1 | length(alt.y)==length(ranged)) {
        yy <- alt.y
      } else {
        warning("alt.y ignored, must be same length as ranged, or else length 1"); alt.y <- NULL
      }
    } else {
      if(is.character(alt.y)) {
        if(typ=="GRanges") {
          cn <- colnames(mcols(ranged)); df <- mcols(ranged)
        } else {
          cn <- colnames(ranged); df <- ranged
        }
        if(!alt.y %in% cn) { stop("alternative y.axis column name ",alt.y," not found in 'ranged'") }
        yy <- df[,alt.y]; rm(df)
      } else {
        warning("invalid value for alt.y, ignoring"); alt.y <- NULL
      }
    }
  }
  if(!is.null(labels)) {
    labels <- paste(labels)
    if(is.character(labels)) {
      if(length(labels)==1 | length(labels)==length(ranged)) {
        if(length(labels)==1) {
          if(typ=="GRanges") {
            cn <- colnames(mcols(ranged)); df <- mcols(ranged)
          } else {
            cn <- colnames(ranged); df <- ranged
          }
          if(!labels %in% cn) { stop("labels column name ",labels," not found in 'ranged'") }
          lab <- df[,labels]; rm(df)
        } else {
          lab <- labels
        }
      } else {
        warning("labels ignored, must be same length as ranged, or else length 1"); labels <- NULL
      }
    } else {
      warning("invalid value for labels, ignoring"); labels <- NULL
    }
  } else {
    lab <- rownames(ranged)
  }
  if(length(chk)>1) {
    warning(length(chk)," chromosomes in 'ranged', only using the first, chr",chk[1])
    ranged <- chrSel(ranged,1)
  }
  if(all(width(ranged)<=1)) { theyAreSnps <- TRUE } else { theyAreSnps <- FALSE }
  scl <- make.divisor(scl)
  xl <- range(c(start(ranged),end(ranged)),na.rm=T)
  xl <- xl + ((diff(xl)*0.1)*c(-1,1))
  xl <- xl/scl
  nr <- nrow(ranged); if(is.null(nr)) { nr <- length(ranged) }
  if(is.null(alt.y)) {
    yl <- c(0,(nr+2))
  } else {
    yl <- range(yy,na.rm=T)
  }
  if(is.numeric(ylim) & length(ylim)==2) {
    ylim <- range(ylim,na.rm=T)
    ydif <- diff(ylim)
    yl <- ylim
  }
  if(is.numeric(xlim) & length(xlim)==2) {
    xlim <- range(xlim,na.rm=T)
    xdif <- diff(xlim)
    xl <- xlim
  }
  if(is.null(alt.y)) {
    YY <- seq(from=yl[1],to=yl[2],length.out=nr+2)[-1]
  } else {
    if(length(yy)==1) { YY <- rep(yy,length(nr)) } else { YY <- yy }
  }
  #print(YY)
  if(!is.null(col)) {
    if(length(col)==1) {
      col <- rep(col,times=nr)
    } else {
      if(length(col)!=nr) { warning("col was not the same length as ranged, using first only"); col <- rep(col[1],nr) }
    }
  }
  if(is.null(col)) {
    if(nr>22) { colz <- rep("black",nr) } else { colz <- get.distinct.cols(nr) }
  } else { colz <- col[1:nr] }
  if(is.null(lab) & do.labs) { lab <- paste(1:nr) } # last resort
  if(!skip.plot.new) {
    position <- c(start(ranged[1,]),end(ranged[1,]))/scl
    Y <- YY[c(1,1)]
    #prv(position,Y)
    TY <- if(theyAreSnps) { "p" } else { "l" }
    if(v.lines) {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col="white", lty=lty, ...)
      abline(v=position,col=colz[1])
    } else {
      plot(x=position, y=Y, xlim=xl, ylim=yl, type=TY, col=colz[1], lty=lty, lwd=lwd, cex=cex, ...)
    }
    st <- 2
  } else {
    st <- 1
  }
  if(nr>1 | st==1) {
    for (cc in st:nr) {
      if(v.lines) {
        abline(v=c(start(ranged[cc,]),end(ranged[cc,]))/scl,col=colz[cc],lty=lty)
      } else {
        if(theyAreSnps) {
          points(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], pch=pch, cex=cex)
        } else {
          lines(x=c(start(ranged[cc,]),end(ranged[cc,]))/scl,y=YY[c(cc,cc)],col=colz[cc], lty=lty, lwd=lwd)
        }
      }
    }
  }
  if(do.labs) {
    for (cc in 1:nr) {
      if(v.lines) { YY <- rep(tail(YY,1),length(YY)) }
      V.scale <- (diff(head(YY,2))*0.5)
      if(length(V.scale)<1 | srt!=90) { V.scale <- 0 }
      text(x=start(ranged[cc,])/scl,y=YY[cc]+V.scale,labels=lab[cc],cex=0.6,pos=pos,offset=0,srt=srt)
    }
  }
}




col2RGB<-function(x){
   suppressMessages(suppressWarnings(require(colorspace)))
  x_RGB<-t(as.numeric(col2rgb(x)))
  return(RGB(x_RGB))
}

rgb2hex <- function(r,g,b) sprintf('#%s',paste(as.hexmode(c(r,g,b)),collapse = ''))


#' Cleavage_rules_fun
#'
#' This is a function will return a list of pre-defined enzyme digestion specificity.
#'
#' @return a table of Cleavage rules
#'
#' @examples
#' Cleavage_rules_fun()
#'
#' @export

Cleavage_rules_fun <- function(){
  c(
  ## Arg-C proteinase
  "arg-c proteinase"="R(?=\\w)",
  ## Asp-N endopeptidase
  "asp-n endopeptidase"="\\w(?=D)",
  ## BNPS-Skatole
  "bnps-skatole"="W(?=\\w)",
  ## Caspase 1
  "caspase1"="(?<=[FWYL]\\w[HAT])D(?=[^PEDQKR])",
  ## Caspase 2
  "caspase2"="(?<=DVA)D(?=[^PEDQKR])",
  ## Caspase 3
  "caspase3"="(?<=DMQ)D(?=[^PEDQKR])",
  ## Caspase 4
  "caspase4"="(?<=LEV)D(?=[^PEDQKR])",
  ## Caspase 5
  "caspase5"="(?<=[LW]EH)D(?=\\w)",
  ## Caspase 6
  "caspase6"="(?<=VE[HI])D(?=[^PEDQKR])",
  ## Caspase 7
  "caspase7"="(?<=DEV)D(?=[^PEDQKR])",
  ## Caspase 8
  "caspase8"="(?<=[IL]ET)D(?=[^PEDQKR])",
  ## Caspase 9
  "caspase9"="(?<=LEH)D(?=\\w)",
  ## Caspase 10
  "caspase10"="(?<=IEA)D(?=\\w)",
  ## Chymotrypsin - high specifity
  "chymotrypsin-high"="([FY](?=[^P]))|(W(?=[^MP]))",
  ## Chymotrypsin - low specifity
  "chymotrypsin-low"="([FLY](?=[^P]))|(W(?=[^MP]))|(M(?=[^PY]))|(H(?=[^DMPW]))",
  ## Clostripain
  "clostripain"="R(?=\\w)",
  ## CNBr
  "cnbr"="M(?=\\w)",
  ## Enterokinase
  "enterokinase"="(?<=[DE][DE][DE])K(?=\\w)",
  ## Factor Xa
  "factor xa"="(?<=[AFGILTVM][DE]G)R(?=\\w)",
  ## Formic acid
  "formic acid"="D(?=\\w)",
  ## Glutamyl endopeptidase
  "glutamyl endopeptidase"="E(?=\\w)",
  ## Granzyme B
  "granzyme-b"="(?<=IEP)D(?=\\w)",
  ## Hydroxylamine
  "hydroxylamine"="N(?=G)",
  ## Iodosobenzoic acid
  "iodosobenzoic acid"="W(?=\\w)",
  ## LysC
  "lysc"="K(?=\\w)",
  ## LysN
  "lysn"="\\w(?=K)",
  ## Neutrophil elastase
  "neutrophil elastase"="[AV](?=\\w)",
  ## NTCB (2-nitro-5-thiocyanobenzoic acid)
  "ntcb"="\\w(?=C)",
  ## Pepsin (pH 1.3)
  "pepsin1.3"="((?<=([^HKR][^P])|(^[^P]))[^R](?=[FL][^P]))|((?<=([^HKR][^P])|(^[^P]))[FL](?=\\w[^P]))",
  ## Pepsin (pH > 2.0)
  "pepsin"="((?<=([^HKR][^P])|(^[^P]))[^R](?=[FLWY][^P]))|((?<=([^HKR][^P])|(^[^P]))[FLWY](?=\\w[^P]))",
  ## Proline endopeptidase
  "proline endopeptidase"="(?<=[HKR])P(?=[^P])",
  ## Proteinase K
  "proteinase k"="[AEFILTVWY](?=\\w)",
  ## Staphylococcal Peptidase I
  "staphylococcal peptidase i"="(?<=[^E])E(?=\\w)",
  ## Thermolysin
  "thermolysin"="[^DE](?=[AFILMV])",
  ## Thrombin
  "thrombin"="((?<=G)R(?=G))|((?<=[AFGILTVM][AFGILTVW]P)R(?=[^DE][^DE]))",
  ## Trypsin
  "trypsin"="([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))")
}

#' grid.ftable
#'
#' This is a function will plot table.
#'
#' @return a table image
#'
#' @examples
#' grid.ftable()
#'
#' @export
grid.ftable <- function(d, padding = unit(4, "mm"), ...) {
  library(gridExtra)
  library(grid)
  nc <- ncol(d)
  nr <- nrow(d)

  ## character table with added row and column names
  extended_matrix <- cbind(c("", rownames(d)),
                           rbind(colnames(d),
                                 as.matrix(d)))

  ## string width and height
  w <- apply(extended_matrix, 2, strwidth, "inch")
  h <- apply(extended_matrix, 2, strheight, "inch")

  widths <- apply(w, 2, max)
  heights <- apply(h, 1, max)

  padding <- convertUnit(padding, unitTo = "in", valueOnly = TRUE)

  x <- cumsum(widths + padding) - 0.5 * padding
  y <- cumsum(heights + padding) - padding

  rg <- rectGrob(x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 width = unit(widths + padding, "in"),
                 height = unit(heights + padding, "in"))

  tg <- textGrob(c(t(extended_matrix)), x = unit(x - widths/2, "in"),
                 y = unit(1, "npc") - unit(rep(y, each = nc + 1), "in"),
                 just = "center")

  g <- gTree(children = gList(rg, tg), ...,
             x = x, y = y, widths = widths, heights = heights)

  grid.draw(g)
  invisible(g)
}

brewer.pal_n<-function(n,colorstyle){
  library(RColorBrewer)
  if (n>9){
    getPalette = colorRampPalette(brewer.pal(8, colorstyle))
    return(getPalette(n))
  }else{
    getPalette = brewer.pal(n, colorstyle)
    return(getPalette)

  }


}


Plot_Ion_image_Png_Cardinal<- function(WKdir, imagefile, png_filename, mz,adduct="M-H", Tolerance= 5, title="",smooth.image = "adaptive",Creat_new_file=TRUE,color="black",contrast.enhance = "suppression"){
  #x11()
  #dev.copy(png,paste(dir,"/",png_filename,'.png',sep=""))
  Tolerance<-round(Tolerance,digits = 5)
  pngfillewrite<-paste(WKdir,"/",png_filename,'.png',sep="")
  #pngfillewrite<-nice_file_create(pngfillewrite)
  tempfilename=tempfile(pattern = "file", tmpdir = tempdir())
  png(paste(tempfilename,".png",sep=""), bg = "transparent")
  par(oma=c(0, 0, 0, 0),tcl = NA,mar=c(2, 0, 0, 0),
      bty="n",pty="s",xaxt="n",
      yaxt="n",
      no.readonly = TRUE,ann=FALSE)
  try(Cardinal::image(imagefile ,mz=mz ,
                      tolerance=Tolerance, units="mz",
                      contrast.enhance = contrast.enhance,
                      smooth.image = smooth.image ,
                      col.regions=intensity.colors_customize(),
                      asp = 1,
                      add=F),
      silent = TRUE

  )
  dev.off()
  try(pngfile<-image_read(paste(tempfilename,".png",sep="")))
  res <- try(pngfile<-image_trim(pngfile),silent = TRUE)
  if (class(res) != "try-error"){
    if (title=="D-mode"){
      pngfile<-image_border(pngfile, "transparent", "0x30")
      pngfile<-image_annotate(pngfile,paste(png_filename,"D-mode","\n",adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 9)
      pngfile<-image_trim(pngfile)
      res2 <- try(pngfileoriginal<-image_read(pngfillewrite),silent = TRUE)
      if ('&'(class(res2) != "try-error",Creat_new_file==FALSE)){
        pngfile<-image_append(c(pngfileoriginal,pngfile))
        image_write(pngfile,pngfillewrite)}
      if ('&'(class(res2) != "try-error",Creat_new_file==TRUE)){
        image_write(pngfile,pngfillewrite)}
      if (class(res2) == "try-error"){
        image_write(pngfile,pngfillewrite)}
    }else{
      pngfile<-image_border(pngfile, "transparent", "30x30")
      pngfile<-image_annotate(pngfile,paste("\n",png_filename,"\n", adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 10,color = color)
      pngfile<-image_trim(pngfile)

      image_write(pngfile,pngfillewrite)
    }
  }
  file.remove(paste(tempfilename,".png",sep=""))
}


simple_ion_image<-function(workdir=WorkingDir(),
                           Image_Type="",
                           plotTolerance=5 ,
                           creat_new_file=TRUE,
                           Denoising=0,
                           color = "black",
                           interpolate =FALSE,
                           ...){

  datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
  listfile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.csv$")
  masslist<-NULL
  for (z in 1:length(datafile)){

    datafile[z]<-gsub(".imzML", "", datafile[z])
    options(warn=-1)
    MALDI_IMAGE <- quiet(importImzMl(paste(file.path(datafile[z]),".imzML",sep="")))
    options(warn=0)
    print(paste("Plot Ion Images from",datafile[z]))
    if (dir.exists(datafile[z])==FALSE){dir.create(datafile[z])}

    for (i in 1:length(listfile)){
      foldername<-str_remove(listfile[i],paste(workdir,"/",sep=""))
      foldername<-str_remove(foldername,".csv")
      outputDir<-paste(datafile[z],"/",foldername,sep="")
      if (dir.exists(outputDir)==FALSE){dir.create(outputDir)}
      masslist<-read.csv(listfile[i],header = TRUE)
      for (j in 1:length(masslist$moleculeNames)){
        if (is.na(masslist$mz[j])==FALSE){
          for (k in 1:length(str_split(masslist$moleculeNames[j],", ")[[1]])){
            Plot_Ion_image_Png(WKdir=outputDir,
                               imagefile=MALDI_IMAGE,
                               png_filename=windows_filename(str_split(masslist$moleculeNames[j],", ")[[1]][k]),
                               mz=masslist$mz[j],
                               adduct=masslist$adduct[j],
                               Tolerance= plotTolerance*masslist$mz[j]/1000000,
                               title="",Creat_new_file=TRUE,Neighbour = Denoising,
                               color = color,interpolate =interpolate )
          }}}
    }
  }
}

simple_ion_image_cardinal<-function(datafile=tk_choose.files(filter =  matrix(c( "imzml file", ".imzML",
                                                                                 "Text", ".txt", "All files", "*"),
                                                                              3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis"),
                                    #workdir=WorkingDir(),
                                    Image_Type="",
                                    plotTolerance=5 ,
                                    creat_new_file=TRUE,
                                    color = "black",
                                    smooth.image = "gaussian",
                                    contrast.enhance = "none",
                                    ...){
  library(stringr)
  library(magick)
  datafile<-gsub(".imzML", "", datafile)
  workdir=base::dirname(datafile[1])
  if (dir.exists(paste0(workdir ,"/Ion images/"))==FALSE){dir.create(paste0(workdir ,"/Ion images/"))}
  setwd(paste0(workdir ,"/Ion images"))
  listfile <- list.files(path=paste0(workdir ,"/Ion images"),full.names = TRUE,pattern = "\\.csv$")
  masslist<-NULL
  for (z in 1:length(datafile)){
    name <-gsub(paste0(base::dirname(datafile[z]),"/"),"",datafile[z])
    folder<-base::dirname(datafile[z])
    MALDI_IMAGE <- Cardinal::readImzML(name, folder )
    print(paste("Plot Ion Images from",name))
    if (dir.exists(paste0(workdir ,"/Ion images/",name))==FALSE){dir.create(paste0(workdir ,"/Ion images/",name))}

    for (i in 1:length(listfile)){
      listfilename<-gsub(base::dirname(listfile[i]),"",listfile[i])
      outputDir<-paste0(workdir ,"/Ion images/",name,"/",listfilename)
      if (dir.exists(outputDir)==FALSE){dir.create(outputDir)}
      masslist<-read.csv(listfile[i],header = TRUE)
      masslist<-masslist[!is.na(masslist$mz),]
      for (j in 1:length(masslist$moleculeNames)){
        if (is.na(masslist$mz[j])==FALSE){
          for (k in 1:length(str_split(masslist$moleculeNames[j],", ")[[1]])){
            Plot_Ion_image_Png_Cardinal(WKdir=outputDir,
                                        imagefile=MALDI_IMAGE,
                                        png_filename=windows_filename(str_split(masslist$moleculeNames[j],", ")[[1]][k]),
                                        mz=masslist$mz[j],
                                        adduct=masslist$adduct[j],
                                        Tolerance= plotTolerance*masslist$mz[j]/1000000,
                                        title="",
                                        Creat_new_file=TRUE,
                                        color = color,
                                        smooth.image = smooth.image,
                                        contrast.enhance = contrast.enhance)
          }}}
    }
  }
}


intensity.colors_customize <- function(n = 100, alpha = 1) {
  col2 <- rainbow(3*n, alpha=alpha)[(2*n):1]
  f <- colorRamp(colors=c("black", "blue", "green", "yellow", "red","#FF00FF","white"))
  alpha <- col2rgb(col2, alpha=TRUE)[[4]]
  col1 <- sapply(seq(from=0, to=1, length.out=n), function(i) do.call(rgb,
                                                                      c(as.list(f(i)), maxColorValue=255, alpha=alpha)))
  col1
}

Pathway_overview_graphite<-function(){
  p_load(graphite,graph )
  humanReactome <- graphite::pathways("hsapiens", "reactome")
  humanmeatbolome <- graphite::pathways("hsapiens", "smpdb")
  metab_url <-
    url("https://romualdi.bio.unipd.it/wp-uploads/2018/04/Terunuma_metabolite_expr.txt")
  mexpr <- read.table(metab_url, header = TRUE, sep = "\t", row.names = NULL,
                      stringsAsFactors = FALSE)

  p <- humanReactome[["Glycolysis"]]
  head(nodes(p))
  head(edges(p))
  head(nodes(p), which = "mixed")
  head(edges(p), which = "mixed")
  pathwayDatabases()
  g <- pathwayGraph(p)
  g <- pathwayGraph(p, which = "mixed")
  pSymbol <- convertIdentifiers(p, "SYMBOL")
  head(nodes(pSymbol))
  reactomeSymbol <- convertIdentifiers(humanReactome[1:5], "SYMBOL")
  cytoscapePlot(convertIdentifiers(reactome$`Unwinding of DNA`, "symbol"), which = "mixed")
  meta_g <- pathwayGraph(p, which = "metabolites")

}

Plot_Ion_image_Png<- function(WKdir, imagefile, png_filename, mz,adduct="M-H", Tolerance= 0.25, title="",Neighbour = 0,Creat_new_file=TRUE,color="black",interpolate =FALSE){
  #x11()


  #dev.copy(png,paste(dir,"/",png_filename,'.png',sep=""))
  Tolerance<-round(Tolerance,digits = 5)
  pngfillewrite<-paste(WKdir,"/",png_filename,'.png',sep="")
  #pngfillewrite<-nice_file_create(pngfillewrite)
  png(paste(WKdir,"/","temp.png",sep=""), bg = "transparent")
  try(plotMsiSlice(imagefile ,mz , tolerance=Tolerance, legend=FALSE,colRamp=colorRamp(c("black", "blue", "green", "yellow", "red","#FF00FF","white")),interpolate =interpolate ),silent = TRUE)
  dev.off()
  try(pngfile<-image_read(paste(WKdir,"/","temp.png",sep="")))
  res <- try(pngfile<-image_trim(pngfile),silent = TRUE)
  if (class(res) != "try-error"){

    kern=Kernel_convolution(Neighbour)
    if (title=="D-mode"){
      pngfile<-image_border(pngfile, "transparent", "0x30")
      pngfile<-image_annotate(pngfile,paste(png_filename,"D-mode","\n",adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 9)
      pngfile<-image_trim(pngfile)
      res2 <- try(pngfileoriginal<-image_read(pngfillewrite),silent = TRUE)
      if ('&'(class(res2) != "try-error",Creat_new_file==FALSE)){
        pngfile<-image_append(c(pngfileoriginal,pngfile))
        image_write(pngfile,pngfillewrite)}
      if ('&'(class(res2) != "try-error",Creat_new_file==TRUE)){
        image_write(pngfile,pngfillewrite)}
      if (class(res2) == "try-error"){
        image_write(pngfile,pngfillewrite)}
    }else{
      pngfile<-image_convolve(pngfile, kern)
      pngfile<-image_border(pngfile, "transparent", "20x20")
      pngfile<-image_annotate(pngfile,paste(png_filename,"\n", adduct,sprintf("%g",mz),"±",Tolerance),gravity = "south",size = 10,color = color)
      pngfile<-image_trim(pngfile)

      image_write(pngfile,pngfillewrite)
    }
  }
  file.remove(paste(WKdir,"/","temp.png",sep=""))
}

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
    layout(title = list(text = paste(plot_data$mz[1])),
           xaxis = list(title = "Region"),
           yaxis = list(title = "Relative Conc.",ticksuffix = "%"),
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

