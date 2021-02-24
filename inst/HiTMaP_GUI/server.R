#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
source("projectdir/Serverproject.R")
source("Proteomics/ServerProteomics.R")
#source("Preprocessing/ServerPreprocessing.R")
if (!require(shiny)) install.packages("shiny")
library(shiny)
if (!require(promises)) install.packages("promises")
library(promises)
library(stringr)
if (!require(DT)) install.packages("DT")
library(DT)
if (!require(data.table)) install.packages("data.table")
library(data.table)
source("FutureTaskProcessor.R")
if (!require(future)) install.packages("future")
library(future)
library(magrittr)
library(shiny)
library(shinyFiles)
library(tidyverse)

DEBUG_CONSOLE_OUTPUT <- "debugConsoleOutput"
NUM_ASYNC_TASKS_RUNNING <- "NUM_ASYNC_TASKS_RUNNING"
ASYNC_DATA <- "ASYNC_DATA"
FAKE_DATA_INFORMATION <- "FAKE_DATA_INFORMATION"
folder_file_DTDT<-"folder_file_DTDT"



Pre_process_par<-c("imzml_file_Pre_processing","RDA_imzml_file_stats","workdir","Segmentation_per_file_pre","tolerance" ,
                   "mzrange_pre", "Segmentation_pre" ,"Segmentation_def_pre" ,"Segmentation_ncomp_pre","peakAlign_pre",
                   "Segmentation_variance_coverage_pre","Smooth_range_pre","rotate_pre","force_Pre_process_pre",
                   "use_Pre_processRDS_pre","reduceBaseline_pre","smoothSignal_pre","peakpick_pre","tolerance_pre")

Proteomics_par<-c("Database_stats","Fastadatabase","Digestion_site","Digestion_site_2","adducts","Modifications_fix","Modifications_var","Load_candidatelist",
                  "imzml_file","RDA_imzml_file_stats","workdir","Segmentation_per_file","tolerance" ,"peakAlign",
                  "mzrange", "Segmentation" ,"Segmentation_def" ,"Segmentation_ncomp",
                  "Segmentation_variance_coverage","Smooth_range","rotate","force_Pre_process",
                  "use_Pre_processRDS","reduceBaseline","smoothSignal","peakpick","tolerance",
                  "threshold","missedCleavages","Substitute_AA_AA","Substitute_AA_mol","Substitute_AA_mol_wwater","Decoy_search","Decoy_mode","Decoy_adducts",
                  "use_previous_candidates","IMS_analysis","FDR_cutoff","peptide_ID_filter","plot_matching_score",
                  "Protein_feature_summary","Formula_with_water","output_candidatelist","Peptide_feature_summary","Region_feature_summary","parallel")


# Define server logic required to draw a histogram
#shiny::shinyServer(function(input, output, session) {
server<-function(input,output,session,WorkingDir_global){ 
    WorkingDir_global<-get("WorkingDir_global", envir = .GlobalEnv)
    
    setwd(WorkingDir_global)
    
    shinyDirChoose(input, 'Projectdir', roots =  c(Root=getwd()))
    
    
    Projectdir <- reactive(input$Projectdir)
    Current_states <- reactiveValues()
    observeEvent(Projectdir(),{
        if (typeof(Projectdir())=="list"){
    Current_states[[folder_file_DTDT]]<<-current_states_filesum(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/"))
        }

    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE)

    observeEvent(Current_states[[folder_file_DTDT]],
        {
            updateSelectizeInput(session, 'imzml_file_Pre_processing', choices = Current_states[[folder_file_DTDT]]$filename[Current_states[[folder_file_DTDT]]$filetype=="IMS data"])
            updateSelectizeInput(session, 'imzml_file', choices = Current_states[[folder_file_DTDT]]$filename[Current_states[[folder_file_DTDT]]$filetype=="IMS data"])
            updateSelectizeInput(session, 'Fastadatabase', choices = Current_states[[folder_file_DTDT]]$filename[Current_states[[folder_file_DTDT]]$filetype=="Database"],
                                 selected = 1)
            
            #updateSelectizeInput(session, 'RDA_imzml_file_stats', choices = Current_states[[folder_file_DTDT]]$filename[Current_states[[folder_file_DTDT]]$filetype %in% c("IMS data")], server = TRUE)
            
            updateSelectizeInput(session, 'ID_folders_stats', choices = Current_states[[folder_file_DTDT]]$filename[Current_states[[folder_file_DTDT]]$filetype %in% c("ID folder")], server = TRUE )
            
        },
        ignoreNULL = TRUE,
        ignoreInit = TRUE
    )
    
    imzmlibd <- reactive(input$imzmlibd)

    output$Projectdir <- renderText({
        if (typeof(Projectdir())=="list"){
            #Projectdirglob
            #paste0(unlist(Projectdir()[[2]]),paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),collapse = '/',sep="")
            paste0("Root/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),sep="/")
            #Projectdir()
        }else{
            "Root"
            #"Root/"
        }
    })
    
    
    
    rotate_pre <- reactive(input$rotate_pre)
    output$rotate_pre <- renderPrint({rotate_pre()})
    
    output$Projectfile <- renderPrint({
        if (typeof(Projectdir())=="list"){
            dir(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/"))            #paste0(unlist(Projectdir()[[2]]),paste0(Projectdir()[[1]],collapse = '/',sep=""),collapse = '/',sep="")

        }else{
            dir(getwd())
        }
    })

    toListen_projectfile <- reactive({
        list(input$imzmlibd,input$database,input$config,input$allfiles,input$Projectdir)
    })
    #
     output$folder_file_DTDT <- DT::renderDataTable({
        DT::datatable(Current_states[[folder_file_DTDT]])
     })


    # output$Projectfile1 <- renderPrint({
    #     (Current_states[[folder_file_DTDT]])
    # })


    observeEvent(toListen_projectfile(),{
        library(stringr)
        if (typeof(input$Projectdir)=="list"){

            current_states_filesum_df<-current_states_filesum(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/"))
           }else{
            current_states_filesum_df<-current_states_filesum(getwd())
           }
            Current_states[[folder_file_DTDT]]<<-current_states_filesum_df

            #Current_states[["imzmlfiles_in_proj"]]<<-current_states_filesum_df$filename[current_states_filesum_df$filetype=="IMS data"]

            
    },
    ignoreNULL = TRUE,
    ignoreInit = TRUE
    )

    

    output$imzmlibd <- renderPrint(imzmlibd())
    
    toListen_upload <- reactive({
        list(input$imzmlibd,input$database,input$config,input$allfiles)
    })

    output$uploadfile <- renderPrint({
        if (typeof(toListen_upload())=="list"){
            toListen_upload()          #paste0(unlist(Projectdir()[[2]]),paste0(Projectdir()[[1]],collapse = '/',sep=""),collapse = '/',sep="")

        }else{
            ""
        }
    })


    observeEvent(toListen_upload(),{
         lapply(toListen_upload(), function(x,Projectdir){
             if (!is.null(x)){
                 file.copy((x$datapath), paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep="")))
                 file.rename(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/",basename(as.character(x$datapath))),paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/",(x$name)))
                 }
         },Projectdir=Projectdir())
},
    ignoreNULL = TRUE,
    ignoreInit = TRUE)

    seesion_task_no<-1
    
    observeEvent(input$Pre_processing_run,if(!is.null(req(input$imzml_file_Pre_processing)))
                 {
                     seesion_task_no <<- seesion_task_no + 1
                     workdir=paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/")
                     print(workdir)
                     if(!dir.exists(paste0(workdir,"/source/"))) dir.create(paste0(workdir,"/source/"))
                     tasktime<-format(Sys.time(), "%Y %b %d %X")
                     input_future=reactiveValuesToList(input)
                     
                     save(list=c("input_future",
                                 "workdir",
                                 "tasktime"
                     ),
                     file=paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"), ascii = T,compress  = FALSE
                     )
                     
                     R_cmd_file<-paste0(workdir,"/source/task_Pre_processing_",gsub(":| ","_",tasktime),"_scource.R")
                     input_future<-input_future[names(input_future) %in% Pre_process_par]
                     fileConn<-file(R_cmd_file)
                     
                     writeLines(
                         c("suppressMessages(suppressWarnings(require(HiTMaP)))",
                           paste0("#load(\"",paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"),"\")"),
                           dput_list(input_future),
                           paste0("workdir=","\"",workdir,"\""),
                           paste0("tasktime=","\"",tasktime,"\""),
                           "HiTMaP:::Preprocessing_segmentation(datafile=input_future$imzml_file_Pre_processing,
                                                      workdir=workdir ,Bypass_Segmentation=T,
                                                      segmentation_num=input_future$Segmentation_per_file_pre,threshold=0.001,
                                                      ppm=input_future$tolerance,
                                                      mzrange=input_future$mzrange_pre,
                                                      Segmentation=input_future$Segmentation_pre,
                                                      Segmentation_def=input_future$Segmentation_def_pre,
                                                      Segmentation_ncomp=input_future$Segmentation_ncomp_pre,
                                                      Segmentation_variance_coverage=input_future$Segmentation_variance_coverage_pre,
                                                      Smooth_range=as.numeric(input_future$Smooth_range_pre),
                                                      rotate_file=input_future$rotate_pre,
                                                      BPPARAM=bpparam(),
                                                      preprocess=list(force_preprocess=input_future$force_Pre_process_pre,
                                                                      use_preprocessRDS=input_future$use_Pre_processRDS_pre,
                                                                      smoothSignal=list(method=input_future$smoothSignal_pre),
                                                                      reduceBaseline=list(method=input_future$reduceBaseline_pre),
                                                                      peakPick=list(method=input_future$peakpick_pre),
                                                                      peakAlign=list(tolerance=input_future$peakAlign_pre, units=\"ppm\")))"),
                                fileConn)
                     # Preprocessing_segmentation(datafile=input_future$imzml_file_Pre_processing,
                     #                            workdir=dir(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/")) ,
                     #                            segmentation_num=input_future$Segmentation_per_file_pre,threshold=0.001,
                     #                            ppm=input_future$tolerance,
                     #                            mzrange=input_future$mzrange_pre,
                     #                            Segmentation=input_future$Segmentation_pre,
                     #                            Segmentation_def=input_future$Segmentation_def_pre,
                     #                            Segmentation_ncomp=input_future$Segmentation_ncomp_pre,
                     #                            Segmentation_variance_coverage=input_future$Segmentation_variance_coverage_pre,
                     #                            Smooth_range=input_future$Smooth_range_pre,
                     #                            rotate_file=input_future$rotate_pre,
                     #                            BPPARAM=bpparam(),
                     #                            preprocess=list(force_preprocess=input_future$force_Pre_process_pre,
                     #                                            use_preprocessRDS=input_future$use_Pre_processRDS_pre,
                     #                                            smoothSignal=list(method=input_future$smoothSignal_pre),
                     #                                            reduceBaseline=list(method=input_future$reduceBaseline_pre),
                     #                                            peakPick=list(method=input_future$peakpick_pre),
                     #                                            peakAlign=list(tolerance=input_future$tolerance_pre, units="ppm")))
                     close(fileConn)
                     
                     
                     
                     # Decoy_adducts, Decoy_mode, Decoy_search, Digestion_site, Digestion_site_2, FDR, 
                     # Fastadatabase, IMS_analysis, Load_candidatelist, 
                     # Modifications_fix, Modifications_var, PCA_run, PMFsearch, Peptide_feature_summary, Pre_processing_run, 
                     # Projectdir, Protein_feature_summary, Proteomics_run, 
                     # RDA_imzml_file_stats, Segmentation, Segmentation_def, Segmentation_def_pre, 
                     # Segmentation_ncomp, Segmentation_ncomp_pre, 
                     # Segmentation_pre, Segmentation_variance_coverage, Segmentation_variance_coverage_pre, 
                     # Smooth_range, Smooth_range_pre, adducts, allfiles, config, creatorLoadbtn, 
                     # database, duration, force_Pre_process, 
                     # force_Pre_process_pre, imzml_file, imzml_file_Pre_processing, imzmlibd, maintabset, missedCleavages, 
                     # mz_tolerance_ppm, mzrange, mzrange_pre, normalize, normalize_pre, peakAlign, peakAlign_pre, 
                     # peptide_ID_filter, reduceBaseline, reduceBaseline_pre, segmentation_run, smoothSignal, 
                     # smoothSignal_pre, start_proc_future, threshold, use_Pre_processRDS, use_Pre_processRDS_pre 
                     # 
                     # Segmentation_ncomp_pre=Segmentation_ncomp_pre
                     # force_Pre_process=force_Pre_process
                     # Smooth_range_pre=Smooth_range_pre
                     # reduceBaseline_pre=reduceBaseline_pre
                     # duration <- input$duration # This variable needs to be created for use in future object. When using fakeDataProcessing(input$duration) an error occurs: 'Warning: Error in : Operation not allowed without an active reactive context.'
                     dataName <- paste0("Pre_processing", "_", seesion_task_no)
                     #fakeDataProcessing() is in FutureTaskProcessor.R

                     #(input_future)<-input_futur()
                     # input_future <- list()
                     # for(i in names(input)){
                     #     input_future[[i]] <- isolate(input[[i]])
                     # }
                     
                     system_fun<-function(x) { system(paste0("Rscript \"",x,"\"")) }
                     
                     startAsyncTask(
                         dataName,
                         future(system_fun(x=R_cmd_file)),
                         callback = function(asyncCallbackResults) {
                             # asyncTaskName = asyncTaskName,
                             # taskResult = taskResult,
                             # submitTime = submitTime,
                             # endTime = endTime,
                             # elapsedTime = elapsedTime,
                             # caughtError = caughtError,
                             # caughtWarning = caughtWarning
                             callback_fun(asyncCallbackResults=asyncCallbackResults)

        
                         },
                         tracklink=R_cmd_file
                     ) #end callback and call to startAsyncTask
                     otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <<-
                         getRunningTasksStatus()
                     asyncDataNumber <<- asyncDataNumber + 1
                     if (asyncDataNumber > 100) {
                         asyncDataNumber <<- 1
                     }
                     
                     library(stringr)
                     
                     if (typeof(input$Projectdir)=="list"){
                         
                         current_states_filesum_df<-current_states_filesum(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/"))
                     }else{
                         current_states_filesum_df<-current_states_filesum(getwd())
                     }
                     Current_states[[folder_file_DTDT]]<<-current_states_filesum_df
                     
                     #Current_states[["imzmlfiles_in_proj"]]<<-current_states_filesum_df$filename[current_states_filesum_df$filetype=="IMS data"]
                     # updateSelectizeInput(session, 'imzml_file_Pre_processing', choices = current_states_filesum_df$filename[current_states_filesum_df$filetype=="IMS data"])
                     # 
                     # updateSelectizeInput(session, 'RDA_imzml_file_stats', choices = current_states_filesum_df$filename[current_states_filesum_df$filetype %in% c("IMS data")], server = TRUE)
                     # 
                     # updateSelectizeInput(session, 'ID_folders_stats', choices = current_states_filesum_df$filename[current_states_filesum_df$filetype %in% c("ID folder")], server = TRUE )
                     # 
                 },
                 ignoreNULL = TRUE,
                 ignoreInit = TRUE)
    
    observeEvent(input$ID_folders_stats,
                 {
                     library(stringr)
                     
                     if (typeof(input$Projectdir)=="list"){
                         
                         workdir<-(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/"))
                         
                     }else{
                         
                         workdir<-(getwd())
                         
                     }
                     
                     stats_png_file <- dir(paste0(workdir,"/",input$ID_folders_stats))
                     
                     stats_png_file <- stats_png_file[str_detect(stats_png_file,regex(".png$|.jpg$",ignore_case = T))]
                     
                     updateSelectizeInput(session, 'Stats_img', choices = stats_png_file, server = TRUE)
                     
                 },
                 ignoreNULL = TRUE,
                 ignoreInit = TRUE)
    
    
    observeEvent(input$PCA_run,if(!is.null(req(input$imzml_file_Pre_processing)))
                 {
                     seesion_task_no <<- seesion_task_no + 1
                     workdir=paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/")
                     print(workdir)
                     if(!dir.exists(paste0(workdir,"/source/"))) dir.create(paste0(workdir,"/source/"))
                     tasktime<-format(Sys.time(), "%Y %b %d %X")
                     input_future=reactiveValuesToList(input)
                     
                     save(list=c("input_future",
                                 "workdir",
                                 "tasktime"
                     ),
                     file=paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"), ascii = T,compress  = FALSE
                     )
                     
                     R_cmd_file<-paste0(workdir,"/source/task_Pre_processing_segmentation_",gsub(":| ","_",tasktime),"_scource.R")
                     input_future<-input_future[names(input_future) %in% Pre_process_par]
                     fileConn<-file(R_cmd_file)
                     writeLines(
                         c("suppressMessages(suppressWarnings(require(HiTMaP)))",
                           paste0("#load(\"",paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"),"\")"),
                           dput_list(input_future),
                           paste0("workdir=","\"",workdir,"\""),
                           paste0("tasktime=","\"",tasktime,"\""),
                           "HiTMaP:::Preprocessing_segmentation(datafile=input_future$imzml_file_Pre_processing,
                                                      workdir=workdir ,Bypass_Segmentation=F,
                                                      segmentation_num=input_future$Segmentation_per_file_pre,threshold=0.001,
                                                      ppm=input_future$tolerance,
                                                      mzrange=input_future$mzrange_pre,
                                                      Segmentation=\"PCA\",
                                                      Segmentation_def=input_future$Segmentation_def_pre,
                                                      Segmentation_ncomp=\"auto-detect\",
                                                      Segmentation_variance_coverage=input_future$Segmentation_variance_coverage_pre,
                                                      Smooth_range=as.numeric(input_future$Smooth_range_pre),
                                                      rotate_file=input_future$rotate_pre,
                                                      BPPARAM=bpparam(),
                                                      preprocess=list(force_preprocess=input_future$force_Pre_process_pre,
                                                                      use_preprocessRDS=TRUE,
                                                                      smoothSignal=list(method=input_future$smoothSignal_pre),
                                                                      reduceBaseline=list(method=input_future$reduceBaseline_pre),
                                                                      peakPick=list(method=input_future$peakpick_pre),
                                                                      peakAlign=list(tolerance=input_future$peakAlign_pre, units=\"ppm\")))"),
                         fileConn)
                     close(fileConn)
                     dataName <- paste0("PCA_segmentation", "_", seesion_task_no)
                     system_fun<-function(x) { system(paste0("Rscript \"",x,"\"")) }
                     
                     startAsyncTask(
                         dataName,
                         future(system_fun(x=R_cmd_file)),
                         callback = function(asyncCallbackResults) {
                             # asyncTaskName = asyncTaskName,
                             # taskResult = taskResult,
                             # submitTime = submitTime,
                             # endTime = endTime,
                             # elapsedTime = elapsedTime,
                             # caughtError = caughtError,
                             # caughtWarning = caughtWarning
                             callback_fun(asyncCallbackResults=asyncCallbackResults)
                             
                             
                         },
                         tracklink=R_cmd_file
                     ) #end callback and call to startAsyncTask
                     otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <<-
                         getRunningTasksStatus()
                     asyncDataNumber <<- asyncDataNumber + 1
                     if (asyncDataNumber > 100) {
                         asyncDataNumber <<- 1
                     }
                     library(stringr)
                     
                     if (typeof(input$Projectdir)=="list"){
                         
                         current_states_filesum_df<-current_states_filesum(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/"))
                     
                         }else{
                          
                         current_states_filesum_df<-current_states_filesum(getwd())
                         }
                     
                     Current_states[[folder_file_DTDT]]<<-current_states_filesum_df
                 },
                 ignoreNULL = TRUE,
                 ignoreInit = TRUE)
    
 
    observeEvent(input$segmentation_run,if(!is.null(req(input$imzml_file_Pre_processing)))
                 {
                     seesion_task_no <<- seesion_task_no + 1
                     workdir=paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/")
                     print(workdir)
                     if(!dir.exists(paste0(workdir,"/source/"))) dir.create(paste0(workdir,"/source/"))
                     tasktime<-format(Sys.time(), "%Y %b %d %X")
                     input_future=reactiveValuesToList(input)
                     
                     save(list=c("input_future",
                                 "workdir",
                                 "tasktime"
                     ),
                     file=paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"), ascii = T,compress  = FALSE
                     )
                     
                     
                     input_future<-input_future[names(input_future) %in% Pre_process_par]
                         
                     R_cmd_file<-paste0(workdir,"/source/task_segmentation_",gsub(":| ","_",tasktime),"_scource.R")
                     
                     fileConn<-file(R_cmd_file)
                     
                     
                     writeLines(
                         
                         c("suppressMessages(suppressWarnings(require(HiTMaP)))",
                           
                           paste0("#load(\"",paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"),"\")"),
                           dput_list(input_future),
                           paste0("workdir=","\"",workdir,"\""),
                           paste0("tasktime=","\"",tasktime,"\""),
                           "HiTMaP:::Preprocessing_segmentation(datafile=input_future$imzml_file_Pre_processing,
                                                      workdir=workdir ,Bypass_Segmentation=F,
                                                      segmentation_num=input_future$Segmentation_per_file_pre,threshold=0.001,
                                                      ppm=input_future$tolerance,
                                                      mzrange=input_future$mzrange_pre,
                                                      Segmentation=input_future$Segmentation_pre,
                                                      Segmentation_def=input_future$Segmentation_def_pre,
                                                      Segmentation_ncomp=input_future$Segmentation_ncomp_pre,
                                                      Segmentation_variance_coverage=input_future$Segmentation_variance_coverage_pre,
                                                      Smooth_range=as.numeric(input_future$Smooth_range_pre),
                                                      rotate_file=input_future$rotate_pre,
                                                      BPPARAM=bpparam(),
                                                      preprocess=list(force_preprocess=input_future$force_Pre_process_pre,
                                                                      use_preprocessRDS=input_future$use_Pre_processRDS_pre,
                                                                      smoothSignal=list(method=input_future$smoothSignal_pre),
                                                                      reduceBaseline=list(method=input_future$reduceBaseline_pre),
                                                                      peakPick=list(method=input_future$peakpick_pre),
                                                                      peakAlign=list(tolerance=input_future$peakAlign_pre, units=\"ppm\")))"),
                         fileConn)
                     
                     close(fileConn)
                     
                     dataName <- paste0(input_future$Segmentation_pre,"_segmentation", "_", seesion_task_no)
                     system_fun<-function(x) { system(paste0("Rscript \"",x,"\"")) }
                     
                     startAsyncTask(
                         dataName,
                         future(system_fun(x=R_cmd_file)),
                         callback = function(asyncCallbackResults) {
                         callback_fun(asyncCallbackResults=asyncCallbackResults)
                         },
                         tracklink=R_cmd_file
                     ) 
                     otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <<-
                         getRunningTasksStatus()
                     asyncDataNumber <<- asyncDataNumber + 1
                     if (asyncDataNumber > 100) {
                         asyncDataNumber <<- 1
                     }
                 },
                 ignoreNULL = TRUE,
                 ignoreInit = TRUE)
    
    observeEvent(input$Proteomics_run,
                 {   
                     task_type<-("Proteomics_run")
                     task_type<-str_split(as.character(task_type),regex("\\$"))[[3]]
                     seesion_task_no <<- seesion_task_no + 1
                     workdir=paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/")
                     print(workdir)
                     if(!dir.exists(paste0(workdir,"/source/"))) dir.create(paste0(workdir,"/source/"))
                     tasktime<-format(Sys.time(), "%Y %b %d %X")
                     input_future=reactiveValuesToList(input)
                     input_future<-input_future[names(input_future) %in% Proteomics_par]

                     #input_future$AA<-ifelse(input_future$Substitute_AA_AA=="",NULL,input_future$Substitute_AA_AA)
                     #input_future$AA_new_formula<-ifelse(input_future$Substitute_AA_mol=="",NULL,input_future$Substitute_AA_mol)
                     #input_future$Formula_with_water<-ifelse(input_future$Substitute_AA_mol_wwater=="",NULL,input_future$Substitute_AA_mol_wwater)
                     input_future$Modifications_fix<-parse_mod(input_future$Modifications_fix)
                     input_future$Modifications_var<-parse_mod(input_future$Modifications_var)
                     #message(names(input_future) )
                     save(list=c("input_future",
                                 "workdir",
                                 "tasktime"
                     ),
                     file=paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"), ascii = T,compress  = FALSE
                     )
                     
                     

                     
                     R_cmd_file<-paste0(workdir,"/source/task_",task_type,"_",gsub(":| ","_",tasktime),"_scource.R")
                     
                     fileConn<-file(R_cmd_file)
                     
                     
                     
                     writeLines(
                         
                         c("suppressMessages(suppressWarnings(require(HiTMaP)))",
                           
                           paste0("#load(\"",paste0(workdir,"/source/input_",gsub(":| ","_",tasktime),".RData"),"\")"),
                           dput_list(input_future),
                           paste0("workdir=","\"",workdir,"\""),
                           paste0("tasktime=","\"",tasktime,"\""),

                         
                         "imaging_identification(
                             datafile=input_future$imzml_file,
                             projectfolder=workdir,
                             threshold=input_future$threshold, 
                             ppm=input_future$tolerance,
                             mode=\"Proteomics\",
                             Digestion_site=c(input_future$Digestion_site,input_future$Digestion_site_2),
                             missedCleavages=input_future$missedCleavages,
                             Fastadatabase=input_future$Fastadatabase,
                             adducts=input_future$adducts,
                             Modifications=list(fixed=input_future$Modifications_fix[['1']],fixmod_position=input_future$Modifications_fix[['2']],
                                                variable=input_future$Modifications_var[['1']],varmod_position=input_future$Modifications_var[['2']]),
                             Substitute_AA=list(AA=input_future$AA,AA_new_formula=input_future$AA_new_formula,Formula_with_water=input_future$Formula_with_water),
                             Decoy_search=TRUE,
                             Decoy_mode = \"isotope\",
                             Decoy_adducts = input_future$Decoy_adducts,
                             mzrange=input_future$mzrange,
                             adjust_score = FALSE,
                             IMS_analysis=TRUE,
                             Load_candidatelist=T,
                             Database_stats=input_future$Database_stats,
                             Bypass_generate_spectrum=FALSE,
                             peptide_ID_filter=input_future$peptide_ID_filter,
                             Protein_feature_summary=input_future$Protein_feature_summary,
                             Peptide_feature_summary=input_future$Peptide_feature_summary,
                             plot_ion_image=FALSE,
                             parallel=future::availableCores() - 2,
                             spectra_segments_per_file=input_future$Segmentation_per_file,
                             Segmentation=input_future$Segmentation,
                             Segmentation_def=input_future$Segmentation_def,
                             Segmentation_ncomp=input_future$Segmentation_ncomp,
                             Segmentation_variance_coverage=input_future$Segmentation_variance_coverage,
                             preprocess=list(force_preprocess=input_future$force_Pre_process,
                                             use_preprocessRDS=input_future$use_Pre_processRDS,
                                             smoothSignal=list(method=input_future$smoothSignal),
                                             reduceBaseline=list(method=input_future$reduceBaseline),
                                             peakPick=list(method=input_future$peakpick),
                                             peakAlign=list(tolerance=input_future$peakAlign, units=\"ppm\")),
                             Smooth_range=input_future$Smooth_range,
                             Virtual_segmentation_rankfile=NULL,
                             Rotate_IMG=input_future$rotate,
                             Region_feature_summary=input_future$Region_feature_summary,
                             Spectrum_validate=TRUE,
                             output_candidatelist=input_future$output_candidatelist,
                             use_previous_candidates=input_future$use_previous_candidates,
                             score_method=\"SQRTP\",
                             plot_cluster_image_grid=FALSE,
                             plot_cluster_image_maxretry=2,
                             plot_cluster_image_overwrite=F,
                             smooth.image=\"gaussian\",
                             componentID_colname=\"Peptide\",
                             ClusterID_colname=\"Protein\",
                             Protein_desc_of_interest=\".\",
                             Protein_desc_of_exclusion=NULL,
                             plot_unique_component=TRUE,
                             FDR_cutoff=0.05,
                             use_top_rank=NULL,
                             plot_matching_score=F,
                             Component_plot_coloure=\"mono\",
                             cluster_color_scale=\"blackwhite\",
                             plot_layout=\"line\",
                             export_Header_table=T,
                             export_footer_table=T,
                             attach_summary_cluster=T,
                             pixel_size_um=25,
                             img_brightness=100,
                             Thread=NULL,
                             cluster_rds_path=NULL,
                         )"),
                         fileConn)
                     
                     close(fileConn)
                     
                     dataName <- paste0(task_type, "_", seesion_task_no)
                     system_fun<-function(x) { system(paste0("Rscript \"",x,"\"")) }
                     
                     startAsyncTask(
                         dataName,
                         future(system_fun(x=R_cmd_file)),
                         callback = function(asyncCallbackResults) {
                             callback_fun(asyncCallbackResults=asyncCallbackResults)
                         },
                         tracklink=R_cmd_file
                     ) 
                     otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <<-
                         getRunningTasksStatus()
                     asyncDataNumber <<- asyncDataNumber + 1
                     if (asyncDataNumber > 100) {
                         asyncDataNumber <<- 1
                     }
                 },
                 ignoreNULL = TRUE,
                 ignoreInit = TRUE)
    
    Stats_img_outputs<-reactive(input$Stats_img)
    
    stats_img_zoom<-reactiveValues(zoom=1)
        
    output$Stats_img_output <- renderImage({

         list(src = paste0(paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/"),"/",input$ID_folders_stats,"/",input$Stats_img),
                 contentType = 'image/png',
                 width = 750 * stats_img_zoom$zoom,
                 alt = "This is alternate text")
        }, deleteFile = F)
    
    observeEvent(input$zoom_p, {
        if (stats_img_zoom$zoom < 1.5) stats_img_zoom$zoom <- stats_img_zoom$zoom + 0.07
    })
    observeEvent(input$zoom_m, {
        if (stats_img_zoom$zoom > 0.5) stats_img_zoom$zoom <- stats_img_zoom$zoom - 0.07
    })

    observeEvent(input$project_bookmark,{
        session$doBookmark()
        },
        ignoreNULL = TRUE,
        ignoreInit = TRUE)
    FAKE_PROCESSED_DATA <- "fakeProcessedData"
    
    otherReactiveValues <-
        reactiveValues() 
    #WARNING- DON'T USE VARIABLES TO INITIALIZE LIST KEYS - the variable name will be used, not the value
    
    asyncDataValues <- reactiveValues()
    
    asyncDataNumber <- 1
    
    otherReactiveValues[[DEBUG_CONSOLE_OUTPUT]] <-
        data.table::data.table(time = format(Sys.time(), "%Y %b %d %X"), message = "Placeholder to be deleted", tracklink = "...")[-1, ]
    
    otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <-
        getRunningTasksStatus()
    
    debugConsole <- function(msg,tracklink="...") {
        time <- format(Sys.time(), "%Y %b %d %X")
        newRow <- data.table::data.table(time = time, message = msg , tracklink=tracklink)
        otherReactiveValues[[DEBUG_CONSOLE_OUTPUT]] <<-
            rbind(newRow,
                  otherReactiveValues[[DEBUG_CONSOLE_OUTPUT]])
        print(paste0(nrow(otherReactiveValues[[DEBUG_CONSOLE_OUTPUT]]), ": ", time, ": ", msg, ": ",tracklink))
        flush.console()
    }
    
    observeEvent(input$Refresh_Running_task,{
        otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <<- getRunningTasksStatus()
    })
    
    callback_fun<-function(asyncCallbackResults){
         asyncTaskName <- asyncCallbackResults[["asyncTaskName"]]
         taskResult <- asyncCallbackResults[["taskResult"]]
         tracklink <- asyncCallbackResults[["tracklink"]]
         debugConsole(
             msg = paste(
                 #paste0("<a href='",  tracklink, "' target='_blank'>",asyncTaskName,"</a>"),
                 paste0("<a href='",  tracklink, "' target='_blank'>",asyncTaskName,"</a>"),
                 switch(typeof(taskResult),
                        "integer"="is Done.",
                        "character"=substr(taskResult,1,30),
                        "list"=paste("returning with data of size",object.size(taskResult))
                 
                 ),
                 paste0("Rscript file: ",  basename(tracklink))
                 
             ),
             tracklink=tracklink
             )
             asyncDataValues[[asyncTaskName]] <<- taskResult
             otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <<- getRunningTasksStatus()
     }
    #need to call processRunningTasks so that the callback to the futureFunction will be hit
    observe({
        invalidateLater(200)
        processRunningTasks(debug = TRUE)
    })
    
    output[[FAKE_DATA_INFORMATION]] <-
        renderText({
            paste(
                "Next item to be updated: ",
                paste0(FAKE_PROCESSED_DATA, "-", asyncDataNumber)
            )
        })
    
    try(output[[ASYNC_DATA]] <-
        renderText({
            myList <- reactiveValuesToList(asyncDataValues)
            if (length(myList) == 0) {
                result <- ""
            } else {
                result <-
                    sapply(1:length(myList), function(i) {
                        paste0("====> ", myList[[i]]$name, ": ", myList[[i]]$test, " ")
                    })
            }
            return(result)
        }))
    
    output[[NUM_ASYNC_TASKS_RUNNING]] <-renderDataTable({
         otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]]

        })
    
    observeEvent(input$start_proc_future,
                 {
                     duration <-
                         input$duration # This variable needs to be created for use in future object. When using fakeDataProcessing(input$duration) an error occurs: 'Warning: Error in : Operation not allowed without an active reactive context.'
                     dataName <-
                         paste0(FAKE_PROCESSED_DATA, "-", asyncDataNumber)
                     #fakeDataProcessing() is in FutureTaskProcessor.R
                     R_cmd_file<-"..."
                     startAsyncTask(
                         dataName,
                         future(fakeDataProcessing(dataName, duration)),
                         callback = function(asyncCallbackResults) {
                             callback_fun(asyncCallbackResults)

                         },
                         tracklink=R_cmd_file
                     ) #end callback and call to startAsyncTask
                     otherReactiveValues[[NUM_ASYNC_TASKS_RUNNING]] <<-
                         getRunningTasksStatus()
                     asyncDataNumber <<- asyncDataNumber + 1
                     if (asyncDataNumber > 100) {
                         asyncDataNumber <<- 1
                     }
                 },
                 ignoreNULL = TRUE,
                 ignoreInit = TRUE)
    
    output[[DEBUG_CONSOLE_OUTPUT]] = renderDataTable({
        DT::datatable(
        {
        data <- otherReactiveValues[[DEBUG_CONSOLE_OUTPUT]]

        if (nrow(data)>0){
            
        
        lapply(1:nrow(data), function(i,data) {
            output[[paste0("downloadData", i)]] <<- downloadHandler(
                filename = function() {
                    if (data$tracklink[i]!="...") {basename(data$tracklink[i])}else{"temp.R"}
                },
                content = function(file) {
                    RCMD<-readLines(data$tracklink[i])
                    #file.copy(data$tracklink[i], file)
                    #tar(tarfile=file, files=data$tracklink[i],compression="none")
                    # message(head(mtcars))
                    # head(RCMD)
                    writeLines(RCMD, con = file)
                }
            )
        },data)
        

        
        #data$tracklink <- NULL
        data %>%
            mutate('Source file' = lapply(1:n(),
                   function(i) paste0('<a href="#" onClick=document.getElementById("downloadData',i, '").click() >Download</a>')
            ))
        } else {
            data.table::data.table(time = format(Sys.time(), "%Y %b %d %X"), message = "Placeholder to be deleted", tracklink = "...")[-1, ]
        }
        },escape =2, options=list(columnDefs = list(list(visible=FALSE, targets=c(3)))))
        
    }, sanitize.text.function = function(x) x)
    
    
    output$hidden_downloads <- renderUI(

            lapply(1:nrow(otherReactiveValues[[DEBUG_CONSOLE_OUTPUT]]), function(i) {
                downloadLink(paste0("downloadData", i), "download", class = "hiddenLink")
            }
            )
        )

    output$projectSavebtn <- downloadHandler(

        filename = function() {
            workdir=paste0(getwd(),"/",paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/")
            rootname<-max(which(str_split(workdir,"/")[[1]] != ""))
            filename=paste(str_split(workdir,"/")[[1]][rootname],"_", Sys.Date(),".tar", sep="")
        
            return(filename[1])
            },
        content = function(filename) {
            
            workdir<-paste0(paste0(unlist(Projectdir()[[1]]),collapse = '/',sep=""),"/")
            #setwd(workdir)
            #message(workdir)
            files_df<-Current_states[[folder_file_DTDT]]
            files_dl<-files_df$filename[files_df$filetype %in% c("ID folder", "Summary")]
            files_dl<-paste0(workdir,files_dl,"/")
            files_dl<-gsub("^/","",files_dl)
            #message(files_dl)
            tar(tarfile=filename, files=files_dl,compression="none")
        },contentType = "application/tar"
    )
    
    observeEvent(input$reset, {
        shinyjs::reset("form")
    })
}


