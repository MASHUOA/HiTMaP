library(magrittr)
library(shiny)
library(shinyFiles)
library(future)
plan(multisession)

folder_file_DTDT<-"folder_file_DTDT"
    #titlePanel("HiTMaP"),
UIproject<-function(){    
    sidebarLayout(
        
    
    
    sidebarPanel(
       # img(src = "tm.bmp", height =72 , width = 260), br(),h5("MALDI-MSI Proteomics annotation"),br(),
       # shinyDirButton("Projectdir", "Chose or creat Project", "Creat and choose a Project folder"),
       # br(),h4("Project"),
       # verbatimTextOutput("Projectdir"),
       # br(),
       # verbatimTextOutput("imzmlibd"),
       # br(),
        fileInput("imzmlibd", "Choose imzml file(s)",
                  multiple = T,
                  accept = c(".imzml",
                             ".ibd")),
        
        fileInput("database", "Choose Fasta database(s)",
                  multiple = T,
                  accept = c(".fasta")),
        
        fileInput("config", "Choose configuration file(s)",
                  multiple = FALSE,
                  accept = c("text/csv",
                             "text/comma-separated-values,text/plain",
                             ".csv")),
        fileInput("allfiles", "Choose and upload all file(s)",
                  multiple = T),
        
        ),
    
        mainPanel(
        h5("project stats"),
        br(),
        #verbatimTextOutput("imzmlibd"),
        br(),
       # verbatimTextOutput("Projectfile"),
        br(),
        DT::dataTableOutput(folder_file_DTDT),
        br()
        ))
    }
   


