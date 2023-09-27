#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
ui <- function(request) {
if (!require(shiny)) install.packages("shiny")
library(shiny)
if (!require(shinyjs)) install.packages("shinyjs")
library(shinyjs)
if (!require(future)) install.packages("future")
library(future)
library(shinythemes)
library(HiTMaP)

if (!require(data.table)) install.packages("data.table")
library(data.table)
library(magrittr)
library(shinyFiles)
#source("projectdir/UIproject.R")
#source("Pre-processing/UIPre-processing.R")
if (!require(DT)) install.packages("DT")
library(DT)
library(stringr)

DEBUG_CONSOLE_OUTPUT <- "debugConsoleOutput"
NUM_ASYNC_TASKS_RUNNING <- "NUM_ASYNC_TASKS_RUNNING"
ASYNC_DATA <- "ASYNC_DATA"
FAKE_DATA_INFORMATION <- "FAKE_DATA_INFORMATION"

folder_file_DTDT<-"folder_file_DTDT"
HiTMaP:::Peptide_modification(retrive_ID=NULL,update_unimod=F)
enzymes_option<-names(HiTMaP::Cleavage_rules_fun())
unimod.modification.df<-merge(unimod.df$modifications,unimod.df$specificity,by.x=c("record_id"),by.y=c("mod_key"),all.x=T)
unimod.modification.df=unimod.modification.df[unimod.modification.df$hidden==0,]
unimod.df$positions_new<-unimod.df$positions
unimod.df$positions_new$position_ext<-str_replace(unimod.df$positions_new$position,"Anywhere|Any N-term|Any C-term","Peptide")
unimod.df$positions_new$position_ext<-str_replace(unimod.df$positions_new$position_ext,"Protein N-term|Protein C-term","Protein")

unimod.modification.df<-merge(unimod.modification.df,unimod.df$positions_new,by.x="position_key",by.y="record_id")
unimod.modification.df<-unimod.modification.df[unimod.modification.df$ex_code_name!="",]
unimod.modification.df<-unimod.modification.df[order(as.numeric(unimod.modification.df$record_id)),]
mod_option<-paste0(unimod.modification.df$record_id,": ",unimod.modification.df$code_name," @ ",unimod.modification.df$position_ext,":",unimod.modification.df$position_key," @ ", unimod.modification.df$one_letter)

adducts_option<-HiTMaP:::Build_adduct_list()$Name

fluidPage(
  theme = shinythemes::shinytheme("lumen"),
  #shinythemes::themeSelector(),
    # tag to have textAreaInput of width 100%
    tags$style(HTML("
      .shiny-input-container:not(.shiny-input-container-inline) {
      width: 100%;}
      .hiddenLink {visibility: hidden;}"
                    )),
    useShinyjs()
    ,
    img(src = "tm.bmp",  width = 360),
  #h2("  MALDI-MSI Proteomics annotation"),
  br(),
  fluidRow(
    column(2,br(),
     shinyDirButton("Projectdir", "Choose or create Project", "Create or choose a Project folder",icon = icon("folder-open"), class = "btn-primary"),br()
      ),
    column(2,br(),
           downloadButton(outputId = "projectSavebtn", label = "Save project data")
    ),
    column(2,br(),
           actionButton("reset", "Reset setting")
    ),
    column(2,br(),
           bookmarkButton(label = "Bookmark setting",title = "Bookmark this project's setting and use the URL for later use.", id = "project_bookmark",icon = icon("save")),br()
    ),
    column(1,br(),
      h5("Project selected:"),br()),
    column(3,br(),h5(textOutput("Projectdir"))),br()),

    tabsetPanel(
        id = "maintabset", type = "tabs",
        tabPanel("Project", wellPanel(
            #UIproject()
            {
      sidebarLayout(
      sidebarPanel(
        # img(src = "tm.bmp", height =72 , width = 260), br(),h5("MALDI-MSI Proteomics annotation"),br(),
 #        shinyDirButton("Projectdir", "Chose or creat Project", "Creat and choose a Project folder"),
 #        br(),h4("Project"),
 #        verbatimTextOutput("Projectdir"),
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
                  accept = c("text/csv/xlsx",
                             "text/comma-separated-values,text/plain",
                             ".csv",".xlsx")),
        fileInput("allfiles", "Choose and upload all file(s)",
                  multiple = T)

      ),

      mainPanel(

        #downloadButton(outputId = "projectSavebtn", label = "Save project data"),

        #h5("Download ID/Summary folder"),br(),

        h3("Project stats"),

        DT::dataTableOutput(folder_file_DTDT)

      ))}
        ),icon = icon("briefcase")),
        tabPanel("Pre-processing", wellPanel(
            #UIPre-processing()
            {
              fluidRow(
                column(3,
                #img(src = "tm.bmp", height =72 , width = 260), br(),h5("MALDI-MSI Proteomics annotation"),br(),

                # verbatimTextOutput("imzmlibd"),
                # br(),
                selectizeInput(
                  'imzml_file_Pre_processing', 'Select IMS file(s)', choices = NULL,
                  multiple = TRUE, options = list(maxItems = 1000)
                ),br(),

                actionButton(inputId="Pre_processing_run", icon=NULL, label="Start Pre-proccessing (optional)"),
                br(),
                selectizeInput(
                  'rotate_pre', 'Select a rotation configuration file (optional)', choices = NULL,
                  multiple = F, options = list(maxItems = 1000)
                ),
                sliderInput("mzrange_pre", label = ("Set m/z Range"), min = 0,
                            max = 4000, value = c(500, 4000)),


                sliderInput("tolerance_pre", label = ("Set tolerance in ppm"), min = 0,
                            max = 50, value = 5),




                checkboxInput("force_Pre_process_pre", label = "Force pre-process data", value = FALSE),

                checkboxInput("use_Pre_processRDS_pre", label = "Use pre-processed RDS", value = TRUE),


                selectInput("normalize_pre", label = ("Normalize method"),
                            choices = list("TIC" = "tic", "RMS" = "rms", "Reference" = "reference", "Disable"="Disable"),
                            selected = "rms"),

                selectInput("smoothSignal_pre", label = ("Smooth Signal method"),
                            choices = list("Gaussian" = "gaussian", "Moving average" = "ma", "Savitzky-Golay" = "sgolay", "Disable"="Disable"),
                            selected = "Disable"),
                selectInput("peakpick_pre", label = ("Peak picking method"),
                            choices = list("Median absolute deviation" = "mad", "Simple" = "simple", "Adaptive" = "adaptive", "Disable"="Disable"),
                            selected = "adaptive"),
                selectInput("reduceBaseline_pre", label = ("Reduce Baseline method"),
                            choices = list("Local minimum" = "locmin", "Median" = "median", "Disable"="Disable"),
                            selected = "Disable"),


                sliderInput("peakAlign_pre", label = ("Set peak alignment tolerance (ppm)"), min = 1,
                            max = 100, value = 8)




              ),
              column(3,
              # selectizeInput(
              #     'RDA_imzml_file_stats', 'Select file(s) for statistical analysis', choices = NULL,
              #     multiple = T, options = list(maxItems = 1000)
              #   ),br(),
              actionButton(inputId="segmentation_run", icon=NULL, label="Start Image segmentation", class = "btn-primary"),
              br(),br(),
              actionButton(inputId="PCA_run", icon=NULL, label="Start PCA analysis (optional)"),
              br(),

              selectInput("Smooth_range_pre", label = ("Denoising level"),
                          choices = list("None" = 0, "Weak R1" = 1, "Medium R2" = 2, "Strong R3" = 3),
                          selected = 1),
              sliderInput("Segmentation_per_file_pre", label = ("Set initial number of segmentation"), min = 1,
                          max = 20, value = 6),
              sliderInput("Segmentation_variance_coverage_pre", label = ("Set variance coverage"), min = 0,
                          max = 1, value = 0.8),



              selectInput("Segmentation_pre", label = ("Segmentation method"),
                          choices = list("Spatial kMeans" = "spatialKMeans", "Spatial Shrunken Centroids" = "spatialShrunkenCentroids", "none" = "none","def_file"="def_file"),
                          selected = 1),

              selectInput("Segmentation_ncomp_pre", label = ("Segmentation nComp"),
                          choices = list("auto-detect" = "auto-detect", "3" = 3, "6" = 6, "9" = 9),
                          selected = 1),


              selectizeInput(
                'Segmentation_def_pre', 'Select one Segmentation label file', choices = NULL,
                multiple = F, options = list(maxItems = 1000)
               )




              ),
              column(6,
                h3("Segmentation Result"),
                br(),
                selectizeInput(
                  'ID_folders_stats', 'Select ID folder for result visualize', choices = NULL,
                  multiple = T, options = list(maxItems = 1)
                ),
                selectizeInput(
                  'Stats_img', 'Select Statistical result', choices = NULL,
                  multiple = T, options = list(maxItems = 1)
                ),
                br(),
                tagList(
                  actionButton(("zoom_m"), "", icon = icon("magnifying-glass-minus"), width = "40px",
                               style = "border-radius: 25px; padding: 0px;"),
                  actionButton(("zoom_p"), "", icon = icon("magnifying-glass-plus"), width = "40px",
                               style = "border-radius: 25px; padding: 0px;")
                ),
                imageOutput("Stats_img_output"),
                #uiOutput("Stats_img_output"),
                br()
              ))
              }
        ),icon = icon("think-peaks")),
         tabPanel("Proteomics annotation", wellPanel(
           {
             fluidRow(
               column(4,
                      actionButton(inputId="Proteomics_run", icon=NULL, label="Run Proteomics pipeline", class = "btn-primary"),br(),br(),
                      h3("Candidates and FDR"),
                      br(),
                      #verbatimTextOutput("imzmlibd"),
                      selectizeInput(
                        'Fastadatabase', 'Select database', choices = NULL,
                        multiple = TRUE, options = list(maxItems = 1)
                      ),

                      # selectInput("mode", label = ("Mode"),
                      #             choices = list("Proteomics" = "Proteomics", "Metabolomics" = "Metabolomics"),
                      #             selected = 1),


                      selectInput(
                        'Digestion_site', 'Enzyme(s)', choices = enzymes_option, selected= "trypsin",
                        multiple = TRUE
                      ),

                      textInput("Digestion_site_2", label = ("Optional digestion rule(s)"), value = ""),
                      sliderInput("missedCleavages", label = ("Missed Cleavages Range"), min = 0,
                                  max = 5, value = c(0, 1)),

                      selectInput(
                        'adducts', 'Adduct(s)', choices =adducts_option, selected= "M+H",
                         multiple=TRUE, selectize=TRUE
                      ),

                      selectizeInput(
                        'Modifications_fix', 'Fixed modification(s)', choices = mod_option,
                        multiple = TRUE, options = list(maxItems = 10)
                      ),

                      selectizeInput(
                        'Modifications_var', 'Variable modification(s)', choices = mod_option,
                        multiple = TRUE, options = list(maxItems = 10)
                      ),
                      textInput("Substitute_AA_AA", label = ("Substituted amino acid"), value = ""),
                      textInput("Substitute_AA_mol", label = ("Formula for Substitution"), value = ""),
                      checkboxInput("Substitute_AA_mol_wwater",label = "Formula without dehydration",value = T),
                      checkboxInput("use_previous_candidates", label = "Load previously generated candidate", value = TRUE),
                      checkboxInput("output_candidatelist", label = "Output generated candidate", value = TRUE),
                      checkboxInput("Database_stats", label = "Output generated candidate", value = FALSE),


                      br()
               ),
               column(4,

                      #img(src = "tm.bmp", height =72 , width = 260), br(),h5("MALDI-MSI Proteomics annotation"),br(),

                      # verbatimTextOutput("imzmlibd"),
                      # br(),
                      selectizeInput(
                        'imzml_file', 'Select multiple file for annotation', choices = NULL,
                        multiple = TRUE, options = list(maxItems = 1000)
                      ),

                     h3("Pre-processing"),br(),
                     selectizeInput(
                       'rotate', 'Select a rotation configuration file (optional)', choices = NULL,
                       multiple = F, options = list(maxItems = 1000)
                     ),
                     sliderInput("mzrange", label = ("Set m/z Range"), min = 0,
                                 max = 4000, value = c(500, 4000)),


                     sliderInput("tolerance", label = ("Set tolerance in ppm"), min = 0,
                                 max = 50, value = 5),




                     checkboxInput("force_Pre_process", label = "Force pre-process data", value = FALSE),

                     checkboxInput("use_Pre_processRDS", label = "Use pre-processed RDS", value = TRUE),


                     selectInput("normalize", label = ("Normalize method"),
                                 choices = list("TIC" = "tic", "RMS" = "rms", "Reference" = "reference"),
                                 selected = "rms"),

                     selectInput("smoothSignal", label = ("Smooth Signal method"),
                                 choices = list("Gaussian" = "gaussian", "Moving average" = "ma", "Savitzky-Golay" = "sgolay"),
                                 selected = "sgolay"),
                     selectInput("peakpick", label = ("Peak picking method"),
                                 choices = list("Median absolute deviation" = "mad", "Simple" = "simple", "Adaptive" = "adaptive"),
                                 selected = "adaptive"),
                     selectInput("reduceBaseline", label = ("Reduce Baseline method"),
                                 choices = list("Local minimum" = "locmin", "Median" = "median"),
                                 selected = "locmin"),


                     sliderInput("peakAlign", label = ("Set peak alignment tolerance (ppm)"), min = 1,
                                 max = 100, value = 8),


                      h3("Segmentation"),br(),

                      selectInput("Smooth_range", label = ("Denoising level"),
                                  choices = list("None" = 0, "Weak R1" = 1, "Medium R2" = 2, "Strong R3" = 3),
                                  selected = 1),
                      sliderInput("Segmentation_per_file", label = ("Set initial number of segmentation"), min = 1,
                                  max = 20, value = 6),
                      sliderInput("Segmentation_variance_coverage", label = ("Set variance coverage"), min = 0,
                                  max = 1, value = 0.8),



                      selectInput("Segmentation", label = ("Segmentation method"),
                                  choices = list("Spatial kMeans" = "spatialKMeans", "Spatial Shrunken Centroids" = "spatialShrunkenCentroids", "none" = "none","def_file"="def_file"),
                                  selected = "spatialKMeans"),

                      selectInput("Segmentation_ncomp", label = ("Segmentation nComp"),
                                  choices = list("auto-detect" = "auto-detect", "3" = 3, "6" = 6, "9" = 9),
                                  selected = "auto-detect"),


                      selectizeInput(
                        'Segmentation_def', 'Select one Segmentation label file', choices = NULL,
                        multiple = F, options = list(maxItems = 1000)
                      )

                     # "missedCleavages","Substitute_AA","Decoy_search","Decoy_mode","Decoy_adducts",
                     # "use_previous_candidates","IMS_analysis","FDR_cutoff","peptide_ID_filter","plot_matching_score",
                     # "Protein_feature_summary","Peptide_feature_summary","Region_feature_summary","parallel"
                     #
               ),
               column(4,
                      br(),h3("IMS analysis"),br(),
                      sliderInput("mz_tolerance_ppm", label = ("Precusor tolerance (ppm)"), min = 1,
                                  max = 100, value = 5),
                      sliderInput("threshold", label = ("Set relative signal threshold"), min = 0,
                                  max = 0.2, value = 0.005),
                      sliderInput("FDR_cutoff", label = ("False discovery rate cut-off"), min = 0,
                                  max = 1, value = 0.05),
                      sliderInput("peptide_ID_filter", label = ("Minimum peptide number for Prtoein ID"), min = 1,
                                  max = 10, value = 2),
                      checkboxInput("Decoy_search", label = "Decoy search", value = TRUE),

                      selectInput("Decoy_mode", label = ("Decoy Mode"),
                                  choices = list("Isotope" = "isotope", "Elements" = "elements", "Adducts"="adducts"),
                                  selected = 1),

                      selectizeInput("Decoy_adducts", label = ("Decoy adducts"),
                                     choices = c("M+ACN+H","M+IsoProp+H","M+DMSO+H","M+Co","M+Ag","M+Cu","M+He","M+Ne","M+Ar","M+Kr","M+Xe","M+Rn"),
                                     multiple = TRUE, options = list(maxItems = 10)),

                      checkboxInput("IMS_analysis", label = "IMS analysis", value = TRUE),
                      checkboxInput("PMFsearch", label = "PMF analysis", value = TRUE),

                      checkboxInput("Protein_feature_summary", label = "Protein feature summary", value = TRUE),
                      checkboxInput("Peptide_feature_summary", label = "Peptide feature summary list", value = TRUE),
                      br()
               ))

           }
         ),icon = icon("fingerprint")),
         tabPanel("Cluster image rendering", wellPanel(
           {
             sidebarLayout(
               sidebarPanel(
                 # img(src = "tm.bmp", height =72 , width = 260), br(),h5("MALDI-MSI Proteomics annotation"),br(),
                 #        shinyDirButton("Projectdir", "Chose or creat Project", "Creat and choose a Project folder"),
                 #        br(),h4("Project"),
                 #        verbatimTextOutput("Projectdir"),
                 # br(),
                 # verbatimTextOutput("imzmlibd"),
                 # br(),
                 selectizeInput(
                   'imzml_file_Pre_processing', 'Select result file for rendering', choices = NULL,
                   multiple = TRUE, options = list(maxItems = 1000)
                 ),

                 h3("Pre-processing"),br(),

                 sliderInput("mzrange", label = ("Set m/z Range"), min = 0,
                             max = 4000, value = c(500, 4000)),


                 selectInput("Smooth_range", label = ("Denoising level"),
                             choices = list("None" = 0, "Weak R1" = 1, "Medium R2" = 2, "Strong R3" = 3),
                             selected = 2),

                 checkboxInput("force_Pre_process", label = "Force pre-process data", value = TRUE),

                 checkboxInput("use_Pre_processRDS", label = "Use pre-processed RDS", value = TRUE),


                 selectInput("normalize", label = ("Normalize method"),
                             choices = list("tic" = "tic", "rms" = "rms", "reference" = "reference"),
                             selected = 2),

                 selectInput("smoothSignal", label = ("Smooth Signal method"),
                             choices = list("gaussian" = "gaussian", "ma" = "ma", "sgolay" = "sgolay"),
                             selected = 3),

                 selectInput("reduceBaseline", label = ("Reduce Baseline method"),
                             choices = list("locmin" = "locmin", "median" = "median"),
                             selected = 1),


                 sliderInput("peakAlign", label = ("Set peakAlign tolerance (ppm)"), min = 1,
                             max = 100, value = 8)



               ),

               mainPanel(
                 h3("Cluster image"),
                 verbatimTextOutput("input_print")

               ))}
         ),icon = icon("layer-group")),
        # tabPanel("Server", wellPanel(
        #     serverCreator()
        # )),
        # tabPanel("Global", wellPanel(
        #     globalCreator()
        # )),
        # tabPanel("Export result", wellPanel(
        #     downloadButton(outputId = "creatorSavebtn", label = "Save shinyWYSIWYG data"),
        #     fileInput(inputId = "creatorLoadbtn", label = "", accept = ".RData", buttonLabel = "Load shinyWYSIWYG data"),
        #     actionButton("input_print_button","input"),
        #     verbatimTextOutput("rotate_pre")
        #
        #     ),icon = icon("file-export
 
 tabPanel("Target selection", wellPanel(
   {
     fluidRow(
       column(3,
              
              selectizeInput(
                'Target_table_file', 'Select Target ID', choices = NULL,
                multiple = F, options = list(maxItems = 1000)
              ),br(),
              
              selectizeInput(
                'ID_table_file', 'Select ID input', choices = NULL,
                multiple = F, options = list(maxItems = 1000)
              ),br(),
              
              selectInput("Analysis_pipeline", label = ("Analysis pipeline"),
                          choices = list("ProteinPilot" = "ProteinPilot", "Fragpipe" = "Fragpipe", "UserTable" = "UserTable"),
                          selected = "ProteinPilot"),
              
              sliderInput("pep_length", label = ("Set peptide length range"), min = 5,
                          max = 50, value = c(7, 40)),
              
              sliderInput("Score_cutoff", label = ("Set score cut-off"), min = 0,
                          max = 100, value = 50),
              
              sliderInput("ppm_cutoff", label = ("Set tolerance in ppm"), min = 0,
                          max = 100, value = 15),
              
              sliderInput("TopN_Feat", label = ("Set TopN Features to report"), min = 0,
                          max = 100, value = 5),
              
              checkboxInput("Peptideatlas_mapping", label = "Use Peptideatlas Data", value = TRUE),
              
              
              
              
       ),column(3,actionButton(inputId="PRM_processing_run", icon=NULL, label="Start PRM candidate selection"),
               
               br())
       )
   }
 ),icon = icon("circle-dot")),
 
        tabPanel("Task manager", wellPanel(
            h5(paste0('Available Cores: ', future::availableCores())),
            column(3,h5('Running task:')),column(3,actionButton('Refresh_Running_task',"Refresh")),
            DT::dataTableOutput(NUM_ASYNC_TASKS_RUNNING),
            numericInput("duration", "Duration", value = 5, min = 0),
            actionButton("start_proc_future", h5("get data using future")),
            h5("Task console (most recent first)"),
            DT::dataTableOutput(DEBUG_CONSOLE_OUTPUT),
            uiOutput("hidden_downloads")
        ),icon = icon("stack-overflow"))
    )
 )
}
