
# code for testing purpose

RUN_imagine_workflow<-function(){
  #Spatial_Quant(ppm=15,adducts=c("M-H","M+Cl"),Quant_list="lipid candidates.csv",cal.mz=T)
  #Spatial_Quant(ppm=15,adducts=c("M-H","M+Cl"),Quant_list="lipid candidates.csv",cal.mz=T,Spectrum_feature_summary = T,Region_feature_summary = T,PMF_analysis = F)
  
  imaging_identification(
    #==============Choose the imzml raw data file(s) to process
    datafile=tk_choose.files(filter =  matrix(c( "imzml file", ".imzML",
                                                 "Text", ".txt", "All files", "*"),
                                              3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis"),
    threshold=0.05, 
    ppm=20,
    Digestion_site="[G]",
    missedCleavages=0:2,
    Fastadatabase="murine_matrisome.fasta",
    adducts=c("M+H","M+NH4","M+Na"),
    #==============TRUE if you want to perform PMF analysis
    PMF_analysis=T,
    #==============TRUE if you want to generate peptide summary in the Summary folder
    Protein_feature_summary=T,
    Peptide_feature_summary=T,
    #==============TRUE if you want to plot peptide in the Ion images folder, make sure there's imzml file in the folder
    plot_ion_image=F,
    #==============Set a number if you want a parallel processing
    parallel=8,
    #==============Set a number (1 to maximum pixels in the data file) if you want to dig more peaks in the raw data
    spectra_segments_per_file=5,
    #==============True if you want the mean spectrum generated based on spatial spectrum profile
    spatialKMeans=T,
    Smooth_range=1
  )
}

imaging_default_parameter<-function(){
  datafile=tk_choose.files(filter = matrix(c( "imzml file", ".imzML",
                                              "Text", ".txt", "All files", "*"),
                                           3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis")
  
  threshold=0.05
  ppm=10
  Digestion_site="[G]"
  missedCleavages=0:2
  Fastadatabase="murine_matrisome.fasta"
  adducts=c("M+H","M+NH4","M+Na")
  #==============TRUE if you want to perform PMF analysis
  PMF_analysis=T
  #==============TRUE if you want to generate peptide summary in the Summary folder
  Protein_feature_summary=T
  Peptide_feature_summary=T
  #==============TRUE if you want to plot peptide in the Ion images folder, make sure there's imzml file in the folder
  plot_ion_image=T
  #==============Set a number if you want a parallel processing
  parallel=8
  #==============Set a number (1 to maximum pixels in the data file) if you want to dig more peaks in the raw data
  spectra_segments_per_file=5
  spatialKMeans=T
  plot_cluster_image=T
}

run_spatial_quant_workflow<-function(){
  imagine_Spatial_Quant(datafile=tk_choose.files(filter =  matrix(c( "imzml file", ".imzML",
                                                                     "Text", ".txt", "All files", "*"),
                                                                  3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis"),
                        threshold=0.00, 
                        ppm=5,
                        #Fastadatabase="murine_matrisome.fasta",
                        Quant_list="Metabolites of Interest.csv",
                        adducts=c("M-H","M+Cl"),
                        cal.mz=F,
                        mzlist_bypass=T,
                        #==============TRUE if you want to plot protein PMF result
                        PMF_analysis=TRUE,
                        #==============TRUE if you want to generate peptide summary in the Summary folder
                        Protein_feature_summary=F,
                        Peptide_feature_summary=F,
                        #PLOT_PMF_Protein=FALSE,
                        #==============TRUE if you want to plot peptide in the Ion images folder, make sure there's imzml file in the folder
                        plot_ion_image=FALSE,
                        #==============Set a number if you want a parallel processing
                        parallel=detectCores()/2,
                        #==============Set a number (1 to maximum pixels in the data file) if you want to dig more peaks in the raw data
                        spectra_segments_per_file=4,
                        spatialKMeans=F,
                        Smooth_range=1,
                        Virtual_segmentation=T,
                        Virtual_segmentation_rankfile="Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy\\radius_rank_bovin.csv",
                        Spectrum_feature_summary=F,
                        Region_feature_summary=F,
                        Region_feature_analysis=F,
                        Cluster_level=="Low")
  
  
  #workflow for aged human lens
  imagine_Spatial_Quant(
    #==============Choose the imzml raw data file(s) to process  make sure the fasta file in the same folder
    datafile=tk_choose.files(filter =  matrix(c( "imzml file", ".imzML",
                                                 "Text", ".txt", "All files", "*"),
                                              3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis"),
    threshold=0.00, 
    ppm=20,
    #Fastadatabase="murine_matrisome.fasta",
    Quant_list="lipid candidates.csv",
    adducts=c("M-H","M+Cl"),
    cal.mz=T,
    mzlist_bypass=F,
    #==============TRUE if you want to plot protein PMF result
    PMF_analysis=TRUE,
    #==============TRUE if you want to generate protein summary in the Summary folder
    Protein_feature_summary=T,
    #==============TRUE if you want to generate protein cluster image in the Summary folder
    plot_cluster_image=T,
    Peptide_feature_summary=T,
    #PLOT_PMF_Protein=FALSE,
    #==============TRUE if you want to plot peptide in the Ion images folder, make sure there's imzml file in the folder
    plot_ion_image=FALSE,
    #==============Set a number if you want a parallel processing
    parallel=detectCores(),
    #==============Set a number (1 to maximum pixels in the data file) if you want to dig more peaks in the raw data
    spectra_segments_per_file=5,
    spatialKMeans=F,
    Smooth_range=1,
    Virtual_segmentation=T,
    Virtual_segmentation_rankfile="Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy\\radius_rank.csv",
    Spectrum_feature_summary=T,
    Region_feature_summary=T,
    Region_feature_analysis=T,
    plot_each_metabolites=F,
    Cluster_level="High",
    Region_feature_analysis_bar_plot=F,
    norm_datafiles=T,
    norm_Type="Median"
  )
  
  imagine_Spatial_Quant(
    #==============Choose the imzml raw data file(s) to process  make sure the fasta file in the same folder
    datafile=tk_choose.files(filter =  matrix(c( "imzml file", ".imzML",
                                                 "Text", ".txt", "All files", "*"),
                                              3, 2, byrow = TRUE),caption  = "Choose single or multiple file(s) for analysis"),
    threshold=0.00, 
    ppm=25,
    #Fastadatabase="murine_matrisome.fasta",
    Quant_list="lipid blast candidates.csv",
    adducts=c("M-H","M+Cl"),
    cal.mz=F,
    mzlist_bypass=F,
    #==============TRUE if you want to plot protein PMF result
    PMF_analysis=TRUE,
    #==============TRUE if you want to generate protein summary in the Summary folder
    Protein_feature_summary=T,
    #==============TRUE if you want to generate protein cluster image in the Summary folder
    plot_cluster_image=T,
    Peptide_feature_summary=T,
    #PLOT_PMF_Protein=FALSE,
    #==============TRUE if you want to plot peptide in the Ion images folder, make sure there's imzml file in the folder
    plot_ion_image=FALSE,
    #==============Set a number if you want a parallel processing
    parallel=detectCores(),
    #==============Set a number (1 to maximum pixels in the data file) if you want to dig more peaks in the raw data
    spectra_segments_per_file=5,
    spatialKMeans=F,
    Smooth_range=1,
    Virtual_segmentation=T,
    Virtual_segmentation_rankfile="Z:\\George skyline results\\maldiimaging\\Maldi_imaging - Copy\\radius_rank.csv",
    Spectrum_feature_summary=T,
    Region_feature_summary=T,
    Region_feature_analysis=T,
    plot_each_metabolites=F,
    Cluster_level="High",
    Region_feature_analysis_bar_plot=F,
    norm_datafiles=T
  )
}

Cardinal<- function(){
  workdir <- WorkingDir()
  datafile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.imzML$")
  
  listfile <- list.files(path=workdir,full.names = TRUE,pattern = "\\.csv$")
  datafile<-gsub(".imzML", "", datafile)
  name <-gsub(paste(workdir,"\\/",sep=""),"",datafile[3])
  
  #imdata <- readImzML(name, folder, as="MSImagingExperiment")
  imdata <- readImzML(name, folder, attach.only=TRUE, as="MSImagingExperiment")
  summarize(data, sum, .by="pixel")
  tmp <- imdata %>%
    smoothSignal() %>%
    reduceBaseline() %>%
    peakPick() %>%
    peakFilter() %>%
    select(x == 1, y == 1)
  process(plot=TRUE,
          par=list(layout=c(1,3)),
          BPPARAM=SerialParam())
  
  plot(data, pixel=1)
  plot(data, coord=list(x=2, y=2))
  
  
  pattern <- factor(c(0, 0, 2, 2, 0, 0, 0, 0, 0, 0, 2, 2, 0,
                      0, 0, 0, 0, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 2, 1, 1, 2,
                      2, 0, 0, 0, 0, 0, 1, 2, 2, 2, 2, 0, 0, 0, 0, 1, 2, 2,
                      2, 2, 2, 0, 0, 0, 0, 2, 2, 2, 2, 2, 2, 2, 0, 0, 0, 2,
                      2, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0),
                    levels=c(0,1,2), labels=c("blue", "black", "red"))
  tmp <- imdata %>%
    smoothSignal() %>%
    reduceBaseline() %>%
    peakPick() %>%
    peakFilter() %>%
    select(x == 1, y == 1) %>%
    process(plot=TRUE,
            par=list(layout=c(1,3)),
            BPPARAM=SerialParam())
}