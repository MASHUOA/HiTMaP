---
title: "HiTMaP"
output:
  html_document: 
    keep_md: yes
    toc: yes
    theme: readable
    number_sections: yes
    df_print: tibble
    highlight: tango
  pdf_document: default
  word_document: default
---
-- An R package of High-resolution Informatics Toolbox for Maldi-imaging Proteomics


## Package installation
This is an tutorial for use of HiTMaP (An R package of High-resolution Informatics Toolbox for Maldi-imaging Proteomics). To access the software use the installation codes as below: 


```r
#install the git package
install.packages("devtools")
Sys.setenv(R_REMOTES_NO_ERRORS_FROM_WARNINGS=T)
library(devtools)
install_github("guoguodigit/HiTMaP",auth_token ="a124a067ed1c84f8fd577c972845573922f1bb0f",force=T)
#Update all dependencies
1
library(HiTMaP)
```


## Proteomics identification on Maldi imaging data file 

Now the HiTMaP is upon running. You could build the candidate list of your target proteome and perform image identification by using the function as below:


```r
#creat candidate list
library(HiTMaP)

imaging_identification(
#==============Choose the imzml raw data file(s) to process  make sure the fasta file in the same folder
               datafile=paste(file.path(path.package(package="HiTMaP")),"/data/Bovin_lens.imzML", sep=""),
               threshold=0.005, 
               ppm=5,
#==============specify the digestion enzyme specificity
               Digestion_site="([KR](?=[^P]))|((?<=W)K(?=P))|((?<=M)R(?=P))",
#==============specify the range of missed Cleavages
               missedCleavages=0:1,
#==============Set the target fasta file
               Fastadatabase="uniprot-bovin.fasta",
#==============Set the possible adducts and fixed modifications
               adducts=c("M+H"),
               Modifications=list(fixed=NULL),
#==============The decoy mode: could be one of the "adducts", "elements" or "isotope"
               Decoy_mode = "isotope",
               use_previous_candidates=F,
               output_candidatelist=T,
#==============Set the parameters for image segmentation
               spectra_segments_per_file=5,
               spatialKMeans=TRUE,
               Smooth_range=1,
               Virtual_segmentation=FALSE,
               Virtual_segmentation_rankfile=NULL,
#==============Set the Score method for hi-resolution isotopic pattern matching
               score_method="SQRTP",
#==============Summarise the protein and peptide features across the project the result can be found at the summary folder
               Protein_feature_summary=TRUE,
               Peptide_feature_summary=TRUE,
               Region_feature_summary=FALSE,
#==============The decoy mode: could be one of the "adducts", "elements" or "isotope"
               plot_cluster_image_grid=FALSE,
               ClusterID_colname="Protein",
               componentID_colname="Peptide",
               Rotate_IMG=NULL,
               )
```


## Project folder and result structure 

In the above function, You have performed proteomics analysis of the sample data file. It is a tryptic Bovin lens MALDI-imaging file which is acquired on an FT-ICR MS.
The function will take the selected data files' root directory as the project folder.
In this example, the project folder will be:


```r
library(HiTMaP)
wd=paste(file.path(path.package(package="HiTMaP")),"/data/", sep="")
wd
```

```
## [1] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data/"
```


After the whole identification process, we will get two types of sub-folders in the project folder:


```r
list.dirs(wd, recursive=FALSE)
```

```
## [1] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data//Bovin_lens ID" 
## [2] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data//Summary folder"
```

1. The one which has an identical name to an input data file contains the identification result of that input:
   + the protein and peptides list of each segmentation region
   + the PMF matching plot of each segmentation
   + the image that indicates segmentations' boundary (applies to either K-mean segmentation (powered by Cardinal) or manually defined segmentation)
   + folders of each region contains the detailed identification process, FDR plots and isotopic pattern matching plots

2. "Summary folder" contains: 
   + the identification summary of protein and peptides across all the data
   + the candidate list of all possible proteins and peptides (if *use_previous_candidates* is set as **TRUE**)
   + the Cluster imaging files of the protein of interest
   
   
## Identification result visulasation and interpretation

Now we could visualize the result by the following functions:

To check the segmentation result over the sample, you need got to each data file ID folder and find the "spatialKMeans_image_plot.png" (if you are using the spatial K-means method for segmentation.)


```r
library(magick)
```

```
## Linking to ImageMagick 6.9.9.14
## Enabled features: cairo, freetype, fftw, ghostscript, lcms, pango, rsvg, webp
## Disabled features: fontconfig, x11
```

```r
p<-image_read(paste0(wd,"/Bovin_lens ID/spatialKMeans_image_plot.png"))
print(p)
```

```
## # A tibble: 1 x 7
##   format width height colorspace matte filesize density
##   <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
## 1 PNG     1024    720 sRGB       FALSE     9656 72x72
```

<img src="README_files/figure-html/VisulazeKmean-1.png" width="1024" />

The pixels in image data now has been categorized into five regions according to the initial setting of segmentation (*spectra_segments_per_file=5*). The rainbow shaped bovine lens segmentation image (on the left panel) shows a unique statistical classification based on the mz features of each region (on the right panel).

The identification will take place on the **mean spectra** of each region. To check the peptide mass fingerprint (PMF) matching quality, 
you could locate the PMF spectrum matching plot of each individual region.


```r
library(magick)
p_pmf<-image_read(paste0(wd,"/Bovin_lens ID/Bovin_lens 3PMF spectrum match.png"))
print(p_pmf)
```

```
## # A tibble: 1 x 7
##   format width height colorspace matte filesize density
##   <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
## 1 PNG     1980   1080 sRGB       FALSE    17148 72x72
```

<img src="README_files/figure-html/unnamed-chunk-1-1.png" width="1980" />

list of Peptides and proteins of each region has also been created so that you may check each individual region's result.


```r
peptide_pmf_result<-read.csv(paste0(wd,"/Bovin_lens ID/Peptide_segment_PMF_RESULT_3.csv"))
head(peptide_pmf_result)
```

```
## # A tibble: 6 x 15
##   Protein    mz Peptide adduct formula isdecoy pepmz charge Intensity
##     <int> <dbl> <fct>   <fct>  <fct>     <int> <dbl>  <int>     <dbl>
## 1   10450 1040. GPSSEC~ M+H    C47H66~       0 1039.      1   960902.
## 2   10450 1451. NTDPIT~ M+H    C66H10~       0 1450.      1    73014.
## 3   10726 1585. ASGCLI~ M+H    C66H11~       0 1584.      1    77588.
## 4   10726 1460. YEYILA~ M+H    C64H99~       0 1459.      1   250624.
## 5   11103 1616. LVKGNY~ M+H    C75H11~       0 1615.      1   252569.
## 6   11203 1464. IDPSAS~ M+H    C61H99~       0 1463.      1    91876.
## # ... with 6 more variables: moleculeNames <fct>, Region <int>,
## #   Score <dbl>, mz_align <dbl>, Rank <int>, desc <fct>
```


```r
protein_pmf_result<-read.csv(paste0(wd,"/Bovin_lens ID/Protein_segment_PMF_RESULT_3.csv"))
head(protein_pmf_result)
```

```
## # A tibble: 6 x 9
##   Protein Proscore isdecoy Intensity Score peptide_count Protein_coverage
##     <int>    <dbl>   <int>     <dbl> <dbl>         <int>            <dbl>
## 1   10450   0.113        0   516958. 0.593             2           0.172 
## 2   10726   0.0749       0   164106. 0.979             2           0.0759
## 3   11103   0.307        0   252569. 2.26              1           0.130 
## 4   11203   0.188        0  1632393. 2.01              2           0.0777
## 5   11209   0.0701       0   116368. 0.964             2           0.0742
## 6   11264   0.0829       0   584037. 0.770             3           0.0965
## # ... with 2 more variables: Intensity_norm <dbl>, desc <fct>
```

**Score** in peptide result table shows the isotopic pattern matching score of the peptide. In Protein result table, it shows the intensity weighted peptide spectrum matching score.

**Proscore** in the protein result table shows the overall estimation of the protein identification Accuracy

A *Peptide_region_file.csv* has also been created to summarise all the IDs in this data file:


```r
Identification_summary_table<-read.csv(paste0(wd,"/Bovin_lens ID/Peptide_region_file.csv"))
head(Identification_summary_table)
```

```
## # A tibble: 6 x 15
##   Protein    mz Peptide adduct formula isdecoy pepmz charge Intensity
##     <int> <dbl> <fct>   <fct>  <fct>     <int> <dbl>  <int>     <dbl>
## 1   10562 1735. IGFEEK~ M+H    C73H11~       0 1734.      1    74414.
## 2    1137 1590. VLGRTG~ M+H    C63H11~       0 1589.      1    75190.
## 3   17261 1509. GLGHGY~ M+H    C69H98~       0 1508.      1  3482147.
## 4    1734 2220. LQRISQ~ M+H    C90H15~       0 2219.      1   239088.
## 5   19661 1375. KFPFTL~ M+H    C67H10~       0 1374.      1   144208.
## 6   19661 1714. SIEHSG~ M+H    C76H11~       0 1713.      1   103318.
## # ... with 6 more variables: moleculeNames <fct>, Region <int>,
## #   Score <dbl>, mz_align <dbl>, Rank <int>, desc <fct>
```

The details of protein/peptide identification process has been save to the folder named by the segmentation:


```r
list.dirs(paste0(wd,"/Bovin_lens ID/"), recursive=FALSE)
```

```
## [1] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data//Bovin_lens ID//1"
## [2] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data//Bovin_lens ID//2"
## [3] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data//Bovin_lens ID//3"
## [4] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data//Bovin_lens ID//4"
## [5] "C:/Users/admgguo484/Documents/R/win-library/3.6/HiTMaP/data//Bovin_lens ID//5"
```
In the identification details folder, you will find a series of FDR file and plots to demonstrate the FDR model and score cutoff threshold:


```r
dir(paste0(wd,"/Bovin_lens ID/1/"), recursive=T)
```

```
##  [1] "FDR.CSV"                                        
##  [2] "FDR.png"                                        
##  [3] "Matching_Score_vs_mz_target-decoy.png"          
##  [4] "Peptide_1st_ID.csv"                             
##  [5] "Peptide_1st_ID_score_rank_SQRTP.csv"            
##  [6] "Peptide_2nd_ID_score_rankSQRTP_Rank_above_3.csv"
##  [7] "Peptide_Score_histogram_target-decoy.png"       
##  [8] "PROTEIN_FDR.CSV"                                
##  [9] "protein_FDR.png"                                
## [10] "Protein_ID_score_rank_SQRTP.csv"                
## [11] "PROTEIN_Score_histogram.png"                    
## [12] "Spectrum.csv"
```

In this folder, you will find the FDR plots for protein and peptide. The software will take the proscore and its FDR model to trim the final identification result.


```r
library(magick)
p_FDR_peptide<-image_read(paste0(wd,"/Bovin_lens ID/3/FDR.png"))
p_FDR_protein<-image_read(paste0(wd,"/Bovin_lens ID/3/protein_FDR.png"))
p_FDR_peptide_his<-image_read(paste0(wd,"/Bovin_lens ID/3/Peptide_Score_histogram_target-decoy.png"))
p_FDR_protein_his<-image_read(paste0(wd,"/Bovin_lens ID/3/PROTEIN_Score_histogram.png"))
p_combined<-image_append(c(p_FDR_peptide,p_FDR_peptide_his,p_FDR_protein,p_FDR_protein_his))
print(p_combined)
```

```
## # A tibble: 1 x 7
##   format width height colorspace matte filesize density
##   <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
## 1 PNG     1920    480 sRGB       FALSE        0 72x72
```

<img src="README_files/figure-html/FDR plot-1.png" width="1920" />



## Identification summary and cluster imaging

In the project summary folder, you will find four files and a sub-folder:

```r
wd_sum=paste(file.path(path.package(package="HiTMaP")),"/data/Summary folder", sep="")
dir(wd_sum)
```

```
## [1] "candidatelist.csv"   "cluster Ion images"  "Peptide_Summary.csv"
## [4] "protein_index.csv"   "Protein_Summary.csv"
```

"candidatelist.csv" and "protein_index.csv" contains the candidates used for this project. They are output after the candidate processing while *output_candidatelist* set as TRUE, and can be used repeatedly while *use_previous_candidates* set as TRUE.

"Peptide_Summary.csv" and "Protein_Summary.csv" contains the table of the project identification summary. You could set the *plot_cluster_image_grid* as TRUE to enable the cluster imaging function. Please be noted that you could indicate *Rotate_IMG* with a CSV file path that indicates the rotation degree of image files. 

**Note**: 90$^\circ$, 180$^\circ$ and 270$^\circ$ are recommended for image rotation. You may find an example CSV file in the library/HiTMaP/data folder.

Now you could visualized the interest proteins and their associated peptides' distribution via cluster imaging function.


```r
p_cluster1<-image_read(paste0(wd,"/Summary folder/cluster Ion images/705_cluster_plot_sum_flex.png"))
print(p_cluster1)
```

```
## # A tibble: 1 x 7
##   format width height colorspace matte filesize density
##   <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
## 1 PNG     9360   3060 sRGB       TRUE    582848 118x118
```

<img src="README_files/figure-html/CLuster imaging-1.png" width="9360" />

```r
#p_cluster1<-image_read(paste0(wd,"/Summary folder/cluster Ion images/705_footer.png"))
##print(p_cluster1)
p_cluster2<-image_read(paste0(wd,"/Summary folder/cluster Ion images/5027_cluster_plot_sum_flex.png"))
print(p_cluster2)
```

```
## # A tibble: 1 x 7
##   format width height colorspace matte filesize density
##   <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
## 1 PNG    15600   3060 sRGB       TRUE    676601 118x118
```

<img src="README_files/figure-html/CLuster imaging-2.png" width="15600" />

```r
#p_cluster2<-image_read(paste0(wd,"/Summary folder/cluster Ion images/5027_footer.png"))
#print(p_cluster2)
p_cluster3<-image_read(paste0(wd,"/Summary folder/cluster Ion images/5479_cluster_plot_sum_flex.png"))
print(p_cluster3)
```

```
## # A tibble: 1 x 7
##   format width height colorspace matte filesize density
##   <chr>  <int>  <int> <chr>      <lgl>    <int> <chr>  
## 1 PNG    15600   3060 sRGB       TRUE    673817 118x118
```

<img src="README_files/figure-html/CLuster imaging-3.png" width="15600" />

```r
#p_cluster3<-image_read(paste0(wd,"/Summary folder/cluster Ion images/5479_footer.png"))
#print(p_cluster3)
```

End of the tutorial, Enjoy~
