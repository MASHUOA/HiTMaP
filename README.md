---
title: "HiTMaP"
output:
  html_document: 
    keep_md: yes
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
##   format width height colorspace matte filesize density
## 1    PNG  1024    720       sRGB FALSE     9656   72x72
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
##   format width height colorspace matte filesize density
## 1    PNG  1980   1080       sRGB FALSE    17148   72x72
```

<img src="README_files/figure-html/unnamed-chunk-1-1.png" width="1980" />

list of Peptides and proteins of each region has also been created so that you may check each individual region's result.


```r
peptide_pmf_result<-read.csv(paste0(wd,"/Bovin_lens ID/Peptide_segment_PMF_RESULT_3.csv"))
head(peptide_pmf_result)
```

```
##   Protein       mz         Peptide adduct         formula isdecoy    pepmz
## 1   10450 1040.451       GPSSECFWK    M+H  C47H66N11O14S1       0 1039.443
## 2   10450 1450.733    NTDPITNIFYPR    M+H   C66H100N17O20       0 1449.725
## 3   10726 1584.816 ASGCLITLDQHNGKK    M+H C66H114N21O22S1       0 1583.809
## 4   10726 1459.706   YEYILATASADSR    M+H    C64H99N16O23       0 1458.699
## 5   11103 1615.812  LVKGNYGFEAEFNK    M+H   C75H111N18O22       0 1614.804
## 6   11203 1463.724   IDPSASRQGYDVR    M+H    C61H99N20O22       0 1462.716
##   charge Intensity   moleculeNames Region     Score mz_align Rank
## 1      1 960901.55       GPSSECFWK      3 0.9018853 1040.453    5
## 2      1  73013.99    NTDPITNIFYPR      3 0.2122419 1450.729   17
## 3      1  77587.51 ASGCLITLDQHNGKK      3 0.3443412 1584.813   23
## 4      1 250623.55   YEYILATASADSR      3 1.5529821 1459.708    6
## 5      1 252568.51  LVKGNYGFEAEFNK      3 2.2638102 1615.810   11
## 6      1  91875.69   IDPSASRQGYDVR      3 0.6023862 1463.727   11
##                                                                                                                 desc
## 1                     tr|A0A3Q1LYG0|A0A3Q1LYG0_BOVIN Uncharacterized protein OS=Bos taurus OX=9913 GN=UTS2 PE=3 SV=1
## 2                     tr|A0A3Q1LYG0|A0A3Q1LYG0_BOVIN Uncharacterized protein OS=Bos taurus OX=9913 GN=UTS2 PE=3 SV=1
## 3         tr|A0A3Q1M3N1|A0A3Q1M3N1_BOVIN DNA excision repair protein ERCC-8 OS=Bos taurus OX=9913 GN=ERCC8 PE=4 SV=1
## 4         tr|A0A3Q1M3N1|A0A3Q1M3N1_BOVIN DNA excision repair protein ERCC-8 OS=Bos taurus OX=9913 GN=ERCC8 PE=4 SV=1
## 5                   tr|A0A3Q1MVY9|A0A3Q1MVY9_BOVIN Ig-like domain-containing protein OS=Bos taurus OX=9913 PE=4 SV=1
## 6 tr|A0A3Q1LTC6|A0A3Q1LTC6_BOVIN Mitogen-activated protein kinase kinase 4 OS=Bos taurus OX=9913 GN=MAP2K4 PE=3 SV=1
```


```r
protein_pmf_result<-read.csv(paste0(wd,"/Bovin_lens ID/Protein_segment_PMF_RESULT_3.csv"))
head(protein_pmf_result)
```

```
##   Protein   Proscore isdecoy Intensity     Score peptide_count
## 1   10450 0.11276258       0  516957.8 0.5926478             2
## 2   10726 0.07492019       0  164105.5 0.9785716             2
## 3   11103 0.30671809       0  252568.5 2.2638102             1
## 4   11203 0.18797348       0 1632392.8 2.0121646             2
## 5   11209 0.07008142       0  116367.6 0.9637757             2
## 6   11264 0.08294995       0  584037.1 0.7702205             3
##   Protein_coverage Intensity_norm
## 1       0.17213115      1.1053730
## 2       0.07588076      1.0089615
## 3       0.12962963      1.0451896
## 4       0.07772021      1.2019852
## 5       0.07419355      0.9800784
## 6       0.09653465      1.1156240
##                                                                                                                 desc
## 1                     tr|A0A3Q1LYG0|A0A3Q1LYG0_BOVIN Uncharacterized protein OS=Bos taurus OX=9913 GN=UTS2 PE=3 SV=1
## 2         tr|A0A3Q1M3N1|A0A3Q1M3N1_BOVIN DNA excision repair protein ERCC-8 OS=Bos taurus OX=9913 GN=ERCC8 PE=4 SV=1
## 3                   tr|A0A3Q1MVY9|A0A3Q1MVY9_BOVIN Ig-like domain-containing protein OS=Bos taurus OX=9913 PE=4 SV=1
## 4 tr|A0A3Q1LTC6|A0A3Q1LTC6_BOVIN Mitogen-activated protein kinase kinase 4 OS=Bos taurus OX=9913 GN=MAP2K4 PE=3 SV=1
## 5                                     tr|E1BC32|E1BC32_BOVIN Tropomodulin 2 OS=Bos taurus OX=9913 GN=TMOD2 PE=4 SV=3
## 6                        tr|A0A3Q1MZJ9|A0A3Q1MZJ9_BOVIN Apolipoprotein A-IV OS=Bos taurus OX=9913 GN=APOA4 PE=3 SV=1
```

*Score* in peptide result table shows the isotopic pattern matching score of the peptide.
*Proscore* in the protein result table shows the overall estimation of identification Accuracy

A *Peptide_region_file.csv* has also been created to summarise all the IDs in this data file:


```r
Identification_summary_table<-read.csv(paste0(wd,"/Bovin_lens ID/Peptide_region_file.csv"))
head(Identification_summary_table)
```

```
##   Protein       mz              Peptide adduct         formula isdecoy
## 1   10562 1734.829      IGFEEKDIAANEENR    M+H   C73H116N21O28       0
## 2    1137 1589.818      VLGRTGSQGQCTQVR    M+H C63H113N24O22S1       0
## 3   17261 1508.739     GLGHGYGGFGGLGFGR    M+H    C69H98N21O18       0
## 4    1734 2220.063 LQRISQSEDEESIVGDGETK    M+H   C90H151N26O39       0
## 5   19661 1374.713          KFPFTLEVYCK    M+H C67H100N13O16S1       0
## 6   19661 1713.846       SIEHSGWNVWSQKR    M+H   C76H113N24O22       0
##      pepmz charge  Intensity        moleculeNames Region     Score
## 1 1733.822      1   74413.62      IGFEEKDIAANEENR      2 2.3401072
## 2 1588.810      1   75190.35      VLGRTGSQGQCTQVR      2 1.7183337
## 3 1507.732      1 3482147.03     GLGHGYGGFGGLGFGR      2 1.4175674
## 4 2219.055      1  239087.96 LQRISQSEDEESIVGDGETK      2 3.1973421
## 5 1373.705      1  144207.75          KFPFTLEVYCK      2 1.7398698
## 6 1712.838      1  103318.09       SIEHSGWNVWSQKR      2 0.4101354
##   mz_align Rank
## 1 1734.824    2
## 2 1589.813    1
## 3 1508.742    5
## 4 2220.064    1
## 5 1374.715   10
## 6 1713.839   22
##                                                                                                                   desc
## 1 tr|G3X6U5|G3X6U5_BOVIN SH3 domain-binding glutamic acid-rich-like protein OS=Bos taurus OX=9913 GN=SH3BGRL PE=3 SV=1
## 2                              sp|Q56JX6|RS28_BOVIN 40S ribosomal protein S28 OS=Bos taurus OX=9913 GN=RPS28 PE=3 SV=1
## 3                               tr|A0A3Q1MB50|A0A3Q1MB50_BOVIN Uncharacterized protein OS=Bos taurus OX=9913 PE=4 SV=1
## 4                                                   sp|Q17Q87|SNN_BOVIN Stannin OS=Bos taurus OX=9913 GN=SNN PE=3 SV=1
## 5                               tr|A0A3Q1LFU2|A0A3Q1LFU2_BOVIN Uncharacterized protein OS=Bos taurus OX=9913 PE=4 SV=1
## 6                               tr|A0A3Q1LFU2|A0A3Q1LFU2_BOVIN Uncharacterized protein OS=Bos taurus OX=9913 PE=4 SV=1
```



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
##   format width height colorspace matte filesize density
## 1    PNG  9360   1560       sRGB  TRUE   173900 118x118
```

<img src="README_files/figure-html/CLuster imaging-1.png" width="9360" />

```r
p_cluster2<-image_read(paste0(wd,"/Summary folder/cluster Ion images/5027_cluster_plot_sum_flex.png"))
print(p_cluster2)
```

```
##   format width height colorspace matte filesize density
## 1    PNG 15600   1560       sRGB  TRUE   299751 118x118
```

<img src="README_files/figure-html/CLuster imaging-2.png" width="15600" />

```r
p_cluster3<-image_read(paste0(wd,"/Summary folder/cluster Ion images/5479_cluster_plot_sum_flex.png"))
print(p_cluster3)
```

```
##   format width height colorspace matte filesize density
## 1    PNG 15600   1560       sRGB  TRUE   304278 118x118
```

<img src="README_files/figure-html/CLuster imaging-3.png" width="15600" />

End of the tutorial, Enjoy~
