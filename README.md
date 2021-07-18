# iCLIP_visualisation_shiny

Shiny app to visualise CLIP data

## Requirements and installation

Installation of the required software and packages should take < 30 min.

Depends on R and RStudio, freely available to install on Windows, macOS or Linux here: 
https://cran.r-project.org/index.html
https://www.rstudio.com/products/rstudio/

The code has been tested with R >= v4.0.0 and <= v4.1.0, and corresponding package versions.

To run the preprocessing and app, packages will need to be installed as follows:

install.packages("shiny") # app only
install.packages("RColorBrewer") # app only
install.packages("tidyverse")
install.packages("htmltools") # app only
install.packages("BiocManager")
BiocManager::install("GenomicRanges") # preprocessing only
BiocManager::install("rtracklayer") # preprocessing only
BiocManager::install("BSgenome.Mmusculus.UCSC.mm10") # app only


## Instructions and test data

### preprocess_app_data.R

This script provides the code used to process the data that goes into the shiny app - both of a GTF file, and of some CLIP crosslink and cluster data. This can be tested on the sample iCLIP crosslink data (iCLIP...chr17_100kb.txt files) and GTF (Mouse_GRCm38.90_chr17_100kb.gtf) provided here:

https://github.com/LouiseMatheson/Process_CLIP_data/blob/main/test_data/

These files encompass chr17:35150000-35250000 (GRCm38), including the Tnf and Lta genes where we see accumulation of CLIP crosslinks in these datasets. 

Data for clusters used in this script, which encompass the same regions, is in the test_data folder of this project. 

Whilst replicates from a single dataset and type of data are included as a sample here, annotations to the code indicate how additional datasets and types of data can be incorporated (including example lines of code which will be ignored since the datasets are not present here).


### iCLIP_shiny_app

This folder contains the code underlying the shiny app used to visualise CLIP data. The folder includes all of the data required to run the app in order to visualise the CLIP crosslinks over the genes within chr17:35150000-35250000 (GRCm38). 

Note that since only a small sample of data is provided to test the app and show how it runs, some of the underlying code is unnecessary, but enables more flexibility if additional datasets are added.

The app can be run by opening the app.R script in RStudio, and clicking on 'Run App'. 

A screenshot of how this should look, displaying Tnf together with the TATTTA motif, is shown below:

![alt text](https://github.com/LouiseMatheson/iCLIP_visualisation_shiny/blob/main/iCLIP_shiny_app/iCLIP_visualisation_app_Tnf.png?raw=true)