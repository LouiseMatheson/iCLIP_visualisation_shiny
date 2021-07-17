# iCLIP_visualisation_shiny

Shiny app to visualise CLIP data

## preprocess_app_data.R

This script provides the code used to process the data that goes into the shiny app - both of a GTF file, and of some CLIP crosslink and cluster data. This can be tested on the sample iCLIP crosslink data (iCLIP...chr17_100kb.txt files) and GTF (Mouse_GRCm38.90_chr17_100kb.gtf) provided here:

https://github.com/LouiseMatheson/Process_CLIP_data/blob/main/test_data/

These files encompass chr17:35150000-35250000 (GRCm38), including the Tnf and Lta genes where we see accumulation of CLIP crosslinks in these datasets. 

Data for clusters used in this script, which encompass the same regions, is in the test_data folder of this project. 


## iCLIP_shiny_app

This folder contains the code underlying the shiny app used to visualise CLIP data. The folder includes all of the data required to run the app in order to visualise the CLIP crosslinks over the genes within chr17:35150000-35250000 (GRCm38). 

Note that since only a small sample of data is provided to test the app and show how it runs, some of the underlying code is unnecessary, but enables more flexibility if additional datasets are added.

The app uses the shiny, RColorBrewer, tidyverse and htmltools packages (available from CRAN with install.packages), and BSgenome.Mmusculus.UCSC.mm10 which is available through Bioconductor (BiocManager::install). 

Once these packages are installed, the app can be run by opening the app.R script in RStudio, and clicking on 'Run App'. 

A screenshot of how this should look, displaying Tnf together with the TATTTA motif, is shown below:

![alt text](https://github.com/LouiseMatheson/iCLIP_visualisation_shiny/blob/main/iCLIP_shiny_app/iCLIP_visualisation_app_Tnf.png?raw=true)