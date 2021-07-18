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

Using the provided test data, preprocessing should take < 2 min to run on a standard desktop computer. The app should launch in < 1 min and, since the dataset is very small, respond to inputs almost immediately.

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


## Adding complete CLIP datasets to the app

Fastq files for the CLIP data generated and/or analysed in our manuscript (https://www.biorxiv.org/content/10.1101/2021.06.03.446738v2) are available via the GEO accessions GSE176313 (ZFP36L1 CLIP) and GSE96074 (ZFP36 CLIP).

Before they can be added to the app, significant crosslink sites and clusters should first be identified using the iCount pipeline (documentation and instructions: https://icount.readthedocs.io/en/latest/).

The required output files from iCount are:

1. tab-delimited text files from iCount peaks, generated using the --scores parameter (ie, including all crosslink sites, together with the score and FDR)
2. BED files output from iCount clusters.

A complete GTF file for the same genome build as used by iCount is also required (eg one option for GRCm38: http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz)

The following modifications then need to be made to the script to incorporate the full data files:
1. Replace sample GTF file name with full GTF
2. Replace filenames in read_tsv for sample iCLIP crosslink data with the full datasets; duplicate as many times as necessary to include all data
3. bind_rows() of all crosslink datasets that have been imported, adding appropriate sample names (optional: find and replace 'all_test_CLIP_data' with 'all_CLIP_data' throughout)
4. Modify the path where the individual chromosome data from all_CLIP_data is saved to be CLIP_data (rather than test_CLIP_data). If additional types of data are imported, these can be separated out and saved into separate subfolders.
5. Modify the code for creating CLIP_datasets to ensure that a) each sample name is assigned to the correct 'dataset' (ie, group of replicates), and appropriate names are provided for each dataset; b) if any additional types of data have been added, these are correctly picked up as a separate 'Type'; c) all samples for which there are replicate datasets are picked up and 'replicates' is set to T for these.
6. If additional data types and therefore subfolders containing the CLIP data have been added, ensure the subfolder paths are defined for each Type of data in CLIPdata_paths. 
7. Replace filenames for sample iCLIP cluster data with the full datasets, and set the 'dataset' to correspond to the Sample name assigned to the crosslinks for that replicate; duplicate bind_rows(add_column(read_tsv(...))) rows as many times as necessary to include all data.
8. For datasets with replicates, duplicate and modify the code for creating GRanges objects and merging based on intersections to perform this for each group of replicate data. Where more than 3 replicates are present, all pairwise intersections should be included. In this case, 'dataset' should correspond to the dataset name assigned to the crosslinks.


Once preprocessing has been completed, no modification to the app.R code should be necessary, and the full datasets should now be detected when it is run. 
