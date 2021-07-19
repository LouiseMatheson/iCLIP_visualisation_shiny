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

Note that the order of loading packages matters - therefore separate sessions should be used to run the preprocessing and the app.

### preprocess_app_data.R

This script provides the code used to process the data that goes into the shiny app - both of a GTF file, and of some CLIP crosslink and cluster data.
It should be run from the parent directory to the 'iCLIP_shiny_app' folder that contains the app, to ensure that all data is saved to the correct locations.

It can be tested on the sample iCLIP crosslink data (iCLIP...chr17_100kb.txt files), cluster data (iCLIP_clusters...txt files) and GTF (Mouse_GRCm38.90_chr17_100kb.gtf) provided here:

https://github.com/LouiseMatheson/Process_CLIP_data/blob/main/test_data/

A completed sample template describing each sample is also required, and is provided for the test data here:
https://github.com/LouiseMatheson/Process_CLIP_data/blob/main/test_data/CLIPdata_template_test.txt

These files encompass chr17:35150000-35250000 (GRCm38), including the Tnf and Lta genes where we see accumulation of CLIP crosslinks in these datasets. 

The script is written to use the test GTF and sample template, so can be run directly; the outputs should overwrite the following files contained with the iCLIP_shiny_app/ folder with identical files:
CLIP_tool_data.RData
test_CLIP_data/merged_CLIP_chr17.Rds
Cluster_data/merged_clusters_chr17.Rds

Note that some warning messages are expected; the following warnings do not indicate a problem:

Problem with `mutate()` column `transcript_support_level`.
i `transcript_support_level = as.numeric(sub(" .*", "", transcript_support_level))`.
i NAs introduced by coercion 

The following named parsers don't match the column names: chromosome 

The following named parsers don't match the column names: chrom 


### iCLIP_shiny_app

This folder contains the code underlying the shiny app used to visualise CLIP data. The folder includes all of the data required to run the app in order to visualise the CLIP crosslinks over the genes within chr17:35150000-35250000 (GRCm38). 

Note that since only a small sample of data is provided to test the app and show how it runs, some of the underlying code is unnecessary, but enables more flexibility if additional datasets are added.

The app can be run by opening the app.R script in RStudio, and clicking on 'Run App'. 

A screenshot of how this should look, displaying Tnf together with the TATTTA motif, is shown below:

![alt text](https://github.com/LouiseMatheson/iCLIP_visualisation_shiny/blob/main/iCLIP_shiny_app/iCLIP_visualisation_app_Tnf.png?raw=true)


## Adding complete CLIP datasets to the app

Fastq files for the CLIP data generated and/or analysed in our manuscript (https://www.biorxiv.org/content/10.1101/2021.06.03.446738v2) are available via the GEO accessions GSE176313 (ZFP36L1 CLIP) and GSE96074 (ZFP36 CLIP).

Before they can be added to the app, significant crosslink sites and clusters should first be identified using the iCount pipeline (documentation and instructions: https://icount.readthedocs.io/en/latest/).

The output files from iCount used are:

1. tab-delimited text files from iCount peaks, generated using the --scores parameter (ie, including all crosslink sites, together with the score and FDR). Alternatively, a BED file with just the significant crosslink sites can be used (but will not accommodate filtering on FDR)
2. BED files output from iCount clusters (optional)

A complete GTF file for the same genome build as used by iCount is also required (eg one option for GRCm38: http://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz)


To incorporate complete datasets, the sample template needs to be completed for each sample. The order of samples/datasets determines the order in which they will be listed for selection, as well as the order they are plotted (the first at the bottom). The following columns should be included (case unimportant):

1. xlink_path - path to the crosslink data. The output from iCount (v2.0.1dev) peaks with all crosslink sites, including FDRs, should be the correct format. Other sources can be used, as long as the required columns are present (order/case unimportant), or BED files in which the score/number of crosslinks is assumed to be in column 5:
* chromosome (or 'chrom'; UCSC or Ensembl format can be used)
* position (alternatives 'hit_position' or 'end' - for BED files or other formats with start/end colnames present instead, uses end to account for 0 indexing)
* strand (can be +/- or 1/-1)
* score (alternatives 'crosslinks' or 'value..unknown_counts.')
* FDR (optional: if not present will be set to 0, thereby disabling filtering on FDR for this sample)

2. sample - a unique, simple name/ID for each sample, typically without spaces/special characters

3. label - the sample name to be used in labelling the individual sample/replicate on plots when replicates are not merged. Note that commas are replaced with newline characters on the plots. Can include eg spaces

4. dataset - defines the replicate samples that belong to the same dataset; when replicates are merged, the number of crosslink sites will be added together for all samples belonging to the same dataset. This column is used for the dataset names in the app for selecting which to display, and on the plots - can include eg spaces. Similar to 'label', commas are replaced by newline characters when plotted. For datasets that do not include replicates (ie the dataset contains only one sample), the 'label' and 'dataset' value should be identical

5. out_path - name of folder where the crosslink data is to be saved. This folder needs to be created within the 'iCLIP_shiny_app' folder, but 'iCLIP_shiny_app/' need not be included in the file path in the sample template (it does not matter if it is!). This must be identical for all datasets of the same 'type' (or for all datasets if no type is specified).

6. type (optional) - if different types of data are included, such that you would not wish to display datasets of different type on the same plot (eg, different RBPs/completely separate projects), 'type' defines the datasets that belong together (therefore must be identical for samples within the same dataset), and is the name displayed for selection of which type of data you wish to display. If a large amount of data is included, this has the advantage that only the relevant data will be loaded, thus making it quicker. If you wish to be able to display all data together, this should be the same for all samples, or can be left blank or the column omitted (in which case all will be designated 'Unspecified'). 

7. cluster_path (optional - the column can be omitted altogether, or left blank for some samples) - file path to the cluster data. BED files such as output from iCount clusters are an appropriate input, if the format is not BED, the column names should include (case unimportant) chromosome (or chrom; UCSC or Ensembl format can be used), start, end and strand (+/- or 1/-1). Paths to cluster data are optional, however if included for one replicate in a dataset, it is assumed there will be data for the other datasets. If clusters are only there for one out of a set of replicates, this will cause an error. 

The cluster data is all saved together into the 'Cluster_data' subdirectory of 'iCLIP_shiny_app'. Even if no cluster data is included, this folder must be retained and preprocessing is still expected to result in very small Rds files being exported to this location.


An example of a completed sample template for a more complex setup is provided here:
https://github.com/LouiseMatheson/Process_CLIP_data/blob/main/test_data/CLIPdata_template_complex_example.txt


Once the sample template is complete and appropriate subfolders created within iCLIP_shiny_app/, the paths to this and to the complete GTF file to be used should be inserted to correctly define CLIP_sample_template and gtf_file at the start of the preprocess_app_data.R script (lines 7 and 6 respectively), and the script can then be run (similar to the test data, this needs to be run from the parent directory to iCLIP_shiny_app/). 


Once preprocessing has been completed, no modification to the app.R code should be necessary, and the full datasets should now be detected when it is run. 
