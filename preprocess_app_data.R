
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)

gtf_filepath <- "test_data/Mouse_GRCm38.90_chr17_100kb.gtf"
CLIP_sample_template <- "test_data/CLIPdata_template_test.txt"

# As a test, GTF file and CLIP dataset subsets covering a 100kb region of chromosome 17 can be used - these can be found here:
# https://github.com/LouiseMatheson/Process_CLIP_data/tree/main/test_data
# The region is chr17:35150000-35250000 (GRCm38) and encompasses several genes, including Tnf and Lta over which we expect to see crosslinks

# A sample template providing sample/dataset names and file paths for each sample is provided in test_data/CLIPdata_template_sample.txt.

# The filepaths to the GTF and the sample template describing all samples to be included need to be provided at the start of the script - 
# if all input files are in the correct format, no other changes to the script should be required.


import(gtf_filepath) %>%
  as_tibble() -> gtf_file

gtf_file %>%
  rename(chromosome = seqnames, feature = type) %>%
  mutate(chromosome = as.character(chromosome)) %>%
  mutate(feature = as.character(feature)) %>%
  mutate(strand = as.character(strand)) %>%
  mutate(transcript_support_level = as.numeric(sub(" .*", "", transcript_support_level))) %>%
  filter(transcript_support_level == 1 | !is.na(ccds_id)) %>%
  # to limit the amount of data we have to load, limiting to only TSL = 1 or transcripts with CCDS ID (can be changed if people find things missing!)
  filter(!feature %in% c("Selenocysteine", "start_codon", "transcript")) %>%
  select(chromosome, feature, start, end, strand, gene_id, gene_name, transcript_id, transcript_name, transcript_biotype) %>%
  mutate(start_oriented = if_else(strand == "+", start-0.5, end+0.5)) %>%
  mutate(end_oriented = if_else(strand == "+", end+0.5, start-0.5)) -> gtf_filtered


gtf_filtered  %>% add_column(feature_type = NA_character_) %>% 
  slice(0) -> gtf_filtered_converted

for(t_id in unique(gtf_filtered$transcript_id)) {
  gtf_filtered %>%
    filter(transcript_id == t_id) -> coord_subset
  if("CDS" %in% coord_subset$feature) {
    coord_subset %>%
      filter(feature != "exon") %>%
      mutate(feature_type = if_else(grepl("utr$", feature), "UTR", "CDS")) -> coord_subset
  } else {
    coord_subset %>%
      mutate(feature_type = if_else(feature == "exon", "UTR", feature)) -> coord_subset
  }
  if(coord_subset$strand[1] == "+") {
    coord_subset %>%
      arrange(start_oriented) %>%
      bind_rows(gtf_filtered_converted) -> gtf_filtered_converted
  } else {
    coord_subset %>%
      arrange(desc(start_oriented)) %>%
      bind_rows(gtf_filtered_converted) -> gtf_filtered_converted
  }
}
gtf_filtered_converted %>%
  filter(!grepl("^JH|GL", chromosome)) -> gtf_filtered_converted

gtf_filtered_converted %>%
  mutate(start = start_oriented, end= end_oriented) %>%
  select(1:9,13) -> gtf_filtered_converted

gtf_filtered_converted %>%
  distinct(gene_name, gene_id, transcript_name, transcript_id) %>%
  arrange(transcript_name) %>%
  mutate(transcript_biotype = gtf_filtered$transcript_biotype[match(transcript_id, gtf_filtered$transcript_id)]) -> 
  gene_names_ids

gtf_filtered_converted %>%
  arrange(transcript_name) %>%
  mutate(feature_type = factor(feature_type, levels = c("CDS", "UTR"))) -> gene_coord

gene_coord %>%
  select(gene_name, chromosome, strand, start, end) %>%
  group_by(gene_name) %>%
  mutate(start = if_else(strand == "+", min(start)+0.5, max(start)-0.5)) %>%
  mutate(end = if_else(strand == "+", max(end)-0.5, min(end)+0.5)) %>%
  ungroup() %>%
  distinct() -> gene_summary

# initially did this by gene_id, but there were duplicated gene names with different IDs, 
# whereas no duplicated IDs with different names. Checked all duplicates (only 16 pairs) - all 
# are where multiple IDs map to the same gene, one or more isoform assigned to each. So better 
# to group by gene name since this includes all isoforms for a given gene. Then can identify gene 
# name from gene_names_ids, and use this to access gene coordinates and summary.


### CLIP data - crosslinks and clusters

read_tsv(CLIP_sample_template) %>%
  rename_with(~ tolower(.x)) %>%
  mutate(out_path = sub("iCLIP_shiny_app/", "", out_path)) -> all_CLIP_sample_data
if(!"type" %in% colnames(all_CLIP_sample_data)) {
  all_CLIP_sample_data %>%
    add_column(type = "Unspecified") -> all_CLIP_sample_data
}
if(!"cluster_path" %in% colnames(all_CLIP_sample_data)) {
  all_CLIP_sample_data %>%
    add_column(cluster_path = NA) -> all_CLIP_sample_data
}

all_CLIP_sample_data %>%
  # some have uppercase in app and code below so simpler to rename than change everything..!
  rename(Sample = sample, Type = type) %>%
  mutate(Type = if_else(is.na(Type), "Unspecified", as.character(Type))) %>%
  group_by(Sample) %>%
  mutate(replicates = sum(all_CLIP_sample_data$dataset == dataset)>1) %>%
  mutate(dataset = if_else(replicates == T, dataset, label)) %>% # force dataset to equal label if there aren't replicates 
  ungroup() -> all_CLIP_sample_data


# could add a check and warning if any 'Type' has more than one out_path assigned (since the data will only be loaded from one place at once, would potentially cause issues)
# But for now assume the template is completed correctly!

# initiate empty list to store all cluster data
all_clusters <- list()

for(outpath in unique(all_CLIP_sample_data$out_path)) {
  all_CLIP_sample_data %>%
    filter(out_path == outpath) -> subset_CLIP_sample_data
  lapply(subset_CLIP_sample_data$xlink_path, function(x) {
    if(grepl("bed$|bed.gz$", x)) {
      read_tsv(x, col_names = c("chromosome", "start", "end", "name", "score", "strand"), col_types = cols(chromosome = col_character())) -> CLIPdata_temp
    } else {
      read_tsv(x, col_types = cols(chromosome = col_character(), chrom = col_character())) -> CLIPdata_temp
      # will give a warning as one won't be found, but missing is just ignored
    }
    if(!"FDR" %in% colnames(CLIPdata_temp)) {
      CLIPdata_temp %>%
        add_column(FDR = 0) -> CLIPdata_temp
    }
    CLIPdata_temp %>%
      rename_with(~ tolower(.x)) %>%
      rename(FDR = fdr) %>%
      rename_with(~ sub("^end$", "position", .x)) %>% # means it can have start/end or position in input data; take end since in BED format this is the actual position
      rename_with(~ sub("^hit_position$", "position", .x)) %>% # legacy from files analysed through iCount prior to Genialis
      rename_with(~ sub("^chrom$", "chromosome", .x)) %>% # accepts either chrom or chromosome in input data
      rename_with(~ sub("^score$", "crosslinks", .x)) %>% # accepts either score or crosslinks in input data
      rename_with(~ sub("^value..unknown_counts.$", "crosslinks", .x)) %>% # legacy from files analysed through iCount prior to Genialis
      select(chromosome, position, strand, crosslinks, FDR) %>%
      mutate(chromosome = sub("chr", "", chromosome)) %>%
      mutate(chromosome = sub("^M$", "MT", chromosome)) %>%
      filter(chromosome %in% gtf_filtered_converted$chromosome) -> CLIPdata_temp
    if(is.numeric(CLIPdata_temp$strand)) {
      CLIPdata_temp$strand <- as.character(factor(CLIPdata_temp$strand, levels = c(-1,1), labels = c("-", "+")))
    }
    return(CLIPdata_temp)
  }) -> subset_CLIP_data
  names(subset_CLIP_data) <- subset_CLIP_sample_data$Sample
  bind_rows(subset_CLIP_data, .id = "Sample") -> subset_CLIP_data
  for(ch in unique(gene_coord$chromosome)) {
    subset_CLIP_data %>%
      filter(chromosome == ch) %>%
      saveRDS(paste0("iCLIP_shiny_app/", sub("/$", "", outpath), "/merged_CLIP_chr", ch, ".Rds"))
  }
  
  lapply(na.omit(subset_CLIP_sample_data$cluster_path), function(x) {
    if(grepl("bed$|bed.gz$", x)) {
      read_tsv(x, col_names = c("chromosome", "start", "end", "name", "score", "strand"), col_types = cols(chromosome = col_character())) %>% 
        mutate(start = start + 1) -> Cluster_temp
    } else {
      read_tsv(x, col_types = cols(chromosome = col_character(), chrom = col_character())) %>%
        rename_with(~ tolower(.x)) -> Cluster_temp
    }
    if(is.numeric(Cluster_temp$strand)) {
      Cluster_temp$strand <- as.character(factor(Cluster_temp$strand, levels = c(-1,1), labels = c("-", "+")))
    }
    Cluster_temp %>%
      rename_with(~ sub("^chrom$", "chromosome", .x)) %>% # accepts either chrom or chromosome in input data
      mutate(chromosome = sub("chr", "", chromosome)) %>%
      mutate(chromosome = sub("^M$", "MT", chromosome)) %>%
      filter(chromosome %in% gtf_filtered_converted$chromosome) %>%
      select(chromosome, start, end, strand)
  }) -> subset_cluster_data
  names(subset_cluster_data) <- subset_CLIP_sample_data$Sample[!is.na(subset_CLIP_sample_data$cluster_path)]
  
  
  for(ds in unique(subset_CLIP_sample_data$dataset[subset_CLIP_sample_data$replicates == T & !is.na(subset_CLIP_sample_data$cluster_path)])) {
    lapply(subset_cluster_data[subset_CLIP_sample_data$Sample[subset_CLIP_sample_data$dataset == ds]], function(x) {
      with(x, GRanges(chromosome, IRanges(start, end), strand = strand))
    }) -> clusters_GR
    
    merged_clusters <- GRanges()
    for(i in 1:(length(clusters_GR)-1))  {
      for(j in (i+1):length(clusters_GR)) {
        merged_clusters <- c(merged_clusters, IRanges::intersect(clusters_GR[[i]], clusters_GR[[j]]))
      }
    }
    merged_clusters %>%
      GenomicRanges::reduce() %>%
      as_tibble() %>%
      rename(chromosome = seqnames) %>%
      select(-width) -> subset_cluster_data[[ds]]
  }
  all_clusters <- c(all_clusters, subset_cluster_data)
}


all_CLIP_sample_data %>%
  select(Sample, label, dataset, Type, replicates) -> CLIP_datasets


all_CLIP_sample_data %>%
  distinct(Type, out_path) %>%
  rename(path = out_path) -> CLIPdata_paths

# Clusters - sorting out sample/dataset annotations, since in the list of cluster data, Sample is used to name 
# those for individual replicates, whilst dataset is used to name merged.

# When samples are merged, we need to be selecting based on dataset names from CLIP_datasets, whereas 
# when not merged we need to use sample names. Need separate columns to define these so that both options 
# are there for those that don't have multiple replicates (otherwise you just don't get clusters for those 
# displayed when replicates are merged)

# Also - if there are no clusters for any dataset, need to provide an empty tibble so it does not fail

if(length(all_clusters) == 0) {
  all_clusters <- tibble(dataset = character(), chromosome = character(), start = double(), end = double(), strand = character(), sample = character())
}
all_clusters %>%
  bind_rows(.id = "dataset") %>%
  mutate(sample = if_else(dataset %in% CLIP_datasets$Sample, dataset, NA_character_)) %>%
  mutate(dataset = if_else(dataset %in% CLIP_datasets$dataset, dataset, CLIP_datasets$dataset[match(dataset, CLIP_datasets$Sample)])) -> all_clusters


# But this isn't quite right because for datasets with replicates, we don't want to include the clusters that come from the individual replicates when they are merged.

all_clusters %>% 
  mutate(dataset = if_else(sample %in% CLIP_datasets$Sample[CLIP_datasets$replicates == T], NA_character_, dataset)) -> all_clusters

for(ch in unique(gene_coord$chromosome)) {
  all_clusters %>%
    filter(chromosome == ch) %>%
    saveRDS(paste0("iCLIP_shiny_app/Cluster_data/merged_clusters_chr", ch, ".Rds"))
}


save(gene_coord, gene_names_ids, gene_summary, CLIP_datasets, CLIPdata_paths, file = "iCLIP_shiny_app/CLIP_tool_data.RData")
