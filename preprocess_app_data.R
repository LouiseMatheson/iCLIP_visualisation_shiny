
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)


# As a test, GTF file and CLIP dataset subsets covering a 100kb region of chromosome 17 can be used - these can be found here:
# https://github.com/LouiseMatheson/Process_CLIP_data/tree/main/test_data
# The region is chr17:35150000-35250000 (GRCm38) and encompasses several genes, including Tnf and Lta over which we expect to see crosslinks
import("Mouse_GRCm38.90_chr17_100kb.gtf") %>%
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


### CLIP data - crosslinks

read_tsv("iCLIP_rep1_chr17_100kb.txt") %>%
  select(1:3,6,8) %>%
  rename(chromosome = chrom, crosslinks = score) %>%
  mutate(chromosome = sub("chr", "", chromosome)) %>%
  mutate(chromosome = sub("^M$", "MT", chromosome)) %>%
  filter(chromosome %in% gtf_filtered_converted$chromosome) -> iCLIP_rep1_data

read_tsv("iCLIP_rep2_chr17_100kb.txt") %>%
  select(1:3,6,8) %>%
  rename(chromosome = chrom, crosslinks = score) %>%
  mutate(chromosome = sub("chr", "", chromosome)) %>%
  mutate(chromosome = sub("^M$", "MT", chromosome)) %>%
  filter(chromosome %in% gtf_filtered_converted$chromosome) -> iCLIP_rep2_data

read_tsv("iCLIP_rep3_chr17_100kb.txt") %>%
  select(1:3,6,8) %>%
  rename(chromosome = chrom, crosslinks = score) %>%
  mutate(chromosome = sub("chr", "", chromosome)) %>%
  mutate(chromosome = sub("^M$", "MT", chromosome)) %>%
  filter(chromosome %in% gtf_filtered_converted$chromosome) -> iCLIP_rep3_data

iCLIP_rep1_data %>%
  add_column(Sample = "Replicate_1") %>%
  bind_rows(add_column(iCLIP_rep2_data, Sample = "Replicate_2")) %>%
  bind_rows(add_column(iCLIP_rep3_data, Sample = "Replicate_3")) -> all_CLIP_data

for(c in unique(all_CLIP_data$chromosome)) {
  all_CLIP_data %>%
    filter(chromosome == c) %>%
    saveRDS(paste0("iCLIP_shiny_app/CLIP_data/merged_CLIP_chr", c, ".Rds"))
}


all_CLIP_data %>%
  distinct(Sample) %>%
  add_column(label = c("iCLIP replicate 1", "iCLIP replicate 2", "iCLIP replicate 3")) %>%
  mutate(dataset = case_when(
    grepl("^CD8.*WT", Sample) ~ "CD8 T cell iCLIP, WT",
    grepl("^CD8.*KO", Sample) ~ "CD8 T cell iCLIP, KO",
    grepl("Darnell_4h", Sample) ~ "Darnell CD4 T cell HITS-CLIP, 4h activation",
    grepl("Darnell_72h", Sample) ~ "Darnell CD4 T cell HITS-CLIP, 72h + reactivation",
    grepl("^m6A_eCLIP", Sample) ~ "iGB m6A eCLIP",
    grepl("^m6A.*input", Sample) ~ "iGB m6A input",
    grepl("^Replicate", Sample) ~ "iCLIP test data",
    TRUE ~ label
  )) %>%
  mutate(Type = if_else(grepl("m6A", dataset), "m6A eCLIP", "ZFP36-family CLIP")) %>%
  mutate(replicates = case_when(
    grepl("^CD8", Sample) ~ T,
    grepl("Darnell_4h", Sample) ~ T,
    grepl("Darnell_72h", Sample) ~ T,
    grepl("m6A_eCLIP", Sample) ~ T,
    grepl("m6A_input", Sample) ~ T,
    grepl("^Replicate", Sample) ~ T,
    TRUE ~ F
  ))  -> CLIP_datasets


## iCLIP clusters

# For this, use the iCLIP_clusters_rep[1-3].txt file in test_data

read_tsv("test_data/iCLIP_clusters_rep1.txt") %>%
  add_column(dataset = "Replicate_1") %>%
  bind_rows(add_column(read_tsv("test_data/iCLIP_clusters_rep2.txt"), dataset = "Replicate_2")) %>%
  bind_rows(add_column(read_tsv("test_data/iCLIP_clusters_rep3.txt"), dataset = "Replicate_3")) %>%
  mutate(chromosome = sub("chr", "", chromosome)) %>%
  mutate(chromosome = sub("^M$", "MT", chromosome)) %>%
  filter(chromosome %in% all_CLIP_data$chromosome) %>%
  select(1:3, 6, 7) -> all_clusters

# Make cluster track for merged replicates - requiring intersect of at least 2 replicates
with(filter(all_clusters, dataset == "Replicate_1"), GRanges(chromosome, IRanges(start, end), strand = strand)) -> Rep1_clusterGR
with(filter(all_clusters, dataset == "Replicate_2"), GRanges(chromosome, IRanges(start, end), strand = strand)) -> Rep2_clusterGR
with(filter(all_clusters, dataset == "Replicate_3"), GRanges(chromosome, IRanges(start, end), strand = strand)) -> Rep3_clusterGR

c(IRanges::intersect(Rep1_clusterGR, Rep2_clusterGR), IRanges::intersect(Rep1_clusterGR, Rep3_clusterGR), IRanges::intersect(Rep2_clusterGR, Rep3_clusterGR)
) %>%
  GenomicRanges::reduce() %>%
  as_tibble() %>%
  rename(chromosome = seqnames) %>%
  select(-width) %>%
  add_column(dataset = "iCLIP test data") %>%
  bind_rows(all_clusters) -> all_clusters

# When samples are merged, we need to be selecting based on dataset names from CLIP_datasets, whereas 
# when not merged we need to use sample names. Need separate columns to define these so that both options 
# are there for those that don't have multiple replicates (otherwise you just don't get clusters for those 
# displayed when replicates are merged)

all_clusters %>%
  mutate(sample = if_else(dataset %in% CLIP_datasets$Sample, dataset, NA_character_)) %>%
  mutate(dataset = if_else(dataset %in% CLIP_datasets$dataset, dataset, CLIP_datasets$dataset[match(dataset, CLIP_datasets$Sample)])) -> all_clusters


# But this isn't quite right because for datasets with replicates, we don't want to include the clusters that come from the individual replicates when they are merged.

all_clusters %>% 
  mutate(dataset = if_else(sample %in% CLIP_datasets$Sample[CLIP_datasets$replicates == T], NA_character_, dataset)) -> all_clusters

for(c in unique(all_clusters$chromosome)) {
  all_clusters %>%
    filter(chromosome == c) %>%
    saveRDS(paste0("iCLIP_shiny_app/CLIP_data/merged_clusters_chr", c, ".Rds"))
}

save(gene_coord, gene_names_ids, gene_summary, CLIP_datasets, file = "iCLIP_shiny_app/CLIP_tool_data.RData")
