#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Reformat mutation table for fragment analysis on UCL cluster      ######
######   - Add CRUK ID                                                    ######
######   - Keep only relevant columns                                     ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-05-06

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")

# -----------------------------------------------------------------------------
# Source required functions & load libraries
# -----------------------------------------------------------------------------

# suppress warning on R build version #
library(fst)
library(data.table) 
library(dplyr) 
library(ggplot2) 
library(cowplot) 
library(RColorBrewer) 
library(tidyr)
library(GenomicRanges)
library(ggpubr)
library(readr)

# -----------------------------------------------------------------------------
# Load required data
# -----------------------------------------------------------------------------

# Read in ctDNA data with CCF z-scores - Black 2025
ctDNA_data <- read_fst("outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst")

# -----------------------------------------------------------------------------
# Create BED file for liftOver
# -----------------------------------------------------------------------------

# Write unique positions only for liftOver
liftover_input <- ctDNA_data %>%
  mutate(
    chrom = paste0("chr", chromosome),
    start = position - 1,
    end   = position,
    hg19_coords = paste(chrom, start, end, sep = "_")
  ) %>%
  select(chrom, start, end, hg19_coords) %>%
  distinct()

write_tsv(
  liftover_input,
  file      = "data/ctDNA_pos_multiple_black_mutations_hg19_unique.bed",
  col_names = FALSE
)

# -----------------------------------------------------------------------------
# Do liftOver at https://genome.ucsc.edu/cgi-bin/hgLiftOver
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Check liftOver results
# -----------------------------------------------------------------------------

# Load liftover output
liftover_output <- read_tsv(
  "data/ctDNA_pos_multiple_black_mutations_hg38_unique.bed",
  col_names = c("chrom_hg38", "start_hg38", "end_hg38", "hg19_coords")
) %>% 
  select(chrom_hg38, start_hg38, end_hg38, hg19_coords)

# Rejoin with metadata
ctDNA_data_hg38 <- ctDNA_data %>% 
  mutate(
    chrom = paste0("chr", chromosome),
    start = position - 1,
    end   = position,
    hg19_coords = paste(chrom, start, end, sep = "_")
  ) %>%
  select(patient, sample, hg19_coords, reference, alternate, ccf_z_score) %>%
  inner_join(liftover_output, by = "hg19_coords")

# Update column names
ctDNA_data_hg38 <- ctDNA_data_hg38 %>% 
  mutate(ref = reference,
         alt = alternate,
         chrom = chrom_hg38,
         start = start_hg38,
         end = end_hg38) %>% 
  select(chrom, start, end, ref, alt, patient, sample, ccf_z_score)

# -----------------------------------------------------------------------------
# Add CRUK patient ID
# -----------------------------------------------------------------------------

# Load conversion table
tracerX_CRUK_conversion <- read_tsv("data/tracerxPublicationKey_170221.txt")

# Add column with CRUK ID
ctDNA_data_hg38_cruk <- ctDNA_data_hg38 %>% 
  mutate(SampleID = ifelse(
    grepl("LTX[0-9]{4}", patient),  # already 4 digits, don't pad
    patient,
    sub("LTX", "LTX0", patient)     # 3 digits, pad to 4
  )) %>% 
  left_join(tracerX_CRUK_conversion, by = "SampleID") %>%
  mutate(tracerx_id = SampleID, cruk_id = PublicationID) %>% 
  select(chrom, start, end, ref, alt, patient, tracerx_id, cruk_id, sample, ccf_z_score)
  
head(ctDNA_data_hg38_cruk)

# -----------------------------------------------------------------------------
# Add column tracking high/low CCF z-score
# -----------------------------------------------------------------------------

# Add column with high/low CCF z-score classification
ctDNA_data_hg38_final <- ctDNA_data_hg38_cruk %>%
  mutate(ccf_z_score_group = ifelse(ccf_z_score > 0, "high", "low"))

table(ctDNA_data_hg38_final$ccf_z_score_group)

# Save for transfer to UCL cluster
write_tsv(
  ctDNA_data_hg38_final,
  "data/ctDNA_data_black_hg38_ccfzscore_highlow_20260507.tsv"
)




















