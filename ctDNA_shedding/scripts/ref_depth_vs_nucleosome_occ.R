#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Explore correlation between depth of reference allele and         ######
######  nucleosome peaks in blood                                         ######
######                                                                    ######
#==============================================================================#
#==============================================================================#


# Author: Hannah Bazin
# Date: 2026-04-23

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
# Make a folder for this analysis run
# -----------------------------------------------------------------------------

date <- gsub("-","",Sys.Date())

analysis_name <- 'nucleosome_mapping_continuous'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )

# -----------------------------------------------------------------------------
# Load required data
# -----------------------------------------------------------------------------

# Read in ctDNA data with CCF z-scores - Black 2025
ctDNA_data_path <- "outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst"
ctDNA_data <- read_fst(ctDNA_data_path)

# Processed ATAC-seq peak scores
blood_atac_peaks <- fread("data/blood_ATACseq_peaks_Corces_2016.csv")

# Read in nucleosome maps
cat("Reading nucleosome maps...\n")

read_nuc_map <- function(path) {
  dt <- data.table::fread(path)
  colnames(dt) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")
  setkey(dt, chr, start, end)
  return(dt)
}

nuc_healthy_25yo <- read_nuc_map("data/GSE114511_25yo_Teo_cfDNA_stable_100bp.bed")
nuc_healthy_70yo <- read_nuc_map("data/GSE114511_70yo_Teo_cfDNA_stable_100bp.bed")
nuc_bcell        <- read_nuc_map("data/GSE36979_Gaffney2012_Bcells_MNase-seq_stable_100bp_hg19.bed")

# -----------------------------------------------------------------------------
# Process ctDNA data
# -----------------------------------------------------------------------------

# Add chr prefix to match nucleosome map format
ctDNA_data$chromosome_chr <- paste0("chr", ctDNA_data$chromosome)

# Normalise for sample depth
ctDNA_data <- ctDNA_data %>%
  group_by(sample) %>%
  mutate(REF_norm = REF / mean(REF, na.rm = TRUE)) %>%
  ungroup()

# Make GRanges objects
mut_gr <- GRanges(seqnames = ctDNA_data$chromosome_chr,
                  ranges = IRanges(start = ctDNA_data$position, 
                                   width = 1))

# -----------------------------------------------------------------------------
# REF depth vs nucleosome occupancy
# -----------------------------------------------------------------------------

plot_ref_vs_nuc <- function(ctDNA_data, mut_gr, nuc_map, map_label) {
  
  nuc_gr <- GRanges(seqnames = nuc_map$chr,
                    ranges = IRanges(start = nuc_map$start, end = nuc_map$end),
                    nuc_occupancy = nuc_map$nuc_occupancy)
  
  hits <- findOverlaps(mut_gr, nuc_gr)
  
  ctDNA_data$nuc_occupancy <- NA_real_
  ctDNA_data$nuc_occupancy[queryHits(hits)] <- nuc_gr$nuc_occupancy[subjectHits(hits)]
  
  ctDNA_matched <- ctDNA_data %>% filter(!is.na(nuc_occupancy))
  
  p <- ggplot(ctDNA_matched, aes(x = nuc_occupancy, y = REF_norm)) +
    geom_point(alpha = 0.1, size = 0.8) +
    geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.8, se = TRUE) +
    stat_cor(method = "spearman", colour = "firebrick") +
    labs(x = "Normalised nucleosome occupancy",
         y = "Normalised REF depth",
         title = map_label,
         subtitle = paste0("n = ", nrow(ctDNA_matched), " mutations")) +
    theme_cowplot(font_size = 14)
  
  return(p)
}

# Call on each dataset
p_70yo        <- plot_ref_vs_nuc(ctDNA_data, mut_gr, nuc_healthy_70yo, "70yo cfDNA (Teo et al. 2018)")
p_25yo        <- plot_ref_vs_nuc(ctDNA_data, mut_gr, nuc_healthy_25yo, "25yo cfDNA (Teo et al. 2018)")
p_bcell       <- plot_ref_vs_nuc(ctDNA_data, mut_gr, nuc_bcell,        "B cell MNase-seq (Gaffney et al. 2012)")

# Combined multi-panel figure
nuc_combined <- plot_grid(
  p_70yo, p_25yo,
  p_bcell, NULL,
  ncol = 2,
  labels = c("A", "B", "C", ""),
  label_size = 16
)

ggsave(paste0(outputs.folder, "REF_norm_vs_nucleosome_maps.pdf"), nuc_combined, width = 10, height = 8)
ggsave(paste0(outputs.folder, "REF_norm_nucleosome_maps.svg"), nuc_combined, width = 10, height = 8)


# -----------------------------------------------------------------------------
# REF depth vs blood ATAC-seq score
# -----------------------------------------------------------------------------

# Fix chromosome naming
blood_atac_peaks[, Chr := paste0("chr", Chr)]

atac_gr <- GRanges(seqnames = blood_atac_peaks$Chr,
                   ranges = IRanges(start = blood_atac_peaks$Start,
                                    end = blood_atac_peaks$End),
                   atac_score = blood_atac_peaks$blood_mean)

hits_atac <- findOverlaps(mut_gr, atac_gr)

ctDNA_data$atac_score <- NA_real_
ctDNA_data$atac_score[queryHits(hits_atac)] <- atac_gr$atac_score[subjectHits(hits_atac)]

ctDNA_matched_atac <- ctDNA_data %>% filter(!is.na(atac_score))

p_atac <- ggplot(ctDNA_matched_atac, aes(x = atac_score, y = REF_norm)) +
  geom_point(alpha = 0.2, size = 0.8) +
  geom_smooth(method = "lm", colour = "firebrick", linewidth = 0.8, se = TRUE) +
  stat_cor(method = "spearman", colour = "firebrick") +
  labs(x = "Blood ATAC-seq peak score (Corces et al. 2016)",
       y = "Normalised REF depth",
       subtitle = paste0("n = ", nrow(ctDNA_matched_atac), " mutations")) +
  theme_cowplot(font_size = 11)

ggsave(paste0(outputs.folder, "ATAC_REF_norm_vs_atac_score.pdf"), p_atac, width = 6, height = 5)
ggsave(paste0(outputs.folder, "ATAC_REF_norm_vs_atac_score.svg"), p_atac, width = 6, height = 5)

p_atac










