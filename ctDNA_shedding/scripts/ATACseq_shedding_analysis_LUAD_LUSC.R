#==============================================================================#
#==============================================================================#
######                                                                    ######
######  ATAC-seq and shedding analysis                                    ######
######                                                                    ######
######  Goal: Test whether chromatin accessibility at mutation sites      ######
######        (TCGA LUAD/LUSC ATAC-seq, Corces et al. 2018) correlates    ######
######        with ctDNA shedding levels (CCF z-score) in TRACERx         ######
######        lung cancer patients                                        ######
######  Input:  - LUAD/LUSC ATAC-seq peak calls (hg38, Corces 2018)       ######
######          - Lifted-over peak calls (hg19)                           ######
######          - TRACERx ctDNA mutation data with CCF z-scores           ######
######          - TRACERx clinical outcomes data for cancer type          ######
######                                                                    ######
######  Output: - Correlation plots of ATAC-seq score vs CCF z-score      ######
#==============================================================================#
#==============================================================================#


# Author: Hannah Bazin
# Date: 2026-04-21

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")
source("scripts/plot_theme_mres_frankell.R")

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
library(svglite)

# -----------------------------------------------------------------------------
# Make a folder for this analysis run
# -----------------------------------------------------------------------------

date <- gsub("-","",Sys.Date())

analysis_name <- 'ATAC_seq_shedding_analysis'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )

# -----------------------------------------------------------------------------
# Load clinical data
# -----------------------------------------------------------------------------

clinical_data <- read_tsv("data/tx842_clinical_outcome_20251211.tsv")

# -----------------------------------------------------------------------------
# Load and reformat ATAC-seq data
# -----------------------------------------------------------------------------

# LUAD hg19
luad_peaks_hg19 <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUAD_peakCalls_hg19.bed",
  header = FALSE,
  sep = "\t",
  col.names = c("chr", "start", "end", "hg38_coords", "score")
)
luad_peaks_hg19$chr <- gsub("chr", "", luad_peaks_hg19$chr)
luad_peaks_hg19 <- subset(luad_peaks_hg19, select = -c(score))
head(luad_peaks_hg19)

# LUAD original
luad_peaks_original <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUAD_peakCalls.txt",
  header = TRUE,
  sep = "\t"
)

# LUSC
lusc_peaks_hg19 <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUSC_peakCalls_hg19.bed",
  header = FALSE,
  sep = "\t",
  col.names = c("chr", "start", "end", "hg38_coords", "score")
)
lusc_peaks_hg19$chr <- gsub("chr", "", lusc_peaks_hg19$chr)
lusc_peaks_hg19 <- subset(lusc_peaks_hg19, select = -c(score))
head(lusc_peaks_hg19)

# LUSC original
lusc_peaks_original <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUSC_peakCalls.txt",
  header = TRUE,
  sep = "\t"
)

# -----------------------------------------------------------------------------
# Load ctDNA data
# -----------------------------------------------------------------------------

# Read in ctDNA data with CCF z-scores
ctDNA_data_path <- "outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst"
ctDNA_data <- read_fst(ctDNA_data_path)

# Add column for cancer type
ctDNA_data <- ctDNA_data %>% left_join(clinical_data, by=c("patient" = "Shorter_ID"))
ctDNA_data %>%
  distinct(patient, histology1_group_central.reviewed) %>%
  count(histology1_group_central.reviewed)

# Filter data per cancer type
ctDNA_LUAD <- ctDNA_data %>% 
  filter(histology1_group_central.reviewed == "LUAD")

ctDNA_LUSC <- ctDNA_data %>% 
  filter(histology1_group_central.reviewed == "LUSC")

# -----------------------------------------------------------------------------
# Add ATAC-seq score from hg38 to hg19 data
# -----------------------------------------------------------------------------

# Add column matching the hg19 dataframe, with +1 to match the hg19 format (0-based vs 1-based format)
luad_peaks_original$hg38_coords <- paste0(luad_peaks_original$seqnames, ":",
                                          luad_peaks_original$start + 1, "-",
                                          luad_peaks_original$end)

# Join on hg38 coordinates
luad_peaks <- left_join(luad_peaks_hg19, luad_peaks_original, by = "hg38_coords")

luad_peaks <- luad_peaks %>%
  select(
    chr,
    start = start.x,
    end = end.x,
    score,
    name,
  )
head(luad_peaks)

# Add column matching the hg19 dataframe, with +1 to match the hg19 format (0-based vs 1-based format)
lusc_peaks_original$hg38_coords <- paste0(lusc_peaks_original$seqnames, ":",
                                          lusc_peaks_original$start + 1, "-",
                                          lusc_peaks_original$end)

# Join on hg38 coordinates
lusc_peaks <- left_join(lusc_peaks_hg19, lusc_peaks_original, by = "hg38_coords")

lusc_peaks <- lusc_peaks %>%
  select(
    chr,
    start = start.x,
    end = end.x,
    score,
    name,
  )
head(lusc_peaks)


# -----------------------------------------------------------------------------
# Find out which mutations fall within ATAC-seq peaks
# -----------------------------------------------------------------------------

# Convert mutations to GRanges object
mutations_LUAD_gr <- GRanges(
  seqnames = ctDNA_LUAD$chromosome,
  ranges = IRanges(start = ctDNA_LUAD$position, end = ctDNA_LUAD$position)
)
mutations_LUSC_gr <- GRanges(
  seqnames = ctDNA_LUSC$chromosome,
  ranges = IRanges(start = ctDNA_LUSC$position, end = ctDNA_LUSC$position)
)

# Convert LUAD peaks to GRanges object
luad_peaks_gr <- GRanges(
  seqnames = luad_peaks$chr,
  ranges = IRanges(start = luad_peaks$start, end = luad_peaks$end),
  score = luad_peaks$score
)

# Find overlaps
overlaps_LUAD <- findOverlaps(mutations_LUAD_gr, luad_peaks_gr)
overlaps_LUAD

# Annotate mutations with ATAC-seq peak information
ctDNA_LUAD$in_luad_peak <- FALSE
ctDNA_LUAD$luad_atac_score <- NA

ctDNA_LUAD$in_luad_peak[queryHits(overlaps_LUAD)] <- TRUE
ctDNA_LUAD$luad_atac_score[queryHits(overlaps_LUAD)] <- luad_peaks_gr$score[subjectHits(overlaps_LUAD)]

table(ctDNA_LUAD$in_luad_peak)
head(ctDNA_LUAD[ctDNA_LUAD$in_luad_peak == TRUE, c("chromosome", "position", "ccf_z_score", "luad_atac_score")])


# Convert LUSC peaks to GRanges object
lusc_peaks_gr <- GRanges(
  seqnames = lusc_peaks$chr,
  ranges = IRanges(start = lusc_peaks$start, end = lusc_peaks$end),
  score = lusc_peaks$score
)

# Find overlaps
overlaps_LUSC <- findOverlaps(mutations_LUSC_gr, lusc_peaks_gr)
overlaps_LUSC

# Annotate mutations with ATAC-seq peak information
ctDNA_LUSC$in_lusc_peak <- FALSE
ctDNA_LUSC$lusc_atac_score <- NA

ctDNA_LUSC$in_lusc_peak[queryHits(overlaps_LUSC)] <- TRUE
ctDNA_LUSC$lusc_atac_score[queryHits(overlaps_LUSC)] <- lusc_peaks_gr$score[subjectHits(overlaps_LUSC)]

table(ctDNA_LUSC$in_lusc_peak)
head(ctDNA_LUSC[ctDNA_LUSC$in_lusc_peak == TRUE, c("chromosome", "position", "ccf_z_score", "lusc_atac_score")])

# -----------------------------------------------------------------------------
# Functions for correlation with shedding
# -----------------------------------------------------------------------------

# Function 1: correlate ATAC-seq score with CCF z-score within peaks
### Does accessibility matter? Among mutations within a peak, does the extent of chromatin accessibility influence shedding?
plot_atac_correlation <- function(ctDNA_data, score_col, peak_col, cancer_type, save_path) {
  
  in_peak <- ctDNA_data[[peak_col]]
  
  # Plot
  p <- ggplot(ctDNA_data[in_peak == TRUE, ],
              aes(x = .data[[score_col]], y = ccf_z_score)) +
    geom_point(alpha = 0.3) +
    stat_cor(method = "spearman", label.x = 0.5, label.y = 7) +
    geom_smooth(method = "lm", colour = "steelblue", fill = "lightblue", alpha = 0.3) +
    labs(
      x = paste(cancer_type, "ATAC-seq peak score (chromatin accessibility)"),
      y = "CCF z-score (ctDNA shedding)",
      title = paste("Chromatin accessibility vs ctDNA shedding -", cancer_type)
    ) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 10))
  print(p)
  ggsave(paste0(save_path, "/atac_correlation_", cancer_type, ".pdf"), p, width = 6, height = 6)
  
  # Correlation
  cor_result <- cor.test(
    ctDNA_data[[score_col]][in_peak == TRUE],
    ctDNA_data$ccf_z_score[in_peak == TRUE],
    method = "spearman"
  )
  print(cor_result)
}

# Function 2: compare CCF z-scores in vs outside peaks
### Do mutations within a chromatin peak have different CCF z-scores than those outside of a peak?
plot_atac_boxplot <- function(ctDNA_data, peak_col, cancer_type, save_path) {
  
  p <- ggplot(ctDNA_data, aes(x = .data[[peak_col]], y = ccf_z_score, 
                              fill = .data[[peak_col]])) +
    geom_boxplot(outlier.alpha = 0.2) +
    scale_fill_manual(values = c("FALSE" = low_col, "TRUE" = high_col)) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.3) +
    labs(
      x = paste("Mutation in", cancer_type, "ATAC-seq peak"),
      y = "CCF z-score") +
    theme_cowplot(font_size = 20) +
    theme(legend.position = "none")
  print(p)
  ggsave(paste0(save_path, "/atac_boxplot_", cancer_type, ".pdf"), p, width = 6, height = 6)
  ggsave(paste0(save_path, "/atac_boxplot_", cancer_type, ".svg"), p, width = 6, height = 6)
  
  # Statistical test
  formula <- as.formula(paste("ccf_z_score ~", peak_col))
  wilcox_result <- wilcox.test(formula, data = ctDNA_data)
  print(wilcox_result)
}


# LUAD
plot_atac_correlation(ctDNA_LUAD, "luad_atac_score", "in_luad_peak", "LUAD", outputs.folder)
plot_atac_boxplot(ctDNA_LUAD, "in_luad_peak", "LUAD", outputs.folder)

# LUSC
plot_atac_correlation(ctDNA_LUSC, "lusc_atac_score", "in_lusc_peak", "LUSC", outputs.folder)
plot_atac_boxplot(ctDNA_LUSC, "in_lusc_peak", "LUSC", outputs.folder)



# -----------------------------------------------------------------------------
# Function: patient-level paired analysis of peak vs non-peak shedding
# -----------------------------------------------------------------------------
analyse_patient_level <- function(ctDNA_data, peak_col, cancer_type, save_path) {
  
  # Paired test across patients
  patient_medians <- ctDNA_data %>%
    group_by(patient, .data[[peak_col]]) %>%
    summarise(median_z = median(ccf_z_score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = all_of(peak_col),
                values_from = median_z,
                names_prefix = "peak_")
  
  test_result <- wilcox.test(patient_medians$peak_TRUE,
                             patient_medians$peak_FALSE,
                             paired = TRUE)
  print(test_result)
  
  # Paired line plot
  patient_long <- patient_medians %>%
    pivot_longer(cols = c(peak_FALSE, peak_TRUE),
                 names_to = "in_peak",
                 values_to = "median_z") %>%
    mutate(in_peak = factor(in_peak,
                            levels = c("peak_FALSE", "peak_TRUE"),
                            labels = c("Outside peak", "In peak")))
  
  p_lines <- ggplot(patient_long, aes(x = in_peak, y = median_z, group = patient)) +
    geom_line(alpha = 0.5, colour = "grey40") +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", paired = TRUE,
                       comparisons = list(c("Outside peak", "In peak"))) +
    labs(
      x = paste("Mutation location (", cancer_type, "ATAC-seq peaks)"),
      y = "Median CCF z-score per patient",
      title = paste("Per-patient shedding: in vs outside open chromatin -", cancer_type),
      subtitle = paste0("n = ", nrow(patient_medians), " patients")
    ) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9))
  
  ggsave(paste0(save_path, "/atac_per_patient_CCF_zscores_lines_", cancer_type, ".pdf"),
         p_lines, width = 6, height = 6)
  
  # Compute n mutations per patient per peak category
  n_labels <- ctDNA_data %>%
    group_by(patient, .data[[peak_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(peak_col),
                values_from = n,
                names_prefix = "n_") %>%
    mutate(label = paste0("n_false=", n_FALSE, ", n_true=", n_TRUE))
  
  # Facetted boxplots
  p_box <- ggplot(ctDNA_data,
                  aes(x = .data[[peak_col]], y = ccf_z_score, fill = .data[[peak_col]])) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    geom_text(data = n_labels,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2,
              size = 2.5, inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    facet_wrap(~ patient, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#E41A1C")) +
    labs(
      x = paste("In", cancer_type, "ATAC-seq peak"),
      y = "CCF z-score",
      title = paste("Per-patient distribution of CCF z-scores -", cancer_type)
    ) +
    theme_cowplot() +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))
  
  ggsave(paste0(save_path, "/atac_per_patient_CCF_zscores_boxplots_", cancer_type, ".pdf"),
         p_box, width = 10, height = 10)
  
  return(list(test = test_result, patient_medians = patient_medians))
}

luad_result <- analyse_patient_level(ctDNA_LUAD, "in_luad_peak", "LUAD", outputs.folder)
lusc_result <- analyse_patient_level(ctDNA_LUSC, "in_lusc_peak", "LUSC", outputs.folder)

# -----------------------------------------------------------------------------
# Save LUAD and LUSC mutation info
# -----------------------------------------------------------------------------

saveRDS(ctDNA_LUAD, "data/ctDNA_LUAD_annotated.rds")
saveRDS(ctDNA_LUSC, "data/ctDNA_LUSC_annotated.rds")








