#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Blood ATAC-seq and ctDNA shedding analysis                        ######
######                                                                    ######
######  Goal: Test whether blood chromatin accessibility (GSE74912,       ######
######        Corces et al. 2016) at mutation sites influences ctDNA      ######
######        shedding (CCF z-score) in TRACERx lung cancer patients      ######
######                                                                    ######
######  Input:  - Blood ATAC-seq counts (GSE74912, hg19)                  ######
######          - ctDNA data with CCF z-scores                            ######
######          - LUAD/LUSC-annotated ctDNA objects (ATAC script 1)       ######
######                                                                    ######
######  Output: - Blood peak vs shedding correlation/boxplots             ######
######          - Peak category analysis (blood/lung/both/neither)        ######
#==============================================================================#
#==============================================================================#


# Author: Hannah Bazin
# Date: 2026-04-22

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
# Load required data
# -----------------------------------------------------------------------------

# Raw ATAC-seq read counts that fell within each peak region
blood_counts <- fread("data/GSE74912_ATACseq_All_Counts.txt.gz")

# Read in ctDNA data with CCF z-scores
ctDNA_data_path <- "outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst"
ctDNA_data <- read_fst(ctDNA_data_path)

# Cancer-type specific ctDNA data
ctDNA_LUAD <- readRDS("data/ctDNA_LUAD_annotated.rds")
ctDNA_LUSC <- readRDS("data/ctDNA_LUSC_annotated.rds")

# Clinical data
clinical_data <- read_tsv("data/tx842_clinical_outcome_20251211.tsv")

# -----------------------------------------------------------------------------
# Process blood ATAC-seq data
# -----------------------------------------------------------------------------

# Separate coordinates from counts
blood_peaks <- blood_counts[, .(Chr, Start, End)]
blood_peaks[, peak_id := paste0(Chr, "_", Start, "_", End)]
blood_peaks[, Chr := gsub("chr", "", Chr)]

counts <- as.matrix(blood_counts[, -(1:3)])
rownames(counts) <- blood_peaks$peak_id

# Remove leukaemic or pre-leukaemic samples
normal_cols <- !grepl("Leuk|LSC|pHSC", colnames(counts))
counts <- counts[, normal_cols]

# Normalise for library size
lib_sizes <- colSums(counts)
log_cpm <- log2(t(t(counts) / lib_sizes) * 1e6 + 1)
blood_mean <- rowMeans(log_cpm)
summary(blood_mean)

# Add blood_mean score to peaks
blood_peaks[, blood_mean := blood_mean]

# -----------------------------------------------------------------------------
# Which mutations fall in a blood chromatin peak?
# -----------------------------------------------------------------------------

# Convert blood peaks to GRanges
blood_peaks_gr <- GRanges(
  seqnames = blood_peaks$Chr,
  ranges = IRanges(start = blood_peaks$Start, end = blood_peaks$End),
  score = blood_peaks$blood_mean
)

# Convert mutations to GRanges
mutations_gr <- GRanges(
  seqnames = ctDNA_data$chromosome,
  ranges = IRanges(start = ctDNA_data$position, end = ctDNA_data$position)
)

# Find overlaps
overlaps_blood <- findOverlaps(mutations_gr, blood_peaks_gr)

# Annotate
ctDNA_data$in_blood_peak <- FALSE
ctDNA_data$blood_atac_score <- NA

ctDNA_data$in_blood_peak[queryHits(overlaps_blood)] <- TRUE
ctDNA_data$blood_atac_score[queryHits(overlaps_blood)] <- blood_peaks_gr$score[subjectHits(overlaps_blood)]

table(ctDNA_data$in_blood_peak)

# Add blood peak annotations to the cancer-type specific objects
blood_annotations <- ctDNA_data %>%
  distinct(patient, chromosome, position, in_blood_peak, blood_atac_score)

ctDNA_LUAD <- ctDNA_LUAD %>%
  left_join(blood_annotations, by = c("patient", "chromosome", "position"))

ctDNA_LUSC <- ctDNA_LUSC %>%
  left_join(blood_annotations, by = c("patient", "chromosome", "position"))

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

# Blood peaks and shedding
plot_atac_correlation(ctDNA_data, "blood_atac_score", "in_blood_peak", "blood", outputs.folder)
plot_atac_boxplot(ctDNA_data, "in_blood_peak", "blood", outputs.folder)

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
  
  # Compute n mutations per patient per peak category
  n_labels <- ctDNA_data %>%
    group_by(patient, .data[[peak_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = all_of(peak_col),
                values_from = n,
                names_prefix = "n_") %>%
    mutate(label = paste0("n_false=", n_FALSE, ", n_true=", n_TRUE))
  
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
      x = paste0("Mutation location (", cancer_type, " ATAC-seq peaks)"),
      y = "Median CCF z-score per patient",
      title = paste("Per-patient shedding: in vs outside open chromatin -", cancer_type),
      subtitle = paste0("n = ", nrow(patient_medians), " patients")
    ) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 11),
          plot.subtitle = element_text(size = 9))
  
  ggsave(paste0(save_path, "/atac_per_patient_CCF_zscores_lines_", cancer_type, ".pdf"),
         p_lines, width = 6, height = 6)
  
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

blood_result <- analyse_patient_level(ctDNA_data, "in_blood_peak", "blood", outputs.folder)


# -----------------------------------------------------------------------------
# Do mutations in peak in blood and not in lung have abnormal shedding?
# -----------------------------------------------------------------------------

# Create peak category: blood only, lung only, both, neither
ctDNA_LUAD <- ctDNA_LUAD %>%
  mutate(peak_category = case_when(
    in_blood_peak & in_luad_peak   ~ "both",
    in_blood_peak & !in_luad_peak  ~ "blood_only",
    !in_blood_peak & in_luad_peak  ~ "lung_only",
    TRUE                            ~ "neither"
  ))

ctDNA_LUSC <- ctDNA_LUSC %>%
  mutate(peak_category = case_when(
    in_blood_peak & in_lusc_peak   ~ "both",
    in_blood_peak & !in_lusc_peak  ~ "blood_only",
    !in_blood_peak & in_lusc_peak  ~ "lung_only",
    TRUE                            ~ "neither"
  ))

# Plot
plot_peak_categories <- function(ctDNA_data, cancer_type, save_path) {
  
  df <- ctDNA_data %>%
    mutate(peak_category = factor(peak_category,
                                  levels = c("neither", "lung_only", "blood_only", "both")))
  
  p <- ggplot(df, aes(x = peak_category, y = ccf_z_score, fill = peak_category)) +
    geom_boxplot(outlier.alpha = 0.2) +
    stat_compare_means(comparisons = list(
      c("neither",    "lung_only"),
      c("neither",    "blood_only"),
      c("neither",    "both"),
      c("lung_only",  "blood_only"),
      c("lung_only",  "both"),
      c("blood_only", "both")
    ), method = "wilcox.test") +
    scale_fill_manual(values = c(
      "neither"   = "#D9D9D9",
      "lung_only" = "#A8D5A2",
      "blood_only" = "#A8C4E0",
      "both"      = "#F4A582"
    )) +
    labs(
      x = "Peak category",
      y = "CCF z-score (ctDNA shedding)",
      title = paste("Shedding by chromatin accessibility context -", cancer_type)
    ) +
    theme_cowplot() +
    theme(legend.position = "none")
  
  print(p)
  ggsave(paste0(save_path, "/peak_category_boxplot_", cancer_type, ".pdf"), p, width = 7, height = 6)
  
  # Count per category
  print(count(df, peak_category))
  
  # Kruskal-Wallis across all groups
  print(kruskal.test(ccf_z_score ~ peak_category, data = df))
}

plot_peak_categories(ctDNA_LUAD, "LUAD", outputs.folder)
plot_peak_categories(ctDNA_LUSC, "LUSC", outputs.folder)


# -----------------------------------------------------------------------------
# Patient-level analysis: Do mutations in peak in blood and not in lung have abnormal shedding?
# -----------------------------------------------------------------------------
analyse_peak_categories_patient_level <- function(ctDNA_data, cancer_type, save_path) {
  
  patient_category_medians <- ctDNA_data %>%
    group_by(patient, peak_category) %>%
    summarise(
      median_z = median(ccf_z_score, na.rm = TRUE),
      n_mutations = n(),
      se = sd(ccf_z_score, na.rm = TRUE) / sqrt(n()),
      ci_low  = median_z - 1.96 * se,
      ci_high = median_z + 1.96 * se,
      .groups = "drop"
    ) %>%
    mutate(peak_category = factor(peak_category,
                                  levels = c("neither", "lung_only", "blood_only", "both")))
  
  
  # Run Kruskal Wallis test (are the medians different)
  kruskal_result <- kruskal.test(median_z ~ peak_category, data = patient_category_medians)
  kruskal_p <- kruskal_result$p.value
  
  # Run Levene's test (are the variances different)
  levene_result <- leveneTest(median_z ~ peak_category, data = patient_category_medians)
  levene_p <- levene_result$`Pr(>F)`[1]
  
  combined_label <- paste0("Levene's test p = ", signif(levene_p, 3),
                          "\nKruskal-Wallis p = ", signif(kruskal_p, 3))
  
  p <- ggplot(patient_category_medians,
              aes(x = peak_category, y = median_z, fill = peak_category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                    position = position_jitter(width = 0.25, seed = 42),
                    size = 0.3, alpha = 0.5) +
    annotate("text", x = 1.2, y = max(patient_category_medians$median_z) * 1.7,
             label = combined_label, size = 4) +
    scale_fill_manual(values = c(
      "neither"    = "#D9D9D9",
      "lung_only"  = "#A8D5A2",
      "blood_only" = "#A8C4E0",
      "both"       = "#F4A582"
    )) +
    labs(
      x = "Peak category",
      y = "Median CCF z-score per patient",
      title = paste("Patient-level shedding by chromatin accessibility context -", cancer_type),
      subtitle = "One median per patient per category"
    ) +
    theme_cowplot() +
    theme(legend.position = "none",
          plot.title = element_text(size = 14))
  
  print(p)
  ggsave(paste0(save_path, "/peak_category_boxplot_patient_level_", cancer_type, ".pdf"),
         p, width = 7, height = 6)
  
  # Statistical tests
  print(kruskal.test(median_z ~ peak_category, data = patient_category_medians))
  print(pairwise.wilcox.test(patient_category_medians$median_z,
                             patient_category_medians$peak_category,
                             p.adjust.method = "BH"))
  
  return(patient_category_medians)
}

patient_category_medians_LUAD <- analyse_peak_categories_patient_level(ctDNA_LUAD, "LUAD", outputs.folder)
patient_category_medians_LUSC <-analyse_peak_categories_patient_level(ctDNA_LUSC, "LUSC", outputs.folder)


# -----------------------------------------------------------------------------
# Save the blood peaks as csv file
# -----------------------------------------------------------------------------

fwrite(blood_peaks, "data/blood_ATACseq_peaks_Corces_2016.csv")





























