#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Combine two A594 cfChromatin replicates by averaging RPKM scores  ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-04-24

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")
source("scripts/plot_theme_mres_frankell.R")

# -----------------------------------------------------------------------------
# Load libraries
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
library(readxl)
library(car)
library(svglite)

# -----------------------------------------------------------------------------
# Make output folder
# -----------------------------------------------------------------------------

date <- gsub("-", "", Sys.Date())

analysis_name <- 'ATAC_seq_shedding_analysis_abbosh'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep = '/')
if (!file.exists(out_dir_general)) dir.create(out_dir_general)

out_dir_logs <- paste(out_dir_general, 'logs', sep = '/')
if (!file.exists(out_dir_logs)) dir.create(out_dir_logs)

outputs.folder <- paste0(out_dir_general, "/", date, "/")
if (!file.exists(outputs.folder)) dir.create(outputs.folder)

# -----------------------------------------------------------------------------
# Load ctDNA data
# -----------------------------------------------------------------------------

ctDNA_data_abbosh <- read_fst("data/ctDNA_data_abbosh_pos_multiple.fst")

# Parse chromosome and position from pure_mutation_id (format: chr:pos:alt)
ctDNA_data_abbosh <- ctDNA_data_abbosh %>%
  mutate(
    chromosome = sub(":.*", "", pure_mutation_id),
    position   = as.numeric(sub("^[^:]+:([^:]+):.*", "\\1", pure_mutation_id))
  )

# -----------------------------------------------------------------------------
# Load and process clinical data — Abbosh et al. 2023, Table_8_patient_data
# -----------------------------------------------------------------------------

clinical_abbosh <- read_excel(
  path = "data/41586_2023_5776_MOESM4_ESM.xlsx",
  sheet = "Table_8_patient_data",
  skip = 22
)

# Recode histology into LUAD, LUSC, Other
clinical_abbosh <- clinical_abbosh %>%
  mutate(histology_group = case_when(
    grepl("adenocarcinoma", Histology, ignore.case = TRUE) ~ "LUAD",
    grepl("squamous",       Histology, ignore.case = TRUE) ~ "LUSC",
    TRUE                                                    ~ "Other"
  ))

cat("Histology breakdown:\n")
print(count(clinical_abbosh, histology_group))

# Join histology onto ctDNA data via PublicationID
ctDNA_data_abbosh <- ctDNA_data_abbosh %>%
  left_join(clinical_abbosh %>% select(PublicationID, Histology, histology_group),
            by = "PublicationID")

cat("Histology in ctDNA data:\n")
ctDNA_data_abbosh %>%
  distinct(PublicationID, histology_group) %>%
  count(histology_group) %>%
  print()

# -----------------------------------------------------------------------------
# Load ATAC-seq peak files
# (same files as Black et al. analysis — hg19 liftover, chr prefix stripped)
# LUAD/LUSC: Corces et al. 2018 (TCGA ATAC-seq)
# Blood: Corces et al. 2016 (GSE74912)
# -----------------------------------------------------------------------------

# LUAD hg19 — Corces et al. 2018
luad_peaks_hg19 <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUAD_peakCalls_hg19.bed",
  header = FALSE, sep = "\t",
  col.names = c("chr", "start", "end", "hg38_coords", "score")
)
luad_peaks_hg19$chr <- gsub("chr", "", luad_peaks_hg19$chr)
luad_peaks_hg19 <- subset(luad_peaks_hg19, select = -c(score))

luad_peaks_original <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUAD_peakCalls.txt",
  header = TRUE, sep = "\t"
)

# LUSC hg19 — Corces et al. 2018
lusc_peaks_hg19 <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUSC_peakCalls_hg19.bed",
  header = FALSE, sep = "\t",
  col.names = c("chr", "start", "end", "hg38_coords", "score")
)
lusc_peaks_hg19$chr <- gsub("chr", "", lusc_peaks_hg19$chr)
lusc_peaks_hg19 <- subset(lusc_peaks_hg19, select = -c(score))

lusc_peaks_original <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUSC_peakCalls.txt",
  header = TRUE, sep = "\t"
)

# Blood ATAC-seq — Corces et al. 2016 (GSE74912)
blood_counts <- fread("data/GSE74912_ATACseq_All_Counts.txt.gz")

# -----------------------------------------------------------------------------
# Add ATAC-seq score from hg38 to hg19 lifted-over peaks
# -----------------------------------------------------------------------------

luad_peaks_original$hg38_coords <- paste0(luad_peaks_original$seqnames, ":",
                                          luad_peaks_original$start + 1, "-",
                                          luad_peaks_original$end)
luad_peaks <- left_join(luad_peaks_hg19, luad_peaks_original, by = "hg38_coords") %>%
  select(chr, start = start.x, end = end.x, score, name)

lusc_peaks_original$hg38_coords <- paste0(lusc_peaks_original$seqnames, ":",
                                          lusc_peaks_original$start + 1, "-",
                                          lusc_peaks_original$end)
lusc_peaks <- left_join(lusc_peaks_hg19, lusc_peaks_original, by = "hg38_coords") %>%
  select(chr, start = start.x, end = end.x, score, name)

# -----------------------------------------------------------------------------
# Process blood ATAC-seq — Corces et al. 2016 (GSE74912)
# -----------------------------------------------------------------------------

blood_peaks <- blood_counts[, .(Chr, Start, End)]
blood_peaks[, peak_id := paste0(Chr, "_", Start, "_", End)]
blood_peaks[, Chr := gsub("chr", "", Chr)]

counts <- as.matrix(blood_counts[, -(1:3)])
rownames(counts) <- blood_peaks$peak_id

# Remove leukaemic or pre-leukaemic samples
normal_cols <- !grepl("Leuk|LSC|pHSC", colnames(counts))
counts <- counts[, normal_cols]

# Normalise for library size
lib_sizes  <- colSums(counts)
log_cpm    <- log2(t(t(counts) / lib_sizes) * 1e6 + 1)
blood_mean <- rowMeans(log_cpm)
blood_peaks[, blood_mean := blood_mean]

# -----------------------------------------------------------------------------
# Convert peaks to GRanges objects
# -----------------------------------------------------------------------------

luad_peaks_gr <- GRanges(
  seqnames = luad_peaks$chr,
  ranges   = IRanges(start = luad_peaks$start, end = luad_peaks$end),
  score    = luad_peaks$score
)

lusc_peaks_gr <- GRanges(
  seqnames = lusc_peaks$chr,
  ranges   = IRanges(start = lusc_peaks$start, end = lusc_peaks$end),
  score    = lusc_peaks$score
)

blood_peaks_gr <- GRanges(
  seqnames = blood_peaks$Chr,
  ranges   = IRanges(start = blood_peaks$Start, end = blood_peaks$End),
  score    = blood_peaks$blood_mean
)

# -----------------------------------------------------------------------------
# Annotate helper function
# -----------------------------------------------------------------------------

annotate_peaks <- function(ctDNA_df, peaks_gr, peak_col, score_col) {
  mutations_gr <- GRanges(
    seqnames = ctDNA_df$chromosome,
    ranges   = IRanges(start = ctDNA_df$position, end = ctDNA_df$position)
  )
  overlaps <- findOverlaps(mutations_gr, peaks_gr)
  
  ctDNA_df[[peak_col]]  <- FALSE
  ctDNA_df[[score_col]] <- NA
  ctDNA_df[[peak_col]][queryHits(overlaps)]  <- TRUE
  ctDNA_df[[score_col]][queryHits(overlaps)] <- peaks_gr$score[subjectHits(overlaps)]
  
  cat(peak_col, "- in peak:", sum(ctDNA_df[[peak_col]]),
      "| outside:", sum(!ctDNA_df[[peak_col]]), "\n")
  
  return(ctDNA_df)
}

# -----------------------------------------------------------------------------
# Annotate full cohort with blood peaks (Corces et al. 2016)
# Blood chromatin is not histology-specific so run on all patients together
# -----------------------------------------------------------------------------

ctDNA_data_abbosh <- annotate_peaks(ctDNA_data_abbosh, blood_peaks_gr,
                                    "in_blood_peak", "blood_atac_score")

# -----------------------------------------------------------------------------
# Split ctDNA data by histology for tumour ATAC analyses
# -----------------------------------------------------------------------------

ctDNA_LUAD_abbosh <- ctDNA_data_abbosh %>% filter(histology_group == "LUAD")
ctDNA_LUSC_abbosh <- ctDNA_data_abbosh %>% filter(histology_group == "LUSC")

cat("LUAD mutations:", nrow(ctDNA_LUAD_abbosh), "\n")
cat("LUSC mutations:", nrow(ctDNA_LUSC_abbosh), "\n")

# Annotate LUAD with LUAD tumour peaks (Corces et al. 2018)
ctDNA_LUAD_abbosh <- annotate_peaks(ctDNA_LUAD_abbosh, luad_peaks_gr,
                                    "in_luad_peak", "luad_atac_score")

# Annotate LUSC with LUSC tumour peaks (Corces et al. 2018)
ctDNA_LUSC_abbosh <- annotate_peaks(ctDNA_LUSC_abbosh, lusc_peaks_gr,
                                    "in_lusc_peak", "lusc_atac_score")

# -----------------------------------------------------------------------------
# Assign peak categories: neither / lung_only / blood_only / both
# -----------------------------------------------------------------------------

ctDNA_LUAD_abbosh <- ctDNA_LUAD_abbosh %>%
  mutate(peak_category = case_when(
    in_blood_peak &  in_luad_peak  ~ "both",
    in_blood_peak & !in_luad_peak  ~ "blood_only",
    !in_blood_peak &  in_luad_peak ~ "lung_only",
    TRUE                           ~ "neither"
  ))

ctDNA_LUSC_abbosh <- ctDNA_LUSC_abbosh %>%
  mutate(peak_category = case_when(
    in_blood_peak &  in_lusc_peak  ~ "both",
    in_blood_peak & !in_lusc_peak  ~ "blood_only",
    !in_blood_peak &  in_lusc_peak ~ "lung_only",
    TRUE                           ~ "neither"
  ))

# -----------------------------------------------------------------------------
# Plotting functions
# -----------------------------------------------------------------------------

# Function 1: correlation — ATAC score vs CCF z-score within peaks
plot_atac_correlation <- function(ctDNA_df, score_col, peak_col, cancer_type, save_path) {
  df <- ctDNA_df[ctDNA_df[[peak_col]] == TRUE, ]
  
  p <- ggplot(df, aes(x = .data[[score_col]], y = ccf_z_score)) +
    geom_point(alpha = 0.3) +
    stat_cor(method = "spearman",
             label.x = 0.5,
             label.y = max(df$ccf_z_score, na.rm = TRUE) * 0.9) +
    geom_smooth(method = "lm", colour = "steelblue", fill = "lightblue", alpha = 0.3) +
    labs(
      x     = paste(cancer_type, "ATAC-seq peak score (chromatin accessibility)"),
      y     = "CCF z-score (ctDNA shedding)",
      title = paste("Chromatin accessibility vs ctDNA shedding —",
                    cancer_type, "— Abbosh et al. 2023")
    ) +
    theme_cowplot() +
    theme(plot.title = element_text(size = 10))
  
  ggsave(paste0(save_path, "atac_correlation_", cancer_type, "_abbosh.pdf"),
         p, width = 6, height = 6)
  
  cor_result <- cor.test(df[[score_col]], df$ccf_z_score, method = "spearman")
  print(cor_result)
}

# Function 2: boxplot — in vs outside peaks
plot_atac_boxplot <- function(ctDNA_df, peak_col, cancer_type, save_path) {
  p <- ggplot(ctDNA_df, aes(x = .data[[peak_col]], y = ccf_z_score,
                            fill = .data[[peak_col]])) +
    geom_boxplot(outlier.alpha = 0.2) +
    scale_fill_manual(values = c("FALSE" = low_col, "TRUE" = high_col)) +
    stat_compare_means(method = "wilcox.test", label = "p.format", label.x = 1.3) +
    labs(
      x     = paste("Mutation in", cancer_type, "ATAC-seq peak"),
      y     = "CCF z-score") +
    theme_cowplot(font_size = 20) +
    theme(legend.position = "none")
  
  ggsave(paste0(save_path, "atac_boxplot_", cancer_type, "_abbosh.pdf"),
         p, width = 6, height = 6)
  ggsave(paste0(save_path, "atac_boxplot_", cancer_type, "_abbosh.svg"),
         p, width = 6, height = 6)
  
  wilcox_result <- wilcox.test(as.formula(paste("ccf_z_score ~", peak_col)),
                               data = ctDNA_df)
  print(wilcox_result)
}

# Function 3: patient-level paired analysis
analyse_patient_level <- function(ctDNA_df, peak_col, cancer_type, save_path) {
  
  patient_medians <- ctDNA_df %>%
    group_by(PublicationID, .data[[peak_col]]) %>%
    summarise(median_z = median(ccf_z_score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from   = all_of(peak_col),
                values_from  = median_z,
                names_prefix = "peak_")
  
  test_result <- wilcox.test(patient_medians$peak_TRUE,
                             patient_medians$peak_FALSE,
                             paired = TRUE)
  print(test_result)
  
  n_labels <- ctDNA_df %>%
    group_by(PublicationID, .data[[peak_col]]) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from   = all_of(peak_col),
                values_from  = n,
                names_prefix = "n_") %>%
    mutate(label = paste0("n_false=", n_FALSE, ", n_true=", n_TRUE))
  
  patient_long <- patient_medians %>%
    pivot_longer(cols = c(peak_FALSE, peak_TRUE),
                 names_to = "in_peak", values_to = "median_z") %>%
    mutate(in_peak = factor(in_peak,
                            levels = c("peak_FALSE", "peak_TRUE"),
                            labels = c("Outside peak", "In peak")))
  
  p_lines <- ggplot(patient_long,
                    aes(x = in_peak, y = median_z, group = PublicationID)) +
    geom_line(alpha = 0.5, colour = "grey40") +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_compare_means(method = "wilcox.test", paired = TRUE,
                       comparisons = list(c("Outside peak", "In peak"))) +
    labs(
      x        = paste0("Mutation location (", cancer_type, " ATAC-seq peaks)"),
      y        = "Median CCF z-score per patient",
      title    = paste("Per-patient shedding: in vs outside open chromatin —",
                       cancer_type, "— Abbosh et al. 2023"),
      subtitle = paste0("n = ", nrow(patient_medians), " patients")
    ) +
    theme_cowplot() +
    theme(plot.title    = element_text(size = 11),
          plot.subtitle = element_text(size = 9))
  
  ggsave(paste0(save_path, "atac_per_patient_lines_", cancer_type, "_abbosh.pdf"),
         p_lines, width = 6, height = 6)
  
  p_box <- ggplot(ctDNA_df,
                  aes(x = .data[[peak_col]], y = ccf_z_score,
                      fill = .data[[peak_col]])) +
    geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3) +
    geom_text(data = n_labels,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.1, vjust = 1.2, size = 2.5, inherit.aes = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    facet_wrap(~ PublicationID, scales = "free_y", ncol = 5) +
    scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "#E41A1C")) +
    labs(
      x     = paste("In", cancer_type, "ATAC-seq peak"),
      y     = "CCF z-score",
      title = paste("Per-patient CCF z-scores —", cancer_type, "— Abbosh et al. 2023")
    ) +
    theme_cowplot() +
    theme(legend.position = "none", strip.text = element_text(size = 8))
  
  ggsave(paste0(save_path, "atac_per_patient_boxplots_", cancer_type, "_abbosh.pdf"),
         p_box, width = 10, height = 10)
  
  return(list(test = test_result, patient_medians = patient_medians))
}

# Function 4: peak category boxplots (neither/lung_only/blood_only/both)
plot_peak_categories <- function(ctDNA_df, cancer_type, save_path) {
  
  df <- ctDNA_df %>%
    mutate(peak_category = factor(peak_category,
                                  levels = c("neither", "lung_only",
                                             "blood_only", "both")))
  
  p <- ggplot(df, aes(x = peak_category, y = ccf_z_score, fill = peak_category)) +
    geom_boxplot(outlier.alpha = 0.2) +
    stat_compare_means(comparisons = list(
      c("neither", "lung_only"),   c("neither", "blood_only"),
      c("neither", "both"),        c("lung_only", "blood_only"),
      c("lung_only", "both"),      c("blood_only", "both")
    ), method = "wilcox.test") +
    scale_fill_manual(values = c(
      "neither"    = "#D9D9D9",
      "lung_only"  = "#A8D5A2",
      "blood_only" = "#A8C4E0",
      "both"       = "#F4A582"
    )) +
    labs(
      x     = "Peak category",
      y     = "CCF z-score (ctDNA shedding)",
      title = paste("Shedding by chromatin accessibility context —",
                    cancer_type, "— Abbosh et al. 2023")
    ) +
    theme_cowplot() +
    theme(legend.position = "none")
  
  ggsave(paste0(save_path, "peak_category_boxplot_", cancer_type, "_abbosh.pdf"),
         p, width = 7, height = 6)
  
  print(count(df, peak_category))
  print(kruskal.test(ccf_z_score ~ peak_category, data = df))
}

# Function 5: patient-level peak category analysis
analyse_peak_categories_patient_level <- function(ctDNA_df, cancer_type, save_path) {
  
  patient_category_medians <- ctDNA_df %>%
    group_by(PublicationID, peak_category) %>%
    summarise(
      median_z    = median(ccf_z_score, na.rm = TRUE),
      n_mutations = n(),
      se          = sd(ccf_z_score, na.rm = TRUE) / sqrt(n()),
      ci_low      = median_z - 1.96 * se,
      ci_high     = median_z + 1.96 * se,
      .groups     = "drop"
    ) %>%
    mutate(peak_category = factor(peak_category,
                                  levels = c("neither", "lung_only",
                                             "blood_only", "both")))
  
  kruskal_result <- kruskal.test(median_z ~ peak_category,
                                 data = patient_category_medians)
  levene_result  <- leveneTest(median_z ~ peak_category,
                               data = patient_category_medians)
  
  combined_label <- paste0("Levene's test p = ",
                           signif(levene_result$`Pr(>F)`[1], 3),
                           "\nKruskal-Wallis p = ",
                           signif(kruskal_result$p.value, 3))
  
  p <- ggplot(patient_category_medians,
              aes(x = peak_category, y = median_z, fill = peak_category)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    geom_pointrange(aes(ymin = ci_low, ymax = ci_high),
                    position = position_jitter(width = 0.25, seed = 42),
                    size = 0.3, alpha = 0.5) +
    annotate("text", x = 1.2,
             y = max(patient_category_medians$median_z, na.rm = TRUE) * 1.7,
             label = combined_label, size = 4) +
    scale_fill_manual(values = c(
      "neither"    = "#D9D9D9",
      "lung_only"  = "#A8D5A2",
      "blood_only" = "#A8C4E0",
      "both"       = "#F4A582"
    )) +
    labs(
      x        = "Peak category",
      y        = "Median CCF z-score per patient",
      title    = paste("Patient-level shedding by chromatin accessibility context —",
                       cancer_type, "— Abbosh et al. 2023"),
      subtitle = "One median per patient per category"
    ) +
    theme_cowplot() +
    theme(legend.position = "none", plot.title = element_text(size = 12))
  
  ggsave(paste0(save_path, "peak_category_patient_level_", cancer_type, "_abbosh.pdf"),
         p, width = 7, height = 6)
  
  print(kruskal.test(median_z ~ peak_category, data = patient_category_medians))
  print(pairwise.wilcox.test(patient_category_medians$median_z,
                             patient_category_medians$peak_category,
                             p.adjust.method = "BH"))
  
  return(patient_category_medians)
}

# -----------------------------------------------------------------------------
# Run analyses — blood ATAC (Corces et al. 2016) — full cohort, all histologies
# -----------------------------------------------------------------------------

plot_atac_correlation(ctDNA_data_abbosh, "blood_atac_score", "in_blood_peak",
                      "blood", outputs.folder)
plot_atac_boxplot(ctDNA_data_abbosh, "in_blood_peak", "blood", outputs.folder)
blood_result <- analyse_patient_level(ctDNA_data_abbosh, "in_blood_peak",
                                      "blood", outputs.folder)

# -----------------------------------------------------------------------------
# Run analyses — LUAD tumour ATAC (Corces et al. 2018)
# -----------------------------------------------------------------------------

plot_atac_correlation(ctDNA_LUAD_abbosh, "luad_atac_score", "in_luad_peak",
                      "LUAD", outputs.folder)
plot_atac_boxplot(ctDNA_LUAD_abbosh, "in_luad_peak", "LUAD", outputs.folder)
luad_patient_result <- analyse_patient_level(ctDNA_LUAD_abbosh, "in_luad_peak",
                                             "LUAD", outputs.folder)
plot_peak_categories(ctDNA_LUAD_abbosh, "LUAD", outputs.folder)
luad_category_medians <- analyse_peak_categories_patient_level(ctDNA_LUAD_abbosh,
                                                               "LUAD", outputs.folder)

# -----------------------------------------------------------------------------
# Run analyses — LUSC tumour ATAC (Corces et al. 2018)
# -----------------------------------------------------------------------------

plot_atac_correlation(ctDNA_LUSC_abbosh, "lusc_atac_score", "in_lusc_peak",
                      "LUSC", outputs.folder)
plot_atac_boxplot(ctDNA_LUSC_abbosh, "in_lusc_peak", "LUSC", outputs.folder)
lusc_patient_result <- analyse_patient_level(ctDNA_LUSC_abbosh, "in_lusc_peak",
                                             "LUSC", outputs.folder)
plot_peak_categories(ctDNA_LUSC_abbosh, "LUSC", outputs.folder)
lusc_category_medians <- analyse_peak_categories_patient_level(ctDNA_LUSC_abbosh,
                                                               "LUSC", outputs.folder)

# -----------------------------------------------------------------------------
# Save annotated objects
# -----------------------------------------------------------------------------

saveRDS(ctDNA_LUAD_abbosh, "data/ctDNA_LUAD_abbosh_annotated.rds")
saveRDS(ctDNA_LUSC_abbosh, "data/ctDNA_LUSC_abbosh_annotated.rds")