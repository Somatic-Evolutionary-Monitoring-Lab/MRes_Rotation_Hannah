#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Map high/low shedding mutations to CML MNase-seq nucleosome map   ######
######  K562 cell line (ENCSR000CXQ)                                      ######
######  Analyses: single base occupancy, 80bp window, dyad distance,      ######
######  nucleosome vs linker binary comparison                            ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-05-13

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
library(rtracklayer)

# -----------------------------------------------------------------------------
# Plotting theme
# -----------------------------------------------------------------------------

source("scripts/plot_theme_mres_frankell.R")

# -----------------------------------------------------------------------------
# Make a folder for this analysis to run
# -----------------------------------------------------------------------------

date <- gsub("-", "", Sys.Date())

analysis_name <- 'nucleosome_mapping_cml_mnase_seq'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if (!file.exists(out_dir_general)) dir.create(out_dir_general)

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if (!file.exists(out_dir_logs)) dir.create(out_dir_logs)

outputs.folder <- paste0(out_dir_general, "/", date, "/")
if (!file.exists(outputs.folder)) dir.create(outputs.folder)

# -----------------------------------------------------------------------------
# Get inputs required for all analyses
# -----------------------------------------------------------------------------

cat("Reading input data...\n")

# Read in ctDNA data with CCF z-scores - BLACK ET AL. 2025
ctDNA_data_black <- read_fst("outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst")

# Read in ctDNA data with CCF z-scores - ABBOSH ET AL. 2023
ctDNA_data_abbosh <- read_fst("data/ctDNA_data_abbosh_pos_multiple.fst")

# Read in CML MNase-seq data - ENCODE ENCSR000CXQ
cat("Reading CML MNase-seq bigwig file...\n")
bw <- BigWigFile("data/ENCFF000VNN.bigWig")

# Read in CML MNase-seq nucleosome dyads - ENCODE ENCSR000CXQ, processed by https://zenodo.org/records/3820875
dyads <- fread(
  "data/K562_dyads_hg19.bed.gz",
  col.names = c("chr", "start", "end")
)

# -----------------------------------------------------------------------------
# Make mutation summary - take mean z-score per mutation per patient
# -----------------------------------------------------------------------------

mutation_summary <- ctDNA_data_black %>%
  group_by(patient, Pos) %>%
  summarise(
    mean_z = mean(ccf_z_score),
    se_z = sd(ccf_z_score) / sqrt(n()),
    ci_lower = mean_z - 1.96 * se_z,
    ci_upper = mean_z + 1.96 * se_z,
    n_samples = n(),
    sig_6samples = dplyr::first(sig_6samples),
    .groups = "drop"
  ) %>%
  tidyr::separate(Pos, into = c("chr_num", "pos", "ref", "alt"),
                  sep = ":", remove = FALSE) %>%
  mutate(
    chr = paste0("chr", chr_num),
    pos = as.numeric(pos)
  ) %>%
  select(-chr_num) %>%
  select(patient, Pos, chr, pos, ref, alt, n_samples, mean_z, se_z,
         ci_lower, ci_upper, sig_6samples)

setDT(mutation_summary)
mutation_summary[, pos_end := pos + 1L]
setkey(mutation_summary, chr, pos, pos_end)

# -----------------------------------------------------------------------------
# Single base nucleosome occupancy
# -----------------------------------------------------------------------------

mut_gr <- GRanges(
  seqnames = mutation_summary$chr,
  ranges = IRanges(
    start = mutation_summary$pos,
    end   = mutation_summary$pos
  )
)

# Look up the mutations' bigWig nucleosome occupancy score
signal <- import(bw, which = mut_gr, as = "NumericList")
mutation_summary$k562_nuc_occ <- unlist(signal)

# Correlation test
cor_singlebase <- cor.test(mutation_summary$k562_nuc_occ, mutation_summary$mean_z,
                           method = "spearman")

# Single base occupancy scatter
ggplot(mutation_summary, aes(x = k562_nuc_occ, y = mean_z)) +
  geom_point(alpha = 0.3, size = 1.5, colour = scatter_col) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = horiz_line_col) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_singlebase$estimate, 3),
                          "\np = ", format(cor_singlebase$p.value, scientific = TRUE, digits = 2)),
           size = 7) +
  labs(x = "MNase-seq nucleosome occupancy",
       y = "Mean CCF z-score") +
  theme_cowplot(font_size = 20)

ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_nuc_occ_singlebase.pdf"),
       width = 7, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_nuc_occ_singlebase.svg"),
       width = 7, height = 6)

# -----------------------------------------------------------------------------
# 80bp window nucleosome occupancy
# -----------------------------------------------------------------------------

mut_gr <- GRanges(
  seqnames = mutation_summary$chr,
  ranges = IRanges(
    start = mutation_summary$pos - 80,
    end   = mutation_summary$pos + 80
  )
)

signal <- import(bw, which = mut_gr, as = "NumericList")
mutation_summary$k562_nuc_occ_80bp <- sapply(signal, mean, na.rm = TRUE)

cor_80bp <- cor.test(mutation_summary$k562_nuc_occ_80bp, mutation_summary$mean_z,
                     method = "spearman")

# -----------------------------------------------------------------------------
# Correlating with nucleosome dyads
# -----------------------------------------------------------------------------

dyads_gr <- GRanges(
  seqnames = dyads$chr,
  ranges = IRanges(start = dyads$start, end = dyads$end)
)

mut_gr <- GRanges(
  seqnames = mutation_summary$chr,
  ranges = IRanges(start = mutation_summary$pos, end = mutation_summary$pos)
)

# For each mutation, find index of nearest dyad in dyad list
nearest_idx <- nearest(mut_gr, dyads_gr)

# Computes the bp distance between each mutation and its nearest dyad
mutation_summary$dist_to_dyad <- distance(mut_gr, dyads_gr[nearest_idx])

# Correlation test
cor_dyad <- cor.test(mutation_summary$dist_to_dyad, mutation_summary$mean_z,
                     method = "spearman")

# Dyad distance scatter
ggplot(mutation_summary %>% filter(dist_to_dyad < 900), aes(x = dist_to_dyad, y = mean_z)) +
  geom_point(alpha = 0.3, size = 1.2, colour = scatter_col) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = horiz_line_col) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_dyad$estimate, 3),
                          "\np = ", format(cor_dyad$p.value, scientific = FALSE, digits = 2)),
           size = 7) +
  labs(x = "Distance to nearest nucleosome dyad (bp)",
       y = "Mean CCF z-score") +
  theme_cowplot(font_size = 20)

ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_dyad_distance.pdf"),
       width = 7, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_dyad_distance.svg"),
       width = 7, height = 6)

# Binary comparison: z-score in nucleosome or linker
mutation_summary[, nucleosome_region := ifelse(dist_to_dyad < 80,
                                               "Nucleosome",
                                               "Linker")]

table(mutation_summary$nucleosome_region)

wilcox_result <- wilcox.test(mean_z ~ nucleosome_region, data = mutation_summary)

# Nucleosome vs Linker boxplot
ggplot(mutation_summary, aes(x = nucleosome_region, y = mean_z, fill = nucleosome_region)) +
  geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = horiz_line_col) +
  scale_fill_manual(values = c("Nucleosome" = high_col, "Linker" = low_col)) +
  annotate("text", x = 1.5, y = max(mutation_summary$mean_z) * 0.9,
           label = paste0("p = ", format(wilcox_result$p.value, scientific = FALSE, digits = 2)),
           size = 7) +
  labs(x = "", y = "Mean CCF z-score") +
  theme_cowplot(font_size = 20) +
  theme(legend.position = "none")

ggsave(paste0(outputs.folder, "ccf_zscore_cml_nucleosome_vs_linker_boxplot.pdf"),
       width = 5, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_cml_nucleosome_vs_linker_boxplot.svg"),
       width = 5, height = 6)


# =============================================================================
# Abbosh et al. 2023 analysis
# =============================================================================

# -----------------------------------------------------------------------------
# Make mutation summary - take mean z-score per mutation per patient
# -----------------------------------------------------------------------------

mutation_summary_abbosh <- ctDNA_data_abbosh %>%
  group_by(PublicationID, pure_mutation_id) %>%
  summarise(
    mean_z       = mean(ccf_z_score, na.rm = TRUE),
    se_z         = sd(ccf_z_score, na.rm = TRUE) / sqrt(n()),
    ci_lower     = mean_z - 1.96 * se_z,
    ci_upper     = mean_z + 1.96 * se_z,
    n_samples    = n(),
    sig_6samples = dplyr::first(sig_6samples),
    .groups      = "drop"
  ) %>%
  tidyr::separate(pure_mutation_id,
                  into = c("chr_num", "pos", "alt"),
                  sep = ":",
                  remove = FALSE) %>%
  mutate(
    chr = paste0("chr", chr_num),
    pos = as.numeric(pos)
  ) %>%
  select(-chr_num) %>%
  select(PublicationID, pure_mutation_id, chr, pos, alt,
         n_samples, mean_z, se_z, ci_lower, ci_upper, sig_6samples)

setDT(mutation_summary_abbosh)
mutation_summary_abbosh[, pos_end := pos + 1L]
setkey(mutation_summary_abbosh, chr, pos, pos_end)

# -----------------------------------------------------------------------------
# Single base nucleosome occupancy
# -----------------------------------------------------------------------------

mut_gr_abbosh <- GRanges(
  seqnames = mutation_summary_abbosh$chr,
  ranges = IRanges(
    start = mutation_summary_abbosh$pos,
    end   = mutation_summary_abbosh$pos
  )
)

# Look up the mutations' bigWig nucleosome occupancy score
signal_abbosh <- import(bw, which = mut_gr_abbosh, as = "NumericList")
mutation_summary_abbosh$k562_nuc_occ <- unlist(signal_abbosh)

# Correlation test
cor_singlebase_abbosh <- cor.test(mutation_summary_abbosh$k562_nuc_occ,
                                  mutation_summary_abbosh$mean_z,
                                  method = "spearman")

# Single base occupancy scatter
ggplot(mutation_summary_abbosh, aes(x = k562_nuc_occ, y = mean_z)) +
  geom_point(alpha = 0.3, size = 1.2, colour = scatter_col) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = horiz_line_col) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_singlebase_abbosh$estimate, 3),
                          "\np = ", format(cor_singlebase_abbosh$p.value, scientific = TRUE, digits = 2)),
           size = 7) +
  labs(x = "MNase-seq nucleosome occupancy",
       y = "Mean CCF z-score") +
  theme_cowplot(font_size = 20)

ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_nuc_occ_singlebase_abbosh.pdf"),
       width = 7, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_nuc_occ_singlebase_abbosh.svg"),
       width = 7, height = 6)

# -----------------------------------------------------------------------------
# 80bp window nucleosome occupancy
# -----------------------------------------------------------------------------

mut_gr_abbosh <- GRanges(
  seqnames = mutation_summary_abbosh$chr,
  ranges = IRanges(
    start = mutation_summary_abbosh$pos - 80,
    end   = mutation_summary_abbosh$pos + 80
  )
)

signal_abbosh <- import(bw, which = mut_gr_abbosh, as = "NumericList")
mutation_summary_abbosh$k562_nuc_occ_80bp <- sapply(signal_abbosh, mean, na.rm = TRUE)

cor_80bp_abbosh <- cor.test(mutation_summary_abbosh$k562_nuc_occ_80bp,
                            mutation_summary_abbosh$mean_z,
                            method = "spearman")

# -----------------------------------------------------------------------------
# Correlating with nucleosome dyads
# -----------------------------------------------------------------------------

mut_gr_abbosh <- GRanges(
  seqnames = mutation_summary_abbosh$chr,
  ranges = IRanges(start = mutation_summary_abbosh$pos,
                   end   = mutation_summary_abbosh$pos)
)

# For each mutation, find index of nearest dyad in dyad list
nearest_idx_abbosh <- nearest(mut_gr_abbosh, dyads_gr)

# Computes the bp distance between each mutation and its nearest dyad
mutation_summary_abbosh$dist_to_dyad <- distance(mut_gr_abbosh, dyads_gr[nearest_idx_abbosh])

# Correlation test
cor_dyad_abbosh <- cor.test(mutation_summary_abbosh$dist_to_dyad,
                            mutation_summary_abbosh$mean_z,
                            method = "spearman")

# Dyad distance scatter
ggplot(mutation_summary_abbosh %>% filter(dist_to_dyad < 900), aes(x = dist_to_dyad, y = mean_z)) +
  geom_point(alpha = 0.3, size = 1.2, colour = scatter_col) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = horiz_line_col) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_dyad_abbosh$estimate, 3),
                          "\np = ", format(cor_dyad_abbosh$p.value, scientific = FALSE, digits = 2)),
           size = 7) +
  labs(x = "Distance to nearest nucleosome dyad (bp)",
       y = "Mean CCF z-score") +
  theme_cowplot(font_size = 20)

ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_dyad_distance_abbosh.pdf"),
       width = 7, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_vs_cml_dyad_distance_abbosh.svg"),
       width = 7, height = 6)

# Binary comparison: z-score in nucleosome or linker
mutation_summary_abbosh[, nucleosome_region := ifelse(dist_to_dyad < 80,
                                                      "Nucleosome",
                                                      "Linker")]

table(mutation_summary_abbosh$nucleosome_region)

wilcox_result_abbosh <- wilcox.test(mean_z ~ nucleosome_region, data = mutation_summary_abbosh)

# Nucleosome vs Linker boxplot
ggplot(mutation_summary_abbosh, aes(x = nucleosome_region, y = mean_z, fill = nucleosome_region)) +
  geom_boxplot(outlier.alpha = 0.2, outlier.size = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = horiz_line_col) +
  scale_fill_manual(values = c("Nucleosome" = high_col, "Linker" = low_col)) +
  annotate("text", x = 1.5, y = max(mutation_summary_abbosh$mean_z) * 0.9,
           label = paste0("p = ", format(wilcox_result_abbosh$p.value, scientific = FALSE, digits = 2)),
           size = 7) +
  labs(x = "", y = "Mean CCF z-score") +
  theme_cowplot(font_size = 20) +
  theme(legend.position = "none")

ggsave(paste0(outputs.folder, "ccf_zscore_cml_nucleosome_vs_linker_boxplot_abbosh.pdf"),
       width = 5, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_cml_nucleosome_vs_linker_boxplot_abbosh.svg"),
       width = 5, height = 6)

