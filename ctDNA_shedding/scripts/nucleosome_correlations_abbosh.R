#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Explore CCF z-score correlations with nucleosomes                 ######
######  from Abbosh et al 2023                                            ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-04-24

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")
source("scripts/plot_theme_mres_frankell.R")

####################################################
#### Source required functions & load libraries ####
####################################################

library(fst)
library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(tidyr)
library(readxl)
library(svglite)

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-", "", Sys.Date())

analysis_name <- 'nucleosome_mapping_continuous_abbosh'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if (!file.exists(out_dir_general)) dir.create(out_dir_general)

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if (!file.exists(out_dir_logs)) dir.create(out_dir_logs)

outputs.folder <- paste0(out_dir_general, "/", date, "/")
if (!file.exists(outputs.folder)) dir.create(outputs.folder)

##############################################
#### Get Inputs required for all analyses ####
##############################################

cat("Reading input data...\n")

# Read in ctDNA data with CCF z-scores
ctDNA_data <- read_fst("data/ctDNA_data_abbosh_pos_multiple.fst")

##############################################
#### Helper functions ####
##############################################

format_p <- function(p) {
  if (p < 0.01) {
    format(p, scientific = TRUE, digits = 2)
  } else {
    format(round(p, 3), scientific = FALSE)
  }
}

# Read in nucleosome maps
cat("Reading nucleosome maps...\n")

read_nuc_map <- function(path) {
  dt <- data.table::fread(path)
  colnames(dt) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")
  setkey(dt, chr, start, end)
  return(dt)
}

nuc_lung_ma      <- read_nuc_map("data/SRA438908_lung_cancer_Ma2017_stable_100bp_hg19.bed")
nuc_lung_snyder  <- read_nuc_map("data/GSE71378_lung_cancer_Snyder_cfDNA_stable_100bp.bed")
nuc_healthy_25yo <- read_nuc_map("data/GSE114511_25yo_Teo_cfDNA_stable_100bp.bed")
nuc_healthy_70yo <- read_nuc_map("data/GSE114511_70yo_Teo_cfDNA_stable_100bp.bed")
nuc_bcell        <- read_nuc_map("data/GSE36979_Gaffney2012_Bcells_MNase-seq_stable_100bp_hg19.bed")

###############################################################
#### Build mutation summary across all patients            ####
###############################################################

cat("Building mutation summary...\n")

mutation_summary_all <- ctDNA_data %>%
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

# Convert to data.table and add pos_end for foverlaps
setDT(mutation_summary_all)
mutation_summary_all[, pos_end := pos + 1L]
setkey(mutation_summary_all, chr, pos, pos_end)

################################################################################################
#### Fast overlap function using foverlaps                                                  ####
################################################################################################

map_nuc_occupancy <- function(mutations, nuc_map, label) {
  cat("Mapping:", label, "\n")
  
  result <- foverlaps(
    mutations,
    nuc_map[, .(chr, start, end, nuc_occupancy)],
    by.x = c("chr", "pos", "pos_end"),
    by.y = c("chr", "start", "end"),
    type = "within",
    nomatch = NA
  )
  
  cat(label, "- Mutations with nucleosome occupancy:", sum(!is.na(result$nuc_occupancy)), "\n")
  cat(label, "- Mutations without (unstable regions):", sum(is.na(result$nuc_occupancy)), "\n")
  cat(label, "- Proportion covered:", round(mean(!is.na(result$nuc_occupancy)) * 100, 1), "%\n")
  
  return(result)
}

################################################################################################
#### Map all nucleosome datasets                                                            ####
################################################################################################

mut_lung_ma      <- map_nuc_occupancy(mutation_summary_all, nuc_lung_ma,      "Ma et al. 2017 lung cancer")
mut_lung_snyder  <- map_nuc_occupancy(mutation_summary_all, nuc_lung_snyder,  "Snyder et al. 2016 lung cancer")
mut_25yo         <- map_nuc_occupancy(mutation_summary_all, nuc_healthy_25yo, "Teo et al. 2018 25yo")
mut_70yo         <- map_nuc_occupancy(mutation_summary_all, nuc_healthy_70yo, "Teo et al. 2018 70yo")
mut_bcell        <- map_nuc_occupancy(mutation_summary_all, nuc_bcell,        "Gaffney et al. 2012 B cell")

################################################################################################
#### Correlation tests                                                                      ####
################################################################################################

run_cor <- function(dt) {
  cor.test(dt$mean_z, dt$nuc_occupancy, method = "spearman", alternative = "two.sided")
}

cor_lung_ma     <- run_cor(mut_lung_ma)
cor_lung_snyder <- run_cor(mut_lung_snyder)
cor_25yo        <- run_cor(mut_25yo)
cor_70yo        <- run_cor(mut_70yo)
cor_bcell       <- run_cor(mut_bcell)

################################################################################################
#### Build plots                                                                            ####
################################################################################################

make_plot <- function(dt, cor_result, colour, title) {
  dt %>%
    filter(!is.na(nuc_occupancy)) %>%
    ggplot(aes(x = nuc_occupancy, y = mean_z)) +
    geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
    geom_smooth(method = "lm", colour = colour, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("rho = ", round(cor_result$estimate, 3),
                            "\np = ", format_p(cor_result$p.value)),
             size = 5) +
    labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score", # Mean CCF z-score across timepoints
         title = title) +
    theme_cowplot(font_size = 14)
}

p_lung_ma     <- make_plot(mut_lung_ma,     cor_lung_ma,     high_col, "Lung cancer cfDNA") # Ma et al 2017
p_lung_snyder <- make_plot(mut_lung_snyder, cor_lung_snyder, high_col, "Lung cancer cfDNA") # Snyder et al 2016
p_25yo        <- make_plot(mut_25yo,        cor_25yo,        mut_col,  "Healthy 25yo cfDNA") # Teo et al 2018
p_70yo        <- make_plot(mut_70yo,        cor_70yo,        mut_col,  "Healthy 70yo cfDNA") # Teo et al 2018
p_bcell       <- make_plot(mut_bcell,       cor_bcell,       mut_col,  "B cell MNase-seq") # Gaffney et al 2012

################################################################################################
#### sig_6samples subset                                                                    ####
################################################################################################

sig_lung_ma     <- mut_lung_ma[sig_6samples == TRUE & !is.na(nuc_occupancy)]
sig_lung_snyder <- mut_lung_snyder[sig_6samples == TRUE & !is.na(nuc_occupancy)]
sig_25yo        <- mut_25yo[sig_6samples == TRUE & !is.na(nuc_occupancy)]
sig_70yo        <- mut_70yo[sig_6samples == TRUE & !is.na(nuc_occupancy)]
sig_bcell       <- mut_bcell[sig_6samples == TRUE & !is.na(nuc_occupancy)]

cor_sig_lung_ma     <- run_cor(sig_lung_ma)
cor_sig_lung_snyder <- run_cor(sig_lung_snyder)
cor_sig_25yo        <- run_cor(sig_25yo)
cor_sig_70yo        <- run_cor(sig_70yo)
cor_sig_bcell       <- run_cor(sig_bcell)

p_sig_lung_ma     <- make_plot(sig_lung_ma,     cor_sig_lung_ma,     high_col, "Lung cancer cfDNA") # Ma et al 2017
p_sig_lung_snyder <- make_plot(sig_lung_snyder, cor_sig_lung_snyder, high_col, "Lung cancer cfDNA") # Snyder et al 2016
p_sig_25yo        <- make_plot(sig_25yo,        cor_sig_25yo,        mut_col,  "Healthy 25yo cfDNA") # Teo et al 2018
p_sig_70yo        <- make_plot(sig_70yo,        cor_sig_70yo,        mut_col,  "Healthy 70yo cfDNA") # Teo et al 2018
p_sig_bcell       <- make_plot(sig_bcell,       cor_sig_bcell,       mut_col,  "B cell MNase-seq") # Gaffney et al 2012

overall_title_sig <- ggdraw() +
  draw_label("CCF z-score vs nucleosome occupancy (sig. mutations) — Abbosh et al. 2023",
             fontface = "bold", size = 14, x = 0.5, hjust = 0.5)

multi_panel_sig <- plot_grid(overall_title_sig,
                             plot_grid(p_sig_lung_ma, p_sig_lung_snyder, p_sig_25yo, p_sig_70yo, p_sig_bcell,
                                       ncol = 3, labels = c("A", "B", "C", "D", "E")),
                             ncol = 1, rel_heights = c(0.06, 1))

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_all_maps_sig6samples_abbosh.pdf"),
       multi_panel_sig, width = 10, height = 6)

################################################################################################
#### Tumour-specific nucleosome occupancy: lung cancer minus healthy 70yo                  ####
################################################################################################

cat("Computing tumour-specific occupancy difference...\n")

# Ma et al. lung - healthy 70yo
nuc_diff_ma <- merge(
  nuc_lung_ma[, .(chr, start, end, occ_lung = nuc_occupancy)],
  nuc_healthy_70yo[, .(chr, start, end, occ_healthy = nuc_occupancy)],
  by = c("chr", "start", "end")
)
nuc_diff_ma[, occ_diff := occ_lung - occ_healthy]
setkey(nuc_diff_ma, chr, start, end)

cat("Ma - Bins in lung map:", nrow(nuc_lung_ma), "\n")
cat("Ma - Bins in healthy 70yo map:", nrow(nuc_healthy_70yo), "\n")
cat("Ma - Shared bins:", nrow(nuc_diff_ma), "\n")
summary(nuc_diff_ma$occ_diff)

# Snyder et al. lung - healthy 70yo
nuc_diff_snyder <- merge(
  nuc_lung_snyder[, .(chr, start, end, occ_lung = nuc_occupancy)],
  nuc_healthy_70yo[, .(chr, start, end, occ_healthy = nuc_occupancy)],
  by = c("chr", "start", "end")
)
nuc_diff_snyder[, occ_diff := occ_lung - occ_healthy]
setkey(nuc_diff_snyder, chr, start, end)

cat("Snyder - Bins in lung map:", nrow(nuc_lung_snyder), "\n")
cat("Snyder - Bins in healthy 70yo map:", nrow(nuc_healthy_70yo), "\n")
cat("Snyder - Shared bins:", nrow(nuc_diff_snyder), "\n")
summary(nuc_diff_snyder$occ_diff)

# Map mutations to diff scores
mut_diff_ma <- foverlaps(
  mutation_summary_all,
  nuc_diff_ma[, .(chr, start, end, occ_diff)],
  by.x = c("chr", "pos", "pos_end"),
  by.y = c("chr", "start", "end"),
  type = "within",
  nomatch = NA
)

mut_diff_snyder <- foverlaps(
  mutation_summary_all,
  nuc_diff_snyder[, .(chr, start, end, occ_diff)],
  by.x = c("chr", "pos", "pos_end"),
  by.y = c("chr", "start", "end"),
  type = "within",
  nomatch = NA
)

cat("Ma - Mutations with diff score:", sum(!is.na(mut_diff_ma$occ_diff)), "\n")
cat("Ma - Mutations without:", sum(is.na(mut_diff_ma$occ_diff)), "\n")
cat("Snyder - Mutations with diff score:", sum(!is.na(mut_diff_snyder$occ_diff)), "\n")
cat("Snyder - Mutations without:", sum(is.na(mut_diff_snyder$occ_diff)), "\n")

# Correlation tests
cor_diff_ma     <- cor.test(mut_diff_ma$mean_z,     mut_diff_ma$occ_diff,     method = "spearman")
cor_diff_snyder <- cor.test(mut_diff_snyder$mean_z, mut_diff_snyder$occ_diff, method = "spearman")

cat("Ma tumour-specific diff: rho =", round(cor_diff_ma$estimate, 3),
    "p =", format(cor_diff_ma$p.value, scientific = TRUE, digits = 2), "\n")
cat("Snyder tumour-specific diff: rho =", round(cor_diff_snyder$estimate, 3),
    "p =", format(cor_diff_snyder$p.value, scientific = TRUE, digits = 2), "\n")

# Plot function for diff
make_diff_plot <- function(dt, cor_result, colour, title) {
  dt %>%
    filter(!is.na(occ_diff)) %>%
    ggplot(aes(x = occ_diff, y = mean_z)) +
    geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
    geom_smooth(method = "lm", colour = colour, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("rho = ", round(cor_result$estimate, 3),
                            "\np = ", format_p(cor_result$p.value)),
             size = 5) +
    labs(x = "Nucleosome occupancy difference",
         y = "Mean CCF z-score", # across timepoints
         title = title) +
    theme_cowplot(font_size = 14)
}

p_diff_ma <- make_diff_plot(mut_diff_ma, cor_diff_ma, low_col, "Tumour-specific")
p_diff_snyder <- make_diff_plot(mut_diff_snyder, cor_diff_snyder, low_col, "Tumour-specific")

multi_diff <- plot_grid(p_diff_ma, p_diff_snyder, ncol = 2, labels = c("A", "B"))

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occ_diff_lung_healthy70_abbosh.pdf"),
       multi_diff, width = 10, height = 6)


# ===================== Final figure with all 7 subplots =====================
row1 <- plot_grid(p_25yo, p_70yo, p_bcell, ncol = 3)
row2 <- plot_grid(p_lung_ma, p_lung_snyder, NULL, ncol = 3)
row3 <- plot_grid(p_diff_ma, p_diff_snyder, NULL, ncol = 3)

multi_panel <- plot_grid(row1, row2, row3, ncol = 1)

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_all_maps_abbosh.pdf"),
       multi_panel, width = 12, height = 16)
ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_all_maps_abbosh.svg"),
       multi_panel, width = 12, height = 16)

cat("Analysis complete.\n")