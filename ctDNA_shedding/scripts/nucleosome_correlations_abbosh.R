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
library(readxl)
library(svglite)

# -----------------------------------------------------------------------------
# Make a folder for this analysis run
# -----------------------------------------------------------------------------

date <- gsub("-","",Sys.Date())

analysis_name <- 'nucleosome_mapping_continuous_abbosh'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )

# -----------------------------------------------------------------------------
# Load data from Abbosh et al. ctDNA
# -----------------------------------------------------------------------------

ctDNA_data <- read_fst("data/ctDNA_data_abbosh_pos_multiple.fst")

# -----------------------------------------------------------------------------
# Load nucleosome maps
# -----------------------------------------------------------------------------

nuc_lung_ma <- data.table::fread("data/SRA438908_lung_cancer_Ma2017_stable_100bp_hg19.bed")
colnames(nuc_lung_ma) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")

nuc_healthy_25yo <- data.table::fread("data/GSE114511_25yo_Teo_cfDNA_stable_100bp.bed")
colnames(nuc_healthy_25yo) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")

nuc_healthy_70yo <- data.table::fread("data/GSE114511_70yo_Teo_cfDNA_stable_100bp.bed")
colnames(nuc_healthy_70yo) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")

nuc_bcell <- data.table::fread("data/GSE36979_Gaffney2012_Bcells_MNase-seq_stable_100bp_hg19.bed")
colnames(nuc_bcell) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")

# -----------------------------------------------------------------------------
# Build mutation summary across all patients
# -----------------------------------------------------------------------------

mutation_summary_abbosh <- ctDNA_data_abbosh_pos_multiple %>%
  group_by(PublicationID, pure_mutation_id) %>%
  summarise(
    mean_z       = mean(ccf_z_score, na.rm = TRUE),
    se_z         = sd(ccf_z_score, na.rm = TRUE) / sqrt(n()),
    ci_lower     = mean_z - 1.96 * se_z,
    ci_upper     = mean_z + 1.96 * se_z,
    n_samples    = n(),
    sig_6samples = first(sig_6samples),
    .groups      = "drop"
  )

# Parse pure_mutation_id (format: chr:pos:alt) into components
mutation_summary_abbosh <- mutation_summary_abbosh %>%
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

# -----------------------------------------------------------------------------
# Helper function: annotate mutations with nucleosome occupancy
# -----------------------------------------------------------------------------

get_nuc_occupancy <- function(mutation_df, nuc_map) {
  mutation_df %>%
    mutate(nuc_occupancy = mapply(function(c, p) {
      idx <- which(nuc_map$chr == c & nuc_map$start <= p & nuc_map$end > p)
      if (length(idx) > 0) nuc_map$nuc_occupancy[idx[1]] else NA
    }, chr, pos))
}

# -----------------------------------------------------------------------------
# Annotate with each nucleosome map
# -----------------------------------------------------------------------------

mutation_summary_abbosh_nuc_lung_ma <- get_nuc_occupancy(mutation_summary_abbosh, nuc_lung_ma)
mutation_summary_abbosh_nuc_25yo    <- get_nuc_occupancy(mutation_summary_abbosh, nuc_healthy_25yo)
mutation_summary_abbosh_nuc_70yo    <- get_nuc_occupancy(mutation_summary_abbosh, nuc_healthy_70yo)
mutation_summary_abbosh_nuc_bcell   <- get_nuc_occupancy(mutation_summary_abbosh, nuc_bcell)

# Coverage summaries
for (nm in c("lung_ma", "25yo", "70yo", "bcell")) {
  df <- get(paste0("mutation_summary_abbosh_nuc_", nm))
  cat(nm, "- covered:", sum(!is.na(df$nuc_occupancy)),
      "| missing:", sum(is.na(df$nuc_occupancy)),
      "| pct:", round(mean(!is.na(df$nuc_occupancy)) * 100, 1), "%\n")
}

# -----------------------------------------------------------------------------
# Spearman correlations — all mutations
# -----------------------------------------------------------------------------

cor_lung_ma <- cor.test(mutation_summary_abbosh_nuc_lung_ma$mean_z,
                        mutation_summary_abbosh_nuc_lung_ma$nuc_occupancy,
                        method = "spearman")
cor_25yo    <- cor.test(mutation_summary_abbosh_nuc_25yo$mean_z,
                        mutation_summary_abbosh_nuc_25yo$nuc_occupancy,
                        method = "spearman")
cor_70yo    <- cor.test(mutation_summary_abbosh_nuc_70yo$mean_z,
                        mutation_summary_abbosh_nuc_70yo$nuc_occupancy,
                        method = "spearman")
cor_bcell   <- cor.test(mutation_summary_abbosh_nuc_bcell$mean_z,
                        mutation_summary_abbosh_nuc_bcell$nuc_occupancy,
                        method = "spearman")

# -----------------------------------------------------------------------------
# Multi-panel plot: all 4 nucleosome maps — all mutations
# -----------------------------------------------------------------------------

x_lim <- range(c(mutation_summary_abbosh_nuc_lung_ma$nuc_occupancy,
                 mutation_summary_abbosh_nuc_25yo$nuc_occupancy,
                 mutation_summary_abbosh_nuc_70yo$nuc_occupancy,
                 mutation_summary_abbosh_nuc_bcell$nuc_occupancy), na.rm = TRUE)
y_lim <- range(mutation_summary_abbosh$mean_z, na.rm = TRUE)

make_nuc_plot <- function(df, cor, colour, title) {
  df %>%
    filter(!is.na(nuc_occupancy)) %>%
    ggplot(aes(x = nuc_occupancy, y = mean_z)) +
    geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
    geom_smooth(method = "lm", colour = colour, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    coord_cartesian(xlim = x_lim, ylim = y_lim) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("rho = ", round(cor$estimate, 3),
                            "\np = ", round(cor$p.value, 3),
                            "\nn = ", sum(!is.na(df$nuc_occupancy))),
             size = 3.5) +
    labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
         title = title) +
    theme_cowplot(font_size = 11)
}

p_lung  <- make_nuc_plot(mutation_summary_abbosh_nuc_lung_ma, cor_lung_ma, "#D7191C",
                         "Lung cancer cfDNA (Ma et al. 2017)")
p_25yo  <- make_nuc_plot(mutation_summary_abbosh_nuc_25yo,    cor_25yo,    "#2C7BB6",
                         "Healthy 25yo cfDNA (Teo et al. 2018)")
p_70yo  <- make_nuc_plot(mutation_summary_abbosh_nuc_70yo,    cor_70yo,    "#FF7F00",
                         "Healthy 70yo cfDNA (Teo et al. 2018)")
p_bcell <- make_nuc_plot(mutation_summary_abbosh_nuc_bcell,   cor_bcell,   "#762A83",
                         "B cell MNase-seq (Gaffney et al. 2012)")

overall_title <- ggdraw() +
  draw_label("CCF z-score vs nucleosome occupancy — Abbosh et al. 2023",
             fontface = "bold", size = 14, x = 0.5, hjust = 0.5)

multi_panel <- plot_grid(overall_title,
                         plot_grid(p_lung, p_25yo, p_70yo, p_bcell,
                                   ncol = 2, labels = c("A", "B", "C", "D")),
                         ncol = 1, rel_heights = c(0.04, 1))

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_all_maps_abbosh.pdf"),
       multi_panel, width = 10, height = 10)

# -----------------------------------------------------------------------------
# Multi-panel plot: sig_6samples mutations only
# -----------------------------------------------------------------------------

# Shared axis limits across sig panels
x_lim_sig <- range(c(
  filter(mutation_summary_abbosh_nuc_lung_ma, sig_6samples)$nuc_occupancy,
  filter(mutation_summary_abbosh_nuc_25yo,    sig_6samples)$nuc_occupancy,
  filter(mutation_summary_abbosh_nuc_70yo,    sig_6samples)$nuc_occupancy,
  filter(mutation_summary_abbosh_nuc_bcell,   sig_6samples)$nuc_occupancy
), na.rm = TRUE)

y_lim_sig <- range(filter(mutation_summary_abbosh, sig_6samples)$mean_z, na.rm = TRUE)

make_nuc_plot_sig <- function(df, colour, title) {
  df_sig <- df %>% filter(sig_6samples == TRUE, !is.na(nuc_occupancy))
  cor    <- cor.test(df_sig$mean_z, df_sig$nuc_occupancy, method = "spearman")
  
  ggplot(df_sig, aes(x = nuc_occupancy, y = mean_z)) +
    geom_point(alpha = 0.5, size = 1.5, colour = "grey40") +
    geom_smooth(method = "lm", colour = colour, se = TRUE) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
    coord_cartesian(xlim = x_lim_sig, ylim = y_lim_sig) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
             label = paste0("rho = ", round(cor$estimate, 3),
                            "\np = ", round(cor$p.value, 3),
                            "\nn = ", nrow(df_sig)), size = 3.5) +
    labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
         title = title) +
    theme_cowplot(font_size = 11)
}

p_sig_lung  <- make_nuc_plot_sig(mutation_summary_abbosh_nuc_lung_ma, "#D7191C",
                                 "Lung cancer cfDNA (Ma et al. 2017)")
p_sig_25yo  <- make_nuc_plot_sig(mutation_summary_abbosh_nuc_25yo,    "#2C7BB6",
                                 "Healthy 25yo cfDNA (Teo et al. 2018)")
p_sig_70yo  <- make_nuc_plot_sig(mutation_summary_abbosh_nuc_70yo,    "#FF7F00",
                                 "Healthy 70yo cfDNA (Teo et al. 2018)")
p_sig_bcell <- make_nuc_plot_sig(mutation_summary_abbosh_nuc_bcell,   "#762A83",
                                 "B cell MNase-seq (Gaffney et al. 2012)")

overall_title_sig <- ggdraw() +
  draw_label("CCF z-score vs nucleosome occupancy (sig. mutations) — Abbosh et al. 2023",
             fontface = "bold", size = 14, x = 0.5, hjust = 0.5)

multi_panel_sig <- plot_grid(overall_title_sig,
                             plot_grid(p_sig_lung, p_sig_25yo, p_sig_70yo, p_sig_bcell,
                                       ncol = 4, labels = c("A", "B", "C", "D")),
                             ncol = 1, rel_heights = c(0.06, 1))

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_all_maps_sig6samples_abbosh.pdf"),
       multi_panel_sig, width = 24, height = 6)

# -----------------------------------------------------------------------------
# Nucleosome occupancy difference: lung cancer vs healthy 70yo
# -----------------------------------------------------------------------------

nuc_diff <- merge(
  nuc_lung_ma[, .(chr, start, end, occ_lung = nuc_occupancy)],
  nuc_healthy_70yo[, .(chr, start, end, occ_healthy = nuc_occupancy)],
  by = c("chr", "start", "end")
)

cat("Bins in lung map:", nrow(nuc_lung_ma), "\n")
cat("Bins in healthy 70yo map:", nrow(nuc_healthy_70yo), "\n")
cat("Shared bins:", nrow(nuc_diff), "\n")

nuc_diff[, occ_diff := occ_lung - occ_healthy]

mutation_summary_abbosh_nuc_diff <- mutation_summary_abbosh %>%
  mutate(occ_diff = mapply(function(c, p) {
    idx <- which(nuc_diff$chr == c & nuc_diff$start <= p & nuc_diff$end > p)
    if (length(idx) > 0) nuc_diff$occ_diff[idx[1]] else NA
  }, chr, pos))

cat("Mutations with diff score:", sum(!is.na(mutation_summary_abbosh_nuc_diff$occ_diff)), "\n")
cat("Mutations without:", sum(is.na(mutation_summary_abbosh_nuc_diff$occ_diff)), "\n")

cor_diff <- cor.test(mutation_summary_abbosh_nuc_diff$mean_z,
                     mutation_summary_abbosh_nuc_diff$occ_diff,
                     method = "spearman")

mutation_summary_abbosh_nuc_diff %>%
  filter(!is.na(occ_diff)) %>%
  ggplot(aes(x = occ_diff, y = mean_z)) +
  geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#4DAF4A", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_diff$estimate, 3),
                          "\np = ", round(cor_diff$p.value, 3)), size = 4) +
  labs(x = "Nucleosome occupancy difference (lung cancer - healthy 70yo)",
       y = "Mean CCF z-score",
       title = "CCF z-score vs tumour-specific nucleosome occupancy change — Abbosh et al. 2023") +
  theme_cowplot()

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occ_diff_lung_healthy70_abbosh.pdf"),
       width = 8, height = 6)










