#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Map high/low shedding mutations to nucleosome position map        ######
######  using different nucleosome maps: lung cancer cfDNA, healthy       ######
######  young and healthy old                                             ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-03-18

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")

####################################################
#### Source required functions & load libraries ####
####################################################

# suppress warning on R build version #
library(fst)
library(data.table) 
library(dplyr) 
library(ggplot2) 
library(cowplot) 
library(RColorBrewer) 
library(tidyr)

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'nucleosome_mapping_continuous'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


##############################################
#### Get Inputs required for all analyses ####
##############################################

# Read in ctDNA data with CCF z-scores
ctDNA_data_path <- "outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst"
ctDNA_data <- read_fst(ctDNA_data_path)

# Read in nucleosome map - LUNG CANCER MA ET AL 2017
nucleosome_map_path_lung_ma <- "data/SRA438908_lung_cancer_Ma2017_stable_100bp_hg19.bed"
nuc_lung_ma <- data.table::fread(nucleosome_map_path_lung_ma)
colnames(nuc_lung_ma) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")
head(nuc_lung_ma)

# Read in nucleosome map - HEALTHY 25YO TEO ET AL. 2018
nucleosome_map_path_healthy_25yo <- "data/GSE114511_25yo_Teo_cfDNA_stable_100bp.bed"
nuc_healthy_25yo <- data.table::fread(nucleosome_map_path_healthy_25yo)
colnames(nuc_healthy_25yo) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")
head(nuc_healthy_25yo)

# Read in nucleosome map - HEALTHY 70YO TEO ET AL. 2018
nucleosome_map_path_healthy_70yo <- "data/GSE114511_70yo_Teo_cfDNA_stable_100bp.bed"
nuc_healthy_70yo <- data.table::fread(nucleosome_map_path_healthy_70yo)
colnames(nuc_healthy_70yo) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")
head(nuc_healthy_70yo)

# Read in clinical features
clinical <- read.delim("data/tx842_clinical_outcome_20251211.tsv")
head(clinical)
colnames(clinical)


################################################################################################
#### What is the average age in this cohort?                                                ####
################################################################################################

my_patients <- unique(ctDNA_data$patient)
my_ages <- clinical %>% filter(Shorter_ID %in% my_patients) %>% pull(age)

summary(my_ages)
mean(my_ages, na.rm = TRUE)
median(my_ages, na.rm = TRUE)

###############################################################
#### Build mutation summary across all patients            ####
###############################################################

mutation_summary_all <- ctDNA_data %>%
  group_by(patient, Pos) %>%
  summarise(
    mean_z = mean(ccf_z_score),
    se_z = sd(ccf_z_score) / sqrt(n()),
    ci_lower = mean_z - 1.96 * se_z,
    ci_upper = mean_z + 1.96 * se_z,
    n_samples = n(),
    sig_6samples = first(sig_6samples),
    .groups = "drop"
  )

# Parse Pos into chr and position
mutation_summary_all <- mutation_summary_all %>%
  tidyr::separate(Pos, into = c("chr_num", "pos", "ref", "alt"), 
                  sep = ":", remove = FALSE) %>%
  mutate(
    chr = paste0("chr", chr_num),
    pos = as.numeric(pos)
  ) %>%
  select(-chr_num) %>% 
  select(patient, Pos, chr, pos, ref, alt, n_samples, mean_z, se_z, ci_lower, ci_upper, sig_6samples)

################################################################################################
#### Does CCF z-score correlate with nucleosome occupancy? LUNG CANCER CFDNA, MA ET AL 2017 ####
################################################################################################

# For each mutation, get the nucleosome occupancy of the overlapping bin
mutation_summary_all_nuc_lung_ma <- mutation_summary_all %>%
  mutate(nuc_occupancy = mapply(function(c, p) {
    idx <- which(nuc_lung_ma$chr == c & nuc_lung_ma$start <= p & nuc_lung_ma$end > p)
    if (length(idx) > 0) nuc_lung_ma$nuc_occupancy[idx[1]] else NA
  }, chr, pos))

# How many mutations have a nucleosome occupancy score?
cat("Mutations with nucleosome occupancy:", sum(!is.na(mutation_summary_all_nuc_lung_ma$nuc_occupancy)), "\n")
cat("Mutations without (unstable regions):", sum(is.na(mutation_summary_all_nuc_lung_ma$nuc_occupancy)), "\n")
cat("Proportion covered:", round(mean(!is.na(mutation_summary_all_nuc_lung_ma$nuc_occupancy)) * 100, 1), "%\n")

# Spearman correlation test between mean_z and nucleosome occupancy
cor_test <- cor.test(mutation_summary_all_nuc_lung_ma$mean_z,
                     mutation_summary_all_nuc_lung_ma$nuc_occupancy,
                     alternative = "two.sided",
                     method = "spearman"
                     )

# Plotting
mutation_summary_all_nuc_lung_ma %>%
  ggplot(aes(x = nuc_occupancy, y = mean_z)) +
  geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#D7191C", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  annotate("text",
           x = Inf, y = Inf,
           hjust = 1.1, vjust = 1.5,
           label = paste0("Spearman rho = ", round(cor_test$estimate, 3),
                          "\np = ", round(cor_test$p.value, 3)),
           size = 4) +
  labs(x = "Normalised nucleosome occupancy",
       y = "Mean CCF z-score",
       title = "CCF z-score vs nucleosome occupancy - Lung cancer cfDNA, Ma et al. 2017",
       subtitle = "n = 5,575 mutations") +
  theme_cowplot()

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_lung_ma_2017.pdf"),
       width = 10, height = 6)

################################################################################################
#### Does CCF z-score correlate with nucleosome occupancy? HEALTHY 25YO, TEO ET AL 2018     ####
################################################################################################

mutation_summary_all_nuc_25yo <- mutation_summary_all %>%
  mutate(nuc_occupancy = mapply(function(c, p) {
    idx <- which(nuc_healthy_25yo$chr == c & nuc_healthy_25yo$start <= p & nuc_healthy_25yo$end > p)
    if (length(idx) > 0) nuc_healthy_25yo$nuc_occupancy[idx[1]] else NA
  }, chr, pos))

cat("25yo - Mutations with nucleosome occupancy:", sum(!is.na(mutation_summary_all_nuc_25yo$nuc_occupancy)), "\n")
cat("25yo - Mutations without (unstable regions):", sum(is.na(mutation_summary_all_nuc_25yo$nuc_occupancy)), "\n")
cat("25yo - Proportion covered:", round(mean(!is.na(mutation_summary_all_nuc_25yo$nuc_occupancy)) * 100, 1), "%\n")

cor_test_25yo <- cor.test(mutation_summary_all_nuc_25yo$mean_z,
                          mutation_summary_all_nuc_25yo$nuc_occupancy,
                          alternative = "two.sided",
                          method = "spearman")


################################################################################################
#### Does CCF z-score correlate with nucleosome occupancy? HEALTHY 70YO, TEO ET AL 2018     ####
################################################################################################

mutation_summary_all_nuc_70yo <- mutation_summary_all %>%
  mutate(nuc_occupancy = mapply(function(c, p) {
    idx <- which(nuc_healthy_70yo$chr == c & nuc_healthy_70yo$start <= p & nuc_healthy_70yo$end > p)
    if (length(idx) > 0) nuc_healthy_70yo$nuc_occupancy[idx[1]] else NA
  }, chr, pos))

cat("70yo - Mutations with nucleosome occupancy:", sum(!is.na(mutation_summary_all_nuc_70yo$nuc_occupancy)), "\n")
cat("70yo - Mutations without (unstable regions):", sum(is.na(mutation_summary_all_nuc_70yo$nuc_occupancy)), "\n")
cat("70yo - Proportion covered:", round(mean(!is.na(mutation_summary_all_nuc_70yo$nuc_occupancy)) * 100, 1), "%\n")

cor_test_70yo <- cor.test(mutation_summary_all_nuc_70yo$mean_z,
                          mutation_summary_all_nuc_70yo$nuc_occupancy,
                          alternative = "two.sided",
                          method = "spearman")


################################################################################################
#### Multi-panel plot: all three nucleosome maps                                            ####
################################################################################################

# Get shared axis limits
all_occ <- c(
  mutation_summary_all_nuc_lung_ma$nuc_occupancy,
  mutation_summary_all_nuc_25yo$nuc_occupancy,
  mutation_summary_all_nuc_70yo$nuc_occupancy
)
all_z <- c(
  mutation_summary_all_nuc_lung_ma$mean_z,
  mutation_summary_all_nuc_25yo$mean_z,
  mutation_summary_all_nuc_70yo$mean_z
)

x_lim <- range(all_occ, na.rm = TRUE)
y_lim <- range(all_z, na.rm = TRUE)

# Build individual plots
p_lung <- mutation_summary_all_nuc_lung_ma %>%
  filter(!is.na(nuc_occupancy)) %>%
  ggplot(aes(x = nuc_occupancy, y = mean_z)) +
  geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#D7191C", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_test$estimate, 3),
                          "\np = ", round(cor_test$p.value, 3)), size = 3.5) +
  labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
       title = "Lung cancer cfDNA (Ma et al. 2017)") +
  theme_cowplot(font_size = 11)

p_25yo <- mutation_summary_all_nuc_25yo %>%
  filter(!is.na(nuc_occupancy)) %>%
  ggplot(aes(x = nuc_occupancy, y = mean_z)) +
  geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#2C7BB6", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_test_25yo$estimate, 3),
                          "\np = ", round(cor_test_25yo$p.value, 3)), size = 3.5) +
  labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
       title = "Healthy 25yo cfDNA (Teo et al. 2018)") +
  theme_cowplot(font_size = 11)

p_70yo <- mutation_summary_all_nuc_70yo %>%
  filter(!is.na(nuc_occupancy)) %>%
  ggplot(aes(x = nuc_occupancy, y = mean_z)) +
  geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#FF7F00", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  coord_cartesian(xlim = x_lim, ylim = y_lim) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_test_70yo$estimate, 3),
                          "\np = ", round(cor_test_70yo$p.value, 3)), size = 3.5) +
  labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
       title = "Healthy 70yo cfDNA (Teo et al. 2018)") +
  theme_cowplot(font_size = 11)

# Combine into one figure
multi_panel <- plot_grid(p_lung, p_25yo, p_70yo, ncol = 3, labels = c("A", "B", "C"))

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_three_maps.pdf"),
       multi_panel, width = 18, height = 6)


################################################################################################
#### Multi-panel plot: all three nucleosome maps for mutations significantly different      ####
#### CCF z-score distribution to 0 across 6 samples (from 3 patients)                       ####
################################################################################################

# Filter to sig_6samples mutations
sig_lung <- mutation_summary_all_nuc_lung_ma %>% filter(sig_6samples == TRUE, !is.na(nuc_occupancy))
sig_25yo <- mutation_summary_all_nuc_25yo %>% filter(sig_6samples == TRUE, !is.na(nuc_occupancy))
sig_70yo <- mutation_summary_all_nuc_70yo %>% filter(sig_6samples == TRUE, !is.na(nuc_occupancy))

# Correlation tests
cor_sig_lung <- cor.test(sig_lung$mean_z, sig_lung$nuc_occupancy, alternative = "two.sided", method = "spearman")
cor_sig_25yo <- cor.test(sig_25yo$mean_z, sig_25yo$nuc_occupancy, alternative = "two.sided", method = "spearman")
cor_sig_70yo <- cor.test(sig_70yo$mean_z, sig_70yo$nuc_occupancy, alternative = "two.sided", method = "spearman")

# Shared axis limits
sig_occ <- c(sig_lung$nuc_occupancy, sig_25yo$nuc_occupancy, sig_70yo$nuc_occupancy)
sig_z <- c(sig_lung$mean_z, sig_25yo$mean_z, sig_70yo$mean_z)
x_lim_sig <- range(sig_occ, na.rm = TRUE)
y_lim_sig <- range(sig_z, na.rm = TRUE)

# Plots
p_sig_lung <- ggplot(sig_lung, aes(x = nuc_occupancy, y = mean_z)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#D7191C", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  coord_cartesian(xlim = x_lim_sig, ylim = y_lim_sig) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_sig_lung$estimate, 3),
                          "\np = ", round(cor_sig_lung$p.value, 3),
                          "\nn = ", nrow(sig_lung)), size = 3.5) +
  labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
       title = "Lung cancer cfDNA (Ma et al. 2017)") +
  theme_cowplot(font_size = 11)

p_sig_25yo <- ggplot(sig_25yo, aes(x = nuc_occupancy, y = mean_z)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#2C7BB6", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  coord_cartesian(xlim = x_lim_sig, ylim = y_lim_sig) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_sig_25yo$estimate, 3),
                          "\np = ", round(cor_sig_25yo$p.value, 3),
                          "\nn = ", nrow(sig_25yo)), size = 3.5) +
  labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
       title = "Healthy 25yo cfDNA (Teo et al. 2018)") +
  theme_cowplot(font_size = 11)

p_sig_70yo <- ggplot(sig_70yo, aes(x = nuc_occupancy, y = mean_z)) +
  geom_point(alpha = 0.5, size = 1.5, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#FF7F00", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  coord_cartesian(xlim = x_lim_sig, ylim = y_lim_sig) +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_sig_70yo$estimate, 3),
                          "\np = ", round(cor_sig_70yo$p.value, 3),
                          "\nn = ", nrow(sig_70yo)), size = 3.5) +
  labs(x = "Normalised nucleosome occupancy", y = "Mean CCF z-score",
       title = "Healthy 70yo cfDNA (Teo et al. 2018)") +
  theme_cowplot(font_size = 11)

multi_panel_sig <- plot_grid(p_sig_lung, p_sig_25yo, p_sig_70yo, ncol = 3, labels = c("A", "B", "C"))

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occupancy_three_maps_sig6samples.pdf"),
       multi_panel_sig, width = 18, height = 6)


################################################################################################
#### Examine correlation between CCF z-score and difference between lung cancer cfDNA       ####
#### and healthy cfDNA (70yo)                                                               ####
################################################################################################

# Join lung and healthy 70yo maps on matching bins
nuc_diff <- merge(
  nuc_lung_ma[, .(chr, start, end, occ_lung = nuc_occupancy)],
  nuc_healthy_70yo[, .(chr, start, end, occ_healthy = nuc_occupancy)],
  by = c("chr", "start", "end")
)

# How many bins are shared?
cat("Bins in lung map:", nrow(nuc_lung_ma), "\n")
cat("Bins in healthy 70yo map:", nrow(nuc_healthy_70yo), "\n")
cat("Shared bins:", nrow(nuc_diff), "\n")

# Compute difference
nuc_diff[, occ_diff := occ_lung - occ_healthy]

summary(nuc_diff$occ_diff)

# For each mutation, get the occupancy difference
mutation_summary_all_nuc_diff <- mutation_summary_all %>%
  mutate(occ_diff = mapply(function(c, p) {
    idx <- which(nuc_diff$chr == c & nuc_diff$start <= p & nuc_diff$end > p)
    if (length(idx) > 0) nuc_diff$occ_diff[idx[1]] else NA
  }, chr, pos))

cat("Mutations with diff score:", sum(!is.na(mutation_summary_all_nuc_diff$occ_diff)), "\n")
cat("Mutations without:", sum(is.na(mutation_summary_all_nuc_diff$occ_diff)), "\n")

# Correlation test
cor_test_diff <- cor.test(mutation_summary_all_nuc_diff$mean_z,
                          mutation_summary_all_nuc_diff$occ_diff,
                          method = "spearman")
cor_test_diff

# Plot
mutation_summary_all_nuc_diff %>%
  filter(!is.na(occ_diff)) %>%
  ggplot(aes(x = occ_diff, y = mean_z)) +
  geom_point(alpha = 0.2, size = 0.8, colour = "grey40") +
  geom_smooth(method = "lm", colour = "#4DAF4A", se = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey60") +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60") +
  annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_test_diff$estimate, 3),
                          "\np = ", round(cor_test_diff$p.value, 3)), size = 4) +
  labs(x = "Nucleosome occupancy difference (lung cancer - healthy 70yo)",
       y = "Mean CCF z-score",
       title = "CCF z-score vs tumour-specific nucleosome occupancy change") +
  theme_cowplot()

ggsave(paste0(outputs.folder, "ccf_zscore_vs_nuc_occ_diff_lung_healthy70.pdf"),
       width = 8, height = 6)







