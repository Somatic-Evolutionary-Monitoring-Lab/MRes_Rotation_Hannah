#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Load & explore CCF z-score distributions from Abbosh et al 2023   ######
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
library(readxl)
library(svglite)
source("scripts/plot_theme_mres_frankell.R")

# -----------------------------------------------------------------------------
# Make a folder for this analysis run
# -----------------------------------------------------------------------------

date <- gsub("-","",Sys.Date())

analysis_name <- 'ccf_analysis_abbosh'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )

# -----------------------------------------------------------------------------
# Load data from Excel sheet
# -----------------------------------------------------------------------------

ctDNA_data_abbosh <- read_excel(
  path = "data/41586_2023_5776_MOESM4_ESM.xlsx",
  sheet = "Table_17_patient_mutation_data",
  skip = 46,
  col_types = "text")

ctDNA_data_abbosh <- ctDNA_data_abbosh %>%
  mutate(
    mutation_ccf       = as.numeric(mutation_ccf),
    mean_tumour_ccf    = as.numeric(mean_tumour_ccf),
    days_post_surgery  = as.numeric(days_post_surgery),
    `mrd_caller_p-value` = as.numeric(`mrd_caller_p-value`),
    high_qual_eclipse  = as.logical(high_qual_eclipse),
    PyCloneClonal_SC   = as.character(PyCloneClonal_SC)
  )

# Create a row for mutation ID without sample ID
ctDNA_data_abbosh <- ctDNA_data_abbosh %>%
  mutate(pure_mutation_id = sub("^[^:]+:", "", mutation_id))


# -----------------------------------------------------------------------------
# Filter for clonal high-quality mutations
# -----------------------------------------------------------------------------

# Keep all clonal, ctDNA positive timepoints after surgery at high ppm
ctDNA_data_abbosh_pos <- ctDNA_data_abbosh %>%
  filter(`mrd_caller_p-value` < 0.01, 
         days_post_surgery >= 1,
         PyCloneClonal_SC == "C",
         high_qual_eclipse == TRUE, # this filters to > 1000 ppm
         !is.na(mutation_ccf)
         )

# Keep all mutations that appear in at least two samples in a patient
ctDNA_data_abbosh_pos_multiple <- ctDNA_data_abbosh_pos %>% 
  group_by(PublicationID, pure_mutation_id) %>% 
  filter(n() >= 2) %>% 
  ungroup()

# -----------------------------------------------------------------------------
# Calculate z-scores
# -----------------------------------------------------------------------------

ctDNA_data_abbosh_pos_multiple <- ctDNA_data_abbosh_pos_multiple %>%
  group_by(tracerx_id, days_post_surgery) %>% # using days_post_surgery as a proxy for sample, assuming no two samples were taken on the same day for the same patient
  mutate(ccf_z_score = (mutation_ccf - mean(mutation_ccf)) / sd(mutation_ccf)) %>% 
  ungroup()

# -----------------------------------------------------------------------------
# Distribution of samples per patient
# -----------------------------------------------------------------------------

p <- ctDNA_data_abbosh_pos_multiple %>%
  distinct(PublicationID, days_post_surgery) %>%
  count(PublicationID) %>%
  arrange(desc(n)) %>%
  mutate(PublicationID = factor(PublicationID, levels = PublicationID)) %>%
  ggplot(aes(x = PublicationID, y = n)) +
  geom_bar(stat = "identity", fill = "#7EB5D6") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8)) +
  labs(x = "Patient", y = "Number of samples",
       title = "Samples per patient (clonal, ctDNA+, post-surgery, high ppm) — Abbosh et al. 2023")

ggsave(paste0(outputs.folder, "samples_per_patient_abbosh.pdf"), p, width = 10, height = 6)
ggsave(paste0(outputs.folder, "samples_per_patient_abbosh.svg"), p, width = 10, height = 6)

# -----------------------------------------------------------------------------
# Z-score histograms for top 10 patients (rows x samples)
# -----------------------------------------------------------------------------

# Select top 10 patients by number of distinct timepoints
top10_patients <- ctDNA_data_abbosh_pos_multiple %>%
  group_by(PublicationID) %>%
  summarise(n_samples = n_distinct(days_post_surgery), .groups = "drop") %>%
  arrange(desc(n_samples)) %>%
  slice_head(n = 10) %>%
  pull(PublicationID)

# Filter to top 10 patients
ctDNA_top10 <- ctDNA_data_abbosh_pos_multiple %>%
  filter(PublicationID %in% top10_patients)

# Build histogram grid
patient_plots <- lapply(top10_patients, function(pat) {
  
  pat_data <- ctDNA_top10 %>% filter(PublicationID == pat)
  timepoints <- pat_data %>% arrange(days_post_surgery) %>% pull(days_post_surgery) %>% unique()
  timepoints <- head(timepoints, 6)  # cap at 6
  
  sample_plots <- lapply(timepoints, function(tp) {
    tp_data <- pat_data %>% filter(days_post_surgery == tp)
    
    samp_label <- paste0(tp, " days")
    
    ggplot(tp_data, aes(x = ccf_z_score)) +
      geom_histogram(binwidth = 0.3, fill = "#7EB5D6", colour = "white", alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
      xlim(-3, 3) +
      theme_cowplot(font_size = 14) +
      labs(x = "CCF z-score", y = "Count", title = samp_label)
  })
  
  # Pad to max 6 panels
  if (length(sample_plots) < 6) {
    empty <- replicate(6 - length(sample_plots),
                       ggplot() + theme_void(),
                       simplify = FALSE)
    sample_plots <- c(sample_plots, empty)
  }
  
  row   <- plot_grid(plotlist = sample_plots, nrow = 1)
  title <- ggdraw() + draw_label(paste0("Patient ", pat),
                                 fontface = "bold", x = 0.01, hjust = 0, size = 16)
  plot_grid(title, row, ncol = 1, rel_heights = c(0.1, 1))
})

# Add overall title
overall_title <- ggdraw() + 
  draw_label("CCF z-score histograms - Abbosh et al. 2023", 
             fontface = "bold", size = 16, x = 0.5, hjust = 0.5)

final_plot <- plot_grid(overall_title,
                        plot_grid(plotlist = patient_plots, ncol = 1),
                        ncol = 1, rel_heights = c(0.02, 1))

ggsave(paste0(outputs.folder, "top10_patients_zscore_histograms_abbosh.pdf"),
       final_plot, width = 18, height = 20)

# -----------------------------------------------------------------------------
# Spearman correlations between two timepoints for patients with exactly 2 samples 
# -----------------------------------------------------------------------------

# Find patients with exactly 2 timepoints
patients_2samples <- ctDNA_data_abbosh_pos_multiple %>%
  group_by(PublicationID) %>%
  summarise(n_samples = n_distinct(days_post_surgery), .groups = "drop") %>%
  filter(n_samples == 2) %>%
  pull(PublicationID)

# Take up to 16 patients
patients_2samples_16 <- head(patients_2samples, 16)

spearman_plots <- lapply(patients_2samples_16, function(pat) {
  
  pat_data <- ctDNA_data_abbosh_pos_multiple %>%
    filter(PublicationID == pat)
  
  two_timepoints <- pat_data %>%
    distinct(days_post_surgery) %>%
    arrange(days_post_surgery) %>%
    pull(days_post_surgery)
  
  tp1 <- two_timepoints[1]
  tp2 <- two_timepoints[2]
  
  wide <- pat_data %>%
    filter(days_post_surgery %in% c(tp1, tp2)) %>%
    select(pure_mutation_id, days_post_surgery, ccf_z_score) %>%
    pivot_wider(names_from = days_post_surgery, values_from = ccf_z_score) %>%
    drop_na()
  
  colnames(wide)[2:3] <- c("z_s1", "z_s2")
  
  spear   <- cor.test(wide$z_s1, wide$z_s2, method = "spearman")
  rho     <- round(spear$estimate, 3)
  p_val   <- signif(spear$p.value, 3)
  p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", p_val))
  sig     <- spear$p.value < 0.05
  
  # Point colour: blue if significant, grey if not
  point_colour <- ifelse(sig, "#0C447C", "grey70")
  title_colour <- ifelse(sig, "black", "grey60")
  
  ggplot(wide, aes(x = z_s1, y = z_s2)) +
    geom_point(size = 1.2, alpha = 0.65, colour = point_colour) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey70") +
    annotate("text",
             x = min(wide$z_s1, na.rm = TRUE),
             y = max(wide$z_s2, na.rm = TRUE),
             label = paste0("rho = ", rho, "\n", p_label),
             hjust = 0, vjust = 1, size = 3,
             colour = point_colour) +
    theme_cowplot(font_size = 10) +
    theme(plot.title = element_text(colour = title_colour)) +
    labs(
      x     = paste0("CCF z-score (", tp1, " days)"),
      y     = paste0("CCF z-score (", tp2, " days)"),
      title = pat
    )
})

overall_title <- ggdraw() +
  draw_label("Spearman correlations between timepoints - Abbosh et al. 2023",
             fontface = "bold", size = 16, x = 0.5, hjust = 0.5)

final_spearman_grid <- plot_grid(overall_title,
                                 plot_grid(plotlist = spearman_plots, ncol = 4, nrow = 4),
                                 ncol = 1, rel_heights = c(0.03, 1))

ggsave(paste0(outputs.folder, "spearman_2sample_patients_4x4_abbosh.pdf"),
       final_spearman_grid, width = 10, height = 10)
ggsave(paste0(outputs.folder, "spearman_2sample_patients_4x4_abbosh.svg"),
       final_spearman_grid, width = 10, height = 10)

# -----------------------------------------------------------------------------
# Violin plot: mutation-level CCF z-scores for patients with >= 6 timepoints
# -----------------------------------------------------------------------------

# Identify patients with at least 6 distinct timepoints
patients_6samples <- ctDNA_data_abbosh_pos_multiple %>%
  group_by(PublicationID) %>%
  summarise(n_samples = n_distinct(days_post_surgery), .groups = "drop") %>%
  filter(n_samples >= 6) %>%
  pull(PublicationID)

cat("Patients with >= 6 timepoints:", length(patients_6samples), "\n")
print(patients_6samples)

# Filter to those patients
ctDNA_abbosh_violin <- ctDNA_data_abbosh_pos_multiple %>%
  filter(PublicationID %in% patients_6samples) %>%
  mutate(PublicationID = as.factor(PublicationID))

# Compute Wilcoxon signed-rank test per patient-mutation
mutation_summary_abbosh <- ctDNA_abbosh_violin %>%
  group_by(PublicationID, pure_mutation_id) %>%
  summarise(
    mean_z      = mean(ccf_z_score),
    se_z        = sd(ccf_z_score) / sqrt(n()),
    ci_lower    = mean_z - 1.96 * se_z,
    ci_upper    = mean_z + 1.96 * se_z,
    mut_p_value = wilcox.test(ccf_z_score, mu = 0)$p.value,
    .groups     = "drop"
  ) %>%
  mutate(sig = mut_p_value < 0.05)

# Plot
pos <- position_jitter(width = 0.25, seed = 42)

ggplot(mutation_summary_abbosh, aes(x = PublicationID, y = mean_z, colour = sig)) +
  
  geom_violin(fill = "grey90", colour = "grey60", width = 0.8) +
  
  geom_point(data = filter(mutation_summary_abbosh, sig == FALSE),
             position = pos, size = 1.2, alpha = 0.6) +
  
  geom_point(data = filter(mutation_summary_abbosh, sig == TRUE),
             position = pos, size = 1.6, alpha = 1) +
  
  geom_errorbar(data = filter(mutation_summary_abbosh, sig == TRUE),
                aes(ymin = ci_lower, ymax = ci_upper),
                position = pos,
                width = 0.02,
                linewidth = 0.4) +
  
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  
  scale_colour_manual(
    values = c("FALSE" = "grey50", "TRUE" = "#E8829A"),
    labels = c("Not significant", "Significant"),
    name   = "Z ≠ 0"
  ) +
  labs(x = "Patient", y = "Mean CCF z-score",
       title = "CCF z-score distribution per patient — Abbosh et al. 2023",
       subtitle = "Patients with ≥ 6 timepoints; Wilcoxon test vs. 0") +
  theme_cowplot() +
  theme(
    legend.position = "top",
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 14),
    axis.text.x     = element_text(angle = 45, hjust = 1)
  )

ggsave(paste0(outputs.folder, "ccf_zscore_violin_abbosh.pdf"),
       width = 8, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_violin_abbosh.svg"),
       width = 8, height = 6)


# Extract significant mutations from >= 6 timepoint patients
sig_mutations_6sample_abbosh <- filter(mutation_summary_abbosh, sig == TRUE)

# Add sig_6samples column to the full dataset
ctDNA_data_abbosh_pos_multiple <- ctDNA_data_abbosh_pos_multiple %>%
  mutate(sig_6samples = paste(PublicationID, pure_mutation_id, sep = "_") %in%
           paste(sig_mutations_6sample_abbosh$PublicationID, 
                 sig_mutations_6sample_abbosh$pure_mutation_id, sep = "_"))

table(ctDNA_data_abbosh_pos_multiple$sig_6samples)


# ---------------------


pos <- position_jitter(width = 0.25, seed = 42)
ggplot(mutation_summary_abbosh, aes(x = PublicationID, y = mean_z, colour = sig)) +
  
  geom_violin(fill = "grey90", colour = "grey60", width = 0.8) +
  
  geom_point(data = filter(mutation_summary_abbosh, sig == FALSE),
             position = pos, size = 1.2, alpha = 0.6) +
  
  geom_point(data = filter(mutation_summary_abbosh, sig == TRUE),
             position = pos, size = 1.6, alpha = 1) +
  
  geom_errorbar(data = filter(mutation_summary_abbosh, sig == TRUE),
                aes(ymin = ci_lower, ymax = ci_upper),
                position = pos,
                width = 0.02,
                linewidth = 0.4) +
  
  geom_hline(yintercept = 0, linetype = "dashed", colour = horiz_line_col) +
  
  scale_colour_manual(
    values = c("FALSE" = wt_col, "TRUE" = mut_col),
    labels = c("Not significant", "Significant"),
    name   = "Z ≠ 0"
  ) +
  labs(x = "Patient", y = "Mean CCF z-score") +
  ggtitle("CCF z-score distribution per patient — Abbosh et al. 2023") +
  theme_cowplot() +
  theme(
    legend.position = "top",
    axis.text       = element_text(size = 11),
    axis.title      = element_text(size = 14),
  )

ggsave(paste0(outputs.folder, "ccf_zscore_violin_abbosh.pdf"), width = 10, height = 6)
ggsave(paste0(outputs.folder, "ccf_zscore_violin_abbosh.svg"), width = 10, height = 6)

# -----------------------------------------------------------------------------
# Save final dataframe
# -----------------------------------------------------------------------------

write_fst(ctDNA_data_abbosh_pos_multiple, "data/ctDNA_data_abbosh_pos_multiple.fst")












