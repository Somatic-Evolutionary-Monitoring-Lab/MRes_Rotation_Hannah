#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Explore CCF z-score distributions across different mutations to   ######
######  identify high / low shedding mutations - multiple patients        ######
######                                                                    ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-03-16

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

analysis_name <- 'ccf_zscores_multi_patient'
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

ctDNA_data_path <- "data/20260307_tracked_mutations_primary_and_met_data_eclipse_annotated.fst"

# Read in table
ctDNA_data <- read_fst(ctDNA_data_path)


#######################################################
#### Filter data to keep clonal high ppm mutations ####
#######################################################

# Keep all clonal, ctDNA positive timepoints after surgery at high ppm
ctDNA_data_pos <- ctDNA_data %>%
  filter(mrd_p_value < 0.01, 
         days_post_surgery >= 1, 
         !is.na(mutation_cluster), 
         !is.na(is_subclonal_sample_p),
         ppm > 1000,
         is_subclonal_sample == FALSE,
         is_trunk == TRUE,
         !is.na(ccf))

# Keep all mutations that appear in at least two samples in a patient
ctDNA_data_pos_multiple <- ctDNA_data_pos %>% 
  group_by(patient, Pos) %>% 
  filter(n() >= 2) %>% 
  ungroup()

n_patients_multiple <- length(unique(ctDNA_data_pos_multiple$patient))


#######################################################################################
#### Analysis 1: Plot for top 10 patients the z-score distributions across samples ####
#######################################################################################


################################################################
#### Calculate z-scores per sample across multiple patients ####
################################################################

# Select top 10 patients by number of samples - later will plot all patients
top10_patients <- ctDNA_data_pos_multiple %>% 
  group_by(patient_name) %>%
  summarise(n_samples = n_distinct(sample)) %>%
  arrange(desc(n_samples)) %>% 
  slice_head(n = 10) %>% 
  pull(patient_name)

# Filter top 10 patients
ctDNA_top10 <- ctDNA_data_pos_multiple %>% 
  filter(patient_name %in% top10_patients)

# Compute z-scores per patient-sample group
ctDNA_top10 <- ctDNA_top10 %>%
  group_by(patient_name, sample) %>% 
  mutate(ccf_z_score = (ccf - mean(ccf)) / sd(ccf)) %>% 
  ungroup()

################################################################
#### Calculate z-scores per sample across multiple patients ####
################################################################

# Select top 10 patients by number of samples
top10_patients <- ctDNA_data_pos_multiple %>%
  group_by(patient) %>%
  summarise(n_samples = n_distinct(sample)) %>%
  arrange(desc(n_samples)) %>%
  slice_head(n = 10) %>%
  pull(patient)

# Filter to top 10 patients
ctDNA_top10 <- ctDNA_data_pos_multiple %>%
  filter(patient %in% top10_patients)

# Compute z-scores per patient-sample group
ctDNA_top10 <- ctDNA_top10 %>%
  group_by(patient, sample) %>%
  mutate(ccf_z_score = (ccf - mean(ccf)) / sd(ccf)) %>%
  ungroup()

##################################################################
#### Plot z-score histograms: rows = patients, cols = samples ####
##################################################################

# Create a list of plots, one per patient
patient_plots <- lapply(top10_patients, function(pat) {
  
  pat_data <- ctDNA_top10 %>% filter(patient_name == pat)
  samples <- unique(pat_data$sample)
  
  # One histogram per sample
  sample_plots <- lapply(samples, function(samp) {
    samp_data <- pat_data %>% filter(sample == samp)
    
    # Shorten sample label for plot title
    samp_label <- gsub(paste0(pat, "_"), "", samp)
    
    ggplot(samp_data, aes(x = ccf_z_score)) +
      geom_histogram(binwidth = 0.3, fill = "#2C7BB6", colour = "white", alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
      geom_vline(xintercept = c(-1.96, 1.96), linetype = "dotted", colour = "red") +
      xlim(-4, 4) +  # fixed x-axis across all panels
      theme_cowplot(font_size = 8) +
      labs(x = "CCF z-score", y = "Count", title = samp_label)
  })
  
  # Pad with empty plots if fewer than 6 samples
  max_samples <- 6
  if (length(sample_plots) < max_samples) {
    empty <- replicate(max_samples - length(sample_plots), 
                       ggplot() + theme_void(), 
                       simplify = FALSE)
    sample_plots <- c(sample_plots, empty)
  }
  
  # Combine samples into one row, with patient label on left
  row <- plot_grid(plotlist = sample_plots, nrow = 1)
  title <- ggdraw() + draw_label(pat, fontface = "bold", x = 0, hjust = 0, size = 10)
  plot_grid(title, row, ncol = 1, rel_heights = c(0.1, 1))
})

# Stack all patient rows
final_plot <- plot_grid(plotlist = patient_plots, ncol = 1)

ggsave(paste0(outputs.folder, "top10_patients_zscore_histograms.pdf"),
       final_plot, width = 18, height = 20)


######################################################################################
#### Analysis 2: Identify consistently low/high shedding mutations across samples ####
######################################################################################

# Compute z-scores for all patients
ctDNA_data_pos_multiple <- ctDNA_data_pos_multiple %>%
  group_by(patient_name, sample) %>%
  mutate(ccf_z_score = (ccf - mean(ccf)) / sd(ccf)) %>%
  ungroup()

# For each patient-mutation, summarise across all samples
mutation_consistency <- ctDNA_data_pos_multiple %>%
  group_by(patient_name, Pos) %>%
  summarise(
    n_samples        = n(),
    n_low            = sum(ccf_z_score < -1.96),
    n_high           = sum(ccf_z_score > 1.96),
    mean_z           = mean(ccf_z_score),
    consistently_low  = (n_low / n_samples) >= 0.5,
    consistently_high = (n_high / n_samples) >= 0.5,
    .groups = "drop"
  )

# Summarise per patient
patient_summary <- mutation_consistency %>%
  group_by(patient_name) %>%
  summarise(
    n_mutations           = n(),
    n_consistently_low    = sum(consistently_low),
    n_consistently_high   = sum(consistently_high),
    pct_consistently_low  = round(100 * n_consistently_low / n_mutations, 1),
    pct_consistently_high = round(100 * n_consistently_high / n_mutations, 1),
    .groups = "drop"
  ) %>%
  arrange(desc(n_consistently_low))

print(patient_summary)


###############################################################
#### Visualise consistent low/high shedding per patient    ####
###############################################################

# Reshape to long format for plotting
patient_summary_long <- patient_summary %>%
  select(patient_name, n_consistently_low, n_consistently_high) %>%
  tidyr::pivot_longer(
    cols = c(n_consistently_low, n_consistently_high),
    names_to = "category",
    values_to = "n_mutations"
  ) %>%
  mutate(category = recode(category,
                           "n_consistently_low"  = "Consistently low (< -1.96)",
                           "n_consistently_high" = "Consistently high (> 1.96)"
  ))

# Order patients by total consistently aberrant mutations
patient_order <- patient_summary %>%
  mutate(total = n_consistently_low + n_consistently_high) %>%
  arrange(desc(total)) %>%
  pull(patient_name)

patient_summary_long$patient_name <- factor(patient_summary_long$patient_name, 
                                            levels = patient_order)

# Plot
ggplot(patient_summary_long, aes(x = patient_name, y = n_mutations, colour = category)) +
  geom_segment(aes(xend = patient_name, y = 0, yend = n_mutations), linewidth = 0.8) +
  geom_point(size = 3) +
  facet_wrap(~ category, ncol = 1) +
  scale_colour_manual(values = c(
    "Consistently low (< -1.96)"  = "#2C7BB6",
    "Consistently high (> 1.96)" = "#D7191C"
  )) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none") +
  labs(x = "Patient", y = "Number of mutations",
       title = "Consistently low/high shedding mutations per patient",
       subtitle = "Mutations with |z| > 1.96 in at least half the samples")

ggsave(paste0(outputs.folder, "consistent_shedding_per_patient.pdf"),
       width = 12, height = 8)





###############################################################
#### Mean CCF z-score with 95% CI per patient per sample   ####
###############################################################

ci_summary <- ctDNA_data_pos_multiple %>%
  group_by(patient_name, sample) %>%
  summarise(
    n         = n(),
    mean_z    = mean(ccf_z_score),
    se        = sd(ccf_z_score) / sqrt(n),
    ci_lower  = mean_z - 1.96 * se,
    ci_upper  = mean_z + 1.96 * se,
    .groups = "drop"
  )

ggplot(ci_summary, aes(x = sample, y = mean_z, colour = patient_name)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ patient_name, scales = "free_x") +
  theme_cowplot(font_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "none") +
  labs(x = "Sample", y = "Mean CCF z-score",
       title = "Mean CCF z-score ± 95% CI per patient per sample")

ggsave(paste0(outputs.folder, "mean_ccf_zscore_CI_per_patient.pdf"),
       width = 16, height = 16)
















