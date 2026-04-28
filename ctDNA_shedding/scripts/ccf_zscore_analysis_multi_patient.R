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

patient_plots <- lapply(top10_patients, function(pat) {
  
  pat_data <- ctDNA_top10 %>% filter(patient == pat)
  samples <- pat_data %>% arrange(days_post_surgery) %>% pull(sample) %>% unique()
  samples <- head(samples, 6)  # cap at 6
  
  sample_plots <- lapply(samples, function(samp) {
    samp_data <- pat_data %>% filter(sample == samp)
    
    samp_label <- paste0(unique(samp_data$days_post_surgery), " days post surgery")
    
    ggplot(samp_data, aes(x = ccf_z_score)) +
      geom_histogram(binwidth = 0.3, fill = "#2C7BB6", colour = "white", alpha = 0.8) +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
      xlim(-3, 3) +
      theme_cowplot(font_size = 12) +
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
                                 fontface = "bold", x = 0.01, hjust = 0, size = 14)
  plot_grid(title, row, ncol = 1, rel_heights = c(0.1, 1))
})

# Add overall title
overall_title <- ggdraw() +
  draw_label("CCF z-score histograms - Black et al. 2025",
             fontface = "bold", size = 16, x = 0.5, hjust = 0.5)

final_plot <- plot_grid(overall_title,
                        plot_grid(plotlist = patient_plots, ncol = 1),
                        ncol = 1, rel_heights = c(0.02, 1))

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
  group_by(patient_name, sample, days_post_surgery) %>%
  summarise(
    n        = n(),
    mean_ccf = mean(ccf),
    se       = sd(ccf) / sqrt(n),
    ci_lower = mean_ccf - 1.96 * se,
    ci_upper = mean_ccf + 1.96 * se,
    .groups  = "drop"
  )

ggplot(ci_summary, aes(x = factor(days_post_surgery), y = mean_ccf, colour = patient_name)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.3) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  facet_wrap(~ patient_name, scales = "free_x") +
  theme_cowplot(font_size = 8) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
        legend.position = "none") +
  labs(x = "Days post surgery", y = "Mean CCF",
       title = "Mean CCF ± 95% CI per patient per sample")

ggsave(paste0(outputs.folder, "mean_ccf_CI_per_patient.pdf"),
       width = 16, height = 16)



##########################################################################
#### Finding mutations that have consistently higher/lower shedding   ####
##########################################################################

# Group mutations per patient
mutation_mean_z <- ctDNA_data_pos_multiple %>% 
  group_by(patient_name, Pos) %>% 
  summarise(
    n_samples = n(),
    mean_z = mean(ccf_z_score),
    .groups = "drop"
  )

# Count number of samples each patient has
patient_sample_counts <- ctDNA_data_pos_multiple %>%
  group_by(patient_name) %>%
  summarise(total_samples = n_distinct(sample), .groups = "drop")

# Compute patient order (those with most outlier mutations first)
patient_order <- mutation_mean_z %>% 
  mutate(is_outlier = case_when(
    mean_z > 1.96 ~ TRUE,
    mean_z < -1.96 ~ TRUE,
    TRUE ~ FALSE
  )) %>% 
  group_by(patient_name) %>% 
  summarise(n_outliers = sum(is_outlier)) %>% 
  arrange(desc(n_outliers)) %>% 
  pull(patient_name)

# Make multipage PDF plots
pdf(paste0(outputs.folder, "mutation_mean_zscore_per_patient.pdf"), 
    width = 12, height = 6)

for (pat in patient_order) {
  
  # Categorise mutations by shedding status
  pat_data <- mutation_mean_z %>%
    filter(patient_name == pat) %>% 
    mutate(shedding_status = case_when(
      mean_z > 1.96  ~ "High",
      mean_z < -1.96 ~ "Low",
      TRUE           ~ "Normal"
    ))
  
  # Plot
  p <- ggplot(pat_data, aes(x = reorder(Pos, -mean_z), y = mean_z, colour = shedding_status)) +
    geom_point(size = 2) + 
    theme_cowplot() +
    geom_hline(yintercept = -1.96, linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept =  1.96, linetype = "dotted", colour = "grey40") +
    geom_hline(yintercept =  0,    linetype = "dashed",  colour = "black") +
    scale_colour_manual(values = c(
      "High"   = "#D7191C",
      "Low"    = "#2C7BB6",
      "Normal" = "grey60"
    )) +
    labs(
      x        = "Mutation",
      y        = "Mean z-score across samples",
      colour   = "Shedding status",
      title    = paste0("CCF z-scores for patient ", pat),
      subtitle = paste0("n mutations = ", nrow(pat_data),
                        " ; n samples = ",
                        patient_sample_counts$total_samples[patient_sample_counts$patient_name == pat])) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
  
  print(p)
}

dev.off()

# Save outlier mutations as csv
outlier_mutations <- mutation_mean_z %>%
  mutate(shedding_status = case_when(
    mean_z > 1.96  ~ "High",
    mean_z < -1.96 ~ "Low",
    TRUE           ~ "Normal"
  )) %>%
  filter(shedding_status != "Normal") %>%
  arrange(patient_name, desc(abs(mean_z)))

write.csv(outlier_mutations, 
          file = paste0(outputs.folder, "outlier_mutations.csv"), 
          row.names = FALSE)

# Save as BED file for later nucleosome mapping
outlier_bed <- outlier_mutations %>%
  tidyr::separate(Pos, into = c("chr_num", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>%
  mutate(
    chr   = paste0("chr", chr_num),
    start = as.numeric(pos) - 1,  # BED is 0-based
    end   = as.numeric(pos)
  ) %>%
  select(chr, start, end, Pos, patient_name, shedding_status, mean_z)

outlier_bed_path <- paste0(outputs.folder, "outlier_mutations.bed")
write.table(outlier_bed, file = outlier_bed_path, 
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



################################################################
#### Plot the violin style plot for three 6-sample patients ####
################################################################

# Convert patient name to factor for plotting
ctDNA_data_pos_multiple_violin <- ctDNA_data_pos_multiple
ctDNA_data_pos_multiple_violin$patient <- as.factor(ctDNA_data_pos_multiple_violin$patient)

# Plot only for patients with 6 samples (for now)
ctDNA_data_pos_multiple_violin_6_samples <- ctDNA_data_pos_multiple_violin %>% 
  filter(patient %in% c("LTX208", "LTX287", "LTX854"))

# Compute the Wilcoxon signed-rank test per mutation to test whether its 
# CCF z-scores are systematically different from 0 across all 6 samples
mutation_summary <- ctDNA_data_pos_multiple_violin_6_samples %>% 
  group_by(patient, Pos) %>% 
  summarise(
    mean_z = mean(ccf_z_score),
    se_z = sd(ccf_z_score) / sqrt(n()),
    ci_lower = mean_z - 1.96 * se_z,
    ci_upper = mean_z + 1.96 * se_z,
    mut_p_value = wilcox.test(ccf_z_score, mu = 0)$p.value,
    .groups = "drop"
  ) %>% 
  mutate(sig = mut_p_value < 0.05)

# Plotting
pos <- position_jitter(width = 0.25, seed = 42)

ggplot(mutation_summary, aes(x = patient, y = mean_z, colour = sig)) +
  
  geom_violin(fill = "grey90", colour = "grey60", width = 0.8) +
  
  # Non-significant points
  geom_point(data = filter(mutation_summary, sig == FALSE),
             position = pos, size = 1.2, alpha = 0.6) +
  
  # Significant points + errorbars
  geom_point(data = filter(mutation_summary, sig == TRUE),
             position = pos, size = 1.6, alpha = 1) +
  
  geom_errorbar(data = filter(mutation_summary, sig == TRUE),
                aes(ymin = ci_lower, ymax = ci_upper),
                position = pos,
                width = 0.02,
                linewidth = 0.4) +
  
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  scale_colour_manual(
    values = c("FALSE" = "grey50", "TRUE" = "#D7191C"),
    labels = c("Not significant", "Significant"),
    name = "Z ≠ 0"
  ) +
  labs(x = "Patient", y = "Mean CCF z-score") +
  theme_cowplot() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 14)
  ) +
  ggtitle("CCF z-score distribution per patient (6-sample patients)")

ggsave(paste0(outputs.folder, "ccf_zscore_violin.pdf"),
       width = 8, height = 6)

# Examine the significant mutations from those 6-sample patients ####
sig_mutations_6sample <- filter(mutation_summary, sig == TRUE)

########################################################
#### Save the final ctDNA dataset with CCF z-scores ####
########################################################

# Add sig_6samples column to the full dataset
ctDNA_data_pos_multiple <- ctDNA_data_pos_multiple %>%
  mutate(sig_6samples = paste(patient, Pos, sep = "_") %in% 
           paste(sig_mutations_6sample$patient, sig_mutations_6sample$Pos, sep = "_"))

table(ctDNA_data_pos_multiple$sig_6samples)

write_fst(ctDNA_data_pos_multiple, paste0(outputs.folder, "ctDNA_data_pos_multiple.fst"))


###########################################################################################
#### Spearman correlations between two samples for subset of patients with two samples ####
###########################################################################################

# Find patients with exactly 2 samples
patients_2samples <- ctDNA_data_pos_multiple %>%
  group_by(patient) %>%
  summarise(n_samples = n_distinct(sample), .groups = "drop") %>%
  filter(n_samples == 2) %>%
  pull(patient)

# Take up to 16 patients
patients_2samples_16 <- head(patients_2samples, 16)

spearman_plots <- lapply(patients_2samples_16, function(pat) {
  
  pat_data <- ctDNA_data_pos_multiple %>%
    filter(patient == pat)
  
  two_samples <- pat_data %>%
    distinct(sample, days_post_surgery) %>%
    arrange(days_post_surgery) %>%
    pull(sample)
  
  s1 <- two_samples[1]
  s2 <- two_samples[2]
  
  wide <- pat_data %>%
    filter(sample %in% c(s1, s2)) %>%
    select(Pos, sample, ccf_z_score) %>%
    pivot_wider(names_from = sample, values_from = ccf_z_score) %>%
    drop_na()
  
  colnames(wide)[2:3] <- c("z_s1", "z_s2")
  
  spear   <- cor.test(wide$z_s1, wide$z_s2, method = "spearman")
  rho     <- round(spear$estimate, 3)
  p_val   <- signif(spear$p.value, 3)
  p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", p_val))
  sig     <- spear$p.value < 0.05
  
  point_colour <- ifelse(sig, "#0C447C", "grey70")
  title_colour <- ifelse(sig, "black", "grey60")
  
  days_s1 <- pat_data %>% filter(sample == s1) %>% pull(days_post_surgery) %>% unique()
  days_s2 <- pat_data %>% filter(sample == s2) %>% pull(days_post_surgery) %>% unique()
  
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
      x     = paste0("CCF z-score (", days_s1, " days)"),
      y     = paste0("CCF z-score (", days_s2, " days)"),
      title = pat
    )
})

overall_title <- ggdraw() +
  draw_label("Spearman correlations between timepoints - Black et al. 2025",
             fontface = "bold", size = 16, x = 0.5, hjust = 0.5)

final_spearman_grid <- plot_grid(overall_title,
                                 plot_grid(plotlist = spearman_plots, ncol = 4, nrow = 4),
                                 ncol = 1, rel_heights = c(0.03, 1))

ggsave(paste0(outputs.folder, "spearman_2sample_patients_4x4.pdf"),
       final_spearman_grid, width = 10, height = 10)









