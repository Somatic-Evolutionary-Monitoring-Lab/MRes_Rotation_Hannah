#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Explore CCF z-score distributions across different mutations to   ######
######  identify high / low shedding mutations - starting with single     ######
######  patient                                                           ######
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
library(svglite)

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'ccf_zscores_single_patient'
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

# Load the LTX to CRUK conversion table
publication_key <- read.delim("data/tracerxPublicationKey_170221.txt", stringsAsFactors = FALSE)

publication_key_min <- publication_key %>%
  mutate(
    ShorterID = paste0(
      "LTX",
      sprintf("%03d", as.numeric(sub("^LTX0*", "", SampleID)))
    )
  )

# Convert patient ID to CRUK ID
ctDNA_data <- ctDNA_data %>%
  left_join(publication_key_min, by = c("patient" = "ShorterID")) %>%
  mutate(patient = PublicationID) %>%
  select(-PublicationID)

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


##################################################
#### How many samples does each patient have? ####
################################################## 

samples_per_patient <- ctDNA_data_pos %>%
  group_by(patient) %>%
  summarise(n_samples = n_distinct(sample)) %>%
  arrange(desc(n_samples))

ggplot(samples_per_patient, aes(x = reorder(patient, -n_samples), y = n_samples)) +
  geom_bar(stat = "identity", fill = "#7EB5D6") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
  labs(x = "Patient", y = "Number of samples", 
       title = "Samples per patient (clonal, ctDNA+, post-surgery, high ppm)")

ggsave(paste0(outputs.folder, "samples_per_patient.pdf"), width = 10, height = 6)
ggsave(paste0(outputs.folder, "samples_per_patient.svg"), width = 10, height = 6)

#############################################################################################
#### Analysis 1: mutation shedding analysis on one single sample from one single patient ####
#############################################################################################

# Keep data only for patient LTX030
ctDNA_LTX030 <- ctDNA_data_pos %>%
  filter(patient == 'LTX030')

# Calculate z-scores for patient LTX030
ccf_values_LTX030 <- ctDNA_LTX030$ccf
z_scores_LTX030 <- (ccf_values_LTX030 - mean(ccf_values_LTX030)) / sd(ccf_values_LTX030)

# Visualise z-score distribution LTX030 (one sample)

# Add z-scores back to dataframe
ctDNA_LTX030$ccf_z_score <- z_scores_LTX030

# Plot 1: Histogram of z-scores
p1 <- ggplot(ctDNA_LTX030, aes(x = ccf_z_score)) +
  geom_histogram(binwidth = 0.3, fill = "#2C7BB6", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", colour = "red") +
  theme_cowplot() +
  labs(x = "CCF z-score", y = "Count",
       title = "Distribution of CCF z-scores (LTX030)",
       subtitle = "Dashed = mean, dotted red = ±2 SD")

# Plot 2: Density plot
p2 <- ggplot(ctDNA_LTX030, aes(x = ccf_z_score)) +
  geom_density(fill = "#2C7BB6", alpha = 0.4) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", colour = "red") +
  theme_cowplot() +
  labs(x = "CCF z-score", y = "Density",
       title = "CCF z-score density (LTX030)")

# Combine and save
combined <- plot_grid(p1, p2, nrow = 2)
ggsave(paste0(outputs.folder, "LTX030_zscore_distribution.pdf"), 
       combined, width = 8, height = 8)


################################################################################################
#### Analysis 2: mutation shedding analysis on multiple samples from one single patient ####
################################################################################################

# Get data for patient LTX063
ctDNA_data_pos_LTX063 <- ctDNA_data_pos_multiple %>%
  filter(patient == 'LTX063')

# Split data by sample
ctDNA_data_pos_LTX063_sample1 <- ctDNA_data_pos_multiple %>%
  filter(sample == 'LTX063_5_26_2015')

ctDNA_data_pos_LTX063_sample2 <- ctDNA_data_pos_multiple %>%
  filter(sample == 'LTX063_6_30_2015')

# Calculate z-scores per sample
ccf_values_LTX063_sample1 <- ctDNA_data_pos_LTX063_sample1$ccf
z_scores_LTX063_sample1 <- (ccf_values_LTX063_sample1 - mean(ccf_values_LTX063_sample1)) / sd(ccf_values_LTX063_sample1)

ccf_values_LTX063_sample2 <- ctDNA_data_pos_LTX063_sample2$ccf
z_scores_LTX063_sample2 <- (ccf_values_LTX063_sample2 - mean(ccf_values_LTX063_sample2)) / sd(ccf_values_LTX063_sample2)

# Add z-scores back to dataframes
ctDNA_data_pos_LTX063_sample1$ccf_z_score <- z_scores_LTX063_sample1
ctDNA_data_pos_LTX063_sample2$ccf_z_score <- z_scores_LTX063_sample2

# Plot histograms
p_s1 <- ggplot(ctDNA_data_pos_LTX063_sample1, aes(x = ccf_z_score)) +
  geom_histogram(binwidth = 0.3, fill = "#2C7BB6", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", colour = "red") +
  theme_cowplot() +
  labs(x = "CCF z-score", y = "Count",
       title = "LTX063 - Sample 1 (26 May 2015)")

p_s2 <- ggplot(ctDNA_data_pos_LTX063_sample2, aes(x = ccf_z_score)) +
  geom_histogram(binwidth = 0.3, fill = "#D7191C", colour = "white", alpha = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "black") +
  geom_vline(xintercept = c(-2, 2), linetype = "dotted", colour = "red") +
  theme_cowplot() +
  labs(x = "CCF z-score", y = "Count",
       title = "LTX063 - Sample 2 (30 Jun 2015)")

combined <- plot_grid(p_s1, p_s2, nrow = 2)
ggsave(paste0(outputs.folder, "LTX063_zscore_distributions_by_sample.pdf"),
       combined, width = 8, height = 8)

#### Are the mutations that have lower CCF the same across both samples?

LTX063_sample1_low <- ctDNA_data_pos_LTX063_sample1 %>% 
  filter(ccf_z_score < -1.96)

LTX063_sample2_low <- ctDNA_data_pos_LTX063_sample2 %>% 
  filter(ccf_z_score < -1.96)

LTX063_low_shedding_shared <- intersect(LTX063_sample1_low$Pos, LTX063_sample2_low$Pos)

## Plot intersection

# Merge the two samples by pos
LTX063_merged <- inner_join(
  ctDNA_data_pos_LTX063_sample1 %>% select(Pos, ccf_z_score),
  ctDNA_data_pos_LTX063_sample2 %>% select(Pos, ccf_z_score),
  by = "Pos",
  suffix = c("_s1", "_s2")
) %>%
  mutate(category = case_when(
    ccf_z_score_s1 < -1.96 & ccf_z_score_s2 < -1.96 ~ "Low in both",
    ccf_z_score_s1 < -1.96 ~ "Low in sample 1 only",
    ccf_z_score_s2 < -1.96 ~ "Low in sample 2 only",
    TRUE ~ "Normal"
  ))

ggplot(LTX063_merged, aes(x = ccf_z_score_s1, y = ccf_z_score_s2, colour = category)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_vline(xintercept = -1.96, linetype = "dotted", colour = "grey40") +
  geom_hline(yintercept = -1.96, linetype = "dotted", colour = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  scale_colour_manual(values = c(
    "Low in both" = "#D7191C",
    "Low in sample 1 only" = "#FDAE61",
    "Low in sample 2 only" = "#2C7BB6",
    "Normal" = "grey80"
  )) +
  theme_cowplot() +
  labs(x = "CCF z-score (Sample 1)", y = "CCF z-score (Sample 2)",
       title = "CCF z-scores across two timepoints (LTX063)",
       colour = "Category")

ggsave(paste0(outputs.folder, "LTX063_zscore_scatter.pdf"), width = 7, height = 6)




####################################################################################
#### Correlation between CCF z-scores from sample 1 and sample 2 in one patient ####
####################################################################################

# Calculate Spearman correlation
spearman_result <- cor.test(LTX063_merged$ccf_z_score_s1, LTX063_merged$ccf_z_score_s2, 
                            method = "spearman")

rho <- round(spearman_result$estimate, 3)
p_val <- signif(spearman_result$p.value, 3)

# Format p-value label
p_label <- ifelse(p_val < 0.001, "p < 0.001", paste0("p = ", p_val))

ggplot(LTX063_merged, aes(x = ccf_z_score_s1, y = ccf_z_score_s2)) +
  geom_point(alpha = 0.7, size = 2, colour = "grey40") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  annotate("text", 
           x = min(LTX063_merged$ccf_z_score_s1, na.rm = TRUE), 
           y = max(LTX063_merged$ccf_z_score_s2, na.rm = TRUE),
           label = paste0("rho = ", rho, "\n", p_label),
           hjust = 0, vjust = 1, size = 4) +
  theme_cowplot() +
  labs(x = "CCF z-score (Sample 1)", y = "CCF z-score (Sample 2)",
       title = "CCF z-scores across two timepoints (LTX063)")

ggsave(paste0(outputs.folder, "LTX063_zscore_scatter_spearman.pdf"), width = 7, height = 6)


























