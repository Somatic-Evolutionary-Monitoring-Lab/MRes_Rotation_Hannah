#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Correlate mutation shedding rates to clinical features            ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-03-17

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

analysis_name <- 'clinical_features'
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

# Read in outlier mutations
outlier_mutations_path <- "outputs/ccf_zscores_multi_patient/20260319/outlier_mutations.csv"
outlier_mutations <- read.csv(outlier_mutations_path)

#########################################################################
#### Check for associations between lung cancer subtype and shedding ####
#########################################################################

# Check the overlap
unique(outlier_mutations$patient_name) %in% clinical$Shorter_ID

# Merge
outlier_clinical <- outlier_mutations %>%
  left_join(clinical, by = c("patient_name" = "Shorter_ID"))

# Show histologies in data by shedding type - but this is skewed by patients with large numbers of mutations
table(outlier_clinical$histology1_group_central.reviewed, outlier_clinical$shedding_status)

# Patient-level: does each patient have mostly high or mostly low shedding mutations?
patient_histology <- outlier_clinical %>%
  group_by(patient_name, histology1_group_central.reviewed, shedding_status) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = shedding_status, values_from = n, values_fill = 0)

print(patient_histology, n = 40)

# Patient level shedding
patient_level <- patient_histology %>%
  mutate(
    total = High + Low,
    pct_high = 100 * High / total,
    dominant = case_when(
      High > Low ~ "Mostly high",
      Low > High ~ "Mostly low",
      TRUE ~ "Equal"
    )
  )
table(patient_level$histology1_group_central.reviewed, patient_level$dominant)

# Patient-level plot
ggplot(patient_level, aes(x = histology1_group_central.reviewed, y = pct_high, 
                          fill = histology1_group_central.reviewed)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  scale_fill_manual(values = c("LUAD" = "#FF7F00", "LUSC" = "#6A3D9A", "Other" = "grey60")) +
  geom_hline(yintercept = 50, linetype = "dashed", colour = "grey40") +
  theme_cowplot() +
  labs(x = "Histology", y = "% of outlier mutations that are high-shedding",
       title = "High-shedding mutation proportion by histology") +
  theme(legend.position = "none")

ggsave(paste0(outputs.folder, "histology_pct_high_shedding.pdf"), width = 7, height = 6)



###################################################################
#### Plot distributions of high/low shed mutations per patient ####
###################################################################

# Order patients by pct_high
patient_level_plot <- patient_level %>%
  arrange(pct_high) %>%
  mutate(patient_name = factor(patient_name, levels = patient_name))

# Reshape for stacked bar
patient_long <- patient_level_plot %>%
  pivot_longer(cols = c(High, Low), names_to = "shedding", values_to = "n_mutations") %>%
  mutate(shedding = factor(shedding, levels = c("Low", "High")))

ggplot(patient_long, aes(x = patient_name, y = n_mutations, fill = shedding)) +
  geom_col() +
  scale_fill_manual(values = c("High" = "#D7191C", "Low" = "#2C7BB6")) +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 5)) +
  labs(x = "Patient", y = "Number of outlier mutations",
       fill = "Shedding status",
       title = "Most patients have exclusively high or low shedding mutations")

ggsave(paste0(outputs.folder, "patient_shedding_profile.pdf"), width = 12, height = 6)


###########################################################
#### Join all clinical features to patient_level once  ####
###########################################################

patient_level <- patient_level %>%
  left_join(clinical %>% select(Shorter_ID, 
                                pTNMStage_v8_lesion1_central.reviewed,
                                SizePath_lesion1_central.reviewed,
                                pTStage_v8_lesion1_central.reviewed,
                                pNStage_lesion1_central.reviewed,
                                smoking_status_group,
                                adjuvant_treatment_YN),
            by = c("patient_name" = "Shorter_ID"))

# Add derived stage grouping
patient_level <- patient_level %>%
  mutate(stage_group = case_when(
    pTNMStage_v8_lesion1_central.reviewed %in% c("1a", "1b", "2a") ~ "Early (I-IIA)",
    pTNMStage_v8_lesion1_central.reviewed %in% c("2b", "3a", "3b") ~ "Late (IIB-IIIB)",
    TRUE ~ NA_character_
  ))


###########################################################
#### Check for associations between shedding and stage ####
###########################################################

# Check distribution of stage vs dominant shedding type
table(patient_level$pTNMStage_v8_lesion1_central.reviewed, patient_level$dominant)

table(patient_level$stage_group, patient_level$dominant)
fisher.test(table(patient_level$stage_group, patient_level$dominant))

#################################################################
#### Check for associations between shedding and tumour size ####
#################################################################

wilcox.test(SizePath_lesion1_central.reviewed ~ dominant, data = patient_level)


####################################################################################
#### Check for associations between shedding and smoking and adjuvant treatment ####
####################################################################################

table(patient_level$smoking_status_group, patient_level$dominant)
table(patient_level$adjuvant_treatment_YN, patient_level$dominant)


table(patient_level$pNStage_lesion1_central.reviewed, patient_level$dominant)


###############################################
#### Summary plot of the clinical features ####
###############################################


# Build summary table of all clinical tests
clinical_summary <- data.frame(
  feature = c("Histology (LUSC vs LUAD)", 
              "Stage (Late vs Early)", 
              "N stage (N1-2 vs N0)",
              "Tumour size (>median vs ≤median)",
              "Smoking (Smoker vs Ex-smoker)",
              "Adjuvant (Yes vs No)"),
  stringsAsFactors = FALSE
)

# Histology
ht <- fisher.test(table(
  patient_level$histology1_group_central.reviewed %in% c("LUSC"),
  patient_level$dominant == "Mostly high"
))

# Stage
st <- fisher.test(table(
  patient_level$stage_group == "Late (IIB-IIIB)",
  patient_level$dominant == "Mostly high"
))

# N stage
nt <- fisher.test(table(
  patient_level$pNStage_lesion1_central.reviewed > 0,
  patient_level$dominant == "Mostly high"
))

# Tumour size (split at median)
median_size <- median(patient_level$SizePath_lesion1_central.reviewed, na.rm = TRUE)
sz <- fisher.test(table(
  patient_level$SizePath_lesion1_central.reviewed > median_size,
  patient_level$dominant == "Mostly high"
))

# Smoking
sm <- fisher.test(table(
  patient_level$smoking_status_group == "Smoker",
  patient_level$dominant == "Mostly high")[c("FALSE", "TRUE"), ])

# Adjuvant
aj <- fisher.test(table(
  patient_level$adjuvant_treatment_YN == "Adjuvant",
  patient_level$dominant == "Mostly high"
))

# Combine results
clinical_summary$OR <- c(ht$estimate, st$estimate, nt$estimate, sz$estimate, sm$estimate, aj$estimate)
clinical_summary$ci_lower <- c(ht$conf.int[1], st$conf.int[1], nt$conf.int[1], sz$conf.int[1], sm$conf.int[1], aj$conf.int[1])
clinical_summary$ci_upper <- c(ht$conf.int[2], st$conf.int[2], nt$conf.int[2], sz$conf.int[2], sm$conf.int[2], aj$conf.int[2])
clinical_summary$p_value <- c(ht$p.value, st$p.value, nt$p.value, sz$p.value, sm$p.value, aj$p.value)

# Order by OR
clinical_summary$feature <- factor(clinical_summary$feature, 
                                   levels = rev(clinical_summary$feature))

# Forest plot
ggplot(clinical_summary, aes(x = OR, y = feature)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2) +
  geom_point(size = 3) +
  geom_text(aes(label = paste0("p = ", round(p_value, 2))), 
            hjust = -0.2, vjust = -0.8, size = 3) +
  scale_x_log10() +
  theme_cowplot() +
  labs(x = "Odds ratio (mostly high vs mostly low shedding)", 
       y = "",
       title = "Clinical features and shedding direction"
       )

ggsave(paste0(outputs.folder, "clinical_summary_forest_plot.pdf"), width = 10, height = 5)



