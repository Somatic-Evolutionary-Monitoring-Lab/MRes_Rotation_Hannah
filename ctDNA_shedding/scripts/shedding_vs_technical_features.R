#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Evaluate whether technical features influence ctDNA shedding      ######
######                                                                    ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-04-01

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

analysis_name <- 'shedding_vs_technical_features'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


##############################################
#### Get inputs required for all analyses ####
##############################################

# Read in ctDNA data with CCF z-scores
ctDNA_data_path <- "outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst"
ctDNA_data_pos_multiple <- read_fst(ctDNA_data_path)


###########################################################################################
#### Correlation between mean_z and technical features                                 ####
###########################################################################################

# Build summary
mutation_correlation_data <- ctDNA_data_pos_multiple %>%
  group_by(patient_name, Pos) %>%
  summarise(
    mean_z            = mean(ccf_z_score, na.rm = TRUE),
    avg_major_cn      = mean(major_cn_mean, na.rm = TRUE),
    avg_multiplicity  = mean(mean_multiplicity, na.rm = TRUE),
    cv_depth          = sd(Depth, na.rm = TRUE) / mean(Depth, na.rm = TRUE),
    avg_depth         = mean(Depth, na.rm = TRUE),
    tnc               = first(tnc_type),
    .groups = "drop"
  )

# Define function to calculate correlation for plot labels
get_corr_label <- function(df, var_x) {
  res <- cor.test(df[[var_x]], df$mean_z, method = "spearman")
  paste0("rho = ", round(res$estimate, 3), " ; p = ", round(res$p.value, 4))
}

#------------------------------------------------------------------------------#
# PLOT A: Correlation with Multiplicity
#------------------------------------------------------------------------------#
p_mult <- ggplot(mutation_correlation_data, aes(x = avg_multiplicity, y = mean_z)) +
  geom_point(alpha = 0.2, colour = "grey30") +
  geom_smooth(method = "lm", colour = "#D7191C", fill = "pink") +
  theme_cowplot() +
  labs(
    title = "Shedding vs multiplicity",
    subtitle = get_corr_label(mutation_correlation_data, "avg_multiplicity"),
    x = "Mean multiplicity",
    y = "Mean CCF z-score"
  )

#------------------------------------------------------------------------------#
# PLOT B: Correlation with Major Copy Number
#------------------------------------------------------------------------------#
p_cn <- ggplot(mutation_correlation_data, aes(x = avg_major_cn, y = mean_z)) +
  geom_point(alpha = 0.2, colour = "grey30") +
  geom_smooth(method = "lm", colour = "#D7191C", fill = "pink") +
  theme_cowplot() +
  labs(
    title = "Shedding vs major copy number",
    subtitle = get_corr_label(mutation_correlation_data, "avg_major_cn"),
    x = "Avg major CN",
    y = "Mean CCF z-score"
  )

#------------------------------------------------------------------------------#
# PLOT C: Correlation with Depth Stability (CV)
#------------------------------------------------------------------------------#
p_cv <- ggplot(mutation_correlation_data, aes(x = cv_depth, y = mean_z)) +
  geom_point(alpha = 0.2, colour = "grey30") +
  geom_smooth(method = "lm", colour = "#D7191C", fill = "pink") +
  theme_cowplot() +
  labs(
    title = "Shedding vs depth stability",
    subtitle = get_corr_label(mutation_correlation_data, "cv_depth"),
    x = "CV of depth (SD/mean)",
    y = "Mean CCF z-score"
  )

#------------------------------------------------------------------------------#
# PLOT D: TNC Motif
#------------------------------------------------------------------------------#
p_tnc <- ggplot(mutation_correlation_data, aes(x = reorder(tnc, mean_z, FUN = median), y = mean_z)) +
  geom_boxplot(outlier.size = 0.5, fill = "grey90") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "black") +
  theme_cowplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6)) +
  labs(title = "Shedding distribution by TNC motif", x = "Trinucleotide context", y = "Mean CCF z-score")

#------------------------------------------------------------------------------#
# Combine and save
#------------------------------------------------------------------------------#
rule_out_2x2_plots <- plot_grid(p_mult, p_cn, p_cv, p_tnc, ncol = 2, align = "hv")

ggsave(paste0(outputs.folder, "shedding_vs_technical_features.pdf"),
       rule_out_2x2_plots, width = 12, height = 10)

#==============================================================================#
# Final correlation table
#==============================================================================#
cor_results <- data.frame(
  Variable = c("Major_CN", "Multiplicity", "CV_Depth", "Avg_Depth"),
  Spearman_Rho = c(
    cor(mutation_correlation_data$avg_major_cn, mutation_correlation_data$mean_z, method = "spearman"),
    cor(mutation_correlation_data$avg_multiplicity, mutation_correlation_data$mean_z, method = "spearman"),
    cor(mutation_correlation_data$cv_depth, mutation_correlation_data$mean_z, method = "spearman", use = "complete.obs"),
    cor(mutation_correlation_data$avg_depth, mutation_correlation_data$mean_z, method = "spearman")
  )
)

write.csv(cor_results, paste0(outputs.folder, "zscore_correlation_table.csv"), row.names = FALSE)
print(cor_results)
