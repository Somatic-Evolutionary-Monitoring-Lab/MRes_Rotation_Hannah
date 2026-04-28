#==============================================================================#
#                                                                              #
#          Do cfChromatin RPKM scores correlate with shedding (CCF z-score)?   #
#                                                                              #
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-04-28

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")

# -----------------------------------------------------------------------------
# Load libraries
# -----------------------------------------------------------------------------

library(data.table)
library(fst)
library(readr)
library(tidyverse)
library(ggplot2)
library(cowplot)

# -----------------------------------------------------------------------------
# Make output folder
# -----------------------------------------------------------------------------

date <- gsub("-", "", Sys.Date())

analysis_name <- 'cfChromatin_vs_shedding'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep = '/')
if (!file.exists(out_dir_general)) dir.create(out_dir_general)

out_dir_logs <- paste(out_dir_general, 'logs', sep = '/')
if (!file.exists(out_dir_logs)) dir.create(out_dir_logs)

outputs.folder <- paste0(out_dir_general, "/", date, "/")
if (!file.exists(outputs.folder)) dir.create(outputs.folder)

# -----------------------------------------------------------------------------
# Load input data
# -----------------------------------------------------------------------------

# Load cfChromatin data (merged in cfChromatin_avg_RPKM_scores.R)
cfChrom_data <- fread("data/A549_cfChromatin_merged_hg19.bedGraph", 
              col.names = c("chr", "start", "end", "RPKM_score"))

# Load ctDNA data with CCF z-scores
ctDNA_data_path <- "outputs/ccf_zscores_multi_patient/20260319/ctDNA_data_pos_multiple.fst"
ctDNA_data <- read_fst(ctDNA_data_path)

# Load clinical data
clinical_data <- read_tsv("data/tx842_clinical_outcome_20251211.tsv")

# -----------------------------------------------------------------------------
# Clean up data
# -----------------------------------------------------------------------------

# Keep only standard chromosomes (not unplaced / unlocalised contigs or mitochondrial chromosome)
standard_chrs <- paste0("chr", c(1:22, "X", "Y"))
cfChrom_data <- cfChrom_data[chr %in% standard_chrs]

# Verify
table(cfChrom_data$chr)

# Add column for cancer type
ctDNA_data <- ctDNA_data %>% left_join(clinical_data, by=c("patient" = "Shorter_ID"))
ctDNA_data %>%
  distinct(patient, histology1_group_central.reviewed) %>%
  count(histology1_group_central.reviewed)

# Convert to data.table
setDT(ctDNA_data)
setDT(cfChrom_data)

# Add chr prefix
ctDNA_data[, chromosome := paste0("chr", chromosome)]

# -----------------------------------------------------------------------------
# Assign cfChromatin RPKM score to each mutation based on 10kb bin
# -----------------------------------------------------------------------------

# foverlaps needs a start and end for the query - for SNVs they are the same
ctDNA_data[, pos_end := position]

# Set keys
setkey(cfChrom_data, chr, start, end)
setkey(ctDNA_data, chromosome, position, pos_end)

# Overlap join
ctDNA_scored <- foverlaps(ctDNA_data, cfChrom_data,
                          by.x = c("chromosome", "position", "pos_end"),
                          by.y = c("chr", "start", "end"),
                          type = "within",
                          nomatch = NA)

# Check
cat("Mutations with RPKM score:", sum(!is.na(ctDNA_scored$RPKM_score)), "\n")
cat("Mutations without score:", sum(is.na(ctDNA_scored$RPKM_score)), "\n")

# -----------------------------------------------------------------------------
# Filter per cancer type
# -----------------------------------------------------------------------------

# Remove mutations without a score
ctDNA_scored_clean <- ctDNA_scored[!is.na(RPKM_score)]

# Filter by cancer type
ctDNA_LUAD <- ctDNA_scored_clean[histology1_group_central.reviewed == "LUAD"]
ctDNA_LUSC <- ctDNA_scored_clean[histology1_group_central.reviewed == "LUSC"]

cat("LUAD mutations scored:", nrow(ctDNA_LUAD), "\n")
cat("LUSC mutations scored:", nrow(ctDNA_LUSC), "\n")

# -----------------------------------------------------------------------------
# Check correlation between RPKM score and CCF z-score and plot
# -----------------------------------------------------------------------------

# Spearman correlation - all samples
cor_all <- cor.test(ctDNA_scored_clean$ccf_z_score, ctDNA_scored_clean$RPKM_score, 
                    method = "spearman")
print(cor_all)

# Spearman correlation - LUAD only
cor_LUAD <- cor.test(ctDNA_LUAD$ccf_z_score, ctDNA_LUAD$RPKM_score, 
                     method = "spearman")
print(cor_LUAD)

# Scatter plot - all samples
p_all <- ggplot(ctDNA_scored_clean, aes(x = RPKM_score, y = ccf_z_score)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", colour = "red") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.05)) +
  annotate("text", 
           x = Inf, y = Inf, 
           hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_all$estimate, 3), 
                          "\np = ", signif(cor_all$p.value, 3)),
           size = 4) +
  labs(x = "cfChromatin RPKM score (A549, 10kb bins)",
       y = "CCF z-score",
       title = "cfChromatin occupancy vs ctDNA shedding (all)") +
  coord_cartesian(clip = "off") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12))

# Scatter plot - LUAD only
p_LUAD <- ggplot(ctDNA_LUAD, aes(x = RPKM_score, y = ccf_z_score)) +
  geom_point(alpha = 0.1, size = 0.5) +
  geom_smooth(method = "lm", colour = "red") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), limits = c(0, 1.05)) +
  annotate("text", 
           x = Inf, y = Inf, 
           hjust = 1.1, vjust = 1.5,
           label = paste0("rho = ", round(cor_LUAD$estimate, 3), 
                          "\np = ", signif(cor_LUAD$p.value, 3)),
           size = 4) +
  labs(x = "cfChromatin RPKM score (A549, 10kb bins)",
       y = "CCF z-score",
       title = "cfChromatin occupancy vs ctDNA shedding (LUAD)") +
  theme_cowplot() +
  theme(plot.title = element_text(size = 12))

# Save plots
ggsave(paste0(outputs.folder, "cfChromatin_vs_shedding_all.pdf"), p_all, width = 6, height = 5)
ggsave(paste0(outputs.folder, "cfChromatin_vs_shedding_LUAD.pdf"), p_LUAD, width = 6, height = 5)









