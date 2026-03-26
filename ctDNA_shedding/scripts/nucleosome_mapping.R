#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Map high/low shedding mutations to nucleosome position map        ######
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

analysis_name <- 'nucleosome_mapping'
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

# Read in outlier mutations (high/low shedding)
outlier_mutations_path <- "outputs/ccf_zscores_multi_patient/20260316/outlier_mutations.csv"
outlier_mutations <- read.csv(outlier_mutations_path)

# Parse Pos into chr and position
outlier_mutations <- outlier_mutations %>% 
  tidyr::separate(Pos, into = c("chr_num", "pos", "ref", "alt"), sep = ":", remove = FALSE) %>% 
  mutate(
    chr = paste0("chr", chr_num),
    pos = as.numeric(pos)
  ) %>%
  select(-chr_num) %>% 
  select(patient_name, Pos, chr, pos, ref, alt, n_samples, mean_z, shedding_status)
head(outlier_mutations)

# Read in nucleosome map
nucleosome_map_path <- "data/SRA438908_lung_cancer_Ma2017_stable_100bp_hg19.bed"
nuc <- data.table::fread(nucleosome_map_path)
colnames(nuc) <- c("chr", "start", "end", "nuc_occupancy", "sd", "rel_deviation")
head(nuc)

##########################################################################
#### Do highly shed mutation have higher nucleosome occupancy scores? ####
##########################################################################

# For each mutation, get the nucleosome occupancy of the overlapping bin (NA if no overlap)
outlier_mutations <- outlier_mutations %>%
  mutate(nuc_occupancy = mapply(function(c, p) {
    idx <- which(nuc$chr == c & nuc$start <= p & nuc$end > p)
    if (length(idx) > 0) nuc$nuc_occupancy[idx[1]] else NA
  }, chr, pos))

# Compare occupancy scores between high and low shedders
outlier_mutations %>%
  filter(!is.na(nuc_occupancy)) %>%
  group_by(shedding_status) %>%
  summarise(
    n = n(),
    mean_occ = mean(nuc_occupancy),
    median_occ = median(nuc_occupancy)
  )

# Wilcoxon test
high_occ <- outlier_mutations %>% filter(shedding_status == "High", !is.na(nuc_occupancy)) %>% pull(nuc_occupancy)
low_occ <- outlier_mutations %>% filter(shedding_status == "Low", !is.na(nuc_occupancy)) %>% pull(nuc_occupancy)
wt <- wilcox.test(high_occ, low_occ)

# Boxplot
outlier_mutations %>%
  filter(!is.na(nuc_occupancy)) %>%
  ggplot(aes(x = shedding_status, y = nuc_occupancy, fill = shedding_status)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.5) +
  scale_fill_manual(values = c("High" = "#D7191C", "Low" = "#2C7BB6")) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey40") +
  theme_cowplot() +
  annotate("text", x = 1.5, y = max(c(high_occ, low_occ)) + 0.1, 
           label = paste0("Wilcoxon p = ", round(wt$p.value, 3)), size = 4) +
  labs(x = "Shedding status", y = "Normalised nucleosome occupancy",
       title = "Nucleosome occupancy at high vs low shedding mutations") +
  theme(legend.position = "none")

ggsave(paste0(outputs.folder, "nuc_occupancy_high_vs_low_shed.pdf"), width = 8, height = 6)

#############################################################################################################
#### Do lowly shed mutations fall in unstable nucleosome regions, i.e. those without an occupancy score? ####
#############################################################################################################

# Fisher test
ft <- fisher.test(table(outlier_mutations$shedding_status, !is.na(outlier_mutations$nuc_occupancy)))

# Prepare data for plot
stable_counts <- outlier_mutations %>%
  mutate(region = ifelse(!is.na(nuc_occupancy), "Stable nucleosome", "Unstable region")) %>%
  group_by(shedding_status, region) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(shedding_status) %>%
  mutate(pct = 100 * n / sum(n))

# Stacked bar plot
ggplot(stable_counts, aes(x = shedding_status, y = pct, fill = region)) +
  geom_col(alpha = 0.8) +
  geom_text(aes(label = paste0(n, " (", round(pct, 1), "%)")), 
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_fill_manual(values = c("Stable nucleosome" = "#4DAF4A", "Unstable region" = "#984EA3")) +
  annotate("text", x = 1.5, y = 105, 
           label = paste0("Fisher p = ", round(ft$p.value, 3)), size = 4) +
  theme_cowplot() +
  labs(x = "Shedding status", y = "Percentage of mutations (%)",
       title = "Mutations in stable vs unstable nucleosome regions",
       fill = "Region") +
  ylim(0, 110)

ggsave(paste0(outputs.folder, "nuc_stable_vs_unstable_regions.pdf"), width = 8, height = 6)






































