#==============================================================================#
#                                                                              #
#          Merge cfChromatin A549 Replicates - Average RPKM Scores            #
#                                                                              #
#  Description: Combines two biological replicates of A549 cfMNase-Seq data   #
#               by averaging RPKM scores per 10kb genomic bin (hg19).         #
#                                                                              #
#  Input:  GSM9265838 - A549 Rep1 cfChromatin, cfMNase-Seq, 30min digestion   #
#          GSM9265839 - A549 Rep2 cfChromatin, cfMNase-Seq, 30min digestion   #
#  Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9265838      #
#          https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM9265839      #
#                                                                              #
#  Output: A549_cfChromatin_merged_hg19.bedGraph                              #
#                                                                              #
#                                                                              #
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-04-28

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")

# -----------------------------------------------------------------------------
# Load libraries
# -----------------------------------------------------------------------------

library(data.table)

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

# Load both replicates
rep1 <- fread("data/GSM9265838_A549_rep1_hg19.bedGraph", 
              col.names = c("chr", "start", "end", "score", "extra"))
rep2 <- fread("data/GSM9265839_A549_rep2_hg19.bedGraph", 
              col.names = c("chr", "start", "end", "score", "extra"))

# Drop the extra column
rep1[, extra := NULL]
rep2[, extra := NULL]

# -----------------------------------------------------------------------------
# Merge two replicates
# -----------------------------------------------------------------------------

# Merge on chr, start, end
merged <- merge(rep1, rep2, by = c("chr", "start", "end"), suffixes = c("_rep1", "_rep2"))

# Average the scores
merged[, RPKM_score_mean := (score_rep1 + score_rep2) / 2]

# Keep only necessary columns
result <- merged[, .(chr, start, end, RPKM_score_mean)]

# Sort
setorder(result, chr, start)

# Save
fwrite(result, "data/A549_cfChromatin_merged_hg19.bedGraph", 
       sep = "\t", col.names = FALSE)

cat("Done! Bins in merged file:", nrow(result), "\n")
cat("Bins only in rep1:", nrow(rep1) - nrow(result), "\n")
cat("Bins only in rep2:", nrow(rep2) - nrow(result), "\n")












