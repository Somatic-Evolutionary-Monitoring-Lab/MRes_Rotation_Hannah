# =============================================================================
# Reformat TCGA ATAC-seq peak calls for liftover (hg38 -> hg19)
# Input:  LUAD and LUSC peak call txt files (hg38, from Corces et al. 2018)
# Output: Trimmed BED files (chr, start, end only) ready for UCSC LiftOver
# =============================================================================

library(dplyr)

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/ctDNA_shedding/")

# -----------------------------------------------------------------------------
# LUAD
# -----------------------------------------------------------------------------

luad_peaks <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUAD_peakCalls.txt",
  header = TRUE,
  sep = "\t"
)

head(luad_peaks)

# Trim to BED format (seqnames, start, end only)
luad_peaks_reformat <- luad_peaks %>%
  select(seqnames, start, end)

# Write out as BED file (no header, no row names)
write.table(
  luad_peaks_reformat,
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUAD_peakCalls_hg38.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

# -----------------------------------------------------------------------------
# LUSC
# -----------------------------------------------------------------------------

lusc_peaks <- read.table(
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUSC_peakCalls.txt",
  header = TRUE,
  sep = "\t"
)

head(lusc_peaks)

# Trim to BED format (seqnames, start, end only)
lusc_peaks_reformat <- lusc_peaks %>%
  select(seqnames, start, end)

# Write out as BED file (no header, no row names)
write.table(
  lusc_peaks_reformat,
  "data/TCGA-ATAC_Cancer_Type-specific_PeakCalls/LUSC_peakCalls_hg38.bed",
  quote = FALSE,
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)