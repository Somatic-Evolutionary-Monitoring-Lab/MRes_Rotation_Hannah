# ============================================================
# GEL Cancer Analysis Plots
# Author: Hannah Bazin
# Date: 2026-03-20
# Data: cancer_analysis_2026-03-20_13-47-35.csv
# ============================================================
# To run this script:
# module load R/4.3.3
# R
# source("~/hbazin_GEL_panel/gel_cancer_plots.R")

# ------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------

library(farver,       lib.loc="~/re_gecip/cancer_pan/R_packages")
library(labeling,     lib.loc="~/re_gecip/cancer_pan/R_packages")
library(withr,        lib.loc="~/re_gecip/cancer_pan/R_packages")
library(RColorBrewer, lib.loc="~/re_gecip/cancer_pan/R_packages")
library(cli,          lib.loc="~/re_gecip/cancer_pan/R_packages")
library(ggplot2,      lib.loc="~/re_gecip/cancer_pan/R_packages")

# Load data
df_raw <- read.csv("~/Downloads/cancer_analysis_2026-03-20_13-47-35.csv")
df <- df_raw

# Handle censored contamination column
# "<1.0" means known to be below 1% - passes 5% threshold
# Numeric values are parsed directly
raw_contam <- df_raw$Tumour.Sample.Cross.Contamination.Percentage
df$contam_passes_qc <- ifelse(
  raw_contam == "<1.0", TRUE,
  !is.na(suppressWarnings(as.numeric(raw_contam))) &
    suppressWarnings(as.numeric(raw_contam)) < 5
)
# For plotting only: numeric values where available, NA for <1%
df$contam_numeric <- suppressWarnings(as.numeric(raw_contam))

cat("Samples with <1.0% contamination:", sum(raw_contam == "<1.0"), "\n")
cat("Samples passing contamination QC:", sum(df$contam_passes_qc), "\n")

# Create output directory
out_dir <- "~/hbazin_GEL_panel"
dir.create(out_dir, showWarnings = FALSE)

# ------------------------------------------------------------
# Plot 1: Number of samples per cancer type
# ------------------------------------------------------------

pdf(file.path(out_dir, "cancer_type_counts.pdf"), width=12, height=6)
ggplot(df, aes(x = reorder(Disease.Type, Disease.Type, function(x) -length(x)))) +
  geom_bar(fill = "steelblue") +
  labs(title = "Number of samples per cancer type in GEL data",
       x = "Cancer type",
       y = "Number of samples") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))
dev.off()

# ------------------------------------------------------------
# Plot 2: TMB distribution per cancer type
# ------------------------------------------------------------

cancer_order <- names(sort(tapply(df$Somatic.Coding.Variants.Per.Mb,
                                  df$Disease.Type,
                                  median, na.rm=TRUE), decreasing=TRUE))
df_plot <- df
df_plot$Disease.Type <- factor(df_plot$Disease.Type, levels=cancer_order)

pdf(file.path(out_dir, "TMB_by_cancer_type.pdf"), width=14, height=7)
ggplot(df_plot, aes(x=Disease.Type, y=Somatic.Coding.Variants.Per.Mb, fill=Disease.Type)) +
  geom_violin(trim=FALSE, alpha=0.7) +
  geom_boxplot(width=0.1, outlier.size=0.5, alpha=0.9, fill="white") +
  scale_y_log10(
    breaks=c(0.01, 0.1, 1, 10, 100, 1000),
    labels=c("0.01", "0.1", "1", "10", "100", "1000")
  ) +
  labs(title="Tumour mutational burden distribution by cancer type",
       x="Cancer type",
       y="Somatic coding variants per Mb (log10)") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        legend.position="none")
dev.off()

# ------------------------------------------------------------
# QC Exploration
# ------------------------------------------------------------

cat("\nTotal samples:", nrow(df), "\n")
cat("Total columns:", ncol(df), "\n")

qc_cols <- c("Tumour.Purity", "Tumour.Autosomal.Coverage",
             "Germline.Autosomal.Coverage",
             "Somatic.Coding.Variants.Per.Mb",
             "Mapping.Rate", "Chimeric.Percentage",
             "Coverage.Homogeneity")

cat("\n--- Missing values in QC columns ---\n")
for (col in qc_cols) {
  n_missing <- sum(is.na(df[[col]]))
  cat(sprintf("%-55s %d missing (%.1f%%)\n", col, n_missing, 100*n_missing/nrow(df)))
}

cat("\n--- Summary statistics ---\n")
print(summary(df[, qc_cols]))

cat("\nSamples with Somatic.Small.Variants.Vcf.Path:",
    sum(!is.na(df$Somatic.Small.Variants.Vcf.Path) &
          df$Somatic.Small.Variants.Vcf.Path != ""), "\n")
cat("Samples with Small.Variants.Tiering.Path:",
    sum(!is.na(df$Small.Variants.Tiering.Path) &
          df$Small.Variants.Tiering.Path != ""), "\n")

# Purity distribution
cat("\nPurity-missing samples - coverage summary:\n")
print(summary(df$Tumour.Autosomal.Coverage[is.na(df$Tumour.Purity)]))
cat("Purity-missing samples - TMB summary:\n")
print(summary(df$Somatic.Coding.Variants.Per.Mb[is.na(df$Tumour.Purity)]))

pdf(file.path(out_dir, "QC_purity.pdf"), width=10, height=5)
par(mfrow=c(1,2))
hist(df$Tumour.Purity, breaks=50, main="Tumour purity distribution",
     xlab="Purity (%)", col="steelblue", border="white")
hist(df$Tumour.Purity[df$Tumour.Purity < 40], breaks=50,
     main="Zoom: purity < 40%",
     xlab="Purity (%)", col="steelblue", border="white")
dev.off()

# Coverage distribution
pdf(file.path(out_dir, "QC_coverage.pdf"), width=10, height=5)
par(mfrow=c(1,2))
hist(df$Tumour.Autosomal.Coverage, breaks=50,
     main="Tumour autosomal coverage",
     xlab="Mean coverage", col="darkgreen", border="white")
abline(v=30, col="red", lty=2, lwd=2)

hist(df$Germline.Autosomal.Coverage, breaks=50,
     main="Germline autosomal coverage",
     xlab="Mean coverage", col="darkgreen", border="white",
     xlim=c(0, max(df$Germline.Autosomal.Coverage, na.rm=TRUE)))
dev.off()

# Contamination distribution
cat("\nContamination: samples censored at <1.0%:", sum(raw_contam == "<1.0"), "\n")
cat("Contamination: samples with numeric value:", sum(!is.na(df$contam_numeric)), "\n")
cat("Contamination summary (numeric values only):\n")
print(summary(df$contam_numeric))

pdf(file.path(out_dir, "QC_contamination.pdf"), width=8, height=5)
hist(df$contam_numeric, breaks=50,
     main=paste0("Tumour cross-contamination\n(",
                 sum(raw_contam == "<1.0"),
                 " samples had contamination below the detection limit of 1%,\nnot shown)"),
     xlab="Contamination (%)", col="tomato", border="white")
abline(v=5, col="red", lty=2, lwd=2)
dev.off()

# ------------------------------------------------------------
# QC filter
# Coverage >= 30x and contamination < 5% applied
# ------------------------------------------------------------

df_pass <- df[
  !is.na(df$Tumour.Autosomal.Coverage) & df$Tumour.Autosomal.Coverage >= 30 &
    df$contam_passes_qc,
]

cat("\n--- QC filter summary ---\n")
cat("Samples before QC:", nrow(df), "\n")
cat("Samples after QC:", nrow(df_pass), "\n")
cat("Samples removed:", nrow(df) - nrow(df_pass), "\n")

cat("\nFails coverage < 30x:", sum(!is.na(df$Tumour.Autosomal.Coverage) &
                                     df$Tumour.Autosomal.Coverage < 30), "\n")
cat("Missing coverage:", sum(is.na(df$Tumour.Autosomal.Coverage)), "\n")
cat("Fails contamination QC:", sum(!df$contam_passes_qc), "\n")

cat("\n--- Post-QC samples per cancer type ---\n")
print(sort(table(df_pass$Disease.Type), decreasing=TRUE))

# Save QC-passing sample list for downstream analysis
write.csv(df_pass[, c("Participant.Id", "Tumour.Sample.Platekey",
                      "Disease.Type", "Somatic.Small.Variants.Vcf.Path",
                      "Small.Variants.Tiering.Path")],
          file.path(out_dir, "samples_post_QC.csv"),
          row.names=FALSE)
cat("\nSaved post-QC sample list to", file.path(out_dir, "samples_post_QC.csv"), "\n")


