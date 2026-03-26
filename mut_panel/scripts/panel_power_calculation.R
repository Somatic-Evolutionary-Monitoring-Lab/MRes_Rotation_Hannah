#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Estimate power depending on bin size and mutation rates           ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-03-24

setwd("/Volumes/RFS/rfs-kh_rfs-rDsHEAv2WP0/hannah/MRes_Rotation_Hannah/mut_panel/")

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
library(pwrss)

#############################################
#### Make a folder for this analysis run ####
#############################################

date <- gsub("-","",Sys.Date())

analysis_name <- 'power_calc'
out_name <- 'outputs'
out_dir_general <- paste(out_name, analysis_name, sep='/')
if( !file.exists(out_dir_general) ) dir.create( out_dir_general )

out_dir_logs <- paste(out_dir_general, 'logs', sep='/')
if( !file.exists(out_dir_logs) ) dir.create( out_dir_logs )

outputs.folder <- paste0( out_dir_general, "/", date, "/" )

if( !file.exists(outputs.folder) ) dir.create( outputs.folder )


##########################################################
#### Plot probability of at least one mutation        ####
##########################################################

# ============================================================
# Background hit probability as a function of bin size
# P(X >= 1) = 1 - exp(-lambda * B)
# lambda = 2 mutations/Mb = 2e-6 mutations/bp (median coding TMB proxy)
# ============================================================

lambda <- 2e-6  # mutations per bp

bin_sizes <- seq(100, 500000, by=100)  # 100bp to 500kb

p_background <- 1 - exp(-lambda * bin_sizes)

pdf(paste0(outputs.folder, "background_hit_probability.pdf"), width=9, height=6)

plot(bin_sizes / 1000, p_background * 100,
     type = "l", lwd = 2, col = "steelblue",
     log = "x",  # log scale on x-axis
     xlab = "Bin size (kb, log scale)",
     ylab = "Background hit probability (%)",
     main = "Expected background mutation probability per bin\n(lambda = 2 mutations/Mb, Poisson model)")

candidate_bins <- c(1000, 5000, 10000, 50000, 100000)
candidate_p <- 1 - exp(-lambda * candidate_bins)

points(candidate_bins / 1000, candidate_p * 100,
       pch = 19, col = "red", cex = 1.2)

label_pos <- c(1, 1, 1, 1, 1)  # all to the right
label_offset <- c(-0.8, -0.8, -0.8, -0.5, -0.5)  # vertical offset

text(candidate_bins / 1000, candidate_p * 100 + label_offset,
     labels = paste0(candidate_bins/1000, "kb (", round(candidate_p*100, 2), "%)"),
     pos = 4, cex = 0.8, col = "red")

abline(h = c(1, 5, 10), lty = 2, col = "grey60")
text(rep(0.12, 3), c(1.4, 5.4, 10.4), 
     labels = c("1%", "5%", "10%"),
     cex = 0.75, col = "grey40")

dev.off()

cat("Plot saved to", paste0(outputs.folder, "background_hit_probability.pdf"), "\n")

##########################################################
#### Power calculation                                ####
##########################################################

# Training set size = 80% of GEL cohort size
N <- round(0.8 * 16341)

# Mutation rate
lambda <- 2e-6

# Approximate human genome size
genome_size <- 3e9

# Significance level before Bonferroni correction
alpha <- 0.05

# Define candidate bin sizes
bin_sizes <- c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000)

f_min_values <- numeric(length(bin_sizes))
p_bg_values <- numeric(length(bin_sizes))

for (i in seq_along(bin_sizes)) {
  B <- bin_sizes[i]
  n_bins <- genome_size / B
  alpha_bonf <- alpha / n_bins
  p_bg <- 1 - exp(-lambda * B)
    
  f_min <- uniroot(function(f) {
    power.z.oneprop(prob = f,
                    null.prob = p_bg,
                    n = N,
                    alpha = alpha_bonf,
                    alternative = "one.sided",
                    verbose = FALSE)$power - 0.8
  }, interval = c(p_bg + 1e-6, 0.99))$root
  
  cat(sprintf("Bin: %6d bp | p_bg: %.4f | min detectable freq: %.3f%% | ~%d patients\n",
              B, p_bg, f_min * 100, round(f_min * N)))
  
  p_bg_values[i] <- p_bg
  f_min_values[i] <- f_min
}

# Plotting
results <- data.frame(
  bin_size_kb = bin_sizes / 1000,
  p_bg = p_bg_values,
  f_min_pct = f_min_values * 100,
  n_patients = round(f_min_values * N)
)

ggplot(results, aes(x = bin_size_kb, y = f_min_pct)) +
  geom_line() +
  geom_point(size = 2) +
  geom_text(aes(label = paste0("n=", n_patients)),
            vjust = -1.2, hjust = 0.9, size = 2.8) +
  scale_x_log10() +
  theme_cowplot() +
  labs(title = "Minimum detectable mutation frequency by bin size",
       subtitle = paste0("80% power, Bonferroni correction, N=", format(N, big.mark=","), ", lambda=", lambda*1e6, "/Mb"),
       x = "Bin size (kb, log scale)",
       y = "Minimum detectable frequency (%)")


ggsave(paste0(outputs.folder, "power_analysis_min_detectable_freq.pdf"),
       width = 9, height = 6)













