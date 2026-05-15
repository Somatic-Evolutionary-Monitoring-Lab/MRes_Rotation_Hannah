# =============================================================================
# Custom ggplot2 theme and colour scheme
# Author: Hannah Bazin
# =============================================================================

# -----------------------------------------------------------------------------
# Colour scheme
# -----------------------------------------------------------------------------

wt_col         <- "#8E9BAE"
mut_col        <- "#E8829A"
scatter_col    <- "#7EB5D6"
high_col       <- "#82CEC9"
low_col        <- "#B8A9E3"
horiz_line_col <- "#D4DCE8"

# -----------------------------------------------------------------------------
# Custom theme
# -----------------------------------------------------------------------------

theme_mres_frankell <- function(base_size = 14) {
  theme_cowplot(font_size = base_size) %+replace%
    theme(
      # Title
      plot.title        = element_text(face = "bold", size = 12, colour = "#2D3436"),
      # Axis text and labels
      axis.text         = element_text(colour = "#2D3436", size = 12),
      axis.title        = element_text(colour = "#2D3436", size = 14),
      # Ticks
      axis.ticks        = element_line(colour = "#2D3436", linewidth = 0.8),
      axis.ticks.length = unit(4, "pt"),
      # Legend
      legend.text       = element_text(colour = "#2D3436"),
      legend.title      = element_text(colour = "#2D3436"),
      # General text
      text              = element_text(colour = "#2D3436")
    )
}