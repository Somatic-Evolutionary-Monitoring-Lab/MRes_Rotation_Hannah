#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Explore CCF distributions across different mutations to identify  ######
######  high / low shedding mutations                                     ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-03-16

setwd("/Volumes/proj-tracerx-lung/tctProjects/frankella/archer_ctdna/")

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

#############################################
#### Make a folder for this analysis run ####
#############################################