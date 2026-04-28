#==============================================================================#
#==============================================================================#
######                                                                    ######
######  Evaluate CCFs of high/low shed mutations in three samples         ######
######                                                                    ######
######                                                                    ######
#==============================================================================#
#==============================================================================#

# Author: Hannah Bazin
# Date: 2026-04-20

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

analysis_name <- 'mets_CCF'
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

# Read in mutation table
mut_table_path <- "../../../Somatic-Evolutionary-Monitoring-Lab/personalis_ctDNA_mets_analysis/inputs/20240724_primary_met_and_peace_data_mutTable_tree_altered.fst"
mut_table <- read_fst(mut_table_path)

#####################################################################################
#### Analysis: do the 3 patients with high/low shed mutations have mets samples? ####
#####################################################################################

print(mut_table[mut_table$patient_id == "LTX208"])
print(mut_table[mut_table$patient_id == "LTX287"])
print(mut_table[mut_table$patient_id == "LTX854"])

































