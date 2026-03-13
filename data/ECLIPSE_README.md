DATE_tracerx_variant_table.tsv – this is all the primary mutation data across the TRACERx cohort excluding the PEACE metastatic mutations.
DATE_peace_variant_table.tsv – this is the PEACE metastatic and primary tumor data combined for LTX208/392 and 474 – it appears it was LTX568 where only primary tumor mutations were tracked but please confirm by looking at the first file. I have that we tracked 641 mutations in LTX474 (this is the interesting case we could publish in the MRD paper), 384 mutations in LTX392 and 485 in LTX208.
 
Whether an MRD call was made and the date data for the above tables are in the sample level databases below:
 
DATE_tracerx_sample_table.tsv
DATE_peace_sample_table.tsv
 
Here are some of the important headers for the variant level table:
 
Tracerx_id – this is the unique identifier for that plasma timepoint for the patient, present in both sample and variant tables (date and day post op data in sample table)
Dao – deep alternate observation, number of deep alternate observations for the variant (deep means supported by at least 5 duplicate reads with same molecular barcode).
Dro – same as above but reference
Ddp – deep depth = total error corrected coverage of position.
Daf – deep allele frequency

Tnc_error_rate – In that library this is the error-rate associated with the trinucleotide context of the variant being evaluated so the rate at which we expect to see error DAOs.
D_group_error_rate – This is used for the MRD caller if we cluster all the TNC error rates into 4 categories this is the clustered error-rate, a better estimate that the individual TNC_error_rates as more data goes into this parameter.
Mnp_error_rate – For DNVs use this column for error rates.
Adjusted_variant_poisson_test  - I use an alpha of <= 0.01 in this column to say whether DAOs associated with a position are real or artefactual, this P-val is basically asking whether the number of DAOs observed deviates from the expected based on error and takes into account the size of the panel.
Failed filters – I recommend filtering the variants based on this column being empty e.g., == “”, or containing “tnc_error_rate” (the latter is only relevant for the MRD caller).
Variant_chooser – If ADX – Archer selected this variant and it may not have been called by our pipeline.
 
In the sample level database there is a column matches_wes – make sure this is TRUE to remove the small number of plasma timepoints that were sample swaps – this is a SNP ID check.


The DATE1_DATE2_tracerx_variant_table_annotated.tsv file contains additional columns with data from the TRACERx exome pipeline including Clonaltiy status, neoantigen status, mutation expression and others. DATE1 indciate the date annotation took place and DATE2 is the date for the ctDNA data itself (ie the same as 'DATE' above). 

In general equivilent columns have the name as in the mutation table but there are a few others:
is_nag - Indicates if the mutation is either a 'high binding' or 'strong' neoantigen in the nag output. Many mutations are NA as only missense mutations are considered in our nag prediction. 
expression_index - This metric quatifies a mutation's expression in regions where it is present accounting for the cellularity and CCF of the mutation in the primary tumour. 
missing_data - indicates patients which have not been included in the latest pipeline uns for the 421 as yet so are missing many annotation. These should be included soon. 
no_longer_called - A minority of mutations (~2%) were called by older verisons of the TRACERx pipeline and hence tracked in our ctDNA panels but are no longer called in the latest version of the pipeline. These will also be NA for many annotations. 

The script used for aming this annotated file is here: /camp/project/proj-tracerx-lung/tctProjects/frankella/archer_ctdna/Scripts/pipeline/overlay_tx_data.R


ECLIPSE columns:


group: Within each clone background noise reads are estimated from within these groups. This is copying the MRD D groups but the number of groups is variable depending on the number of mutations / overall amount of background niose to ensure that background niose is accurately estimated. These are groups are determined for each clone using kmeans clustering of the background noise and the group number proportional to the total amount of noise summed across all variants. 

supporting_reads_no_niose: supporting deep reads (varcount/daos) but with supporting reads assigned to noise removed

vaf_no_background: deep vaf (daf) but using supporting_reads_no_niose as varcount

mutCPN: the average num of copies across all tumour cells in the primary (will be < 1 (the CCF) in subclonal mutations with 1 copy)

multiplicity: the number of mutation copies in mutated cells (ie cannot be < 1 and is in reality an integer if one 1 cn state for each mutation)

int_multiplicity: multiplicity rounded to integer

purity: The % tumour content of the sample (can be used as ctDNA fraction adjusted for CN, background noise, 
LOD variants and outlier variants). 

ccf: calculated ccf of each mutation using the NEJM 2017 equations between vaf, multiplicity, cn & purity

clone_purity_mutation: This is the ccf of the mutation * purity - should be directly preportional the number of cells habouring this mutation in the body 

mulit_modal_p: p value for if the clone is multimodal indicting that a large number of copy number alterations have occured since we measured them in the primary tumour. We will only be sensitive to this at relatively high purities (>1% most likely)

subsequecent_cin: is mulit_modal_p significant (< 0.01)

int_multiplicity_preCIN: the old copy number (before CIN) found in the primary tumour - NA if no subsequent cin

subsequent_amplfiication: has this mutation been amplified (> CN 2) in a subsequent CIN event (can’t tell exact CN in these cases) - NA if no subsequent cin

cn_change: change in CN from when measured at op to when subsequent CIN detected in ctDNA

new_clone_ccf: in some cases you can determine the CCF of the new clone with the subsequent CIN if there are 
completely deleted mutations in the new clone which are still present in the parent clones
clone_purity_mutation_LOD: a measure of the limit of detection for clone_purity_mutation (ie cn_ajd_vaf if varcount is 1)

ccf_LOD = A quantification of thw limit of detection for ccf in this mutaiton (ie ccf if varcount was 1)

clonal_purity_mut_LOD = A quantification of thw limit of detection for clonal purity in this mutaiton (ie ccf if varcount was 1)

matched_lod: FALSE flags mutations which have varcount of 0 because of a low LOD (eg low depth) compared to the variants with read support. These are exlcuded from some calculationes eg clone ccf calculations

zscore: Z score for good quality mutations (not hard filtered and matched_lod = TRUE) in clone

clone_normsd: normalied sd (sd / mean) for the clone after excluding bad quality mutations and outliers

is_outlier: identifier outliers - these are identified and excluded until the normalied SD fo the clone is in an acceptable range (based on high fraction clonal clones)

poor_quality_clone: Is the clone trustworthy with a coherent vaf distribution? TRUE if too many outliers have to be excluded to get to an acceptable normalised SD

clone_ccf: ccf of a clone (mean ccf of quality mutations)

clone_purity: clone_ccf * purity - hould be directly preportional the number of cells from this clone in the body 

clone_present_p: p value for a if a clone is present. Uses possion test comparing background noise (lambda) to the total number of observations. Only includes low niose variants. Should ese p = 0.01 or can do 
multiple hypothesis correction across a specific cohort. 

power_clone_ccf: The equivalent ccf if we observed enough reads to detect a clone a p = 0.01 using poisson test in this clone (based on depth background noise etc in the mutations in this clone)

power_sample_ccf: The equivalent ccf if we observed enough reads to detect a 5 mutation subclone in a given sample at p = 0.01 using poisson test (based on depth background noise etc in the mutations int his clone)

power_purity: The equivalent purity that would be detectable at p = 0.01 in this sample using the clonal 
mutations (similar to MDAF but useful for tool)

is_subclonal_sample_p:  Does this clone have significantly lower CCFs than the clonal cluster? If not we can consider it probably clonal in this sample (/pseudo clonal) - perhaps should also add an absolute clone_ccf theshold as well (0.9 or 0.95) to exclude clones with a lot of uncertainty in ccf and few mutations. 

is_subclonal_sample: Is the clone subclonal in this sample ie with is_subclonal_sample_p < 0.01









