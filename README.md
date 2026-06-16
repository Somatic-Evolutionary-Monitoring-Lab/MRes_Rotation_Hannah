# MRes Rotation — Determinants of ctDNA shedding & Tumour-naive panel design

This repository contains the analysis code for my MRes rotation project in the 
Frankell Lab (Early Cancer Institute, University of Cambridge). The project has 
two parts:

1. **Determinants of ctDNA shedding** — investigating why some mutations are 
   consistently over- or under-represented in plasma in the TRACERx NSCLC cohort
2. **Tumour-naive ctDNA panel design (ECLIPSE 2.0)** — designing a somatic
   mutation enrichment panel using Genomics England whole-genome sequencing data
   from ~13,500 tumour samples across 21 cancer types

---

## Part 1: Determinants of ctDNA shedding

Cancer cell fraction (CCF) z-scores were computed per mutation from ECLIPSE 
output to identify mutations that are consistently over- or under-represented 
in plasma relative to the sample mean. We investigated whether chromatin 
architecture and fragment-level properties of ctDNA explain this variability.

### `ctDNA_shedding/scripts/`
R scripts for chromatin architecture analyses, including:
- **Nucleosome positioning** — correlation of CCF z-scores with nucleosome 
  occupancy from NucPosDB (cfDNA and MNase-seq datasets)
- **Myeloid chromatin occupancy** — CML MNase-seq (K562, ENCSR000CXQ) 
  nucleosome occupancy at mutation loci
- **Chromatin accessibility** — ATAC-seq analyses using tumour 
  (LUAD/LUSC, Corces 2018) and blood (GSE74912, Corces 2016) peak calls
- **Clinical and technical features** — correlations with patient-level 
  shedding covariates
- **CCF z-score analysis** — computation and visualisation across the 
  Abbosh 2023 and Black 2025 TRACERx datasets

### `01.personalis/01.shedding/`
Python scripts for fragment-level analysis of ctDNA, including:
- **`shedding_fragmentomics_analysis.py`** — computes per-mutation and 
  patient-level fragment length and fragment end distance metrics for 
  mutant vs wild-type fragments
- **`shedding_fragmentomics_plots.ipynb`** — visualisation of fragmentomics 
  results
- **`submit_shedding_fragmentomics.sh`** — job submission script for the 
  UCL HPC cluster

---

## Part 2: Tumour-naive ctDNA panel design (ECLIPSE 2.0)

Somatic mutations from ~13,500 tumour samples across 21 cancer types were 
aggregated into 500bp genomic bins using somAgg WGS data from the 100,000 
Genomes Project (Genomics England, release 12). A panel of ~20,000 bins 
spanning 10Mb was selected by mutation density in a training set and validated 
on an independent test set by computing mutation enrichment relative to genomic
background.

> The code for this project cannot be made available as it resides within 
> the Genomics England Research Environment (GEL RE).

---

## Data availability

Input data and results are not included in this repository. Fragment-level BED 
files are stored on the UCL cluster (contact Woody for access). Mutation tables 
are derived from the TRACERx cohort (Abbosh et al. 2023, Black et al. 2025) 
and are not publicly available (contact Katie for access).
