import pandas as pd
import numpy as np
import os

# -----------------------------------------------------------------------------
# Load mutation table
# -----------------------------------------------------------------------------

DATA_DIR   = "/SAN/colcc/tracerx_personalis_pipeline/hannah/MRes_Rotation_Hannah/01.personalis/00.data"
OUTPUT_DIR = "/SAN/colcc/tracerx_personalis_pipeline/hannah/MRes_Rotation_Hannah/01.personalis/01.shedding/00.fragmentomics"
BED_DIR     = os.path.join(DATA_DIR, "00.bed_files")
MUT_DIR     = os.path.join(DATA_DIR, "01.mutation_table")

# Set to false to run on all data, true to run on a single patient
TEST_RUN = False

ctDNA_mutations = pd.read_csv(
    os.path.join(MUT_DIR, "ctDNA_data_black_hg38_ccfzscore_highlow_20260507.tsv"),
    sep="\t"
)

print("Loaded ctDNA_mutations")

# -----------------------------------------------------------------------------
# Compute mutation- and patient-level fragmentomics characteristics
# -----------------------------------------------------------------------------

# If test run, subset only to one patient
if TEST_RUN:
    ctDNA_mutations = ctDNA_mutations[ctDNA_mutations["cruk_id"] == "CRUK0080"]

# Initialise patient level results
patient_level_results = []

for cruk_id, patient_mutations in ctDNA_mutations.groupby("cruk_id"):
    
    # Initialise counters
    n_mutations_processed = 0
    n_mutations_no_overlap = 0
    n_mutations_no_mut_frags = 0
    n_mutations_no_wt_frags = 0

    bed_path = os.path.join(BED_DIR, f"cfDNA_{cruk_id}.fragments.bed")

    # Skip if no BED file for this patient
    if not os.path.exists(bed_path):
        print(f"No BED file for {cruk_id}, skipping")
        continue
    
    print(f"Loading {cruk_id}: {len(patient_mutations)} mutations")
    
    # Load BED file with column names
    bed = pd.read_csv(bed_path, sep="\t", header=None,
                      names=["chrom", "frag_start", "frag_end", "frag_length",
                             "cigar", "mapq", "start_motif", "end_motif",
                             "is_mutant", "mut_position", "ref_bed", "alt_bed", "sample_id"],
                             dtype=str)
    
    # Convert columns to correct types after loading
    bed["frag_start"]  = bed["frag_start"].astype(int)
    bed["frag_end"]    = bed["frag_end"].astype(int)
    bed["frag_length"] = bed["frag_length"].astype(int)
    bed["is_mutant"]   = bed["is_mutant"] == "True"


    # For each patient, this will be a table with summary of fragment characteristics per wt/mut status
    mutation_level_results = []

    for _, mutation in patient_mutations.iterrows():

        # Filter bed to relevant chromosome
        bed_chr = bed[bed["chrom"] == mutation["chrom"]]

        # Find all fragments overlapping this mutation position
        # Because bed is 0-based, the mutation position is "end" for SNVs (which these are)
        overlapping = bed_chr[
            (bed_chr["frag_start"] < mutation["end"]) &
            (bed_chr["frag_end"] >= mutation["end"])
        ]

        # Skip if no overlapping fragments found
        if len(overlapping) == 0:
            n_mutations_no_overlap += 1
            continue

        # Compute distances from fragment ends to mutation position
        overlapping = overlapping.copy()
        overlapping["distance_left"] = mutation["end"] - overlapping["frag_start"]
        overlapping["distance_right"] = overlapping["frag_end"] - mutation["end"]
        # Take the min distance as the measure for "fragment end proximity"
        overlapping["min_distance"] = overlapping[["distance_left", "distance_right"]].min(axis=1)

        # Split into wt and mutant fragments
        wt_frags  = overlapping[overlapping["is_mutant"] == False]
        mut_frags = overlapping[overlapping["is_mutant"] == True]
        
        # Skip if no wt fragments found
        if len(wt_frags) == 0:
            n_mutations_no_wt_frags += 1
            continue
        
        # Skip if no mutant fragments found
        if len(mut_frags) == 0:
            n_mutations_no_mut_frags += 1
            continue
        
        # Increase counter for processed mutations (that had overlapping fragments and wt fragments and mut fragments)
        n_mutations_processed += 1

        # Save summary for this mutation
        mutation_level_results.append({
            "cruk_id":                   cruk_id,
            "tracerx_id":                mutation["patient"],
            "sample":                    mutation["sample"],
            "chrom":                     mutation["chrom"],
            "pos":                       mutation["end"],
            "ref":                       mutation["ref"],
            "alt":                       mutation["alt"],
            "ccf_z_score":               mutation["ccf_z_score"],
            "ccf_z_score_group":         mutation["ccf_z_score_group"],
            "nb_wt":                     len(wt_frags),
            "nb_mut":                    len(mut_frags),
            "median_frag_length_wt":     np.median(wt_frags["frag_length"])     if len(wt_frags)  > 0 else np.nan,
            "median_frag_length_mut":    np.median(mut_frags["frag_length"])    if len(mut_frags) > 0 else np.nan,
            "median_distance_left_wt":   np.median(wt_frags["distance_left"])   if len(wt_frags)  > 0 else np.nan,
            "median_distance_left_mut":  np.median(mut_frags["distance_left"])  if len(mut_frags) > 0 else np.nan,
            "median_distance_right_wt":  np.median(wt_frags["distance_right"])  if len(wt_frags)  > 0 else np.nan,
            "median_distance_right_mut": np.median(mut_frags["distance_right"]) if len(mut_frags) > 0 else np.nan,
            "median_min_distance_wt":    np.median(wt_frags["min_distance"])    if len(wt_frags)  > 0 else np.nan,
            "median_min_distance_mut":   np.median(mut_frags["min_distance"])   if len(mut_frags) > 0 else np.nan
        })

    # Save the mutation table for this patient
    mutation_level_results_df = pd.DataFrame(mutation_level_results)
    mutation_level_results_df.to_csv(
        os.path.join(OUTPUT_DIR, f"{cruk_id}_per_mut_fragmentomics.tsv"),
        sep="\t", index=False
    )
    print(f"Saved {cruk_id}: {n_mutations_processed} mutations processed, "
        f"{n_mutations_no_overlap} skipped (no overlapping fragments), "
        f"{n_mutations_no_wt_frags} skipped (no WT fragments), "
        f"{n_mutations_no_mut_frags} skipped (no mutant fragments)")
    
    # Compute patient level fragmentomics characteristics
    # Filter by CCF group
    high_muts = mutation_level_results_df[mutation_level_results_df["ccf_z_score_group"] == "high"]
    low_muts  = mutation_level_results_df[mutation_level_results_df["ccf_z_score_group"] == "low"]

    # Compute patient level fragmentomics characteristics
    patient_level_results.append({
        "cruk_id":                       cruk_id,
        "n_mutations_high":              len(high_muts),
        "n_mutations_low":               len(low_muts),
        "median_frag_length_wt_high":    np.median(high_muts["median_frag_length_wt"])    if len(high_muts) > 0 else np.nan,
        "median_frag_length_mut_high":   np.median(high_muts["median_frag_length_mut"])   if len(high_muts) > 0 else np.nan,
        "median_frag_length_wt_low":     np.median(low_muts["median_frag_length_wt"])     if len(low_muts)  > 0 else np.nan,
        "median_frag_length_mut_low":    np.median(low_muts["median_frag_length_mut"])    if len(low_muts)  > 0 else np.nan,
        "delta_frag_length_high":        np.median(high_muts["median_frag_length_mut"])   - np.median(high_muts["median_frag_length_wt"])  if len(high_muts) > 0 else np.nan,
        "delta_frag_length_low":         np.median(low_muts["median_frag_length_mut"])    - np.median(low_muts["median_frag_length_wt"])   if len(low_muts)  > 0 else np.nan,
        "median_dist_left_wt_high":      np.median(high_muts["median_distance_left_wt"])  if len(high_muts) > 0 else np.nan,
        "median_dist_left_mut_high":     np.median(high_muts["median_distance_left_mut"]) if len(high_muts) > 0 else np.nan,
        "median_dist_left_wt_low":       np.median(low_muts["median_distance_left_wt"])   if len(low_muts)  > 0 else np.nan,
        "median_dist_left_mut_low":      np.median(low_muts["median_distance_left_mut"])  if len(low_muts)  > 0 else np.nan,
        "delta_dist_left_high":          np.median(high_muts["median_distance_left_mut"]) - np.median(high_muts["median_distance_left_wt"]) if len(high_muts) > 0 else np.nan,
        "delta_dist_left_low":           np.median(low_muts["median_distance_left_mut"])  - np.median(low_muts["median_distance_left_wt"])  if len(low_muts)  > 0 else np.nan,
        "median_dist_right_wt_high":     np.median(high_muts["median_distance_right_wt"]) if len(high_muts) > 0 else np.nan,
        "median_dist_right_mut_high":    np.median(high_muts["median_distance_right_mut"])if len(high_muts) > 0 else np.nan,
        "median_dist_right_wt_low":      np.median(low_muts["median_distance_right_wt"])  if len(low_muts)  > 0 else np.nan,
        "median_dist_right_mut_low":     np.median(low_muts["median_distance_right_mut"]) if len(low_muts)  > 0 else np.nan,
        "delta_dist_right_high":         np.median(high_muts["median_distance_right_mut"])- np.median(high_muts["median_distance_right_wt"]) if len(high_muts) > 0 else np.nan,
        "delta_dist_right_low":          np.median(low_muts["median_distance_right_mut"]) - np.median(low_muts["median_distance_right_wt"])  if len(low_muts)  > 0 else np.nan,
    })

# Save the patient-level characteristics
patient_level_results_df = pd.DataFrame(patient_level_results)
patient_level_results_df.to_csv(
    os.path.join(OUTPUT_DIR, "patient_level_fragmentomics.tsv"),
    sep="\t", index=False
)
print(f"Saved patient level results: {len(patient_level_results_df)} patients")

