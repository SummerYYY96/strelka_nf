import os
import glob
import pandas as pd

def find_matching_file(directory, id):
    # Use glob to search for files containing the ID
    print(directory)
    print(id)
    pattern = os.path.join(directory, f"*{id}*.bam")
    matching_files = glob.glob(pattern)
    print(matching_files)
    
    # Check if any files were found
    if matching_files:
        return matching_files[0]  # Return the first match (assuming there's a unique match)
    else:
        return None

# Load the CSV file
csv_file_path = "/gpfs/data/molecpathlab/development/chip_strelka_runs/data/AACR_953_age_final_65_plus.csv"
df = pd.read_csv(csv_file_path)
# Randomly select five rows
sampled_df = df.sample(n=5, random_state=42)  # random_state is set for reproducibility
# Create the "Sample" column
sampled_df["Sample"] = sampled_df["Tumor"] + "_" + sampled_df["Normal"]

# Base directory path
base_path = "output/alignments/deduplicated"

# Find matching BAM and BAI files
sampled_df["TumorBam"] = sampled_df.apply(lambda row: find_matching_file(os.path.join(row['Run'],base_path), row['Tumor']), axis=1)
sampled_df["TumorBai"] = sampled_df["TumorBam"].apply(lambda x: x + ".bai" if x else None)

sampled_df["NormalBam"] = sampled_df.apply(lambda row: find_matching_file(os.path.join(row['Run'],base_path), row['Normal']), axis=1)
sampled_df["NormalBai"] = sampled_df["NormalBam"].apply(lambda x: x + ".bai" if x else None)

# Reorder columns for output
output_df = sampled_df[["Sample", "Tumor", "Normal", "TumorBam", "TumorBai", "NormalBam", "NormalBai"]]

# Save the result to a TSV file
output_file_path = "/gpfs/data/molecpathlab/development/chip_strelka_runs/data/sample.pairs.tsv"
output_df.to_csv(output_file_path, sep="\t", index=False)
sample_df_path = "/gpfs/data/molecpathlab/development/chip_strelka_runs/data/sample_df.csv"
sampled_df.to_csv(sample_df_path, index=False)