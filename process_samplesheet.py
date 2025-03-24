import os
import glob
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-f','--filepath', type=str,
                   help='file path to all all the 65+ patients with Tumor and Normal Sample IDs')
parser.add_argument("-o", "--output", help = "Output Path")

args = parser.parse_args()
csv_file_path = args.filepath
output_dir = args.output

def find_matching_file(directory, id):
    # Search for files containing the ID
    pattern = os.path.join(directory, f"*{id}*.dd.ra.rc.bam")
    matching_files = glob.glob(pattern)
    
    # Check if any files were found
    if matching_files:
        return matching_files[0]  # Return the first match
    else:
        return None

# Load the CSV file
df = pd.read_csv(csv_file_path)
sampled_df = df
sampled_df["Sample"] = sampled_df["Tumor"] + "_" + sampled_df["Normal"]

base_path = "output/alignments/recalibrated"

# Find matching BAM and BAI files
sampled_df["TumorBam"] = sampled_df.apply(lambda row: find_matching_file(os.path.join(row['Run'],base_path), row['Tumor']), axis=1)
sampled_df["TumorBai"] = sampled_df["TumorBam"].apply(lambda x: x + ".bai" if x else None)

sampled_df["NormalBam"] = sampled_df.apply(lambda row: find_matching_file(os.path.join(row['Run'],base_path), row['Normal']), axis=1)
sampled_df["NormalBai"] = sampled_df["NormalBam"].apply(lambda x: x + ".bai" if x else None)

# Get samplesheet columns
output_df = sampled_df[["Sample", "Tumor", "Normal", "TumorBam", "TumorBai", "NormalBam", "NormalBai"]]

output_file_path = os.path.join(output_dir,"sample.pairs.tsv")
output_df.to_csv(output_file_path, sep="\t", index=False)
sample_df_path = os.path.join(output_dir,"sample_df_all.csv")
sampled_df.to_csv(sample_df_path, index=False)