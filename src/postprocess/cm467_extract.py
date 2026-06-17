import pandas as pd
import os

# Define paths
dataset_path = "/home/amity/Documents/experiments/data/reference/s2_dataset_extracted/Dataset1/Dataset1.txt"
output_path = "/home/amity/Documents/experiments/results/CM467_extracted_data.tsv"

# Create results directory if it doesn't exist
os.makedirs(os.path.dirname(output_path), exist_ok=True)

print(f"Reading dataset from {dataset_path}...")
# Read the tab-separated file
# We use low_memory=False to avoid type inference warnings on large files
df = pd.read_csv(dataset_path, sep='\t', low_memory=False)

# Check if Intensity CM467 exists
if 'Intensity CM467' not in df.columns:
    print("Error: 'Intensity CM467' column not found in Dataset1.txt")
    print("Available columns:", df.columns.tolist())
else:
    # Filter for rows where Intensity CM467 is not NaN
    cm467_df = df[df['Intensity CM467'].notna()].copy()
    
    # Select relevant columns
    # Sequence, Intensity CM467, Proteins, Gene names, Protein names
    relevant_cols = ['Sequence', 'Intensity CM467', 'Proteins', 'Gene names', 'Protein names', 'Length', 'Mass', 'Score']
    
    # Ensure all relevant columns exist
    existing_cols = [col for col in relevant_cols if col in cm467_df.columns]
    cm467_df = cm467_df[existing_cols]
    
    print(f"Extracted {len(cm467_df)} rows for CM467.")
    
    # Save to TSV
    cm467_df.to_csv(output_path, sep='\t', index=False)
    print(f"Saved extracted data to {output_path}")

    # Display some stats
    print("\nExtraction Summary:")
    print(f"Total peptides found for CM467: {len(cm467_df)}")
    if 'Length' in cm467_df.columns:
        print("Peptide length distribution:")
        print(cm467_df['Length'].value_counts().sort_index())
