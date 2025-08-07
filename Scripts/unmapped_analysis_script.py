import pandas as pd
import pysam

def get_mapping_info(sam_file):
    """Extract mapping information from a SAM file."""
    try:
        samfile = pysam.AlignmentFile(sam_file, "r")
    except FileNotFoundError:
        print(f"Error: File {sam_file} not found.")
        return {}
    
    mapping_info = {}
    for read in samfile.fetch():
        gene_id = read.reference_name
        if gene_id not in mapping_info:
            mapping_info[gene_id] = 0
        if not read.is_unmapped:
            mapping_info[gene_id] += 1
    return mapping_info

def load_and_merge_sam_files(files):
    """Load and merge SAM file data into a single DataFrame."""
    dataframes = []
    for file in files:
        mapping_info = get_mapping_info(file)
        df = pd.DataFrame(list(mapping_info.items()), columns=['gene_id', f'mapped_reads_{file.split(".")[0]}'])
        dataframes.append(df)
    
    merged_df = dataframes[0]
    for df in dataframes[1:]:
        merged_df = merged_df.merge(df, on='gene_id', how='outer')
    
    merged_df.fillna(0, inplace=True)
    return merged_df

def filter_genes_by_threshold(df, threshold):
    """Filter genes by total mapped reads threshold."""
    return df[df['total_mapped_reads'] >= threshold]

# List of SAM files
sam_files = ["data/6A.sam", "data/6U.sam", "data/EA.sam", "data/EU.sam"]

# Load and merge SAM file data
merged_df = load_and_merge_sam_files(sam_files)

# Summarize total mapped reads
merged_df['total_mapped_reads'] = merged_df.iloc[:, 1:].sum(axis=1)

# Define regex patterns for gene ID extraction
v1_pattern = r'^(\d+)\|'
v2_pattern = r'^(Vocar\d+)'
gene_name_pattern = r'^(Vocar\.\d+s\d+\.1)'

# Extract V1, V2, and gene names into new columns
merged_df['V1_ID'] = merged_df['gene_id'].str.extract(v1_pattern)
merged_df['V2_ID'] = merged_df['gene_id'].str.extract(v2_pattern)
merged_df['Gene_Name'] = merged_df['gene_id'].str.extract(gene_name_pattern)

# Process version 1 genes
version1_df = merged_df.drop(columns=['V2_ID', 'Gene_Name']).dropna(subset=['V1_ID'])
v1_threshold = version1_df['total_mapped_reads'].quantile(0.25)
filtered_v1_genes = filter_genes_by_threshold(version1_df, v1_threshold)

# Process version 2 genes
version2_df = merged_df.drop(columns=['V1_ID', 'Gene_Name']).dropna(subset=['V2_ID'])
v2_threshold = version2_df['total_mapped_reads'].quantile(0.25)
filtered_v2_genes = filter_genes_by_threshold(version2_df, v2_threshold)

# Output filtered DataFrames
print(filtered_v1_genes)
print(filtered_v2_genes)

# Save filtered genes to CSV
filtered_v1_genes.to_csv('filtered_v1_genes.csv', index=False)
filtered_v2_genes.to_csv('filtered_v2_genes.csv', index=False)
