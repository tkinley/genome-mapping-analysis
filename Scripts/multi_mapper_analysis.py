import pysam
import pandas as pd

def count_multi_mappers(sam_file):
    mapping_info = {}
    with pysam.AlignmentFile(sam_file, "r") as samfile:
        for read in samfile.fetch():
            if read.is_unmapped:
                continue
            gene_id = read.reference_name
            if gene_id not in mapping_info:
                mapping_info[gene_id] = {"multi_mapped": 0, "unique_mapped": 0}
            if read.mapping_quality == 0:
                mapping_info[gene_id]["multi_mapped"] += 1
            else:
                mapping_info[gene_id]["unique_mapped"] += 1
    return pd.DataFrame.from_dict(mapping_info, orient="index").reset_index().rename(columns={"index": "gene_id"})

# Count multi-mappers and unique-mappers in each SAM file
df_6A = count_multi_mappers("data/6A.sam")
df_6U = count_multi_mappers("data/6U.sam")
df_EA = count_multi_mappers("data/EA.sam")
df_EU = count_multi_mappers("data/EU.sam")

# Combine all DataFrames
merged_df = df_6A.set_index('gene_id').add(df_6U.set_index('gene_id'), fill_value=0) \
                 .add(df_EA.set_index('gene_id'), fill_value=0) \
                 .add(df_EU.set_index('gene_id'), fill_value=0) \
                 .reset_index()

# Summarize the total counts
merged_df['total_multi_mapped'] = merged_df['multi_mapped'] + merged_df['unique_mapped']
merged_df['total_unique_mapped'] = merged_df['unique_mapped']

# Print the first few rows of the DataFrame
print(merged_df.head())

# Save the merged DataFrame to a CSV file
merged_df.to_csv("results/merged_mapping_info.csv", index=False)
