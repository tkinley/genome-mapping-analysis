import pysam
import pandas as pd

def count_multi_mappers(sam_file):
    samfile = pysam.AlignmentFile(sam_file, "r")
    mapping_info = {'gene_id': [], 'multi_mapped': [], 'unique_mapped': []}
    for read in samfile.fetch():
        if read.is_unmapped:
            continue
        gene_id = read.reference_name
        if read.mapping_quality == 0:
            mapping_info['gene_id'].append(gene_id)
            mapping_info['multi_mapped'].append(1)
            mapping_info['unique_mapped'].append(0)
        else:
            mapping_info['gene_id'].append(gene_id)
            mapping_info['multi_mapped'].append(0)
            mapping_info['unique_mapped'].append(1)
    return pd.DataFrame(mapping_info).drop_duplicates(subset=['gene_id'])


def main():
    """Run multi-mapper and unique-mapper analysis on SAM files."""
    # Count multi-mappers and unique-mappers in each SAM file
    df_6A = count_multi_mappers("data/6A.sam")
    df_6U = count_multi_mappers("data/6U.sam")
    df_EA = count_multi_mappers("data/EA.sam")
    df_EU = count_multi_mappers("data/EU.sam")

    # Combine all DataFrames
    merged_df = (
        df_6A.set_index('gene_id')
        .add(df_6U.set_index('gene_id'), fill_value=0)
        .add(df_EA.set_index('gene_id'), fill_value=0)
        .add(df_EU.set_index('gene_id'), fill_value=0)
        .reset_index()
    )

    # Summarize the total counts
    merged_df['total_multi_mapped'] = (
        merged_df['multi_mapped'] + merged_df['unique_mapped']
    )
    merged_df['total_unique_mapped'] = merged_df['unique_mapped']

    # Print the first few rows of the DataFrame
    print(merged_df.head())

    # Save the merged DataFrame to a CSV file
    merged_df.to_csv("results/merged_mapping_info.csv", index=False)


if __name__ == "__main__":
    main()
