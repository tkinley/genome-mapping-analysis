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
