# Genome Mapping Analysis

This repository contains scripts for analyzing genome mapping data from SAM files. The scripts extract mapping information, summarize mapped reads across genes, and identify multi-mapped reads.

## Structure

- `Scripts/`: Contains the Python scripts for data processing and analysis.
- `data/`: Directory to store input SAM files.
- `results/`: Directory to store output CSV files.

## Scripts

### 1. Multi-Mapper Analysis

The `multi_mapper_analysis.py` script identifies and counts multi-mapped and unique-mapped reads in SAM files.

### 2. Unmapped Analysis

The `unmapped_analysis_script.py` script merges data from multiple SAM files, summarizes mapped reads per gene, and filters genes using read-count thresholds.

## Usage

1. Install the required packages:
    ```bash
    pip install -r requirements.txt
    ```

2. Run the multi-mapper analysis script:
    ```bash
    python Scripts/multi_mapper_analysis.py
    ```

3. Run the unmapped analysis script:
    ```bash
    python Scripts/unmapped_analysis_script.py
    ```

## Requirements

- pandas
- pysam

## Example Data

Place your SAM files in the `data/` directory. 
