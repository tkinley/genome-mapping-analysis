# Genome Mapping Analysis

This repository contains scripts for analyzing genome mapping data from SAM files. The scripts extract mapping information, filter genes based on mapped reads, and identify multi-mapped reads.

## Structure

- `scripts/`: Contains the Python scripts for data processing and analysis.
- `data/`: Directory to store input SAM files.
- `results/`: Directory to store output CSV files.

## Scripts

### 1. Gene Filtering

The `gene_filtering.py` script loads SAM files, merges mapping data, and filters genes based on mapped read thresholds.

### 2. Multi-Mapper Analysis

The `multi_mapper_analysis.py` script identifies and counts multi-mapped and unique-mapped reads in SAM files.

## Usage

1. Install the required packages:
    ```bash
    pip install -r requirements.txt
    ```

2. Run the gene filtering script:
    ```bash
    python scripts/gene_filtering.py
    ```

3. Run the multi-mapper analysis script:
    ```bash
    python scripts/multi_mapper_analysis.py
    ```

## Requirements

- pandas
- pysam

## Example Data

Place your SAM files in the `data/` directory. 
