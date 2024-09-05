# Preprocessing Pipeline

This repository contains the preprocessing pipeline for RNA-seq data, which includes the following steps:

1. **Quality Control**: Assessing the quality of raw reads using FastQC.
2. **Trimming**: Removing low-quality bases and adapter sequences using fastp.
3. **Alignment**: Mapping cleaned reads to the reference genome using STAR.
4. **Quantification**: Counting reads per gene using featureCounts.
5. **Quality Assessment**: Generating summary reports using MultiQC and Qualimap.

The final output of the preprocessing pipeline is the gene expression for downstream analysis, which can be found in Horse_AF_project\RNA-seq\data\count

## Structure
- `Snakefile/`: Contains the Snakefile used for processing of reads 
- `qc-raw-report/`: Contains MultiQC reports for raw data.
- `qc-clean-report/`: Contains MultiQC reports for cleaned data.
- `mapping-report/`: Contains MultiQC reports for mapping.
- `count-report/`: Contains MultiQC reports for quantification.
- `qc-qualimap/`: Contains Qualimap reports.