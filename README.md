# Horse Metformin Project

This study investigates the effect of metformin on atrial fibrillation (AF) in a cohort of retired racehorses. We induced AF in 20 horses and monitored them over a four-month period. Half of these horses were treated with metformin, while the other half received a placebo. Additionally, a sham group of four healthy horses in normal sinus rhythm was included as controls.

## Experimental Groups:
- **Placebo Group**: Horses with AF treated with placebo for four months (n = 10)
- **Metformin Group**: Horses with AF treated with metformin for four months (n = 10)
- **Sham Group**: Sham-operated horses in normal sinus rhythm with no treatment (n = 4)

For proteomics, baseline biopsies from the right atrium were collected before inducing AF, allowing for a time-course comparison. All other analyses focus on samples obtained from the left and right atrium after four months of disease.

The project is divided into the following main analysis workflows:

### Proteomics
1. **Terminal Samples**: Analyzes terminal samples for differential protein abundance, enrichment, and pathway interactions.
2. **Timecourse Samples**: Evaluates dynamic changes in protein abundance across different time points in horses treated with metformin compared to controls.

### RNA-seq
Performs transcriptome-wide differential expression analysis and gene enrichment using RNA-seq data.

### Integration of Omics
Combines multi-omics data (proteomics and transcriptomics) for joint dimensionality reduction analysis (MCIA) and functional enrichment of high-loading features.

## Directory Structure
- Horse_Metformin_project.Rproj
- README.md
- Integration of Omics
  - analysis
    - 01_Venn_Diagrams
    - 02_Enrichment_MCIA_Loadings
  - data
- Proteomics
  - Terminal
    - analysis
    - data
  - Timecourse
    - analysis
    - data
- RNA-seq
  - analysis
  - data
  - preprocessing
## Directory Descriptions

### RNA-seq
Contains all RNA sequencing data and analysis scripts.

- `analysis/`: Scripts for differential gene expression analysis, gene enrichment, and other downstream analyses.
  - `01_dge/`: Scripts for differential gene expression analysis.
  - `02_gene-enrichment/`: Gene Ontology (GO) and KEGG pathway enrichment analysis.
  - `03_mitoXplorer/`: Mitochondrial-focused transcriptomics analysis.
  - `04_ChordDiagram/`: Chord diagrams for visualizing gene connections.

- `data/`: Includes count matrices with gene expression data, gene annotations, and metadata, as well as RNA-seq data formatted for analysis with mitoXplorer.

- `preprocessing/`: Scripts and outputs for initial data preprocessing, quality control, and cleaning steps before downstream analysis (including Snakefiles).

### Proteomics
The proteomics analyses are divided into terminal samples and time-course data.

- `analysis/`: Contains scripts for differential protein abundance analysis, enrichment, and pathway interaction analysis.
  - `01_dge/`: Scripts for differential abundance analysis.
  - `02_gene-enrichment/`: Gene Ontology (GO) and KEGG pathway enrichment analysis.
  - `03_string/`: Analyzes protein-protein interactions using STRING and exports networks to Cytoscape for visualization.
  - `04_matrisome/`: Characterizes extracellular matrix and matrisome-related proteins.

- `data/`: Includes count matrices of raw protein abundances, protein annotations (including STRING and matrisome annotations), and metadata.

### Integration of Omics
This folder contains scripts and data for integrating multiple omics platforms using MCIA and conducting enrichment analysis based on the identified loading features.

- `analysis/`: Contains scripts for generating Venn diagrams, enrichment analysis of MCIA loadings, and visualization outputs.

- `data/`: Includes loading results from MCIA performed in the OmicsAnalyst Platform.

## Running the Analysis
The analysis scripts are available in the relevant analysis subdirectories as Rmarkdown files (.Rmd). For visualization of the results, HTML files (.html) document the steps of data analysis.

The project uses `pacman` for package management to simplify dependency loading.
