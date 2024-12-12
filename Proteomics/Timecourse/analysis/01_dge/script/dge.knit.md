---
title: "Differential Abundance Analysis TimeCourse"
author: "Simon Haugaard & Ali Altıntaş"
date: "2024-12-06"
output: 
  html_document:
    toc: true
    toc_float: true
    toc_collapsed: true
    toc_depth: 3
    number_sections: true
---




# Required R libraries

``` r
if (!require("pacman")) install.packages("pacman")
```

```
## Indlæser krævet pakke: pacman
```

``` r
pacman::p_load("edgeR", "readr", "readxl", "biomaRt", "magrittr", "tibble", "stringr", 
               "ggplot2", "data.table", "patchwork", "openxlsx", "dplyr", "missForest", 
               "RColorBrewer", "limma", "DEqMS", "preprocessCore", "DEP", 
               "SummarizedExperiment", "Metrics", "fdrtool", "aamisc", "sva")

# install aamisc package for MDS and Volcano plots
# Commented out for future use:
# pacman::p_load("qvalue", "rain", "limma", "devtools")
# url <- "https://cran.r-project.org/src/contrib/Archive/HarmonicRegression/HarmonicRegression_1.0.tar.gz"
# pkgFile <- "HarmonicRegression_1.0.tar.gz"
# download.file(url = url, destfile = pkgFile)
# install.packages(pkgs=pkgFile, type="source", repos=NULL)
# file.remove(pkgFile)
# pacman::p_load_gh("altintasali/aamisc")

# Colours for publication
publication_colors <- c("Control" = "#285291", "Metformin" = "#7C1516", "Sham" = "#9B9B9B")
```

# Read data and filtrer

``` r
# Load count data from a local file and metadata through Excel
# Define file paths
count_file <- "../../../data/count/report.unique_genes_matrix.tsv"
meta_file <- "../../../data/metadata/meta_proteomics.xlsx"
geneinfo_file <- "../../../data/gene_annotation/horse_gene_annotation.tsv.gz"

# Read count data
count <- readr::read_delim(count_file)
```

```
## Rows: 2967 Columns: 74
## ── Column specification ────────────────────────────────────────────────────────
## Delimiter: "\t"
## chr  (1): Genes
## dbl (73): D:\Mass_spectrometry\Raw_data\Joakim\Simon Horse second run2\1A_RA...
## 
## ℹ Use `spec()` to retrieve the full column specification for this data.
## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
```

``` r
# Read metadata
meta <- readxl::read_excel(meta_file)

# Remove unnecessary or blank columns from MS machine
count <- count %>% select(-matches("D:\\\\Mass_spectrometry\\\\Raw_data\\\\Joakim\\\\Simon Horse second run2\\\\10A_RA11_1_24830.d"))

# Read and clean gene annotation file
geneinfo <- fread(geneinfo_file)
setnames(geneinfo, old = names(geneinfo), new = c("ENSEMBL", "ENSEMBLv", "Description_detailed", 
                                                  "Chr", "Start", "End", "Strand", "GENENAME", "ENTREZID", "Description"))
geneinfo <- geneinfo %>% select(-ENSEMBLv, -Description_detailed) %>% distinct(ENSEMBL, .keep_all = TRUE)

# Merge gene annotation with count data
annot <- merge(count[, "Genes", drop = FALSE], geneinfo, by.x = "Genes", by.y = "GENENAME", all.x = TRUE)

# Clean count data
count <- count %>% remove_rownames() %>% column_to_rownames(var="Genes")
cleaned_names <- gsub("D:\\\\Mass_spectrometry\\\\Raw_data\\\\Joakim\\\\Simon Horse second run2\\\\|_.*", "", colnames(count))
colnames(count) <- cleaned_names
meta$`Sample ID` <- cleaned_names

#subset data to only RA terminal and baseline samples
meta <- meta[meta$Region %in% c("RA"), ]
count <- count[, colnames(count) %in% meta$`Sample ID`]


# Remove low count samples from Count and Meta
samples_to_remove <- c("25C", "12A", "22A", "19A")
count <- count[, !colnames(count) %in% samples_to_remove]
meta <- meta %>% filter(!`Sample ID` %in% samples_to_remove)

# Calculate the number of proteins in each sample
num_proteins <- colSums(!is.na(count))
num_proteins <- num_proteins[meta$`Sample ID`]  # Ensure matching order with meta

# Define color scheme for groups
KUalt <- c("#7C1516", "#285291", "#434343", "#999999")  # Custom color palette
group_colors <- setNames(KUalt[1:length(unique(meta$Group))], unique(meta$Group))
bar_colors <- group_colors[meta$Group]

# Plot the number of proteins in each sample before filtering
par(mfrow = c(1, 1))
barplot(num_proteins, 
        main = "Number of Proteins in Each Sample\nBefore Filtering",
        xlab = "Sample",
        ylab = "Number of Proteins",
        col = bar_colors,
        border = "black",
        ylim = c(0, 2500),
        las = 2)

# Print the total number of proteins before filtering
cat("The total number of proteins before filtering was", max(num_proteins), "\n")
```

```
## The total number of proteins before filtering was 2610
```

``` r
# Identify genes with less than three valid values in any condition
invalid_genes <- unique(unlist(lapply(unique(meta$Condition), function(condition) {
  condition_columns <- meta$`Sample ID`[meta$Condition == condition]
  condition_count <- count[, condition_columns, drop = FALSE]
  rownames(condition_count)[rowSums(!is.na(condition_count)) < 3]
})))

# Remove invalid genes from the original count table
filtered_count <- count[!rownames(count) %in% invalid_genes, ]

# Calculate the number of proteins in each sample after filtering
num_proteins_after <- colSums(!is.na(filtered_count))
num_proteins_after <- num_proteins_after[meta$`Sample ID`]  # Ensure matching order with meta

# Plot settings
par(mfrow = c(1, 2))

# First plot: Number of proteins in each sample before filtering
barplot(num_proteins, 
        main = "Number of Proteins in Each Sample\nBefore Filtering",
        xlab = "Sample",
        ylab = "Number of Proteins",
        col = bar_colors,
        border = "black",
        ylim = c(0, 2500),
        las = 2)

# Second plot: Number of proteins in each sample after filtering
barplot(num_proteins_after, 
        main = "Number of Proteins in Each Sample\nAfter Filtering",
        xlab = "Sample",
        ylab = "Number of Proteins",
        col = bar_colors,
        border = "black",
        ylim = c(0, 2500),
        las = 2)
```

<img src="dge_files/figure-html/read_data-1.png" width="672" />

``` r
# Print the total and average number of proteins before and after filtering
cat("The total number of proteins before filtering was", max(num_proteins), "\n")
```

```
## The total number of proteins before filtering was 2610
```

``` r
cat("The average number of proteins before filtering was", mean(num_proteins), "\n")
```

```
## The average number of proteins before filtering was 2168.578
```

``` r
cat("The total number of proteins after filtering was", max(num_proteins_after), "\n")
```

```
## The total number of proteins after filtering was 1288
```

``` r
cat("The average number of proteins after filtering was", mean(num_proteins_after), "\n")
```

```
## The average number of proteins after filtering was 1258.933
```
# DEP-package  https://bioconductor.org/packages/release/bioc/vignettes/DEP/inst/doc/DEP.html
## Generate a SummarizedExperiment object

``` r
# Generate proxy columns similar to DEP's sample data to make the structure compatible.
proxy_column_names <- c("Protein.IDs", "Majority.protein.IDs", "Protein.names", "Gene.names", 
                        "Fasta.headers", "Peptides", "Razor...unique.peptides", "Unique.peptides", 
                        "Only.identified.by.site", "Reverse", "Potential.contaminant")
proxy_df <- data.frame(matrix(NA, nrow=nrow(filtered_count), ncol=length(proxy_column_names)))
colnames(proxy_df) <- proxy_column_names

# Assign row names to specific columns
rownames_to_columns <- rownames(filtered_count)
proxy_df[1:4] <- lapply(1:4, function(i) rownames_to_columns)

# Merge and rename sample columns
merged_df <- cbind(proxy_df, filtered_count)
colnames(merged_df)[12:ncol(merged_df)] <- paste0("LFQ.intensity.", colnames(filtered_count))

# Reorder columns
lfq_columns <- grep("^LFQ.intensity", colnames(merged_df), value = TRUE)
new_order <- c(proxy_column_names[1:8], lfq_columns, proxy_column_names[9:11])
data <- merged_df[, new_order]

print(colnames(data))
```

```
##  [1] "Protein.IDs"             "Majority.protein.IDs"   
##  [3] "Protein.names"           "Gene.names"             
##  [5] "Fasta.headers"           "Peptides"               
##  [7] "Razor...unique.peptides" "Unique.peptides"        
##  [9] "LFQ.intensity.1A"        "LFQ.intensity.1B"       
## [11] "LFQ.intensity.2A"        "LFQ.intensity.2B"       
## [13] "LFQ.intensity.3A"        "LFQ.intensity.3B"       
## [15] "LFQ.intensity.4A"        "LFQ.intensity.4B"       
## [17] "LFQ.intensity.5A"        "LFQ.intensity.5B"       
## [19] "LFQ.intensity.6A"        "LFQ.intensity.6B"       
## [21] "LFQ.intensity.7A"        "LFQ.intensity.7B"       
## [23] "LFQ.intensity.8A"        "LFQ.intensity.8B"       
## [25] "LFQ.intensity.9A"        "LFQ.intensity.9B"       
## [27] "LFQ.intensity.10A"       "LFQ.intensity.10B"      
## [29] "LFQ.intensity.11A"       "LFQ.intensity.11B"      
## [31] "LFQ.intensity.12B"       "LFQ.intensity.13A"      
## [33] "LFQ.intensity.13B"       "LFQ.intensity.14A"      
## [35] "LFQ.intensity.14B"       "LFQ.intensity.15A"      
## [37] "LFQ.intensity.15B"       "LFQ.intensity.16A"      
## [39] "LFQ.intensity.16B"       "LFQ.intensity.17A"      
## [41] "LFQ.intensity.17B"       "LFQ.intensity.18A"      
## [43] "LFQ.intensity.18B"       "LFQ.intensity.19B"      
## [45] "LFQ.intensity.20A"       "LFQ.intensity.20B"      
## [47] "LFQ.intensity.22B"       "LFQ.intensity.23A"      
## [49] "LFQ.intensity.23B"       "LFQ.intensity.24A"      
## [51] "LFQ.intensity.24B"       "LFQ.intensity.25A"      
## [53] "LFQ.intensity.25B"       "Only.identified.by.site"
## [55] "Reverse"                 "Potential.contaminant"
```

``` r
#Are there any duplicated gene names?
data$Gene.names %>% duplicated() %>% any()
```

```
## [1] FALSE
```

``` r
# Make unique names using the annotation in the "Gene.names" column as primary names and the annotation in "Protein.IDs" as name for those that do not have an gene name.
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers

# Prepare the metadata for DEP analysis
experimental_design <- data.frame(
  label = colnames(filtered_count),
  condition = meta$Condition,
  replicate = NA  # Initialize replicate column with NA
)

# Assign replicate numbers within each condition
experimental_design <- experimental_design %>%
  group_by(condition) %>%
  mutate(replicate = row_number()) %>%
  ungroup()

# Set row names
row.names(experimental_design) <- experimental_design$label
```

```
## Warning: Setting row names on a tibble is deprecated.
```

``` r
# Ensure the row names of experimental_design are unique
experimental_design <- as.data.frame(experimental_design)
rownames(experimental_design) <- make.unique(rownames(experimental_design))

# Print the experimental design to verify
print(experimental_design)
```

```
##    label             condition replicate
## 1     1A RA_Metformin_Baseline         1
## 2     1B  RA_Metformin_4months         1
## 3     2A   RA_Control_Baseline         1
## 4     2B    RA_Control_4months         1
## 5     3A RA_Metformin_Baseline         2
## 6     3B  RA_Metformin_4months         2
## 7     4A      RA_Sham_Baseline         1
## 8     4B       RA_Sham_4months         1
## 9     5A RA_Metformin_Baseline         3
## 10    5B  RA_Metformin_4months         3
## 11    6A   RA_Control_Baseline         2
## 12    6B    RA_Control_4months         2
## 13    7A RA_Metformin_Baseline         4
## 14    7B  RA_Metformin_4months         4
## 15    8A   RA_Control_Baseline         3
## 16    8B    RA_Control_4months         3
## 17    9A RA_Metformin_Baseline         5
## 18    9B  RA_Metformin_4months         5
## 19   10A   RA_Control_Baseline         4
## 20   10B    RA_Control_4months         4
## 21   11A RA_Metformin_Baseline         6
## 22   11B  RA_Metformin_4months         6
## 23   12B    RA_Control_4months         5
## 24   13A RA_Metformin_Baseline         7
## 25   13B  RA_Metformin_4months         7
## 26   14A   RA_Control_Baseline         5
## 27   14B    RA_Control_4months         6
## 28   15A   RA_Control_Baseline         6
## 29   15B    RA_Control_4months         7
## 30   16A   RA_Control_Baseline         7
## 31   16B    RA_Control_4months         8
## 32   17A RA_Metformin_Baseline         8
## 33   17B  RA_Metformin_4months         8
## 34   18A      RA_Sham_Baseline         2
## 35   18B       RA_Sham_4months         2
## 36   19B  RA_Metformin_4months         9
## 37   20A   RA_Control_Baseline         8
## 38   20B    RA_Control_4months         9
## 39   22B    RA_Control_4months        10
## 40   23A RA_Metformin_Baseline         9
## 41   23B  RA_Metformin_4months        10
## 42   24A      RA_Sham_Baseline         3
## 43   24B       RA_Sham_4months         3
## 44   25A      RA_Sham_Baseline         4
## 45   25B       RA_Sham_4months         4
```

``` r
# Create the SummarizedExperiment Object
# This will allow DEP to use the data for further analysis and visualization
data_se <- make_se(data_unique, LFQ_columns, experimental_design)
```
## Protein Coverage & filtrer on missing values

``` r
# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)
```

<img src="dge_files/figure-html/Protein Coverage-1.png" width="672" />

``` r
plot_numbers(data_se)
```

<img src="dge_files/figure-html/Protein Coverage-2.png" width="672" />

``` r
plot_coverage(data_se)
```

<img src="dge_files/figure-html/Protein Coverage-3.png" width="672" />

``` r
# Filter for proteins that are identified in all replicates of at least one condition 
# Commented out as we have already performed filtrering manually in the first step, which unlike the terminal samples, will be used here. 
#data_filt <- filter_missval(data_se, thr = 0)
#plot_frequency(data_filt)
```
## Normalizaiton

``` r
#Check un-normalized data
plot_normalization(data_se)
```

<img src="dge_files/figure-html/Normalization-1.png" width="672" />

``` r
#VSN-normalization
data_se_norm <- normalize_vsn(data_se)

#Plot VSN-norm data
meanSdPlot(data_se_norm, rank = TRUE)
```

<img src="dge_files/figure-html/Normalization-2.png" width="672" />

``` r
plot_normalization(data_se_norm)
```

<img src="dge_files/figure-html/Normalization-3.png" width="672" />

``` r
###Verdict: After variance stabilisation, the median (a reasonable estimator of the standard deviation of feature level data conditional on the mean) is approximately a horizontal line.
```
### Impute data for missing values

``` r
# Plot a heatmap of proteins with missing values
plot_missval(data_se_norm)
```

<img src="dge_files/figure-html/Imputation for Missing Values-1.png" width="672" />

``` r
# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_se_norm)
```

<img src="dge_files/figure-html/Imputation for Missing Values-2.png" width="672" />

``` r
# Typically, missing values in proteomics data are "Missing Not At Random" (MNAR),
# meaning proteins with missing values often have lower intensities and are close to the detection limit.
# For MNAR data, appropriate imputation methods include left-censored imputation techniques such as:
# - Quantile Regression-based Left-Censored (QRILC)
# - Random draws from a left-shifted distribution ("MinProb" or "man")

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp_MinProb <- impute(data_se_norm, fun = "MinProb", q = 0.01) 
```

```
## Imputing along margin 2 (samples/columns).
```

```
## [1] 0.4919121
```

``` r
# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute( data_se_norm, fun = "man", shift = 1.8, scale = 0.3) 


# Impute missing data using Quantile Regression Imputation of Left-Censored Data 
data_imp_qrilc <- impute( data_se_norm, fun = "QRILC") 
```

```
## Imputing along margin 2 (samples/columns).
```

``` r
# Plot intensity distributions before and after imputation
# This plot helps assess the effectiveness of each imputation method by visualizing the changes in intensity distributions.
plot_imputation( data_se_norm, data_imp_MinProb, data_imp_man, data_imp_qrilc)
```

<img src="dge_files/figure-html/Imputation for Missing Values-3.png" width="672" />

``` r
### Verdict: After visual inspection, it appears that the MinProb handles the MNAR nature of the data best. and will be used downstream. 
```

## Differential abundance analysis using DEP

``` r
# Note: This analysis is intended as an extra quality control (QC) step.
# The main differential abundance analysis is performed using the `limma` package due to the more complex study design.
#Set contrasts
con <- c("RA_Control_4months_vs_RA_Sham_4months", "RA_Metformin_4months_vs_RA_Control_4months")

# Perform differential analysis using the specified contrasts - best on non-imputed data, honestly
data_diff <- test_diff(data_imp_qrilc, type = "manual", test = con)
```

```
## Tested contrasts: RA_Control_4months_vs_RA_Sham_4months, RA_Metformin_4months_vs_RA_Control_4months
```

``` r
# Denote significant proteins based on user defined cutoffs
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(0))

# Generate a long data.frame
df_long <- get_df_long(dep)
df_wide <- get_df_wide(dep)
```

# MDS-plot

``` r
# Transform raw count data and remove NA values
logCounts <- log2(count + 1)  # Add 1 to avoid log(0) issues
logCounts_noNA <- na.omit(logCounts)

# Create DGEList object and define the design matrix
annot_reordered <- annot[match(rownames(logCounts_noNA), annot$Genes), ]
d <- DGEList(counts = logCounts_noNA, genes = annot_reordered, samples = meta)
design <- model.matrix(~0 + Condition, data = d$samples)

# 1. MDS Plot for Raw Data
mds_raw <- plotMDS(logCounts_noNA, plot = FALSE)
mds_plot_raw <- aamisc::ggMDS(mds = mds_raw, meta = d$samples, color.by = "Group", shape.by = "Timepoint", legend.position = "right", text.by = "Horse", text.size = 1.5) +
  scale_color_manual(values = publication_colors) +
  labs(title = "MDS: Raw Data", x = "MDS Dimension 1", y = "MDS Dimension 2") +
  theme(legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5))

# 2. MDS Plot for Data Adjusted by Blocking for Horse
logCounts_blocked <- removeBatchEffect(logCounts_noNA, batch = as.factor(d$samples$Horse), design = design)
```

```
## Coefficients not estimable: batch22 batch23
```

```
## Warning: Partial NA coefficients for 769 probe(s)
```

``` r
mds_blocked <- plotMDS(logCounts_blocked, plot = FALSE)
mds_plot_blocked <- aamisc::ggMDS(mds = mds_blocked, meta = d$samples, color.by = "Group", shape.by = "Timepoint", legend.position = "right", text.by = "Horse", text.size = 1.5) +
  scale_color_manual(values = publication_colors) +
  labs(title = "MDS: Blocked for Horse", x = "MDS Dimension 1", y = "MDS Dimension 2") +
  theme(legend.title = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0.5))


# Display the plots
print(mds_plot_raw)
```

<img src="dge_files/figure-html/unnamed-chunk-1-1.png" width="672" />

``` r
print(mds_plot_blocked)
```

<img src="dge_files/figure-html/unnamed-chunk-1-2.png" width="672" />


# Checking histograms and save/load imputations

``` r
# Histograms will be generated for four datasets:
# 1. `normalized_unimputed`: Normalized but not imputed data.
# 2. `imputed_data_MinProb`: Data imputed using the "MinProb" method.
# 3. `imputed_data_qrilc`: Data imputed using the "QRILC" method.
# 4. `unfiltered_count`: Raw log-transformed counts after initial filtering.

# Extract the assay data from data_se_norm (normalized_unimputed)
normalized_unimputed <- as.data.frame(assay(data_se_norm))
colnames(normalized_unimputed) <- colnames(filtered_count)
rownames(normalized_unimputed) <- rownames(data_se_norm)

# Extract the assay data from data_imp_MinProb (imputed)
imputed_data_MinProb <- as.data.frame(assay(data_imp_MinProb))
colnames(imputed_data_MinProb) <- colnames(filtered_count)
rownames(imputed_data_MinProb) <- rownames(data_se_norm)

# Extract the assay data from data_imp_qrilc
imputed_data_qrilc <- as.data.frame(assay(data_imp_qrilc))
colnames(imputed_data_qrilc) <- colnames(filtered_count)
rownames(imputed_data_qrilc) <- rownames(data_se_norm)

# Extract the assay data from data_imp_man
imputed_data_man <- as.data.frame(assay(data_imp_man))
colnames(imputed_data_man) <- colnames(filtered_count)
rownames(imputed_data_man) <- rownames(data_se_norm)

# Plot histogram for normalized data
par(mfrow = c(1, 4))
hist(as.vector(as.matrix(normalized_unimputed)), breaks = 50, main = "Normalized Data", xlab = "Intensity")
hist(as.vector(as.matrix(imputed_data_MinProb)), breaks = 50, main = "MinProb Imputed Data", xlab = "Intensity")
hist(as.vector(as.matrix(imputed_data_qrilc)), breaks = 50, main = "QRILC Imputed Data", xlab = "Intensity")
hist(as.vector(as.matrix(filtered_count)), breaks = 50, main = "Unfiltered Data", xlab = "Intensity")
```

<img src="dge_files/figure-html/Intensity-Histograms-1.png" width="960" />

``` r
# Define the output folder path
output_folder <- "../../../../Timecourse/analysis/01_dge/output/"

# Function to save a matrix with row names as a separate column
save_matrix <- function(matrix, file_path) {
  # Create a data frame from the matrix, including row names as a column
  matrix_df <- as.data.frame(matrix)
  matrix_df$GeneName <- rownames(matrix)
  fwrite(matrix_df, file = file_path, sep = "\t", row.names = FALSE)
}

# Function to load a matrix and restore row names
load_matrix <- function(file_path) {
  matrix_df <- fread(file_path, data.table = FALSE)
  rownames(matrix_df) <- matrix_df$GeneName
  matrix_df$GeneName <- NULL
  return(as.matrix(matrix_df))  # Convert back to a matrix
}

# Save the matrices
 # save_matrix(normalized_unimputed, paste0(output_folder, "normalized_unimputed.tsv"))
 # save_matrix(imputed_data_MinProb, paste0(output_folder, "imputed_data_MinProb.tsv"))
 # save_matrix(imputed_data_qrilc, paste0(output_folder, "imputed_data_qrilc.tsv"))

# Load the matrices in the future
normalized_unimputed <- load_matrix(paste0(output_folder, "normalized_unimputed.tsv"))
imputed_data_MinProb <- load_matrix(paste0(output_folder, "imputed_data_MinProb.tsv"))
imputed_data_qrilc <- load_matrix(paste0(output_folder, "imputed_data_qrilc.tsv"))
```

# Differential Abundance Analysis
## Limma without SVA

``` r
# Define the design matrix
design <- model.matrix(~ 0 + Condition, d$samples)
colnames(design) <- gsub("Condition", "", colnames(design))  # Clean column names

# Define contrasts for comparisons of interest
con <- makeContrasts(
  Metformin_vs_AF_RA = RA_Metformin_4months - RA_Control_4months,
  AF_vs_Sham_RA = RA_Control_4months - RA_Sham_4months,
  Terminal_vs_Baseline_Control = RA_Control_4months - RA_Control_Baseline,
  Terminal_vs_Baseline_Metformin = RA_Metformin_4months - RA_Metformin_Baseline,
  Diff_Treatment = (RA_Metformin_4months - RA_Metformin_Baseline) - (RA_Control_4months - RA_Control_Baseline),
  Diff_Disease = (RA_Control_4months - RA_Control_Baseline) - (RA_Sham_4months - RA_Sham_Baseline),
  Baseline_Difference_Metf_vs_AF = RA_Metformin_Baseline - RA_Control_Baseline,
  Baseline_Difference_AF_vs_Sham = RA_Control_Baseline - RA_Sham_Baseline,
  levels = design)
con
```

```
##                        Contrasts
## Levels                  Metformin_vs_AF_RA AF_vs_Sham_RA
##   RA_Control_4months                    -1             1
##   RA_Control_Baseline                    0             0
##   RA_Metformin_4months                   1             0
##   RA_Metformin_Baseline                  0             0
##   RA_Sham_4months                        0            -1
##   RA_Sham_Baseline                       0             0
##                        Contrasts
## Levels                  Terminal_vs_Baseline_Control
##   RA_Control_4months                               1
##   RA_Control_Baseline                             -1
##   RA_Metformin_4months                             0
##   RA_Metformin_Baseline                            0
##   RA_Sham_4months                                  0
##   RA_Sham_Baseline                                 0
##                        Contrasts
## Levels                  Terminal_vs_Baseline_Metformin Diff_Treatment
##   RA_Control_4months                                 0             -1
##   RA_Control_Baseline                                0              1
##   RA_Metformin_4months                               1              1
##   RA_Metformin_Baseline                             -1             -1
##   RA_Sham_4months                                    0              0
##   RA_Sham_Baseline                                   0              0
##                        Contrasts
## Levels                  Diff_Disease Baseline_Difference_Metf_vs_AF
##   RA_Control_4months               1                              0
##   RA_Control_Baseline             -1                             -1
##   RA_Metformin_4months             0                              0
##   RA_Metformin_Baseline            0                              1
##   RA_Sham_4months                 -1                              0
##   RA_Sham_Baseline                 1                              0
##                        Contrasts
## Levels                  Baseline_Difference_AF_vs_Sham
##   RA_Control_4months                                 0
##   RA_Control_Baseline                                1
##   RA_Metformin_4months                               0
##   RA_Metformin_Baseline                              0
##   RA_Sham_4months                                    0
##   RA_Sham_Baseline                                  -1
```

``` r
# Estimate correlation due to repeated measures (blocking for Horse)
corfit <- duplicateCorrelation(imputed_data_MinProb, design, block = as.factor(d$samples$Horse))
cat("Estimated consensus correlation:", corfit$consensus.correlation, "\n")
```

```
## Estimated consensus correlation: -7.09684e-05
```

``` r
# Fit the linear model
fit <- lmFit(imputed_data_MinProb, design, block = as.factor(d$samples$Horse), correlation = corfit$consensus.correlation)
rownames(fit$coefficients) <- rownames(imputed_data_MinProb)  # Retain row names for mapping

# Apply contrasts and empirical Bayes moderation
fit.contrast <- contrasts.fit(fit, contrasts = con)
fit.contrast <- eBayes(fit.contrast, robust = TRUE, trend = TRUE)

# Extract Differential Gene Expression (DGE) results
res <- list()  # List to store DGE results for each contrast
for (contrast in colnames(con)) {
  res_tmp <- topTable(fit.contrast, coef = contrast, adjust.method = "BH", number = Inf)
  res_tmp <- res_tmp[!is.na(res_tmp$t), ]  # Remove rows with NA values
  res_tmp$Contrast <- contrast  # Add contrast identifier
  res[[contrast]] <- res_tmp  # Store results
}

# Combine results across all contrasts into a single data frame
res_all <- do.call(rbind, res)

# Map Gene Names for easier interpretation (optional)
res_all$GeneName <- sapply(seq_len(nrow(res_all)), function(i) {
  gsub(paste0("^", res_all$Contrast[i], "\\."), "", rownames(res_all)[i])
})

# Split results for easier access/output
res_split <- split(res_all, res_all$Contrast)

# Visualize P-value Distributions
for (contrast in names(res_split)) {
  p_values <- res_split[[contrast]]$P.Value
  ggplot(data = data.frame(P.Value = p_values), aes(x = P.Value)) +
    geom_histogram(bins = 50, color = "black", fill = "lightblue") +
    theme_minimal() +
    labs(
      title = paste("P-value Histogram for", contrast),
      x = "P-value",
      y = "Frequency"
    ) +
    xlim(0, 1)  # Set x-axis limits between 0 and 1
  print(last_plot())
}
```

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-1.png" width="960" />

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-2.png" width="960" />

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-3.png" width="960" />

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-4.png" width="960" />

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-5.png" width="960" />

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-6.png" width="960" />

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-7.png" width="960" />

```
## Warning: Removed 2 rows containing missing values or values outside the scale range
## (`geom_bar()`).
```

<img src="dge_files/figure-html/limma-8.png" width="960" />

``` r
# VERDICT: Rising hill in P-value histogram detected (e.g. metformin_vs_AF_RA)
# Comment: A rising hill in the P-value histogram may indicate a systematic batch effect
# or a missing covariate that was not included in the initial design matrix.

cat("Detected potential systematic effects (rising hill in P-value histogram). Proceeding with SVA for adjustment...\n")
```

```
## Detected potential systematic effects (rising hill in P-value histogram). Proceeding with SVA for adjustment...
```



## Limma with SVA

``` r
# In this workflow, SVA is used to remove hidden variance, and a blocking factor 
# accounts for variance introduced by repeated measures (e.g., two samples from the same horse).

# Extract biological condition information
condition <- d$samples$Condition  # Specify the biological condition

# Apply SVA to Adjust for Hidden Confounders
## Full model including condition
mod <- model.matrix(~ condition, data = d$samples) 

## Null model excluding biological information
mod0 <- model.matrix(~ 1, data = d$samples)

# Define the design matrix for downstream analysis
# Exclude intercept for simpler contrast definitions in lmFit
design <- model.matrix(~ 0 + condition, data = d$samples)

# Estimate the number of surrogate variables
n.sv <- num.sv(imputed_data_MinProb, mod, method = "be")  
cat("Number of surrogate variables estimated:", n.sv, "\n")
```

```
## Number of surrogate variables estimated: 8
```

``` r
# Run SVA to estimate surrogate variables
svobj <- sva(imputed_data_MinProb, mod, mod0, n.sv = n.sv)
```

```
## Number of significant surrogate variables is:  8 
## Iteration (out of 5 ):1  2  3  4  5
```

``` r
# Test on WSVA
wsva(y = imputed_data_MinProb, design = design, plot = TRUE, n.sv = n.sv) # Also Finds 8 surrogate variables
```

<img src="dge_files/figure-html/unnamed-chunk-2-1.png" width="672" />

```
##           SV1         SV2        SV3       SV4         SV5         SV6
## 1A  1.0170971 -1.25398260 -0.9203357 1.0796055 -0.28919025  1.42365548
## 1B  1.0477800 -1.15699226 -0.8269477 1.0292767 -0.18927732  0.93656216
## 2A  1.0138614 -0.74111269 -1.0663759 0.9162007 -1.01699050  0.86344124
## 2B  0.9938475 -0.55409943 -1.0264045 0.9235017 -0.17017534  0.91867023
## 3A  1.0220489 -0.79966915 -1.0765506 0.9760710 -1.27608863  1.09940105
## 3B  1.0486587 -1.13690070 -0.6826656 1.0145288 -0.61974330  1.21337964
## 4A  0.7765506 -1.22937714 -1.9986527 1.0869874 -0.61151088  2.01232509
## 4B  0.9935679  0.08634185 -0.8409351 0.9950601  1.23510682  1.13202941
## 5A  0.9957720 -1.11736613 -1.4273800 0.8777079 -0.46094988  0.71250112
## 5B  1.0240529 -1.25239242 -0.7689843 0.9897433 -1.51414156  1.68892633
## 6A  0.9801450 -1.44371663 -1.2319859 0.9931967  1.13280493  1.42494601
## 6B  1.0052125  0.01367230 -0.7513574 0.9655632  0.94848867  1.01249110
## 7A  1.0121371 -1.02776188 -1.2018788 0.8372954 -0.39227267  0.83664904
## 7B  1.0227216 -0.76315104 -0.6757432 1.1485396 -1.26893145  1.04262740
## 8A  1.0086994 -1.10643406 -1.3044271 0.8436565  0.36158918  0.49067052
## 8B  1.0481232 -1.06492210 -0.6596542 0.9719714 -1.38142960  1.29359950
## 9A  0.9985343 -1.00863983 -0.9861130 1.1407475 -0.86201836  1.08799497
## 9B  1.0273786 -0.38602186 -0.8672133 0.9311572 -0.05531169  1.05366938
## 10A 0.9377540 -0.85171884 -1.4214371 0.9683597 -0.50113491  0.70550115
## 10B 1.0397173 -0.80535205 -0.5899910 1.1394358 -0.45218734  0.02964132
## 11A 0.9381384 -1.29597918 -1.2976243 1.1803860 -0.04067446  0.84757452
## 11B 1.0190957  0.15784179 -0.6915113 1.0094660  0.79128385  1.41576385
## 12B 1.0358260 -1.30360297 -0.8526275 0.9434360 -0.13422102  1.12202863
## 13A 0.6465615 -1.22073173 -0.3816647 0.7667680 -0.59365706  0.44273495
## 13B 1.0277462 -1.03901278 -0.9698071 0.7755931 -0.27681644  0.44058829
## 14A 0.9792192 -0.41420217 -1.0444255 1.0878228 -2.00808294  1.16509832
## 14B 1.0363395 -0.96986792 -0.6367784 1.0519027 -1.15709445  1.18129456
## 15A 0.9893579 -1.06473446 -1.3320707 0.9050185 -0.20695364 -0.43085291
## 15B 1.0503491 -1.03425833 -0.9304579 0.7524322 -0.87677150  0.72535834
## 16A 1.0467508 -1.02350603 -1.0690229 0.8085938 -0.33435635  1.00222155
## 16B 0.9963057 -1.01015719 -0.6807928 1.0258378 -1.36718873  1.11541806
## 17A 1.0241838 -1.23083609 -0.9690620 1.0818548  0.98048799  0.28556982
## 17B 1.0463102 -1.07849707 -0.9510019 0.7700829 -0.28378214  0.44473441
## 18A 0.9786495 -1.37755031 -1.1858105 1.0547257  1.13685062  0.77031853
## 18B 1.0259749 -0.76372410 -0.7575374 0.8626013 -1.45203393  1.22311680
## 19B 1.0505361 -1.04844919 -0.6286679 1.0301848 -0.81588322  0.92572698
## 20A 0.8644761 -1.41980396 -1.4937535 1.3555504  0.71738575 -0.02176214
## 20B 1.0489730 -0.55313120 -0.5728718 0.9319941 -1.75363502  1.33958966
## 22B 1.0149074 -0.61949081 -0.6556935 1.1300950 -0.53087078  0.64588798
## 23A 0.9655400 -1.43785292 -1.0810095 1.0205243  0.47367901  1.05865300
## 23B 1.0349637 -0.48238891 -0.8247340 0.9344748 -0.54559751  0.51751728
## 24A 0.9725134 -0.92612665 -1.4037629 1.1014202 -2.44094656 -0.35101176
## 24B 1.0378459 -0.73545244 -0.5533829 1.1485575 -1.28601980  0.90640602
## 25A 0.9987605 -1.51315841 -0.7877455 1.0186798  1.89279125  1.22964497
## 25B 1.0384414 -0.18332997 -0.6614969 1.0756166 -0.52628328  0.46557312
##           SV7        SV8
## 1A  1.0417751  0.9537164
## 1B  0.9292773  1.1499795
## 2A  1.0191331  0.8983325
## 2B  1.1403034  0.9568046
## 3A  1.0189757  0.8043103
## 3B  0.9964634  1.0224089
## 4A  0.7936713  1.3710771
## 4B  1.0079163  0.9323738
## 5A  1.0262370  0.7834478
## 5B  0.9985203  0.6868392
## 6A  0.7431116  0.6553652
## 6B  0.9198022  0.7823037
## 7A  1.0070442  1.0874378
## 7B  0.8902253  1.3024706
## 8A  1.0590529  1.3073453
## 8B  0.8815998  0.7764818
## 9A  0.9039270  0.6669610
## 9B  0.9731131  0.7400826
## 10A 1.3104102  1.8409287
## 10B 0.8232952  1.5216074
## 11A 1.1531298  1.0487371
## 11B 0.9428981  0.6016632
## 12B 0.9763376  0.9522037
## 13A 0.9581426  0.8811937
## 13B 1.0951320  0.9271929
## 14A 1.1434317  0.9603831
## 14B 0.9764873  1.0953167
## 15A 0.8995054  1.0960068
## 15B 0.9328239  0.6756996
## 16A 1.0068353  0.6781485
## 16B 1.0832678  1.0285193
## 17A 0.8555823  1.0587031
## 17B 0.9493854  0.8953877
## 18A 0.8829000  0.8301960
## 18B 1.0480114  0.7487345
## 19B 0.9152127  1.0281717
## 20A 1.2034111 -0.1153156
## 20B 0.9906489  0.6316683
## 22B 1.0260031  1.5225095
## 23A 1.1521094  1.0418286
## 23B 0.9142789  1.0946919
## 24A 1.0536291  0.4252301
## 24B 0.8921456  1.2039588
## 25A 1.3011811  1.0356674
## 25B 0.8437989  1.1576940
```

``` r
# Visualize Surrogate Variables to Ensure No Biological Information is Captured
# Prepare SV data for plotting
sv_data <- as.data.frame(svobj$sv)
colnames(sv_data) <- paste0("SV", seq_len(ncol(svobj$sv)))
sv_data <- cbind(d$samples, sv_data)

# Generate pairwise scatter plots of surrogate variables
sv_plots <- list()  # Store plots
sv_cols <- colnames(sv_data)[grep("^SV", colnames(sv_data))]

for (i in seq_along(sv_cols)) {
  for (j in seq_along(sv_cols)) {
    if (i < j) {  # Plot each pair only once
      sv_plots[[paste0("SV", i, "_SV", j)]] <- ggplot(sv_data, aes_string(x = sv_cols[i], y = sv_cols[j], 
                                                                          color = "Group", shape = "Timepoint")) +
        geom_point(size = 3, alpha = 0.8) +
        theme_minimal() +
        labs(title = paste("Surrogate Variables:", sv_cols[i], "vs", sv_cols[j]),
             x = paste("Surrogate Variable", i),
             y = paste("Surrogate Variable", j))
    }
  }
}
```

```
## Warning: `aes_string()` was deprecated in ggplot2 3.0.0.
## ℹ Please use tidy evaluation idioms with `aes()`.
## ℹ See also `vignette("ggplot2-in-packages")` for more information.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

``` r
# Print all scatter plots for manual inspection
for (plot_name in names(sv_plots)) {
  print(sv_plots[[plot_name]])
}
```

<img src="dge_files/figure-html/unnamed-chunk-2-2.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-3.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-4.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-5.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-6.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-7.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-8.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-9.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-10.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-11.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-12.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-13.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-14.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-15.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-16.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-17.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-18.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-19.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-20.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-21.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-22.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-23.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-24.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-25.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-26.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-27.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-28.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-29.png" width="672" />

``` r
cat("None of the surrogate variables captured variability related to disease group.\n")
```

```
## None of the surrogate variables captured variability related to disease group.
```

``` r
# Incorporate Surrogate Variables into the Design Matrix for downstream analysis
modSv <- cbind(design, svobj$sv)

# Run duplicateCorrelation to estimate correlation between repeated samples (blocking by horse)
corfit <- duplicateCorrelation(imputed_data_MinProb, modSv, block = as.factor(d$samples$Horse))
cat("Consensus correlation for repeated measures:", corfit$consensus.correlation, "\n") # 0.1
```

```
## Consensus correlation for repeated measures: 0.100393
```

``` r
# Fit the linear model with limma, accounting for repeated measures
fit <- lmFit(imputed_data_MinProb, modSv, block = as.factor(d$samples$Horse), correlation = corfit$consensus.correlation)
cat("Consensus correlation for repeated measures:", corfit$consensus.correlation, "\n") # 0.09
```

```
## Consensus correlation for repeated measures: 0.100393
```

``` r
# Ensure column names in the design matrix are syntactically valid
colnames(modSv) <- make.names(colnames(modSv))

# Define contrasts with the updated column names
con <- makeContrasts(
  Metformin_vs_AF_RA = conditionRA_Metformin_4months - conditionRA_Control_4months,
  AF_vs_Sham_RA = conditionRA_Control_4months - conditionRA_Sham_4months,
  Terminal_vs_Baseline_Control = conditionRA_Control_4months - conditionRA_Control_Baseline,
  Terminal_vs_Baseline_Metformin = conditionRA_Metformin_4months - conditionRA_Metformin_Baseline,
  Diff_Treatment = (conditionRA_Metformin_4months - conditionRA_Metformin_Baseline) - (conditionRA_Control_4months - conditionRA_Control_Baseline),
  Diff_Disease = (conditionRA_Control_4months - conditionRA_Control_Baseline) - (conditionRA_Sham_4months - conditionRA_Sham_Baseline),
  Baseline_Difference_Metf_vs_AF = conditionRA_Metformin_Baseline - conditionRA_Control_Baseline,
  Baseline_Difference_AF_vs_Sham = conditionRA_Control_Baseline - conditionRA_Sham_Baseline,
  levels = modSv)
con
```

```
##                                 Contrasts
## Levels                           Metformin_vs_AF_RA AF_vs_Sham_RA
##   conditionRA_Control_4months                    -1             1
##   conditionRA_Control_Baseline                    0             0
##   conditionRA_Metformin_4months                   1             0
##   conditionRA_Metformin_Baseline                  0             0
##   conditionRA_Sham_4months                        0            -1
##   conditionRA_Sham_Baseline                       0             0
##   X                                               0             0
##   X                                               0             0
##   X                                               0             0
##   X                                               0             0
##   X                                               0             0
##   X                                               0             0
##   X                                               0             0
##   X                                               0             0
##                                 Contrasts
## Levels                           Terminal_vs_Baseline_Control
##   conditionRA_Control_4months                               1
##   conditionRA_Control_Baseline                             -1
##   conditionRA_Metformin_4months                             0
##   conditionRA_Metformin_Baseline                            0
##   conditionRA_Sham_4months                                  0
##   conditionRA_Sham_Baseline                                 0
##   X                                                         0
##   X                                                         0
##   X                                                         0
##   X                                                         0
##   X                                                         0
##   X                                                         0
##   X                                                         0
##   X                                                         0
##                                 Contrasts
## Levels                           Terminal_vs_Baseline_Metformin Diff_Treatment
##   conditionRA_Control_4months                                 0             -1
##   conditionRA_Control_Baseline                                0              1
##   conditionRA_Metformin_4months                               1              1
##   conditionRA_Metformin_Baseline                             -1             -1
##   conditionRA_Sham_4months                                    0              0
##   conditionRA_Sham_Baseline                                   0              0
##   X                                                           0              0
##   X                                                           0              0
##   X                                                           0              0
##   X                                                           0              0
##   X                                                           0              0
##   X                                                           0              0
##   X                                                           0              0
##   X                                                           0              0
##                                 Contrasts
## Levels                           Diff_Disease Baseline_Difference_Metf_vs_AF
##   conditionRA_Control_4months               1                              0
##   conditionRA_Control_Baseline             -1                             -1
##   conditionRA_Metformin_4months             0                              0
##   conditionRA_Metformin_Baseline            0                              1
##   conditionRA_Sham_4months                 -1                              0
##   conditionRA_Sham_Baseline                 1                              0
##   X                                         0                              0
##   X                                         0                              0
##   X                                         0                              0
##   X                                         0                              0
##   X                                         0                              0
##   X                                         0                              0
##   X                                         0                              0
##   X                                         0                              0
##                                 Contrasts
## Levels                           Baseline_Difference_AF_vs_Sham
##   conditionRA_Control_4months                                 0
##   conditionRA_Control_Baseline                                1
##   conditionRA_Metformin_4months                               0
##   conditionRA_Metformin_Baseline                              0
##   conditionRA_Sham_4months                                    0
##   conditionRA_Sham_Baseline                                  -1
##   X                                                           0
##   X                                                           0
##   X                                                           0
##   X                                                           0
##   X                                                           0
##   X                                                           0
##   X                                                           0
##   X                                                           0
```

``` r
# Apply contrasts and run eBayes
fit <- contrasts.fit(fit, con)
```

```
## Warning in contrasts.fit(fit, con): row names of contrasts don't match col
## names of coefficients
```

``` r
fit <- eBayes(fit, robust = TRUE, trend = TRUE)

# Extract DGE results using the BH method for FDR correction
res <- list()  # List to store DGE results
for (i in colnames(con)) {
  res_tmp <- topTable(fit, coef = i, adjust.method = "BH", number = Inf)  # Get top table results
  res_tmp <- res_tmp[!is.na(res_tmp$t), ]  # Remove rows with NA values
  res_tmp$Contrast <- i #TODO: replaced this part >>> rep(i, nrow(res_tmp))  # Store the contrast name
  res[[i]] <- res_tmp  # Add to the results list
  
  # Print the number of differentially expressed genes based on adjusted p-values
  n_adj_pval <- nrow(res_tmp[res_tmp$adj.P.Val < 0.05, ])
  print(paste('Number of differentially expressed genes for', i, 'based on adjusted p-value (BH) =', n_adj_pval))
}
```

```
## [1] "Number of differentially expressed genes for Metformin_vs_AF_RA based on adjusted p-value (BH) = 0"
## [1] "Number of differentially expressed genes for AF_vs_Sham_RA based on adjusted p-value (BH) = 8"
## [1] "Number of differentially expressed genes for Terminal_vs_Baseline_Control based on adjusted p-value (BH) = 359"
## [1] "Number of differentially expressed genes for Terminal_vs_Baseline_Metformin based on adjusted p-value (BH) = 394"
## [1] "Number of differentially expressed genes for Diff_Treatment based on adjusted p-value (BH) = 1"
## [1] "Number of differentially expressed genes for Diff_Disease based on adjusted p-value (BH) = 10"
## [1] "Number of differentially expressed genes for Baseline_Difference_Metf_vs_AF based on adjusted p-value (BH) = 1"
## [1] "Number of differentially expressed genes for Baseline_Difference_AF_vs_Sham based on adjusted p-value (BH) = 42"
```

``` r
# Combine all results into a single data frame
res_all <- do.call(rbind, res)

# Map Gene Names Manually
res_all$GeneName <- sapply(seq_len(nrow(res_all)), function(i) {
  gsub(paste0("^", res_all$Contrast[i], "\\."), "", rownames(res_all)[i])
})

# Split results by contrast for easier output
res_split <- split(res_all, res_all$Contrast)

# Optionally, save the results to files
 openxlsx::write.xlsx(x = res_split, file = "../../../../Timecourse/analysis/01_dge/output/dge_results.xlsx", asTable = TRUE)
 data.table::fwrite(x = res_all, file = "../../../../Timecourse/analysis/01_dge/output/dge_results.gz", sep = "\t")

# Visualize Results with P-value Histograms
for (contrast in names(res_split)) {
  p_values <- res_split[[contrast]]$P.Value
  ggplot(data = data.frame(P.Value = p_values), aes(x = P.Value)) +
    geom_histogram(breaks = seq(0, 1, by = 0.05), color = "black", fill = "lightblue") +
    theme_minimal() +
    labs(title = paste("P-value Histogram for", contrast),
         x = "P-value",
         y = "Frequency") +
    xlim(0, 1)  # Set x-axis limits
  print(last_plot())
}
```

<img src="dge_files/figure-html/unnamed-chunk-2-30.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-31.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-32.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-33.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-34.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-35.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-36.png" width="672" /><img src="dge_files/figure-html/unnamed-chunk-2-37.png" width="672" />


### Volcano plots

``` r
volcano_plots <- list()
for (i in names(res)){
  volcano_plots[[i]] <- ggVolcano(x = res[[i]], 
                                  fdr = 0.05,
                                  fdr.column = "adj.P.Val", 
                                  pvalue.column = "P.Value", 
                                  logFC = 0, 
                                  logFC.column = "logFC", 
                                  text.size = 2) + 
    theme_bw(base_size = 10) + 
    ggtitle(i)
}
```

```
## Warning in max(x[get(fdr.column) <= .fdr][, get(pvalue.column)]): no
## non-missing arguments to max; returning -Inf
```

``` r
patchwork::wrap_plots(volcano_plots, ncol = 3)
```

<img src="dge_files/figure-html/volcano_plots-1.png" width="960" />

## Volcano Plot for Publicaiton

``` r
library(ggrepel)

# Create a named vector for ENSEMBL to Genes mapping
ensembl_to_Genes <- setNames(annot_reordered$Genes, annot_reordered$ENSEMBL)

# Map ENSEMBL IDs to Gene Names in the `res` list
for (contrast_name in names(res)) {
  # Ensure the dataframe has ENSEMBL IDs as rownames
  if (!"ENSEMBL" %in% colnames(res[[contrast_name]])) {
    res[[contrast_name]]$ENSEMBL <- rownames(res[[contrast_name]])
  }
  
  # Map Gene Names using the annotation
  res[[contrast_name]]$Genes <- ensembl_to_Genes[res[[contrast_name]]$ENSEMBL]
  
  # Replace NA values in Genes with ENSEMBL IDs (to ensure plotting works even if some gene names are missing)
  res[[contrast_name]]$Genes[is.na(res[[contrast_name]]$Genes)] <- res[[contrast_name]]$ENSEMBL[is.na(res[[contrast_name]]$Genes)]
}

# Load the volcano plot helper function
source("volcano_helpers.R")

# Create lists for volcano plots with and without labels
volcano_plots_no_labels <- list()
volcano_plots_with_labels <- list()

# Iterate over each contrast and create volcano plots
for (contrast_name in names(res)) {
  # Ensure the Genes column is present for labeling
  if (!"Genes" %in% colnames(res[[contrast_name]])) {
    res[[contrast_name]]$Genes <- sapply(rownames(res[[contrast_name]]), function(x) gsub(".*_", "", x))
  }
  
  # Generate volcano plots using the helper function
  volcano_plots <- create_custom_volcano_plot(
    df = res[[contrast_name]],
    logFC_col = "logFC",
    pvalue_col = "P.Value",
    adj_pvalue_col = "adj.P.Val",
    contrast_name = contrast_name,
    fc_cutoff = 0,  # Set fold-change cutoff for significance
    pvalue_cutoff = 0.05,  # Set p-value cutoff
    save_plot = TRUE, 
    output_path = "../output/",  # Adjust output path if needed
    show_labels = TRUE  # Generate both labeled and unlabeled plots
  )
  
  # Store the plots
  volcano_plots_no_labels[[contrast_name]] <- volcano_plots$No_Labels
  volcano_plots_with_labels[[contrast_name]] <- volcano_plots$With_Labels
}
```

```
## Warning: ggrepel: 297 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 323 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

``` r
# Combine and display volcano plots without labels
combined_volcano_no_labels <- patchwork::wrap_plots(volcano_plots_no_labels, ncol = 3)

# Combine and display volcano plots with labels
combined_volcano_with_labels <- patchwork::wrap_plots(volcano_plots_with_labels, ncol = 3)

# Display the combined volcano plots
print(combined_volcano_no_labels)
```

<img src="dge_files/figure-html/Volcano Plots for Publicaiton-1.png" width="672" />

``` r
print(combined_volcano_with_labels)
```

```
## Warning: ggrepel: 354 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 390 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

```
## Warning: ggrepel: 38 unlabeled data points (too many overlaps). Consider
## increasing max.overlaps
```

<img src="dge_files/figure-html/Volcano Plots for Publicaiton-2.png" width="672" />

``` r
# Print individual volcano plots with labels for key contrasts
print(volcano_plots_with_labels[["Diff_Treatment"]])
```

<img src="dge_files/figure-html/Volcano Plots for Publicaiton-3.png" width="672" />

``` r
print(volcano_plots_with_labels[["Diff_Disease"]])
```

<img src="dge_files/figure-html/Volcano Plots for Publicaiton-4.png" width="672" />

``` r
# ggsave("../output/volcano_Diff_Treatment_for_figure.png", (volcano_plots_with_labels[["Diff_Treatment"]]), 
#        dpi = 600, width = 4, height = 3, units = "in")
```
# Plotting One 

# Understanding Treatment Direction and Plotting Protein Abundance

``` r
# Understanding LogFC for "Diff_Treatment"
# Positive LogFC indicates:
# 1) Metformin group shows a relatively higher level compared to the control group.
# 2) This can occur if:
#    a) Metformin increases and control decreases.
#    b) Both increase, but metformin shows a greater increase.
#    c) Both decrease, but metformin shows a smaller decrease.

# Negative LogFC indicates:
# 1) Metformin group shows a relatively lower level compared to the control group.
# 2) This can occur if:
#    a) Metformin decreases and control increases.
#    b) Both decrease, but metformin shows a greater decrease.
#    c) Both increase, but metformin shows a smaller increase.

# Load ggpubr for enhanced visualization
library(ggpubr)

# Step 1: Define the function to create a violin plot with individual points and mean line
plot_protein_counts <- function(protein_of_interest) {
  # Check if the protein is present in the dataset
  if (protein_of_interest %in% rownames(data_se_norm)) {
    # Extract the normalized counts for the protein
    protein_counts <- assay(data_se_norm)[protein_of_interest,]
    
    # Create a data frame for plotting
    protein_df <- data.frame(
      SampleID = colnames(data_se_norm), # Extract sample IDs
      Count = protein_counts,            # Protein counts
      Condition = meta$Condition         # Experimental conditions
    )
    
    # Generate violin plot with jitter points and mean line
    p <- ggplot(protein_df, aes(x = Condition, y = Count, fill = Condition)) +
      geom_violin(trim = FALSE, alpha = 0.5) +                        # Violin plot
      geom_jitter(width = 0.2, size = 2, alpha = 0.7) +               # Add individual points
      stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "red") +  # Mean line
      labs(
        title = paste("Violin plot of normalized counts for", protein_of_interest),  # Plot title
        x = "Condition", y = "Normalized Abundance"                                  # Axis labels
      ) +
      theme_pubr() +                                                  # Apply a clean theme
      scale_fill_brewer(palette = "Set3")                             # Color palette

    print(p)  # Display the plot
  } else {
    # Output message if protein is not found
    cat("The protein", protein_of_interest, "is not present in the dataset.\n")
  }
}

# Step 2: Visualize selected proteins with significant treatment effects to understand, why it is upregulated and downregulated (see comment in beginning of chunk)

# A) Electron Transport Chain Genes - negatively regulated (blue) in "Diff_Treatment"
plot_protein_counts("NDUFA6")   # Associated with metformin in studies related to Aortic Aneurysms, increases less in metformin group
```

<img src="dge_files/figure-html/unnamed-chunk-3-1.png" width="672" />

``` r
plot_protein_counts("ATP5F1C")  # Important candidate gene for treatment effects, small decrease in metformin group
```

<img src="dge_files/figure-html/unnamed-chunk-3-2.png" width="672" />

``` r
plot_protein_counts("KARS1")    # Downregulated in response to treatment 
```

<img src="dge_files/figure-html/unnamed-chunk-3-3.png" width="672" />

``` r
# B) Proteasome Proteins - negatively regulated (blue) in "Diff_Treatment"
plot_protein_counts("PSMC5")    # Protein linked to proteasome function
```

<img src="dge_files/figure-html/unnamed-chunk-3-4.png" width="672" />

``` r
plot_protein_counts("PSMC4")    # Another candidate for proteasome-associated mechanisms
```

<img src="dge_files/figure-html/unnamed-chunk-3-5.png" width="672" />

``` r
plot_protein_counts("PSMD11")   # Related to proteasome degradation processes
```

<img src="dge_files/figure-html/unnamed-chunk-3-6.png" width="672" />

``` r
# C) Heat Shock Proteins (Hsp90 Family) - negatively regulated (blue) in "Diff_Treatment"
plot_protein_counts("SUGT1")    # Chaperone activity protein
```

```
## Warning: Removed 1 row containing non-finite outside the scale range
## (`stat_ydensity()`).
```

```
## Warning: Removed 1 row containing non-finite outside the scale range
## (`stat_summary()`).
```

```
## Warning: Removed 1 row containing missing values or values outside the scale range
## (`geom_point()`).
```

<img src="dge_files/figure-html/unnamed-chunk-3-7.png" width="672" />

``` r
plot_protein_counts("HSP90AA1") # Hsp90 protein linked to stress response
```

<img src="dge_files/figure-html/unnamed-chunk-3-8.png" width="672" />

``` r
plot_protein_counts("DYNC1H1")  # Motor protein involved in cellular transport
```

<img src="dge_files/figure-html/unnamed-chunk-3-9.png" width="672" />

``` r
# D) Detoxification of Reactive Oxygen Species (ROS) - positively regulated (red) in "Diff_Treatment"
plot_protein_counts("TXNRD2")   # Key enzyme in ROS detoxification
```

<img src="dge_files/figure-html/unnamed-chunk-3-10.png" width="672" />

``` r
plot_protein_counts("TXN")      # Thioredoxin, involved in redox regulation
```

<img src="dge_files/figure-html/unnamed-chunk-3-11.png" width="672" />

``` r
plot_protein_counts("SOD3")     # Superoxide dismutase, a primary ROS scavenger
```

<img src="dge_files/figure-html/unnamed-chunk-3-12.png" width="672" />

``` r
# Additional Proteins of Interest
plot_protein_counts("DDAH1")    # Dimethylarginine Dimethylaminohydrolase 1
```

<img src="dge_files/figure-html/unnamed-chunk-3-13.png" width="672" />

``` r
plot_protein_counts("COQ8A")    # Coenzyme Q8 homolog
```

<img src="dge_files/figure-html/unnamed-chunk-3-14.png" width="672" />

``` r
plot_protein_counts("RICTOR")   # Component of mTOR complex
```

```
## The protein RICTOR is not present in the dataset.
```

``` r
plot_protein_counts("YWHAE")    # 14-3-3 protein epsilon, signaling protein
```

<img src="dge_files/figure-html/unnamed-chunk-3-15.png" width="672" />

``` r
plot_protein_counts("TXNDC5")   # Protein disulfide isomerase, ROS response
```

<img src="dge_files/figure-html/unnamed-chunk-3-16.png" width="672" />

``` r
plot_protein_counts("GNAI2")
```

<img src="dge_files/figure-html/unnamed-chunk-3-17.png" width="672" />

### Focus on GNAI2

``` r
library(ggpubr)

# Define the function to create a violin plot with individual points and mean line
plot_protein_counts <- function(protein_of_interest) {
  # Check if the protein is present in the dataset
  if (protein_of_interest %in% rownames(data_se_norm)) {
    # Extract the normalized counts for the protein
    protein_counts <- assay(data_se_norm)[protein_of_interest,]
    
    # Create a data frame for plotting
    protein_df <- data.frame(
      SampleID = colnames(data_se_norm), # Extract sample IDs
      Count = protein_counts,            # Protein counts
      Condition = meta$Condition         # Experimental conditions
    )
    
    # Generate violin plot with jitter points and mean line
    p <- ggplot(protein_df, aes(x = Condition, y = Count, fill = Condition)) +
      geom_violin(trim = FALSE, alpha = 0.7, color = "black") +       # Violin plot with border
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.6, color = "black") +  # Individual points
      stat_summary(fun = mean, geom = "crossbar", width = 0.5, fatten = 2, color = "darkred") +  # Mean line
      labs(
        title = paste("Normalized Protein Counts for", protein_of_interest),  # Plot title
        x = "Condition", y = "Normalized Abundance"                           # Axis labels
      ) +
      theme_pubr() +                                                          # Apply a clean, publication-ready theme
      scale_fill_brewer(palette = "Set3") +                                   # Use a color palette with enough colors
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),     # Centered and bold title
        axis.title = element_text(size = 12, face = "bold"),                  # Bold axis titles
        axis.text = element_text(size = 10, color = "black"),                 # Custom axis text
        legend.position = "none"                                              # Remove legend
      )
    
    # Print the plot
    print(p)
  } else {
    # Output message if the protein is not found
    cat("The protein", protein_of_interest, "is not present in the dataset.\n")
  }
}

# Plot the protein counts for "GNAI2"
plot_protein_counts("GNAI2")
```

<img src="dge_files/figure-html/unnamed-chunk-4-1.png" width="672" />

# Session Info

``` r
sessionInfo()
```

```
## R version 4.4.2 (2024-10-31 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 10 x64 (build 19045)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=Danish_Denmark.utf8  LC_CTYPE=Danish_Denmark.utf8   
## [3] LC_MONETARY=Danish_Denmark.utf8 LC_NUMERIC=C                   
## [5] LC_TIME=Danish_Denmark.utf8    
## 
## time zone: Europe/Copenhagen
## tzcode source: internal
## 
## attached base packages:
## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] ggpubr_0.6.0                ggrepel_0.9.6              
##  [3] sva_3.54.0                  BiocParallel_1.40.0        
##  [5] genefilter_1.88.0           mgcv_1.9-1                 
##  [7] nlme_3.1-166                aamisc_0.1.5               
##  [9] fdrtool_1.2.18              Metrics_0.1.4              
## [11] SummarizedExperiment_1.36.0 Biobase_2.66.0             
## [13] GenomicRanges_1.58.0        GenomeInfoDb_1.42.0        
## [15] IRanges_2.40.0              S4Vectors_0.44.0           
## [17] BiocGenerics_0.52.0         MatrixGenerics_1.18.0      
## [19] DEP_1.28.0                  preprocessCore_1.68.0      
## [21] DEqMS_1.24.0                matrixStats_1.4.1          
## [23] RColorBrewer_1.1-3          missForest_1.5             
## [25] dplyr_1.1.4                 openxlsx_4.2.7.1           
## [27] patchwork_1.3.0             data.table_1.16.2          
## [29] ggplot2_3.5.1               stringr_1.5.1              
## [31] tibble_3.2.1                magrittr_2.0.3             
## [33] biomaRt_2.62.0              readxl_1.4.3               
## [35] readr_2.1.5                 edgeR_4.4.0                
## [37] limma_3.62.1                pacman_0.5.1               
## 
## loaded via a namespace (and not attached):
##   [1] splines_4.4.2               later_1.3.2                
##   [3] norm_1.0-11.1               filelock_1.0.3             
##   [5] R.oo_1.27.0                 cellranger_1.1.0           
##   [7] XML_3.99-0.17               lifecycle_1.0.4            
##   [9] httr2_1.0.6                 rstatix_0.7.2              
##  [11] doParallel_1.0.17           vroom_1.6.5                
##  [13] lattice_0.22-6              MASS_7.3-61                
##  [15] MultiAssayExperiment_1.32.0 backports_1.5.0            
##  [17] sass_0.4.9                  rmarkdown_2.29             
##  [19] jquerylib_0.1.4             yaml_2.3.10                
##  [21] httpuv_1.6.15               doRNG_1.8.6                
##  [23] zip_2.3.1                   MsCoreUtils_1.18.0         
##  [25] DBI_1.2.3                   abind_1.4-8                
##  [27] zlibbioc_1.52.0             R.utils_2.12.3             
##  [29] purrr_1.0.2                 HarmonicRegression_1.0     
##  [31] AnnotationFilter_1.30.0     itertools_0.1-3            
##  [33] rappdirs_0.3.3              sandwich_3.1-1             
##  [35] circlize_0.4.16             GenomeInfoDbData_1.2.13    
##  [37] annotate_1.84.0             MSnbase_2.32.0             
##  [39] ncdf4_1.23                  codetools_0.2-20           
##  [41] DelayedArray_0.32.0         DT_0.33                    
##  [43] xml2_1.3.6                  tidyselect_1.2.1           
##  [45] gmm_1.8                     shape_1.4.6.1              
##  [47] UCSC.utils_1.2.0            farver_2.1.2               
##  [49] gmp_0.7-5                   BiocFileCache_2.14.0       
##  [51] jsonlite_1.8.9              multtest_2.62.0            
##  [53] GetoptLong_1.0.5            Formula_1.2-5              
##  [55] survival_3.7-0              iterators_1.0.14           
##  [57] systemfonts_1.1.0           foreach_1.5.2              
##  [59] tools_4.4.2                 progress_1.2.3             
##  [61] ragg_1.3.3                  Rcpp_1.0.13-1              
##  [63] glue_1.8.0                  gridExtra_2.3              
##  [65] SparseArray_1.6.0           xfun_0.49                  
##  [67] qvalue_2.38.0               shinydashboard_0.7.2       
##  [69] withr_3.0.2                 BiocManager_1.30.25        
##  [71] fastmap_1.2.0               fansi_1.0.6                
##  [73] digest_0.6.37               mime_0.12                  
##  [75] R6_2.5.1                    textshaping_0.4.0          
##  [77] imputeLCMD_2.1              colorspace_2.1-1           
##  [79] RSQLite_2.3.8               R.methodsS3_1.8.2          
##  [81] hexbin_1.28.5               utf8_1.2.4                 
##  [83] tidyr_1.3.1                 generics_0.1.3             
##  [85] prettyunits_1.2.0           PSMatch_1.10.0             
##  [87] httr_1.4.7                  htmlwidgets_1.6.4          
##  [89] S4Arrays_1.6.0              pkgconfig_2.0.3            
##  [91] NISTunits_1.0.1             gtable_0.3.6               
##  [93] blob_1.2.4                  ComplexHeatmap_2.22.0      
##  [95] impute_1.80.0               XVector_0.46.0             
##  [97] htmltools_0.5.8.1           carData_3.0-5              
##  [99] rain_1.40.0                 MALDIquant_1.22.3          
## [101] ProtGenerics_1.38.0         clue_0.3-66                
## [103] scales_1.3.0                tmvtnorm_1.6               
## [105] png_0.1-8                   knitr_1.49                 
## [107] rstudioapi_0.17.1           tzdb_0.4.0                 
## [109] reshape2_1.4.4              rjson_0.2.23               
## [111] curl_6.0.1                  cachem_1.1.0               
## [113] zoo_1.8-12                  GlobalOptions_0.1.2        
## [115] parallel_4.4.2              AnnotationDbi_1.68.0       
## [117] mzID_1.44.0                 vsn_3.74.0                 
## [119] pillar_1.9.0                grid_4.4.2                 
## [121] vctrs_0.6.5                 pcaMethods_1.98.0          
## [123] promises_1.3.0              randomForest_4.7-1.2       
## [125] car_3.1-3                   dbplyr_2.5.0               
## [127] xtable_1.8-4                cluster_2.1.6              
## [129] evaluate_1.0.1              mvtnorm_1.3-2              
## [131] cli_3.6.3                   locfit_1.5-9.10            
## [133] compiler_4.4.2              rlang_1.1.4                
## [135] crayon_1.5.3                rngtools_1.5.2             
## [137] ggsignif_0.6.4              labeling_0.4.3             
## [139] QFeatures_1.16.0            affy_1.84.0                
## [141] plyr_1.8.9                  stringi_1.8.4              
## [143] assertthat_0.2.1            munsell_0.5.1              
## [145] Biostrings_2.74.0           lazyeval_0.2.2             
## [147] Matrix_1.7-1                hms_1.1.3                  
## [149] bit64_4.5.2                 KEGGREST_1.46.0            
## [151] statmod_1.5.0               shiny_1.9.1                
## [153] mzR_2.40.0                  broom_1.0.7                
## [155] igraph_2.1.1                memoise_2.0.1              
## [157] affyio_1.76.0               bslib_0.8.0                
## [159] bit_4.5.0
```
