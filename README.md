#  NEWMODEL
## Overview
NEWMODEL is a tool that integrates single-cell transcriptomic data with eQTL reference panels and performs cell-type-specific TWAS via a Bayesian multi-task regression framework, taking individual-level genotype and bulk expression data as input and outputting gene-trait association statistics across cell types.
## Data Preparation
The detailed description of the data that need to be prepared for implementing the model is stated as below:
1.  **sc_data (single-cell expression matrix)**: A matrix with genes as rows and cells as columns from single-cell RNA-seq data, used to construct the cell-type reference atlas.
2.  **bulk_data (bulk tissue expression matrix)**: A matrix with genes as rows and samples as columns from bulk RNA-seq data, used for deconvolution to obtain cell-type-specific expression.
3.  **Genotype files in genotype_dir**: One CSV file per gene `{gene_name}_sum.csv` with samples as rows and SNPs as columns, used for mr.mash regression modeling.
4.  **Cell-type expression files in zcellgene_dir**: One CSV file per gene `{gene_name}.csv` with samples as rows and cell types as columns, containing cell-type-specific expression values deconvoluted by BayesPrism.
5.  **GWAS summary statistics (gwas_data)**: Must contain at least four columns: chromosome (CHROM), base position (POS), effect size (BETA), and standard error (SE), used to calculate Z-scores.
6.  **genechrom_vector (gene chromosome information)**: An RDS file storing a vector indicating the chromosome of each gene, used to match with GWAS data.
7.  **LD reference panel (V matrix)**: A SNP correlation matrix computed from genotype data using the `big_cor()` function, used as the denominator for Z-score calculation.
8.  **Output directory (output_dir/mrresult)**: Used to save mr.mash model fits`{gene_name}.RData`, containing `fit$mu1` (effect size matrix).
## Dependency
NEWMODEL is built based on R package "BayesPrism" and "mr.mash.alpha", we sincerely appreciate the authors of these packages.
