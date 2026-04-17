# Genealogical Analysis of Grapevine Cultivars R scripts
R pipeline for analysis of SNP data (heterozygosity, PCA, IBD, clustering) based on the SNPRelate package (Zheng et al., 2012).
#### To run this scripts, you will need R installed along with the following packages:
#### install.packages("BiocManager")
#### BiocManager::install("SNPRelate")
#### install.packages("writexl")
## 01_heterozygosity_calc 
Converts VCF data into GDS format, performs LD pruning, and calculates individual observed (Ho) and expected (He) heterozygosity.
## 02_filter_gds
Cleans the GDS dataset based on the results from 01_heterozygosity_calc (e.g., removing samples with abnormal heterozygosity or high inbreeding). If any samples are removed, it automatically performs a new LD pruning for the remaining dataset.
## 03_pca_analysis
Performs Principal Component Analysis (PCA) on the cleaned GDS dataset. It extracts the genetic variance for the first two components and generates a publication-ready, high-resolution (600 DPI) TIFF plot using ggplot2 and ggrepel for non-overlapping labels.
## 04_ibd_kinship
Calculates Identity-By-Descent (IBD) and kinship coefficients using the KING-robust method, which is ideal for populations with structure/stratification. It automatically classifies the relationships between sample pairs (e.g., Clones, Parent-offspring, Full sibs, Unrelated) based on standard theoretical thresholds and exports the results to an Excel file.

