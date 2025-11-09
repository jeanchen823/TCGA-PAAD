# TCGA-PAAD
This study identifies CXCL10 as a prognostic biomarker in pancreatic cancer through TCGA data mining, linking its expression to immune-rich tumor microenvironments and improved patient survival.
Bioinformatics Methods Used

# 1. ESTIMATE Algorithm
Computes ImmuneScore, StromalScore, and ESTIMATEScore
Estimates the proportion of immune and stromal components in tumor microenvironment
Higher scores indicate more immune/stromal content
# 2. CIBERSORT
Deconvolutes bulk RNA-seq data to estimate 22 immune cell type proportions
Uses leukocyte signature matrix (LM22) as reference
Provides tumor-infiltrating immune cell (TIC) profiles
# 3. Differential Expression Analysis
limma package for identifying differentially expressed genes (DEGs)
Comparison between high-score vs low-score groups
Threshold: |log2FC| > 1, FDR < 0.05
# 4. Survival Analysis
Kaplan-Meier curves with log-rank test
Univariate COX regression for survival-associated genes
Compares high vs low expression groups
# 5. Protein-Protein Interaction (PPI) Network
Constructed using STRING database
Visualized with Cytoscape software
Interaction confidence threshold > 0.95
# 6. Functional Enrichment Analysis
GO (Gene Ontology) - biological process, molecular function, cellular component
KEGG (Kyoto Encyclopedia of Genes and Genomes) - pathway enrichment
Uses clusterProfiler R package
# 7. Gene Set Enrichment Analysis (GSEA)
Analyzes Hallmark and C7 immunologic gene sets from MSigDB
Identifies enriched pathways in high/low BTK expression groups
Threshold: NOM p < 0.05, FDR q < 0.06
# 8. Statistical Tests
Wilcoxon rank sum test - for two-group comparisons
Kruskal-Wallis test - for multiple group comparisons
Pearson correlation - for correlation analysis
# 9. Data Visualization
Heatmaps (pheatmap package)
Venn diagrams - for intersection analysis
Violin plots and scatter plots - for distribution and correlation
# 10. Data Source
TCGA database -  samples (tumor,normal)
RNA-seq data (level 3) and clinical information
