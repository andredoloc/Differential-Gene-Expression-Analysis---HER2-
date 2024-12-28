# Differential-Gene-Expression-Analysis - HER2+

This project is the final assignment for the Bio Principles and Cellular Organisation course (ANAT40040). 
It focuses on analyzing differential gene expression, pathway enrichment, and survival outcomes in HER2-positive breast cancer samples using RNA-seq data.


## **Code Explanation**

### **1. Setup and Libraries**
- The script starts by loading the R libraries for data manipulation, visualization, and bioinformatics analyses:
  - **`DESeq2`**: for differential expression analysis.
  - **`clusterProfiler` and `org.Hs.eg.db`**: for pathway enrichment.
  - **`ggplot2`, `pheatmap`**: for creating plots and heatmaps.
  - **`glmnet` and `survival`**: for survival analysis using Lasso-Cox regression.

### **2. Data Preprocessing**
1. **Loading Data**
   - RNA-seq (`data_mrna_seq_v2_rsem.txt`), clinical (`data_clinical_patient.txt`), and CNA (`data_cna.txt`) datasets are imported for analysis.
2. **Cleaning and Formatting**
   - Missing gene names in RNA-seq data are filled using CNA data.
   - Duplicate rows are resolved by retaining rows with the highest expression values.
   - Patient IDs in all datasets are reformatted for consistency.
3. **Subset Common Patients**
   - Only samples present in all three datasets (RNA-seq, CNA, clinical) are retained.

### **3. Differential Expression Analysis**
- Uses **`DESeq2`** to compare HER2-amplified vs. non-amplified samples:
  - Genes with low counts are removed before analysis.
  - Outputs:
    - `log2FoldChange`: magnitude of expression differences.
    - `padj`: adjusted p-values for statistical significance.
  - Top genes are selected by:
    - Statistical significance (`padj`).
    - Biological relevance (absolute `log2FoldChange` > 1).

### **4. Pathway Enrichment Analysis**
- Identifies biological pathways associated with differentially expressed genes using **KEGG**.
  - **Underexpressed Genes**: immune-related and signaling pathways.
  - **Overexpressed Genes**: DNA replication, cell cycle.
- Generates:
  - **Dotplots**: visualise pathway significance and gene counts.
  - **Treeplots**: show relationships between enriched pathways.

### **5. Principal Component Analysis (PCA)**
- Visualizes variability in gene expression across samples:
  - PCA identifies clusters based on HER2 amplification.
  - The percentage of variance explained by the first two components is annotated on the plot.

### **6. Heatmap Visualization**
- Heatmap 1: Top 10 statistically significant genes.
- Heatmap 2: Top 10 genes by `log2FoldChange`.
- Annotations indicate HER2 status for each sample.

### **7. Survival Analysis**
- Uses Lasso-Cox regression to identify genes influencing survival:
  - Survival data combines time (`OS_MONTHS`) and status (`OS_STATUS`).
  - Lasso regularization selects the most informative genes.
  - Outputs:
    - Genes with positive coefficients (higher risk).
    - Genes with negative coefficients (lower risk).
  - Risk scores stratify patients into high- and low-risk groups.
  - Kaplan-Meier curves visualize survival probabilities for each group.
