library(DESeq2)
library(dplyr)
library(tidyverse)
library(ggupset)
library(clusterProfiler)
library(org.Hs.eg.db)  # Annotation for human genes 
library(enrichplot)
library(glmnet)
library(pathview)
library(survival)
library(pheatmap)
library(plotly)
library(survminer)
library(ggplot2)
set.seed(0)


# RNA-seq data
rna_seq <- read.table("data_mrna_seq_v2_rsem.txt", header = TRUE, sep = "\t")

# Clinical Data
clinical_data <- read.table("data_clinical_patient.txt", header = TRUE, sep = "\t")

# CNA Data
cna_data <- read.table("data_cna.txt", header = TRUE, sep = "\t")

dim(rna_seq)     
dim(clinical_data) 
dim(cna_data)  

# What we want to do next is fill in some of the blanks in the `Hugo_Symbol` column of the `rna_seq` data.
# We're gonna use `cna_data`, since it doesn't have any blank spaces.

blank_values <- rna_seq$Hugo_Symbol == ""
sum(blank_values)

blank_values2 <- cna_data$Hugo_Symbol == ""
sum(blank_values2)

# We match Entrez_Gene_Id from rna_seq (blanks) to cna_data...
matched_symbols <- match(rna_seq$Entrez_Gene_Id[blank_values], cna_data$Entrez_Gene_Id)

# ...and replace only non-NA matches
valid_matches <- !is.na(matched_symbols)  # Identify valid indices
rna_seq$Hugo_Symbol[blank_values][valid_matches] <- cna_data$Hugo_Symbol[matched_symbols[valid_matches]]

# We replace NAs created with blanks...
rna_seq$Hugo_Symbol[is.na(rna_seq$Hugo_Symbol)] <- ""

# ...and verify result
sum(rna_seq$Hugo_Symbol == "")

# We managed to fill in only 2 blanks, but that's okay.

# Next we remove the unannotated rows, and we also deal with duplicates. 
# For the latter, we keep the row with the highest counts, as it's more useful for the downstream analyses. 
# Another option would be to sum all the counts in the duplicated rows and consider them as an only gene. 
# But that is not necessarily the best approach, as genes might have isoforms that slightly differ from each other, 
# so putting them all together wouldn't make too much sense from a biological standpoint.

rna_seq <- rna_seq[!(is.na(rna_seq$Hugo_Symbol) | rna_seq$Hugo_Symbol == ""), ]

rna_seq <- rna_seq %>%
  group_by(Hugo_Symbol) %>%
  slice_max(rowSums(across(starts_with("TCGA")))) %>%
  ungroup()

sum(is.na(rownames(rna_seq)))       # Should be both 0
sum(duplicated(rownames(rna_seq))) 

# Next we reformat the IDs in the `rna_seq` data and in `cna` data, to ensure consistency with the `clinical_data`.
colnames(rna_seq)[3:ncol(rna_seq)] <- gsub("\\.", "-", colnames(rna_seq)[3:ncol(rna_seq)])
colnames(rna_seq)[3:ncol(rna_seq)] <- gsub("-01$", "", colnames(rna_seq)[3:ncol(rna_seq)])

colnames(cna_data)[3:ncol(cna_data)] <- gsub("\\.", "-", colnames(cna_data)[3:ncol(cna_data)])
colnames(cna_data)[3:ncol(cna_data)] <- gsub("-01$", "", colnames(cna_data)[3:ncol(cna_data)])

# and we verify the reformatted IDs
print(head(colnames(rna_seq)[3:ncol(rna_seq)]))
print(head(colnames(cna_data)[3:ncol(cna_data)]))

# Removing columns with ALL NA values from the datasets
rna_seq <- rna_seq[, colSums(is.na(rna_seq)) < nrow(rna_seq)]
cna_data <- cna_data[, colSums(is.na(cna_data)) < nrow(cna_data)]
clinical_data <- clinical_data[, colSums(is.na(clinical_data)) < nrow(clinical_data)]

# Verifying successful deletion
sum(colSums(is.na(rna_seq)) == nrow(rna_seq))  
sum(colSums(is.na(cna_data)) == nrow(cna_data)) 
sum(colSums(is.na(clinical_data)) == nrow(clinical_data))

# Removing genes with all zero counts from the RNA dataset
rna_seq <- rna_seq[rowSums(rna_seq[, 3:ncol(rna_seq)]) > 0, ]
sum(rowSums(rna_seq[, 3:ncol(rna_seq)]) == 0)

# First we extract patient IDs
rna_seq_patient_ids <- colnames(rna_seq)[3:ncol(rna_seq)]
cna_patient_ids <- colnames(cna_data)[3:ncol(cna_data)]
clinical_patient_ids <- clinical_data$PATIENT_ID

# Then find common patient IDs across all three datasets
common_patient_ids <- Reduce(intersect, list(rna_seq_patient_ids, cna_patient_ids, clinical_patient_ids))

print(length(common_patient_ids))
print(head(common_patient_ids))

# We find columns corresponding to common IDs
rna_seq_columns <- which(colnames(rna_seq) %in% common_patient_ids)

# Subset RNA-seq data
rna_seq <- rna_seq[, c(1, 2, rna_seq_columns)]  # We wnat to retain the first two columns 
# (Hugo_Symbol, Entrez_Gene_Id) 

# Same for CNA data...
cna_columns <- which(colnames(cna_data) %in% common_patient_ids)
cna_data <- cna_data[, c(1, 2, cna_columns)]  

# ...and for clinical data
clinical_data <- clinical_data[clinical_data$PATIENT_ID %in% common_patient_ids, ]

# Check dimensions
print(dim(rna_seq))     
print(dim(cna_data))     
print(dim(clinical_data)) 

# Checking patient IDs alignment, both lines should return TRUE
all(colnames(rna_seq)[3:ncol(rna_seq)] == colnames(cna_data)[3:ncol(cna_data)])
all(colnames(rna_seq)[3:ncol(rna_seq)] %in% clinical_data$PATIENT_ID)          

dim(rna_seq[, -c(1, 2)])   # Dimensions of numeric matrix
length(rna_seq$Hugo_Symbol)

dim(rna_seq)       # 20232  1070
dim(clinical_data) # 1068   37
dim(cna_data)      # 25128  1070

# Getting ERBB2 row from the CNA data subset
erbb2_cna <- as.numeric(cna_data[cna_data$Hugo_Symbol == "ERBB2", 
                                 3:ncol(cna_data)])

metadata <- data.frame(
  row.names = colnames(cna_data)[3:ncol(cna_data)],  # Set patient IDs as row names
  ERBB2_Status = ifelse(erbb2_cna > 0, 1, 0)         # HER2 status: 1 = Amplified, 
                                                     # 0 = Not Amplified
)
metadata$ERBB2_Status <- factor(metadata$ERBB2_Status)

# Before normalising we need integer gene counts, so we round the `rna_seq` values.

count_matrix <- round(as.matrix(rna_seq[, -c(1,2)])) 
rownames(count_matrix) <- rna_seq$Hugo_Symbol

count_matrix[is.na(count_matrix)] = 0  # Impute NAs with zeros
count_matrix[count_matrix<0] = 0

smallestGroupSize <- 3
keep <- rowSums(count_matrix >= 10) >= smallestGroupSize
count_matrix <- count_matrix[keep,]

# Check if all column names of count_matrix are present in metadata...
all(colnames(count_matrix) %in% rownames(metadata))
# ...and viceversa
all(rownames(metadata) %in% colnames(count_matrix)) 

dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata, 
  design = ~ ERBB2_Status 
)

dds <- DESeq(dds)

res <- results(dds, contrast = c("ERBB2_Status", "1", "0"))
res <- res[!is.na(rownames(res)), ]

# `results()` compares the two statuses and outputs:
# - log2FoldChange: the magnitude of expression changes between the two groups
# - pvalue: significance of differential expression for each gene
# - padj: adjusted p-values (multiple testing correction)

# Top 10 genes with adjusted p-values (STATISTICALLY SIGNIFICANT)
top10_genes <- res[order(res$padj, na.last = NA)[1:10], ] # orders genes based on their adjusted p-value
# in ascending order
# Also removes genes with NA values for padj 

top10_genes <- as.data.frame(top10_genes)
rownames(top10_genes) #<- rownames(res)[order(res$padj, na.last = NA)[1:10]]

top10_genes

# Some checks
sum(is.na(rownames(res)))  
sum(is.na(res$log2FoldChange))  

max_abs_log2FoldChange <- max(abs(res$log2FoldChange), na.rm = TRUE)
max_abs_log2FoldChange

# BIOLOGICALLY SIGNIFICANT GENES
res_sig <- res[res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]

top10_genes_fc <- res_sig[order(abs(res_sig$log2FoldChange), decreasing = TRUE)[1:10], ]

top10_genes_fc <- as.data.frame(top10_genes_fc)

rownames(top10_genes_fc) 
top10_genes_fc

# PATHWAY ENRICHMENT
Results <- res[!is.na(res$padj) & res$padj < 0.05, ]

# Overexpressed and underexpressed genes
DE_over <- rownames(Results[Results$log2FoldChange > 0, ])
DE_under <- rownames(Results[Results$log2FoldChange < 0, ])

# We map overexpressed and underexpressed genes to Entrez IDs
entrez_under <- bitr(DE_under, 
                     fromType = "SYMBOL", 
                     toType = "ENTREZID", 
                     OrgDb = org.Hs.eg.db)

entrez_over <- bitr(DE_over, 
                    fromType = "SYMBOL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db)

# KEGG enrichment for overexpressed genes

kegg_over <- enrichKEGG(
  gene = entrez_over$ENTREZID,
  organism = "human",           
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff  = 0.05
)

# KEGG enrichment for underexpressed genes

kegg_under <- enrichKEGG(
  gene = entrez_under$ENTREZID,
  organism = "human",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff  = 0.05
)

# Dotplots
png("kegg_under_dotplot.png", width = 1000, height = 800)
dotplot(kegg_under, showCategory=15) + 
  ggtitle("Kegg Pathway Enrichment Under Expressed")
dev.off()

png("kegg_over_dotplot.png", width = 1000, height = 800)
dotplot(kegg_over, showCategory=15) + 
  ggtitle("Kegg Pathway Enrichment Over Expressed")
dev.off()

# Pairwise similarity and treeplot
png("kegg_under_treeplot.png", width = 1000, height = 800)
kegg_enrichment_under_pw <- pairwise_termsim(kegg_under)
treeplot(kegg_enrichment_under_pw) + 
  ggtitle("KEGG Pathway Enrichment (Underexpressed Genes)") +
  theme(legend.position = "bottom", legend.direction = "vertical")
dev.off()

png("kegg_over_treeplot.png", width = 1000, height = 800)
kegg_enrichment_over_pw <- pairwise_termsim(kegg_over)
treeplot(kegg_enrichment_over_pw) + 
  ggtitle("KEGG Pathway Enrichment (Overexpressed Genes)") +
  theme(legend.position = "bottom", legend.direction = "vertical")
dev.off()

dds <- dds[rowSums(counts(dds)) > 10, ] # We keep genes with more than 10 counts

# Variance stabilising transformation
vst_data <- vst(dds, blind = FALSE)

vst_matrix <- assay(vst_data)  # Extract the VST-transformed matrix

# Plot the distribution of VST values
hist(vst_matrix, breaks = 50, main = "Distribution of VST Expression Values",
     xlab = "VST Expression Values", col = "skyblue")

# The skewed distribution indicates that most genes have relatively low expression values, with fewer genes exhibiting higher expression.
# This reflects the biological reality that only a subset of genes are highly expressed or activated in the analysed samples.

top_DE <- rownames(res)[order(res$padj, na.last = NA)][1:10]
vst_DE <- assay(vst_data)[top_DE,]#[1:10],]

if (!all(rownames(metadata) == colnames(vst_matrix))) {
  stop("Error: Row names of metadata do not match column names of vst_matrix!")
} else {
  message("Row names of metadata match column names of vst_matrix.")
}

annotation_colors <- list(ERBB2_Status = c("1" = "#E74C3C", "0" = "#27AE60"))
annotation_col = data.frame(ERBB2_Status = metadata$ERBB2_Status)
rownames(annotation_col) <- colnames(vst_DE)

png("heatmap_top_DE_genes.png", width = 2000, height = 1200)

pheatmap(
  vst_DE,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Heatmap of Top 20 DE Genes"
)

dev.off()

# This heatmap visualises the expression levels of the top 10 differentially expressed (DE) genes across samples, with annotations indicating ERBB2 status. 
# Each row represents a gene, and each column corresponds to a sample (patient). 
# The color gradient corresponds to normalised expression values, with red indicating higher expression relative to the mean, blue indicating lower expression, and yellow representing expression levels near the mean. 
# Samples and genes are hierarchically clustered to identify patterns of co-expression and similarity among samples. 
# For example, ERBB2 and GRB7, which are located within the amplicon on chromosome 17 and are known drivers of HER2-positive breast cancer, show significant overexpression in HER2-amplified samples.

top_DE_fc <- rownames(results)[order(abs(Results$log2FoldChange), decreasing = TRUE)][1:10]
vst_DE_fc <- assay(vst_data)[top_DE_fc, ]

if (!all(rownames(metadata) == colnames(vst_matrix))) {
  stop("Error: Row names of metadata do not match column names of vst_matrix!")
} else {
  message("Row names of metadata match column names of vst_matrix.")
}

annotation_colors <- list(ERBB2_Status = c("1" = "#E74C3C", "0" = "#27AE60"))
annotation_col <- data.frame(ERBB2_Status = metadata$ERBB2_Status)
rownames(annotation_col) <- colnames(vst_DE_fc)

png("heatmap_top_log2FC_genes.png", width = 2000, height = 1200)

pheatmap(
  vst_DE_fc,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "row",
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  show_colnames = FALSE,
  show_rownames = TRUE,
  main = "Heatmap of Top 10 Genes by log2FoldChange"
)

dev.off()

# PCA
vst <- vst(dds)
pca_data <- plotPCA(vst, intgroup=c("ERBB2_Status"), returnData = TRUE) 

# Percentage variance explained by PC1 and PC2
percentVar <- round(100 * attr(pca_data, "percentVar"))

ggplot(pca_data, aes(x = PC1, y = PC2, color = ERBB2_Status)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust point size and transparency
  theme_minimal(base_size = 14) +     # Use a clean minimal theme
  scale_color_manual(values = c("1" = "#E74C3C", "0" = "#27AE60")) +  # Custom color palette
  labs(
    title = "PCA Plot of VST-Transformed Data",
    subtitle = "Differentiating HER2 Amplified vs Not Amplified Tumors",
    x = paste0("PC1: ", percentVar[1], "% Variance"),
    y = paste0("PC2: ", percentVar[2], "% Variance"),
    color = "HER2 Status"
  ) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "right"
  )

# SURVIVAL MODEL

DE_genes <- rownames(res_sig)

# Subset VST values for significant DEGs
vst_DE_genes <- vst_matrix[DE_genes, ]

all(rownames(metadata) == colnames(vst_DE_genes)) # Should return TRUE

# Align clinical data with VST data
clinical_data_filtered <- clinical_data[clinical_data$PATIENT_ID %in% rownames(metadata), ]

clinical_data_filtered$OS_STATUS <- ifelse(
  clinical_data_filtered$OS_STATUS == "1:DECEASED", 1, 0)

# Create a survival data object
survival_data <- data.frame(
  SurvivalTime = as.numeric(clinical_data_filtered$OS_MONTHS),   # Overall survival time in months
  SurvivalStatus = as.numeric(clinical_data_filtered$OS_STATUS), # 1 = death, 0 = alive
  row.names = clinical_data_filtered$PATIENT_ID
)

survival_data$SurvivalTime[which(survival_data$SurvivalTime <= 0)] <- 1e-10 # so the Cox model can work
# Ensure the survival data matches the VST sample order
vst_DE_genes <- vst_DE_genes[, rownames(survival_data)]

# Convert VST matrix to predictors (x)
x <- t(vst_DE_genes)  # Samples as rows, genes as columns

# Survival response (y): Combine survival time and status into a Surv object
y <- Surv(survival_data$SurvivalTime, survival_data$SurvivalStatus)

# We check for NAs in SurvivalTime and SurvivalStatus
sum(is.na(survival_data$SurvivalTime))  # Count missing survival times
sum(is.na(survival_data$SurvivalStatus))  # Count missing survival statuses

# Fit Lasso-Cox model 
lasso_cox_model <- glmnet(
  x = x, 
  y = y, 
  family = "cox",    # Cox proportional hazards model
  alpha = 1          # Lasso regularization
)

# Cross-validation to select the best lambda
cv_fit <- cv.glmnet(
  x = x, 
  y = y, 
  family = "cox", 
  alpha = 1
)

# Best lambda
best_lambda <- cv_fit$lambda.min

# Extract coefficients at the best lambda
coef <- coef(lasso_cox_model, s = best_lambda)
selected_genes <- rownames(coef)[which(coef != 0)]
selected_genes

# Filter genes with non-zero coefficients
non_zero_coef <- coef[coef != 0, , drop = FALSE]

# Separate positive and negative coefficients
positive_genes <- non_zero_coef[non_zero_coef > 0, , drop = FALSE]
negative_genes <- non_zero_coef[non_zero_coef < 0, , drop = FALSE]

positive_genes_df <- data.frame(Gene = rownames(positive_genes), Coefficient = positive_genes[, 1])
negative_genes_df <- data.frame(Gene = rownames(negative_genes), Coefficient = negative_genes[, 1])

cat("Positive Coefficients:\n")
print(rownames(positive_genes_df))

cat("Negative Coefficients:\n")
print(rownames(negative_genes_df))

# Predict risk scores
risk_scores <- predict(lasso_cox_model, newx = x, s = best_lambda, type = "link")

# Stratify into high- and low-risk groups
risk_group <- ifelse(risk_scores > median(risk_scores), "High Risk", "Low Risk")

# Add risk group to survival data
survival_data$RiskGroup <- risk_group

# Kaplan-Meier curve
fit <- survfit(Surv(SurvivalTime, SurvivalStatus) ~ RiskGroup, data = survival_data)

plot(fit, col = c("#E74C3C", "#27AE60"), main = "Survival Curve by Risk Group", xlab = "Time (months)", ylab = "Survival Probability")

legend("topright", 
       legend = c("High Risk", "Low Risk"), 
       col = c("#E74C3C", "#27AE60"), 
       lty = 1)


