# Project
#   Differential expression analysis: Nucleus accumbens SST interneuron nuclei
#   Saline vs 7.5mg/kg cocaine (1 hour) - mouse
#
# Description
#   Script 03: QC + filtering + normalization + limma-voom DE + figures and tables

# Load libraries

library(SummarizedExperiment)
library(edgeR)
library(limma)
library(ggplot2)
library(pheatmap)

# Create output folders for figures and tables (if they don't exist)

dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

# Load prepared RSE object (includes: counts, group, assigned_gene_prop, lib_size)
rse <- readRDS("output/rse/rse_prepared.rds")

# Explore the object
rse

# Extract metadata for each sample
meta <- as.data.frame(colData(rse))

# Confirm size of samples per group
table(meta$group)


# Quality control (QC)

# QC- assigned gene proportion

# This metric approximates the proportion of reads assigned to genes by featureCounts
#prints a numeric summary (min, cuartiles, median, mean, max)
#Helps detecting outliers or samples with low quality
summary(meta$assigned_gene_prop)

# Boxplot by group (assigned_gene_prop)
p_qc_geneprop <- ggplot(meta, aes(x = group, y = assigned_gene_prop)) +
  geom_boxplot() +
  theme_bw(base_size = 14) +
  labs(x = "Group", y = "Assigned gene proportion")

ggsave(
  filename = "output/figures/qc_assigned_gene_prop_boxplot.png",
  plot = p_qc_geneprop,
  width = 6,
  height = 4,
  dpi = 300
)

# QC-library size (from counts)
#Total number of reads counted for each sample (library size)
#is a key QC metric that can indicate potential issues with sequencing depth or sample quality.

summary(meta$lib_size)

#log10 transformation of library size for better visualization
p_qc_libsize <- ggplot(meta, aes(x = group, y = log_lib_size)) +
  geom_boxplot() +
  theme_bw(base_size = 14) +
  labs(x = "Group", y = "Log10(Library size + 1)")

ggsave(
  filename = "output/figures/qc_log_library_size_boxplot.png",
  plot = p_qc_libsize,
  width = 6,
  height = 4,
  dpi = 300
)


# Filtering genes and normalization with edgeR

# Create DGEList object from raw counts
# DGEList contains counts, gene annotations (rowData) and samples (colData)
dge <- DGEList(
  counts = assay(rse, "counts"),
  genes = rowData(rse),
  samples = meta
)

# Filter lowly expressed genes (robust approach)
#returns a logical vector for each gene
keep <- filterByExpr(dge, group = meta$group)
table(keep) #Number of genes retained after filtering

# Save number of genes before/after filtering
n_genes_before <- nrow(dge$counts)
dge <- dge[keep, , keep.lib.sizes = FALSE]
n_genes_after <- nrow(dge$counts)

#print summary of filtering results
cat("Genes before filtering:", n_genes_before, "\n")
cat("Genes after filtering: ", n_genes_after, "\n")
cat(
  "Percent retained:     ",
  round(n_genes_after / n_genes_before * 100, 2),
  "%\n",
  sep = ""
)

# TMM normalization
#Adjust for differences in libraru size/composition to compare between samples
dge <- calcNormFactors(dge)

#Converts counts to CPM and then log transforms
# Compute logCPM for exploratory plots
logCPM <- cpm(dge, log = TRUE)

# Density plot of logCPM by group
df_density <- data.frame(
  value = as.vector(logCPM),
  group = rep(meta$group, each = nrow(logCPM))
)

#Density plot of logCPM values for all samples
p_density <- ggplot(df_density, aes(x = value, color = group)) +
  geom_density() +
  theme_bw(base_size = 14) +
  labs(x = "logCPM", y = "Density")

ggsave(
  filename = "output/figures/qc_logcpm_density.png",
  plot = p_density,
  width = 6,
  height = 4,
  dpi = 300
)

# Sample-to-sample correlation heatmap (logCPM) uses al the filtered genes
cor_mat <- cor(logCPM)

#creates  a table to anotate the columns of the heatmap (group for sample)
annotation_col <- data.frame(group = meta$group)
row.names(annotation_col) <- colnames(logCPM)

#heatmap of sample correlations

pheatmap(
  cor_mat,
  annotation_col = annotation_col,
  filename = "output/figures/qc_sample_correlation_heatmap.png",
  width = 7,
  height = 6
)


# Differential expression with limma-voom

# Model matrix:
# We include assigned_gene_prop as a technical covariate to adjust for variability
design <- model.matrix(~ assigned_gene_prop + group, data = meta)
colnames(design)

# Voom transformation (also saves a mean-variance trend plot)
png(
  "output/figures/voom_mean_variance_trend.png",
  width = 1200,
  height = 900,
  res = 150
)
vGene <- voomWithQualityWeights(dge, design, plot = TRUE)
dev.off()

# Fit model and empirical Bayes moderation
#lmFit() fits a linear model for each gene
#eBayes() applies empirical Bayes moderation of variances to improve statistical inference
fit <- eBayes(lmFit(vGene, design))

# Extract DE results: effect of cocaine vs saline
# Since group baseline is saline, groupcocaine corresponds to cocaine - saline
de_results <- topTable(
  fit,
  coef = "groupcocaine",
  number = Inf,
  sort.by = "P"
)

# Save results
write.csv(
  de_results,
  "output/tables/de_results_cocaine_vs_saline.csv",
  row.names = FALSE
)

# Resumen nominal (comparabilidad con el paper; exploratorio)
n_nominal <- sum(de_results$P.Value < 0.05, na.rm = TRUE)
cat("Genes con p nominal < 0.05:", n_nominal, "\n")

# Guardar top genes por p nominal
top50_nominal <- de_results[order(de_results$P.Value), ][1:50, ]
write.csv(top50_nominal, "output/tables/top50_nominal_p.csv", row.names = FALSE)


#Count significant genes based on FDR and logFC thresholds

alpha <- 0.05
lfc_cut <- 1

# Conteo con FDR (principal)
sig_fdr <- de_results[
  de_results$adj.P.Val < alpha & abs(de_results$logFC) > lfc_cut,
]
n_up_fdr <- sum(sig_fdr$logFC > 0)
n_down_fdr <- sum(sig_fdr$logFC < 0)

cat("DEGs (FDR <", alpha, "and |logFC| >", lfc_cut, ")\n")
cat("Up in cocaine:", n_up_fdr, "\n")
cat("Up in saline: ", n_down_fdr, "\n")

# Conteo exploratorio con p nominal
sig_nom <- de_results[
  de_results$P.Value < alpha & abs(de_results$logFC) > lfc_cut,
]
n_up_nom <- sum(sig_nom$logFC > 0)
n_down_nom <- sum(sig_nom$logFC < 0)

cat("Candidatos (p nominal <", alpha, "and |logFC| >", lfc_cut, ")\n")
cat("Up in cocaine:", n_up_nom, "\n")
cat("Up in saline: ", n_down_nom, "\n")

#data frame for up/down counts plot
df_counts_nom <- data.frame(
  direction = c("Up in cocaine", "Up in saline"),
  n = c(n_up_nom, n_down_nom)
)

p_updown_nom <- ggplot(df_counts_nom, aes(x = direction, y = n)) +
  geom_col() +
  theme_bw(base_size = 14) +
  labs(x = NULL, y = "Number of candidate genes (p < 0.05 nominal)")

ggsave(
  filename = "output/figures/de_up_down_counts_nominal.png",
  plot = p_updown_nom,
  width = 6,
  height = 4,
  dpi = 300
)

# Exploratory plots: PCA and MDS (using voom-transformed expression)

# PCA using voom expression values

#vGene$E is the matrix of voom-transformed expression values
pca <- prcomp(t(vGene$E), scale. = TRUE)
#percentVar calculates the percentage of variance explained by each principal component
percentVar <- round(100 * (pca$sdev^2 / sum(pca$sdev^2)), 1)

#create dataframe with cordinates of PC1/PC2 for each sample and group for coloring
pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  group = meta$group
)

#PCA plot with axis labels showing %variance
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  theme_bw(base_size = 14) +
  labs(
    x = paste0("PC1: ", percentVar[1], "%"),
    y = paste0("PC2: ", percentVar[2], "%")
  )

#Save PCA plot
ggsave(
  filename = "output/figures/pca_plot.png",
  plot = p_pca,
  width = 6,
  height = 5,
  dpi = 300
)

# MDS plot (limma)
png("output/figures/mds_plot.png", width = 1200, height = 900, res = 150)
plotMDS(vGene$E, labels = meta$group, col = as.integer(meta$group))
legend(
  "topright",
  legend = levels(meta$group),
  col = 1:length(levels(meta$group)),
  pch = 16
)
dev.off()


# Volcano plot and MA plot

# Volcano plot
volcano_df <- de_results
volcano_df$neglog10p <- -log10(volcano_df$P.Value)

#clasificates genes based on FDR and logFC
volcano_df$significance <- "Not significant"
volcano_df$significance[
  volcano_df$adj.P.Val < alpha & volcano_df$logFC > lfc_cut
] <- "Up in cocaine"
volcano_df$significance[
  volcano_df$adj.P.Val < alpha & volcano_df$logFC < -lfc_cut
] <- "Up in saline"

p_volcano <- ggplot(
  volcano_df,
  aes(x = logFC, y = neglog10p, color = significance)
) +
  geom_point(alpha = 0.6, size = 1) +
  theme_bw(base_size = 14) +
  labs(x = "log2 Fold Change (cocaine vs saline)", y = "-log10(P-value)")

ggsave(
  filename = "output/figures/volcano_plot.png",
  plot = p_volcano,
  width = 7,
  height = 5,
  dpi = 300
)

# MA plot (limma)
png("output/figures/ma_plot.png", width = 1200, height = 900, res = 150)
plotMA(fit, coef = "groupcocaine")
dev.off()


# Heatmap (top 50 most significant genes)

top_n <- 50
top_ids <- rownames(de_results)[order(de_results$adj.P.Val)][1:top_n]

exprs_heatmap <- vGene$E[top_ids, ]

# Annotate columns by group
annotation_col2 <- data.frame(group = meta$group)
row.names(annotation_col2) <- colnames(exprs_heatmap)

pheatmap(
  exprs_heatmap,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = annotation_col2,
  filename = "output/figures/heatmap_top50.png",
  width = 7,
  height = 9
)

# Save top tables for report writing
top20 <- de_results[order(de_results$adj.P.Val), ][1:20, ]
write.csv(top20, "output/tables/top20_genes.csv", row.names = FALSE)
