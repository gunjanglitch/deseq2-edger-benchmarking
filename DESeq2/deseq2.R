#loading libraries
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)
library(apeglm)

#setting color palette 
dblue = "#1f78b4"
llight = "#6baed6"
ddark = "#08306b"

#setting the working directory
setwd("C:/Users/Asus/Desktop/benchmarking/")

#load data
count_gene = read.csv("brownbear_counts.csv", row.names = 1, check.names = FALSE )
meta_gene = read.csv("brownbear_metadata.csv", row.names = 1)

head(count_gene)
head(meta_gene)
dim(count_gene)
dim(meta_gene)
counts = count_gene[, rownames(meta_gene)]

#creating DESeq2 dataset
differexp=DESeqDataSetFromMatrix(countData = round(count_gene), colData = meta_gene, design = ~condition)

#factor levels
differexp$condition = relevel(differexp$condition, ref = "Active")

#filtering low count gene
keptgene = rowSums(counts(differexp) >= 10) >= 2
differexp = differexp[keptgene, ]

#library size inspection
write.csv(colSums(counts(differexp)), "deseq2_librarysizes.csv")

#normalization
differexp = estimateSizeFactors(differexp)
norm_count = counts(differexp, normalized=TRUE)
write.csv(norm_count, "deseq2_normalized_counts.csv")

#dispersion estimation and GLM fitting
differexp = estimateDispersions(differexp)
differexp = nbinomWaldTest(differexp)

# === PROOF: Dispersion shrinkage for intermediate-count genes ===
png("deseq2_dispersion_trend.png", 900, 700)
plotDispEsts(differexp, main="DESeq2 Dispersion Estimation")
dev.off()

mean_counts <- rowMeans(counts(differexp))
deseq_disp <- mcols(differexp)$dispGeneEst
edgeR_tag <- read.csv("edgeR_tagwise_dispersion.csv", row.names=1)[,1]

length(names(deseq_disp))
length(mean_counts)
length(deseq_disp)
length(edgeR_tag[names(deseq_disp)])

df_disp <- data.frame(
  gene = names(deseq_disp),
  mean = mean_counts,
  deseq = deseq_disp,
  edger = edgeR_tag[names(deseq_disp)]
)

# intermediate-count genes (the “edge” zone)
df_mid <- df_disp[df_disp$mean > 10 & df_disp$mean < 1000, ]
write.csv(df_mid, "deseq2_vs_edgeR_midgene_dispersion.csv")

#==============================================================

#differentiation expression analysis
signi_genes = results(differexp)
signigenes_ordered = signi_genes[order(signi_genes$padj), ]
write.csv(as.data.frame(signigenes_ordered), "deseq2_fullresult.csv")

#log2 fold change shinkage
resLFC = lfcShrink(differexp, coef = "condition_Hibernating_vs_Active", type = "apeglm")
write.csv(as.data.frame(resLFC), "deseq2_lfc_shrunken.csv")

#rlog transformation (for PCA & Heatmaps)
rld = rlog(differexp, blind = TRUE)
write.csv(assay(rld), "dese2_rlog_matrix.csv")

#sample distance heatmap
sampledist = dist(t(assay(rld)))
sampledistmatrix = as.matrix(sampledist)
pheatmap(sampledistmatrix, clustering_distance_rows = sampledist, clustering_distance_cols = sampledist, color = colorRampPalette(c(llight, dblue, ddark))(100), main = "Sample Distance Matrix (rlog)", filename = "deseq2_sample_distance_)matrix.png")

#PCA plot 
pca = plotPCA(rld, intgroup = "condition")

pca_custom = pca + scale_color_manual(values = c("Active" = llight, "Hibernating" = ddark)) + ggtitle("PCA plot (rlog, DESeq2)") + theme_minimal(base_size = 14)
pca_custom
ggsave("deseq2_pca_plot.png", pca_custom, width = 7, height = 6)

#MA Plot
png("deseq2_ma_plot.png", width = 900, height = 700)
plotMA(signi_genes, ylim =c(-5, 5), colSig = dblue, alpha = 0.05, main = "MA Plot (DESeq2)")
dev.off()

#volcano plot
EnhancedVolcano(signi_genes, lab = rownames(signi_genes), x = "log2FoldChange", y = "pvalue", title = "Volcano Plot (DESeq2)", pCutoff = 0.05, FCcutoff = 1, col = c("grey50", llight, dblue, ddark))
ggsave("deseq2_volcano_plot.png", width = 9, height = 7)

#Heatmap of Top 20 DEGs
top_genes = head(rownames(signigenes_ordered), 20)
head(top_genes)
heat_data = assay(rld)[top_genes, ]
pheatmap(heat_data, color = colorRampPalette(c(llight, dblue, ddark))(100), scale = "row", main = "Top 20 Differentially Expressed Genes", filename = "deseq2_heatmap_top20.png")

#DEG Counts
deg_counts = sum(signi_genes$padj < 0.05, na.rm = TRUE)
up = sum(signi_genes$log2FoldChange > 0 & signi_genes$padj < 0.05, na.rm = TRUE)
down = sum(signi_genes$log2FoldChange < 0 & signi_genes$padj < 0.05, na.rm = TRUE)
write.csv(data.frame(total_DEGs = deg_counts, upregulated = up, downregulated = down), "deseq2_DEG_counts.csv")

#saving significant/up/down genes
res_df = as.data.frame(signi_genes)

sig_df = res_df[which(res_df$padj < 0.05), ]
write.csv(sig_df, "deseq2_significant_genes.csv")

up_df = res_df[which(res_df$log2FoldChange > 0 & res_df$padj < 0.05), ]
write.csv(up_df, "deseq2_upregulated_genes.csv")

down_df = res_df[which(res_df$log2FoldChange < 0 & res_df$padj < 0.05), ]
write.csv(down_df, "deseq2_downregulated.csv")

#final output
write.csv(as.data.frame(signigenes_ordered), "deseq2_final_results_table.csv")

sessionInfo()
