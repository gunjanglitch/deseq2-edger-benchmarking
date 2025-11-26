#loading libraries
library(edgeR)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)

#setting color palette
rmain = "#e31a1c"
rlight = "#fc9272"
rdark = "#99000d"

#setting working directory
setwd("C:/Users/Asus/Desktop/benchmarking/")

#load data
count_edger = read.csv("brownbear_counts.csv", row.names = 1, check.names = FALSE)
meta_edger = read.csv("brownbear_metadata.csv", row.names = 1)

#create DGEList object
group = factor(meta_edger$condition)
dge = DGEList(counts = count_edger, group = group)

#set factor level
dge$samples$group = relevel(dge$samples$group, ref = "Active")

#filter low expressed genes (CPM based)
keep = rowSums(cpm(dge) > 1) >=2
dge = dge[keep, ,keep.lib.sizes = FALSE]

#calculate library size
dge =  calcNormFactors(dge)
write.csv(dge$samples, "edgeR_librarysize_normcounts.csv")

#TMM normalization
dge = calcNormFactors(dge, method = "TMM")

#estimate dispersion
design =  model.matrix(~ group)
dge = estimateDisp(dge, design)
write.csv(dge$common.dispersion, "edgeR_common_dispersion.csv")
write.csv(dge$trended.dispersion, "edgeR_trended_dispersion.csv")
write.csv(dge$tagwise.dispersion, "edgeR_tagwise_dispersion.csv")

#build GLM fit
fit = glmQLFit(dge, design)

#quasi-likelihood F-test
qlf = glmQLFTest(fit, coef = 2)
res = topTags(qlf, n = Inf)$table
write.csv(res, "edgeR_quasi_result.csv")

#compute logCPM matrix
logCPM = cpm(dge, log = TRUE)
write.csv(logCPM, "edgeR_logCPM_matrix.csv")

#MDS Plot 
png("edgeR_MDS_plot.png", width = 900, height = 700)
plotMDS(dge, col=ifelse(group == "Active", rlight, rdark), main = "MDS Plot (edgeR LOGCPM)")
dev.off()

#sample dsiatnce visualization
dist_mat = dist(t(logCPM))
dist_mat_matrix = as.matrix(dist_mat)
pheatmap(dist_mat_matrix, color = colorRampPalette(c(rlight, rmain, rdark))(100), main = "Sample Disatnce Matrix (logCPM)", filename = "edgeR_sample_distance_matrix.png")

#BCV Plot
png("edgeR_BCV_plot.png", width = 900, height = 700)
plotBCV(dge,
        xlab = "Average logCPM",
        ylab = "BCV",
        pch = 16,
        cex = 0.7,
        col.tagwise = rmain,
        col.trend = "red",
        col.common = "blue",
        main = "edgeR_BCV_Plot ")
dev.off()

#MA Plot
png("edgeR_MA_plot.png", width = 900, height = 700)
plotMD(qlf, col=c("grey70", rmain), main = "MA Plot (edgeR)")
dev.off()

#volcano plot
EnhancedVolcano(res, lab = rownames(res), x="logFC", y="PValue", title = "Volcano Plot (edgeR)", pCutoff = 0.05, FCcutoff = 1.0, col = c("grey50", rlight, rmain, rdark))
ggsave("edgeR_volcano_plot.png", width = 9, height = 7)

#heatmap of top 20 DEGs
top_genes = rownames(res)[1:20]
heat_data = logCPM[top_genes, ]
pheatmap(heat_data, scale = "row", color = colorRampPalette(c(rlight, rmain, rdark))(100),main = "Top 20 DEGs (edgeR)", filename = "edgeR_heatmap_top20.png")

#DEG Counts
deg_count = sum(res$FDR < 0.05)
up = sum(res$logFC > 0 & res$FDR < 0.05)
down = sum(res$logFC < 0 & res$FDR < 0.05)
write.csv(data.frame(total_DEGs = deg_count, upregulated = up, downregulated = down), "edgeR_DEG_counts.csv")

#saving significant/up/down genes
res_df = as.data.frame(res)
sig_df = res_df[which(res_df$FDR < 0.05), ]
write.csv(sig_df, "edgeR_significant_genes.csv")

up_df = res_df[which(res_df$logFC > 0 & res_df$FDR < 0.05), ]
write.csv(up_df, "edgeR_upregulates_genes.csv")

down_df =  res_df[which(res_df$logFC < 0 & res$FDR < 0.05), ]
write.csv(down_df, "edgeR_downregulated_genes.csv")

#summary table
write.csv(res, "edgeR_final_result_table.csv")
