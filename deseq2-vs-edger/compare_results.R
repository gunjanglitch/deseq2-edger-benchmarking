#loading libraries
library(tidyverse)
library(pheatmap)
library(ggpubr)
library(cowplot)
library(grid)
library(gridExtra)
library(reshape2)
library(EnhancedVolcano)
library(VennDiagram)
library(UpSetR)
library(pROC)
library(PRROC)
library(DESeq2)
library(edgeR)
library(ggVennDiagram)
library(ggplot2)
#installing packages
install.packages("VennDiagram")
install.packages("UpSetR")
install.packages("pROC")
install.packages("PRROC")
install.packages("ggVennDiagram")

#setting working directory
setwd("C:/Users/Asus/Desktop/benchmarking/")

#setting files paths
deseq2_fp = file.path("deseq2_fullresult.csv")
deseq2_shrink =  file.path("deseq2_lfc_shrunken.csv")
deseq2_rlog = file.path("deseq2_rlog_matrix.csv")

edgeR_fp = file.path("edgeR_final_result_table.csv")
edgeR_logcpm = file.path("edgeR_logCPM_matrix.csv")

deseq_time_fp
edger_time_fp

#color palette
cols = list(
  deseq2 = "#1f78b4",
  edger = "#e31a1c",
  both = "#33a02c",
  neutral = "#b2b2b2"
)

#load files
deseq = read.csv(deseq2_fp, row.names = 1)
deseq_shr = read.csv(deseq2_shrink, row.names = 1)
deseq_rlog = read.csv(deseq2_rlog, row.names = 1)
edger = read.csv(edgeR_fp, row.names = 1)
edger_logcpm = read.csv(edgeR_logcpm, row.names = 1)

deseq$gene = rownames(deseq)
edger$gene = rownames(edger)
head(deseq$gene)
head(edger$gene)
#harmonize gene names and merge tables
#add gene coloum
add_gene = function(df){
  if (is.null(df)) return(NULL)
  df$gene = rownames(df)
  df
}

deseq = add_gene(deseq)
edger = add_gene(edger)

#standardised names
names(deseq) = make.names(names(deseq))
names(edger) = make.names(names(edger))

#select useful columns
deseq_sel =  deseq %>% select(gene, log2FoldChange, pvalue, padj, everything())
edger_sel = edger %>% select(gene, logFC, PValue, FDR, everything())

#merge by gene
merged = full_join(
  deseq_sel %>% rename_with(~ paste0("deseq_", .)),
  edger_sel %>% rename_with(~ paste0("edger_", .)),
  by =c("deseq_gene" = "edger_gene")
)

#fix gene column
merged$gene = merged$deseq_gene
merged$gene[is.na(merged$gene)] <- merged$edger_gene[is.na(merged$gene)]
merged <- merged %>% select(gene, everything())

#DEG Sets and counts
#deseq_deg = deseq %>% filter(padj < 0.05) %>% pull(gene)
dseq_deg <- rownames(deseq)[deseq$padj < 0.05 & abs(deseq$log2FoldChange) > 1]
head(dseq_deg)
edger_deg = edger %>% filter(FDR < 0.05) %>% pull(gene)
head(edger_deg)
summary_deg = data.frame(
  method = c("DESeq2", "edgeR"),
  total = c(length(dseq_deg), length(edger_deg))
)

write.csv(summary_deg, "deseq_edger_compare_deg_counts.csv")

#overlap, jaccard, venn, upset
common = intersect(dseq_deg, edger_deg)
head(common)
union_all = union(dseq_deg, edger_deg)

jaccard = length(common)/length(union_all)

#png("clean_compare_venn.png", 600, 600)
#draw.pairwise.venn(
  #area1 = length(deseq_deg),
  #area2 = length(edger_deg),
  #cross.area = length(common),
  #category = c("DESeq2", "edgeR"),
  #fill = c("#4C72B0", "#DD8452"),
  #alpha = c(0.6, 0.6),
  #lty = "solid",
  #cat.cex = 3,
  #cex = 3,
  #cat.pos = c(-20, 20),
  #cat.dist = c(0.08, 0.08),
  #lwd = 2,
  #main = "Overlap of DEG: DESeq2 vs edgeR", main.cex = 3
#)
#dev.off()



   
#upset
upset_list = list(DESeq2 = dseq_deg, edgeR = edger_deg)
png("compare_upset.pdf")
fromList(upset_list) %>% upset()
dev.off()

#fold-change scatter + correlations
compare = merged %>%
  filter(!is.na(deseq_log2FoldChange), !is.na(edger_logFC)) %>% select(gene, deseq_log2FoldChange, edger_logFC)

#scatter plot
ggplot(compare, aes(deseq_log2FoldChange, edger_logFC)) + geom_point(alpha = 0.5) + geom_smooth(method = "lm") + theme_minimal() + labs(x = "DESeq2 Log2FC", y ="edgeR LogFC") + ggtitle("Log2FC Comparision") 
ggsave("compare_scatter_logFC.png", width = 6, height = 6)

#correlations
cor_df = data.frame(
  pearson = cor(compare$deseq_log2FoldChange, compare$edger_logFC, method = "pearson"),
  spearman = cor(compare$deseq_log2FoldChange, compare$edger_logFC, method = "spearman")
)
write.csv(cor_df, "compare_fc_correlation.csv")

#volcano - DESeq2
png("deseq2_compare_volcano.png", 900, 800)
EnhancedVolcano(deseq, lab = deseq$gene, x="log2FoldChange", y="pvalue", pCutoff=0.05, FCcutoff=1)
dev.off()

#volcano - edgeR
png("edger_compare_volcano.png", 900, 800)
EnhancedVolcano(edger, lab = edger$gene, x="logFC", y="PValue", pCutoff = 0.05, FCcutoff = 1)
dev.off()

#MA Plot - DESeq2
png("deseq2_compare_ma.png", 800, 700)
plot(log10(deseq$baseMean+1), deseq$log2FoldChange, pch=20, cex=0.5, xlab="log10 baseMean", ylab="log2FC")
dev.off()

#MA Plot - edgeR
head(rownames(logCPM))
head(edger$gene)
sum(edger$gene %in% rownames(logCPM))
edger = read.csv("edgeR_final_result_table.csv", row.names = 1)
head(rownames(edger))
edger$gene = rownames(edger)
mean_lcpm <- rowMeans(logCPM, na.rm = TRUE)
edger$mean_lcpm <- mean_lcpm[edger$gene]
ed_clean <- edger %>% filter(!is.na(mean_lcpm))
png("edger_compare_MA.png", 800, 700)
plot(ed_clean$mean_lcpm,
     ed_clean$logFC,
     pch = 20, cex = 0.5, col = cols$edger,
     xlab = "mean logCPM",
     ylab = "logFC",
     main = "edgeR MA Plot")
dev.off()

#method specific DE Genes
deseq_only = setdiff(dseq_deg, edger_deg)
head(deseq_only)
edger_only = setdiff(edger_deg, dseq_deg)
head(edger_only)
write.csv(deseq_only, "compare_deseq_only.csv", row.names = FALSE)
write.csv(edger_only, "compare_edger_only.csv", row.names = FALSE)
write.csv(common, "compare_common_deg.csv", row.names = FALSE)

#PCA & MDS + Heatmap
meta <- read.csv("brownbear_metadata.csv")
rownames(meta) <- meta$sample
meta
list.files()
counts = read.csv()
# rlog transform
rld <- rlog(dds)

# PCA using DESeq2 built-in function
pcaData <- plotPCA(rld, intgroup = c("condition", "individual"), returnData = TRUE)

# % variance
percentVar <- round(100 * attr(pcaData, "percentVar"))

# PCA plot
p <- ggplot(pcaData, aes(PC1, PC2, color = condition, shape = individual, label = name)) +
  geom_point(size = 4) +
  geom_text(vjust = -1, size = 3) +
  theme_minimal(base_size = 14) +
  xlab(paste0("PC1 (", percentVar[1], "%)")) +
  ylab(paste0("PC2 (", percentVar[2], "%)")) +
  ggtitle("PCA of Brown Bear RNA-seq (rlog)") +
  theme(plot.title = element_text(hjust=0.5))

# Save
ggsave("compare_brownbear_pca_rlog.png", plot = p, width = 6, height = 5)
ls()
head(countdata)
colnames(countdata)
rm(list = ls())
ls()
