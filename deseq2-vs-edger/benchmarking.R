library(dplyr)
library(ggplot2)
library(VennDiagram)
library(pheatmap)
setwd("C:/Users/Asus/Desktop/benchmarking/")
deseq = read.csv("deseq2_fullresult.csv", row.names = 1, check.names = FALSE)
edger = read.csv("edgeR_final_result_table.csv", row.names = 1, check.names = FALSE)
rlog = read.csv("deseq2_rlog_matrix.csv", row.names = 1, check.names = FALSE)
logcpm = read.csv("edgeR_logCPM_matrix.csv", row.names = 1, check.names = FALSE)
#total deg
n_deseq = sum(deseq$padj < 0.05, na.rm = TRUE)
n_deseq
n_edger = sum(edger$FDR < 0.05, na.rm = TRUE)
n_edger
#upregulated/downregulated deseq
up_deseq = sum(deseq$log2FoldChange > 0 & deseq$padj < 0.05, na.rm = TRUE)
up_deseq
down_deseq = sum(deseq$log2FoldChange < 0 & deseq$padj < 0.05, na.rm = TRUE)
down_deseq
#upregulated/downregulated edger
up_edger = sum(edger$logFC > 0 & edger$FDR < 0.05)
up_edger
down_edger = sum(edger$logFC < 0 & edger$FDR < 0.05)
down_edger
#overlap deg
#total deg in deseq
deseq_deg = rownames(deseq[deseq$padj < 0.05, ])
deseq_deg
length(deseq_deg)
#total deg in edger
edger_deg = rownames(edger[edger$FDR < 0.05, ])
edger_deg
length(edger_deg)
#common in deseq & edger
common_deg = intersect(deseq_deg, edger_deg)
common_deg
length(common_deg)
#unique gene in deseq
unique_deseq = setdiff(deseq_deg, edger_deg)
length(unique_deseq)
#unique gene in edger
unique_edger = setdiff(edger_deg, deseq_deg)
length(unique_edger)
#--------------intermediate steps-----------------
length(rownames(deseq))
length(rownames(edger))
deseq$gene = rownames(deseq)
edger$gene = rownames(edger)
merged_fc = merge(
  deseq[, c("gene", "log2FoldChange")],
  edger[, c("gene", "logFC")],
  by = "gene"
)
#log2fc correlation
fc_cor = cor(merged_fc$log2FoldChange, merged_fc$logFC, use = "pairwise.complete.obs")
fc_cor

#log2fc correlation plot
ggplot(merged_fc, aes(x = log2FoldChange, y = logFC)) + geom_point(alpha = 0.4) + geom_smooth(method = "lm", se = FALSE, color = "red") + labs(title = "Log2 Fold Change Correlation (DESeq2 vs edgeR)", x = "DESeq2 log2FoldChange", y="edgeR lofFC") + theme_bw()

ggsave("compare_logfc_cor_plot.png", width=2000, height=1500, units = "px")

#adjusted p-value correlation
merged_p = merge(deseq[, c("gene", "padj")],
                 edger[, c("gene", "FDR")],
                 by = "gene")
#correlation with spearman
p_cor_spearman = cor(merged_p$padj, merged_p$FDR, method = "spearman", use = "pairwise.complete.obs")
p_cor_spearman
#correlation with kendall
p_cor_kendall = cor(merged_p$padj, merged_p$FDR, method = "kendall", use = "pairwise.complete.obs")
p_cor_kendall
#correlation with pearson
p_cor_pearson = cor(merged_p$padj, merged_p$FDR, method = "pearson", use = "pairwise.complete.obs")
p_cor_pearson

ggplot(merged_p, aes(x = -log10(padj), y =-log10(FDR))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(
    title = "Adjusted p-value Correlation (DESeq2 vs edgeR)",
    x="-log10(DESeq2 padj)",
    y="-log10(edgeR FDR)"
  ) +
  theme_bw()

ggsave("compare_padj_cor_spearman.png", width = 6, height = 5)

#Sign ConCordance
sign_deseq = ifelse(deseq$log2FoldChange > 0, "up", "down")
sign_deseq
sign_edger = ifelse(edger$logFC > 0, "up", "down")
sign_concordance = mean(sign_deseq == sign_edger) * 100
sign_concordance

write.csv(
  data.frame(sign_concordance_percent = sign_concordance),
  "sign_concordance.csv",
  row.names = FALSE
)
png("compare_sign_concordance.png", width = 800, height = 600)
barplot(
  sign_concordance,
  names.arg = "Sign Concordance (%)",
  col = "skyblue",
  main = "Sign Concordance between DESeq2 and edgeR",
  ylab = "Percentage (%)",
  ylim = c(0, 100)
)
dev.off()

df <- data.frame(metric = "Sign Concordance (%)",
                 value = sign_concordance)
ggplot(df, aes(x = metric, y = value)) +
  geom_col(fill = "skyblue") +
  ylim(0, 100) +
  labs(
    title = "Sign Concordance between DESeq2 and edgeR",
    y = "Percentage (%)",
    x = ""
  ) +
  theme_bw()

#Top 100 Overlap
top100_deseq = rownames(deseq[order(deseq$padj), ])[1:100]
top100_edger = rownames(edger[order(edger$FDR), ])[1:100]
top100_common = intersect(top100_deseq, top100_edger)
top100_common
length(top100_common)
top100_deseq
length(top100_deseq)
top100_edger
length(top100_edger)

df <- data.frame(
  gene = c(top100_deseq, top100_edger, top100_common),
  source = c(
    rep("DESeq2", length(top100_deseq)),
    rep("edgeR", length(top100_edger)),
    rep("Common", length(top100_common))
  )
)
write.csv(df, "top_100_overlap.csv", row.names = FALSE)


#Normalized counts correlation (rlog vs logCPM)
common_genes = intersect(rownames(rlog), rownames(logcpm))
common_genes
norm_cor = cor(rlog[common_genes, 1], logcpm[common_genes, 1])
norm_cor

#dispersion trend proxy (variance correlation)

var_rlog = apply(rlog[common_genes, ], 1, var)
var_cpm = apply(logcpm[common_genes, ], 1, var)
disp_cor = cor(var_rlog, var_cpm)
disp_cor

sessionInfo()
