##Volcano plot of DEGs
filename <- "exp-profile/mtx/qc/DEG_analysis/myeloid_vs_lymphoid/results/Myeloid_vs_Lymphoid.degs.txt"

res <- read.table(filename, sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
head(res)

library(org.Hs.eg.db)
res.df <- as.data.frame(res)
#res.df$symbols <- mapIds(org.Hs.eg.db, keys = rownames(res.df), keytype = 'ENSEMBL', column = "SYMBOL")
head(res.df)

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages('BiocManager')

#BiocManager::install('EnhancedVolcano')


library(EnhancedVolcano)

res.sub <- subset(res.df, res.df$type == 'Known_lncRNA' | res.df$type == 'Novel_lncRNA')
head(res.sub)

select <- c('KIAA0513','MCOLN1', 'ETS2', 'NAMPT', 'SAT1', 'SLC9B2', 'SIT1', 'MT-ND1', 'MT-ND3', 'NEFL','STAMBPL1' )
EnhancedVolcano(res.sub, x = "log2FoldChange", y = "padj", lab = res.sub$symbols,
                pCutoff = 1e-4, FCcutoff = 2)

ggsave(filename = 'exp-profile/mtx/qc/DEG_analysis/myeloid_vs_lymphoid/results/Myel_vs_Lymph.lncRNA-deg.pdf', width = 17, height = 15, units = "cm")
