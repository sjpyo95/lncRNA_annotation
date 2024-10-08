##GSEA analysis
#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db") #org.Mm.eg.db for mouse

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("clusterProfiler")

#if (!requireNamespace('BiocManager', quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("AnnotationDbi")

#install.packages("https://cran.r-project.org/src/contrib/rlang_1.0.6.tar.gz", repo=NULL,type='source')
library(clusterProfiler)
library(org.Hs.eg.db)

library(AnnotationDbi)
library(DOSE)
library(enrichplot)
deg <- read.table(file = 'exp-profile/mtx/qc/DEG_analysis/myeloid/DEG_pcg_correlation.mtx.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
dim(deg)
deg <- unique(deg)
rownames(deg) <- deg$symbol
res <- deg[,c(-1:-3)]
dim(res)

d.cor <- cor(t(res), method = 'pearson')

dim(d.cor)
d.cor <- as.data.frame(d.cor)
head(d.cor[,1])
res <- subset(d.cor, abs(d.cor$`GK-AS1`) >= 0.7)[1]
dim(res)
org <- res$`GK-AS1`
names(org) <- rownames(res)

genelist <- na.omit(org)

#res <- t(res)
genelist <- sort(genelist, decreasing = TRUE)
#myel <- deg[deg$log2FoldChange >= 2 & deg$padj < 0.05 & deg$type == 'PCG' & deg$baseMean > 50,]
#lymph <- deg[deg$log2FoldChange <= -2 & deg$padj < 0.05 & deg$type == 'PCG' & deg$baseMean > 50,]

head(genelist)

gse <- gseGO(geneList = genelist,
             ont = 'BP',
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             pvalueCutoff = 0.05)
p <- cnetplot(gse, categorySize = "pvalue", foldChange = genelist, showCategory = 5)
p
pdf(file = "exp-profile/mtx/qc/DEG_analysis/myeloid/GK-AS1.pcg_correlation.pdf", width = 25, height = 15)
p
dev.off()
library(ggridges)
library(ggplot2)
library(viridis)
library(hrbrthemes)
gse <- as.data.frame(gse)
head(summary(gse))
ggplot(gse, aes(x = , y = Description, fill = cut)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none")
