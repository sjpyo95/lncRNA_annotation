#Clustering based on lncRNA expression across cell-types and find signatures
library(ggplot2)
library("ape")
tpm_mtx <- read.table('dataset/expression/mtxs/batchQC/all/imm.qc.info.celltypeAvg.tpm.genetype_edit.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)

tpm_mtx <- subset(tpm_mtx, type != 'PCG')
head(tpm_mtx)
dim(tpm_mtx)
mtx <- tpm_mtx[,c(-1,-2,-3)]

srt_mtx <- subset(mtx, apply(mtx,1,function(x) any(x >= 5)))#1TPM 이상이 하나 이상의 cell-type에서 나타나는 gene만 가져오기
dim(srt_mtx)

dd <- dist(t(srt_mtx), method = "euclidean")

hc <- hclust(dd, method = 'ward.D2')
hcd <- as.dendrogram(hc)

plot(hcd, type = 'rectangle', ylab = 'Expression')
