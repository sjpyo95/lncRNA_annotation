##Gene expression correlation
library(dplyr)
library(corrplot)
library(ComplexHeatmap)
library(reshape)
library('magick')
library(circlize)
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "cornsilk", "red"))
col_fun(seq(-3, 3))

filename <- 'dataset/expression/mtxs/deg_results/all/correlation/deg_hallmarks/mtx/Plasmablasts_deg_lncRNAs.hallmarks.mtx.txt'
mtx <- read.table(file = filename, sep = '\t', header = TRUE, check.names = FALSE)
grtmp <- strsplit(filename,'/')[[1]]
grtmp2 <- strsplit(grtmp[length(grtmp)],'_deg')[[1]]
gr <- grtmp2[1]
gr

head(mtx)
dim(mtx)
mtx <- distinct(mtx)

mtx <- subset(mtx, abs(logFC) >= 3)
mtx <- mtx[order(mtx$pathway),]
datExpr <- mtx[,-c(1:5)]
rownames(datExpr) <- mtx$ID
summary(datExpr)
datExpr <- datExpr[apply(datExpr>=1,1,any),]
lncs <- subset(mtx, pathway == 'DEG')$ID
pcgs <- subset(mtx, pathway != 'DEG')$ID
length(lncs)
length(pcgs)

head(datExpr)
dim(datExpr)
#anno <- data.frame(ID = rownames(datExpr), pathway = datExpr$pathway)
#head(anno)
#dim(anno)
#lanno <- list()
#str(lanno)

#for (i in c(1:nrow(anno))) {
#  pathway <- anno$pathway[i]
#  id <- anno$ID[i]
#  if (!(pathway %in% names(lanno))) {
#    lanno[[pathway]] <- id
#  } else {
#    lanno[[pathway]] <- c(lanno[[pathway]], id)
#  }
#}

#str(lanno)


cormat <- round(cor(t(datExpr)),2)
#head(cormat)
cormat <- reorder_cormat(cormat)
dim(cormat)
cormat <- as.data.frame(cormat)
#colnames(cormat)
lncvspcg <- cormat[,(rownames(cormat) %in% lncs)]
lncvspcg <- lncvspcg[(colnames(cormat) %in% pcgs),]
#head(lncvspcg)
dim(lncvspcg)

pdf(file = paste("results/correlation_plots/", gr, '.degs_lnc_vs_hallmarks.cor.heatmap.pdf'), width = 10, height = 8)

#ha = rowAnnotation(foo = anno_empty(border = FALSE))
heat <- Heatmap(lncvspcg, name = 'Preason cor.', column_title = paste(gr, ' DEG lncRNAs vs hallmarks correlation', sep = ''), col = col_fun, use_raster = TRUE)

heat
draw(heat)

dev.off()

write.table(lncvspcg, file = paste('dataset/expression/mtxs/deg_results/all/correlation/deg_hallmarks/cor_table/', gr,'.degs_lnc_vs_hallmark.cor.txt', sep = ''), append = FALSE, sep = '\t', eol = '\n', na = 'NA', row.names = TRUE, col.names = TRUE, quote = FALSE)

