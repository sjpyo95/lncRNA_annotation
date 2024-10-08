if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)
library(cluster)
install.packages('seriation')
library(seriation)
#update.packages()
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

save_pheatmap_pdf <- function(x, filename, width=20, height=20) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

col_fun = colorRamp2(c(0, 3, 6, 12, 15), c("cornflowerblue", 'mediumspringgreen', 'yellow','orangered',"indianred4"))
col_fun(seq(-3, 3)) 

table <- read.table('exp-profile/mtx/qc/deg/group1/lncRNAs/deg.lncRNAs.top10s.mtx.txt', sep = '\t', header = TRUE, check.names = FALSE)
head(table)
dim(table)


#table <- subset(table, table$type != 'PCG')

exptables <- table[,c(-1,-2,-3,-4,-5,-6)]
head(exptables)
exptables <- log2(exptables+1)
rownames(exptables) <- table$symbol
head(exptables)
summary(exptables)

heat <- Heatmap(exptables, name = 'log2 TPM+1', col = col_fun, column_order = colnames(exptables))
heat
save_pheatmap_pdf(heat, 'results/figure3/tcell.degs_lncRNAs.heatmap.pdf')
