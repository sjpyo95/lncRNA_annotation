install.packages("devtools")
devtools::install_github("zhangyuqing/sva-devel")
library('sva')

count_matrix <- read.table('exp-profile/mtx/raw/imm.RC.txt', sep = '\t', header = TRUE, check.names = FALSE)
meta <- read.table('exp-profile/groupinfo/imm33.groups.txt', sep = '\t', header = TRUE, check.names = FALSE)
dim(count_matrix)
head(count_matrix)
dim(meta)
head(meta)
batch <- meta$batch
head(batch)
all(meta$samples %in% colnames(count_matrix))
all(meta$samples == colnames(count_matrix))
group <- meta$group1

group
count_matrix$sum <- rowSums(count_matrix[,c(-1,-2)])
filt.matrix <- subset(count_matrix, sum != 0)
dim(filt.matrix)
filt.matrix$sum <- NULL
head(filt.matrix)
adjusted_counts <- ComBat_seq(filt.matrix[,c(-1,-2)], batch=batch)
head(adjusted_counts)
dim(adjusted_counts)

adjusted_counts <- data.frame(adjusted_counts, check.names = FALSE)
adjusted_counts <- cbind(ID = filt.matrix$ID, symbol = filt.matrix$symbol, adjusted_counts)
row.names(adjusted_counts) <- NULL
head(adjusted_counts)
write.table(adjusted_counts, file = 'exp-profile/mtx/qc/rc/imm.QC.RC.txt', append = FALSE, sep = '\t', eol = '\n', na = 'NA', row.names = FALSE, col.names = TRUE, quote = FALSE)
