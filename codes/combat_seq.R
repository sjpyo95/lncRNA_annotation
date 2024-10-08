# install.packages("devtools")
#devtools::install_github("zhangyuqing/sva-devel")
library('sva')
count_matrix <- read.table('exp-profile/mtx/qc/rc/imm.QC.srt.RC.txt', sep = '\t', header = TRUE, check.names = FALSE)

coldata <- read.table('exp-profile/groupinfo/imm128samples.groups.txt', sep = '\t', header = TRUE, check.names = FALSE)
batch <- coldata$batch
batch
group <- coldata$group1
head(group)
cov1 <- rep(c(0,1), 4)
cov2 <- c(0,0,1,1,0,0,1,1)
covar_mat <- cbind(cov1, cov2)
covar_mat
head(count_matrix)
adjusted <- ComBat_seq(count_matrix, batch=batch, group=group)
head(adjusted)
adjusted <- data.frame(adjusted, check.names = FALSE)
adjusted <- cbind(ID = count_matrix$ID, symbol = count_matrix$symbol, adjusted)
row.names(adjusted_counts) <- NULL
head(adjusted_counts)
write.table(adjusted, file = 'exp-profile/mtx/qc/rc/imm.batchQC.MR.txt', append = FALSE, sep = '\t', eol = '\n', na = 'NA', row.names = FALSE, col.names = TRUE, quote = FALSE)
