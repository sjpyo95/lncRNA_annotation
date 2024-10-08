###Count expressing genes in each cell-types
library(ggplot2)
tpm_mtx <- read.table('exp-profile/mtx/qc/info/cv/10TPM_count.txt', sep = '\t', header = TRUE, check.names = FALSE)
head(tpm_mtx)

celltypes <- unique(tpm_mtx$sample)
celltypes
tpm_mtx <- tpm_mtx[grep('lncRNA', tpm_mtx$type),]
tpm_mtx <- subset(tpm_mtx, genetype != 'lncRNA')
tpm_mtx$sample <- factor(tpm_mtx$sample, levels = celltypes)
tpm_mtx$type <- factor(tpm_mtx$type, levels = c('Novel_lncRNA', 'Known_lncRNA'))
tpm_mtx$genetype <- factor(tpm_mtx$genetype, levels = c('lincRNA','antisense', 'intervening', 'flanking_region'))
tpm_mtx
ggplot(tpm_mtx, aes(x = sample, y = count, color = type, fill = type, alpha = genetype)) + geom_bar(stat = 'identity')+theme_classic()+theme(legend.title = element_blank(),axis.text = element_text(colour = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.0,), legend.key.size = unit(0.3, 'cm'), axis.title = element_text(size = 10), legend.text = element_text(size = 10), plot.title = element_text(hjust = 0.5))+ xlab('Immune cell-types') + ylab('Count') + ggtitle('33 Immune cell-types highly expressing lncRNAs count (10>=TPM)')

ggsave(filename = 'exp-profile/mtx/qc/stat/graph/10TPM_hi_expressing_lncRNA_count.pdf', width = 20, height = 17, units = "cm")


mtx <- read.table('dataset/expression/mtxs/batchQC/all/imm.qc.info.celltypeAvg.tpm.genetype_edit.txt', sep = '\t', header = TRUE, check.names = FALSE)
