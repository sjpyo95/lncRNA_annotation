## Cell-type-specific expression patterns (boxplot)

library(ggplot2)

degs <- read.table('dataset/expression/mtxs/batchQC/degs.all.1tpmhigher.txt', sep = '\t', header = TRUE, check.names = FALSE)
head(degs)
summary(degs)
degs$expression <- log10(degs$expression)

degs$type <- factor(degs$type, levels = c('PCG','Known_lncRNA', 'Novel_lncRNA'))
ggplot(degs, aes(x = degs$`cell-type`, y= expression, fill = type))+geom_boxplot()+
  theme_classic()+ theme(legend.title = element_blank(),
                         axis.text.x = element_text(colour = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.0,),
                         legend.key.size = unit(0.3, 'cm'),
                         axis.title = element_text(size = 15),
                         legend.text = element_text(size = 10))+
  xlab('Cell-types') + ylab('log10(TPM)') + ggtitle('Human Immune cells expression (TPM ≥ 1)') + theme(plot.title = element_text(hjust = 0.5))
