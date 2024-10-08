###lncRNA genomic context 별 stats
library(ggplot2)
library(scales)

count <- read.table('dataset/expression/mtxs/batchQC/info/stat/', sep = '\t', header = TRUE, check.names = FALSE)
head(count)
count <- subset(count, sample == 'Total')
count$genetype <- factor(count$genetype, levels = c('lincRNA-novel', 'lincRNA-known',
                                                    'antisense-novel',
                                                    'antisense-known',
                                                    'intervening-novel',
                                                    'intervening-known',
                                                    'flanking_region-novel',
                                                    'flanking_region-known'))

ggplot(count, aes(x="", y=count, fill=genetype))+
  geom_bar(width = 1, stat = "identity") + ggtitle('Immune-lncRNA genomic loci') + coord_polar('y', start = 0)

ggsave(filename = 'results/sub-results/gene-types.pie.pdf', width = 20, height = 17, units = "cm")
