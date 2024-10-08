df <- read.table('exp-profile/mtx/qc/info/ubiq.10up.top5.mtx.txt', sep = '\t', header = TRUE,check.names = FALSE)
head(df)
ubi <- c('LNC-MIR3687-1', 'MALAT1', 'PSMA3-AS1', 'SNHG1', 'NEAT1', 'FGD5-AS1')
df <- subset(df, symbol %in% ubi)
dim(df)
df_melt <- melt(df)
head(df_melt)
df_melt$logval <- log2(df_melt$value)
ggplot(data = df_melt, aes(x=variable, y = logval, group = symbol))+
  geom_line(aes(color=symbol, size = 1.3)) + geom_point(aes(color=symbol, size = 3)) + theme_classic() + 
  theme(legend.title = element_blank(),
                                     axis.text.x = element_text(colour = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.0,),
                                     legend.key.size = unit(0.3, 'cm'),
                                     axis.title = element_text(size = 15),
                                     legend.text = element_text(size = 10))+
  xlab('Samples') + ylab('log2count') + ggtitle('Human Immune cells expression (count > 0)') + theme(plot.title = element_text(hjust = 0.5))
