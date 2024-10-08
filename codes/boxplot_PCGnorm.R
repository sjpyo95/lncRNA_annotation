library(ggplot2)

library(reshape2)

filename <-'exp-profile/mtx/qc/info/imm_15group.PCG_lncRNA.TPM.txt'
df <- read.table(filename, sep = '\t', header = TRUE,check.names = FALSE, row.names = 1)
head(df)
unique(df$type)
unique(subset(df, type == 'Known_lncRNA')$genetype)
summary(subset(df, type == 'Ig'))
#df[df=='Ig'] <- 'PCG'
df <- subset(df,type=='PCG' | type=='Known_lncRNA' | type == 'Novel_lncRNA')
head(df)
srt_df <- df[,c(-1,-3)]
head(srt_df)
dim(srt_df)
df_melt <- melt(srt_df)
head(df_melt)
dim(df_melt)
dim(subset(df_melt, variable == 'Plasmablasts' & type == 'Known_lncRNA'))
df_melt <- subset(df_melt, value >= 10)
summary(df_melt)

df_melt$logexp <- log(df_melt$value)


head(df_melt)

nrow(df_melt)

celltypes <- colnames(df)
celltypes

new_mdf <- data.frame(type = character(), variable = character(), value = numeric(), logexp = numeric(), normlog = numeric())
new_mdf
for (celltype in celltypes)  {
  celldf <- subset(df_melt, variable == celltype)
  cellpcg <- subset(celldf, type == 'PCG', select = logexp)
  cmed <- median(cellpcg$logexp)
  celldf$normlog <- celldf$logexp/cmed
  new_mdf <- rbind(new_mdf, celldf)
}
head(new_mdf)
#Q <- quantile(new_mdf$normlog, probs=c(.1, .90), na.rm = FALSE)
#iqr <- IQR(new_mdf$norm)
#up <-  Q[2]+1.5*iqr # Upper Range  
#low<- Q[1]-1.5*iqr # Lower Range
#new_mdf <- subset(new_mdf, normlog > (Q[1] - 1.5*iqr) & normlog < (Q[2]+1.5*iqr))
new_mdf$type <- factor(new_mdf$type, levels = c('PCG', 'Known_lncRNA', 'Novel_lncRNA'))
p <- ggplot(new_mdf, aes(x = variable, y=normlog, fill = type))+geom_boxplot()+
  theme_classic()+ theme(legend.title = element_blank(),
                         axis.text.x = element_text(colour = 'black', size = 10, angle = 90, hjust = 1, vjust = 0.0,),
                         legend.key.size = unit(0.3, 'cm'),
                         axis.title = element_text(size = 15),
                         legend.text = element_text(size = 10))+
  xlab('Immune cell-types') + ylab('log(TPM)') + ggtitle('Immune cells gene expression pattern (TPM ≥ 10)') + theme(plot.title = element_text(hjust = 0.5)) #+ coord_cartesian(ylim = c(0.8,2))
p

#pdfName <- paste(paste(strsplit(filename,".txt"), collapse=''),'.all.boxplot.pdf',sep = '')
ggsave(filename = 'exp-profile/mtx/qc/stat/graph/PCG_k_n_lncRNA.lognorm.10TPM.pdf', width = 30, height = 20, units = "cm")

