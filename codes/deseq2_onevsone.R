#one vs one
dds <- DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~group4)
dds
dds <- DESeq(dds)
group <- coldata$group4

grouptype <- group[!duplicated(group)]
grouptype

grouptype <- group[!duplicated(group)]
#for (i in c(1:(length(grouptype)))) {
#  others <- grouptype[-i]
#  grouptype[i]
#  newdir <- paste(grouptype[i],'/',sep = '')
#  if (dir.exists(paste('dataset/expression/degs/one_vs_one/', newdir, sep = ''))) next
#  print(grouptype[i])
#  for (j in c(1:length(others))) {
#    
#    res <- results(dds, contrast = c('group1', grouptype[i], others[j]))
#    res.filt <- subset(res, padj <= 0.05)
#    res.filt <- subset(res.filt, log2FoldChange >= 1)
#    filename <- paste(grouptype[i],'_', others[j],'.deg.txt', sep = '')


#     dir.create(paste('dataset/expression/degs/one_vs_one/', newdir,'/', sep = ''),showWarnings = FALSE)
#    write.table(res.filt, file = paste('dataset/expression/degs/one_vs_one/',newdir,filename, sep = ''), append = FALSE, sep = '\t', eol = '\n', na = 'NA', row.names = TRUE, col.names = TRUE, quote = FALSE)
#  }
#}

for (i in c(1:(length(grouptype)))) {
  others <- grouptype[-i]
  newdir <- paste(grouptype[i],'/',sep = '')
  if (dir.exists(paste('dataset/expression/degs/one_vs_one/group3/', newdir, sep = ''))) next
  print(paste('dataset/expression/degs/one_vs_one/', newdir,'/', sep = ''))
  for (j in c(1:length(others))) {
    
    res <- results(dds, contrast = c('group3', grouptype[i], others[j]))
    res.filt <- subset(res, padj <= 0.05)
    res.filt <- subset(res.filt, log2FoldChange >= 1)
    filename <- paste(grouptype[i],'_', others[j],'.deg.txt', sep = '')
    
    
    dir.create(paste('dataset/expression/degs/one_vs_one/group3/', newdir,'/', sep = ''),showWarnings = FALSE)
    write.table(res.filt, file = paste('dataset/expression/degs/one_vs_one/group3/',newdir,filename, sep = ''), append = FALSE, sep = '\t', eol = '\n', na = 'NA', row.names = TRUE, col.names = TRUE, quote = FALSE)
  }
}



res <- results(dds, contrast = c('group4', grouptype[], others[j]))