##PCA test##
install.packages("FactoMineR")
install.packages("ggplot2")
install.packages("devtools")
library(devtools)
install_github("protViz/quantable")
library(FactoMineR)
library(ggplot2)
library(quantable)
library(matrixStats)
tpm_mtx <- read.table('exp-profile/mtx/qc/info/imm.PCG.TPM.txt', sep = '\t', header = TRUE, row.names = 1, check.names = FALSE)
head(tpm_mtx)
dim(tpm_mtx)
tpm_mtx <- subset(tpm_mtx, type=='Known_lncRNA' | type == 'Novel_lncRNA')
mtx <- tpm_mtx[,c(-1:-3)]

head(mtx)

srt_mtx <- subset(mtx, apply(mtx[,-1],1,function(x) any(x >= 10)))#10TPM 이상이 하나 이상의 cell-type에서 나타나는 gene만 가져오기
dim(srt_mtx)

std <-transform(srt_mtx, SD=apply(srt_mtx,1, sd, na.rm = TRUE))
head(std)
dim(std)

#cv<- CV(srt_mtx[,-1], top = 0, na.rm = TRUE)
#cv<- as.data.frame(cv)

cv_mtx <- cbind(srt_mtx, cv)

hist(cv)

hivar_mtx <- subset(cv_mtx, cv >= 100)
head(hivar_mtx)
dim(hivar_mtx)
tail(hivar_mtx)
#PCA

res <- PCA(t(hivar_mtx[,c(-1,-35)]))
summary(res)
pc <- as.data.frame(res$ind$coord[,1:2])
head(pc)


pcPlot <- ggplot(pc, aes(x = Dim.1, y = Dim.2, color = rownames(pc)))+ geom_point(size = 5.5)+
  xlab(paste('PC1 (',as.character(round(res$eig[1,2], 1)), ')%',sep = ''))+ ylab(paste('PC2 (',as.character(round(res$eig[2,2],1)), ')%', sep = ''))+
  #scale_color_manual(values = c("9999FF", '3333FF',  "9933FF",  '6600FF', 'FF9966', )) + 
  ggtitle('LncRNAs with high CV between immune cells (CV >= 200)')+#xlim(-10,90)+ylim(-10,90)+
  theme_classic()#+ geom_text(aes(y=Dim.2+1, label=rownames(pc)), size=5)

pcPlot

ggsave(filename = 'exp-profile/mtx/qc/stat/graph/200_higher_cv_lncRNAs.PCA.pdf', width = 20, height = 14, units = "cm")

write.table(hivar_mtx, file = 'exp-profile/mtx/qc/info/imm.cv200.TPM.txt', append = TRUE, sep = '\t', eol = '\n', na = 'NA', row.names = TRUE, col.names = TRUE, quote = FALSE)


#srt_mtx[1,c(-1,-35)]
#mean(as.numeric(srt_mtx[1,c(-1,-35)]))
#exps <- srt_mtx[2,c(-1,-35)]
#mean <- mean(as.numeric(srt_mtx[2,c(-1,-35)]))
#mapply('/', exps, mean)
#variable genes selection
#srt_mtx$cv <- 0
#srt_mtx[1,c(-1,-35)]
#library(matrixStats)
#for (r in 1:nrow(srt_mtx)) {
#  exps <- as.numeric(srt_mtx[r,c(-1,-35)])
#  meanval <- 
#  cv_val <- sd(as.numeric(srt_mtx[r,c(-1,-35)]))/mean(as.numeric(srt_mtx[r,c(-1,-35)])) * 100
#  srt_mtx$cv[r] <- cv_val
#}

#srt_mtx <- srt_mtx[order(-srt_mtx$cv),]
#srt_mtx$var <- rowVars(as.matrix(srt_mtx))
#head(srt_mtx)
#srt_mtx <- srt_mtx[order(-srt_mtx$var),]
#head(srt_mtx)
#dim(srt_mtx)
#srt_mtx$var <- NULL