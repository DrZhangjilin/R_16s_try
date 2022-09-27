setwd("G:/R/Rdata4_new")
library(Hmisc)
genus <- read.table("otu_new_1.xls",row.names = 1,header=T,comment.char='',sep='\t')
genus_corr <- rcorr(t(genus), type = 'spearman')
r <- genus_corr$r
r[abs(r) < 0.8] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'genus_corr_C.matrix.csv', col.names = NA, sep = ',', quote = FALSE)

genus <- read.table("otu_new_2.xls",row.names = 1,header=T,comment.char='',sep='\t')
genus_corr <- rcorr(t(genus), type = 'spearman')
r <- genus_corr$r
r[abs(r) < 0.8] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'genus_corr_DL.matrix.csv', col.names = NA, sep = ',', quote = FALSE)

genus <- read.table("otu_new_3.xls",row.names = 1,header=T,comment.char='',sep='\t')
genus_corr <- rcorr(t(genus), type = 'spearman')
r <- genus_corr$r
r[abs(r) < 0.8] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'genus_corr_NL.matrix.csv', col.names = NA, sep = ',', quote = FALSE)

genus <- read.table("otu_new_4.xls",row.names = 1,header=T,comment.char='',sep='\t')
genus_corr <- rcorr(t(genus), type = 'spearman')
r <- genus_corr$r
r[abs(r) < 0.8] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'genus_corr_NP.matrix.csv', col.names = NA, sep = ',', quote = FALSE)

genus <- read.table("otu_new_5.xls",row.names = 1,header=T,comment.char='',sep='\t')
genus_corr <- rcorr(t(genus), type = 'spearman')
r <- genus_corr$r
r[abs(r) < 0.8] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'genus_corr_NPDL.matrix.csv', col.names = NA, sep = ',', quote = FALSE)

genus <- read.table("otu_new_6.xls",row.names = 1,header=T,comment.char='',sep='\t')
genus_corr <- rcorr(t(genus), type = 'spearman')
r <- genus_corr$r
r[abs(r) < 0.8] <- 0
p <- genus_corr$P
p <- p.adjust(p, method = 'BH')    #可选 p 值校正，这里使用 BH 法校正 p 值
p[p>=0.05] <- -1
p[p<0.05 & p>=0] <- 1
p[p==-1] <- 0
z <- r * p
diag(z) <- 0    #将相关矩阵中对角线中的值（代表了自相关）转为 0
head(z)[1:6,1:6]
write.table(data.frame(z, check.names = FALSE), 'genus_corr_NPNL.matrix.csv', col.names = NA, sep = ',', quote = FALSE)
