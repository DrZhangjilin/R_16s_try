setwd("G:/R/Rdata4_new")
library(vegan)

otu_or <- read.table("otu_table2.xls",skip=1,row.names = 1,header=T,comment.char='',sep='\t')
otu <- otu_or[,-ncol(otu_or)]

#事先过滤一些低丰度或低频的类群
#otu_absolute <- otu[which(rowSums(otu) >= 100), ]  #只保留绝对丰度总和高于 100 的otu

#otu1 <- otu_absolute #二者结合时释放
#otu1 <- otu
#otu1[otu1>0] <- 1
#otu_sample <- otu[which(rowSums(otu1) >= 18), ]    #只保留在 18 个及以上样本中出现的otu

#上述两种筛选可以单独使用也可以同时使用

#下面按照相对丰度筛选
#otu_relative <-(decostand(otu,'total',2)/36)#按列标准化otu,36个(24)样本，6(4)个处理，6 组重复，除以36是为了使总样本和为100%
otu_relative <-(decostand(otu,'total',2)/36)#单个处理
otu_relative <- otu_relative[which(rowSums(otu_relative) >= 0.0005), ] #只保留相对丰度总和高于 万分之5 的otu

otu_relative[otu_relative>0] <- 1
otu_sample <- otu_relative[which(rowSums(otu_relative) >= 19), ] #只保留在 18 个及以上样本中出现的otu

otu_sa = subset(otu_or, rownames(otu_or) %in% c(rownames(otu_sample)))#从原otu表中筛选符合条件的原始otu数据

write.table(otu_sa,"G:/R/Rdata4_new/2021.12.20_group_network/otu_all_0.0001.xls",sep="\t", quote=FALSE,col.names = NA)