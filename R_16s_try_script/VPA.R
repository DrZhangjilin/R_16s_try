#setwd("G:/R/Rdata2")
setwd("G:/R/Rdata5")
library(openxlsx)
library(vegan)
library(adespatial)
otu <- read.table("otu_table2.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')
otu <- otu[,-ncol(otu)]
#删除多余的列
otu=t(otu)
#转置
#a <- read.xlsx("cor.xlsx",sheet=2,rowNames=T)
#b <- read.xlsx("cor.xlsx",sheet=3,rowNames=T)
#c <- read.xlsx("cor.xlsx",sheet=4,rowNames=T)
#d <- read.xlsx("cor.xlsx",sheet=5,rowNames=T)
#e <- read.xlsx("cor.xlsx",sheet=6,rowNames=T)
env.data <- read.table("envi.xls", row.names = 1, fill = T, header=T, sep="\t")
#otu<-otu[match(row.names(a),row.names(otu)),]
#print(rownames(a) == rownames(otu))
#print(rownames(b) == rownames(otu))
#print(rownames(c) == rownames(otu))
otu<-otu[match(row.names(env.data),row.names(otu)),]
print(rownames(env.data) == rownames(otu))

##by log
#a <- log1p(a)
#b <- log1p(b)
#c <- log1p(c)
env.data <- log1p(env.data)
##delete NA
#a <- na.omit(a)
#b <- na.omit(b)
#c <- na.omit(c)
env.data  <- na.omit(env.data )

#f <- cbind(a,b,c)
#环境变量的前向选择
#spe.f <- rda(otu,f)
#R2a.f <- RsquareAdj(spe.f)$adj.r.squared
#forward.sel(otu,f,adjR2thresh=R2a.f,nperm=9999)

#选择环境变量（可依据方差膨胀因子分析和前向分析）
env<- env.data[,-c(8,11,12,13,14,15,16,18,19)]
a <- env[,c(4,6,8,9,10,11)]#C、N
b <- env[,c(3,5,7)]#微生物
c <- env[,1:2]#环境

#ITS
#a <- a[,c(1,2,5,8,11)]
#b <- b[,-c(2,4)]


#f <- cbind(a,b,c)
#vpa <- varpart(otu,f,e,transfo='hel')
vpa <- varpart(otu,a,b,c,transfo='hel')

#vpa <- varpart(otu,a[,1:4],a[,5:11],transfo='hel')
vpa
plot(vpa,digits=2,Xnames=c("土壤C、N","微生物","环境") , id.size=1.5,bg = 2:4,
     cutoff=0)
#plot(vpa,digits=2,Xnames=c("环境因子","处理") , id.size=0.7,bg = 2:4,
     #cutoff=0)
#检验
anova(rda(otu,a),permutations=how(nperm=999))
anova(rda(otu,b),permutations=how(nperm=999))
anova(rda(otu,c),permutations=how(nperm=999))
anova(rda(otu,f),permutations=how(nperm=999))
anova(rda(otu,e),permutations=how(nperm=999))

