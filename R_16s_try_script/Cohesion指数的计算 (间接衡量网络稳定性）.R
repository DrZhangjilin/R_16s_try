setwd("G:/R/Rdata4_new/1zipi_new")
library(ggplot2)
library(reshape2)
otu<-read.table('otu_table2.xls',row.names = 1,skip = 1,header = T,sep ="\t",comment.char = "")

C <- otu[,c("C4","C12","C17","C23","C30","C34")]
DL <- otu[,c("C1","C8","C14","C19","C27","C33")]
NL <- otu[,c("C2","C7","C15","C21","C25","C32")]
NP <- otu[,c("C3","C9","C13","C20","C26","C31")]
NPDL <- otu[,c("C6","C10","C16","C24","C29","C35")]
NPNL <- otu[,c("C5","C11","C18","C22","C28","C36")]

C1 <- read.table("G:/R/Rdata4_new/1zipi_new/node_C.csv",header = T,row.names=1,sep=",")
DL1 <- read.table("G:/R/Rdata4_new/1zipi_new/node_DL.csv",header = T,row.names=1,sep=",")
NL1 <- read.table("G:/R/Rdata4_new/1zipi_new/node_NL.csv",header = T,row.names=1,sep=",")
NP1 <- read.table("G:/R/Rdata4_new/1zipi_new/node_NP.csv",header = T,row.names=1,sep=",")
NPDL1 <- read.table("G:/R/Rdata4_new/1zipi_new/node_NPDL.csv",header = T,row.names=1,sep=",")
NPNL1 <- read.table("G:/R/Rdata4_new/1zipi_new/node_NPNL.csv",header = T,row.names=1,sep=",")

C2 = subset(C, rownames(C) %in% c(rownames(C1)))
DL2 = subset(DL, rownames(DL) %in% c(rownames(DL1)))
NL2 = subset(NL, rownames(NL) %in% c(rownames(NL1)))
NP2 = subset(NP, rownames(NP) %in% c(rownames(NP1)))
NPDL2 = subset(NPDL, rownames(NPDL) %in% c(rownames(NPDL1)))
NPNL2 = subset(NPNL, rownames(NPNL) %in% c(rownames(NPNL1)))

otu1 <- t(C2)
otu2 <- t(DL2)
otu3 <- t(NL2)
otu4 <- t(NP2)
otu5 <- t(NPDL2)
otu6 <- t(NPNL2)
##四个参数可调
##保留物种的阈值。如共有100个样本，0.1为只保留出现在大于10个样本中的物种。
pers.cutoff <- 0.00
#迭代次数。推荐>=200.但是速度会慢。
iter <- 999
##tax.shuffle=T,采用我文章中介绍的第一种随机化方法；tax.shuffle=F,采用第二种随机化方法
tax.shuffle <- T
##是否提供自己的相关性矩阵，代替零模型结果
use.custom.cors <- F

#读入数据，行为样本，列为物种
#b <- read.table("otu.txt",header = T,row.names = 1,sep='\t')
b <- otu1
##三个函数
#1.统计0的个数
zero <- function(vec){
  num.zero <- length(which(vec==0))
  return(num.zero)
}

#2.挑出来所有相关性的负值并求平均值。
neg.mean <- function(vector){
  neg.vals <- vector[which(vector<0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals)==0) n.mean <- 0
  return(n.mean)
}

#3.挑出来所有相关性的正值并求平均值。
pos.mean <- function(vector){
  pos.vals <- vector[which(vector>0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals)==0) p.mean <- 0
  return(p.mean)
}


#读入自己的相关矩阵
if(use.custom.cors==T){
  custom.cor.mat <- read.csv("your_path_here.csv",header = T,row.names=1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2]==dim(custom.cor.mat)[2])
}

#去掉都是0的行和列
c <- as.matrix(b)
c <- c[rowSums(c)>0,colSums(c)>0]

#计算样本的总序列数
rowsums.orig <- rowSums(c)

#确定物种个数的阈值
zero.cutoff <- ceiling(pers.cutoff*dim(c)[1])

#根据阈值筛选物种
d <- c[,apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
#移除丰度都是0的样本
d <- d[rowSums(d)>0,]

#如果自己提供了相关性表，也要根据阈值筛选
if(use.custom.cors==T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c,2,zero) < (dim(c)[1]-zero.cutoff), apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
}

#计算丰度百分比
rel.d <- d/rowsums.orig

#计算实际相关矩阵
cor.mat.true <- cor(rel.d)

#计算零模型
med.tax.cors <- vector()
set.seed(999)
if(use.custom.cors==F){
  if(tax.shuffle){#第一种零模型计算方法，打乱列
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- matrix(numeric(0),dim(rel.d)[1],dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #对每个物种
        for(j in 1:dim(rel.d)[2]){
          #先替换所有的值
          perm.rel.d[,j] <- sample(rel.d[,j])
        }
        
        #focal物种维持原状不随机
        perm.rel.d[,which.taxon] <- rel.d[,which.taxon]
        
        #算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  } else{#第二种零模型。打乱样本内所有丰度，打乱行
  for(which.taxon in 1:dim(rel.d)[2]){
    perm.cor.vec.mat <- vector()
    
    for(i in 1:iter){
      perm.rel.d <- rel.d
      #对每个物种
      for(j in 1:dim(rel.d)[1]){
        #先取出大于0的值
        which.replace <- which(rel.d[j,]>0)
        #focal taxon 去掉，不随机
        which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
        #对样本内，不是0的，且不是focal taxon的物种进行随机化
        perm.rel.d[j,which.replace.nonfocal] <- sample(rel.d[j,which.replace.nonfocal])
      }
      
      #计算相关性
      cor.mat.null <- cor(perm.rel.d)
      
      #保存每次对于focal物种迭代的结果
      perm.cor.vec.mat <- cblind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
      
    }
    med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
    #运行情况
    if(which.taxon %% 20==0){print(which.taxon)}
    }
  }
}
  
#实际相关-零模型相关
ifelse(use.custom.cors==T,{
  obs.exp.cors.mat <- custom.cor.mat.sub
},{
  obs.exp.cors.mat <- cor.mat.true-med.tax.cors
})

#计算正负相关的连通性
connectedness.pos <- apply(obs.exp.cors.mat,2,pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat,2,neg.mean)

#根据定义，计算cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

###输出
#output <- list(connectedness.neg,connectedness.pos,cohesion.neg,cohesion.pos)
#names(output) <- c("Negative Connectedness","Positive Connectedness","Negative Cohesion","Positive Cohesion")
output <- list(cohesion.neg,cohesion.pos)
names(output) <- c("Negative Cohesion","Positive Cohesion")
output_table <- as.data.frame(output)

print(output_table)

write.table(output_table,"Cohesion_C.csv",sep = ',', row.names = T, quote = FALSE)

###########################################################################################################

#读入数据，行为样本，列为物种
#b <- read.table("otu.txt",header = T,row.names = 1,sep='\t')
b <- otu2
##三个函数
#1.统计0的个数
zero <- function(vec){
  num.zero <- length(which(vec==0))
  return(num.zero)
}

#2.挑出来所有相关性的负值并求平均值。
neg.mean <- function(vector){
  neg.vals <- vector[which(vector<0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals)==0) n.mean <- 0
  return(n.mean)
}

#3.挑出来所有相关性的正值并求平均值。
pos.mean <- function(vector){
  pos.vals <- vector[which(vector>0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals)==0) p.mean <- 0
  return(p.mean)
}


#读入自己的相关矩阵
if(use.custom.cors==T){
  custom.cor.mat <- read.csv("your_path_here.csv",header = T,row.names=1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2]==dim(custom.cor.mat)[2])
}

#去掉都是0的行和列
c <- as.matrix(b)
c <- c[rowSums(c)>0,colSums(c)>0]

#计算样本的总序列数
rowsums.orig <- rowSums(c)

#确定物种个数的阈值
zero.cutoff <- ceiling(pers.cutoff*dim(c)[1])

#根据阈值筛选物种
d <- c[,apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
#移除丰度都是0的样本
d <- d[rowSums(d)>0,]

#如果自己提供了相关性表，也要根据阈值筛选
if(use.custom.cors==T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c,2,zero) < (dim(c)[1]-zero.cutoff), apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
}

#计算丰度百分比
rel.d <- d/rowsums.orig

#计算实际相关矩阵
cor.mat.true <- cor(rel.d)

#计算零模型
med.tax.cors <- vector()
set.seed(999)
if(use.custom.cors==F){
  if(tax.shuffle){#第一种零模型计算方法，打乱列
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- matrix(numeric(0),dim(rel.d)[1],dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #对每个物种
        for(j in 1:dim(rel.d)[2]){
          #先替换所有的值
          perm.rel.d[,j] <- sample(rel.d[,j])
        }
        
        #focal物种维持原状不随机
        perm.rel.d[,which.taxon] <- rel.d[,which.taxon]
        
        #算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  } else{#第二种零模型。打乱样本内所有丰度，打乱行
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- rel.d
        #对每个物种
        for(j in 1:dim(rel.d)[1]){
          #先取出大于0的值
          which.replace <- which(rel.d[j,]>0)
          #focal taxon 去掉，不随机
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          #对样本内，不是0的，且不是focal taxon的物种进行随机化
          perm.rel.d[j,which.replace.nonfocal] <- sample(rel.d[j,which.replace.nonfocal])
        }
        
        #计算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cblind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  }
}

#实际相关-零模型相关
ifelse(use.custom.cors==T,{
  obs.exp.cors.mat <- custom.cor.mat.sub
},{
  obs.exp.cors.mat <- cor.mat.true-med.tax.cors
})

#计算正负相关的连通性
connectedness.pos <- apply(obs.exp.cors.mat,2,pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat,2,neg.mean)

#根据定义，计算cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

###输出
#output <- list(connectedness.neg,connectedness.pos,cohesion.neg,cohesion.pos)
#names(output) <- c("Negative Connectedness","Positive Connectedness","Negative Cohesion","Positive Cohesion")
output <- list(cohesion.neg,cohesion.pos)
names(output) <- c("Negative Cohesion","Positive Cohesion")
output_table <- as.data.frame(output)

print(output_table)

write.table(output_table,"Cohesion_DL.csv",sep = ',', row.names = T, quote = FALSE)

###############################################################################################################################

#读入数据，行为样本，列为物种
#b <- read.table("otu.txt",header = T,row.names = 1,sep='\t')
b <- otu3
##三个函数
#1.统计0的个数
zero <- function(vec){
  num.zero <- length(which(vec==0))
  return(num.zero)
}

#2.挑出来所有相关性的负值并求平均值。
neg.mean <- function(vector){
  neg.vals <- vector[which(vector<0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals)==0) n.mean <- 0
  return(n.mean)
}

#3.挑出来所有相关性的正值并求平均值。
pos.mean <- function(vector){
  pos.vals <- vector[which(vector>0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals)==0) p.mean <- 0
  return(p.mean)
}


#读入自己的相关矩阵
if(use.custom.cors==T){
  custom.cor.mat <- read.csv("your_path_here.csv",header = T,row.names=1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2]==dim(custom.cor.mat)[2])
}

#去掉都是0的行和列
c <- as.matrix(b)
c <- c[rowSums(c)>0,colSums(c)>0]

#计算样本的总序列数
rowsums.orig <- rowSums(c)

#确定物种个数的阈值
zero.cutoff <- ceiling(pers.cutoff*dim(c)[1])

#根据阈值筛选物种
d <- c[,apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
#移除丰度都是0的样本
d <- d[rowSums(d)>0,]

#如果自己提供了相关性表，也要根据阈值筛选
if(use.custom.cors==T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c,2,zero) < (dim(c)[1]-zero.cutoff), apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
}

#计算丰度百分比
rel.d <- d/rowsums.orig

#计算实际相关矩阵
cor.mat.true <- cor(rel.d)

#计算零模型
med.tax.cors <- vector()
set.seed(999)
if(use.custom.cors==F){
  if(tax.shuffle){#第一种零模型计算方法，打乱列
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- matrix(numeric(0),dim(rel.d)[1],dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #对每个物种
        for(j in 1:dim(rel.d)[2]){
          #先替换所有的值
          perm.rel.d[,j] <- sample(rel.d[,j])
        }
        
        #focal物种维持原状不随机
        perm.rel.d[,which.taxon] <- rel.d[,which.taxon]
        
        #算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  } else{#第二种零模型。打乱样本内所有丰度，打乱行
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- rel.d
        #对每个物种
        for(j in 1:dim(rel.d)[1]){
          #先取出大于0的值
          which.replace <- which(rel.d[j,]>0)
          #focal taxon 去掉，不随机
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          #对样本内，不是0的，且不是focal taxon的物种进行随机化
          perm.rel.d[j,which.replace.nonfocal] <- sample(rel.d[j,which.replace.nonfocal])
        }
        
        #计算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cblind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  }
}

#实际相关-零模型相关
ifelse(use.custom.cors==T,{
  obs.exp.cors.mat <- custom.cor.mat.sub
},{
  obs.exp.cors.mat <- cor.mat.true-med.tax.cors
})

#计算正负相关的连通性
connectedness.pos <- apply(obs.exp.cors.mat,2,pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat,2,neg.mean)

#根据定义，计算cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

###输出
#output <- list(connectedness.neg,connectedness.pos,cohesion.neg,cohesion.pos)
#names(output) <- c("Negative Connectedness","Positive Connectedness","Negative Cohesion","Positive Cohesion")
output <- list(cohesion.neg,cohesion.pos)
names(output) <- c("Negative Cohesion","Positive Cohesion")
output_table <- as.data.frame(output)

print(output_table)

write.table(output_table,"Cohesion_NL.csv",sep = ',', row.names = T, quote = FALSE)

############################################################################################################################################

#读入数据，行为样本，列为物种
#b <- read.table("otu.txt",header = T,row.names = 1,sep='\t')
b <- otu4
##三个函数
#1.统计0的个数
zero <- function(vec){
  num.zero <- length(which(vec==0))
  return(num.zero)
}

#2.挑出来所有相关性的负值并求平均值。
neg.mean <- function(vector){
  neg.vals <- vector[which(vector<0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals)==0) n.mean <- 0
  return(n.mean)
}

#3.挑出来所有相关性的正值并求平均值。
pos.mean <- function(vector){
  pos.vals <- vector[which(vector>0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals)==0) p.mean <- 0
  return(p.mean)
}


#读入自己的相关矩阵
if(use.custom.cors==T){
  custom.cor.mat <- read.csv("your_path_here.csv",header = T,row.names=1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2]==dim(custom.cor.mat)[2])
}

#去掉都是0的行和列
c <- as.matrix(b)
c <- c[rowSums(c)>0,colSums(c)>0]

#计算样本的总序列数
rowsums.orig <- rowSums(c)

#确定物种个数的阈值
zero.cutoff <- ceiling(pers.cutoff*dim(c)[1])

#根据阈值筛选物种
d <- c[,apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
#移除丰度都是0的样本
d <- d[rowSums(d)>0,]

#如果自己提供了相关性表，也要根据阈值筛选
if(use.custom.cors==T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c,2,zero) < (dim(c)[1]-zero.cutoff), apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
}

#计算丰度百分比
rel.d <- d/rowsums.orig

#计算实际相关矩阵
cor.mat.true <- cor(rel.d)

#计算零模型
med.tax.cors <- vector()
set.seed(999)
if(use.custom.cors==F){
  if(tax.shuffle){#第一种零模型计算方法，打乱列
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- matrix(numeric(0),dim(rel.d)[1],dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #对每个物种
        for(j in 1:dim(rel.d)[2]){
          #先替换所有的值
          perm.rel.d[,j] <- sample(rel.d[,j])
        }
        
        #focal物种维持原状不随机
        perm.rel.d[,which.taxon] <- rel.d[,which.taxon]
        
        #算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  } else{#第二种零模型。打乱样本内所有丰度，打乱行
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- rel.d
        #对每个物种
        for(j in 1:dim(rel.d)[1]){
          #先取出大于0的值
          which.replace <- which(rel.d[j,]>0)
          #focal taxon 去掉，不随机
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          #对样本内，不是0的，且不是focal taxon的物种进行随机化
          perm.rel.d[j,which.replace.nonfocal] <- sample(rel.d[j,which.replace.nonfocal])
        }
        
        #计算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cblind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  }
}

#实际相关-零模型相关
ifelse(use.custom.cors==T,{
  obs.exp.cors.mat <- custom.cor.mat.sub
},{
  obs.exp.cors.mat <- cor.mat.true-med.tax.cors
})

#计算正负相关的连通性
connectedness.pos <- apply(obs.exp.cors.mat,2,pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat,2,neg.mean)

#根据定义，计算cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

###输出
#output <- list(connectedness.neg,connectedness.pos,cohesion.neg,cohesion.pos)
#names(output) <- c("Negative Connectedness","Positive Connectedness","Negative Cohesion","Positive Cohesion")
output <- list(cohesion.neg,cohesion.pos)
names(output) <- c("Negative Cohesion","Positive Cohesion")
output_table <- as.data.frame(output)

print(output_table)

write.table(output_table,"Cohesion_NP.csv",sep = ',', row.names = T, quote = FALSE)
#####################################################################################################################################

#读入数据，行为样本，列为物种
#b <- read.table("otu.txt",header = T,row.names = 1,sep='\t')
b <- otu5
##三个函数
#1.统计0的个数
zero <- function(vec){
  num.zero <- length(which(vec==0))
  return(num.zero)
}

#2.挑出来所有相关性的负值并求平均值。
neg.mean <- function(vector){
  neg.vals <- vector[which(vector<0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals)==0) n.mean <- 0
  return(n.mean)
}

#3.挑出来所有相关性的正值并求平均值。
pos.mean <- function(vector){
  pos.vals <- vector[which(vector>0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals)==0) p.mean <- 0
  return(p.mean)
}


#读入自己的相关矩阵
if(use.custom.cors==T){
  custom.cor.mat <- read.csv("your_path_here.csv",header = T,row.names=1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2]==dim(custom.cor.mat)[2])
}

#去掉都是0的行和列
c <- as.matrix(b)
c <- c[rowSums(c)>0,colSums(c)>0]

#计算样本的总序列数
rowsums.orig <- rowSums(c)

#确定物种个数的阈值
zero.cutoff <- ceiling(pers.cutoff*dim(c)[1])

#根据阈值筛选物种
d <- c[,apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
#移除丰度都是0的样本
d <- d[rowSums(d)>0,]

#如果自己提供了相关性表，也要根据阈值筛选
if(use.custom.cors==T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c,2,zero) < (dim(c)[1]-zero.cutoff), apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
}

#计算丰度百分比
rel.d <- d/rowsums.orig

#计算实际相关矩阵
cor.mat.true <- cor(rel.d)

#计算零模型
med.tax.cors <- vector()
set.seed(999)
if(use.custom.cors==F){
  if(tax.shuffle){#第一种零模型计算方法，打乱列
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- matrix(numeric(0),dim(rel.d)[1],dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #对每个物种
        for(j in 1:dim(rel.d)[2]){
          #先替换所有的值
          perm.rel.d[,j] <- sample(rel.d[,j])
        }
        
        #focal物种维持原状不随机
        perm.rel.d[,which.taxon] <- rel.d[,which.taxon]
        
        #算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  } else{#第二种零模型。打乱样本内所有丰度，打乱行
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- rel.d
        #对每个物种
        for(j in 1:dim(rel.d)[1]){
          #先取出大于0的值
          which.replace <- which(rel.d[j,]>0)
          #focal taxon 去掉，不随机
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          #对样本内，不是0的，且不是focal taxon的物种进行随机化
          perm.rel.d[j,which.replace.nonfocal] <- sample(rel.d[j,which.replace.nonfocal])
        }
        
        #计算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cblind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  }
}

#实际相关-零模型相关
ifelse(use.custom.cors==T,{
  obs.exp.cors.mat <- custom.cor.mat.sub
},{
  obs.exp.cors.mat <- cor.mat.true-med.tax.cors
})

#计算正负相关的连通性
connectedness.pos <- apply(obs.exp.cors.mat,2,pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat,2,neg.mean)

#根据定义，计算cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

###输出
#output <- list(connectedness.neg,connectedness.pos,cohesion.neg,cohesion.pos)
#names(output) <- c("Negative Connectedness","Positive Connectedness","Negative Cohesion","Positive Cohesion")
output <- list(cohesion.neg,cohesion.pos)
names(output) <- c("Negative Cohesion","Positive Cohesion")
output_table <- as.data.frame(output)

print(output_table)

write.table(output_table,"Cohesion_NPDL.csv",sep = ',', row.names = T, quote = FALSE)

################################################################################################################################

#读入数据，行为样本，列为物种
#b <- read.table("otu.txt",header = T,row.names = 1,sep='\t')
b <- otu6
##三个函数
#1.统计0的个数
zero <- function(vec){
  num.zero <- length(which(vec==0))
  return(num.zero)
}

#2.挑出来所有相关性的负值并求平均值。
neg.mean <- function(vector){
  neg.vals <- vector[which(vector<0)]
  n.mean <- mean(neg.vals)
  if(length(neg.vals)==0) n.mean <- 0
  return(n.mean)
}

#3.挑出来所有相关性的正值并求平均值。
pos.mean <- function(vector){
  pos.vals <- vector[which(vector>0)]
  p.mean <- mean(pos.vals)
  if(length(pos.vals)==0) p.mean <- 0
  return(p.mean)
}


#读入自己的相关矩阵
if(use.custom.cors==T){
  custom.cor.mat <- read.csv("your_path_here.csv",header = T,row.names=1)
  custom.cor.mat <- as.matrix(custom.cor.mat)
  #check that correlation matrix and abundance matrix have the same dimension
  print(dim(b)[2]==dim(custom.cor.mat)[2])
}

#去掉都是0的行和列
c <- as.matrix(b)
c <- c[rowSums(c)>0,colSums(c)>0]

#计算样本的总序列数
rowsums.orig <- rowSums(c)

#确定物种个数的阈值
zero.cutoff <- ceiling(pers.cutoff*dim(c)[1])

#根据阈值筛选物种
d <- c[,apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
#移除丰度都是0的样本
d <- d[rowSums(d)>0,]

#如果自己提供了相关性表，也要根据阈值筛选
if(use.custom.cors==T){
  custom.cor.mat.sub <- custom.cor.mat[apply(c,2,zero) < (dim(c)[1]-zero.cutoff), apply(c,2,zero)<(dim(c)[1]-zero.cutoff)]
}

#计算丰度百分比
rel.d <- d/rowsums.orig

#计算实际相关矩阵
cor.mat.true <- cor(rel.d)

#计算零模型
med.tax.cors <- vector()
set.seed(999)
if(use.custom.cors==F){
  if(tax.shuffle){#第一种零模型计算方法，打乱列
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- matrix(numeric(0),dim(rel.d)[1],dim(rel.d)[2])
        rownames(perm.rel.d) <- rownames(rel.d)
        colnames(perm.rel.d) <- colnames(rel.d)
        
        #对每个物种
        for(j in 1:dim(rel.d)[2]){
          #先替换所有的值
          perm.rel.d[,j] <- sample(rel.d[,j])
        }
        
        #focal物种维持原状不随机
        perm.rel.d[,which.taxon] <- rel.d[,which.taxon]
        
        #算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cbind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  } else{#第二种零模型。打乱样本内所有丰度，打乱行
    for(which.taxon in 1:dim(rel.d)[2]){
      perm.cor.vec.mat <- vector()
      
      for(i in 1:iter){
        perm.rel.d <- rel.d
        #对每个物种
        for(j in 1:dim(rel.d)[1]){
          #先取出大于0的值
          which.replace <- which(rel.d[j,]>0)
          #focal taxon 去掉，不随机
          which.replace.nonfocal <- which.replace[!(which.replace %in% which.taxon)]
          #对样本内，不是0的，且不是focal taxon的物种进行随机化
          perm.rel.d[j,which.replace.nonfocal] <- sample(rel.d[j,which.replace.nonfocal])
        }
        
        #计算相关性
        cor.mat.null <- cor(perm.rel.d)
        
        #保存每次对于focal物种迭代的结果
        perm.cor.vec.mat <- cblind(perm.cor.vec.mat,cor.mat.null[,which.taxon])
        
      }
      med.tax.cors <- cbind(med.tax.cors,apply(perm.cor.vec.mat,1,median))
      #运行情况
      if(which.taxon %% 20==0){print(which.taxon)}
    }
  }
}

#实际相关-零模型相关
ifelse(use.custom.cors==T,{
  obs.exp.cors.mat <- custom.cor.mat.sub
},{
  obs.exp.cors.mat <- cor.mat.true-med.tax.cors
})

#计算正负相关的连通性
connectedness.pos <- apply(obs.exp.cors.mat,2,pos.mean)
connectedness.neg <- apply(obs.exp.cors.mat,2,neg.mean)

#根据定义，计算cohesion
cohesion.pos <- rel.d %*% connectedness.pos
cohesion.neg <- rel.d %*% connectedness.neg

###输出
#output <- list(connectedness.neg,connectedness.pos,cohesion.neg,cohesion.pos)
#names(output) <- c("Negative Connectedness","Positive Connectedness","Negative Cohesion","Positive Cohesion")
output <- list(cohesion.neg,cohesion.pos)
names(output) <- c("Negative Cohesion","Positive Cohesion")
output_table <- as.data.frame(output)

print(output_table)

write.table(output_table,"Cohesion_NPNL.csv",sep = ',', row.names = T, quote = FALSE)

##############################################################################################

C1 <- read.table("G:/R/Rdata4_new/1zipi_new/Cohesion_C.csv",header = T,row.names=1,sep=",")
C2 <- read.table("G:/R/Rdata4_new/1zipi_new/Cohesion_DL.csv",header = T,row.names=1,sep=",")
C3 <- read.table("G:/R/Rdata4_new/1zipi_new/Cohesion_NL.csv",header = T,row.names=1,sep=",")
C4 <- read.table("G:/R/Rdata4_new/1zipi_new/Cohesion_NP.csv",header = T,row.names=1,sep=",")
C5 <- read.table("G:/R/Rdata4_new/1zipi_new/Cohesion_NPDL.csv",header = T,row.names=1,sep=",")
C6 <- read.table("G:/R/Rdata4_new/1zipi_new/Cohesion_NPNL.csv",header = T,row.names=1,sep=",")

Cohesion <- rbind(C1,C2,C3,C4,C5,C6)
Cohesion$Treat <- c(rep("C",6),rep("DL",6),rep("NL",6),rep("NP",6),rep("NPDL",6),rep("NPNL",6))
Cohesion$`|Negative.Cohesion|` <- abs(Cohesion$Negative.Cohesion)
Cohesion$Total.Cohesion <- Cohesion$Positive.Cohesion+Cohesion$`|Negative.Cohesion|`
Cohesion$`Negative/Positive.Cohesion` <- Cohesion$`|Negative.Cohesion|` / Cohesion$Positive.Cohesion
colnames(Cohesion) <- c("Negative Cohesion","Positive Cohesion","Treat","|Negative Cohesion|","Total Cohesion","Negative/Positive Cohesion")
write.table(Cohesion,"Cohesion_all.csv",sep = ',', row.names = T, quote = FALSE)

dat <- read.table("G:/R/Rdata4_new/1zipi_new/Cohesion_all.csv",header = T,row.names=1,sep=",")
dat$Sample <- rownames(dat)
dat <- dat[,-1]
dat1 <- melt(dat)

ONE_LSD <- function(data,group,compare,value){
  library(agricolae)
  
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    # sub_dat <- subset(data,group == i)
    sub_dat <- data[data[,group]==i,]
    # fit <- aov(value ~ compare,sub_dat)
    fit <- aov(sub_dat[,value] ~ sub_dat[,compare] )
    out <- LSD.test(fit,'sub_dat[, compare]',p.adj='BH')#进行了p值校正
    #out$groups就可获取多重比较字母列表
    out$groups$type <- i
    out$groups$compare <- rownames(out$groups)
    
    a <- rbind(a,merge(out$means[,1:2], out$groups,by='sub_dat[, value]'))
  }
  names(a) <- c('mean','std','lsd',group,compare)
  return(a)
}

df <- ONE_LSD(dat1,'variable','Treat','value')
head(df)

p <- ggplot(dat1)+
geom_jitter(aes(x=Treat,y=value,color=Treat),width=0.25,shape=20,size=2,show.legend = F)+#散点
stat_boxplot(aes(x=Treat,y=value),geom="errorbar",size=0.6,width=0.25)+
geom_boxplot(aes(x=Treat,y=value,fill=Treat))+#color非填充，fill填充
#outlier.colour="red", outlier.shape=8, 
#outlier.size=4,outlier.stroke =0.5,
#outlier.alpha = 0.5,outlier.fill = "red")+
#scale_y_continuous(expand = c(0.1,0))+
geom_text(data=df,aes(x=Treat,y=mean+1.3*std,label=lsd,vjust=1),size=10,family="serif")+
  facet_wrap(.~variable,scales = "free_y")+
  stat_summary(aes(x=Treat,y=value),fun="mean",geom="point",shape=23,size=1,fill="white")+
  labs(x='',y='Cohesion')+
  #ggprism::theme_prism()+
  #theme(axis.text.x = element_text(angle = 45))+
  #ggprism函数中的theme
  theme_bw()+
  guides(fill=guide_legend(title = "Treat"))+
  theme(text=element_text(size=28,  family="serif",face = "bold"),
        axis.text.x = element_text(size=20,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=20,  family="serif",colour = "black"))
p
ggsave("./Positive Cohesion.pdf", p, width = 400, height = 400, units = "mm")