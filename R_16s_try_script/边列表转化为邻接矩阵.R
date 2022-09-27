setwd("F:/R/Rdata2")
library(igraph)
dat <- read.csv("F:/R/Rdata2/result（舍去万分之一）/C_Gephi_edge.csv")
dat <- as.data.frame(dat)
ed <- dat[,1:2]
ed <- as.matrix(ed)#graph.edgelist（）只能识别矩阵
edgelist <- graph.edgelist(ed,directed=FALSE)#无向，directed=FALSE
adjacency <- get.adjacency(edgelist)
adj_matrix <- as.matrix(adjacency)#把S4类结果保存成矩阵
adj <- data.frame(adj_matrix,check.names=FALSE)#矩阵保存为数据框方便输出
write.table(adj,"network_adj_matrix_unwC.txt",sep="\t",col.names = NA,quote=FALSE)

