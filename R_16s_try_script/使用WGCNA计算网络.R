setwd("G:/R/Rdata4_new/2021.12.22WGCNA")
library(WGCNA)
library("multtest")
library(reshape2)
library(tidyverse)

#calNetwork函数共包含了6个参数，
#第一个参数data = data为数据框；
#第二个参数filter = FALSE为逻辑参数；
#第三个参数n = 3只有在filter = TRUE时才有意义或者才需要定义；
#第四个参数method = “spearman”指计算选用的方法，共包含了三种方法，分别是：pearson, kendall和spearman；
#第五个参数cutoff_cor = 0.6
#第六个参数cutoff_p = 0.05，指过滤的相关系数和显著性p值，默认只有相关系数的绝对值大于0.6且显著性小于0.05的边才会被保留。


#函数使用：
#1）在所有函数的参数不修改的默认的情况下，该函数使用的网络推断方法为spearman相关分析；
#2）默认筛选OTU或者物种的规则为OTU或物种所在的样方，不为0的个数大于一半以上；例如如果有50个样方和100个物种，物种1到100只有在25个样方中的值均不为0时才会进入函数计算，否则将被过滤掉；
#3）默认filter为FALSE指以2）中默认的方式进行网络计算，当然如果选择filter = TRUE则意味着需要自己定义不为0的样方个数，这里需要对n进行定义，如n = 3；
#4）默认使用的相关性计算方法为spearman，默认使用的cutoff相关性和p值分别为0.6和0.05。

calNetwork <- function(data = data,filter = FALSE,
                       n = 3,method = "spearman",
                       cutoff_cor = 0.6, cutoff_p = 0.05){
  require(WGCNA)
  require("multtest")
  require(reshape2)
  require(tidyverse)
  
  cutoff_cor <- cutoff_cor
  cutoff_p <- cutoff_p
  if(filter){
    dat2 <- data[,colSums(data!= "0") > n]
  } else{
    dat2 <- data[,colSums(data!= "0") > round(nrow(data)/2,0)]
  }
  set.seed(123)
  dat_net <- corAndPvalue(dat2,method = method)#calculate correlation and significance
  dat_cor <- dat_net$cor#select correlation
  dat_p <- dat_net$p#select p values
  
  dat_cor[upper.tri(dat_cor,diag=TRUE)]=NA
  dat_p[upper.tri(dat_p,diag=TRUE)]=NA
  dat_cor_melt <- melt(dat_cor)
  dat_p_melt <- melt(dat_p)
  colnames(dat_cor_melt) <- c("from","to","correlation")
  colnames(dat_p_melt) <- c("from","to","p")
  
  dat_cor_melt %>%
    left_join(dat_p_melt, by = c("from","to")) %>%
    dplyr::filter(p!="NA") %>%
    arrange(p) -> dat_net
  #adjust p values
  procs <- c("ABH")#select benjamini
  p.adjust <- mt.rawp2adjp(dat_net$p,procs)
  
  dat_net2 <- data.frame(dat_net,p.adjust$adjp) %>%
    dplyr::filter(abs(correlation)> cutoff_cor,ABH < cutoff_p) %>%
    dplyr::select(from:correlation,ABH)
  
  colnames(dat_net2) <- c("from","to","correlation","adjust.p")
  return(dat_net2)
}
########################################

bac1 <- read.table("otu_new_1.xls",row.names = 1,header=T,comment.char='',sep='\t')
bac1 <- as.data.frame(t(bac1))
C_edge_list <- calNetwork(bac1,cutoff_cor = 0.8,cutoff_p = 0.05,filter = FALSE)
C_edge_list$Weight <- abs(C_edge_list$correlation)
C_edge_list <- C_edge_list[,-4]
colnames(C_edge_list) <- c("Source", "Target","correlation","weight")

c( as.character(C_edge_list$from), as.character(C_edge_list$to)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> C_node_list
colnames(C_node_list) <- c("ID", "degree")

write.table(C_edge_list, 'C_edge_list.csv', row.names=F, sep = ',', quote = FALSE)
write.table(as.data.frame(C_node_list), 'C_node_list.csv', row.names=F, sep = ',', quote = FALSE)

