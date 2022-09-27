#if (!requireNamespace("devtools", quietly=TRUE))
  #install.packages("devtools")
#library(devtools)
#if (!requireNamespace("amplicon", quietly=TRUE))
  #install_github("microbiota/amplicon")
#suppressWarnings(suppressMessages(library(amplicon)))

setwd("G:/R/Rdata3")
library(amplicon)
# 读取特征ASV表，多样性分析输入抽平标准化的表
otu <- read.table("otu_table.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')

map<-read.table("mapping.xls",row.names = 1,header = T,sep='\t',comment.char='',check.names=F)

# 基于特征表、原数据、距离类型、分组列名、是否添加置信椭圆，是否添加样本标签
p = beta_cpcoa(otu, map, dis="bray", groupID= 'Treat', ellipse=T, label=T)+
  theme(text=element_text(size=15,face="bold"))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = 'black'))
p#最简单出图，无法细致修改
