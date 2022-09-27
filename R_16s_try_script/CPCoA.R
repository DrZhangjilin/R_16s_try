#if (!requireNamespace("devtools", quietly=TRUE))
  #install.packages("devtools")
#library(devtools)
#if (!requireNamespace("amplicon", quietly=TRUE))
  #install_github("microbiota/amplicon")
#suppressWarnings(suppressMessages(library(amplicon)))

setwd("G:/R/Rdata3")
library(amplicon)
# ��ȡ����ASV���������Է��������ƽ��׼���ı�
otu <- read.table("otu_table.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')

map<-read.table("mapping.xls",row.names = 1,header = T,sep='\t',comment.char='',check.names=F)

# ������������ԭ���ݡ��������͡������������Ƿ�����������Բ���Ƿ�����������ǩ
p = beta_cpcoa(otu, map, dis="bray", groupID= 'Treat', ellipse=T, label=T)+
  theme(text=element_text(size=15,face="bold"))+
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'transparent', color = 'black'))
p#��򵥳�ͼ���޷�ϸ���޸�