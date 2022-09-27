setwd("G:/R/Rdata2")
library(tidyverse)
library(vegan)
library(linkET)
library(ggraph)
library(openxlsx)


otu1<- read.table("G:/R/Rdata2/otu_table2.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')
otu1 <- otu1 %>% select(-ncol(otu1)) %>% 
  t() %>% decostand("hellinger")#转置后数据标准化

otu3 <- read.table("G:/R/Rdata2/functional_table.xls",row.names = 1,header=T,comment.char='',sep='\t')
otu3 <- otu3 %>% 
  t() %>% decostand("hellinger")#转置后数据标准化

env.data <-read.xlsx('env_20220512.xlsx',sheet = 8,rowNames = T)
env <- env.data %>% 
  #select(SM, pH, SOC, TN, LC, RC, `LF-C`, `HF-C`) %>% 
  select(SM, pH, TN, LC, `LF-C`, MBC, MBN, MR, SIR) %>%
  log1p() %>% #数据log化
  na.omit()#去除NA

otu1<-otu1[match(row.names(env),row.names(otu1)),]
otu3<-otu3[match(row.names(env),row.names(otu3)),]
print(rownames(env) == rownames(otu1))
print(rownames(env) == rownames(otu3))

spe <- list(Bacteria=otu1,Function=otu3)
Bacteria <- dim(spe[[1]])[2]#dim()取出行列数
Function <- dim(spe[[2]])[2]

mantel <- mantel_test(as.data.frame(spe), env,
                      mantel_fun="mantel.partial",#偏mantel分析
                      #spec_dist = dist_func(method="bray"),
                      #env_dist = dist_func(method="euclidean"),
                      seed=123,
                      env_ctrl=TRUE,
                      spec_select = list(Bacteria=1:Bacteria,#分割数据，前面为细菌数据，后面为功能数据
                                         Function=(Bacteria+1):(Bacteria+Function)) )%>%
  mutate(rd= cut(abs(r), breaks = c(0, 0.1, 0.3,Inf),
                 labels = c("0 - 0.1", "0.1 - 0.3", "> 0.3"),right=FALSE),#定义Mantel的R值范围标签，便于出图
         pd = cut(p, breaks = c(-Inf, 0.01 ,0.05,Inf),
                  labels = c("< 0.01", "0.01 - 0.05","> 0.05"),right=FALSE),#定义Mantel检验的p值范围标签，便于出图
         col=cut(r,breaks=c(-Inf,0,Inf),
                 labels=c("r<0","r>0")))

write.table(mantel,"偏mantel0630_新.xls",sep="\t", quote=FALSE)

p <- qcorrplot(correlate(env), type = "upper",diag=FALSE) +#绘制理化数据热图
  geom_square() +#定义成方块状
  geom_mark(sep='\n',size=4,sig_thres=0.05,fontface="bold",family="serif")+
  geom_couple(aes(colour = pd, size = rd), #lty线型
              curvature=0.1,
              data = mantel,
              nudge_x = 0.2,
              label.size = 8,
              label.colour = "black",
              label.fontface = "bold") +#定义连线
  scale_linetype_manual(values=c("solid","dashed"))+#自定义线型
  scale_size_manual(values = c(0.5, 1.5,3))+
  #scale_colour_manual(values = c("#D95F02","grey"))+
  scale_colour_manual(values = c("#FF0000","#FF7F00","grey"))+
  #scale_linetype_manual(values = c("dotted", "solid")) +
  scale_fill_gradient2(midpoint=0,low="blue",mid="white",high="red",space='Lab')+
  #geom_diag_label(angle=45,size=2.5,color='black',fontface="bold")+
  #remove_y_axis()+
  #theme(axis.text.x=element_text(colour = "black",size = 8,face="bold"),
  #  axis.text.y=element_text(colour = "black",size = 8,face="bold")) +
  guides(size = guide_legend(title = "Mantel's r",#定义图例
                             override.aes = list(colour="grey35"),
                             order = 2),
         colour = guide_legend(title = "Mantle's P",
                               override.aes = list(size=3),
                               order = 1),
         fill = guide_colorbar(title = "Pearson's r", order = 3))+
  theme(text=element_text(size=16,  family="serif",face = "bold"),
        axis.text.x = element_text(size=16,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=16,  family="serif",colour = "black"))
p

ggsave("./偏Mantel test_0630_新.pdf", p, width = 200, height = 200, units = "mm")
