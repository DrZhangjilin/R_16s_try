setwd("G:/R/Rdata4_new/1zipi_new")
library(openxlsx)
library(ggplot2)
library(tidyverse)
library(reshape2)
#otu<-read.xlsx('otu_table2_0.8.xlsx',sheet = 1,rowNames = T)
otu<-read.table('otu_table2.xls',row.names = 1,skip = 1,header = T,sep ="\t",comment.char = "")
#otu <- t(otu)##物种丰度矩阵，行是样本，列是物种变量
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

#在 R 中，可使用 spaa 包的函数 niche.overlap() 计算生态位重叠指数
library(spaa)

#详情 ?niche.overlap
#有多种指数可供选择，通过 method 参数指定，这里以 levins 的生态位重叠指数为例
#niche_overlap <- niche.overlap(otu2, method = 'pianka')
#niche_overlap
#上述结果默认以 dist 类型存储，可转换为 matrix 类型输出到本地
#niche_overlap <- as.matrix(niche_overlap)
#write.table(niche_overlap, 'C_niche_overlap.csv', sep = ',', row.names = T,col.names = NA, quote = FALSE)

####生态位重叠计算
otu1[,colSums(otu1) > 0] %>%
  niche.overlap(method = "pianka") %>%
  data.matrix() %>%
  as_tibble(rownames = "otuid") ->bac_overlap

bac_overlap[upper.tri(bac_overlap,diag=FALSE)]=NA
bac_overlap %>%
  reshape2::melt(id.vars = "otuid") %>%
  filter(value!="NA") %>%
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> bac_overlap2

####基于重叠程度计算核心物种
#--------calculate niche overlap for each otu/species-------------
otu1[,colSums(otu1) > 0] %>%
  as_tibble() %>%
  colnames() -> bac_names

calculate_scores <- function(data, predict = bac_names) {
  x <- list()
  i = 0
  while (i < length(predict)) {
    i = i+1
    fam = predict[[i]]
    x[[i]] <- data %>%
      filter(str_detect(paired,fam)) %>%
      summarize(sum(value))
  }
  dt <- matrix(unlist(x),ncol = 1,nrow=length(predict))
  rownames(dt) <- predict
  colnames(dt) <- "Total_niche_overlap"
  dt %>%
    as_tibble(rownames = "otuid") %>%
    arrange(desc(Total_niche_overlap))
}

core1 <- calculate_scores(bac_overlap2)
##############################################################################################################
###########################################################################################################################

otu2[,colSums(otu2) > 0] %>%
  niche.overlap(method = "pianka") %>%
  data.matrix() %>%
  as_tibble(rownames = "otuid") ->bac_overlap

bac_overlap[upper.tri(bac_overlap,diag=FALSE)]=NA
bac_overlap %>%
  reshape2::melt(id.vars = "otuid") %>%
  filter(value!="NA") %>%
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> bac_overlap2

#--------calculate niche overlap for each otu/species-------------
otu2[,colSums(otu2) > 0] %>%
  as_tibble() %>%
  colnames() -> bac_names

core2 <- calculate_scores(bac_overlap2)
##############################################################################################################
###########################################################################################################################

otu3[,colSums(otu3) > 0] %>%
  niche.overlap(method = "pianka") %>%
  data.matrix() %>%
  as_tibble(rownames = "otuid") ->bac_overlap

bac_overlap[upper.tri(bac_overlap,diag=FALSE)]=NA
bac_overlap %>%
  reshape2::melt(id.vars = "otuid") %>%
  filter(value!="NA") %>%
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> bac_overlap2


#--------calculate niche overlap for each otu/species-------------
otu3[,colSums(otu3) > 0] %>%
  as_tibble() %>%
  colnames() -> bac_names
core3 <- calculate_scores(bac_overlap2)

##############################################################################################################
###########################################################################################################################

otu4[,colSums(otu4) > 0] %>%
  niche.overlap(method = "pianka") %>%
  data.matrix() %>%
  as_tibble(rownames = "otuid") ->bac_overlap

bac_overlap[upper.tri(bac_overlap,diag=FALSE)]=NA
bac_overlap %>%
  reshape2::melt(id.vars = "otuid") %>%
  filter(value!="NA") %>%
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> bac_overlap2


#--------calculate niche overlap for each otu/species-------------
otu4[,colSums(otu4) > 0] %>%
  as_tibble() %>%
  colnames() -> bac_names
core4 <- calculate_scores(bac_overlap2)
##############################################################################################################
###########################################################################################################################

otu5[,colSums(otu5) > 0] %>%
  niche.overlap(method = "pianka") %>%
  data.matrix() %>%
  as_tibble(rownames = "otuid") ->bac_overlap

bac_overlap[upper.tri(bac_overlap,diag=FALSE)]=NA
bac_overlap %>%
  reshape2::melt(id.vars = "otuid") %>%
  filter(value!="NA") %>%
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> bac_overlap2

#--------calculate niche overlap for each otu/species-------------
otu5[,colSums(otu5) > 0] %>%
  as_tibble() %>%
  colnames() -> bac_names
core5 <- calculate_scores(bac_overlap2)

##############################################################################################################
###########################################################################################################################

otu6[,colSums(otu6) > 0] %>%
  niche.overlap(method = "pianka") %>%
  data.matrix() %>%
  as_tibble(rownames = "otuid") ->bac_overlap

bac_overlap[upper.tri(bac_overlap,diag=FALSE)]=NA
bac_overlap %>%
  reshape2::melt(id.vars = "otuid") %>%
  filter(value!="NA") %>%
  mutate(paired = str_c(otuid,variable,sep = "_") ) -> bac_overlap2

#--------calculate niche overlap for each otu/species-------------
otu6[,colSums(otu6) > 0] %>%
  as_tibble() %>%
  colnames() -> bac_names
core6 <- calculate_scores(bac_overlap2)

################################################################################
################################################################################
################################################################################
total <- rbind(core1,core2,core3,core4,core5,core6)
total$Treat <- c(rep('C', 156),
               rep('DL', 129),
               rep('NL', 150),
               rep('NP', 183),
               rep('NPDL',156),
               rep('NPNL',152))
write.table(total, 'total_niche_overlap.xls', sep = '\t', col.names = NA, quote = FALSE)
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

dat <- melt(total)##宽数据变长数据
df <- ONE_LSD(dat,'variable','Treat','value')
head(df)

p = ggplot(dat)+
  geom_jitter(aes(x=Treat,y=value,color=Treat),width=0.25,shape=20,size=2,show.legend = F)+#鏁ｇ偣
  stat_boxplot(aes(x=Treat,y=value),geom="errorbar",size=0.6,width=0.25)+
  geom_boxplot(aes(x=Treat,y=value,fill=Treat))+
  #outlier.colour="red", outlier.shape=8, 
  #outlier.size=4,outlier.stroke =0.5,
  #outlier.alpha = 0.5,outlier.fill = "red")+
  #scale_y_continuous(expand = c(0.1,0))+
  geom_text(data=df,aes(x=Treat,y=200,label=lsd,vjust=1),size=10,family="serif")+
  #stat_summary(aes(x=group,y=observed),fun="mean",geom="point",shape=23,size=1,fill="white")+
  labs(x='',y='Niche overlap')+
  #ggprism::theme_prism()+
  #theme(axis.text.x = element_text(angle = 45))+
  theme_bw()+
  guides(fill=guide_legend(title = "Treat"))+
  theme(text=element_text(size=48,  family="serif",face = "bold"),
        axis.text.x = element_text(size=36,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=36,  family="serif",colour = "black"))
p
ggsave("./生态位重叠2.pdf", p, width = 350, height = 400, units = "mm")
