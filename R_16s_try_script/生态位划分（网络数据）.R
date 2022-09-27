setwd("G:/R/Rdata4_new/1zipi_new")
library(openxlsx)
library(ggplot2)
#otu<-read.xlsx('otu_table2_0.8.xlsx',sheet = 1,rowNames = T)
otu<-read.table('otu_table2.xls',row.names = 1,skip = 1,header = T,sep ="\t",comment.char = "")
#otu <- t(otu)#物种丰度矩阵，行是样本，列是物种变量
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
#在 R 中，可使用 spaa 包的函数 niche.width() 计算生态位宽度指数，详情加载 spaa 包后 ?niche.width
#library(spaa)
#这里以 Levins 生态位宽度指数为例，Shannon 生态位宽度指数可通过 method 参数修改
#niche_width <- niche.width(otu, method = 'levins')
#niche_width

#画图展示分布
#boxplot(unlist(niche_width), ylab = 'niche breadth index')

#使用 EcolUtils 包函数 spec.gen() 实现 Specialist/Generalist species 的划分
#详情加载 EcolUtils 包后 ?spec.gen
library(EcolUtils)

#niche.width.method = 'levins'，基于 Levins（1968）的公式计算生态位宽度指数；若要计算 Shannon 生态位宽度指数可修改此参数
#n = 1000，随机化重排 1000 次
#probs = c(0.025, 0.975)，计算双侧 95% 置信区间为准划分
#set.seed(123)
#spec_gen <- spec.gen(otu, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genC <- spec.gen(otu1, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genDL <- spec.gen(otu2, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNL <- spec.gen(otu3, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNP <- spec.gen(otu4, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNPDL <- spec.gen(otu5, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNPNL <- spec.gen(otu6, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))

#tail(spec_gen)


#输出 Specialist/Generalist species 的划分
#write.table(spec_gen, 'spec_gen_0.822.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genC, 'spec_gen_C.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genDL, 'spec_gen_DL.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNL, 'spec_gen_NL.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNP, 'spec_gen_NP.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNPDL, 'spec_gen_NPDL.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNPNL, 'spec_gen_NPNL.xls', sep = '\t', col.names = NA, quote = FALSE)

#spec_gen$group <- rep("All",196)
spec_genC$group <- rep("C",156)
spec_genDL$group <- rep("DL",129)
spec_genNL$group <- rep("NL",150)
spec_genNP$group <- rep("NP",183)
spec_genNPDL$group <- rep("NPDL",156)
spec_genNPNL$group <- rep("NPNL",152)

#all <- rbind(spec_gen,spec_genC,spec_genDL,spec_genNL,spec_genNP,spec_genNPDL,spec_genNPNL)
#write.table(all, 'spec_gen_0.8all22.xls', sep = '\t', col.names = NA, quote = FALSE)

#all<-read.table('spec_gen_0.822.xls',header = T,sep = '\t')
C<-read.table('spec_gen_C.xls',header = T,sep = '\t')
DL<-read.table('spec_gen_DL.xls',header = T,sep = '\t')
NL<-read.table('spec_gen_NL.xls',header = T,sep = '\t')
NP<-read.table('spec_gen_NP.xls',header = T,sep = '\t')
NPDL<-read.table('spec_gen_NPDL.xls',header = T,sep = '\t')
NPNL<-read.table('spec_gen_NPNL.xls',header = T,sep = '\t')

#a <- unlist(all)
#sum(a=="SPECIALIST")
#sum(a=="GENERALIST")
#sum(a=="NON SIGNIFICANT")

b <- unlist(C)
sum(b=="SPECIALIST")
sum(b=="GENERALIST")
sum(b=="NON SIGNIFICANT")

c <- unlist(DL)
sum(c=="SPECIALIST")
sum(c=="GENERALIST")
sum(c=="NON SIGNIFICANT")

d <- unlist(NL)
sum(d=="SPECIALIST")
sum(d=="GENERALIST")
sum(d=="NON SIGNIFICANT")

e <- unlist(NP)
sum(e=="SPECIALIST")
sum(e=="GENERALIST")
sum(e=="NON SIGNIFICANT")

f <- unlist(NPDL)
sum(f=="SPECIALIST")
sum(f=="GENERALIST")
sum(f=="NON SIGNIFICANT")

g <- unlist(NPNL)
sum(g=="SPECIALIST")
sum(g=="GENERALIST")
sum(g=="NON SIGNIFICANT")

shengtaiwei<-read.xlsx('shengtaiwei.xlsx',sheet = 1,rowNames = F)
library(reshape2)
shengtaiwei <- shengtaiwei[,-1]
data <- melt(shengtaiwei)

#data$percent <- (data$value)/196
#sum(data$percent)
sum(data$value)
#data <- data[-(1:3),]
yanse = colorRampPalette(colors = c("#984EA3","#BFBFBF","#A65628"))(3)#设置颜色，

p1 <- ggplot(data,aes(x = group,y = value,fill = species)) +
  geom_bar(position="fill",stat = "identity",width = 0.6) +
  scale_fill_manual(values = yanse) +
  labs(x='',y='Relative Abundance')+
  scale_x_discrete(limits = c("C","DL","NL",'NP','NPDL','NPNL'))+
  guides(fill=guide_legend(reverse = TRUE,title = "Classification"))+
  theme_bw()+
  theme(text=element_text(size=28,  family="serif",face = "bold"),
        axis.text.x = element_text(size=20,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=20,  family="serif",colour = "black"))+#serif在R中表示新
  scale_y_continuous(expand = c(0,0) )
p1
#ggsave("./生态位分化.pdf", p1, width = 350, height = 400, units = "mm")

C1 <- C[,1]
DL1 <- DL[,1]
NL1 <- NL[,1]
NP1 <- NP[,1]
NPDL1 <- NPDL[,1]
NPNL1 <- NPNL[,1]

C2 <- C[,2]
DL2 <- DL[,2]
NL2 <- NL[,2]
NP2 <- NP[,2]
NPDL2 <- NPDL[,2]
NPNL2 <- NPNL[,2]

C2 <- C[,3]
DL2 <- DL[,3]
NL2 <- NL[,3]
NP2 <- NP[,3]
NPDL2 <- NPDL[,3]
NPNL2 <- NPNL[,3]

dat1 <- data.frame(Names = c(C1,DL1,NL1,NP1,NPDL1,NPNL1),
                   Niche.breadth = c(C2,DL2,NL2,NP2,NPDL2,NPNL2),
                   group = c(rep('C', 156),
                               rep('DL', 129),
                               rep('NL', 150),
                               rep('NP', 183),
                               rep('NPDL',156),
                               rep('NPNL',152)))
dat1$properties <- rep("Niche breadth",length(dat1$Niche.breadth))#要修正离群值
#rownames(box) <-1: length(box$observed)#较大的离群值以去除离群值后的平均值代替，a[x,y] <- n

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

df <- ONE_LSD(dat1,'properties','group','Niche.breadth')

head(df)

p = ggplot(data = dat1, aes(x=group, y=Niche.breadth)) + 
  geom_violin(aes(fill = group)) +#小提琴图
  #stat_boxplot(geom="errorbar",size=0.6,width=0.2)+
  geom_boxplot(width = 0.2, fill = "white", notchwidth=1)+
  #stat_summary(fun=mean, geom="point", shape=1, size=3, color="red") +
  #geom_jitter(width=0.25,shape=20,size=2.5)+#散点
  stat_summary(fun="mean",geom="point",shape=23,size=2,fill="red")+
  #ggplot(dat1)+
  #geom_jitter(aes(x=group,y=Niche.breadth,color=group),width=0.25,shape=20,size=2,show.legend = F)+#散点
 # stat_boxplot(aes(x=group,y=Niche.breadth),geom="errorbar",size=0.6,width=0.25)+
 # geom_boxplot(aes(x=group,y=Niche.breadth,fill=group))+#color非填充，fill填充
  #outlier.colour="red", outlier.shape=8, 
  #outlier.size=4,outlier.stroke =0.5,
  #outlier.alpha = 0.5,outlier.fill = "red")+
  #scale_y_continuous(expand = c(0.1,0))+
  geom_text(data=df,aes(x=group,y=6,label=lsd,vjust=1),size=10,family="serif")+
  #stat_summary(aes(x=group,y=observed),fun="mean",geom="point",shape=23,size=1,fill="white")+
  labs(x='',y='Niche breadth')+
  #ggprism::theme_prism()+
  #theme(axis.text.x = element_text(angle = 45))+
  #ggprism函数中的theme
  theme_bw()+
  guides(fill=guide_legend(title = "Treat"))+
  theme(text=element_text(size=28,  family="serif",face = "bold"),
        axis.text.x = element_text(size=20,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=20,  family="serif",colour = "black"))
p
#ggsave("./Niche breadth.pdf", p, width = 400, height = 400, units = "mm")
