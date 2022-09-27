setwd("G:/R/Rdata2")
library(openxlsx)
library(ggplot2)
#otu<-read.xlsx('otu_table2_0.8.xlsx',sheet = 1,rowNames = T)
otu<-read.xlsx('otu_new.xlsx',sheet = 1,rowNames = T)
otu <- t(otu)#物种丰度矩阵，行是样本，列是物种变量
C <- otu[c("C4","C12","C17","C23","C30","C34"),]
DL <- otu[c("C1","C8","C14","C19","C27","C33"),]
NL <- otu[c("C2","C7","C15","C21","C25","C32"),]
NP <- otu[c("C3","C9","C13","C20","C26","C31"),]
NPDL <- otu[c("C6","C10","C16","C24","C29","C35"),]
NPNL <- otu[c("C5","C11","C18","C22","C28","C36"),]

#�? R 中，可使�? spaa 包的函数 niche.width() 计算生态位宽度指数，详情加�? spaa 包后 ?niche.width
library(spaa)
#这里�? Levins 生态位宽度指数为例，Shannon 生态位宽度指数可通过 method 参数修改
niche_width <- niche.width(otu, method = 'levins')
#niche_width

#画图展示分布
#boxplot(unlist(niche_width), ylab = 'niche breadth index')

#使用 EcolUtils 包函�? spec.gen() 实现 Specialist/Generalist species 的划�?
#详情加载 EcolUtils 包后 ?spec.gen
library(EcolUtils)

#niche.width.method = 'levins'，基�? Levins�?1968）的公式计算生态位宽度指数；若要计�? Shannon 生态位宽度指数可修改此参数
#n = 1000，随机化重排 1000 �?
#probs = c(0.025, 0.975)，计算双�? 95% 置信区间为准划分
set.seed(123)
spec_gen <- spec.gen(otu, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genC <- spec.gen(C, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genDL <- spec.gen(DL, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNL <- spec.gen(NL, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNP <- spec.gen(NP, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNPDL <- spec.gen(NPDL, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))
set.seed(123)
spec_genNPNL <- spec.gen(NPNL, niche.width.method = 'levins', perm.method = 'quasiswap', n = 1000, probs = c(0.025, 0.975))

tail(spec_gen)


#输出 Specialist/Generalist species 的划�?
write.table(spec_gen, 'spec_gen_0.822.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genC, 'spec_gen_0.8C22.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genDL, 'spec_gen_0.8DL22.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNL, 'spec_gen_0.8NL22.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNP, 'spec_gen_0.8NP22.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNPDL, 'spec_gen_0.8NPDL22.xls', sep = '\t', col.names = NA, quote = FALSE)
write.table(spec_genNPNL, 'spec_gen_0.8NPNL22.xls', sep = '\t', col.names = NA, quote = FALSE)

spec_gen$group <- rep("All",196)
spec_genC$group <- rep("C",196)
spec_genDL$group <- rep("DL",196)
spec_genNL$group <- rep("NL",196)
spec_genNP$group <- rep("NP",196)
spec_genNPDL$group <- rep("NPDL",196)
spec_genNPNL$group <- rep("NPNL",196)

all <- rbind(spec_gen,spec_genC,spec_genDL,spec_genNL,spec_genNP,spec_genNPDL,spec_genNPNL)
write.table(all, 'spec_gen_0.8all22.xls', sep = '\t', col.names = NA, quote = FALSE)

all<-read.table('spec_gen_0.822.xls',header = T,sep = '\t')
C<-read.table('spec_gen_0.8C22.xls',header = T,sep = '\t')
DL<-read.table('spec_gen_0.8DL22.xls',header = T,sep = '\t')
NL<-read.table('spec_gen_0.8NL22.xls',header = T,sep = '\t')
NP<-read.table('spec_gen_0.8NP22.xls',header = T,sep = '\t')
NPDL<-read.table('spec_gen_0.8NPDL22.xls',header = T,sep = '\t')
NPNL<-read.table('spec_gen_0.8NPNL22.xls',header = T,sep = '\t')

a <- unlist(all)#allΪdata.frame
sum(a=="SPECIALIST")
sum(a=="GENERALIST")
sum(a=="NON SIGNIFICANT")

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

shengtaiwei<-read.xlsx('shengtaiwei.xlsx',sheet = 4,rowNames = F)
library(reshape2)
data <- melt(shengtaiwei)

data$percent <- (data$value)/196
sum(data$percent)
#data <- data[-(1:3),]
yanse = colorRampPalette(colors = c("#984EA3","#BFBFBF","#A65628"))(3)#设置颜色�?

p <- ggplot(data,aes(x = group,y = percent,fill = species)) +
  geom_bar(position="fill",stat = "identity",width = 0.6) +
  scale_fill_manual(values = yanse) +
  labs(x='',y='Relative Abundance')+
  scale_x_discrete(limits = c("C","DL","NL",'NP','NPDL','NPNL'))+
  guides(fill=guide_legend(reverse = TRUE,title = "Classification"))+
  theme_bw()+
  theme(text=element_text(size=18,  family="serif",face = "bold"),
        axis.text.x = element_text(size=15,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=15,  family="serif",colour = "black"))+#serif在R中表示新
  scale_y_continuous(expand = c(0,0) )
p

box_all <- read.table('spec_gen_0.8all22.xls',header = T,sep = '\t')
box <- box_all[,c(2,7)]
box <- box[-(1:196),]
box$properties <- rep("Niche breadth",length(box$observed))#要修正离群�?
rownames(box) <-1: length(box$observed)#较大的离群值以去除离群值后的平均值代替，a[x,y] <- n
#矫正了DL/NP/NPNL/NPDL四组的较大离群值（2,1,1,1�?

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
    out <- LSD.test(fit,'sub_dat[, compare]',p.adj='BH')#进行了p值校�?
    #out$groups就可获取多重比较字母列表
    out$groups$type <- i
    out$groups$compare <- rownames(out$groups)
    
    a <- rbind(a,merge(out$means[,1:2], out$groups,by='sub_dat[, value]'))
  }
  names(a) <- c('mean','std','lsd',group,compare)
  return(a)
}

df <- ONE_LSD(box,'properties','group','observed')

#head(df)

p = ggplot(box)+
  geom_jitter(aes(x=group,y=observed,color=group),width=0.25,shape=20,size=2,show.legend = F)+#散点
  stat_boxplot(aes(x=group,y=observed),geom="errorbar",size=0.6,width=0.25)+
  geom_boxplot(aes(x=group,y=observed,fill=group))+#color非填充，fill填充
  #outlier.colour="red", outlier.shape=8, 
  #outlier.size=4,outlier.stroke =0.5,
  #outlier.alpha = 0.5,outlier.fill = "red")+
  #scale_y_continuous(expand = c(0.1,0))+
  geom_text(data=df,aes(x=group,y=1.85,label=lsd,vjust=1),size=10,family="serif")+
  #stat_summary(aes(x=group,y=observed),fun="mean",geom="point",shape=23,size=1,fill="white")+
  labs(x='',y='Niche breadth')+
  #ggprism::theme_prism()+
  #theme(axis.text.x = element_text(angle = 45))+
  #ggprism函数中的theme
  theme_bw()+
  guides(fill=guide_legend(title = "Treat"))+
  theme(text=element_text(size=22,  family="serif",face = "bold"),
        axis.text.x = element_text(size=18,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=15,  family="serif",colour = "black"))
p
