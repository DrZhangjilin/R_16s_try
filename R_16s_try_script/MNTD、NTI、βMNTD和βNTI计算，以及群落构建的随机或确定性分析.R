setwd("G:/R/Rdata2")
#读取 OTU 丰度表
#otu <- read.table("otu99.xls",header=T,skip=1,row.names=1,comment.char='',sep='\t')
#otu <- data.frame(t(otu))

#读取 OTU 的进化树文件，需借助 ape 包
#library(ape)

#tree <- read.tree('rep_phylo.tre')

##使用 picante 包计算 MNTD 和 NTI
#library(picante)

#修剪系统发育树以仅包含群落数据集中存在的物种或具有非缺失性状数据的物种
#tree <- prune.sample(otu, tree)

#从进化树中获取 OTU 之间的系统发育距离
#dis <- cophenetic(tree)

#计算 MNTD，生成 MNTD 的零分布以及获取 NTI，详情 ?ses.mntd
#null.model 用于指定零模型方法；runs 指定获得 MNTD 零分布的置换次数，本示例以 999 次随机为例
#set.seed(123)
#mntd <- ses.mntd(otu, dis, abundance.weighted = TRUE, null.model = 'taxa.labels', runs = 999)
#mntd

#对 mntd.obs.z 取负值就是 NTI
#mntd$NTI <- mntd$mntd.obs.z * -1
#mntd

#输出
#write.csv(mntd, 'mntd.csv', quote = FALSE)

library(NST)
library(picante)
library(ape)

#读取 OTU 丰度表
otu <- read.table("otu_table2_0.8.xls",header=T,row.names=1,comment.char='',sep='\t',skip=1)
otu <- otu[,-ncol(otu)]
otu <- data.frame(t(otu))


#读取样本分组文件
group <- read.table('map1.xls', row.names = 1,header=T,sep='\t')

#读取 OTU 的进化树文件
tree <- read.tree('rep_phylo.tre')

#修剪系统发育树以仅包含群落数据集中存在的物种或具有非缺失性状数据的物种
tree <- prune.sample(otu, tree)

#函数 pNST() 本是用来计算 pNST 指数的，通过添加参数 SES=TRUE 即可同时计算 betaMNTD 和 betaNTI，详情 ?pNST
#phylo.shuffle=TRUE 意为在系统发育树中随机置换以随机化物种之间的系统发育关系
#taxo.null.model，随机化丰度矩阵的零模型算法，默认为 NULL 意为只进行系统发育洗牌。若有需要指定具体的零模型，请更改此处
#rand 指定获得 betaMNTD 零分布的置换次数，本示例以 1000 次随机为例
#nworker 指定多线程，本示例 4 线程运行
#此外，pNST() 还可通过参数 RC=TRUE 计算 Raup-Crick（不知此时是否需要修改 taxo.null.model？），但是本示例没有再执行它，细节部分还请参阅函数帮助
set.seed(123)
pnst <- pNST(comm = otu, tree = tree, group = group, phylo.shuffle = TRUE, taxo.null.model = NULL,
             pd.wd = tempdir(), abundance.weighted = TRUE, rand = 1000, nworker = 6, SES = TRUE, RC = FALSE)


#本篇主要是计算 betaMNTD 和 betaNTI，暂且不管其它的（比方说 pNST），因此在结果项里面，主要看这个就行了
#names(pnst)
betaMNTD <- pnst$index.pair
head(betaMNTD)

#输出两两样本的 betaMNTD 和 betaNTI
write.csv(betaMNTD, 'G:/R/Rdata4_new/betaMNTD_0.8.csv', quote = FALSE, row.names= FALSE)

#其它的，如计算随机或确定因素的贡献率
#提取不同分组内的样本对
###################################################
C <- rownames(subset(group, Treat=='C'))
betaMNTD_C <- subset(betaMNTD, name1 %in% C & name2 %in% C)
DL <- rownames(subset(group, Treat=='DL'))
betaMNTD_DL <- subset(betaMNTD, name1 %in% DL & name2 %in% DL)
NL <- rownames(subset(group, Treat=='NL'))
betaMNTD_NL <- subset(betaMNTD, name1 %in% NL & name2 %in% NL)
NP <- rownames(subset(group, Treat=='NP'))
betaMNTD_NP <- subset(betaMNTD, name1 %in% NP & name2 %in% NP)
NPDL <- rownames(subset(group, Treat=='NPDL'))
betaMNTD_NPDL <- subset(betaMNTD, name1 %in% NPDL & name2 %in% NPDL)
NPNL <- rownames(subset(group, Treat=='NPNL'))
betaMNTD_NPNL <- subset(betaMNTD, name1 %in% NPNL & name2 %in% NPNL)


#|βNTI|<2的数量/总数量=随机因素贡献率
RC <- nrow(betaMNTD_C[which(abs(betaMNTD_C$bNTI.wt)<2), ])/nrow(betaMNTD_C)  #C 组
RDL <- nrow(betaMNTD_DL[which(abs(betaMNTD_DL$bNTI.wt)<2), ])/nrow(betaMNTD_DL)  #DL 组
RNL <- nrow(betaMNTD_NL[which(abs(betaMNTD_NL$bNTI.wt)<2), ])/nrow(betaMNTD_NL)
RNP <- nrow(betaMNTD_NP[which(abs(betaMNTD_NP$bNTI.wt)<2), ])/nrow(betaMNTD_NP)
RNPDL <- nrow(betaMNTD_NPDL[which(abs(betaMNTD_NPDL$bNTI.wt)<2), ])/nrow(betaMNTD_NPDL)
RNPNL <- nrow(betaMNTD_NPNL[which(abs(betaMNTD_NPNL$bNTI.wt)<2), ])/nrow(betaMNTD_NPNL)

#|βNTI|>2的数量/总数量=确定因素贡献率
DC <- nrow(betaMNTD_C[which(abs(betaMNTD_C$bNTI.wt)>2), ])/nrow(betaMNTD_C)  #C 组
DDL <- nrow(betaMNTD_DL[which(abs(betaMNTD_DL$bNTI.wt)>2), ])/nrow(betaMNTD_DL)  #DL 组
DNL <- nrow(betaMNTD_NL[which(abs(betaMNTD_NL$bNTI.wt)>2), ])/nrow(betaMNTD_NL)
DNP <- nrow(betaMNTD_NP[which(abs(betaMNTD_NP$bNTI.wt)>2), ])/nrow(betaMNTD_NP)
DNPDL <- nrow(betaMNTD_NPDL[which(abs(betaMNTD_NPDL$bNTI.wt)>2), ])/nrow(betaMNTD_NPDL)
DNPNL <- nrow(betaMNTD_NPNL[which(abs(betaMNTD_NPNL$bNTI.wt)>2), ])/nrow(betaMNTD_NPNL)

a <- data.frame(RC,RDL,RNL,RNP,RNPDL,RNPNL,DC,DDL,DNL,DNP,DNPDL,DNPNL)
write.csv(a,"G:/R/Rdata4_new/随机、确定贡献率（0.8）,csv",quote=F,row.names=F)

#再如比较两组 betaNTI 是否存在显著差异（这里以非参数的 wilcox test 为例）
library(ggpubr)
betaMNTD_C$group <- 'C'
betaMNTD_DL$group <- 'DL'
betaMNTD_NL$group <- 'NL'
betaMNTD_NP$group <- 'NP'
betaMNTD_NPDL$group <- 'NPDL'
betaMNTD_NPNL$group <- 'NPNL'
betaMNTD_group <- rbind(betaMNTD_C, betaMNTD_DL,betaMNTD_NL,
                        betaMNTD_NP,betaMNTD_NPDL,betaMNTD_NPNL)
betaMNTD_group$group <- as.factor(betaMNTD_group$group)#将node_properties列转化成因子

p <- ggboxplot(data = betaMNTD_group, x = 'group', y = 'bNTI.wt', fill = 'group') +
  stat_boxplot(geom="errorbar",size=0.6,width=0.2)+
  geom_boxplot(outlier.colour="red", 
               outlier.shape=8,
               outlier.size=4,
               outlier.stroke =0.5,
               outlier.alpha = 0.5,
               outlier.fill = "red",aes(fill=group))+
  geom_jitter(width=0.25,shape=20,size=2.5)+stat_summary(fun="mean",geom="point",shape=23,size=3,fill="white")+
  #stat_compare_means(method = 'wilcox.test',
                     #comparisons = list(c('C', 'DL'),c('C','NL'),c('C','NP')))+
  labs(x='',y = 'betaNTI')+
  scale_y_continuous(breaks = seq(-2,7,2))+#限定纵轴的范围
  stat_compare_means(method = "anova",label.y=7)+ # Add global p-value
  
  geom_hline(yintercept =c(2,-2), color = 'DimGrey', size = 0.5,linetype="dashed")+
  guides(fill=guide_legend(title = "Treat"))+
  theme_bw()+
  theme(text=element_text(size=22,  family="serif",face = "bold"),
        axis.text.x = element_text(size=18,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=15,  family="serif",colour = "black"))

p

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

betaMNTD_group$factor <-  c(rep('betaNTI', 90))
df <- ONE_LSD(betaMNTD_group,'factor','group','bNTI.wt')
head(df)

write.table(betaMNTD_group,"G:/R/Rdata4_new/betaMNTD_group_0.8.xls",sep="\t", quote=FALSE,col.names = NA)

p1 <- ggplot(data = betaMNTD_group, aes(x=group, y=bNTI.wt)) + 
  geom_violin(aes(fill = group)) +#小提琴图
  #stat_boxplot(geom="errorbar",size=0.6,width=0.2)+
  geom_boxplot(width = 0.2, fill = "white", notchwidth=1)+
  #stat_summary(fun=mean, geom="point", shape=1, size=3, color="red") +
  #geom_jitter(width=0.25,shape=20,size=2.5)+#散点
  stat_summary(fun="mean",geom="point",shape=23,size=2,fill="red")+
  labs(x='',y = 'betaNTI')+
  scale_y_continuous(breaks = seq(-2,7,2))+#限定纵轴的范围
  stat_compare_means(method = "anova",label.y=-1, family="serif",size=8)+ # Add global p-value
  geom_hline(yintercept =c(2,-2), color = 'DimGrey', size = 0.5,linetype="dashed")+
  geom_text(data=df,aes(x=group,y=6.5,label=lsd,),size=8,family="serif")+
  guides(fill=guide_legend(title = "Treat"))+
  theme_bw()+
  theme(text=element_text(size=22,  family="serif",face = "bold"),
        axis.text.x = element_text(size=18,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=15,  family="serif",colour = "black"))

p1
