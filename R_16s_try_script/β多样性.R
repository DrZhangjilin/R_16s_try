library(openxlsx)
library(ggpubr)

otu<-read.xlsx('otu_table2_0.8.xlsx',sheet = 1,rowNames = T)
#otu<-otu[,-ncol(otu)]
otu<-t(otu)

#计算与 Beta 多样性有关的群落相异指数，例如使用 vegan 包计算 Bray-curtis 距离，详情加载 vegan 包后 ?vegdist
dis <- vegan::vegdist(otu, method = 'bray')

#以矩阵形式输出
dis <- as.matrix(dis)
write.table(dis, 'Bray-curtis_0.8.txt', sep = '\t', col.names = NA, quote = FALSE)
#这里直接通过上述提供的物种丰度表计算获得Bray-Curtis相异矩阵。
#如果您希望计算其它类型的相似性或相异指数，如微生物群落中另一常见的Unifrac相异指数，也是可以的，具体使用哪种相似性或相异指数根据实际情况选择。
#此外，当然如果您在一开始就准备好了这个样本间成对的相似性或相异矩阵文件（而不再通过分类群的丰度表从头计算），也可以直接使用它。

#读取 Bray-curtis 距离矩阵

dis <- read.delim('Bray-curtis_0.8.txt', row.names = 1)




dis <- read.delim('G:/R/Rdata2/weighted_unifrac_biaozhun.txt', row.names = 1)
#setwd("G:/R/Rdata2")
#读取样本分组信息
group <- read.table('map1.xls', header=T,row.names=1,comment.char='',sep='\t')
rownames(group)->group$samples

##例如，比较 Env1、Env2、Env3 三组之间，群落的 Beta 多样性差异
#根据分组获得组内距离矩阵
env1 <- subset(group, Treat == 'C')$samples
dis_env1 <- dis[env1,env1]

env2 <- subset(group, Treat == 'LA')$samples
dis_env2 <- dis[env2,env2]

env3 <- subset(group, Treat == 'LR')$samples
dis_env3 <- dis[env3,env3]

env4 <- subset(group, Treat == 'PR')$samples
dis_env4 <- dis[env4,env4]

env5 <- subset(group, Treat == 'LAPR')$samples
dis_env5 <- dis[env5,env5]

env6 <- subset(group, Treat == 'LRPR')$samples
dis_env6 <- dis[env6,env6]


#将矩阵转化为向量，以便用于作图和统计
dis_env1 <- as.vector(as.dist(dis_env1))
dis_env2 <- as.vector(as.dist(dis_env2))
dis_env3 <- as.vector(as.dist(dis_env3))
dis_env4 <- as.vector(as.dist(dis_env4))
dis_env5 <- as.vector(as.dist(dis_env5))
dis_env6 <- as.vector(as.dist(dis_env6))
#构建作图数据集
dat <- data.frame(
  dis = c(dis_env1, dis_env2, dis_env3, dis_env4,dis_env5,dis_env6),
  group = factor(c(
    rep('C', length(dis_env1)),
    rep('LA', length(dis_env2)),
    rep('LR', length(dis_env3)),
    rep('PR', length(dis_env4)),
    rep('LAPR', length(dis_env5)),
    rep('LRPR', length(dis_env6))
  ), levels = c('C', 'LA', 'LR', 'PR','LAPR','LRPR'))
)
write.table(dat, 'unifrac_diversity_all_biaozhun.xls', sep = '\t', col.names = NA, quote = FALSE)

#使用 ggplot2 绘制各组内 Bray-curtis 距离指数分布的箱线图
library(ggplot2)
yanse <-c("#999999","#F781BF","#A65628","#FFFF33","#FF7F00","#984EA3",
          "#377EB8","#74D944","#E41A1C","#DA5724","#CE50CA",
          "#D3D93E","#C0717C","#CBD588","#D7C1B1","#5F7FC7","#673770",
          "#3F4921","#CD9BCD","#38333E","#689030","#AD6F3B","#4DAF4A")
#p <- ggplot(dat, aes(group, dis)) +
 # geom_boxplot(aes(fill = group), width = 0.6) +
  #scale_fill_manual(values =yanse) +
 # theme(panel.grid = element_blank(), panel.background = element_blank(),
 #       axis.line = element_line(colour = 'black'), legend.position = 'none') +
 # labs(x = NULL, y = 'Unifrac dissimilarity\n')+
 # theme_bw()+
 # guides(fill=guide_legend(title = "Treat"))+
 # theme(text=element_text(size=20,  family="serif",face = "bold"),
  #      axis.text.x = element_text(size=15,  family="serif",colour = "black"))#serif在R中表示新罗马字体
p <- ggplot(dat)+
  geom_jitter(aes(x=group,y=dis,color=group),width=0.25,shape=20,size=2,show.legend = F)+#散点
  stat_boxplot(aes(x=group,y=dis),geom="errorbar",size=0.6,width=0.25)+
  geom_boxplot(aes(x=group,y=dis,fill=group))+#color非填充，fill填充
  #stat_summary(aes(x=group,y=dis),fun="mean",geom="point",shape=23,size=1,fill="white")+
  #scale_fill_manual(values =group)+
  #outlier.colour="red", outlier.shape=8,
  #outlier.size=4,outlier.stroke =0.5,
  #outlier.alpha = 0.5,outlier.fill = "red")+
  #facet_wrap(.~node_properties,scales = "free_y")+#分面
  #scale_y_continuous(expand = c(0.1,0))+
  #geom_text(data=df,aes(x=treat,y=0,label=lsd,vjust=1),size=10,family="serif")+
  #stat_summary(aes(x=treat,y=v),fun="mean",geom="point",shape=23,size=1,fill="white")+
  labs(x='',y="UniFrac dissimilarity\n")+
  #ggprism::theme_prism()+
  #theme(axis.text.x = element_text(angle = 45))+
  #ggprism函数中的theme
  #stat_compare_means( aes(x=group, y=dis),method = "kruskal.test",label.x="DL",label.y=0.14, family="serif",size=8)+ # Add global p-value
  theme_bw()+
  #guides(fill=guide_legend(title = "Treat"))+
  guides(fill="none")+
  theme(text=element_text(size=10,  family="serif",face = "bold"),
        axis.text.x = element_text(size=8,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=8,  family="serif",colour = "black"))
p
ggsave('./系统发育多样性.pdf',p,width=100,height=85,units="mm")

p1 <- ggplot(data = dat, aes(x=group, y=dis)) +
  geom_violin(aes(fill = group)) +#小提琴图
  #stat_boxplot(geom="errorbar",size=0.6,width=0.2)+
  geom_boxplot(width = 0.2, fill = "white", notchwidth=1)+
  scale_fill_manual(values =yanse)+
  #stat_summary(fun=mean, geom="point", shape=1, size=3, color="red") +
  #geom_jitter(width=0.25,shape=20,size=2.5)+#散点
  stat_summary(fun="mean",geom="point",shape=23,size=2,fill="red")+
  labs(x='',y = 'Unifrac dissimilarity\n')+
  #scale_y_continuous(breaks = seq(-2,7,2))+#限定纵轴的范围
  stat_compare_means(method = "kruskal.test",label.y=0.14, family="serif",size=8)+ # Add global p-value
  #stat_compare_means(label = "p.signif", method = "wilcox.test",comparisons = list(
  #  c('C', 'NP')
   # c('C', 'DL')
    #c('C', 'WN'),
    #c('W', 'N'),
    #c('W', 'WN'),
   # c('C', 'NPNL')
 # ))+
  #geom_text(data=df,aes(x=group,y=6.5,label=lsd,),size=8,family="serif")+
  guides(fill=guide_legend(title = "Treat"))+
  theme_bw()+
  theme(text=element_text(size=22,  family="serif",face = "bold"),
        axis.text.x = element_text(size=18,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=15,  family="serif",colour = "black"))

p1

#四组的整体差异分析，使用 Kruskal-Wallis Test 执行，详情 ?kruskal.test
kruskal.test(dis~group, data = dat)

#如果整体显著再进行两两分组的比较，使用 Wilcoxon 秩和检验执行双侧检验，详情 ?wilcox.test
wilcox.test(dis_env1, dis_env2, alternative = 'two.sided')
wilcox.test(dis_env1, dis_env3, alternative = 'two.sided')
wilcox.test(dis_env1, dis_env4, alternative = 'two.sided')
wilcox.test(dis_env1, dis_env5, alternative = 'two.sided')
wilcox.test(dis_env1, dis_env6, alternative = 'two.sided')


#考虑到 Wilcoxon 秩和检验体现了中位数的差异，因此计算三组数据的中位数以评估 Beta 多样性的高低水平
median(dis_env1)
median(dis_env2)
median(dis_env3)
median(dis_env4)

