setwd("G:/R/Rdata2")
#setwd("G:/R/Rdata3")
#导入作图包
library(vegan)
library(ggplot2)
library(openxlsx)
library(ggpubr)


#读取otu数据文件
otu <- read.table("otu_table2.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')
otu <- otu[,-ncol(otu)]
otu1<-data.frame(t(otu))
#根据物种组成计算样方距离，结果以 dist 数据类型存储
#bray_dis <- vegdist(otu1, method = 'bray')
bray_dis <- vegdist(otu1, method = 'euclidean')
#样方排序坐标
pcoa <- cmdscale(bray_dis, k =2, eig = TRUE)
site <- data.frame(pcoa$point)[1:2]
site$name <- rownames(site)

#读取分组文件
mapping <- read.table("mapping.xls",row.names = 1,header=T,comment.char='',sep='\t')
#site$group <- c(rep('A', 4), rep('B', 4), rep('C', 4), rep('D', 4))
#将分组文件和数据文件以行名合并
merged=merge(site,mapping,by="row.names",all.x=TRUE)
#使用 wascores() 被动添加物种得分（坐标），丰度加权平均方法
species <- wascores(pcoa$points, otu1)
species<-data.frame(species)

pcoa_exp <- pcoa$eig/sum(pcoa$eig)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 2), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

p <- ggplot(data = merged, aes(X1, X2, shape = Treat, color = Treat)) +
  theme_classic()+
  geom_point(aes(color = Treat),size=2) +
  stat_ellipse(aes(fill = Treat), geom = 'polygon', level = 0.68, alpha = 0.2, show.legend = FALSE) +   #添加置信椭圆，注意不是聚类
  #scale_shape_manual(values = c(3, 8, 10, 19)) +
  #scale_color_manual(values =c('red3', 'orange3', 'green3','blue')) + #点的颜色
  #scale_fill_manual(values =c('red3', 'orange3', 'green3','blue')) + #椭圆填充的颜色
  #theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'),
  #前者是图内背景线的颜色；后者是框的颜色 
  #plot.title = element_text(hjust = 0.5)) + #hjust指总title在横线上的位置
  geom_vline(xintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #过坐标零点的两条线
  geom_hline(yintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #过坐标零点的两条线
  #geom_text(data = merged, aes(label = name), color = "darkorchid", size = 3) +    #丰度物种标签
  labs(x = pcoa1, y = pcoa2)  #大标题title = 'PCoA丰度物种'
p
#定义经典背景
p1 <- p+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'))
p1
