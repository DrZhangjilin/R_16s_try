setwd("G:/R/Rdata2")
setwd("G:/R/Rdata3")
library(ggplot2)
library(vegan)
library(openxlsx)

otu <- read.table("otu_table2.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')#读取默认路径下的任意一个物种数据框
otu <- otu[,-ncol(otu)]
otu <- data.frame(t(otu))#如果行为otu，列为样方数据，则需进行转置
bray_dis <- vegdist(otu, method = 'bray')

nmds_dis <- metaMDS(bray_dis, k = 2, try=100,trymax=100)
nmds_dis$stress#查看stress值

###坐标轴保留4位小数，scale_y_continuous(labels=scaleFUN)######字母f之前的数字表示保留的小数位
scaleFUN <- function(x) sprintf("%.4f", x)
stress <- scaleFUN(nmds_dis$stress)
stress

nmds_dis_site <- data.frame(nmds_dis$points)#样方得分
nmds_dis_species <- wascores(nmds_dis$points, otu)#物种得分

#ordiplot(nmds_dis, type = 'none', main = paste('样方,Stress =', round(nmds_dis$stress, 4)))
#points(nmds_dis,pch =c(rep('+',4),rep('⊿',4),rep('Ж',4),rep('Ф',4)), cex = 1, col = c(rep('red', 4), rep('orange', 4), rep('green3', 4),rep('blue',4)))
#legend("bottomright",                                 #图例位置为右上角
 #      legend=c("C","W","N","WN"),        #图例内容
  #     col=c("red","orange","green3","blue"),                 #图例颜色
   #    pch =c(rep('+',1),rep('⊿',1),rep('Ж',1),rep('Ф',1)),bty="n",ncol=1,pt.cex=1)

#添加分组信息
nmds_dis_site$name <- rownames(nmds_dis_site)
mapping <- read.table("mapping.xls",row.names = 1,header=T,comment.char='',sep='\t')#读取分组文件
nmds_dis_site <- nmds_dis_site[match(row.names(mapping),row.names(nmds_dis_site)),]
nmds_dis_site$Treat <- mapping$Treat

p <- ggplot(data = nmds_dis_site, aes(MDS1, MDS2, shape = Treat , color = Treat)) + #取MDS1和MDS2来绘图
  theme_classic()+ #定义经典背景
  geom_point(aes(color = Treat)) + #颜色按分组填充
  stat_ellipse(aes(fill = Treat), geom = 'polygon', level = 0.68, alpha = 0.2, show.legend = FALSE) +       #添加置信椭圆，注意不是聚类
  #geom_polygon(data=nmds_dis_site,aes(fill=Treat,alpha=0.2),show.legend=F)+
  #annotate("text",x=min(nmds_dis_site$MDS1),y=min(nmds_dis_site$MDS2),hjust=0,vjust=0,label=paste("Stress:",nmds_dis$stress, 4))+
  #scale_shape_manual(values = c(3, 8, 10, 19,21,25)) +
  #scale_color_manual(values = c('red3', 'orange3', 'green3','blue','yellow','gray')) +
  #scale_fill_manual(values = c('red3', 'orange3', 'green3','blue','yellow','gray')) +
  #theme(legend.position = 'none') + #去掉图例
  geom_vline(xintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #中间竖线
  geom_hline(yintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #中间横线
  labs(x = 'NMDS1', y = 'NMDS2') +
  annotate('text', label = paste('Stress =', stress), x = 0.28, y = 0.3, size = 4, colour = 'black') #标注stress函数值
#annotate('text', label = 'C', x = -0.1, y = 0.1, size = 5, colour = 'red3') + #定义大样点文字A
#annotate('text', label = 'W', x = -0.03, y = 0.2, size = 5, colour = 'orange3') + #定义大样点文字B
#annotate('text', label = 'N', x = 0.1, y = -0.2, size = 5, colour = 'green3')+ #定义大样点文字C
#annotate('text', label = 'WN', x = 0.4, y = -0.075, size = 5, colour = 'blue') #定义大样点文字C

p #出图
p1 <- p+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'))
#经典方框背景
p1
