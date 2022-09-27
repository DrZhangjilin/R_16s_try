#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("agricolae")
#install.packages("vegan")

setwd("G:/R/Rdata2")#Rdata3为ITS的文件夹，Rdata2为16S的文件夹

#安装需要的包，如已安装，可省略这一步
library(ggplot2)
library(ggpubr)
library(agricolae)
library(vegan)
#加载需要的分析包，如果还没有安装这些包，请使用

###ITS和16S的otu文件名一样
otu <- read.table("otu_table.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')

#'E:/lesson1/feature-table.taxonomy.txt'为文件路径，注意斜线方向
#row.names = 1指定第一列为行名
#skip=1跳过第一行不读
#header=T指定第一个有效行为列名
#sep='\t'表示指定制表符为分隔符
#comment.char=''表示设置注释符号为空字符‘’，这样#后面的内容就不会被省略

otu1 <- otu[,-ncol(otu)]
#删除多余的列
otu2=t(otu1)
#转置

shannon=diversity(otu2,"shannon",base=2)
#计算香农指数
simpson=diversity(otu2,"simpson")
#计算辛普森指数
alpha=data.frame(shannon,simpson)
#合并数据
write.table(alpha,"alpha-summary3.tsv",sep = '\t',quote=F)
#储存结果，根据此结果进行SPSS单因素方差分析


#############手动修改显著性标签


#map<-read.table("mapping.xls",row.names = 1,header = T,sep='\t',comment.char='',check.names=F)
#读取分组表格
#row.names = 1表示指定第一列为行名
#header = T表示指定第一个有效行为列名
#sep='\t'表示指定制表符为分隔符
#comment.char=''表示设置注释符号为空字符‘’，这样#后面的内容就不会被省略
#check.names=F表示读取过程中不对行名和列名做任何修改
#group<-map['Treat']
#提取需要的分组，'Group1'为表中分组列名



#####手动修改保存的上述合并过的文件并读入（修改字母标记顺序）
alpha<-read.table('alpha-summary3.xls',header = T,row.names = 1,sep = '\t')


#读取alpha多样性表
#alpha2<-alpha[match(rownames(group),rownames(alpha)),]
#print(rownames(group) == rownames(alpha2))
#重排alpha的行的顺序，使其与group的样本id（行名）顺序一致
#data<-data.frame(group,alpha)
#合并两个表格



###坐标轴保留3位小数，scale_y_continuous(labels=scaleFUN)######字母f之前的数字表示保留的小数位
scaleFUN <- function(x) sprintf("%.3f", x)

#下述p1、p2代表16S，p3、p4代表ITS
data <- alpha
p3 <- ggplot(data,aes(x=Treat,y=shannon))+ stat_boxplot(geom="errorbar",size=0.6,width=0.2)+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4,outlier.stroke =0.5,outlier.alpha = 0.5,outlier.fill = "red",aes(fill=Treat))+geom_jitter(width=0.25,shape=20,size=2.5)+stat_summary(fun="mean",geom="point",shape=23,size=3,fill="white")+
  scale_y_continuous(labels=scaleFUN)
p4 <- ggplot(data,aes(x=Treat,y=simpson))+ stat_boxplot(geom="errorbar",size=0.6,width=0.2)+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4,outlier.stroke =0.5,outlier.alpha = 0.5,outlier.fill = "red",aes(fill=Treat))+geom_jitter(width=0.25,shape=20,size=2.5)+stat_summary(fun="mean",geom="point",shape=23,size=3,fill="white")+
  scale_y_continuous(labels=scaleFUN)
p3
p4
#带散点和平均值
#data = data指定数据表格
#x=Group1指定作为x轴的数据列名
#y=shannon指定作为y轴的数据列名
#geom_boxplot()表示画箱线图


#依据spss单因素方差分析结果，手动整理数据，把显著性字母对应标记
library(openxlsx)
a <- read.xlsx("alpha_data1.xlsx",sheet=1)
p3<-p3+geom_text(data=a,aes(x=treat,y=5.2,label=Shannon))#label表示要标记的显著性字母的一列
p3
p4<-p4+geom_text(data=a,aes(x=treat,y=1,label=Simpson))
p4


#去掉图例（为了合并后好看）
p1 <- p1+theme(legend.position = 'none')
p1
p2 <- p2+theme(legend.position = 'none')
p2

#mycompare=list(c('A','B'),c('A','C'),c('B','C'))
#指定多重比较的分组对
#p1<-p+stat_compare_means(comparisons=mycompare,label = "p.signif",method = 'wilcox')
#p1
#添加星号，添加显著性标记的第一种方法,使用wilcoxon非参数检验方法

#anova <- aov(shannon~Description,data = data)
#plotdata<-duncan.test(anova,"Description",console = TRUE, alpha = 0.05)
#plotdata<-data.frame(id=rownames(plotdata$groups),plotdata$groups)
#p2<-p+geom_text(data = plotdata,aes(x=id,y=7.6,label=groups))
#p2
#添加显著性字母，添加显著性标记的第二种方法，注意此方法为参数检验，要求alpha多样性指数服从正太分布



#坐标轴标签和主题标签
p3 <- p3+labs(x="",y="",title="ITS")
p3
p4 <- p4+labs(x="",y="",title="")
p4

#更改字体字号和大小
p3=p3+theme(text = element_text(size = 15,face = "bold"))
p4=p4+theme(text = element_text(size = 15,face = "bold"))
p3
p4

#经典背景
p3 <- p3+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())
p3
p4 <- p4+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())
p4


#拼图准备（上层图像去掉横坐标，左侧图像去掉图例，右侧图像去掉纵坐标轴标签，下层图像去掉主题标签）
p1 <- p1+theme(axis.text.x = element_blank())
p3 <- p3+theme(axis.text.x = element_blank())

p1
p2
p3
p4

####调整单幅图的边界，减少拼图后的缝隙！顺序（上，右，下，左）
p6 <- p1+theme( plot.margin=unit(c(0.5,0.5,-0.5,0.5),"cm" ))
p6
p7 <- p2+theme( plot.margin=unit(c(-0.5,0.5,0.5,0.5),"cm" ))
p7
p8 <- p3+theme( plot.margin=unit(c(0.5,0,-0.5,-0.5),"cm" ))
p8
p9 <- p4+theme( plot.margin=unit(c(-0.5,0,0.5,-0.5),"cm" ))
p9

#合并图像
library(gridExtra)
p10 <- grid.arrange(p6,p8,p7,p9,nrow=2,ncol=2,widths=c(2,2.15))

#保存方案1
dev.copy(pdf,"whatever.pdf")
dev.off()
#保存方案2
pdf("filename.pdf", width = 10, height = 9) # Open a new pdf file
grid.arrange(p6,p8,p7,p9,nrow=2,ncol=2,widths=c(2,2.15))# Write the grid.arrange in the file
dev.off()
#！！！！！保存方案3（不适用）
ggsave(file="alpha多样性比较（4*4）.tiff",plot=arrangeGrob(p1,p3,p2,p4,nrow=2,ncol=2),units="px",width=1000,height=900)


#ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
#p6 <- arrangeGrob(p1,p3,p2,p4,nrow=2,ncol=2)
#ggsave(file="alpha多样性比较（4*4）.tiff",plot=p6,units="px",width=1000,height=900)

#library(grid)
#grid.draw(p5)
#ggsave(file="alpha多样性比较（4*4）.tiff",plot=p5,units="px",width=1000,height=900)

#修改字号字形
#ggsave(plot = p,'E:/lesson1/shannon_boxplot2.png',height = 7,width = 6,dpi = 300)
#plot = p指定刚刚叠加修饰好的p
# ‘E:/practice1/shannon_boxplot.pdf’表示存储路径和文件名字，
#注意文件的后缀名，它将决定图片格式
#建议储存为pdf格式，具体原因，请关注后面的讲座