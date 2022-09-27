#install.packages("ggplot2")
#install.packages("ggpubr")
#install.packages("agricolae")
#install.packages("vegan")

setwd("G:/R/Rdata2")#Rdata3ΪITS���ļ��У�Rdata2Ϊ16S���ļ���

#��װ��Ҫ�İ������Ѱ�װ����ʡ����һ��
library(ggplot2)
library(ggpubr)
library(agricolae)
library(vegan)
#������Ҫ�ķ������������û�а�װ��Щ������ʹ��

###ITS��16S��otu�ļ���һ��
otu <- read.table("otu_table.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')

#'E:/lesson1/feature-table.taxonomy.txt'Ϊ�ļ�·����ע��б�߷���
#row.names = 1ָ����һ��Ϊ����
#skip=1������һ�в���
#header=Tָ����һ����Ч��Ϊ����
#sep='\t'��ʾָ���Ʊ���Ϊ�ָ���
#comment.char=''��ʾ����ע�ͷ���Ϊ���ַ�����������#��������ݾͲ��ᱻʡ��

otu1 <- otu[,-ncol(otu)]
#ɾ���������
otu2=t(otu1)
#ת��

shannon=diversity(otu2,"shannon",base=2)
#������ũָ��
simpson=diversity(otu2,"simpson")
#��������ɭָ��
alpha=data.frame(shannon,simpson)
#�ϲ�����
write.table(alpha,"alpha-summary3.tsv",sep = '\t',quote=F)
#�����������ݴ˽������SPSS�����ط������


#############�ֶ��޸������Ա�ǩ


#map<-read.table("mapping.xls",row.names = 1,header = T,sep='\t',comment.char='',check.names=F)
#��ȡ�������
#row.names = 1��ʾָ����һ��Ϊ����
#header = T��ʾָ����һ����Ч��Ϊ����
#sep='\t'��ʾָ���Ʊ���Ϊ�ָ���
#comment.char=''��ʾ����ע�ͷ���Ϊ���ַ�����������#��������ݾͲ��ᱻʡ��
#check.names=F��ʾ��ȡ�����в����������������κ��޸�
#group<-map['Treat']
#��ȡ��Ҫ�ķ��飬'Group1'Ϊ���з�������



#####�ֶ��޸ı���������ϲ������ļ������루�޸���ĸ���˳��
alpha<-read.table('alpha-summary3.xls',header = T,row.names = 1,sep = '\t')


#��ȡalpha�����Ա�
#alpha2<-alpha[match(rownames(group),rownames(alpha)),]
#print(rownames(group) == rownames(alpha2))
#����alpha���е�˳��ʹ����group������id��������˳��һ��
#data<-data.frame(group,alpha)
#�ϲ���������



###�����ᱣ��3λС����scale_y_continuous(labels=scaleFUN)######��ĸf֮ǰ�����ֱ�ʾ������С��λ
scaleFUN <- function(x) sprintf("%.3f", x)

#����p1��p2����16S��p3��p4����ITS
data <- alpha
p3 <- ggplot(data,aes(x=Treat,y=shannon))+ stat_boxplot(geom="errorbar",size=0.6,width=0.2)+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4,outlier.stroke =0.5,outlier.alpha = 0.5,outlier.fill = "red",aes(fill=Treat))+geom_jitter(width=0.25,shape=20,size=2.5)+stat_summary(fun="mean",geom="point",shape=23,size=3,fill="white")+
  scale_y_continuous(labels=scaleFUN)
p4 <- ggplot(data,aes(x=Treat,y=simpson))+ stat_boxplot(geom="errorbar",size=0.6,width=0.2)+geom_boxplot(outlier.colour="red", outlier.shape=8, outlier.size=4,outlier.stroke =0.5,outlier.alpha = 0.5,outlier.fill = "red",aes(fill=Treat))+geom_jitter(width=0.25,shape=20,size=2.5)+stat_summary(fun="mean",geom="point",shape=23,size=3,fill="white")+
  scale_y_continuous(labels=scaleFUN)
p3
p4
#��ɢ���ƽ��ֵ
#data = dataָ�����ݱ���
#x=Group1ָ����Ϊx�����������
#y=shannonָ����Ϊy�����������
#geom_boxplot()��ʾ������ͼ


#����spss�����ط������������ֶ��������ݣ�����������ĸ��Ӧ���
library(openxlsx)
a <- read.xlsx("alpha_data1.xlsx",sheet=1)
p3<-p3+geom_text(data=a,aes(x=treat,y=5.2,label=Shannon))#label��ʾҪ��ǵ���������ĸ��һ��
p3
p4<-p4+geom_text(data=a,aes(x=treat,y=1,label=Simpson))
p4


#ȥ��ͼ����Ϊ�˺ϲ���ÿ���
p1 <- p1+theme(legend.position = 'none')
p1
p2 <- p2+theme(legend.position = 'none')
p2

#mycompare=list(c('A','B'),c('A','C'),c('B','C'))
#ָ�����رȽϵķ����
#p1<-p+stat_compare_means(comparisons=mycompare,label = "p.signif",method = 'wilcox')
#p1
#�����Ǻţ����������Ա�ǵĵ�һ�ַ���,ʹ��wilcoxon�ǲ������鷽��

#anova <- aov(shannon~Description,data = data)
#plotdata<-duncan.test(anova,"Description",console = TRUE, alpha = 0.05)
#plotdata<-data.frame(id=rownames(plotdata$groups),plotdata$groups)
#p2<-p+geom_text(data = plotdata,aes(x=id,y=7.6,label=groups))
#p2
#������������ĸ�����������Ա�ǵĵڶ��ַ�����ע��˷���Ϊ�������飬Ҫ��alpha������ָ��������̫�ֲ�



#�������ǩ�������ǩ
p3 <- p3+labs(x="",y="",title="ITS")
p3
p4 <- p4+labs(x="",y="",title="")
p4

#���������ֺźʹ�С
p3=p3+theme(text = element_text(size = 15,face = "bold"))
p4=p4+theme(text = element_text(size = 15,face = "bold"))
p3
p4

#���䱳��
p3 <- p3+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())
p3
p4 <- p4+theme(panel.grid = element_blank(), panel.background = element_rect(fill = 'transparent', color = 'black'), legend.title = element_blank(), legend.key = element_blank())
p4


#ƴͼ׼�����ϲ�ͼ��ȥ�������꣬���ͼ��ȥ��ͼ�����Ҳ�ͼ��ȥ�����������ǩ���²�ͼ��ȥ�������ǩ��
p1 <- p1+theme(axis.text.x = element_blank())
p3 <- p3+theme(axis.text.x = element_blank())

p1
p2
p3
p4

####��������ͼ�ı߽磬����ƴͼ��ķ�϶��˳���ϣ��ң��£���
p6 <- p1+theme( plot.margin=unit(c(0.5,0.5,-0.5,0.5),"cm" ))
p6
p7 <- p2+theme( plot.margin=unit(c(-0.5,0.5,0.5,0.5),"cm" ))
p7
p8 <- p3+theme( plot.margin=unit(c(0.5,0,-0.5,-0.5),"cm" ))
p8
p9 <- p4+theme( plot.margin=unit(c(-0.5,0,0.5,-0.5),"cm" ))
p9

#�ϲ�ͼ��
library(gridExtra)
p10 <- grid.arrange(p6,p8,p7,p9,nrow=2,ncol=2,widths=c(2,2.15))

#���淽��1
dev.copy(pdf,"whatever.pdf")
dev.off()
#���淽��2
pdf("filename.pdf", width = 10, height = 9) # Open a new pdf file
grid.arrange(p6,p8,p7,p9,nrow=2,ncol=2,widths=c(2,2.15))# Write the grid.arrange in the file
dev.off()
#�������������淽��3�������ã�
ggsave(file="alpha�����ԱȽϣ�4*4��.tiff",plot=arrangeGrob(p1,p3,p2,p4,nrow=2,ncol=2),units="px",width=1000,height=900)


#ggsave <- ggplot2::ggsave; body(ggsave) <- body(ggplot2::ggsave)[-2]
#p6 <- arrangeGrob(p1,p3,p2,p4,nrow=2,ncol=2)
#ggsave(file="alpha�����ԱȽϣ�4*4��.tiff",plot=p6,units="px",width=1000,height=900)

#library(grid)
#grid.draw(p5)
#ggsave(file="alpha�����ԱȽϣ�4*4��.tiff",plot=p5,units="px",width=1000,height=900)

#�޸��ֺ�����
#ggsave(plot = p,'E:/lesson1/shannon_boxplot2.png',height = 7,width = 6,dpi = 300)
#plot = pָ���ոյ������κõ�p
# ��E:/practice1/shannon_boxplot.pdf����ʾ�洢·�����ļ����֣�
#ע���ļ��ĺ�׺������������ͼƬ��ʽ
#���鴢��Ϊpdf��ʽ������ԭ�����ע����Ľ���