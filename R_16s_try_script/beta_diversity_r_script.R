if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
BiocManager::install("phyloseq")
#> chooseBioCmirror()镜像修改与选择
#安装phyloseq
library("phyloseq")
library("ggplot2")
library("ggalt")
#加载必要的包
setwd("G:/R/Rdata2")
#设置工作目录，即数据存储的目录
qiimedata <- import_qiime(otufilename = "otu_biaozhun.xls", mapfilename = "map1.xls", treefilename="rep_phylo.tre", refseqfilename = "rep_set_aligned_pfiltered.fasta")
#qiimedata <- import_qiime(otufilename = "otu_table2.xls", mapfilename = "map.txt")
#读取数据,读取数据，参数都是文件名，注意加后缀
#otufilename指定out表格，mapfilename指定map文件（分组数据）， treefilename指定有根进化树文件，refseqfilename指定代表序列文件
otu<-qiimedata@otu_table@.Data
#从导入的数据中提取otu表格
#otu<-otu_table(qiimedata)#可选的otu表格提取方法
otu <- t(otu)
mapping <- read.table("map1.xls",row.names = 1,header=T,comment.char='',sep='\t')#读取分组文件
otu <- otu[match(row.names(mapping),row.names(otu)),]

#sum_of_otus<-colSums(otu)
#计算各个otu检测到的总序列数
#selected_otu<-names(sum_of_otus)[sum_of_otus>18]
#获取总序列数大于18的otu id

#sub_qiimedata <- prune_taxa(selected_otu, qiimedata)
#筛选总序列数大于10的otu phyloseq数据
#sub_qiimedata=subset_taxa(sub_qiimedata,Kingdom=="Bacteria")#根据注释分类进行筛选otu的方法
set.seed(123)
#weighted_unifrac <- distance(sub_qiimedata, method='wunifrac')
weighted_unifrac <- distance(qiimedata, method='wunifrac')
#计算样本间加权UniFrac矩阵，必须需要进化树文件
#unweighted_unifrac <- distance(sub_qiimedata, method='unifrac')
#计算样本间非加权UniFrac矩阵，必须需要进化树文件
#bray_curtis <- distance(sub_qiimedata, method='bray')
#计算样本间Bray-Curtis距离矩阵，可以不需要进化树文件
write.table(as.matrix(weighted_unifrac),"weighted_unifrac_biaozhun.txt",sep = '\t',quote = FALSE,col.names = NA)
#write.table(as.matrix(unweighted_unifrac),"unweighted_unifrac1.txt",sep = '\t',quote = FALSE,col.names = NA)
#write.table(as.matrix(bray_curtis),"bray_curtis1.txt",sep = '\t',quote = FALSE,col.names = NA)
#保存三个距离矩阵


#nmds_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "NMDS",try=100,trymax=100)
#nmds_of_bray_curtis1<-ordinate(physeq=sub_qiimedata,distance = 'wunifrac',method = "NMDS")
#nmds_of_bray_curtis2<-ordinate(physeq=sub_qiimedata,distance = 'unifrac',method = "NMDS")
#基于Bray-Curtis距离矩阵的NMDS排序分析
set.seed(123)
nmds_of_bray_curtis<-ordinate(physeq=qiimedata,distance = 'wunifrac',method = "NMDS")
p<-plot_ordination(qiimedata, nmds_of_bray_curtis,
                   type="samples", color="Treat",shape = "Treat")+
  geom_point(size=4)+
  geom_vline(xintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #过坐标零点的两条线
  geom_hline(yintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #过坐标零点的两条线
  #stat_ellipse(level=0.68)+#置信椭圆（一般样本量多的时候用）
  geom_encircle(aes(fill=Treat), alpha = 0.2, show.legend = F) +#样本边包围的多边形
  annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis$stress,4)),
           x=0,y=0.06,size=6,colour='black',family="serif")+
  theme_bw()+#白底灰线背景
  theme(legend.key = element_blank())+
  theme(text=element_text(size=16,  family="serif",face = "bold"),
        axis.text.x = element_text(size=15,  family="serif",colour = "black"),
        axis.text.y = element_text(size=15,  family="serif",colour = "black"))#serif在R中表示新罗马字体
p
#p1<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis, type="samples", color="Treat",shape = "Treat")
#p1
#p2<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis1, type="samples", color="Treat",shape = "Treat")
#p2
#p3<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis2, type="samples", color="Treat",shape = "Treat")
#p3
#将NMDS排序分析结果可视化

#p1<-p1 + geom_point(size=3)+ stat_ellipse(level=0.68)+theme(text = element_text(size = 15))+annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis$stress,4)),x=-0.3,y=0.3,size=6,colour='black')
#p1
#p2<-p2 + geom_point(size=3)+ stat_ellipse(level=0.68)+theme(text = element_text(size = 15))+annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis1$stress,4)),x=-0.05,y=0.05,size=6,colour='black')
#p2
#p3<-p3+ geom_point(size=3)+ stat_ellipse(level=0.68)+theme(text = element_text(size = 15))+annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis2$stress,4)),x=0.07,y=0.15,size=6,colour='black')
#p3
#对图片进行适当修饰， stat_ellipse()加椭圆， ggtitle()加标题


#ggsave(plot = p,"nmds_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#保存图片

set.seed(123)
pcoa_of_bray_curtis<-ordinate(physeq=qiimedata,distance = 'wunifrac',method = "PCoA")

#以下过程为提取1轴、2轴的解释度
a <- pcoa_of_bray_curtis[["values"]]
a <- data.frame(a)
pcoa_exp <- a$Eigenvalues/sum(a$Eigenvalues)
pcoa1 <- paste('PCoA axis1 :', round(100*pcoa_exp[1], 1), '%')
pcoa2 <- paste('PCoA axis2 :', round(100*pcoa_exp[2], 2), '%')

order_y=c("C","LA","LR","PR","LAPR","LRPR")
qiimedata@sam_data$Treat=factor(qiimedata@sam_data$Treat,
                                levels = order_y)

p<-plot_ordination(qiimedata, pcoa_of_bray_curtis,
                   type="samples", color="Treat",shape = "Treat")+
  geom_point(size=2)+
  geom_vline(xintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #过坐标零点的两条线
  geom_hline(yintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #过坐标零点的两条线
  #stat_ellipse(level=0.68)+
  geom_encircle(s_shape=1,expand=0.01,aes(fill=Treat), alpha = 0.2, show.legend = F)+
  labs(x = pcoa1, y = pcoa2)+
  theme_bw()+#白底灰线背景
  theme(legend.key = element_blank())+
  #coord_fixed()+
  theme(text=element_text(size=10,  family="serif",face = "bold"),
        axis.text.x = element_text(size=8,  family="serif",colour = "black"),
        axis.text.y = element_text(size=8,  family="serif",colour = "black"))#serif在R中表示新罗马字体
p
ggsave("./PCoA_20220702.pdf", p, width = 100, height = 80, units = "mm")


#pcoa_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "PCoA")
#基于Bray-Curtis距离矩阵的PCoA排序分析
#p4<-plot_ordination(sub_qiimedata, pcoa_of_bray_curtis, type="samples", color="Description",shape = "Description")
#p4
#将PCoA排序分析结果可视化
#p5<-p4+ scale_colour_manual(values=c("#B22222","#FF00FF","#0000CD","#00FF00","#FFC125","8B7D7B")) +stat_ellipse(level=0.85)+ geom_point(size=3) +theme(text = element_text(size = 15))
#p5
#对图片进行适当修饰
#用scale_colour_manual(values=c())自定义颜色，可查颜色的16进制对照表
#ggsave(plot = p,"pcoa_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#保存图片

#p2 <- p1+theme_bw(base_line_size = 1.3,base_rect_size = 1.3)#设置经典背景
#p2

p1 <- p1+ theme(legend.key = element_blank(),#去掉图例背景色
                panel.background = element_rect(fill = "transparent",colour = NA),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                plot.background = element_rect(fill = "transparent",colour = NA),
                axis.line=element_line(colour="black"))#加上坐标轴
p1
#p2 <- p2+ theme(legend.key = element_blank(),#去掉图例背景色
               # panel.background = element_rect(fill = "transparent",colour = NA),
               # panel.grid.minor = element_blank(),
               # panel.grid.major = element_blank(),
               # plot.background = element_rect(fill = "transparent",colour = NA),
               # axis.line=element_line(colour="black"))#加上坐标轴
#p2
#p3 <- p3+ theme(legend.key = element_blank(),#去掉图例背景色
               # panel.background = element_rect(fill = "transparent",colour = NA),
               # panel.grid.minor = element_blank(),
               # panel.grid.major = element_blank(),
               # plot.background = element_rect(fill = "transparent",colour = NA),
               # axis.line=element_line(colour="black"))#加上坐标轴
#p3
