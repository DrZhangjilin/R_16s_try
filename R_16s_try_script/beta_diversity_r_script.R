if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.13")
BiocManager::install("phyloseq")
#> chooseBioCmirror()�����޸���ѡ��
#��װphyloseq
library("phyloseq")
library("ggplot2")
library("ggalt")
#���ر�Ҫ�İ�
setwd("G:/R/Rdata2")
#���ù���Ŀ¼�������ݴ洢��Ŀ¼
qiimedata <- import_qiime(otufilename = "otu_biaozhun.xls", mapfilename = "map1.xls", treefilename="rep_phylo.tre", refseqfilename = "rep_set_aligned_pfiltered.fasta")
#qiimedata <- import_qiime(otufilename = "otu_table2.xls", mapfilename = "map.txt")
#��ȡ����,��ȡ���ݣ����������ļ�����ע��Ӻ�׺
#otufilenameָ��out����mapfilenameָ��map�ļ����������ݣ��� treefilenameָ���и��������ļ���refseqfilenameָ�����������ļ�
otu<-qiimedata@otu_table@.Data
#�ӵ������������ȡotu����
#otu<-otu_table(qiimedata)#��ѡ��otu������ȡ����
otu <- t(otu)
mapping <- read.table("map1.xls",row.names = 1,header=T,comment.char='',sep='\t')#��ȡ�����ļ�
otu <- otu[match(row.names(mapping),row.names(otu)),]

#sum_of_otus<-colSums(otu)
#�������otu��⵽����������
#selected_otu<-names(sum_of_otus)[sum_of_otus>18]
#��ȡ������������18��otu id

#sub_qiimedata <- prune_taxa(selected_otu, qiimedata)
#ɸѡ������������10��otu phyloseq����
#sub_qiimedata=subset_taxa(sub_qiimedata,Kingdom=="Bacteria")#����ע�ͷ������ɸѡotu�ķ���
set.seed(123)
#weighted_unifrac <- distance(sub_qiimedata, method='wunifrac')
weighted_unifrac <- distance(qiimedata, method='wunifrac')
#�����������ȨUniFrac���󣬱�����Ҫ�������ļ�
#unweighted_unifrac <- distance(sub_qiimedata, method='unifrac')
#����������Ǽ�ȨUniFrac���󣬱�����Ҫ�������ļ�
#bray_curtis <- distance(sub_qiimedata, method='bray')
#����������Bray-Curtis������󣬿��Բ���Ҫ�������ļ�
write.table(as.matrix(weighted_unifrac),"weighted_unifrac_biaozhun.txt",sep = '\t',quote = FALSE,col.names = NA)
#write.table(as.matrix(unweighted_unifrac),"unweighted_unifrac1.txt",sep = '\t',quote = FALSE,col.names = NA)
#write.table(as.matrix(bray_curtis),"bray_curtis1.txt",sep = '\t',quote = FALSE,col.names = NA)
#���������������


#nmds_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "NMDS",try=100,trymax=100)
#nmds_of_bray_curtis1<-ordinate(physeq=sub_qiimedata,distance = 'wunifrac',method = "NMDS")
#nmds_of_bray_curtis2<-ordinate(physeq=sub_qiimedata,distance = 'unifrac',method = "NMDS")
#����Bray-Curtis��������NMDS�������
set.seed(123)
nmds_of_bray_curtis<-ordinate(physeq=qiimedata,distance = 'wunifrac',method = "NMDS")
p<-plot_ordination(qiimedata, nmds_of_bray_curtis,
                   type="samples", color="Treat",shape = "Treat")+
  geom_point(size=4)+
  geom_vline(xintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #����������������
  geom_hline(yintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #����������������
  #stat_ellipse(level=0.68)+#������Բ��һ�����������ʱ���ã�
  geom_encircle(aes(fill=Treat), alpha = 0.2, show.legend = F) +#�����߰�Χ�Ķ����
  annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis$stress,4)),
           x=0,y=0.06,size=6,colour='black',family="serif")+
  theme_bw()+#�׵׻��߱���
  theme(legend.key = element_blank())+
  theme(text=element_text(size=16,  family="serif",face = "bold"),
        axis.text.x = element_text(size=15,  family="serif",colour = "black"),
        axis.text.y = element_text(size=15,  family="serif",colour = "black"))#serif��R�б�ʾ����������
p
#p1<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis, type="samples", color="Treat",shape = "Treat")
#p1
#p2<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis1, type="samples", color="Treat",shape = "Treat")
#p2
#p3<-plot_ordination(sub_qiimedata, nmds_of_bray_curtis2, type="samples", color="Treat",shape = "Treat")
#p3
#��NMDS�������������ӻ�

#p1<-p1 + geom_point(size=3)+ stat_ellipse(level=0.68)+theme(text = element_text(size = 15))+annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis$stress,4)),x=-0.3,y=0.3,size=6,colour='black')
#p1
#p2<-p2 + geom_point(size=3)+ stat_ellipse(level=0.68)+theme(text = element_text(size = 15))+annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis1$stress,4)),x=-0.05,y=0.05,size=6,colour='black')
#p2
#p3<-p3+ geom_point(size=3)+ stat_ellipse(level=0.68)+theme(text = element_text(size = 15))+annotate('text',label=paste('Stress=',round(nmds_of_bray_curtis2$stress,4)),x=0.07,y=0.15,size=6,colour='black')
#p3
#��ͼƬ�����ʵ����Σ� stat_ellipse()����Բ�� ggtitle()�ӱ���


#ggsave(plot = p,"nmds_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#����ͼƬ

set.seed(123)
pcoa_of_bray_curtis<-ordinate(physeq=qiimedata,distance = 'wunifrac',method = "PCoA")

#���¹���Ϊ��ȡ1�ᡢ2��Ľ��Ͷ�
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
  geom_vline(xintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #����������������
  geom_hline(yintercept = 0, color = 'DimGrey', size = 0.5,linetype="dashed") + #����������������
  #stat_ellipse(level=0.68)+
  geom_encircle(s_shape=1,expand=0.01,aes(fill=Treat), alpha = 0.2, show.legend = F)+
  labs(x = pcoa1, y = pcoa2)+
  theme_bw()+#�׵׻��߱���
  theme(legend.key = element_blank())+
  #coord_fixed()+
  theme(text=element_text(size=10,  family="serif",face = "bold"),
        axis.text.x = element_text(size=8,  family="serif",colour = "black"),
        axis.text.y = element_text(size=8,  family="serif",colour = "black"))#serif��R�б�ʾ����������
p
ggsave("./PCoA_20220702.pdf", p, width = 100, height = 80, units = "mm")


#pcoa_of_bray_curtis<-ordinate(physeq=sub_qiimedata,distance = 'bray',method = "PCoA")
#����Bray-Curtis��������PCoA�������
#p4<-plot_ordination(sub_qiimedata, pcoa_of_bray_curtis, type="samples", color="Description",shape = "Description")
#p4
#��PCoA�������������ӻ�
#p5<-p4+ scale_colour_manual(values=c("#B22222","#FF00FF","#0000CD","#00FF00","#FFC125","8B7D7B")) +stat_ellipse(level=0.85)+ geom_point(size=3) +theme(text = element_text(size = 15))
#p5
#��ͼƬ�����ʵ�����
#��scale_colour_manual(values=c())�Զ�����ɫ���ɲ���ɫ��16���ƶ��ձ�
#ggsave(plot = p,"pcoa_of_bary_curtis.pdf",dpi = 300,width = 7,height = 6)
#����ͼƬ

#p2 <- p1+theme_bw(base_line_size = 1.3,base_rect_size = 1.3)#���þ��䱳��
#p2

p1 <- p1+ theme(legend.key = element_blank(),#ȥ��ͼ������ɫ
                panel.background = element_rect(fill = "transparent",colour = NA),
                panel.grid.minor = element_blank(),
                panel.grid.major = element_blank(),
                plot.background = element_rect(fill = "transparent",colour = NA),
                axis.line=element_line(colour="black"))#����������
p1
#p2 <- p2+ theme(legend.key = element_blank(),#ȥ��ͼ������ɫ
               # panel.background = element_rect(fill = "transparent",colour = NA),
               # panel.grid.minor = element_blank(),
               # panel.grid.major = element_blank(),
               # plot.background = element_rect(fill = "transparent",colour = NA),
               # axis.line=element_line(colour="black"))#����������
#p2
#p3 <- p3+ theme(legend.key = element_blank(),#ȥ��ͼ������ɫ
               # panel.background = element_rect(fill = "transparent",colour = NA),
               # panel.grid.minor = element_blank(),
               # panel.grid.major = element_blank(),
               # plot.background = element_rect(fill = "transparent",colour = NA),
               # axis.line=element_line(colour="black"))#����������
#p3