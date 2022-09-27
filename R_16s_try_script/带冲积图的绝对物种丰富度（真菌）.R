setwd("G:/R/Rdata2")
rm(list=ls())
library(tidyverse)#数据整理与数据转换包，用了一些更好用更易懂的函数
library(ggprism)
library(ggsci)
library(vegan)
library(ggalluvial)
#library(RColorBrewer)
#library(ggsci)
otu <- read.table("otu_table2.xls",row.names = 1,skip=1,header=T,comment.char='',sep='\t')#细菌
otu=otu[,-ncol(otu)]
head(otu, n = 3)
#释放即为相对丰度
otu <-(decostand(otu,'total',2)/6)#按列标准化otu,36个样本，6个处理，6 组重复，除以6是为了使和为100%
colSums(otu)#若想后面做成相对丰度的差异比较，可把这两行释放出来即可
tax <- read.table('taxonomy.xls',row.names = 1,sep='\t',comment.char='',header=T)
head(tax, n = 3)
metadata<- read.table('map1.xls',sep='\t',header = T)
metadata=metadata[,c(1:2)]
head(metadata, n = 3)
dat <- merge(x=otu,y=tax,by='row.names')
head(dat, n = 3)
dat =dplyr::rename(dat,OTUID = Row.names)
head(dat, n = 3)
##按Phylum水平分组汇总(根据自己需求更改要展示的物种水平)
#aa<-aggregate(dat[,2:(ncol(otu)+1)],by=list(dat$Class),FUN=sum)
aa<-aggregate(dat[,2:(ncol(otu)+1)],by=list(dat$Phylum),FUN=sum)
head(aa)

#######################三种排序方法，任选其一
#1
# aa<- mutate(aa,v=apply(aa[,c(2:ncol(aa))],1,sum))
# cc<- arrange(aa,desc(v))
# head(cc)
# cc<-select(cc,-v)
# head(cc)
# row.names(cc)=cc$Phylum
# head(cc)
# cc<-select(cc,-Phylum)
# head(cc)
#2
# row.names(aa)=aa$Phylum
# head(aa)
# aa<-select(aa,-Phylum)
# head(aa)
# cc<-aa[order(rowSums(aa),decreasing=T),]
#3
row.names(aa)=aa$Group.1
head(aa)
aa<-dplyr::select(aa,-Group.1)
head(aa, n = 3)
#根据行求和结果对数据排序
order<-sort(rowSums(aa[,2:ncol(aa)]),index.return=TRUE,decreasing=T)
#根据列求和结果对表格排序
cc<-aa[order$ix,]
head(cc, n = 3)
##只展示排名前10的物种，之后的算作Others(根据需求改数字)
dd<-rbind(colSums(cc[10:as.numeric(length(rownames(cc))),]),cc[9:1,])
head(dd, n = 3)
rownames(dd)[1]<-"Others"
head(dd, n = 3)
##再与metadata合并
bb<-merge(t(dd),dplyr::select(metadata,SampleID,Group),
          by.x = "row.names",by.y ="SampleID")
head(bb, n = 3)
##宽数据变长数据
kk<-tidyr::gather(bb,Phylum,Abundance,-c(Group,Row.names))
#kk<-tidyr::gather(bb,Species,Abundance,-c(Group,Row.names))
#将未注释到的Unassigned也改为Others,你也可以不改，有以下两种方式
kk$Phylum<-ifelse(kk$Phylum=='Unassigned','Others',kk$Phylum)#1
#kk$Phylum<-ifelse(kk$Species=='Unassigned','Others',kk$Species)#1
kk[kk$Phylum=='Unassigned','Phylum']='Others'               #2
##根据Group,Phylum分组运算
hh <- kk %>%
  group_by(Group,Phylum) %>%
  dplyr :: summarise(Abundance=sum(Abundance))
#hh <- kk %>%
#group_by(Group,Species) %>%
#dplyr :: summarise(Abundance=sum(Abundance))
#head(hh, n = 3)

##更改因子向量的levels
hh$Phylum = factor(hh$Phylum,order = T,levels = row.names(dd))
#hh$Species = factor(hh$Species,order = T,levels = row.names(dd))

yanse <-c('#BFBFBF', '#FF00FF','#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#B22222','gray')#要确保颜色够用，否则会报错
##按物种丰度排序好的堆积柱形图

p1 <- ggplot(hh,aes(x = Group,y = Abundance*100,#绝对丰度不乘以100
                    fill = Phylum,alluvium = Phylum, stratum = Phylum)) +
  geom_alluvium(aes(fill = Phylum),alpha = .5,width = 0.6) +
  geom_stratum(aes(fill = Phylum),width = 0.6)+
  scale_y_continuous(expand = c(0,0))
p2 <- p1 + ylab(label = "Relative Abundance(%)") + xlab(label = "")#Absolute abundance，绝对丰度；Relative Abundance相对丰度
p3 <- p2 +
  scale_fill_manual(values = yanse) +
  #scale_fill_aaas()+#ggsci配色,只提供10个
  #scale_fill_brewer(palette="Set1")+#只提供9个
  theme_bw()+ theme(panel.grid=element_blank()) +
  theme(panel.border = element_blank()) +
  theme(panel.background=element_rect(fill='transparent', color='black'),plot.margin = unit(c(3,5,1,1),"mm"))
p4 <- p3 + theme(axis.text.x=element_text(colour="black",size=28,face = "bold",vjust = 0.5,hjust = 0,angle = 315)) +
  theme(axis.text.y=element_text(colour = "black",size = 28)) +
  theme(axis.line = element_line(colour = "black"))+
  theme(axis.title.y = element_text(size = 28,face = "bold",margin = unit(c(0,1,0,1),"lines")))
p5 <- p4 + theme(legend.text = element_text(colour = "black",size = 28)) +
  theme(legend.title = element_text(size = 28,colour = "black")) +
  theme(text = element_text(family = "serif"))
p5




###进行处理间各物种非参数检验的多组比较
#数据整理与转换
head(bb,n = 3)
cc =dplyr::select(bb,Row.names,Group,everything(),-Others)
head(cc,n = 3)
dat <- gather(cc,Genus,v,-(Row.names:Group))
head(dat,n = 3)


##################保存图片
#ggsave("./p1.pdf", p1, width = 350, height = 200, units = "mm")
#ggsave("./p2.pdf", p2, width = 350, height = 200, units = "mm")

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


df3 <- ONE_LSD(dat,'Genus','Group','v')
head(df3)
p = ggplot(dat)+
  #geom_jitter(aes(x=Group,y=v,color=Group),width=0.25,shape=20,size=2,show.legend = F)+#散点
  stat_boxplot(aes(x=Group,y=v),geom="errorbar",size=0.6,width=0.25)+
  geom_boxplot(aes(x=Group,y=v,fill=Group))+#color非填充，fill填充
  #outlier.colour="red", outlier.shape=8,
  #outlier.size=4,outlier.stroke =0.5,
  #outlier.alpha = 0.5,outlier.fill = "red")+
  geom_text(data=df3,aes(x=Group,y=mean+1.3*std,label=lsd))+
  facet_wrap(.~Genus,scales = "free_y")+
  #stat_summary(aes(x=Group,y=v),fun="mean",geom="point",shape=23,size=1,fill="white")+
  labs(x='',y='Relative Abundance')+##Absolute abundance，绝对丰度；Relative Abundance相对丰度
  #ggprism::theme_prism()+
  #theme(axis.text.x = element_text(angle = 45))+
  #ggprism函数中的theme
  theme_bw()+
  guides(fill=guide_legend(title = "Treat"))+
  theme(text=element_text(size=16,  family="serif",face = "bold"),
        axis.text.x = element_text(size=10,  family="serif",colour = "black"))#serif在R中表示新罗马字体
p
ggsave("./差异分析_phylum0702.pdf", p5, width = 400, height = 400, units = "mm")
ggsave("./堆叠柱状图（属）.pdf", p5, width = 350, height = 400, units = "mm")
#pdf10*16
