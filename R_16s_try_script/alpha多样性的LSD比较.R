library(ggplot2)
library(ggprism)
library(plyr)
library(reshape2)
library(tidyverse)
library(vegan)
setwd("G:/R/Rdata2")

dat <- read.table('./alpha_div.xls',row.names = 1,header = T,stringsAsFactors = F)#读入α多样性数据
head(dat, n = 3)

rname <- row.names(dat)
rname <- t(rname)
rname <- t(rname)
dat$Row.names <- rname

dat <- melt(dat,id.vars = -c(1:5),variable.name = 'diversity')#宽数据变长数据
head(dat, n = 3)

dat$diversity <- as.factor(dat$diversity)#将α列转化成因子
names(dat)[4] <- 'v'#给value重新赋列名
head(dat, n = 3)


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

df <- ONE_LSD(dat,'diversity','Group','v')
head(df)


PMCMR_compare1 <- function(data,group,compare,value){
  library(multcompView)
  library(multcomp)
  library(PMCMRplus)
  library(PMCMR)
  a <- data.frame(stringsAsFactors = F)
  type <- unique(data[,group])
  for (i in type)
  {
    g1=compare
    sub_dat <- data[data[,group]==i,]
    names(sub_dat)[names(sub_dat)==compare] <- 'g1'
    names(sub_dat)[names(sub_dat)==value] <- 'value'
    sub_dat$g1 <- factor(sub_dat$g1)
    options(warn = -1)

    k <- PMCMRplus::kwAllPairsNemenyiTest(value ~ g1,data=sub_dat)
    n <- as.data.frame(k$p.value)
    h <- n %>%
      mutate(compare=rownames(n)) %>%
      gather(group,p,-compare,na.rm = TRUE) %>%
      unite(compare,group,col="G",sep="-")
    dif <- h$p
    names(dif) <- h$G
    dif
    difL <- multcompLetters(dif)
    K.labels <- data.frame(difL['Letters'], stringsAsFactors = FALSE)
    K.labels$compare = rownames(K.labels)
    K.labels$type <- i

    mean_sd <- merge(aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=sd),
                     aggregate(sub_dat[['value']],by=list(sub_dat[,'g1']),FUN=mean),by="Group.1"
    )
    names(mean_sd) <- c('compare','std','mean')
    a <- rbind(a,merge(mean_sd,K.labels,by='compare'))
  }
  names(a) <- c(compare,'std','mean','Letters',group)
  return(a)
}

##################################
df <- PMCMR_compare1(dat,'diversity','Group','v')
df

order_x=c("PD_whole_tree","chao1","observed_species","shannon","simpson")
df$diversity = factor(df$diversity,
                        levels = order_x )

order_y=c("C","LA","LR","PR","LAPR","LRPR")
dat$Group = factor(dat$Group,
                      levels = order_y )


p = ggplot(dat)+
  geom_jitter(aes(x=Group,y=v,color=Group),width=0.25,shape=20,size=1,show.legend = F)+#散点
  stat_boxplot(aes(x=Group,y=v),geom="errorbar",size=0.6,width=0.25)+
  geom_boxplot(aes(x=Group,y=v,fill=Group))+#color非填充，fill填充
  #outlier.colour="red", outlier.shape=8,
  #outlier.size=4,outlier.stroke =0.5,
  #outlier.alpha = 0.5,outlier.fill = "red")+
  facet_wrap(.~diversity,scales = "free_y",ncol = 5)+#分面
  #scale_y_continuous(expand = c(0.1,0))+
  geom_text(data=df,aes(x=Group,y=mean+1.4*std,label=lsd),size=3,family="serif")+
  #stat_summary(aes(x=treat,y=v),fun="mean",geom="point",shape=23,size=1,fill="white")+
  labs(x='',y='Alpha diversity')+
  #ggprism::theme_prism()+
  #theme(axis.text.x = element_text(angle = 45))+
  #ggprism函数中的theme
  theme_bw()+
  #guides(fill=guide_legend(title = "Treat"))+
  guides(fill="none")+
  theme(text=element_text(size=10,  family="serif",face = "bold"),
        axis.text.x = element_text(size=8,  family="serif",colour = "black",angle=-45),
        axis.text.y = element_text(size=8,  family="serif",colour = "black"))
p
ggsave('./alpha_div_20220702.pdf',p,width=200,height=50,units="mm")
