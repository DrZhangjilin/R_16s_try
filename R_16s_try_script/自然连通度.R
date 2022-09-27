setwd("G:/R/Rdata4_new")

library(igraph)
library(ggplot2)
library(ggpmisc)
a=read.csv("G:/R/Rdata4_new/edge_C.csv")
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
###remove node###
nc.attack <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), decreasing=T)
  sapply(1:round(vcount(ig))*0.50,#随机抽取百分之80的点
         function(i) {
           ind <- hubord[1:i]
           tmp <- delete_vertices(ig, V(ig)$name[ind])
           natcon(tmp)
         })
}
nc=nc.attack(ig)
#nc=read.table("natural connectivity_nodesC.txt",sep="\t")
#nc=as.matrix(nc)
plot(seq(0,0.8,len=length(nc)), nc, type='l', ylim=c(0,max(nc)),xlab="Proportion of removed nodes", ylab="natural connectivity")
hist(diff(nc), breaks=30, col='red', ylim=c(0,100), main="Fragility rates")
write.table(nc,file="natural connectivity_nodesC1.txt",sep="\t", quote=FALSE)

a=read.csv("G:/R/Rdata4_new/edge_DL.csv")
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
###remove node###
nc.attack <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), decreasing=T)
  sapply(1:round(vcount(ig))*0.50,#随机抽取百分之80的点
         function(i) {
           ind <- hubord[1:i]
           tmp <- delete_vertices(ig, V(ig)$name[ind])
           natcon(tmp)
         })
}
nc=nc.attack(ig)
#nc=read.table("natural connectivity_nodes.txt",sep="\t")
#nc=as.matrix(nc)
plot(seq(0,0.8,len=length(nc)), nc, type='l', ylim=c(0,max(nc)),xlab="Proportion of removed nodes", ylab="natural connectivity")
hist(diff(nc), breaks=30, col='red', ylim=c(0,100), main="Fragility rates")
write.table(nc,file="natural connectivity_nodesDL1.txt",sep="\t", quote=FALSE)

a=read.csv("G:/R/Rdata4_new/edge_NL.csv")
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
###remove node###
nc.attack <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), decreasing=T)
  sapply(1:round(vcount(ig))*0.50,#随机抽取百分之80的点
         function(i) {
           ind <- hubord[1:i]
           tmp <- delete_vertices(ig, V(ig)$name[ind])
           natcon(tmp)
         })
}
nc=nc.attack(ig)
#nc=read.table("natural connectivity_nodes.txt",sep="\t")
#nc=as.matrix(nc)
plot(seq(0,0.8,len=length(nc)), nc, type='l', ylim=c(0,max(nc)),xlab="Proportion of removed nodes", ylab="natural connectivity")
hist(diff(nc), breaks=30, col='red', ylim=c(0,100), main="Fragility rates")
write.table(nc,file="natural connectivity_nodesNL1.txt",sep="\t", quote=FALSE)

a=read.csv("G:/R/Rdata4_new/edge_NP.csv")
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
###remove node###
nc.attack <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=T)
  # hubord <- order(rank(betweenness(ig)), decreasing=T)
  sapply(1:round(vcount(ig))*0.50,#随机抽取百分之80的点
         function(i) {
           ind <- hubord[1:i]
           tmp <- delete_vertices(ig, V(ig)$name[ind])
           natcon(tmp)
         })
}
nc=nc.attack(ig)
#nc=read.table("natural connectivity_nodes.txt",sep="\t")
#nc=as.matrix(nc)
plot(seq(0,0.8,len=length(nc)), nc, type='l', ylim=c(0,max(nc)),xlab="Proportion of removed nodes", ylab="natural connectivity")
hist(diff(nc), breaks=30, col='red', ylim=c(0,100), main="Fragility rates")
write.table(nc,file="natural connectivity_nodesNP1.txt",sep="\t", quote=FALSE)

a=read.csv("G:/R/Rdata4_new/edge_NPDL.csv")
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
###remove node###
nc.attack <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=T)
  # hubord <- order(rank(betweenness(ig)), decreasing=T)
  sapply(1:round(vcount(ig))*0.50,#随机抽取百分之80的点
         function(i) {
           ind <- hubord[1:i]
           tmp <- delete_vertices(ig, V(ig)$name[ind])
           natcon(tmp)
         })
}
nc=nc.attack(ig)
#nc=read.table("natural connectivity_nodes.txt",sep="\t")
#nc=as.matrix(nc)
plot(seq(0,0.8,len=length(nc)), nc, type='l', ylim=c(0,max(nc)),xlab="Proportion of removed nodes", ylab="natural connectivity")
hist(diff(nc), breaks=30, col='red', ylim=c(0,100), main="Fragility rates")
write.table(nc,file="natural connectivity_nodesNPDL1.txt",sep="\t", quote=FALSE)

a=read.csv("G:/R/Rdata4_new/edge_NPNL.csv")
b=a[,c(1:2)]
b=t(b)
b=t(b)
ig=graph_from_edgelist(b,directed=FALSE)
natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
###remove node###
nc.attack <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=T)
  #hubord <- order(rank(betweenness(ig)), decreasing=T)
  sapply(1:round(vcount(ig))*0.50,#随机抽取百分之80的点
         function(i) {
           ind <- hubord[1:i]
           tmp <- delete_vertices(ig, V(ig)$name[ind])
           natcon(tmp)
         })
}
nc=nc.attack(ig)
#nc=read.table("natural connectivity_nodes.txt",sep="\t")
#nc=as.matrix(nc)
plot(seq(0,0.8,len=length(nc)), nc, type='l', ylim=c(0,max(nc)),xlab="Proportion of removed nodes", ylab="natural connectivity")
hist(diff(nc), breaks=30, col='red', ylim=c(0,100), main="Fragility rates")
write.table(nc,file="natural connectivity_nodesNPNL1.txt",sep="\t", quote=FALSE)


#################################################################################################################################################
#################################################################################################################################################

ncC=read.table("natural connectivity_nodesC1.txt",sep="\t")
ncC=as.matrix(ncC)
ncDL=read.table("natural connectivity_nodesDL1.txt",sep="\t")
ncDL=as.matrix(ncDL)
ncNL=read.table("natural connectivity_nodesNL1.txt",sep="\t")
ncNL=as.matrix(ncNL)
ncNP=read.table("natural connectivity_nodesNP1.txt",sep="\t")
ncNP=as.matrix(ncNP)
ncNPDL=read.table("natural connectivity_nodesNPDL1.txt",sep="\t")
ncNPDL=as.matrix(ncNPDL)
ncNPNL=read.table("natural connectivity_nodesNPNL1.txt",sep="\t")
ncNPNL=as.matrix(ncNPNL)
#有问题
dat1 <- data.frame(natural_connectivity = c(ncC,ncDL,ncNL,ncNP,ncNPDL,ncNPNL),
                   remove_node = c(1:156,1:129,1:150,1:183,1:156,1:152),
                   network = c(rep('C', 156),
                               rep('DL', 129),
                               rep('NL', 150),
                               rep('NP', 183),
                               rep('NPDL',156),
                               rep('NPNL',152)))
p <- ggplot(dat1, aes(remove_node, natural_connectivity, color = network)) +
  #geom_smooth(se = F,span=0.2)+#拟合的曲线，se指定是否存在置信区间,span=指定曲线光滑程度
  geom_smooth(method = "lm",formula = y~x,se=TRUE,show.legend = F)+
  stat_poly_eq(aes(label = paste(..eq.label.., ..rr.label.., stat(p.value.label), sep = "~`,`~")),
               formula = y~x, parse = TRUE, label.x.npc = 'right', label.y.npc = 'top', size = 10,
               family="serif",fontface="bold")+
  geom_point(size=2)+
  #geom_line(size=2)+#折线图
  facet_wrap(.~network,scales = "free_y",ncol = 3)+#分面
  labs(x="Remove nodes",y="Natural connectivity")+
  guides(color=guide_legend(title = "Treat",override.aes = list(size=2.5)))+#override.aes = list(size=2.5)修改图例中点的大小
  theme_bw()+#白底灰线背景
  theme(legend.key = element_blank())+
  theme(text=element_text(size=48,  family="serif",face = "bold"),
        axis.text.x = element_text(size=36,  family="serif",colour = "black"),
        axis.text.y = element_text(size=36,  family="serif",colour = "black"))#serif在R中表示新罗马字体
p
#ggsave("./Natural connectivity0.5he.pdf", p, width = 800, height = 400, units = "mm")
ggsave("./连通性9.pdf", p, width = 800, height = 400, units = "mm")
#pdf7*8