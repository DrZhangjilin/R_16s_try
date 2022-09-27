setwd("G:/R/Rdata4_new")
##igraph 包计算网络模块
library(igraph)

#输入数据示例，邻接矩阵
#这是一个微生物互作网络，数值“1”表示微生物 OTU 之间存在互作，“0”表示无互作
adjacency_unweight <- read.delim('1network_adj_matrix_NPNL.txt', row.names = 1, sep = '\t', check.names = FALSE)
head(adjacency_unweight)[1:6]    #邻接矩阵类型的网络文件

#邻接矩阵 -> igraph 的邻接列表，获得非含权的无向网络
igraph <- graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
igraph    #igraph 的邻接列表

#计算节点度
V(igraph)$degree <- degree(igraph)

#模块划分，详情 ?cluster_fast_greedy，有多种模型
set.seed(123)
V(igraph)$modularity <- membership(cluster_fast_greedy(igraph))

#输出各节点（微生物 OTU）名称、节点度、及其所划分的模块的列表
nodes_list <- data.frame(
    nodes_id = V(igraph)$name, 
	degree = V(igraph)$degree, 
	modularity = V(igraph)$modularity
)
head(nodes_list)    #节点列表，包含节点名称、节点度、及其所划分的模块

write.table(nodes_list, '1nodes_listNPNL.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##计算模块内连通度（Zi）和模块间连通度（Pi）
source('zi_pi.r')
#######################################################################################################
#上述的邻接矩阵类型的网络文件
adjacency_unweight <- read.delim('1network_adj_matrix_NPNL.txt', row.names = 1, sep = '\t', check.names = FALSE)

#节点属性列表，包含节点所划分的模块
nodes_list <- read.delim('1nodes_listNPNL.txt', row.names = 1, sep = '\t', check.names = FALSE)

#两个文件的节点顺序要一致
nodes_list <- nodes_list[rownames(adjacency_unweight), ]

#计算模块内连通度（Zi）和模块间连通度（Pi）
#指定邻接矩阵、节点列表、节点列表中节点度和模块度的列名称
zi_pi <- zi.pi(nodes_list, adjacency_unweight, degree = 'degree', modularity_class = 'modularity')
head(zi_pi)

write.table(zi_pi, '1zi_pi_resultNPNL.txt', sep = '\t', row.names = FALSE, quote = FALSE)

##可再根据阈值对节点划分为 4 种类型，并作图展示其分布
library(ggplot2)
library(ggrepel)

zi_pi <- na.omit(zi_pi)   #NA 值最好去掉，不要当 0 处理
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Peripherals'
zi_pi[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Connectors'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities < 0.62),'type'] <- 'Module hubs'
zi_pi[which(zi_pi$within_module_connectivities > 2.5 & zi_pi$among_module_connectivities > 0.62),'type'] <- 'Network hubs'

zi_pi$nodes_id[which(zi_pi$within_module_connectivities < 2.5 & zi_pi$among_module_connectivities < 0.62)] <- ""

p <- ggplot(zi_pi, aes(among_module_connectivities, within_module_connectivities)) +
  geom_point(aes(color = type), alpha = 0.5, size = 6) +
  scale_color_manual(values = c('gray','red','blue','purple'), 
    limits = c('Peripherals', 'Connectors', 'Module hubs', 'Network hubs'))+
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), 
        panel.background = element_blank(), legend.key = element_blank()) +
 # geom_text_repel(aes(color=factor(type), label=nodes_id),max.overlaps = 40,size=2)+
  geom_label_repel(aes(fill=type,label=nodes_id),
                   show.legend=FALSE,fontface="bold", color="white", 
                   box.padding=unit(0.35, "lines"), 
                   segment.colour = "grey50",size=2,max.overlaps = 200)+#加标签
  labs(x = 'Among-module connectivities', y = 'Within-module connectivities', color = '') +
  geom_vline(xintercept = 0.62) +
  geom_hline(yintercept = 2.5)+
  theme_bw()+#白底灰线背景
  theme(legend.key = element_blank())+
  theme(text=element_text(size=28,  family="serif",face = "bold"),
        axis.text.x = element_text(size=18,  family="serif",colour = "black"),
        axis.text.y = element_text(size=18,  family="serif",colour = "black"))#serif在R中表示新罗马字体
  #scale_x_continuous(limits = c(0,1))+ 
  #scale_y_continuous(limits = c(-3,4))
p
ggsave('./1NPNL_zipi.pdf',p,width=250,height=200,units="mm")
