setwd("G:/R/Rdata4_new/2021.12.20_group_network")
library(igraph)
C<-read.table("G:/R/Rdata4_new/2021.12.20_group_network/edge_C.csv",header = T,sep=",")
DL<-read.table("G:/R/Rdata4_new/2021.12.20_group_network/edge_DL.csv",header = T,sep=",")
NL<-read.table("G:/R/Rdata4_new/2021.12.20_group_network/edge_NL.csv",header = T,sep=",")
NP<-read.table("G:/R/Rdata4_new/2021.12.20_group_network/edge_NP.csv",header = T,sep=",")
NPDL<-read.table("G:/R/Rdata4_new/2021.12.20_group_network/edge_NPDL.csv",header = T,sep=",")
NPNL<-read.table("G:/R/Rdata4_new/2021.12.20_group_network/edge_NPNL.csv",header = T,sep=",")

f <- function(x){
  dat <- as.data.frame(x)
  ed <- dat[,1:2]
  ed <- as.matrix(ed)
  edgelist <- graph.edgelist(ed,directed=FALSE)
  adjacency <- get.adjacency(edgelist)
  adj_matrix <- as.matrix(adjacency)
  adj <- data.frame(adj_matrix,check.names=FALSE)
  return(adj)
}

set.seed(123)
adj_C <- f(C)
set.seed(123)
adj_DL <- f(DL)
set.seed(123)
adj_NL <- f(NL)
set.seed(123)
adj_NP <- f(NP)
set.seed(123)
adj_NPDL <- f(NPDL)
set.seed(123)
adj_NPNL <- f(NPNL)

write.table(adj_C,"1network_adj_matrix_C.txt",sep="\t",col.names = NA,quote=FALSE)
write.table(adj_DL,"1network_adj_matrix_DL.txt",sep="\t",col.names = NA,quote=FALSE)
write.table(adj_NL,"1network_adj_matrix_NL.txt",sep="\t",col.names = NA,quote=FALSE)
write.table(adj_NP,"1network_adj_matrix_NP.txt",sep="\t",col.names = NA,quote=FALSE)
write.table(adj_NPDL,"1network_adj_matrix_NPDL.txt",sep="\t",col.names = NA,quote=FALSE)
write.table(adj_NPNL,"1network_adj_matrix_NPNL.txt",sep="\t",col.names = NA,quote=FALSE)

adjacency_unweight_C <- read.delim('1network_adj_matrix_C.txt', row.names = 1, sep = '\t', check.names = FALSE)
adjacency_unweight_DL <- read.delim('1network_adj_matrix_DL.txt', row.names = 1, sep = '\t', check.names = FALSE)
adjacency_unweight_NL <- read.delim('1network_adj_matrix_NL.txt', row.names = 1, sep = '\t', check.names = FALSE)
adjacency_unweight_NP <- read.delim('1network_adj_matrix_NP.txt', row.names = 1, sep = '\t', check.names = FALSE)
adjacency_unweight_NPDL <- read.delim('1network_adj_matrix_NPDL.txt', row.names = 1, sep = '\t', check.names = FALSE)
adjacency_unweight_NPNL <- read.delim('1network_adj_matrix_NPNL.txt', row.names = 1, sep = '\t', check.names = FALSE)

#函数，得出网络属性
net <- function(adjacency_unweight){
  igraph = graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
  nodes_num <- length(V(igraph))
  edges_num <- length(E(igraph))
  average_degree <- mean(degree(igraph))
  nodes_connectivity <- vertex.connectivity(igraph)
  edges_connectivity <- edge.connectivity(igraph)
  average_path_length <- average.path.length(igraph, directed = FALSE)
  graph_diameter <- diameter(igraph, directed = FALSE)
  graph_density <- graph.density(igraph)
  clustering_coefficient <- transitivity(igraph)
  betweenness_centralization <- centralization.betweenness(igraph)$centralization
  degree_centralization <- centralization.degree(igraph)$centralization
  fc <- cluster_fast_greedy(igraph)
  modularity <- modularity(igraph, membership(fc))
  assortativity <- assortativity_degree(igraph, directed = FALSE)
  # betweenness <- betweenness(igraph,v=V(igraph),directed = F)
  igraph_character <- data.frame(
    nodes_num,    #节点数量（number of nodes）
    edges_num,    #边数量（number of edges）
    average_degree,    #平均度（average degree)
    nodes_connectivity,    #节点连通度（vertex connectivity）
    edges_connectivity,    #边连通度（edges connectivity）
    average_path_length,    #平均路径长度（average path length）
    graph_diameter,    #网络直径（diameter）
    graph_density,    #图密度（density）
    clustering_coefficient,    #聚类系数（clustering coefficient）
    betweenness_centralization,    #介数中心性（betweenness centralization)
    degree_centralization,    #度中心性
    modularity,    #模块性（modularity）
    assortativity#同配性
    #betweenness #介数
  )
  return(igraph_character)
}


set.seed(123)
igraph_character_C <- net(adjacency_unweight_C)
set.seed(123)
igraph_character_DL <- net(adjacency_unweight_DL)
set.seed(123)
igraph_character_NL <- net(adjacency_unweight_NL)
set.seed(123)
igraph_character_NP <- net(adjacency_unweight_NP)
set.seed(123)
igraph_character_NPDL <- net(adjacency_unweight_NPDL)
set.seed(123)
igraph_character_NPNL <- net(adjacency_unweight_NPNL)

#转置
hebing <- function(igraph_character){
  a <- data.frame(t(igraph_character))
  return(a)
} 

t_C <- hebing(igraph_character_C)
t_DL <- hebing(igraph_character_DL)
t_NL <- hebing(igraph_character_NL)
t_NP <- hebing(igraph_character_NP)
t_NPDL <- hebing(igraph_character_NPDL)
t_NPNL <- hebing(igraph_character_NPNL)

all <- data.frame(t_C,t_DL,t_NL,t_NP,t_NPDL,t_NPNL)
names(all) <- c("C","DL","NL","NP","NPDL","NPNL") #直接改名字

write.table(all, '1igraph_character_all_new.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)

#函数，得出度和中心性的列表附加介数
deg <- fuction(adjacency_unweight){
  igraph = graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
  degree <- degree(igraph)
  return(degree)
}

cen <- function(adjacency_unweight){
  igraph = graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
  centrality <- edge.betweenness(igraph)
  return(centrality)
}

# cen <- function(adjacency_unweight){
#   igraph = graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
#   eigenvector_centrality <- evcent(igraph)$vector
#   return(eigenvector_centrality)
# }


bet <- function(adjacency_unweight){
  igraph = graph_from_adjacency_matrix(as.matrix(adjacency_unweight), mode = 'undirected', weighted = NULL, diag = FALSE)
  betweenness <- betweenness(igraph,v=V(igraph),directed = F)
  return(betweenness)
}

#############################################
set.seed(123)
deg_C <- data.frame(deg(adjacency_unweight_C))
names(deg_C) <- c("Degree")
set.seed(123)
deg_DL <- data.frame(deg(adjacency_unweight_DL))
names(deg_DL) <- c("Degree")
set.seed(123)
deg_NL <- data.frame(deg(adjacency_unweight_NL))
names(deg_NL) <- c("Degree")
set.seed(123)
deg_NP <- data.frame(deg(adjacency_unweight_NP))
names(deg_NP) <- c("Degree")
set.seed(123)
deg_NPDL <- data.frame(deg(adjacency_unweight_NPDL))
names(deg_NPDL) <- c("Degree")
set.seed(123)
deg_NPNL <- data.frame(deg(adjacency_unweight_NPNL))
names(deg_NPNL) <- c("Degree")

write.table(deg_C, '1Degree_C.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(deg_DL, '1Degree_DL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(deg_NL, '1Degree_NL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(deg_NP, '1Degree_NP.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(deg_NPDL, '1Degree_NPDL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(deg_NPNL, '1Degree_NPNL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)

set.seed(123)
cen_C <- data.frame(cen(adjacency_unweight_C))
names(cen_C) <- c("Centrality")
set.seed(123)
cen_DL <- data.frame(cen(adjacency_unweight_DL))
names(cen_DL) <- c("Centrality")
set.seed(123)
cen_NL <- data.frame(cen(adjacency_unweight_NL))
names(cen_NL) <- c("Centrality")
set.seed(123)
cen_NP <- data.frame(cen(adjacency_unweight_NP))
names(cen_NP) <- c("Centrality")
cen_NPDL <- data.frame(cen(adjacency_unweight_NPDL))
names(cen_NPDL) <- c("Centrality")
cen_NPNL <- data.frame(cen(adjacency_unweight_NPNL))
names(cen_NPNL) <- c("Centrality")

write.table(cen_C, '1Centrality_C.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(cen_DL, '1Centrality_DL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(cen_NL, '1Centrality_NL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(cen_NP, '1Centrality_NP.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(cen_NPDL, '1Centrality_NPDL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(cen_NPNL, '1Centrality_NPNL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)

set.seed(123)
bet_C <- data.frame(bet(adjacency_unweight_C))
names(bet_C) <- c("Betweenness")
set.seed(123)
bet_DL <- data.frame(bet(adjacency_unweight_DL))
names(bet_DL) <- c("Betweenness")
set.seed(123)
bet_NL <- data.frame(bet(adjacency_unweight_NL))
names(bet_NL) <- c("Betweenness")
set.seed(123)
bet_NP <- data.frame(bet(adjacency_unweight_NP))
names(bet_NP) <- c("Betweenness")
set.seed(123)
bet_NPDL <- data.frame(bet(adjacency_unweight_NPDL))
names(bet_NPDL) <- c("Betweenness")
set.seed(123)
bet_NPNL <- data.frame(bet(adjacency_unweight_NPNL))
names(bet_NPNL) <- c("Betweenness")

write.table(bet_C, '1Betweenness_C.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(bet_DL, '1Betweenness_DL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(bet_NL, '1Betweenness_NL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(bet_NP, '1Betweenness_NP.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(bet_NPDL, '1Betweenness_NPDL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
write.table(bet_NPNL, '1Betweenness_NPNL.xls', sep = '\t', row.names = T, quote = FALSE,col.names = NA)
