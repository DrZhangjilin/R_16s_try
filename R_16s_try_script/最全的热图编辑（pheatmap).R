#setwd("G:/R/Rdata2")
setwd("G:/R/Rdata3")
library(pheatmap) # 加载包
library(ggplot2) # 加载包
data <- read.table("otu_table_L2_3.xls",  # 读取的数据文件名称，这里文件是放在工作目录下
                   header=T, # 数据集第一行为变量名
                   row.names=1, # 第一列为行名
                   sep="\t") # 指定分隔符号
dim(data) # 查看变量有多少行多少列
head(data) # 查看数据集前六行

p <- pheatmap(data)#绘制基本款热图
p <- pheatmap(data, scale="row")#设置标准化方向scale，对其横向标准化


#设置边框为白色，去掉横向、纵向聚类。
p <- pheatmap(data, scale="row",
              border="white", # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F)


#去掉横纵坐标中的id
p <- pheatmap(data,scale="row",
              border="white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = F, #去掉横、纵坐标id
              show_colnames = F)


# 不显示右上角图例
p <- pheatmap(data,scale="row",
              border="white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = F, #去掉横、纵坐标id
              show_colnames = F,
              legend = F) # 去掉图例


#设置图例范围
p <- pheatmap(data,scale="row",
              border="white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = F, #去掉横、纵坐标id
              show_colnames = F,
              legend = T, # 添加图例
              legend_breaks=c(-1.5,0,1.5)) # 设置图例范围
# 也可以设置legend_breaks=c(-2,0,2)试试


#设置图中字的大小，使用fondsize参数来设置。
p <- pheatmap(data,scale="row",
              border="white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = F, # 去掉图例
              fontsize = 8)  # 设置字体大小
# 也可以设置其他大小试试。


#分别指定横向和纵向字体大小，使用fontsize_row和fontsize_col参数设置。
p <- pheatmap(data,scale="row",
              border="white",  # 设置边框为白色
              cluster_cols = F, # 去掉横向、纵向聚类
              cluster_rows = F,
              show_rownames = T, #显示横、纵坐标id
              show_colnames = T,
              legend = T, # 显示图例
              fontsize_row = 12, # 分别设置横向和纵向字体大小
              fontsize_col = 16)  


#设置聚类的距离类型，使用clustering_distance_rows参数指定，分为如下几类：correlation，euclidean，maximum，manhattan，canberra，binary，minkowski。
p <- pheatmap(data,scale="row",
              border="white", # 设置边框为白色
              cluster_cols = T, # 显示横向、纵向聚类
              cluster_rows = T,
              clustering_distance_rows = "correlation", # 设置聚类的距离类型
              treeheight_col = 50, # 分别设置横、纵向聚类树高
              treeheight_row = 45)
#调整聚类的方法，使用clustering_method参数指定，可选有'ward', 'ward.D', 'ward.D2', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.
p <- pheatmap(data,scale="row",
              border="white", # 设置边框为白色
              cluster_cols = T, # 显示横向、纵向聚类
              cluster_rows = T,
              clustering_distance_rows = "euclidean", # 设置聚类的距离类型
              clustering_method="complete", # 设置聚类方法
              treeheight_col = 50, # 分别设置横、纵向聚类树高
              treeheight_row = 45)


#设置分组标签的角度，可以使用参数angle_col指定，可选有270、0、45、90、315等。
p <- pheatmap(data,scale="row",
              angle_col = 45, # 设置显示角度
              clustering_distance_rows = "euclidean",
              clustering_method="complete",
              border="white",  
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)


#给图形增加标题，可以使用main参数指定。
p <- pheatmap(data, scale="row", border="white",
              main="Treat", # 设置图形标题
              angle_col = 315,
              clustering_distance_rows = "euclidean",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)


#分别调整热图方块宽度和高度，可以使用cellwidth和cellheight参数指定。
p <- pheatmap(data, scale="row", border="white",
              cellwidth = 50,cellheight = 7.5, # 设置热图方块宽度和高度
              main="Treat", angle_col = 0,
              clustering_distance_rows = "euclidean",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)


#根据热图聚类对其进行区块儿划分，可以使用cutree_cols和cutree_rows参数指定
p <- pheatmap(data,scale="row",
              cutree_cols = 6, cutree_rows =5, # 列划为6块，行为5块
              main="Treat",angle_col = 0,
              clustering_distance_rows = "euclidean",
              clustering_method="complete",border="white",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)
#注：此处我去掉了热图块儿的大小；cellwidth = 8,cellheight = 6

#在上图基础上增加边缘线。
p <- pheatmap(data,scale="row", border="#8B0A50",
              cutree_cols = 6, cutree_rows =5,
              main="Treat",angle_col = 0,
              clustering_distance_rows = "euclidean",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)


###热图上是否展示数值，大小和颜色，大小以及数值展示类型，可以使用display_numbers、fontsize_number、number_color、number_format等参数设置。
#使用display_numbers参数指定是否显示数值。
p <- pheatmap(data,scale="row",border="#8B0A50",
              display_numbers = T, # 热图上显示数值
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)

###使用fontsize_number参数指定数值的显示大小。
p <- pheatmap(data,scale="row",border="#8B0A50",
              fontsize_number = 10, display_numbers = T,
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)

###使用number_color参数指定数值的颜色。
p <- pheatmap(data,scale="row",border="#8B0A50",
              number_color="red",
              fontsize_number = 10,display_numbers = T,
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)

###使用number_format参数指定数值显示类型，下图显示为科学计数法。
p <- pheatmap(data,scale="row", border="#8B0A50", number_color="red",
              number_format="%.2e",
              fontsize_number = 10,display_numbers = T,
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)

####对热图方块儿进行标记；display_numbers，如果该值大于1，则为+，否则为-
p <- pheatmap(data,scale="row",
              number_color="red",number_format="%.2e",
              border="#8B0A50",
              fontsize_number = 16,
              display_numbers = matrix(ifelse(data > 1, "+", "-"), nrow(data)),
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)
# 也可以设置display_numbers = matrix(ifelse(data > 2, "++", "-"), nrow(data))试试




#换色卡
color = colorRampPalette(c("navy","white","firebrick3"))(100)


#构建横（纵）向分组信息
map <- read.csv("otu_table_L2_3_map.csv",  # 读取的数据文件名称，这里文件是放在工作目录下
                header=T, # 数据集第一行为变量名
                row.names=1, # 第一列为行名
                stringsAsFactors = F)
#annotation_col=xx标注纵向分组信息，annotation_row=xx标注横向分组信息。
#annotation_names_row 是否显示行注释的名称。annotation_names_col是否显示列注释的名称。
p <- pheatmap(data, color=color,
              annotation_row=map,annotation_names_row=F,
              scale="row", border="#8B0A50",
              cutree_cols = 6, cutree_rows =5,
              main="Treat",angle_col = 0,
              clustering_distance_rows = "euclidean",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20,
              legend_breaks=c(-1.5,0,1.5))







####构建纵向和横向分组信息(有6个处理，分别是：盐、干旱和热应激)；以及时间：0-3day，对3类基因21个基因进行分组，分别是："WRKY", "AP2", "YABBY"。
####构建纵向分组信息
annotation_col = data.frame(Deal_with = factor(rep(c("Salt", "Drought","Heat"), 3)),
                            Day=factor(rep(c("Day1", "Day2","Day3"), 3)))
rownames(annotation_col)
colnames(data)
rownames(annotation_col) <- colnames(data)
head(annotation_col)
pheatmap(data, annotation_col = annotation_col)
p <- pheatmap(data, scale="row",
              annotation_col = annotation_col,
              number_color="red",number_format="%.2e",
              border="#8B0A50",fontsize_number = 8,
              display_numbers = matrix(ifelse(data > 2, "++", "-"), nrow(data)),
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)

####构建横向分组信息
annotation_row = data.frame(GeneClass = factor(rep(c("WRKY", "AP2", "YABBY"),7)))
rownames(annotation_row) <- rownames(data)
head(annotation_row)
pheatmap(data, annotation_row =annotation_row)
p <- pheatmap(data, scale="row",
              number_color="red", 
              annotation_row =annotation_row,
              number_format="%.2e",
              border="#8B0A50",
              fontsize_number = 16, 
              display_numbers = matrix(ifelse(data > 2, "++", "-"), nrow(data)),
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)

####共同组合上述二者
p <- pheatmap(data,
              annotation_col = annotation_col,
              annotation_row = annotation_row)
p <- pheatmap(data,scale="row",
              number_color="red",
              annotation_col = annotation_col,
              annotation_row = annotation_row,
              number_format="%.2e",border="#8B0A50",
              fontsize_number = 15,#热图数字大小_numbers = matrix(ifelse(data > 2, "++", "-"),nrow(data)),
              cutree_cols = 3,cutree_rows =4,
              main="Gene1",angle_col = 0,
              clustering_distance_rows = "minkowski",
              clustering_method="complete",
              cluster_cols = T,treeheight_col = 20,
              cluster_rows = T,treeheight_row = 20)
p