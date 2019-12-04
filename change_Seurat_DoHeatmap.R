### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-12-04
### Email: jieandze1314@gmail.com
### Title: 改造Seurat包装好的DoHeatmap函数
### ---------------
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)

# install.packages('Seurat')
library(Seurat)
library(stringr)   
library(dplyr)  

############################
# 原来Seurat做的热图
############################
load('sce_out_for_heatmap_all.Rdata')
sce
sce$seurat_clusters
top10 <- sce.markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
DoHeatmap(sce,top10$gene,size=3)

############################
# 接下来，进行DoHeatmap的改造
############################
# 第1步，得到表达矩阵
cts <- GetAssayData(sce, slot = "counts")
cts[1:4,1:4]
cts <- log10(cts + 1)

# 第2步，得到小的top10表达矩阵
# 将cluster排序，保证之后热图的顺序
new_cluster <- sort(sce$seurat_clusters)
head(new_cluster)
cts <- as.matrix(cts[top10$gene, names(new_cluster)])

# 第3步，做一个图例
library(pheatmap)
ac=data.frame(cluster=new_cluster)
rownames(ac)=colnames(mat)
pheatmap(cts,show_colnames =F,show_rownames = T,
         cluster_rows = F,
         cluster_cols = F,
         #gaps_col = 1:7,
         annotation_col=ac
         )

############################
# 再使用ComplexHeatmap画图
############################
#BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
# 列标题的颜色框
color <- rainbow(9)
names(color) <- levels(new_cluster)
top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color), # 设置填充色
                       labels = levels(new_cluster), 
                       labels_gp = gpar(cex = 0.5, col = "white"))) 

Heatmap(mat,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        show_column_names = FALSE,
        show_row_names = TRUE,
        column_split = new_cluster,
        heatmap_legend_param = list(
          title = "log10(count+1)",
          title_position = "leftcenter-rot"
        ),
        top_annotation = top_anno,
        column_title = NULL)







