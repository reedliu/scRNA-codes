### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2020-01-03
### Email: jieandze1314@gmail.com
### Title: 改造inferCNV的热图
### ---------------

rm(list = ls())
options(stringsAsFactors = F)

#--- 获得热图数据 ------
# 比如要画：infercnv.preliminary.png(https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2020-01-03-024640.png)，就需要用到preliminary.infercnv_obj
library(infercnv)
infercnv_obj = readRDS('output_dir/preliminary.infercnv_obj')
# 提取：处理后的表达矩阵
expr <- infercnv_obj@expr.data
dim(expr)
expr[1:4,1:4]

# 提取：表达矩阵样本中正常细胞的位置
(normal_loc <- infercnv_obj@reference_grouped_cell_indices)
normal_loc <- normal_loc$WT
# 提取：表达矩阵样本中肿瘤细胞的位置
(tumor_loc <- infercnv_obj@observation_grouped_cell_indices)
tumor_loc <- tumor_loc$`1`


#--- 如何绘制热图 ------
# 通过观察：infercnv.preliminary.png，发现纵坐标是样本，并进行了聚类；横坐标是染色体位置（chr1...），记录了各个基因的位置
# 图片的解释：https://github.com/broadinstitute/inferCNV/wiki/Interpreting-the-figure

## Step1: 将基因名与染色体位置对应
gn <- rownames(expr)
length(gn)
head(gn)
# 加载基因位置信息文件(这个存储的是排序后的)
geneFile <- read.table('infercnv-heatmap/表达矩阵数据/geneFile.txt')
geneFile[1:4,1:4]
length(geneFile$V1);length(gn) #下面需要根据gn对geneFile取小

sub_geneFile <-  geneFile[geneFile$V1%in%gn,]
dim(sub_geneFile)

expr[1:4,1:4]
head(sub_geneFile)
identical(rownames(expr),sub_geneFile$V1)

## Step2: 拆分矩阵
# 整体分成两部分：上面的热图是正常细胞，下面是肿瘤细胞，并且我们知道了各自的位置，就能先把各自的小表达矩阵提取出来

norm_expr <- expr[,normal_loc]
norm_expr$chr <- as.factor(sub_geneFile$V2) #最后加上一列：对应的chr信息
table(norm_expr$chr) #这个信息就与横坐标的间隔对应，chr间隔越大表示其中包含的基因越多

tumor_expr <- expr[,tumor_loc]
tumor_expr$chr <- as.factor(sub_geneFile$V2)
dim(tumor_expr)


## Step3-1: 开始画图=》第一次尝试【先画tumor的】
library(ComplexHeatmap)
# 列标题的颜色框
color <- rainbow(length(unique(tumor_expr$chr)))
# 设定不同chr出现的顺序
new_cluster <- tumor_expr$chr
names(color) <- levels(new_cluster)

top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color), # 设置填充色
                       labels = levels(new_cluster), 
                       labels_gp = gpar(cex = 0.5, col = "white"))) 

# 需要画：行是样本，列是基因
pdf("test1-heatmap.pdf")
Heatmap(t(as.matrix(tumor_expr[,-ncol(tumor_expr)])),
        cluster_rows = T,
        cluster_columns = F,
        show_column_names = F,
        show_row_names = F,
        column_split = new_cluster,
        heatmap_legend_param = list(
          title = "Expression",
          title_position = "leftcenter-rot"
        ),
        top_annotation = top_color,
        column_title = NULL)
dev.off()

## Step3-2: 画图=》第二次尝试【先画tumor的】

# 调整配色
library("RColorBrewer")
# 来自函数：infercnv::plot_cnv
get_group_color_palette <- function () {
  return(colorRampPalette(RColorBrewer::brewer.pal(12, "Set3")))
}
color <- get_group_color_palette()(length(unique(tumor_expr$chr)))
top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color), # 设置填充色
                       labels = levels(new_cluster), 
                       labels_gp = gpar(cex = 0.9, col = "black"))) 

n <- t(tumor_expr[,-ncol(tumor_expr)])
n[n>1.6]=1.6
n[n<0.4]=0.4

pdf("test2-heatmap.pdf",width = 15,height = 10)
ht_tumor = Heatmap(as.matrix(n),
                  cluster_rows = T,
                  cluster_columns = F,
                  show_column_names = F,
                  show_row_names = F,
                  column_split = new_cluster,
                  heatmap_legend_param = list(
                    title = "Modified Expression",
                    title_position = "leftcenter-rot", # 图例标题位置
                    at=c(0.4,1.6), #图例范围
                    legend_height = unit(3, "cm") #图例长度
                  ),
                  top_annotation = top_color,
                  row_title = "Observations (Cells)",
                  row_title_side = c("right"),
                  column_title = "Genomic Region",
                  column_title_side = c("bottom"))
draw(ht_tumor, heatmap_legend_side = "left") # 图例位置

dev.off()

## Step3-3: 画图=》第三次尝试【组合normal和tumor】
m <- t(norm_expr[,-ncol(norm_expr)])
m[m>1.6]=1.6
m[m<0.4]=0.4

pdf("test3-heatmap.pdf",width = 25,height = 30)
ht_normal = Heatmap(as.matrix(m),
                   cluster_rows = T,
                   cluster_columns = F,
                   show_column_names = F,
                   show_row_names = F,
                   column_split = new_cluster,
                   row_title = "References (Cells)",
                   row_title_side = c("right"),
                   row_title_rot = 90,
                   row_title_gp = gpar(fontsize = 25),
                   column_title = NULL, 
                   heatmap_legend_param = list(
                     title = "Modified Expression",
                     title_position = "leftcenter-rot", # 图例标题位置
                     title_gp = gpar(fontsize = 20),# 图例标题大小
                     at=c(0.4,1.6), #图例范围
                     legend_height = unit(6, "cm")),#图例长度
                   width = 20, height = 5
                   ) 

ht_tumor = Heatmap(as.matrix(n),
                   cluster_rows = T,
                   cluster_columns = F,
                   show_column_names = F,
                   show_row_names = F,
                   column_split = new_cluster,
                   show_heatmap_legend=F,
                   top_annotation = top_color,
                   row_title = "Observations (Cells)",
                   row_title_side = c("right"),
                   row_title_rot = 90,
                   row_title_gp = gpar(fontsize = 25),
                   column_title = "Genomic Region",
                   column_title_side = c("bottom"),
                   column_title_gp = gpar(fontsize = 25),
                   width = 20, height = 10,
                   heatmap_height = 15)

# 设置图例位置，并竖直排列
draw(ht_normal %v% ht_tumor)

dev.off()


