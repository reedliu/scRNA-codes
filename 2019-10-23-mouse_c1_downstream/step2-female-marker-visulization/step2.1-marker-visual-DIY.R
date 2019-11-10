### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-DIY-marker基因可视化
### ---------------

#############################
# 准备表达矩阵和分群信息
#############################
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")
load('../female_rpkm.Rdata')

# 加载之前HCPC分群结果（代码在step1.1-tSNE-DIY.R）
cluster <- read.csv('../step1-female-RPKM-tSNE/step1.1-D-female_clustering.csv')
female_clustering=cluster[,2];names(female_clustering)=cluster[,1]
table(female_clustering)

#############################
# 准备marker基因矩阵
#############################
## 作者选择的14个marker基因
markerGenes <- c(
  "Nr2f1",
  "Nr2f2",
  "Maf",
  "Foxl2",
  "Rspo1",
  "Lgr5",
  "Bmp2",
  "Runx1",
  "Amhr2",
  "Kitl",
  "Fst",
  "Esr2",
  "Amh",
  "Ptges"
)
## 提取marker基因小表达矩阵
gene_subset <- as.matrix(log(females[rownames(females) %in% markerGenes,]+1))
dim(gene_subset)
gene_subset[1:4,1:4]
## 然后对这个小表达矩阵再细分，根据4个cluster的列名，也即是前面female_clustering=cluster[,2];names(female_clustering)=cluster[,1]这一步的目的
cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C1"])]
cl2_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C2"])]
cl3_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C3"])]
cl4_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C4"])]
## 重新再组合起来
heatmap_gene_subset <- cbind(
  cl1_gene_subset, 
  cl2_gene_subset,
  cl3_gene_subset,
  cl4_gene_subset
)
## 根据marker基因的顺序，重新排列这个矩阵
# 之前是这样
rownames(heatmap_gene_subset);markerGenes
# 得到marker基因在heatmap_gene_subset中的位置
match(markerGenes,rownames(heatmap_gene_subset))
# 然后就能提取出和marker基因顺序一样的heatmap_gene_subset
heatmap_gene_subset <- heatmap_gene_subset[match(markerGenes,rownames(heatmap_gene_subset)),]
# 之后是这样
rownames(heatmap_gene_subset);markerGenes

## 修改表达矩阵的列名，得到6个时间点信息
heatmap_female_stages <- sapply(strsplit(colnames(heatmap_gene_subset), "_"), `[`, 1)
table(heatmap_female_stages)

#############################
# 用包装好的pheatmap画热图
#############################
# https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2019-10-23-025853.png
# 看到H图中，列被分成了4栏，那么这个就是根据colbreaks来划分的。colbreaks的意思就是从哪里到哪里这是一块。当有多个分组又想画一个分割线的话，这个参数就很有用
colbreaks <- c(
  ncol(cl1_gene_subset),
  ncol(cl1_gene_subset)+ncol(cl2_gene_subset), 
  ncol(cl1_gene_subset)+ncol(cl2_gene_subset)+ncol(cl3_gene_subset)
)

# 然后就是上色，6个时间点和4个群使用自定义的颜色
cluster_color <- c(
  C1="#560047",
  C2="#a53bad", 
  C3="#eb6bac", 
  C4="#ffa8a0"
)

stage_color=c(
  E10.5="#2754b5", 
  E11.5="#8a00b0", 
  E12.5="#d20e0f", 
  E13.5="#f77f05", 
  E16.5="#f9db21",
  P6="#43f14b"
)

# 开始画热图
library(pheatmap)
png("step2.1-A-female_marker_heatmap.png")
plot_heatmap_2(
  heatmap_gene_subset, 
  female_clustering, 
  heatmap_female_stages, 
  rowbreaks, 
  colbreaks,
  cluster_color,
  stage_color
)
dev.off()

#############################
# 用包装好的ggboxplot画小提琴图
#############################
pdf("step2.1-B-markers-violin.pdf", width=10, height=22)
require(gridExtra)

# 每个基因的小提琴图都有4个cluster，对它们用不同的颜色
female_clusterPalette <- c(
  "#560047", 
  "#a53bad", 
  "#eb6bac", 
  "#ffa8a0"
)
# 每个基因做一个小提琴图，并用for循环保存在p这个列表中
p <- list()
for (genes in markerGenes) {
  p[[genes]] <- violin_gene_exp(
    genes, 
    females, 
    female_clustering, 
    female_clusterPalette
  )
}
# 最后组合起来，每列显示3张
do.call(grid.arrange,c(p, ncol=3))
dev.off()


#############################
# 用包装好的geom_point+geom_density_2d画等高线图
#############################
pdf("step2.1-C-markers-tSNE-density.pdf", width=16, height=28)
require(gridExtra)

load('../step1-female-RPKM-tSNE/step1.1-C-female_tsne.Rdata')
p <- list()
for (genes in markerGenes) {
  p[[genes]] <- tsne_gene_exp(
    female_tsne,
    genes, 
    females
  )
}

do.call(grid.arrange,c(p, ncol=4))
dev.off()
















