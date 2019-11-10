### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-11-10
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-聚类后绘制热图
### ---------------

rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

load('step5.2-for_heatmap.Rdata')
load('step5.2-gene_clustering.Rdata')
# 3个谱系（2个特有+1个共有）
table(cellLin)
# 4种细胞类型
table(cellType)
# 17个基因分类结果
table(clustering[,1])

#############################
# 热图准备之--配置基因分组颜色
#############################
gene_cluster_palette <- c(
  '#a6cee3',
  '#1f78b4',
  '#b2df8a',
  '#33a02c',
  '#fb9a99',
  '#e31a1c',
  '#fdbf6f',
  '#ff7f00',
  '#cab2d6',
  '#6a3d9a',
  '#ffff99',
  '#b15928', 
  '#49beaa', 
  '#611c35', 
  '#2708a0',
  '#fccde5',
  '#bc80bd'
)
gene_cluster_colors <- gene_cluster_palette[1:max(clusters)]
names(gene_cluster_colors) <- 1:max(clusters)

#############################
# 热图准备之--配置行、列注释信息
#############################
# 行注释：每个基因属于哪个组
annotation_row <- data.frame(clustering=clustering)

# 列注释：三种信息cell lineage, cell cluster type, cell stage
annotation_col <- data.frame(
  cellLineages=cellLin,
  cellType=cellType,
  Stages=sapply(strsplit(colnames(data_heatmap), "_"), `[`, 1)
)
rownames(annotation_col) <- colnames(data_heatmap)

# 3个cell lineages颜色
cellLinCol <- c(
  "#3b3561", 
  "#c8c8c8", 
  "#ff6663"
)
names(cellLinCol) <- unique(cellLin)

# 4个cell clusters颜色
cellTypeCol <- c(
  C2="#a53bad", 
  C1="#560047", 
  C3="#eb6bac", 
  C4="#ffa8a0"
)
names(cellTypeCol) <- unique(cellType)

#############################
# 热图准备之--把三种注释信息颜色放在一起
#############################
# 3个谱系（2个特有+1个共有）+4种细胞类型 +17个基因分类结果
annotation_colors <- list(
  cellType=cellTypeCol,
  cellLineages=cellLinCol,
  clustering=gene_cluster_colors,
  Stages=c(
    E10.5="#2754b5", 
    E11.5="#8a00b0", 
    E12.5="#d20e0f", 
    E13.5="#f77f05", 
    E16.5="#f9db21",
    P6="#43f14b"
  )
)

# 调画板
cold <- colorRampPalette(c('#f7fcf0','#41b6c4','#253494','#081d58','#081d58'))
warm <- colorRampPalette(c('#ffffb2','#fecc5c','#e31a1c','#800026','#800026'))
mypalette <- c(rev(cold(21)), warm(20))
breaksList = seq(-2.2, 2.5, by = 0.2)

#############################
# 热图
#############################
library(pheatmap)
tiff(file="step5.3-A-female_heatmap_DE_genes_k_17_pval_005.tiff", 
     res = 300, height = 21, width = 18, units = 'cm')
gene_clustering <- pheatmap(
  data_heatmap, 
  scale="row",
  gaps_col=length(cellType_L2),
  show_colnames=FALSE, 
  show_rownames=FALSE, 
  cluster_cols=FALSE,
  clustering_method="ward.D",
  annotation_row=annotation_row,
  annotation_col=annotation_col,
  annotation_colors=annotation_colors,
  cutree_rows=17, 
  annotation_names_row=FALSE,
  color=mypalette
)
dev.off()

