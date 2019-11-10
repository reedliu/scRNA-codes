### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-上一步的差异基因聚类
### ---------------

rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

#############################
# 加载数据
#############################
# RPKM进行可视化，count矩阵进行差异分析
load('../female_rpkm.Rdata')
load('../female_count.Rdata')
females[1:4,1:4]
female_count[1:4,1:4]
# 谱系推断结果(细胞谱系归一化成为百分比)
load('../step4-psudotime/step4.3-female-psudotime-percent.Rdata')

# 6个发育时期获取
head(colnames(female_count))
female_stages <- sapply(strsplit(colnames(female_count), "_"), `[`, 1)
names(female_stages) <- colnames(female_count)
table(female_stages)
# 4个cluster获取
cluster <- read.csv('../step1-female-RPKM-tSNE/step1.1-D-female_clustering.csv')
female_clustering=cluster[,2];names(female_clustering)=cluster[,1]
table(female_clustering)

# 上一步的差异基因
load('step5.1-lineage_sig_gene.Rdata')
head(female_lineage1_sig_gene_pseudoT)
head(female_lineage2_sig_gene_pseudoT)

#############################
# 对两个谱系的差异基因操作
#############################
## 各自提取两个谱系中差异显著的基因
female_lineage1_clustering <- female_lineage1_sig_gene_pseudoT[female_lineage1_sig_gene_pseudoT$qval<0.05,]
female_lineage2_clustering <- female_lineage2_sig_gene_pseudoT[female_lineage2_sig_gene_pseudoT$qval<0.05,]

## 取两个谱系全部的基因，并进行去重复
# 需要理解unique()的用法: unique(c(1,2,1,2,3))
gene_list <- unique(rownames(female_lineage1_clustering), 
                    rownames(female_lineage2_clustering))
length(gene_list)

# 提取无重复HVGs的RPKM小表达矩阵
de_matrix <- log(females[rownames(females) %in% gene_list,]+1)
dim(females)
dim(de_matrix)

# 有了基因，就要分配到不同谱系的细胞

#############################
# 对两个谱系的细胞操作
#############################
## 第一个谱系
# 得到第一个谱系的细胞百分比
L1_lineage <- female_pseudotime[!is.na(female_pseudotime[,1]),1]
# 对百分比进行升序排序
L1_ordered_lineage <- L1_lineage[order(L1_lineage, 
                                       decreasing = FALSE)]
# 根据第一个谱系的排序后的细胞名称，得到属于它的表达矩阵
L1_cells <- de_matrix[,names(L1_ordered_lineage)]
(L1_cells[1:4,1:4])

## 同理得到L2
if(T){
  ## 第二个谱系
  # 得到第二个谱系的细胞百分比
  L2_lineage <- female_pseudotime[!is.na(female_pseudotime[,2]),2]
  # 对百分比进行升序排序（细胞数量从少到多）
  L2_ordered_lineage <- L2_lineage[order(L2_lineage, 
                                         decreasing = FALSE)]
  # 根据第二个谱系的细胞名称，得到属于它的表达矩阵
  L2_cells <- de_matrix[,names(L2_ordered_lineage)]
}

## 提取细胞名
L1_lineage_cells <- names(L1_ordered_lineage)
length(L1_lineage_cells)
L2_lineage_cells <- names(L2_ordered_lineage)
length(L2_lineage_cells)

# 看到L1_lineage有423个，L2_lineage有294个，而总共563个细胞。那么我们想知道，哪些细胞是两个谱系分化之前共有的，哪些是特有的
comp_list <- comparelists(L1_lineage_cells, L2_lineage_cells)
common_cells <- comp_list$intersect
L1_spe_cells <- L1_lineage_cells[!L1_lineage_cells %in% comp_list$intersect]
L2_spe_cells <- L2_lineage_cells[!L2_lineage_cells %in% comp_list$intersect]
# 共有154个，L1特有269个，L2特有140个
length(common_cells);length(L1_spe_cells);length(L2_spe_cells)

## 将细胞ID和谱系信息（包含common、L1、L2）对应起来
# 将L1特有和共有的标记成L1_cellLin
L1_cellLin <- c(
  rep_along("common cells", common_cells), 
  rep_along("L1 cells", L1_spe_cells)
)
names(L1_cellLin) <- c(common_cells, L1_spe_cells)
# 将L1_cellLin按照之前得到的L1表达矩阵列名重新排序
L1_cellLin <- L1_cellLin[match(colnames(L1_cells),names(L1_cellLin) )]

# 同理得到L2
if(T){
  L2_cellLin <- c(
    rep_along("common cells", common_cells), 
    rep_along("L2 cells", L2_spe_cells)
  )
  names(L2_cellLin) <- c(common_cells, L2_spe_cells)
  # 将L1_cellLin按照之前得到的L1表达矩阵列名重新排序
  L2_cellLin <- L2_cellLin[match(colnames(L2_cells),names(L2_cellLin) )]
}

## 将表达矩阵细胞名和cluster分群对应起来，并重新命名细胞ID
# Get the cell cluster of the cells for each cell lineage
cellType_L1 <- female_clustering[colnames(L1_cells)]
colnames(L1_cells) <- paste(colnames(L1_cells), "L1", sep="_")
names(L1_cellLin) <- colnames(L1_cells)

cellType_L2 <- female_clustering[colnames(L2_cells)]
colnames(L2_cells) <- paste(colnames(L2_cells), "L2", sep="_")
names(L2_cellLin) <- colnames(L2_cells)

## 细胞谱系信息合并
cellLin <- c(
  L2_cellLin,
  L1_cellLin
)

## 细胞cluster信息合并
cellType <- c(
  cellType_L2,
  cellType_L1
)

#############################
# 基因分类
#############################
## 表达量局部加权回归散点平滑法（locally weighted scatterplot smoothing，LOESS）
L2_cells_smooth <- smooth_gene_exp(
  L2_cells, 
  L2_ordered_lineage, 
  span=0.4
)
L1_cells_smooth <- smooth_gene_exp(
  L1_cells, 
  L1_ordered_lineage, 
  span=0.4
)
# 合并局部加权回归表达矩阵
data_heatmap <- data.frame(
  L2_cells_smooth,
  L1_cells_smooth
)

# 利用pheatmap函数进行层次聚类，只是为了调用它的算法而已，不是真的作图
set.seed(123)
gene_clustering <- pheatmap::pheatmap(
  data_heatmap, 
  scale="row", 
  clustering_method="ward.D",
  silent=TRUE
)

# 挑出17个基因cluster
clusters <- cutree(gene_clustering$tree_row, k = 17)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
head(clustering)
table(clustering[,1])


write.csv(clustering, 
          file="step5.2-gene_clustering_kmeans_k17_scaled.csv")

save(clusters,clustering,file = 'step5.2-gene_clustering.Rdata')
save(data_heatmap,cellLin,cellType,cellType_L2,
     file = 'step5.2-for_heatmap.Rdata')

