### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-23
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-Seurat-tSNE分群
### ---------------
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

#############################
# 使用原始count矩阵
#############################
# 下载：https://raw.githubusercontent.com/IStevant/XX-XY-mouse-gonad-scRNA-seq/master/data/female_count.Robj
load(file="../female_count.Robj")
load('../female_rpkm.Rdata')
# 直接对细胞和基因过滤
female_count <- female_count[rownames(female_count) %in% rownames(females),!colnames(female_count) %in% grep("rep",colnames(female_count), value=TRUE)]
female_count[1:3,1:3]
save(female_count,file = '../female_count.Rdata')

#############################
# 对细胞操作=》细胞发育时期的获取
#############################
load('../female_count.Rdata')
head(colnames(female_count))
female_stages <- sapply(strsplit(colnames(female_count), "_"), `[`, 1)
names(female_stages) <- colnames(female_count)
table(female_stages)

#############################
# 使用Seurat V3
#############################
library(Seurat)
packageVersion('Seurat')
load('../female_count.Rdata')
## 构建对象
sce_female <- CreateSeuratObject(counts = female_count, 
                           project = "sce_female", 
                           min.cells = 1, min.features = 0)
sce_female

## 添加样本注释信息
sce_female <- AddMetaData(object = sce_female, 
                      metadata = apply(female_count, 2, sum), 
                      col.name = 'nUMI_raw')
sce_female <- AddMetaData(object = sce_female, 
                      metadata = female_stages, 
                      col.name = 'female_stages')

## 数据归一化
sce_female <- NormalizeData(sce_female)
sce_female[["RNA"]]@data[1:3,1:3]

## 找差异基因HVGs
sce_female <- FindVariableFeatures(sce_female, 
                             selection.method = "vst", 
                             nfeatures = 2000)

# HVGs可视化
VariableFeaturePlot(sce_female)
# 保存这2000个HVGs
seurat3_HVGs <- VariableFeatures(sce_female)
save(seurat3_HVGs,file = 'step1.2-A-seurat3_HVGs.Rdata')
# 检查与之前得到的HVGs重合度
load('step1.1-A-females_hvg_matrix.Rdata')
load('step1.2-A-seurat3_HVGs.Rdata')
length(intersect(rownames(females_data),seurat3_HVGs))
# 结果和之前822个HVGs有434个重合

## 数据标准化
# 默认只对FindVariableFeatures得到的HVGs进行操作
sce_female <- ScaleData(object = sce_female, 
                    vars.to.regress = c('nUMI_raw'), 
                    model.use = 'linear', 
                    use.umi = FALSE)
## PCA降维
sce_female <- RunPCA(sce_female, 
                     features = VariableFeatures(object = sce_female))
## 降维后聚类
# 这里可以多选一些PCs
sce_female <- FindNeighbors(sce_female, dims = 1:20)
sce_female <- FindClusters(sce_female, resolution = 0.3)

## 继续tSNE非线性降维
ElbowPlot(sce_female)
sce_female_tsne <- RunTSNE(sce_female, dims = 1:9)
save(sce_female_tsne,file = 'step1.2-B-seurat3-female-tsne.Rdata')

DimPlot(object = sce_female_tsne, reduction = "tsne",
        group.by = 'female_stages')
ggsave('step1.2-C-tSNE_stages.pdf')

DimPlot(sce_female_tsne, reduction = "tsne")
ggsave('step1.2-D-tSNE_clusters.pdf')


## 加载DIY的tSNE结果比较一下
cluster1 <- read.csv('step1.1-D-female_clustering.csv')
load('step1.2-B-seurat3-female-tsne.Rdata')
cluster2 <- as.data.frame(Idents(sce_female_tsne))
# 把它们放在一起比较，前提条件是它们的行名相同
identical(cluster1[,1],rownames(cluster2))
table(cluster1[,2],cluster2[,1])






