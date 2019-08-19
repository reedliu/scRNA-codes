### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-07-18
### Email: jieandze1314@gmail.com
### Blog: www.jieandze1314.com
### CAAS/AGIS/SDAU 
### Update Log: 2019-07-18 Seurat V3.0 pbmc3k pre-process
### Update Log: 2019-08-19 scale、PCA
### From: https://satijalab.org/seurat/v3.0/pbmc3k_tutorial.html
### ---------------

install.packages('Seurat')
## 安装早期版本
if(F){
    install.packages('devtools')
    devtools::install_version(package = 'Seurat', version = package_version('2.3.0'))
}
## 降级R包版本
if(F){
    remove.packages('Seurat')
    pkgs = c( 'mixtools', 'lars', 'dtw', 'doSNOW', 'hdf5r' )
    #pkgs=c('jackstraw','slingshot')
    BiocManager::install(pkgs,ask = F,update = F)
    packageurl <- "https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_3.0.0.tar.gz"
    # 以后只需要修改这个版本号即可
    install.packages(packageurl, repos=NULL, type="source")
    library(Seurat)
}

###########################################
## 0-加载数据集
###########################################
rm(list = ls())
options(stringsAsFactors = F)

library(Seurat)
packageVersion("seurat") # 检查版本信息
library(dplyr)

# 设置Seurat对象
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")

# 这个矩阵长什么样？和我们平常见到的bulk转录组矩阵一样吗？
# 找几个基因，看看前30个细胞的情况。其中存在大量的.来代表0
pbmc.data[c("CD3D", "TCL1A", "MS4A1"), 1:30]

# 这里Seurat采用了一种聪明的再现方式（sparse），比原来用0表示的矩阵（dense）大大减小了空间占用，可以对比一下
if(F){
    dense.size <- object.size(as.matrix(pbmc.data))
    dense.size
    
    sparse.size <- object.size(pbmc.data)
    sparse.size
    
    dense.size/sparse.size #空间缩小23.8倍
}

# 用原始数据（非标准化）先构建一个对象
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

#访问这个对象(两个中括号是对assay的快速访问,pbmc[[]]就如同object@meta.data)
str(pbmc)
pbmc[["RNA"]]@counts
pbmc[[]]
pbmc@version

###########################################
## 1-预处理
###########################################
## 1.1质控筛选
# [[符号可以在pbmc这个对象中添加metadata，利用这个方法还可以挑出其他的feature
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
# 检查下metadata
head(pbmc@meta.data, 3)
# 绘制多个feature指标
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# FeatureScatter：绘制两组feature指标的相关性，最后组合在一起（CombinePlots）
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
# 最后进行一个过滤操作
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

###########################################
## 2-归一化
###########################################
## 2.1 数据归一化
pbmc <- NormalizeData(pbmc)
pbmc[["RNA"]]@data[1:3,1:3] #normalized
# 当然其中都是默认参数，它等于
if(F){
    tmp=pbmc
    tmp <- NormalizeData(tmp, normalization.method = "LogNormalize", scale.factor = 100000)
    tmp@assays$RNA[,]
}

## 2.2 找差异基因HVGs
# 有三种算法：vst、mean.var.plot、dispersion；默认选择2000个HVG；默认vst
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
# 可视化（with and without labels）
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

###########################################
## 3-标准化
###########################################
## 3.1 数据标准化
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)
if(F){
    #默认只对FindVariableFeatures得到的HVGs进行操作
    pbmc <- ScaleData(pbmc) 
}
pbmc[["RNA"]]@scale.data[1:3,1:3] #scaled
pbmc[["RNA"]]@counts[1:3,1:3] #raw

###########################################
## 4-线性降维--PCA
###########################################
## 4.1 线性降维--PCA
# 默认选择之前鉴定的差异基因（2000个）作为input；默认得到50个成分
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
# 结果保存在reductions 这个接口中

## 4.2 做个检查
print(pbmc[["pca"]], dims = 1:3, nfeatures = 5)
# 或者 print(pbmc@reductions$pca, dims = 1:3, nfeatures = 5)
VizDimLoadings(pbmc, dims = 1:3, reduction = "pca")
# 两个成分散点图
DimPlot(pbmc, reduction = "pca")
# 每个主成分的基因热图
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

## 4.3 选择合适数量的主成分来代表整个数据集
ElbowPlot(pbmc)
# 推荐开始不确定时要多选一些主成分

###########################################
## 5-聚类
###########################################
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)
head(Idents(pbmc), 5)
table(Idents(pbmc))
save(pbmc,file = "pbmc_5_cluster.RData")

###########################################
## 6-非线性降维tSNE、UMAP
###########################################
# umap
pbmc_umap <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc_umap, reduction = "umap")
# tsne
pbmc_tsne <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc_tsne, reduction = "tsne")
# 对计算过程很费时的结果保存一下
save(pbmc_umap, file = 'pbmc_umap.Rdata')

###########################################
## 7-找差异表达基因（cluster biomarker）
###########################################
pbmc = pbmc_tsne
# 找到cluster1中的marker基因
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
# 找到相对于cluster0、cluster3，在cluster5中的marker
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)

# 一步到位的办法
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# 这一步过滤好好理解(进行了分类、排序、挑前2个)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# 如果要使用其他的差异检验方法
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)

# 可视化
VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ", "PPBP", "CD8A"))

# DoHeatmap 对给定细胞和基因绘制热图，例如对每个cluster绘制前20个marker
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

###########################################
## 8-赋予每个cluster细胞类型
###########################################
# 重点在于找到marker对应的细胞类型（结合生物背景知识）
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()










