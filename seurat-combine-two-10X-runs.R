### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-07-19
### Email: jieandze1314@gmail.com
### CAAS/AGIS/SDAU 
### Update Log: 2019-07-19  Seurat3.0 - Combining Two 10X Runs
### Ref: https://satijalab.org/seurat/v3.0/merge_vignette.html
### ---------------

###########################
# 准备
###########################
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
if (!requireNamespace("Seurat"))
  BiocManager::install("Seurat")

rm(list = ls()) 
options(warn=-1) # 在全局变量中关闭warning
suppressMessages(library(Seurat))

###########################
# 原文背景
###########################
# 第一个病人的tumor表达矩阵,原文得到了2243个治疗之前（Pre）+ 5188个治疗之后（Ar）的细胞

###########################
# 加载cellranger filter的矩阵
###########################
tumor.pre.data <- Read10X(data.dir = "pre-hg38/")
tumor.pre <- CreateSeuratObject(counts = tumor.pre.data, project = "TumorPRE")
tumor.pre
# An object of class Seurat 
# 33694 features across 2099 samples within 1 assay 
# Active assay: RNA (33694 features)

tumor.ar.data <- Read10X(data.dir = "ar-hg38/")
tumor.ar <- CreateSeuratObject(counts = tumor.ar.data, project = "TumorAR")
tumor.ar
# An object of class Seurat 
# 33694 features across 4800 samples within 1 assay 
# Active assay: RNA (33694 features)

# 组合两个对象
tumor.combined <- merge(tumor.pre, y = tumor.ar, add.cell.ids = c("TumorPre", "TumorAr"), project = "p1tumor")
tumor.combined
# An object of class Seurat 
# 33694 features across 6899 samples within 1 assay 
# Active assay: RNA (33694 features)

###########################
# 走一波标准流程之得到对象
###########################
# 将10X的sparse松散数据变成标准的dense matrix
tumor_data <- as.matrix(GetAssayData(object = tumor.combined, slot = "counts"))
dim(tumor_data) # 33694 genes and 6899 cells from filtered cellranger count matrix

# 标准化
normTumor <- log2(1 + sweep(tumor_data, 2, median(colSums(tumor_data))/colSums(tumor_data), '*')) # Normalization

cellTypes <- sapply(colnames(normTumor), function(x) unlist(strsplit(x, "_"))[1]) 
cellTypes <-ifelse(cellTypes == 'TumorPre', 'Tumor_Before', 'Tumor_AcquiredResistance')
table(cellTypes) 
# 原文：Tumor_AcquiredResistance   Tumor_Before 
#        5188                       2243 


# 表达矩阵的质量控制
#简单看看每个细胞表达基因的数量，和每个基因在多少个细胞里面表达。
fivenum(apply(normTumor,1,function(x) sum(x>0) ))
# hg38-RP11-34P13.3 hg38-CH17-212P11.6        hg38-ODF3L2 
# 0                  0                  7 
# hg38-SLC39A13        hg38-MT-CO2 
# 275               6899
boxplot(apply(normTumor,1,function(x) sum(x>0) ))
fivenum(apply(normTumor,2,function(x) sum(x>0) ))
# TumorPre_AGAATAGCAAGTTAAG TumorPre_TGACGGCTCGGTGTTA 
# 404                      1205 
# TumorAr_CTCTAATTCCCTTGTG  TumorAr_AGTGAGGCAAGAGTCG 
# 1555                      2154 
# TumorAr_GATCGTAGTCATATGC 
# 5886 
hist(apply(normTumor,2,function(x) sum(x>0) ))

# 创建seurat对象
tumor <- CreateSeuratObject(normTumor, min.cells = 1, min.features = 0, project = '10x_Tumor')
tumor
# 原文21,861 genes and 7,431 cells，这里得到22,072 genes and 6,899 cells

# 添加metadata
# 3.0版本可以直接使用 object$name <- vector，当然也可以用AddMetaData
tumor <- AddMetaData(object = tumor, 
                     metadata = apply(tumor_data, 2, sum), 
                     col.name = 'nUMI_raw')
tumor <- AddMetaData(object = tumor, metadata = cellTypes, col.name = 'cellTypes')

###########################
# 走一波标准流程之质控
###########################
# 绘图的分组就会调用metadata的信息，例如添加的cellTypes就是免疫治疗前后的分组
features=c("nFeature_RNA", "nUMI_raw")
VlnPlot(object = tumor, 
        features = features, 
        group.by = 'cellTypes', ncol = 2)

# 3.0版本将GenePlot替换为FeatureScatter
# 版本2 GenePlot(object = sce, gene1 = "nUMI", gene2 = "nGene")
FeatureScatter(tumor,feature1 = "nUMI_raw",feature2 = "nFeature_RNA")

# 可以看看高表达量基因是哪些
# 3.0版本要将sce@raw.data替换成GetAssayData(object = , assay= ,slot = )
tail(sort(Matrix::rowSums(GetAssayData(tumor,assay = "RNA"))))
## 散点图可视化任意两个基因的一些属性（通常是细胞的度量）
# 这里选取两个基因。
tmp=names(sort(Matrix::rowSums(GetAssayData(tumor,assay = "RNA")),decreasing = T))
# 本文结果是
# hg38-RPS2 hg38-MT-CO1 hg38-MT-CO3 hg38-MT-CO2  hg38-RPS18 hg38-MALAT1
# 33834.55    34637.13    34879.63    35755.88    35997.47    46330.78 

# 原文结果是：
# RPS2   MT-CO1   MT-CO3   MT-CO2    RPS18   MALAT1 
# 34406.98 35104.22 35400.39 36274.50 36712.57 47996.87 

# 版本2 GenePlot(object = sce, gene1 = tmp[1], gene2 = tmp[2])
FeatureScatter(object = tumor, feature1 = tmp[1], feature2 = tmp[2])

# 散点图可视化任意两个细胞的一些属性（通常是基因的度量）
# 这里选取两个细胞

# 3.0版本将CellPlot替换成CellScatter，sce@cell.names换为colnames
# 版本2 CellPlot(sce,sce@cell.names[3],sce@cell.names[4],do.ident = FALSE)
CellScatter(tumor, colnames(tumor)[3],colnames(tumor)[4])

###########################
# 走一波标准流程之聚类可视化
###########################
# This process consists of data normalization and variable feature selection, data scaling, a PCA on variable features, construction of a shared-nearest-neighbors graph, and clustering using a modularity optimizer. Finally, we use a t-SNE to visualize our clusters in a two-dimensional space.

start_time <- Sys.time()
# Cluster tumor
tumor <- ScaleData(object = tumor, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)
# 3.0版本将FindVariableGenes换为FindVariableFeatures，另外将原来的cutoff进行整合，x轴统一归到mean.cutoff中，y轴归到dispersion.cutoff中
# 版本2 tumor <- FindVariableGenes(object = tumor, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
tumor <- FindVariableFeatures(object = tumor, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf))

tumor <- RunPCA(object = tumor, pc.genes = VariableFeatures(tumor))

tumor <- RunTSNE(object = tumor, dims.use = 1:10)

DimPlot(tumor, group.by = 'cellTypes', colors.use = c('#EF8A62', '#67A9CF'))

end_time <- Sys.time()
end_time - start_time








