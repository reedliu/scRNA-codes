### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-18
### Email: jieandze1314@gmail.com
### Title: 10X scRNA免疫治疗-Tumor-走Seurat流程并探索marker基因表达
### ---------------
rm(list = ls()) 
options(warn=-1) 
suppressMessages(library(Seurat))
packageVersion('Seurat')
#############################
# 读取表达矩阵
#############################
# 下载地址：ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE117nnn/GSE117988/suppl/GSE117988_raw.expMatrix_Tumor.csv.gz
start_time <- Sys.time()
raw_dataTumor <- read.csv('./GSE117988_raw.expMatrix_Tumor.csv.gz', header = TRUE, row.names = 1)
end_time <- Sys.time()
end_time - start_time

#############################
# 走标准SeuratV2流程
#############################
# 注意：之前PBMC包含了四个时间点的样本,这里的Tumor包含了2个样本（治疗之前Pre和复发AR）
# 参考https://www.researchgate.net/publication/328016998_Supplementary_Material_6/data/5bb2eac9299bf13e605a0a74/41467-2018-6300-MOESM6-ESM.txt
##-- Tumor
# Load the data & Normalization according to Zheng et al., 2017
dim(raw_dataTumor) # 7,431 cells and 21,861 genes - already filtered

## step1: 归一化
dataTumor <- log2(1 + sweep(raw_dataTumor, 2, 
                            median(colSums(raw_dataTumor))/colSums(raw_dataTumor), '*')) # Normalization
head(colnames(dataTumor))

## step2: 自定义划分细胞类型
cellTypes <- sapply(colnames(dataTumor), function(x) ExtractField(x, 2, '[.]'))
cellTypes <-ifelse(cellTypes == '1', 'Tumor_Before', 'Tumor_AcquiredResistance')
table(cellTypes)

## step3: 表达矩阵质控
# 第一点：基因在多少细胞表达 
fivenum(apply(dataTumor,1,function(x) sum(x>0) ))
# 第二点：细胞中有多少表达的基因
fivenum(apply(dataTumor,2,function(x) sum(x>0) ))

## step4: 创建Seurat对象
tumor <- CreateSeuratObject(dataTumor, 
                           min.cells = 1, min.features = 0, project = '10x_Tumor')
tumor 

## step5: 添加metadata (nUMI 和 细胞类型)
tumor <- AddMetaData(object = tumor, metadata = apply(raw_dataTumor, 2, sum), col.name = 'nUMI_raw')
tumor <- AddMetaData(object = tumor, metadata = cellTypes, col.name = 'cellTypes')

## step6:  聚类标准流程
start_time <- Sys.time()
tumor <- ScaleData(object = tumor, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)
# 如果要使用多线程
# tumor <- ScaleData(object = tumor, vars.to.regress = c('nUMI_raw'), model.use = 'linear',
#                    use.umi = FALSE,do.par=T,num.cores =16)
tumor <- FindVariableGenes(object = tumor, mean.function = ExpMean, 
                           dispersion.function = LogVMR, x.low.cutoff = 0.0125, 
                           x.high.cutoff = 3, y.cutoff = 0.5)
tumor <- RunPCA(object = tumor, pc.genes = tumor@var.genes)
tumor <- RunTSNE(object = tumor, dims.use = 1:10, perplexity = 25)
end_time <- Sys.time()
end_time - start_time

TSNEPlot(tumor, group.by = 'cellTypes', colors.use = c('#EF8A62', '#67A9CF'))
save(tumor,file = 'patient1.Tumor.V2.output.Rdata')

#############################
# 基因可视化
#############################
rm(list = ls()) 
options(warn=-1) 
start_time <- Sys.time()
load('patient1.Tumor.V2.output.Rdata')
end_time <- Sys.time()
end_time - start_time
TSNEPlot(tumor, group.by = 'cellTypes', colors.use = c('#EF8A62', '#67A9CF'))
# 取出log归一化后的表达矩阵
count_matrix=tumor@data
count_matrix[1:4,1:4]
# 取出细胞分群信息
cluster=tumor@meta.data$cellTypes
table(cluster)
# 提取基因信息
allGenes = row.names(tumor@raw.data)
allGenes[grep('HLA',allGenes)]
# 对HLA-A操作
FeaturePlot(object = tumor, 
            features.plot ='HLA-A', 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne") 
# 看HLA-A表达量
table(count_matrix['HLA-A',]>0, cluster)
tmp <- chisq.test(table(count_matrix['HLA-A',]>0, cluster))
tmp$p.value

# 对HLA-B操作
FeaturePlot(object = tumor, 
            features.plot ='HLA-B', 
            cols.use = c("grey", "blue"), 
            reduction.use = "tsne") 

table(count_matrix['HLA-B',]>0, cluster)


# 利用循环找相似表达模式的基因（在一个群表达，另外一个群不表达，找p值极显著的）
HLA_genes <- allGenes[grep('HLA',allGenes)]

HLA_result <- c()
for (gene in HLA_genes) {
  tmp <- chisq.test(table(count_matrix[gene,]>0, cluster))
  if (tmp$p.value<0.01) {
    HLA_result[gene] <- gene
  }
}
names(HLA_result)

