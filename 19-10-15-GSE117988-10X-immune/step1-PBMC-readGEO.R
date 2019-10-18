### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-15
### Email: jieandze1314@gmail.com
### Title: 10X scRNA免疫治疗-PBMC-读取GEO并走Seurat流程
### ---------------
rm(list = ls()) 
options(warn=-1) 
suppressMessages(library(Seurat))

#############################
# 读取表达矩阵
#############################
start_time <- Sys.time()
raw_dataPBMC <- read.csv('./GSE117988_raw.expMatrix_PBMC.csv.gz', header = TRUE, row.names = 1)
end_time <- Sys.time()
end_time - start_time

#############################
# 复现作者的分群结果
#############################
# 参考https://www.researchgate.net/publication/328016998_Supplementary_Material_6/data/5bb2eac9299bf13e605a0a74/41467-2018-6300-MOESM6-ESM.txt
##-- PBMC
# Load the data & Normalization according to Zheng et al., 2017
dim(raw_dataPBMC) 
## step1: 归一化
dataPBMC <- log2(1 + sweep(raw_dataPBMC, 2, 
                           median(colSums(raw_dataPBMC))/colSums(raw_dataPBMC), '*')) # Normalization

head(colnames(dataPBMC))

## step2: 自定义划分时间点
# 作者利用的是Seurat V2的ExtractField函数
# for Seurat V2
timePoints <- sapply(colnames(dataPBMC), 
                     function(x) ExtractField(x, 2, '[.]'))
# 但如果使用V3，就要用常规方法strsplit，并且注意点号是正则匹配符，需要用\\来转义
timePoints <- sapply(colnames(dataPBMC), function(x) unlist(strsplit(x, "\\."))[2]) 
table(timePoints)
timePoints <-ifelse(timePoints == '1', 'PBMC_Pre', 
                    ifelse(timePoints == '2', 'PBMC_EarlyD27',
                           ifelse(timePoints == '3', 'PBMC_RespD376', 'PBMC_ARD614')))

## step3: 表达矩阵质控
# 第一点：基因在多少细胞表达 
fivenum(apply(dataPBMC,1,function(x) sum(x>0) ))
# 第二点：细胞中有多少表达的基因
fivenum(apply(dataPBMC,2,function(x) sum(x>0) ))

## step4: 创建Seurat对象
PBMC <- CreateSeuratObject(dataPBMC, 
                           min.cells = 1, min.features = 0, project = '10x_PBMC')
PBMC # 17,712 genes and 12,874 cells

## step5: 添加metadata (nUMI and timePoints)
PBMC <- AddMetaData(object = PBMC, metadata = apply(raw_dataPBMC, 2, sum), col.name = 'nUMI_raw')
PBMC <- AddMetaData(object = PBMC, metadata = timePoints, col.name = 'TimePoints')

## step6:  聚类标准流程
# Seurat V2
PBMC <- ScaleData(object = PBMC, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)
PBMC <- FindVariableGenes(object = PBMC, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes)
PBMC <- FindClusters(object = PBMC, reduction.type = "pca", dims.use = 1:10, resolution = 1, k.param = 35, save.SNN = TRUE) # 13 clusters
PBMC <- RunTSNE(object = PBMC, dims.use = 1:10)
TSNEPlot(PBMC, colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'))


# Seurat V3
# detach("package:Seurat", unload=TRUE)
# .libPaths()
# library(Seurat, lib.loc="/var/R-3.6.0/library")
# packageVersion('Seurat')
PBMC_V3 <- ScaleData(object = PBMC, vars.to.regress = c('nUMI_raw'), model.use = 'linear', use.umi = FALSE)

PBMC_V3 <- FindVariableFeatures(object = PBMC_V3, mean.function = ExpMean, dispersion.function = LogVMR, mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf))

PBMC_V3 <- RunPCA(object = PBMC_V3, pc.genes = VariableFeatures(PBMC_V3))

PBMC_V3 <- FindNeighbors(PBMC_V3, reduction = "pca", dims = 1:10,
                      k.param = 35)
PBMC_V3 <- FindClusters(object = PBMC_V3, 
                     resolution = 0.9, verbose=F) 

PBMC_V3 <- RunTSNE(object = PBMC_V3, dims.use = 1:10,
                reduction.name='RNA_snn_res.0.9')

DimPlot(PBMC_V3, cols = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'))
# TSNEPlot(PBMC_V3, colors = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'))

# save(PBMC_V3,file = 'patient1.PBMC.V3.output.Rdata')
save(PBMC,PBMC_V3,file = 'patient1.PBMC.V2V3.output.Rdata')





