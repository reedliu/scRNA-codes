### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-18
### Email: jieandze1314@gmail.com
### Title: 10X scRNA免疫治疗-以分群为例，一定要使用R包吗？
### ---------------
rm(list = ls()) 
options(warn=-1) 
load('patient1.PBMC_RespD376_for_DEG.Rdata')

count_matrix[1:4,1:4]
table(cluster)

#############################
# 使用常规函数，不使用包
#############################
# 根据分群因子，赋予不同颜色
color <- cluster
levels(color) <- rainbow(2)
table(color)
if(T){
  choosed_count <- count_matrix
  # 表达矩阵过滤
  choosed_count <- choosed_count[apply(choosed_count, 1, sd)>0,]
  choosed_count <- choosed_count[names(head(sort(apply(choosed_count, 1, sd),decreasing = T),1000)),]
  # 然后进行PCA分析
  pca_out <- prcomp(t(choosed_count),scale. = T)
  library(ggfortify)
  autoplot(pca_out, col=color) +theme_classic()+ggtitle('PCA plot')
  
  str(pca_out)
  pca_out$x[1:3,1:3]
  library(Rtsne)
  tsne_out <- Rtsne(pca_out$x[,1:6], perplexity = 10,
                    pca = F, max_iter = 2000,
                    verbose = T)
  tsnes_cord <- tsne_out$Y
  colnames(tsnes_cord) <- c('tSNE1','tSNE2')
  ggplot(tsnes_cord, aes(x=tSNE1, y = tSNE2)) + geom_point(col=color) + theme_classic()+ggtitle('tSNE plot')
}

#############################
# 使用Seurat V2
#############################
PBMC <- CreateSeuratObject(count_matrix, 
                           min.cells = 1, min.features = 0, 
                           project = '10x_PBMC')
PBMC
PBMC <- AddMetaData(object = PBMC, metadata = apply(count_matrix, 2, sum), col.name = 'nUMI_raw')
PBMC <- AddMetaData(object = PBMC, metadata = cluster, col.name = 'cluster')
VlnPlot(PBMC, features.plot  = c("nUMI_raw", "nGene"), 
        group.by = 'cluster',nCol = 2) #怀疑之前的PCA结果是由于基因和UMI数量导致的
# 下面进行校正
PBMC <- ScaleData(PBMC, vars.to.regress =  c("nUMI_raw", "nGene"),
                  model.use = 'linear',
                  use.umi = F)
# 这个FindVariableGenes其实和前面的找top1000 sd基因目的一样，就是为了增加分析的有效性，认为高变化的基因存储的信息更值得关注
PBMC <- FindVariableGenes(object = PBMC, mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                          y.cutoff = 0.5)
length(PBMC@var.genes)
PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes)
PBMC <- RunTSNE(object = PBMC, dims.use = 1:20)
TSNEPlot(PBMC, group.by='cluster')

#############################
# 使用Seurat V3
#############################
PBMC_V3 <- CreateSeuratObject(counts = count_matrix,
                                 min.cells = 1, 
                                 min.features = 0, 
                                 project = "10x_PBMC_V3")
PBMC_V3 <- AddMetaData(object = PBMC_V3, metadata = apply(count_matrix, 2, sum), col.name = 'nUMI_raw')
PBMC_V3 <- AddMetaData(object = PBMC_V3, metadata = cluster, col.name = 'cluster')
VlnPlot(PBMC_V3, features =  c("nUMI_raw", "nFeature_RNA"), ncol = 2,group.by = 'cluster')
PBMC_V3 <- ScaleData(object = PBMC_V3, vars.to.regress =  c("nUMI_raw", "nFeature_RNA"), 
                     model.use = 'linear', use.umi = FALSE)

PBMC_V3 <- FindVariableFeatures(object = PBMC_V3, mean.function = ExpMean, 
                                dispersion.function = LogVMR, 
                                mean.cutoff = c(0.0125,3), 
                                dispersion.cutoff = c(0.5,Inf))
length(VariableFeatures(object = PBMC_V3))
PBMC_V3 <- RunPCA(PBMC_V3, features = VariableFeatures(object = PBMC_V3))
pbmc_tsne <- RunTSNE(PBMC_V3, dims = 1:20)
DimPlot(pbmc_tsne, reduction = "tsne",group.by='cluster')


#############################
# 使用monocle
#############################
library(monocle) 
# 1.表达矩阵
expr_matrix <- as.matrix(count_matrix)
# 2.细胞信息
sample_ann <- data.frame(cells=names(count_matrix),  
                         cellType=cluster)
rownames(sample_ann)<- names(count_matrix)
# 3.基因信息
gene_ann <- as.data.frame(rownames(count_matrix))
rownames(gene_ann)<- rownames(count_matrix)
colnames(gene_ann)<- "genes"
# 然后转换为AnnotatedDataFrame对象
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
# 最后构建CDS对象
sc_cds <- newCellDataSet(
  expr_matrix, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=0.5)

cds=sc_cds
cds <- detectGenes(cds, min_expr = 1)
expressed_genes <- row.names(subset(cds@featureData@data,
                                    num_cells_expressed >= 1))
length(expressed_genes)
cds <- cds[expressed_genes,]

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds) # 挑有差异的
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) # 挑表达量不太低的
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)  # 准备聚类基因名单
# 进行降维
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
# 进行聚类
cds <- clusterCells(cds, num_clusters = 4) 
# Distance cutoff calculated to 1.812595  
plot_cell_clusters(cds, 1, 2, color = "cellType")

#############################
# 使用scater
#############################
suppressMessages(library(scater))
## 创建 scater 要求的对象
# 在导入对象之前，最好是将表达量数据存为矩阵
sce <- SingleCellExperiment(
  assays = list(counts = as.matrix(count_matrix)), 
  colData = data.frame(cluster)
)
# 预处理
exprs(sce) <- log2(calculateCPM(sce ) + 1)
# 降维
sce <- runPCA(sce)
plotReducedDim(sce, use_dimred = "PCA", 
               colour_by= "cluster")
set.seed(1000)
sce <- runTSNE(sce, perplexity=10)
plotTSNE(sce, 
         colour_by= "cluster")
































