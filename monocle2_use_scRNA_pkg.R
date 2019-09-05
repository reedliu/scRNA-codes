### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-09-02
### Email: jieandze1314@gmail.com
### Blog: www.jieandze1314.com
### Update Log: 2019-09-02 用scRNA包学monocle2
### ---------------

rm(list = ls()) 
Sys.setenv(R_MAX_NUM_DLLS=999)
options(stringsAsFactors = F)
options(warn=-1) 

##########################
# 关于monocle包降级
##########################
if(F){
  packageVersion("monocle")
  remove.packages('monocle')
  pkgs = c( 'mixtools', 'lars', 'dtw', 'doSNOW', 'hdf5r' ) 
  #pkgs=c('jackstraw','slingshot')
  BiocManager::install(pkgs,ask = F,update = F)
  BiocManager::install("monocle")
  library(monocle)
  packageVersion("monocle")
}
library(monocle)

##########################
# 加载scRNA包
##########################
library(scRNAseq)
# 加载测试数据 
data(fluidigm)

##########################
# monocle的CDS对象需要3个要素
##########################
# 第一个：RSEM表达矩阵（ct = count）
assay(fluidigm)  <-  assays(fluidigm)$rsem_counts
ct <- floor(assays(fluidigm)$rsem_counts)
ct[1:4,1:4] 
# 第二个：临床信息
sample_ann <- as.data.frame(colData(fluidigm))
# 第三个：基因注释信息（必须包含一列是gene_short_name）
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
# 然后转换为AnnotatedDataFrame对象
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)

# 最后构建CDS对象
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

##########################
# 质控过滤
##########################
cds=sc_cds
cds ## 原始数据有： 26255 features, 130 samples 

# 设置一个基因表达量的过滤阈值，结果会在cds@featureData@data中新增一列num_cells_expressed，记录这个基因在多少细胞中有表达
cds <- detectGenes(cds, min_expr = 0.1)
print(head(cds@featureData@data))
# 在monocle版本2.12.0中，取消了fData函数（此前在2.10版本中还存在），不过在monocle3中又加了回来
# 如果遇到不能使用fData的情况，就可以采用备选方案：cds@featureData@data

## 然后进行基因过滤
expressed_genes <- row.names(subset(cds@featureData@data,
                                    num_cells_expressed >= 5))
length(expressed_genes)
cds <- cds[expressed_genes,]
cds
# 过滤基因后剩下：13385 features, 130 samples 

## 还可以进行细胞层面过滤
# 依然是：如果不支持使用pData()函数，可以使用cds@phenoData@data来获得各种细胞注释信息
print(head(cds@phenoData@data))
# 比如我们看一下细胞注释的第一个NREADS信息
tmp=pData(cds)
fivenum(tmp[,1])
## [1]    91616   232899   892209  8130850 14477100

# 如果要过滤细胞，其实也是利用subset函数，不过这里不会对细胞过滤
valid_cells <- row.names(cds@phenoData@data)
cds <- cds[,valid_cells]
cds 

##########################
# 聚类
##########################
## 不使用marker基因聚类
# step1：判断使用哪些基因进行细胞分群判断
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds) # 挑有差异的
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) # 挑表达量不太低的
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id) # 准备聚类基因名单
plot_ordering_genes(cds) 
# step2:然后选一下主成分
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'
# step3:然后选一下主成分，这里选前6个成分（大概在第一个拐点处）
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
cds <- clusterCells(cds, num_clusters = 4) 
# step4：可视化
plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")
table(cds@phenoData@data$Biological_Condition)

# 如果要去除一些干扰因素，例如批次效应
if(F){
  cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                         reduction_method = 'tSNE',
                         residualModelFormulaStr = "~Biological_Condition + num_genes_expressed",
                         verbose = T)
  cds <- clusterCells(cds, num_clusters = 4)
  plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")
  # 可以看到，去掉本来的生物学意义后，最后细胞是会被打散的。
  # 所以residualModelFormulaStr这个东西的目的就是磨平它参数包含的差异
  
  # 如果去除生物意义以外的效应
  cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                         reduction_method = 'tSNE',
                         residualModelFormulaStr = "~NREADS + num_genes_expressed",
                         verbose = T)
  cds <- clusterCells(cds, num_clusters = 4)
  plot_cell_clusters(cds, 1, 2, color = "Biological_Condition")
}

##########################
# 差异分析
##########################
start=Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~Biological_Condition")
end=Sys.time()
end-start
# 然后得到差异基因
sig_genes <- subset(diff_test_res, qval < 0.1)
head(sig_genes[,c("gene_short_name", "pval", "qval")] )
# 作图（注意要将基因名变成character）
cg=as.character(head(sig_genes$gene_short_name))
plot_genes_jitter(cds[cg,], 
                  grouping = "Biological_Condition", ncol= 2)
# 加个分组颜色
plot_genes_jitter(cds[cg,],
                  grouping = "Biological_Condition",
                  color_by = "Biological_Condition",
                  nrow= 3,
                  ncol = NULL )

# 常规的箱线图
boxplot(log10(cds@assayData$exprs["A1BG",]+1) ~ cds@phenoData@data$Biological_Condition)

##########################
# 推断发育轨迹
##########################
# step1：选合适基因
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
cds <- setOrderingFilter(cds, ordering_genes)
plot_ordering_genes(cds)
# step2: 降维
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')
# step3: 细胞排序
cds <- orderCells(cds)
# 最后可视化
plot_cell_trajectory(cds, color_by = "Biological_Condition")  
# 基因在不同细胞中的表达量变化进行绘图
plot_genes_in_pseudotime(cds[cg,],
                         color_by = "Biological_Condition")


