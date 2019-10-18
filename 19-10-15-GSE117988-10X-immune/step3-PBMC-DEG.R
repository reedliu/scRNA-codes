### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-17
### Email: jieandze1314@gmail.com
### Title: 10X scRNA免疫治疗-PBMC-差异分析
### ---------------
rm(list = ls()) 
options(warn=-1)
load('patient1.PBMC.V2.output.Rdata')

#############################
# 准备数据=> SubsetData()取子集
#############################
PBMC_RespD376 = SubsetData(PBMC,TimePoints =='PBMC_RespD376')
table(PBMC_RespD376@ident)
# 提取出来第4、10群
PBMC_RespD376_for_DEG = SubsetData(PBMC_RespD376,
                                   PBMC_RespD376@ident %in% c(4,10))


#############################
# 利用monocle V2构建对象 
#############################
# 需要三样东西：表达矩阵、细胞信息、基因信息
count_matrix=PBMC_RespD376_for_DEG@data
dim(count_matrix)
# 细胞分群信息
cluster=PBMC_RespD376_for_DEG@ident
table(cluster)
# 基因信息
gene_annotation <- as.data.frame(rownames(count_matrix))
# 开始monocle
library(monocle) 
packageVersion('monocle')
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
  lowerDetectionLimit=1)


#############################
# monocle V2质控过滤
#############################
cds=sc_cds
cds <- detectGenes(cds, min_expr = 0.1)
# 结果保存在cds@featureData@data
print(head(cds@featureData@data))
# 基因过滤
expressed_genes <- row.names(subset(cds@featureData@data,
                                    num_cells_expressed >= 5))
length(expressed_genes)
cds <- cds[expressed_genes,]

#############################
# monocle V2 聚类
#############################
# step1：dispersionTable() 目的是判断使用哪些基因进行细胞分群
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
disp_table <- dispersionTable(cds) # 挑有差异的
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) # 挑表达量不太低的
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)  # 准备聚类基因名单
plot_ordering_genes(cds) 
# 图中黑色的点就是被标记出来一会要进行聚类的基因

# [可省略]step2：plot_pc_variance_explained()  选一下主成分
plot_pc_variance_explained(cds, return_all = F) # norm_method='log'

# step3: 聚类
# 进行降维
cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                       reduction_method = 'tSNE', verbose = T)
# 进行聚类
cds <- clusterCells(cds, num_clusters = 4) 
# Distance cutoff calculated to 1.812595  
plot_cell_clusters(cds, 1, 2, color = "cellType")


#############################
# monocle V2 差异分析
#############################
# 这个过程比较慢！
start=Sys.time()
diff_test_res <- differentialGeneTest(cds,
                                      fullModelFormulaStr = "~cellType")
end=Sys.time()
end-start
# 得到差异基因
sig_genes <- subset(diff_test_res, qval < 0.1)
nrow(sig_genes)
head(sig_genes[,c("genes", "pval", "qval")] )

#############################
# 热图可视化
#############################
htmapGenes=c(
  'GAPDH','CD52','TRAC','IL32','ACTB','ACTG1','COTL1',
  'GZMA','GZMB','GZMH','GNLY'
)
htmapGenes %in% rownames(sig_genes)
# 最原始
library(pheatmap)
dat=count_matrix[htmapGenes,]
pheatmap(dat)
# 修改
n=t(scale(t(dat)))
n[n>2]=2 
n[n< -1]= -1
ac=data.frame(group=cluster)
rownames(ac)=colnames(n)

pheatmap(n,annotation_col = ac,
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F, 
         cluster_rows = F)





