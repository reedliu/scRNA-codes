### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-差异分析-monocle
### ---------------

#############################
# 准备表达矩阵和分群信息
#############################
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

#差异分析一般都要求使用count矩阵
load('../female_count.Rdata')
# 6个发育时期获取
head(colnames(female_count))
female_stages <- sapply(strsplit(colnames(female_count), "_"), `[`, 1)
names(female_stages) <- colnames(female_count)
table(female_stages)
# 4个cluster获取
cluster <- read.csv('../step1-female-RPKM-tSNE/step1.1-D-female_clustering.csv')
female_clustering=cluster[,2];names(female_clustering)=cluster[,1]
table(female_clustering)

#############################
# 利用monocle V2构建对象 
# 需要三样东西：表达矩阵、细胞信息、基因信息
#############################
## 表达矩阵
dim(female_count)
## 细胞分群信息（包括6个stage和4个cluster）
table(female_stages)
table(female_clustering)
## 基因信息
gene_annotation <- as.data.frame(rownames(count_matrix))

## 开始monocle
library(monocle) 
packageVersion('monocle')
# 直接使用作者包装的函数，代码更简洁
DE_female <- prepare_for_DE (
  female_count, 
  female_clustering, 
  female_stages
)

# 想了解详细的步骤，如下
if(F){
  # 1.表达矩阵
  expr_matrix <- as.matrix(count_matrix)
  # 2.细胞信息
  sample_ann <- data.frame(cells=names(count_matrix), 
                           stages=stage, 
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
  
  sc_cds <- detectGenes(sc_cds, min_expr = 5)
  sc_cds <- sc_cds[fData(sc_cds)$num_cells_expressed > 10, ]
  sc_cds <- estimateSizeFactors(sc_cds)
  sc_cds <- estimateDispersions(sc_cds)
}

#############################
# 利用monocle V2进行差异分析
#############################
# 其实也是包装了differentialGeneTest函数
start_time <- Sys.time()
female_DE_genes <- findDEgenes(
  DE_female, 
  qvalue=0.05
)
end_time <- Sys.time()
end_time - start_time
# 得到了"4435 significantly DE genes (FDR<0.05)."
save(female_DE_genes,file = 'step3.1-A-DEG-dataframe.Rdata')

#############################
# 利用monocle V2进行差异分析
#############################
load('../female_rpkm.Rdata')
# 上面得到的4435个差异基因，我们要知道这些基因在哪个cluster
# 作者先得到了每个差异基因在不同cluster的平均表达量，然后找平均表达量最大的那个cluster，就认为这个基因属于这个cluster
de_clusters <- get_up_reg_clusters(
  females, 
  female_clustering, 
  female_DE_genes
)

# write.csv(
#   de_clusters, 
#   quote = FALSE, 
#   file= "step3.1-B-DEG-in-cluster.csv"
# )

save(de_clusters,file='step3.1-DEG-monocle_summary.Rdata')


