### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-monocle找谱系中变化显著的基因
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


#############################
# Monocle找不同谱系之间高变化基因
#############################
# 在第一个谱系中变化更剧烈的
female_lineage1_sig_gene_pseudoT <- get_var_genes_pseudotime(
  females, 
  female_count, 
  female_pseudotime, 
  lineageNb=1, 
  female_clustering
)

dim(female_lineage1_sig_gene_pseudoT)
# 从中找到差异显著的基因，根据qval<0.05过滤
# 从12612个基因里面挑选出2861个差异显著的基因
female_lineage1_sig_gene_pseudoT <- female_lineage1_sig_gene_pseudoT[female_lineage1_sig_gene_pseudoT$qval<0.05,]
dim(female_lineage1_sig_gene_pseudoT)

write.csv(female_lineage1_sig_gene_pseudoT, 
          file= "step5.1-lineage1_pseudotime_DE_genes.csv")


## 在第二个谱系中变化更剧烈的
female_lineage2_sig_gene_pseudoT <- get_var_genes_pseudotime(
  females, 
  female_count, 
  female_pseudotime, 
  lineageNb=2, 
  female_clustering
)

dim(female_lineage2_sig_gene_pseudoT)
# 从11937个基因里面挑选出2182个差异显著的基因
female_lineage2_sig_gene_pseudoT <- female_lineage2_sig_gene_pseudoT[female_lineage2_sig_gene_pseudoT$qval<0.05,]
dim(female_lineage2_sig_gene_pseudoT)

write.csv(female_lineage2_sig_gene_pseudoT, 
          file="step5.1-lineage2_pseudotime_DE_genes.csv")

save(female_lineage1_sig_gene_pseudoT,
     female_lineage2_sig_gene_pseudoT,
     file = 'step5.1-lineage_sig_gene.Rdata')







