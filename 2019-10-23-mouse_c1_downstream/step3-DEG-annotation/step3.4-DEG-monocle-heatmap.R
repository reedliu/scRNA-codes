### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-差异分析-monocle结果画热图
### ---------------

#############################
# 准备表达矩阵和分群信息
#############################
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

# 加载RPKM表达矩阵
load('../female_rpkm.Rdata')
# 加载monocle差异分析的结果
load('step3.1-DEG-monocle_summary.Rdata')
# 6个发育时期获取
head(colnames(females))
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(females) <- colnames(females)
table(female_stages)
# 4个cluster获取
cluster <- read.csv('../step1-female-RPKM-tSNE/step1.1-D-female_clustering.csv')
female_clustering=cluster[,2];names(female_clustering)=cluster[,1]
table(female_clustering)

#############################
# 获得差异基因
#############################
e_genes <- de_clusters
gene_names <- subset(de_genes, qval<0.0001)
gene_names <- gene_names$genes

gene_names <- get_top_up_reg_clusters(de_clusters, 20)
gene_subset <- as.matrix(log(females[rownames(females) %in% gene_names,]+1))

cl1_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C1"])]
cl2_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C2"])]
cl3_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C3"])]
cl4_gene_subset <- gene_subset[, colnames(gene_subset) %in% names(female_clustering[female_clustering=="C4"])]

heatmap_gene_subset <- cbind(
  cl1_gene_subset, 
  cl2_gene_subset,
  cl3_gene_subset,
  cl4_gene_subset
)

markerGenes <- c(
  "Nr2f1",
  "Nr2f2",
  "Maf",
  "Foxl2",
  "Rspo1",
  "Lgr5",
  "Bmp2",
  "Runx1",
  "Amhr2",
  "Kitl",
  "Fst",
  "Esr2",
  "Amh",
  "Ptges"
)

heatmap_gene_subset <- heatmap_gene_subset[order(match(rownames(heatmap_gene_subset), markerGenes)),]
heatmap_female_stages <- sapply(strsplit(colnames(heatmap_gene_subset), "_"), `[`, 1)
table(heatmap_female_stages)
female_stages=heatmap_female_stages
table(female_clustering)
rowbreaks <- c(6, 15)

colbreaks <- c(
  ncol(cl1_gene_subset),
  ncol(cl1_gene_subset)+ncol(cl2_gene_subset), 
  ncol(cl1_gene_subset)+ncol(cl2_gene_subset)+ncol(cl3_gene_subset)
)

cluster_color <- c(
  C1="#560047",
  C2="#a53bad", 
  C3="#eb6bac", 
  C4="#ffa8a0"
)

stage_color=c(
  E10.5="#2754b5", 
  E11.5="#8a00b0", 
  E12.5="#d20e0f", 
  E13.5="#f77f05", 
  E16.5="#f9db21",
  P6="#43f14b"
)

library(pheatmap)
png("step3.4-A-monocle_DEG_heatmap.png")
plot_heatmap_2(
  heatmap_gene_subset, 
  female_clustering, 
  female_stages, 
  rowbreaks, 
  colbreaks,
  cluster_color,
  stage_color
)
dev.off()








