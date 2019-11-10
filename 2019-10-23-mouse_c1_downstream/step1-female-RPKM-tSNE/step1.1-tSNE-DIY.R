### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-23
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-DIY-tSNE分群
### ---------------
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)

source("../analysis_functions.R")
#############################
# 对细胞操作=》细胞发育时期的获取
#############################
load('../female_rpkm.Rdata')
dim(females)

head(colnames(females))
## 取下划线分隔的第一部分
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
# 或者
female_stages <- sapply(strsplit(colnames(females), "_"), 
                        function(x)x[1])
# 再或者
female_stages <- stringr::str_split(colnames(females),'_', simplify = T)[,1]

names(female_stages) <- colnames(females)
table(female_stages)

#############################
# 对基因操作=》基因过滤与统计
#############################
## 去掉在所有细胞都不表达的基因
(dim(females))
females <- females[rowSums(females)>0,]
(dim(females))

## 探索mean,sd,mad,cv
mean_per_gene <- apply(females, 1, mean, na.rm = TRUE) 
sd_per_gene <- apply(females, 1, sd, na.rm = TRUE) 
mad_per_gene <- apply(females, 1, mad, na.rm = TRUE) 
cv = sd_per_gene/mean_per_gene
library(matrixStats)
var_per_gene <- rowVars(as.matrix(females))
cv2=var_per_gene/mean_per_gene^2
# 存储统计结果
cv_per_gene <- data.frame(mean = mean_per_gene,
                          sd = sd_per_gene,
                          mad=mad_per_gene,
                          var=var_per_gene,
                          cv=cv,
                          cv2=cv2)
rownames(cv_per_gene) <- rownames(females)
head(cv_per_gene)
# 根据表达量过滤统计结果
cv_per_gene=cv_per_gene[cv_per_gene$mean>1,]
# 简易的可视化
with(cv_per_gene,plot(log10(mean),log10(cv2)))
# 更加复杂的（其实就是求每列之间的相关性）
library(psych)
pairs.panels(cv_per_gene, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
)
## 使用作者的DIY函数
females_data <- getMostVarGenes(females, fitThr=2)
females_data <- log(females_data+1)
dim(females_data)
females_data[1:4,1:4]
save(females_data,file = 'step1.1-A-females_hvg_matrix.Rdata')

#############################
# 6个发育时期RtSNE分析
#############################
## 针对上面的822个HVGs进行操作
female_sub_pca <- FactoMineR::PCA(
  t(females_data), 
  ncp = ncol(females_data), 
  graph=FALSE
)

## 然后挑选最显著的主成分，作为tSNE的输入
significant_pcs <- jackstraw::permutationPA(
  female_sub_pca$ind$coord, 
  B = 100, 
  threshold = 0.05, 
  verbose = TRUE, 
  seed = NULL
)$r
significant_pcs

## 使用上面jackstraw挑出的显著主成分进行tSNE
# 6个时期给定6个颜色
female_stagePalette <- 	c(
  "#2754b5", 
  "#8a00b0", 
  "#d20e0f", 
  "#f77f05", 
  "#f9db21",
  "#43f14b"
)
female_tsne <- run_plot_tSNE(
  pca=female_sub_pca,
  pc=significant_pcs,
  iter=5000,
  conditions=female_stages,
  colours=female_stagePalette
)
ggsave('step1.1-B-female_tSNE_by_stage.pdf')
save(female_tsne,file='step1.1-C-female_tsne.Rdata')

#############################
# 根据PCA结果进行层次聚类
# Hierarchical Clustering On Principle Components (HCPC)
#############################
## 使用9个显著主成分重新跑PCA
res.pca <- FactoMineR::PCA(
  t(females_data), 
  ncp = significant_pcs, 
  graph=FALSE
)
## 作者根据经验认为分成4群比较好解释，于是设置4
res.hcpc <- FactoMineR::HCPC(
  res.pca, 
  graph = FALSE,
  min=4
)
plot(res.hcpc, choice ="tree", cex = 0.6)
## 得到分群结果
female_clustering <- res.hcpc$data.clust$clust
table(female_clustering)
## 重新命名
female_clustering <- paste("C", female_clustering, sep="")
names(female_clustering) <- rownames(res.hcpc$data.clust)
# 将C1和C2调换位置[这些都是基于作者自己的理解]
female_clustering[female_clustering=="C1"] <- "C11"
female_clustering[female_clustering=="C2"] <- "C22"
female_clustering[female_clustering=="C22"] <- "C1"
female_clustering[female_clustering=="C11"] <- "C2"
table(female_clustering)
write.csv(female_clustering, file="step1.1-D-female_clustering.csv")

## 还是基于之前tSNE坐标，对4个cluster可视化
# 为4种cluster设置颜色
female_clusterPalette <- c(
  "#560047", 
  "#a53bad", 
  "#eb6bac", 
  "#ffa8a0"
)
# 还是利用原来tSNE坐标
head(female_t_sne)

# 作者DIY的函数
female_t_sne_new_clusters <- plot_tSNE(
  tsne=female_t_sne, 
  conditions=female_clustering, 
  colours= female_clusterPalette
)
ggsave('step1.1-E-tSNE_cluster.pdf')









