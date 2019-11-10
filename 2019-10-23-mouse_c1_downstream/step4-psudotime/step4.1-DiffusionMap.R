### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-谱系推断-DiffusionMap
### ---------------
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

#############################
# 准备表达矩阵、HVGs、分群信息
#############################
# 表达矩阵
load('../female_rpkm.Rdata')
dim(females)

# HVGs
load('../step1-female-RPKM-tSNE/step1.1-A-females_hvg_matrix.Rdata')
dim(females_data)

# 6个发育时期获取
head(colnames(females))
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(female_stages) <- colnames(females)
table(female_stages)

# 4个cluster获取
cluster <- read.csv('../step1-female-RPKM-tSNE/step1.1-D-female_clustering.csv')
female_clustering=cluster[,2];names(female_clustering)=cluster[,1]
table(female_clustering)

#############################
# 进行DiffusionMap
#############################
# 包装的代码很简单
female_dm <- run_diffMap(
  females_data, 
  female_clustering,
  sigma=15
)

save(female_dm,females_data,female_clustering,female_stages,
     file = 'step4.1-diffusionMap_output.Rdata')
#############################
# 作图探索
#############################
# 画出特征值，这个很像PCA的碎石图screeplots或者elbowplot，也是看拐点
plot_eigenVal(
  dm=female_dm
)

# 探索4个分群
female_clusterPalette <- c(
  C1="#560047",
  C2="#a53bad", 
  C3="#eb6bac", 
  C4="#ffa8a0"
)

plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_clustering, 
  colour=female_clusterPalette
)

# 探索6个发育时间
female_stagePalette=c(
  E10.5="#2754b5", 
  E11.5="#8a00b0", 
  E12.5="#d20e0f", 
  E13.5="#f77f05", 
  E16.5="#f9db21",
  P6="#43f14b"
)

plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_stages, 
  colour=female_stagePalette
)

# 可以看到，这个函数主要就是选取了前3个DC成分，做了三维空间的映射，然后把点的颜色分别按照cluster和stage两种不同的属性上色




