### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-差异分析-ROTS
### ---------------

#############################
# 准备表达矩阵和分群信息
#############################
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

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
# 使用ROTS
#############################
# ROTS包（Reproducibility-optimized test statistic），对每个亚群和其他几个亚群共同体进行比较
# ROTS可以使用RPKM值；它的速度可能会很慢
# 之前写过：https://github.com/reedliu/academic-me/blob/64d883639974750205e1b0ec4977d61c4e6211ca/content/post/cnposts/053-scRNA-learning-20.md 
# 差异分析重点就在：表达矩阵和分组信息

library(ROTS)
library(plyr)
load('../female_rpkm.Rdata')

## 配置最重要的分组信息
# 当前分组是这样
table(female_clustering)
# 我们只需要用数字来表示，于是可以用substring(x, start, stop)获取
groups<-female_clustering
groups<-as.numeric(substring(as.character(groups),2,2))
table(groups)
# 首先针对第1群和其他群比较(把其他亚群定义为234)
groups[groups!=1]<-234

## 然后配置表达矩阵
RPKM.full=females
ROTS_input<-RPKM.full[rowMeans(RPKM.full)>=1,]
ROTS_input<-as.matrix(log2(ROTS_input+1))

# 第1群
start_time <- Sys.time()
results_pop1 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop1<-data.frame(summary(results_pop1, fdr=1))
end_time <- Sys.time()
(end_time - start_time)

# 然后第2群
groups[groups!=2]<-134
start_time <- Sys.time()
results_pop2 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop2<-data.frame(summary(results_pop2, fdr=1))
end_time <- Sys.time()
(end_time - start_time)

# 然后第3群
groups[groups!=3]<-124
start_time <- Sys.time()
results_pop3 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop3<-data.frame(summary(results_pop3, fdr=1))
end_time <- Sys.time()
(end_time - start_time)

# 然后第4群
groups[groups!=4]<-123
start_time <- Sys.time()
results_pop4 = ROTS(data = ROTS_input, groups = groups , B = 1000 , K = 500 , seed = 1234)
summary_pop3<-data.frame(summary(results_pop4, fdr=1))
end_time <- Sys.time()
(end_time - start_time)

# 都得到以后，共同保存
save(summary_pop1,summary_pop2,summary_pop3,summary_pop4,
     file = 'step3.2-DEG-ROTS_summary.Rdata')








