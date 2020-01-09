### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2020-01-03
### Email: jieandze1314@gmail.com
### Title: 改造inferCNV的热图
### ---------------

rm(list = ls())
options(stringsAsFactors = F)

#-- 1.Genes x Cells的表达矩阵(matrix) ----
expFile <- read.table('infercnv-heatmap/表达矩阵数据/expFile.txt')
expFile[1:4,1:4]

#-- 2.基因位置信息文件 ----
geneFile <- read.table('infercnv-heatmap/表达矩阵数据/geneFile.txt')
# 第一列对应第一个文件的行名，其余三列则是基因的位置
geneFile[1:4,1:4]
# 表达矩阵的行名，要和基因坐标文件的基因名顺序一致
identical(rownames(expFile),geneFile$V1)

#-- 3.样本注释信息文件 ----
groupFile <- read.table('infercnv-heatmap/表达矩阵数据/groupFiles.txt')
head(groupFile)
# 基因名不能有重复
dim(groupFile)
length(unique(groupFile$V1))
# 第二列是细胞的分组
table(groupFile$V2)

# 需要将groupFile和expFile的样本对应上
head(colnames(expFile))
head(groupFile$V1)
library(stringr)
colnames(expFile) <- str_replace(colnames(expFile),"\\.","-")
identical(groupFile$V1,colnames(expFile))

save(expFile,geneFile,groupFile,file = "step1-prepare.Rdata")



