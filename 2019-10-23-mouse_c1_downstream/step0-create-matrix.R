### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-23
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-创建表达矩阵
### ---------------
rm(list = ls()) 
options(warn=-1) 

#############################
# 下载并加载作者数据
#############################
# 作者做好的定量RPKM结果：https://github.com/IStevant/XX-XY-mouse-gonad-scRNA-seq/raw/master/data/female_rpkm.Robj
# 一会要用的基因列表：https://github.com/IStevant/XX-XY-mouse-gonad-scRNA-seq/blob/master/data/prot_coding.csv
# 作者给出的DIY-Function代码：https://raw.githubusercontent.com/IStevant/XX-XY-mouse-gonad-scRNA-seq/master/scripts/analysis_functions.R

load(file="female_rpkm.Robj")

## 去掉重复细胞
#（例如：同一个细胞建库两次，这里作者用“rep”进行了标记） 
grep("rep",colnames(female_rpkm))
colnames(female_rpkm)[256:257]
female_rpkm <- female_rpkm[,!colnames(female_rpkm) %in% grep("rep",colnames(female_rpkm), value=TRUE)]

## 只保留编码基因(去掉类似：X5430419D17Rik、BC003331等)
prot_coding_genes <- read.csv(file="prot_coding.csv", row.names=1)
females <- female_rpkm[rownames(female_rpkm) %in% as.vector(prot_coding_genes$x),]
save(females,file = 'female_rpkm.Rdata')
