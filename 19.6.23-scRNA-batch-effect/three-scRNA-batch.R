### ---------------
###
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-06-23
### Email: jieandze1314@gmail.com
### Blog: https://www.jieandze1314.com/
### CAAS/AGIS/SDAU 
### Update Log: 2019-06-23  First version 
### ---------------

#########################
## 准备数据 GSE81076
#########################
rm(list = ls())  
options(stringsAsFactors = F)
gse81076.df <- read.table("GSE81076_D2_3_7_10_17.txt.gz", sep='\t', 
                          header=TRUE, row.names=1)
dim(gse81076.df)

# sub使用：https://en.wikibooks.org/wiki/R_Programming/Text_Processing#Detecting_the_presence_of_a_substring
# \\1返回^(D[0-9]+).*的小括号中内容
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse81076.df))
table(donor.names)

plate.id <- sub("^D[0-9]+(.*)_.*", "\\1", colnames(gse81076.df))
table(plate.id)

gene.symb <- gsub("__chr.*$", "", rownames(gse81076.df))

is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)

# id转换
library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.symb, keytype="SYMBOL", column="ENSEMBL")
# 除了标准id，还有很多ERCC需要更新到gene.ids
gene.ids[is.spike] <- gene.symb[is.spike]
# 去重复和NA（这里采用反选的方法）
keep <- !is.na(gene.ids) & !duplicated(gene.ids)
gse81076.df <- gse81076.df[keep,]
rownames(gse81076.df) <- gene.ids[keep]
summary(keep)

# 创建对象存储count和metadata
library(SingleCellExperiment)
sce.gse81076 <- SingleCellExperiment(list(counts=as.matrix(gse81076.df)),
                                     colData=DataFrame(Donor=donor.names, Plate=plate.id),
                                     rowData=DataFrame(Symbol=gene.symb[keep]))
# 重新设定一下ERCC位置
isSpike(sce.gse81076, "ERCC") <- grepl("^ERCC-", rownames(gse81076.df)) 
sce.gse81076  

# 质控和标准化
library(scater)
sce.gse81076 <- calculateQCMetrics(sce.gse81076, compact=TRUE)
QC <- sce.gse81076$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE))
# 去掉低质量的样本
discard <- low.lib | low.genes | high.spike
sce.gse81076 <- sce.gse81076[,!discard]
summary(discard)

# 标准化前戏--细胞聚类 
library(scran)
clusters <- quickCluster(sce.gse81076, min.mean=0.1)
table(clusters)
# 计算size factors
sce.gse81076 <- computeSumFactors(sce.gse81076, min.mean=0.1, clusters=clusters)
# 另外对ERCC也单独计算size factor(
sce.gse81076 <- computeSpikeFactors(sce.gse81076, general.use=FALSE)
# 最后计算标准化的log表达量，用作下游分析
sce.gse81076 <- normalize(sce.gse81076)

# 鉴定高变异基因HVGs
block <- paste0(sce.gse81076$Plate, "_", sce.gse81076$Donor)
fit <- trendVar(sce.gse81076, block=block, parametric=TRUE) 
dec <- decomposeVar(sce.gse81076, fit)
plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
OBis.spike <- isSpike(sce.gse81076)
points(dec$mean[is.spike], dec$total[is.spike], col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)
# 挑出来HVGs
dec.gse81076 <- dec
dec.gse81076$Symbol <- rowData(sce.gse81076)$Symbol
dec.gse81076 <- dec.gse81076[order(dec.gse81076$bio, decreasing=TRUE),]
head(dec.gse81076,2)
save(dec.gse81076,sce.gse81076,file = "gse81076.Rdata")

#########################
## 准备数据 GSE85241
#########################
rm(list = ls())  
options(stringsAsFactors = F)
gse85241.df <- read.table("GSE85241_cellsystems_dataset_4donors_updated.csv.gz", sep='\t', header=TRUE, row.names=1)
dim(gse85241.df)

# 提取meta信息（donor、plate）
donor.names <- sub("^(D[0-9]+).*", "\\1", colnames(gse85241.df))
table(donor.names)

plate.id <- sub("^D[0-9]+\\.([0-9]+)_.*", "\\1", colnames(gse85241.df))
table(plate.id)

# 提取基因、ERCC信息
gene.symb <- gsub("__chr.*$", "", rownames(gse85241.df))
is.spike <- grepl("^ERCC-", gene.symb)
table(is.spike)

# 基因转换
library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.symb, keytype="SYMBOL", column="ENSEMBL")
gene.ids[is.spike] <- gene.symb[is.spike]

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
gse85241.df <- gse85241.df[keep,]
rownames(gse85241.df) <- gene.ids[keep]
summary(keep)

# 创建单细胞对象
sce.gse85241 <- SingleCellExperiment(list(counts=as.matrix(gse85241.df)),
                                     colData=DataFrame(Donor=donor.names, Plate=plate.id),
                                     rowData=DataFrame(Symbol=gene.symb[keep]))
isSpike(sce.gse85241, "ERCC") <- grepl("^ERCC-", rownames(gse85241.df)) 
sce.gse85241  

# 质控和标准化
sce.gse85241 <- calculateQCMetrics(sce.gse85241, compact=TRUE)
QC <- sce.gse85241$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3)
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE))

# 去掉低质量细胞
discard <- low.lib | low.genes | high.spike
sce.gse85241 <- sce.gse85241[,!discard]
summary(discard)

# 聚类
clusters <- quickCluster(sce.gse85241, min.mean=0.1, method="igraph")
table(clusters)
# 标准化
sce.gse85241 <- computeSumFactors(sce.gse85241, min.mean=0.1, clusters=clusters)
summary(sizeFactors(sce.gse85241))
sce.gse85241 <- computeSpikeFactors(sce.gse85241, general.use=FALSE)
summary(sizeFactors(sce.gse85241, "ERCC"))
sce.gse85241 <- normalize(sce.gse85241)

# 鉴定HVGs
block <- paste0(sce.gse85241$Plate, "_", sce.gse85241$Donor)
fit <- trendVar(sce.gse85241, block=block, parametric=TRUE) 
dec <- decomposeVar(sce.gse85241, fit)
plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
is.spike <- isSpike(sce.gse85241)
points(dec$mean[is.spike], dec$total[is.spike], col="red", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

dec.gse85241 <- dec
dec.gse85241$Symbol <- rowData(sce.gse85241)$Symbol
dec.gse85241 <- dec.gse85241[order(dec.gse85241$bio, decreasing=TRUE),]
head(dec.gse85241,2)
save(dec.gse85241,sce.gse85241,file = "gse85241.Rdata")

#########################
## 准备数据 E-MTAB-5061
#########################
rm(list = ls())  
options(stringsAsFactors = F)
# 样本信息
header <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt", 
                     nrow=1, sep="\t", comment.char="")
header[1,1:4]
ncells <- ncol(header) - 1L
# 基因和count信息
col.types <- vector("list", ncells*2 + 2)
col.types[1:2] <- "character"
col.types[2+ncells + seq_len(ncells)] <- "integer"
e5601.df <- read.table("pancreas_refseq_rpkms_counts_3514sc.txt", 
                       sep="\t", colClasses=col.types)

gene.data <- e5601.df[,1:2]
e5601.df <- e5601.df[,-(1:2)]
colnames(e5601.df) <- as.character(header[1,-1])
dim(e5601.df)

# ERCC
is.spike <- grepl("^ERCC-", gene.data[,2])
table(is.spike)
# 基因转换
library(org.Hs.eg.db)
gene.ids <- mapIds(org.Hs.eg.db, keys=gene.data[,1], keytype="SYMBOL", column="ENSEMBL")
gene.ids[is.spike] <- gene.data[is.spike,2]

keep <- !is.na(gene.ids) & !duplicated(gene.ids)
e5601.df <- e5601.df[keep,]
rownames(e5601.df) <- gene.ids[keep]
summary(keep)

# metadata操作
metadata <- read.table("E-MTAB-5061.sdrf.txt", header=TRUE, 
                       sep="\t", check.names=FALSE)
m <- match(colnames(e5601.df), metadata$`Assay Name`)
stopifnot(all(!is.na(m)))
metadata <- metadata[m,]
donor.id <- metadata[["Characteristics[individual]"]]
table(donor.id)

# 单细胞对象
sce.e5601 <- SingleCellExperiment(list(counts=as.matrix(e5601.df)),
                                  colData=DataFrame(Donor=donor.id),
                                  rowData=DataFrame(Symbol=gene.data[keep,1]))
isSpike(sce.e5601, "ERCC") <- grepl("^ERCC-", rownames(e5601.df)) 
sce.e5601  
# 质控
sce.e5601 <- calculateQCMetrics(sce.e5601, compact=TRUE)
QC <- sce.e5601$scater_qc
low.lib <- isOutlier(QC$all$log10_total_counts, type="lower", nmad=3)
low.genes <- isOutlier(QC$all$log10_total_features_by_counts, type="lower", nmad=3) 
high.spike <- isOutlier(QC$feature_control_ERCC$pct_counts, type="higher", nmad=3)
low.spike <- isOutlier(QC$feature_control_ERCC$log10_total_counts, type="lower", nmad=2)
data.frame(LowLib=sum(low.lib), LowNgenes=sum(low.genes), 
           HighSpike=sum(high.spike, na.rm=TRUE), LowSpike=sum(low.spike))
# 舍弃低质量细胞
discard <- low.lib | low.genes | high.spike | low.spike
sce.e5601 <- sce.e5601[,!discard]
summary(discard)

# 聚类
clusters <- quickCluster(sce.e5601, min.mean=1, method="igraph")
table(clusters)


# 标准化
sce.e5601 <- computeSumFactors(sce.e5601, min.mean=1, clusters=clusters)
sce.e5601 <- computeSpikeFactors(sce.e5601, general.use=FALSE)
sce.e5601 <- normalize(sce.e5601)

# 可视化
donors <- sort(unique(sce.e5601$Donor))
# is.spike <- isSpike(sce.e5601)
par(mfrow=c(ceiling(length(donors)/2), 2), 
    mar=c(4.1, 4.1, 2.1, 0.1))
collected <- list() 
for (x in unique(sce.e5601$Donor)) {
  current <- sce.e5601[,sce.e5601$Donor==x]
  if (ncol(current)<2L) { next }
  current <- normalize(current)
  fit <- trendVar(current, parametric=TRUE) 
  dec <- decomposeVar(current, fit)
  plot(dec$mean, dec$total, xlab="Mean log-expression",
       ylab="Variance of log-expression", pch=16, main=x)
  points(fit$mean, fit$var, col="red", pch=16)
  curve(fit$trend(x), col="dodgerblue", add=TRUE)
  collected[[x]] <- dec
}

dec.e5601 <- do.call(combineVar, collected)
dec.e5601$Symbol <- rowData(sce.e5601)$Symbol
dec.e5601 <- dec.e5601[order(dec.e5601$bio, decreasing=TRUE),]
head(dec.e5601,3)
save(dec.e5601,sce.e5601,file = "e5601.Rdata")

#########################
## remove batch effetcs
#########################
rm(list = ls())  
options(stringsAsFactors = F)
load("gse81076.Rdata")
load("gse85241.Rdata")
load("e5601.Rdata")
# 选择前1000基因
top.e5601 <- rownames(dec.e5601)[seq_len(1000)]
top.gse85241 <- rownames(dec.gse85241)[seq_len(1000)]
top.gse81076 <- rownames(dec.gse81076)[seq_len(1000)]
# https://www.r-bloggers.com/intersect-for-multiple-vectors-in-r/
chosen <- Reduce(intersect, list(top.e5601, top.gse85241, top.gse81076))

# 添加gene symbol
symb <- mapIds(org.Hs.eg.db, keys=chosen, keytype="ENSEMBL", column="SYMBOL")
DataFrame(ID=chosen, Symbol=symb)

# 对每个数据集取前1000个HVGs
original <- list(logcounts(sce.e5601)[chosen,],
                 logcounts(sce.gse81076)[chosen,],
                 logcounts(sce.gse85241)[chosen,])
corrected <- do.call(mnnCorrect, c(original, list(k=20, sigma=0.1)))
str(corrected$corrected)

corrected$pairs

# 检验校正的作用
# omat是原始矩阵，mat是校正后的
omat <- do.call(cbind, original)
mat <- do.call(cbind, corrected$corrected)
# 将mat列名去掉
colnames(mat) <- NULL
sce <- SingleCellExperiment(list(original=omat, corrected=mat))
# 用lapply对三个列表进行循环操作，求列数，为了给rep设置一个重复值
colData(sce)$Batch <- rep(c("e5601", "gse81076", "gse85241"),
                          lapply(corrected$corrected, ncol))

sce

# 做个t-sne图
osce <- runTSNE(sce, exprs_values="original", set.seed=100)
ot <- plotTSNE(osce, colour_by="Batch") + ggtitle("Original")
csce <- runTSNE(sce, exprs_values="corrected", set.seed=100)
ct <- plotTSNE(csce, colour_by="Batch") + ggtitle("Corrected")
multiplot(ot, ct, cols=2)

# 利用4个marker基因检测
ct.gcg <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000115263") + ggtitle("Alpha cells (GCG)")
ct.ins <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000254647") + ggtitle("Beta cells (INS)")
ct.sst <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000157005") + ggtitle("Delta cells (SST)")
ct.ppy <- plotTSNE(csce, by_exprs_values="corrected", 
                   colour_by="ENSG00000108849") + ggtitle("PP cells (PPY)")
multiplot(ct.gcg, ct.ins, ct.sst, ct.ppy, cols=2)




