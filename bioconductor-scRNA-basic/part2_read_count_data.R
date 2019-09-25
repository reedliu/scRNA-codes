### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date From: 2019-09-13
### Email: jieandze1314@gmail.com
### Source From： https://bioconductor.org/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/reads.html
### ---------------

#############################
# 首先下载、解压数据
#############################
BiocManager::install("BiocFileCache")
library(BiocFileCache)
# browseVignettes('BiocFileCache') 查看帮助文档
# 首先自定义一个存放数据的目录
bfc <- BiocFileCache("raw_data", ask = FALSE)
# file.path可以是本地路径或者网络链接（直接下载到BiocFileCache指定的目录）
lun.zip <- bfcrpath(bfc, 
                    file.path("https://www.ebi.ac.uk/arrayexpress/files",
                              "E-MTAB-5522/E-MTAB-5522.processed.1.zip"))
lun.sdrf <- bfcrpath(bfc, 
                     file.path("https://www.ebi.ac.uk/arrayexpress/files",
                               "E-MTAB-5522/E-MTAB-5522.sdrf.txt"))
unzip(lun.zip, exdir=tempdir())


#############################
# 数据设置
#############################
# ------ 加载表达矩阵 ---------
plate1 <- read.delim(file.path(tempdir(), "counts_Calero_20160113.tsv"), 
                     header=TRUE, row.names=1, check.names=FALSE)
plate1[1:2,1:2]
dim(plate1)
plate2 <- read.delim(file.path(tempdir(), "counts_Calero_20160325.tsv"), 
                     header=TRUE, row.names=1, check.names=FALSE)
# 提取基因长度
gene.lengths <- plate1$Length 
plate1 <- as.matrix(plate1[,-1]) 
plate2 <- as.matrix(plate2[,-1])
# 把两个细胞板的表达矩阵按列整合起来，但首先要检测它们的行名是否一致
stopifnot(identical(rownames(plate1), rownames(plate2)))
all.counts <- cbind(plate1, plate2)

# ------创建SCE对象 ------
suppressMessages(library(SingleCellExperiment))
sce <- SingleCellExperiment(list(counts=all.counts))
rowData(sce)$GeneLength <- gene.lengths
sce

# ------提取ERCC------
isSpike(sce, "ERCC") <- grepl("^ERCC", rownames(sce))
summary(isSpike(sce, "ERCC"))

# 除了ERCC，还有SIRV（舍去）
is.sirv <- grepl("^SIRV", rownames(sce))
summary(is.sirv)
sce <- sce[!is.sirv,] 

# ------ 添加整合细胞注释信息 ---------
metadata <- read.delim(lun.sdrf, check.names=FALSE, header=TRUE)
metadata[1:3,1:3]

m <- match(colnames(sce), metadata[["Source Name"]]) 
# 这样m就存储了一系列的位置，将表达矩阵的列和细胞注释信息的行顺序一一对应起来
stopifnot(all(!is.na(m))) # 检查是否完整（包含所有的细胞）
metadata <- metadata[m,]
# 这样就保证了我们只取到和表达矩阵相关的细胞注释，其他无用信息可以去除
head(colnames(metadata))

## 先添加细胞来源
colData(sce)$Plate <- factor(metadata[["Factor Value[block]"]])
## 然后添加细胞表型
pheno <- metadata[["Factor Value[phenotype]"]]
levels(pheno) <- c("induced", "control")
colData(sce)$Oncogene <- pheno
## 最后看看新增的细胞表型数据
table(colData(sce)$Oncogene, colData(sce)$Plate)
# save(sce,file = 'bioinfoplanet_sce.Rdata')
# load('bioinfoplanet_sce.Rdata')

# ------ 添加整合基因注释信息 ---------
library(org.Mm.eg.db)
symb <- mapIds(org.Mm.eg.db, keys=rownames(sce), 
               keytype="ENSEMBL", column="SYMBOL")
rowData(sce)$ENSEMBL <- rownames(sce)
rowData(sce)$SYMBOL <- symb
head(rowData(sce))

## scater去重、填NA
library(scater)
rownames(sce) <- uniquifyFeatureNames(rowData(sce)$ENSEMBL, 
                                      rowData(sce)$SYMBOL)
head(rownames(sce))

## 添加基因所在染色体位置信息
if(!require("TxDb.Mmusculus.UCSC.mm10.ensGene")) 
  BiocManager::install("TxDb.Mmusculus.UCSC.mm10.ensGene",
                       update = F,ask = F)

# 可以参考：http://bioinfoblog.it/tag/txdb/
library(TxDb.Mmusculus.UCSC.mm10.ensGene)
# TxDb storing transcript annotations（包括CDS、Exon、Transcript的start、end、chr-location、strand）
columns(TxDb.Mmusculus.UCSC.mm10.ensGene)
# 返回一个GRanges对象，使用transcripts()查看其中的内容
head(transcripts(TxDb.Mmusculus.UCSC.mm10.ensGene,columns=c('CDSCHROM')))
# 补充：如果是想看org.db其中的内容
if(F){
  keytypes(org.Hs.eg.db)
  head(mappedkeys(org.Hs.egENSEMBL))
}

location <- mapIds(TxDb.Mmusculus.UCSC.mm10.ensGene, 
                   keys=rowData(sce)$ENSEMBL, 
                   column="CDSCHROM", keytype="GENEID")
rowData(sce)$CHR <- location
summary(location=="chrM")

#############################
# 细胞质控
#############################
# -----首先得到一些质控指标---------
# 之前ERCC已经使用isSpike()添加到sce中，这里再将另一个影响因素Mt添加进来
mito <- which(rowData(sce)$CHR=="chrM")
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=mito))
head(colnames(colData(sce)), 10)

## 用scater的multiplot()作图
# 首先得到细胞类型和批次的组合信息
sce$PlateOnco <- paste0(sce$Oncogene, ".", sce$Plate)
multiplot(
  plotColData(sce, y="total_counts", x="PlateOnco"),
  plotColData(sce, y="total_features_by_counts", x="PlateOnco"),
  plotColData(sce, y="pct_counts_ERCC", x="PlateOnco"),
  plotColData(sce, y="pct_counts_Mt", x="PlateOnco"),
  cols=2)
# 两两指标之间查看
par(mfrow=c(1,3))
plot(sce$total_features_by_counts, sce$total_counts/1e6, xlab="Number of expressed genes",
     ylab="Library size (millions)")
plot(sce$total_features_by_counts, sce$pct_counts_ERCC, xlab="Number of expressed genes",
     ylab="ERCC proportion (%)")
plot(sce$total_features_by_counts, sce$pct_counts_Mt, xlab="Number of expressed genes",
     ylab="Mitochondrial proportion (%)")

# -----为每个质控指标鉴定离群点---------
# 一般也就是去掉低文库大小、低表达基因、高spike-in（返回与细胞数量等长的逻辑值）
libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", 
                          log=TRUE, batch=sce$PlateOnco)
feature.drop <- isOutlier(sce$total_features_by_counts, 
                          nmads=3, type="lower", 
                          log=TRUE, batch=sce$PlateOnco)
spike.drop <- isOutlier(sce$pct_counts_ERCC, nmads=3, type="higher",
                        batch=sce$PlateOnco)
# 可以查看isOutlier设置的阈值
attr(libsize.drop, "thresholds")
attr(spike.drop, "thresholds")
# 然后过滤
keep <- !(libsize.drop | feature.drop | spike.drop)
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop),
           BySpike=sum(spike.drop), Remaining=sum(keep))

sce$PassQC <- keep
# 保存一份QC后、过滤前的数据，方便日后加载
saveRDS(sce, file="416B_preQC.rds")
sce <- sce[,keep]
dim(sce)

# -----（可选）检查抛弃的离群点---------
# 初步做个散点图
if(F){
  # 计算丢弃和保留的细胞平均表达量
  library(SingleCellExperiment)
  sce.full.416b <- readRDS("416B_preQC.rds")
  
  library(scater)
  lost <- calcAverage(counts(sce.full.416b)[,!sce.full.416b$PassQC])
  kept <- calcAverage(counts(sce.full.416b)[,sce.full.416b$PassQC])
  
  # 在上面得到的平均值中，将每个数都与平均值中（除0以外）最小的数进行比较，取最大的那个值作为最终的平均值
  capped.lost <- pmax(lost, min(lost[lost>0]))
  capped.kept <- pmax(kept, min(kept[kept>0]))
  
  plot(capped.lost, capped.kept, xlab="Average count (discarded)", 
       ylab="Average count (retained)", log="xy", pch=16)
  is.spike <- isSpike(sce.full.416b)
  points(capped.lost[is.spike], capped.kept[is.spike], col="red", pch=16)
  is.mito <- rowData(sce.full.416b)$is_feature_control_Mt
  points(capped.lost[is.mito], capped.kept[is.mito], col="dodgerblue", pch=16)
}
# 再仔细检查一下logFC
if(F){
  library(edgeR)
  # 重点在于DGEList中group的设置：数小的是control，也就是logFC的分母
  y <- cbind(lost, kept)
  y <- DGEList(y, group=c(2,1))
  design <- model.matrix(~group, data=y$samples)
  predlfc <- predFC(y, design)[,2]
  info <- data.frame(logFC=predlfc, Lost=lost, Kept=kept, 
                     row.names=rownames(sce.full.416b))
  head(info[order(info$logFC, decreasing=TRUE),], 20)
}
# 其他一些质控方法，例如检测PCA
if(F){
  sce.tmp <- runPCA(sce.full.416b, use_coldata=TRUE, 
                    detect_outliers=TRUE)
  table(sce.tmp$outlier)
}

## Date：2019-09-22
#############################
# 细胞周期推断
#############################
set.seed(100)
library(scran)
mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                package="scran"))
system.time(assignments <- cyclone(sce, mm.pairs, 
                       gene.names=rowData(sce)$ENSEMBL))

# 具体规则是：如果一个细胞在G1中得分大于0.5，并且它高于G2/M的得分，那么这个细胞就被划分到G1期；如果细胞在G2/M中得分大于0.5，并且高于G1的得分，那么它就划为G2/M期；如果细胞的G1、G2/M得分都不大于0.5，那么它就划为S期。
plot(assignments$score$G1, assignments$score$G2M, 
     xlab="G1 score", ylab="G2/M score", pch=16)

sce$phases <- assignments$phases
table(sce$phases)
# 热图
if(F){
  library(pheatmap)
  # 取差异前100基因
  cg=names(tail(sort(apply(assay(sce),1,sd)),50))
  # 矩阵归一化
  n=t(scale(t(assay(sce)[cg,])))
  # 原来的样本注释信息 df中包含了 g、plate  、n_g、all信息，现在新增phases信息
  df=data.frame(plate=sce$Plate,cellcycle=assignments$phases )
  rownames(df)=colnames(n)
  pheatmap(n,show_colnames =F,show_rownames = F,
           annotation_col=df)
  dev.off()
}
save(sce,assignments,file = '416B_cell_cycle.Rdata')

#############################
# 基因层面检查
#############################
# -----检查高表达基因---------
fontsize <- theme(axis.text=element_text(size=12), 
                  axis.title=element_text(size=16))
plotHighestExprs(sce, n=50) + fontsize

# -----过滤低丰度基因之平均表达量---------
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey80", 
     xlab=expression(Log[10]~"average count"))
demo.keep <- ave.counts >= 1
filtered.sce <- sce[demo.keep,]
summary(demo.keep)

# -----过滤低丰度基因之细胞数量---------
num.cells <- nexprs(sce, byrow=TRUE)
smoothScatter(log10(ave.counts), num.cells, 
              ylab="Number of cells", 
              xlab=expression(Log[10]~"average count"))
to.keep <- num.cells > 0
sce <- sce[to.keep,]
summary(to.keep)

#############################
# 对细胞文库差异进行normalization
#############################
# -----去卷积法根据内源基因计算size factor---------
sce <- computeSumFactors(sce)
summary(sizeFactors(sce))
# 看看计算的size factor与lib size的关系
if(F){
  plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
       xlab="Library size (millions)", ylab="Size factor",
       col=c("red", "black")[sce$Oncogene], pch=16)
  # 右下角添加图例，对数据中两种细胞类型进行了区分
  legend("bottomright", col=c("red", "black"), pch=16, cex=1.2,
         legend=levels(sce$Oncogene))
}

# -----根据spike-in转录本单独计算size factor---------
sce <- computeSpikeFactors(sce, type="ERCC", general.use=FALSE)
# -----将size factor应用到normalization这一步---------
sce <- normalize(sce)

#############################
# 模拟基因表达中的技术噪音
#############################
# -----trendVar得到技术噪音方差--------
var.fit <- trendVar(sce, parametric=TRUE, block=sce$Plate,
                    loess.args=list(span=0.3))
# -----decomposeVar得到生物因素方差--------
var.out <- decomposeVar(sce, var.fit)
head(var.out,2)

# -----可视化 Mean-Variance plot--------
plot(var.out$mean, var.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
     ylab="Variance of log-expression")
curve(var.fit$trend(x), col="dodgerblue", lwd=2, add=TRUE)
cur.spike <- isSpike(sce)
points(var.out$mean[cur.spike], var.out$total[cur.spike], col="red", pch=16)

# -----检验生物因素方差结果是否由离群点造成--------
chosen.genes <- order(var.out$bio, decreasing=TRUE)[1:10]
plotExpression(sce, features=rownames(var.out)[chosen.genes]) + fontsize

## Date：2019-09-23
#############################
# 去除批次效应
#############################
library(limma)
# -----修改原来的实验比较分组--------
# 谁在前表示谁作为对照
sce$Oncogene <- relevel(sce$Oncogene,ref = 'control')
# -----检查好分组信息后去除批次--------
assay(sce, "corrected") <- removeBatchEffect(
  logcounts(sce), 
  design=model.matrix(~sce$Oncogene), 
  batch=sce$Plate
  )
assayNames(sce)

#############################
# PCA去除表达量中的技术噪音
#############################
sce <- denoisePCA(sce, technical=var.out, assay.type="corrected")
dim(reducedDim(sce, "PCA")) 

#############################
# 低维空间可视化
#############################
# -----PCA--------
# 根据细胞生物学分类
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
               colour_by="Oncogene") + fontsize
# 根据plate信息画PCA=》检验批次效应去除效果
plotReducedDim(sce, use_dimred="PCA", ncomponents=3, 
               colour_by="Plate") + fontsize

# -----tSNE--------
# 比较多个Perplexity参数
if(T){
  set.seed(100)
  out5 <- plotTSNE(sce, 
                   run_args=list(use_dimred="PCA", perplexity=5),
                   colour_by="Oncogene") + 
    fontsize + ggtitle("Perplexity = 5")
  
  set.seed(100)
  out10 <- plotTSNE(sce, 
                    run_args=list(use_dimred="PCA", perplexity=10),
                    colour_by="Oncogene") + 
    fontsize + ggtitle("Perplexity = 10")
  
  set.seed(100)
  out20 <- plotTSNE(sce, 
                    run_args=list(use_dimred="PCA", perplexity=20),
                    colour_by="Oncogene") + 
    fontsize + ggtitle("Perplexity = 20")
  
  multiplot(out5, out10, out20, cols=3)
  
}

# 好像perplexity=20效果还不错，因此下面继续跑runTSNE时就可以设置perplexity=20
set.seed(100)
sce <- runTSNE(sce, use_dimred="PCA", perplexity=20)
reducedDimNames(sce)

#############################
# 自己进行细胞聚类
#############################
# ----降维结果+层次聚类=》cluster--------
pcs <- reducedDim(sce, "PCA")
my.dist <- dist(pcs)
my.tree <- hclust(my.dist, method="ward.D2")
library(dynamicTreeCut)
my.clusters <- unname(cutreeDynamic(my.tree, 
                                    distM=as.matrix(my.dist), 
                                    minClusterSize=10, verbose=0))
# 检查与批次之间关系（结果无关）
table(my.clusters, sce$Plate)
# 检查与生物分类之间关系（结果有关）
table(my.clusters, sce$Oncogene)

# 得到cluster后，用tsne结果进行映射
sce$cluster <- factor(my.clusters)
plotTSNE(sce, colour_by="cluster") + fontsize

# 检查分群结果
library(cluster)
clust.col <- scater:::.get_palette("tableau10medium") # hidden scater colours
sil <- silhouette(my.clusters, dist = my.dist)
sil.cols <- clust.col[ifelse(sil[,3] > 0, sil[,1], sil[,2])]
sil.cols <- sil.cols[order(-sil[,1], sil[,3])]
plot(sil, main = paste(length(unique(my.clusters)), "clusters"), 
     border=sil.cols, col=sil.cols, do.col.sort=FALSE) 

# ----每个cluster检测marker基因--------
markers <- findMarkers(sce, my.clusters, block=sce$Plate)
# 看一下第1个cluster的top10差异基因
marker.set <- markers[["1"]]
head(marker.set, 10)

write.table(marker.set, file="416B_marker_1.tsv", sep="\t", 
            quote=FALSE, col.names=NA)

# 热图
top.markers <- rownames(marker.set)[marker.set$Top <= 10]
plotHeatmap(sce, features=top.markers, columns=order(sce$cluster), 
            colour_columns_by=c("cluster", "Plate", "Oncogene"),
            cluster_cols=FALSE, center=TRUE, symmetric=TRUE, 
            zlim=c(-5, 5),
            show_colnames = F) 


#############################
# 最后保存一下数据
#############################
saveRDS(file="416B_data.rds", sce)


