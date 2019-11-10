### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-Seurat-marker基因可视化
### ---------------

#############################
# 准备工作
#############################
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)

# 要使用Seurat，最后加载原始表达矩阵
load('../female_count.Rdata')

# 发育时期获取
head(colnames(female_count))
female_stages <- sapply(strsplit(colnames(female_count), "_"), `[`, 1)
names(female_stages) <- colnames(female_count)
table(female_stages)

# 加载之前HCPC分群结果（代码在step1.1-tSNE-DIY.R）
cluster <- read.csv('../step1-female-RPKM-tSNE/step1.1-D-female_clustering.csv')
female_clustering=cluster[,2];names(female_clustering)=cluster[,1]
table(female_clustering)

## 作者选择的14个marker基因
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

#############################
# 使用Seurat V3
#############################
library(Seurat)
packageVersion('Seurat')
## 构建对象
sce_female <- CreateSeuratObject(counts = female_count, 
                                 project = "sce_female", 
                                 min.cells = 1, min.features = 0)
sce_female

## 添加样本注释信息
sce_female <- AddMetaData(object = sce_female, 
                          metadata = apply(female_count, 2, sum), 
                          col.name = 'nUMI_raw')
sce_female <- AddMetaData(object = sce_female, 
                          metadata = female_stages, 
                          col.name = 'female_stages')

## 数据归一化
sce_female <- NormalizeData(sce_female)
sce_female[["RNA"]]@data[1:3,1:3]

## 找差异基因HVGs
sce_female <- FindVariableFeatures(sce_female, 
                                   selection.method = "vst", 
                                   nfeatures = 2000)
## 数据标准化
sce_female <- ScaleData(object = sce_female, 
                        vars.to.regress = c('nUMI_raw'), 
                        model.use = 'linear', 
                        use.umi = FALSE)
## PCA降维
sce_female <- RunPCA(sce_female, 
                     features = VariableFeatures(object = sce_female))
## 降维后聚类
sce_female <- FindNeighbors(sce_female, dims = 1:20)
sce_female <- FindClusters(sce_female, resolution = 0.3)

## 继续tSNE非线性降维
sce_female_tsne <- RunTSNE(sce_female, dims = 1:9)
DimPlot(object = sce_female_tsne, reduction = "tsne")

# 小提琴图
pdf('step2.2-A-seurat3_VlnPlot.pdf', width=10, height=15)
VlnPlot(object = sce_female_tsne, features =  markerGenes , 
        pt.size = 0.2,ncol = 4)
dev.off()

# 基因表达量热图
pdf('step2.2-B-seurat3_FeaturePlot.pdf', width=10, height=15)
FeaturePlot(object = sce_female_tsne, features = markerGenes ,
            pt.size = 0.2,ncol = 3)
dev.off()

#############################
# 比较step2.1-作者代码和step2.2-Seurat的结果
#############################
# https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2019-10-24-025222.png
# 结果类似，那么如果我们自己做小提琴图呢？

# 其实就需要表达矩阵和分类信息
# 找代码：https://rpkgs.datanovia.com/ggpubr/reference/ggviolin.html

# 就画其中的Nr2f2基因
## 分类信息在此
group <- Seurat::Idents(sce_female)
## 表达矩阵在此
nr2f2 <- as.numeric(log(female_count['Nr2f2',]+1))

## boxplot画一个箱线图
boxplot(nr2f2~group)

## ggboxplot画一个箱线图
df <- data.frame(expr=nr2f2,
                 group=group)
female_clusterPalette <- c(
  "#560047", 
  "#a53bad", 
  "#eb6bac", 
  "#ffa8a0"
)
my_comparisons <- list( c("0", "1"), c("1", "2"), c("2", "3") )
ggboxplot(df, x = "group", y = "expr",
          color = "group", palette = female_clusterPalette)+ 
  stat_compare_means(comparisons = my_comparisons)

## ggviolin再画一个小提琴图
ggviolin(df, "group", "expr", fill = "group",
         palette = female_clusterPalette,
         add = "boxplot", add.params = list(fill = "white"))+
  stat_compare_means(comparisons = my_comparisons)

## ggstatsplot再画一个小提琴图
# https://github.com/IndrajeetPatil/ggstatsplot
# 基本上都是X轴分组，y轴连续变量

library(ggstatsplot)
ggstatsplot::ggbetweenstats(
  data = df,
  x = group,
  y = expr,
  messages = FALSE
) + # further modification outside of ggstatsplot
  ggplot2::coord_cartesian(ylim = c(3, 8)) +
  ggplot2::scale_y_continuous(breaks = seq(3, 8, by = 1))







