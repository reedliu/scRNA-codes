### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-09-10
### Email: jieandze1314@gmail.com
### ---------------

# 目的：利用scRNA-seq包的表达矩阵测试几个R包寻找HVGs，然后画upsetR图看看不同方法的HVG的重合情况。

rm(list = ls())  
options(stringsAsFactors = F)

#############################
# 关于测试数据scRNAseq
#############################
# 包中内置了Pollen et al. 2014 的数据集（https://www.nature.com/articles/nbt.2967），到19年8月为止，已经有446引用量了。只不过原文完整的数据是 23730 features， 301 samples，这个包中只选取了4种细胞类型：pluripotent stem cells 分化而成的 neural progenitor cells (NPC，神经前体细胞) ，还有 GW16（radial glia，放射状胶质细胞） 、GW21（newborn neuron，新生儿神经元） 、GW21+3(maturing neuron，成熟神经元) 

library(scRNAseq)
## ----- Load Example Data -----
data(fluidigm)

# 得到RSEM矩阵
assay(fluidigm)  <- assays(fluidigm)$rsem_counts
ct <- floor(assays(fluidigm)$rsem_counts)
ct[1:4,1:4]

# 样本注释信息
pheno_data <- as.data.frame(colData(fluidigm))
table(pheno_data$Coverage_Type)
table(pheno_data$Biological_Condition)

#############################
# 利用Seurat V3
#############################
suppressMessages(library(Seurat))
packageVersion("Seurat") # 检查版本信息

seurat_sce <- CreateSeuratObject(counts = ct,
                          meta.data = pheno_data,
                          min.cells = 5, 
                          min.features = 2000, 
                          project = "seurat_sce")
seurat_sce <- NormalizeData(seurat_sce)
# 默认取前2000个
seurat_sce <- FindVariableFeatures(seurat_sce, selection.method = "vst", nfeatures=2000)
VariableFeaturePlot(seurat_sce)
seurat_hvg <- VariableFeatures(seurat_sce)
length(seurat_hvg)
head(seurat_hvg)

#############################
# 利用Monocle
#############################
library(monocle)
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
sample_ann <- pheno_data
# 然后转换为AnnotatedDataFrame对象
pd <- new("AnnotatedDataFrame",
          data=sample_ann)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
monocle_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
monocle_cds
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds)
disp_table <- dispersionTable(monocle_cds) 
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1) 
monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id) 
plot_ordering_genes(monocle_cds) 

monocle_hvg <- unsup_clustering_genes[order(unsup_clustering_genes$mean_expression, decreasing=TRUE),][,1]
length(monocle_hvg)
# 也取前2000个
monocle_hvg <- monocle_hvg[1:2000]
head(monocle_hvg)

#############################
# 利用scran
#############################
library(SingleCellExperiment)
library(scran)
scran_sce <- SingleCellExperiment(list(counts=ct))
scran_sce <- computeSumFactors(scran_sce) 
scran_sce <- normalize(scran_sce)
fit <- trendVar(scran_sce, parametric=TRUE,use.spikes=FALSE) 
dec <- decomposeVar(scran_sce, fit)

plot(dec$mean, dec$total, xlab="Mean log-expression", 
     ylab="Variance of log-expression", pch=16)
curve(fit$trend(x), col="dodgerblue", add=TRUE)

scran_df <- dec
scran_df <- scran_df[order(scran_df$bio, decreasing=TRUE),]
scran_hvg <- rownames(scran_df)[scran_df$FDR<0.1]
length(scran_hvg)
head(scran_hvg)

#############################
# 利用M3Drop
#############################
library(M3DExampleData)
# 需要提供表达矩阵(expr_mat)=》normalized or raw (not log-transformed) 
HVG <-M3Drop::BrenneckeGetVariableGenes(expr_mat=ct, spikes=NA, suppress.plot=FALSE, fdr=0.1, minBiolDisp=0.5, fitMeanQuantile=0.8)
M3Drop_hvg <- rownames(HVG)
length(M3Drop_hvg)
head(M3Drop_hvg)

#############################
# 自定义函数
#############################
# 来自：(Brennecke et al 2013 method) Extract genes with a squared coefficient of variation >2 times the fit regression 
if(F){
  getMostVarGenes <- function(
    data=data,				# RPKM matrix
    fitThr=1.5, 			# Threshold above the fit to select the HGV
    minMeanForFit=1			# Minimum mean gene expression level
  ){
    # data=females;fitThr=2;minMeanForFit=1	
    # Remove genes expressed in no cells
    data_no0 <- as.matrix(
      data[rowSums(data)>0,]
    )
    # Compute the mean expression of each genes
    meanGeneExp <- rowMeans(data_no0)
    names(meanGeneExp)<- rownames(data_no0)
    
    # Compute the squared coefficient of variation
    varGenes <- rowVars(data_no0)
    cv2 <- varGenes / meanGeneExp^2
    
    # Select the genes which the mean expression is above the expression threshold minMeanForFit
    useForFit <- meanGeneExp >= minMeanForFit
    
    # Compute the model of the CV2 as a function of the mean expression using GLMGAM
    fit <- glmgam.fit( cbind( a0 = 1, 
                              a1tilde = 1/meanGeneExp[useForFit] ), 
                       cv2[useForFit] )
    a0 <- unname( fit$coefficients["a0"] )
    a1 <- unname( fit$coefficients["a1tilde"])
    
    # Get the highly variable gene counts and names
    fit_genes <- names(meanGeneExp[useForFit])
    cv2_fit_genes <- cv2[useForFit]
    fitModel <- fit$fitted.values
    names(fitModel) <- fit_genes
    HVGenes <- fitModel[cv2_fit_genes>fitModel*fitThr]
    print(length(HVGenes))
    
    # Plot the result
    plot_meanGeneExp <- log10(meanGeneExp)
    plot_cv2 <- log10(cv2)
    plotData <-  data.frame(
      x=plot_meanGeneExp[useForFit],
      y=plot_cv2[useForFit],
      fit=log10(fit$fitted.values),
      HVGenes=log10((fit$fitted.values*fitThr))
    )
    p <- ggplot(plotData, aes(x,y)) +
      geom_point(size=0.1) +
      geom_line(aes(y=fit), color="red") +
      geom_line(aes(y=HVGenes), color="blue") +
      theme_bw() +
      labs(x = "Mean expression (log10)", y="CV2 (log10)")+
      ggtitle(paste(length(HVGenes), " selected genes", sep="")) +
      theme(
        axis.text=element_text(size=16),
        axis.title=element_text(size=16),
        legend.text = element_text(size =16),
        legend.title = element_text(size =16 ,face="bold"),
        legend.position= "none",
        plot.title = element_text(size=18, face="bold", hjust = 0.5),
        aspect.ratio=1
      )+
      scale_color_manual(
        values=c("#595959","#5a9ca9")
      )
    print(p)
    
    # Return the RPKM matrix containing only the HVG
    HVG <- data_no0[rownames(data_no0) %in% names(HVGenes),]
    return(HVG)
  }
}
library(statmod)
diy_hvg <- rownames(getMostVarGenes(ct))
length(diy_hvg)
head(diy_hvg)

#############################
# upsetR
#############################
require(UpSetR)
input <- fromList(list(seurat=seurat_hvg, monocle=monocle_hvg, 
                       scran=scran_hvg, M3Drop=M3Drop_hvg, diy=diy_hvg))
upset(input)









