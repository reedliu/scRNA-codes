### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-差异分析-ROTS结果画热图
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
# 加载ROTS 的差异分析结果
load('step3.2-DEG-ROTS_summary.Rdata')
head(summary_pop1)
head(summary_pop2)
head(summary_pop3)
head(summary_pop4)

# 6个发育时期获取
head(colnames(females))
female_stages <- sapply(strsplit(colnames(females), "_"), `[`, 1)
names(females) <- colnames(females)
table(female_stages)
# 4个cluster获取
cluster <- read.csv('../step1-female-RPKM-tSNE/step1.1-D-female_clustering.csv')
female_clustering=cluster[,2];names(female_clustering)=cluster[,1]
table(female_clustering)

RPKM.full=females


#############################
# 每个亚群挑选top18的差异基因绘制热图
#############################
population_subset<-c(rownames(summary_pop1[summary_pop1$ROTS.statistic<0,])[1:18],
                     rownames(summary_pop2[summary_pop2$ROTS.statistic<0,])[1:18],
                     rownames(summary_pop3[summary_pop3$ROTS.statistic<0,])[1:18],
                     rownames(summary_pop4[summary_pop4$ROTS.statistic<0,])[1:18])

RPKM_heatmap<-RPKM.full[population_subset,]
# 将列名重新按照4个cluster排序
RPKM_heatmap<-RPKM_heatmap[,order(female_clustering)]
RPKM_heatmap<-log2(RPKM_heatmap+1)

popul.col<-sort(female_clustering);table(popul.col)
popul.col<-replace(popul.col, popul.col=='C1',"#1C86EE" )
popul.col<-replace(popul.col, popul.col=='C2',"#00EE00" )
popul.col<-replace(popul.col, popul.col=='C3',"#FF9912" )
popul.col<-replace(popul.col, popul.col=='C4',"#FF3E96" )
library(gplots)

pdf("step3.5-A-ROST-DEG-heatmap.pdf")
heatmap.2(as.matrix(RPKM_heatmap),
          ColSideColors = as.character(popul.col), tracecol = NA, 
          dendrogram = "none",col=bluered, labCol = FALSE, 
          scale="none", key = TRUE, symkey = F, symm=F,  key.xlab = "", 
          key.ylab = "", density.info = "density", 
          key.title = "log2(RPKM+1)", keysize = 1.2, denscol="black", Colv=FALSE)
dev.off()




















