### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-谱系发育基因可视化
### ---------------
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")

#############################
# 载入之前结果
#############################
# 最初我们知道细胞有6个时期（就是取样的6个时间点）；然后进行聚类发现这些细胞能分成4个cluster（意思就是虽然是一个时间点取的细胞，依然可能属于不同类型）；后来进行谱系推断，又增加了一个细胞属性（就是不同的发育轨迹）
load('../female_rpkm.Rdata')
load(file = 'step4.1-diffusionMap_output.Rdata')
load(file = 'step4.2-female_pseudotime.Rdata')

#############################
# 对谱系推断结果进行归一化
#############################
# 目的就是让两条轨迹可以比较，采用的方法就是每条轨迹的每个值分别除以各自的最大值
## 第一条
pseudotime_lin <- female_pseudotime[,"curve1"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin1_percent <- (pseudotime_lin*100)/max_pseudotime
## 第二条
pseudotime_lin <- female_pseudotime[,"curve2"]
max_pseudotime <- max(pseudotime_lin, na.rm = TRUE)
pseudotime_lin2_percent <- (pseudotime_lin*100)/max_pseudotime
# 现在female_pseudotime中的两条轨迹结果都在0-100之间了
female_pseudotime[,"curve1"] <- pseudotime_lin1_percent
female_pseudotime[,"curve2"] <- pseudotime_lin2_percent

save(female_pseudotime,file='step4.3-female-psudotime-percent.Rdata')
# 不同于stage和cluster两种细胞属性，这个发育谱系属性是连续型的 
# 既然细胞是按照一定顺序排列的，那么就会有一些基因表达量会跟着这个连续变量进行变化
# 而过去只有离散型的分类变量，因此只能先通过差异分析得到结果，然后对结果去注释

#############################
# 可视化
#############################
# 作者包装的代码非常复杂

## 给一个颜色
source('../analysis_colors.R')

## 包装好的函数
# 做第一个谱系的某个基因
plot_smoothed_gene_per_lineage(
  rpkm_matrix=females, # RPKM表达矩阵
  pseudotime=female_pseudotime,  #谱系推断结果
  lin=c(1), # 对第一个谱系操作
  gene="Amhr2",  #画Amhr2基因变化
  stages=female_stages, # 发育时间点分类
  clusters=female_clustering, # cluster分类
  stage_colors=female_stagePalette,
  cluster_colors=female_clusterPalette,
  lineage_colors=female_clusterPalette2
)

# 做某个基因在两个谱系中的变化
plot_smoothed_gene_per_lineage(
  rpkm_matrix=females, 
  pseudotime=female_pseudotime, 
  lin=c(1,2),
  gene="Amhr2", 
  stages=female_stages, 
  clusters=female_clustering, 
  stage_colors=female_stagePalette,
  cluster_colors=female_clusterPalette,
  lineage_colors=female_clusterPalette2
)

# 如果要进行批量作图
gene_list <- c(
  "Sall4",
  "Sox11",
  "Gata4",
  "Lgr5",
  "Runx1",
  "Foxl2",
  "Hey2",
  "Wnt5a",
  "Pdgfra",
  "Nr2f2",
  "Sfrp1",
  "Ifitm3",
  "Ptch1",
  "Wnt4",
  "Rspo1",
  "Cdkn1b",
  "Gli1",
  "Tcf21",
  "Nr0b1",
  "Nr0b2",
  "Nr5a1",
  "Nr6a1"
)

plot_smoothed_genes <- function(genes, lin){
  female_clusterPalette2 <- c("#ff6663", "#3b3561")
  for (gene in genes){
    plot_smoothed_gene_per_lineage(
      rpkm_matrix=females, 
      pseudotime=female_pseudotime, 
      lin=lin,
      gene=gene, 
      stages=female_stages, 
      clusters=female_clustering, 
      stage_colors=female_stagePalette,
      cluster_colors=female_clusterPalette,
      lineage_colors=female_clusterPalette2
    )
  }
}

pdf("step4.3-interesting_genes_in_lineage.pdf", width=4, height=4)
plot_smoothed_genes(gene_list, 1) # plot only lineage 1
plot_smoothed_genes(gene_list, 2) # plot only lineage 2
plot_smoothed_genes(gene_list, c(1,2)) # plot the two moleages in the same graph to see the divergence
dev.off()









