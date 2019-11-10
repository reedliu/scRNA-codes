### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-11-10
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-聚类后的基因进一步注释
### ---------------
rm(list=ls())
options(stringsAsFactors = F)
library(clusterProfiler)

#############################
# 作者重新进行整合分组（17组变7组）
#############################
# 下载作者做好的整合数据：https://raw.githubusercontent.com/IStevant/XX-XY-mouse-gonad-scRNA-seq/master/data/female_lineages_DE_gene_pseudotime_clustered_annotated.csv
dyn_genes <- read.csv(file="../female_lineages_DE_gene_pseudotime_clustered_annotated.csv")
gene_names <- dyn_genes$Genes

#基因ID转换
entrez_genes <- bitr(gene_names, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")

# 提取有对应Entrez ID的变化基因，drop=FALSE确保返回数据框
gene_clusters <- dyn_genes[dyn_genes$Genes %in% entrez_genes$SYMBOL,,drop=FALSE]

# 准备进行GoSemSim
de_gene_clusters <- data.frame(
  ENTREZID=entrez_genes[!duplicated(entrez_genes$SYMBOL),"ENTREZID"],
  Gene_Clusters=gene_clusters$Gene.categories
)

# 进行富集分析
formula_res <- compareCluster(
  ENTREZID~Gene_Clusters, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont		   = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Run simplified GO enrochment analysis
lineage1_ego <- simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
)

write.csv(formula_res@compareClusterResult, 
          file="step5.4-A-female_compared_GO_term_DE_cluster.csv")
write.csv(lineage1_ego@compareClusterResult, 
          file="step5.4-B-female_compared_symplified_GO_term_DE_cluster.csv")

pdf(file="step5.4-C-female_GO_term_DE_genes_clusters.pdf", width=11, height=8)
dotplot(formula_res, showCategory=3)+ theme(aspect.ratio=0.8)
dotplot(lineage1_ego, showCategory=3)+ theme(aspect.ratio=2)
dev.off()
