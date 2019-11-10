### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-差异分析-注释
### ---------------

rm(list=ls())
options(stringsAsFactors = F)
library(clusterProfiler)

#############################
# 先对monocle结果进行注释
#############################
load(file = 'step3.1-DEG-monocle_summary.Rdata')
## 得到差异基因名
de_genes <- subset(de_clusters, qval<0.05)
length(de_genes$genes)

##然后基因ID转换
entrez_genes <- bitr(de_genes$genes, fromType="SYMBOL", 
                     toType="ENTREZID", 
                     OrgDb="org.Mm.eg.db")
# 作者剔除掉一个基因名
entrez_genes <- entrez_genes[!entrez_genes$ENTREZID %in% "101055843",]
length(entrez_genes$SYMBOL)

# 因为有一些de_genes的SYMBOL没有对应的ENTREZID，因此看到少了100多个基因
# 然后，把存在ENTREZID的那些基因的基因名和cluster信息提取出来
de_gene_clusters <- de_genes[de_genes$genes %in% entrez_genes$SYMBOL,
                             c("genes", "cluster")]

# 保持de_gene_clusters$genes的顺序不变，将他的symbol变成entrez ID
de_gene_clusters <- data.frame(
  ENTREZID=entrez_genes$ENTREZID[entrez_genes$SYMBOL %in% de_gene_clusters$genes],
  cluster=de_gene_clusters$cluster
)

## 【必须】将差异基因对应到每个cluster
list_de_gene_clusters <- split(de_gene_clusters$ENTREZID, 
                               de_gene_clusters$cluster)

## 然后可以将4个cluster同时进行GO注释，并且放在一起（耗时！）
start_time <- Sys.time()
formula_res <- compareCluster(
  ENTREZID~cluster, 
  data=de_gene_clusters, 
  fun="enrichGO", 
  OrgDb="org.Mm.eg.db",
  ont		   = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.01,
  qvalueCutoff  = 0.05
)
end_time <- Sys.time()
(end_time - start_time)
# 可视化
pdf('step3.3-A-DEG_GO_each_cluster.pdf',width = 11,height = 6)
dotplot(formula_res, showCategory=5)
dev.off()

## 简化GO
# The simplified version of enriched result is more clear and give us a more comprehensive view of the whole story
start_time <- Sys.time()
lineage1_ego <- simplify(
  formula_res, 
  cutoff=0.5, 
  by="p.adjust", 
  select_fun=min
)
end_time <- Sys.time()
(end_time - start_time)
# Time difference of 1.588762 mins

pdf('step3.3-B-DEG_GO_each_cluster_simplified.pdf',width = 11,height = 6)
dotplot(lineage1_ego, showCategory=5)
dev.off()


# 保存结果
write.csv(formula_res@compareClusterResult, 
          file="step3.3-C-DEG_GO_each_cluster.csv")
write.csv(lineage1_ego@compareClusterResult, 
          file="step3.3-D-DEG_GO_each_cluster_simplified.csv")

#############################
# 以第一个GO:0140014为例，进行探索
#############################
## 首先获得GO:0140014中的基因
library(org.Mm.eg.db)
# 结果有3百多万个
go2gene <- toTable(org.Mm.egGO2ALLEGS)

# 其实看到里面的基因ID存在重复，一个go_id会对应多个同样的gene_id，那么就要
## 去重复，获取unique gene id
uni_gene <- unique(go2gene$gene_id[go2gene$go_id == 'GO:0140014'])
length(uni_gene)

## 然后拿到第一群的差异基因
c1_genes <- list_de_gene_clusters[['C1']]
length(c1_genes)
# 这里的1284比之前的1236多个几十个，说明存在几十个基因没有GO注释

## 找Cluster1在GO:0140014中的基因，其实就是找c1_genes和go2gene的交叉
overlap_genes <- intersect(c1_genes,uni_gene)
length(overlap_genes)








