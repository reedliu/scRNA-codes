### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2020-01-03
### Email: jieandze1314@gmail.com
### Title: 改造inferCNV的热图
### ---------------

rm(list = ls())
options(stringsAsFactors = F)

load('step1-prepare.Rdata')
library(infercnv)
#--- 再次检查----------
expFile[1:4,1:4]
head(groupFile)
identical(groupFile$V1,colnames(expFile))# 与exp列名对应
head(geneFile)
sum(duplicated(geneFile$V1)) #没有重复基因名
identical(geneFile$V1,rownames(expFile))# 与exp行名对应

# 需要将文件保存为txt文件，然后读入，expFile可以直接读取
write.table(groupFile,file = "groupFile.txt",sep = '\t',quote = F,col.names = F,row.names = F)
write.table(geneFile,file = "geneFile.txt",sep = '\t',quote = F,col.names = F,row.names = F)

#---- 两步走 ------------
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=expFile,
                                    annotations_file="groupFile.txt",
                                    delim="\t",
                                    gene_order_file= "geneFile.txt",
                                    ref_group_names=c("WT"))  
if(!dir.exists("output_dir")){
  dir.create("output_dir")
}
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir="output_dir", 
                             cluster_by_groups=TRUE,#聚类（是否根据细胞注释文件的分组对肿瘤细胞进行分群）
                             denoise=TRUE, #去噪
                             HMM=F) # 是否基于HMM预测CNV
## Tips：
## 设置HMM=TRUE 的计算时间会长于HMM=FALSE,因此可以先设置HMM=FALSE快速查看结果


#--- 部分有用的结果------
# 全部的解释在：https://github.com/broadinstitute/inferCNV/wiki/Output-Files
# 当前版本（2019-3-13）作者也认为输出的结果太冗余，未来会进行整合

# infercnv.preliminary.png : 初步的inferCNV展示结果（未经去噪或HMM预测）
# 
# infercnv.png : 最终inferCNV产生的去噪后的热图.
# 
# infercnv.references.txt : 正常细胞矩阵.
# 
# infercnv.observations.txt : 肿瘤细胞矩阵.
# 
# infercnv.observation_groupings.txt : 肿瘤细胞聚类后的分组关系.
# 
# infercnv.observations_dendrogram.txt : NEWICK格式，展示肿瘤细胞间的层次关系.






