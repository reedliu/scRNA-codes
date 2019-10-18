### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-15
### Email: jieandze1314@gmail.com
### Title: 10X scRNA免疫治疗-PBMC-细胞亚群的生物学命名
### ---------------

#############################
# 加载PBMC聚类结果
#############################
rm(list = ls()) 
options(warn=-1) 
suppressMessages(library(Seurat))
start_time <- Sys.time()
load('./patient1.PBMC.V2V3.output.Rdata')
end_time <- Sys.time()
end_time - start_time

colP<-c('green4', 
        'pink', 
        '#FF7F00', 
        'orchid', 
        '#99c9fb', 
        'dodgerblue2', 
        'grey30', 
        'yellow', 
        'grey60', 
        'grey', 
        'red', 
        '#FB9A99', 
        'black'
)

TSNEPlot(PBMC, 
         colors.use =  colP,
         do.label = T)
#############################
# 将细胞名称对应上去
#############################
## Seurat V3
a=read.table('celltype-patient1-PBMC.txt')
new.cluster.ids <- as.character(a[,2])
names(new.cluster.ids) <- levels(PBMC_V3)
PBMC_V3 <- RenameIdents(PBMC_V3, new.cluster.ids)

DimPlot(PBMC_V3, reduction = "tsne", label = TRUE, pt.size = 0.5, cols = colP) + NoLegend()

## Seurat V2
# # 先与表中第一列对应
# match(as.numeric(as.character(PBMC@ident)),a[,1])
# 
# # 第一点需要注意的是：为什么先用as.character后用as.numeric，而不是直接用as.numeric？
# # 原因就是PBMC@ident存储的是因子型变量，直接取只会得到它们的位置信息，而不是真实的分群信息
# head(as.numeric(PBMC@ident))
# head(as.character(PBMC@ident))
# head(as.numeric(as.character(PBMC@ident)))

# 第二点需要注意的是：match函数的规则是，A要在B中找到对应位置，那么就是 match(A,B)
labels=a[match(as.numeric(as.character(PBMC@ident)),a[,1]),2]
table(labels)
table(PBMC@ident)

PBMC@meta.data$labels=labels

TSNEPlot(PBMC, group.by = 'labels',
         colors.use =  colP,
         do.label = T)
# 修改颜色顺序
colP=colP[match(levels(as.factor(labels)),a[,2])]
TSNEPlot(PBMC, group.by = 'labels',
         colors.use =  colP,
         do.label = T)

#############################
# 按时间拆分分群结果
#############################
TimePoints = PBMC@meta.data$TimePoints
table(TimePoints)

# 取子集（以PBMC_RespD376时间点为例）
PBMC_RespD376 = SubsetData(PBMC,TimePoints =='PBMC_RespD376')
TSNEPlot(PBMC_RespD376, 
         colors.use = c('green4', 'pink', '#FF7F00', 'orchid', '#99c9fb', 'dodgerblue2', 'grey30', 'yellow', 'grey60', 'grey', 'red', '#FB9A99', 'black'),
         do.label = T)
# ggsave('PBMC_RespD376_PBMC_tSNE.pdf')

count_matrix=PBMC_RespD376_for_DEG@data
count_matrix[1:4,1:4]
cluster=PBMC_RespD376_for_DEG@ident
table(cluster)
save(count_matrix,cluster,
     file = 'patient1.PBMC_RespD376_for_DEG.Rdata')











