### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-08-22
### Email: jieandze1314@gmail.com
### Blog: www.jieandze1314.com
### Update Log: 2019-08-22 Monocle-3
### From: https://cole-trapnell-lab.github.io/monocle3/monocle3_docs/
### ---------------

#############################
# 配置镜像及安装包
#############################
file.edit('~/.Rprofile')
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")) #对应清华源
options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/") #对应中科大源

bioc_pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
               'limma', 'S4Vectors', 'SingleCellExperiment',
               'SummarizedExperiment')
for (bpkg in bioc_pkgs){
  if (! require(bpkg,character.only=T) ) {
    BiocManager::install(bpkg,ask = F,update = F)
    require(bpkg,character.only=T)
  }
}
# (可选)如果要在细胞聚类时设定分辨率参数(resolution)，就要安装一个python包louvain
# 而这个python包安装又需要virtualenv
if(F){
  install.packages('reticulate')
  library(reticulate)
  # https://cran.r-project.org/web/packages/reticulate/vignettes/python_packages.html
  ## 先使用virtualenv方法：
  # Rstudio terminal :
  # /usr/bin/python2.7 -m pip install --upgrade --user virtualenv
  virtualenv_create("r-reticulate")
  virtualenv_install("r-reticulate", "louvain")
  louvain <- import("louvain")
  
  ## 再使用conda方法（推荐）：
  conda_create("r-reticulate")
  conda_install(envname = "r-reticulate", packages="louvain")
  use_python("~/miniconda3/envs/r-reticulate/lib/python3.7/site-packages/")
  # RETICULATE_PYTHON="~/miniconda3/envs/r-reticulate/bin/python3"
  py_config()
}


# 现在安装3版本，还是要通过github
devtools::install_github('cole-trapnell-lab/monocle3')
suppressMessages(library(monocle3))

#############################
# 创建对象
#############################
rm(list = ls())  
options(warn=-1)  
# 下载数据(如果已经有了，可以跳过)
expression_matrix <- readRDS(file = 'cao_l2_expression.rds')
cell_metadata <- readRDS(file = 'cao_l2_colData.rds')
gene_annotation <- readRDS(file = 'cao_l2_rowData.rds')

dim(expression_matrix)
expression_matrix[1:3,1:3]
cell_metadata[1:3,1:3]
head(gene_annotation)

# 创建CDS(Cell Data Set)对象
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds

#############################
# 数据预处理（归一化、标准化、降维）
#############################

start=Sys.time()
cds <- preprocess_cds(cds, num_dim = 100)
end=Sys.time()
(end-start)
# 检查一下使用的主成分（PCs）是否能够抓取最主要的基因表达变化信息
plot_pc_variance_explained(cds)

#############################
# 继续在preprocess_cds（PCA）的基础上降维UMAP/tSNE
#############################
start=Sys.time()
cds <- reduce_dimension(cds, preprocess_method = "PCA",reduction_method = c("UMAP"))
end=Sys.time()
(end-start)

plot_cells(cds)
# 根据细胞类型上色
plot_cells(cds, color_cells_by="cao_cell_type")
# 根据不同基因表达上色
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))

# 换种降维方法tSNE
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type")

#############################
# 批次效应
#############################
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)
# 假入要对批次进行校正
cds = preprocess_cds(cds, num_dim = 100, residual_model_formula_str = "~ plate")
# 校正后再降维
cds = reduce_dimension(cds)
# 再次检查批次效应
plot_cells(cds, color_cells_by="plate", label_cell_groups=FALSE)

#############################
# 细胞聚类
#############################
cds = cluster_cells(cds, resolution=c(10^seq(-6,-1)))
plot_cells(cds)
# 将小的cluster聚合成一个个大的partition
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")
# 用细胞类型对cluster进行编号
plot_cells(cds, color_cells_by="cao_cell_type")
# 简化label
plot_cells(cds, color_cells_by="cao_cell_type", label_groups_by_cluster=FALSE)

save(cds, file="L2_clustered_cds.Rdata")
#############################
# 找Marker基因
#############################
marker_test_res = top_markers(cds, group_cells_by="partition", reference_cells=1000, cores=8)
marker_test_res[1:4,1:4]
# 过滤得到unique gene id
library(tidyverse)
top_specific_markers = marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)
top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
# 对每组的marker基因可视化
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="maximal_on_diag",
                    max.size=3)

# 如果要在每组多看几个marker，只需要修改一下top_n
if(T){
  top_specific_markers = marker_test_res %>%
    filter(fraction_expressing >= 0.10) %>%
    group_by(cell_group) %>%
    top_n(3, pseudo_R2)
  
  top_specific_marker_ids = unique(top_specific_markers %>% pull(gene_id))
  
  plot_genes_by_group(cds,
                      top_specific_marker_ids,
                      group_cells_by="partition",
                      ordering_type="cluster_row_col",
                      max.size=3)
}

#############################
# 对细胞类型进行注释
#############################
if(F){
  # 先将partitions的分组由因子型转为字符型
  colData(cds)$assigned_cell_type = as.character(partitions(cds))
  # 再对字符型重新定义
  colData(cds)$assigned_cell_type = dplyr::recode(colData(cds)$assigned_cell_type,
                                                  "1"="Body wall muscle",
                                                  "2"="Germline",
                                                  "3"="Unclassified neurons",
                                                  "4"="Seam cells",
                                                  "5"="Coelomocytes",
                                                  "6"="Pharyngeal epithelia",
                                                  "7"="Vulval precursors",
                                                  "8"="Non-seam hypodermis",
                                                  "9"="Intestinal/rectal muscle",
                                                  "10"="Touch receptor neurons",
                                                  "11"="Pharyngeal neurons",
                                                  "12"="Am/PH sheath cells",
                                                  "13"="NA",
                                                  "14"="Unclassified neurons",
                                                  "15"="flp-1(+) interneurons",
                                                  "16"="Canal associated neurons",
                                                  "17"="Pharyngeal gland",
                                                  "18"="Other interneurons",
                                                  "19"="Ciliated sensory neurons",
                                                  "20"="Ciliated sensory neurons",
                                                  "21"="Ciliated sensory neurons",
                                                  "22"="Ciliated sensory neurons",
                                                  "23"="Ciliated sensory neurons",
                                                  "24"="Ciliated sensory neurons",
                                                  "25"="Oxygen sensory neurons",
                                                  "26"="Ciliated sensory neurons",
                                                  "27"="Unclassified neurons",
                                                  "28"="Pharyngeal gland",
                                                  "29"="Ciliated sensory neurons",
                                                  "30"="Ciliated sensory neurons",
                                                  "31"="Ciliated sensory neurons",
                                                  "32"="Ciliated sensory neurons",
                                                  "33"="Pharyngeal muscle",
                                                  "34"="Failed QC")
  # 另外，想从中取子集、过滤的话
  cds[,colData(cds)$assigned_cell_type != "Failed QC"]
}
plot_cells(cds, group_cells_by="partition", color_cells_by="assigned_cell_type")

#############################
# 自动化注释--Garnett
#############################
## 预处理
# step-1：首先根据上面得到的细胞类型assigned_cell_type，找top_marker
assigned_type_marker_test_res = top_markers(cds,
                                            group_cells_by="assigned_cell_type",
                                            reference_cells=1000,
                                            cores=8)
# step-2：过滤（阈值自定义）
garnett_markers = assigned_type_marker_test_res %>%
  filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
  group_by(cell_group) %>%
  top_n(5, marker_score)
# step-3：去重复
garnett_markers = garnett_markers %>% group_by(gene_short_name) %>%
  filter(n() == 1)
# step-4：生成marker文件
generate_garnett_marker_file(garnett_markers, file="./marker_file.txt")

## 配置
## Install the monocle3 branch of garnett
devtools::install_github("cole-trapnell-lab/garnett", ref="monocle3")
library(garnett)
# install gene database for worm
BiocManager::install("org.Ce.eg.db")
## 自己构建
if(F){
  colData(cds)$garnett_cluster = clusters(cds)
  worm_classifier <- train_cell_classifier(cds = cds,
                                           marker_file = "./marker_file.txt",
                                           db=org.Ce.eg.db::org.Ce.eg.db,
                                           cds_gene_id_type = "ENSEMBL",
                                           num_unknown = 50,
                                           marker_file_gene_id_type = "SYMBOL",
                                           cores=8)
  cds = classify_cells(cds, worm_classifier,
                       db = org.Ce.eg.db::org.Ce.eg.db,
                       cluster_extend = TRUE,
                       cds_gene_id_type = "ENSEMBL")
  plot_cells(cds,
             group_cells_by="partition",
             color_cells_by="cluster_ext_type")
}
## 使用别人的（https://cole-trapnell-lab.github.io/garnett/classifiers/ceWhole）
load('ceWhole')
cds = classify_cells(cds, ceWhole,
                     db = org.Ce.eg.db::org.Ce.eg.db,
                     cluster_extend = TRUE,
                     cds_gene_id_type = "ENSEMBL")


#############################
# 构建发育轨迹
#############################
rm(list = ls())  
options(warn=-1)  
# step-0:载入数据,创建对象
expression_matrix <- readRDS(file = './monocle-2-trajectories/packer_embryo_expression.rds')
cell_metadata <- readRDS(file = './monocle-2-trajectories/packer_embryo_colData.rds')
gene_annotation <- readRDS(file = './monocle-2-trajectories/packer_embryo_rowData.rds')
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
# step-1:预处理
cds <- preprocess_cds(cds, num_dim = 100, 
                      residual_model_formula_str = "~ bg.300.loading + bg.400.loading + bg.500.1.loading 
                      + bg.500.2.loading + bg.r17.loading + bg.b01.loading + bg.b02.loading")
# step-2:降维（UMAP）+可视化
cds <- reduce_dimension(cds)
plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell.type")
# 看不同基因在不同分支的表达量
if(F){
  ciliated_genes = c("che-1",
                     "hlh-17",
                     "nhr-6",
                     "dmd-6",
                     "ceh-36",
                     "ham-1")
  
  plot_cells(cds,
             genes=ciliated_genes,
             label_cell_groups=FALSE,
             show_trajectory_graph=FALSE)
}
# step-3:聚类(得到的每个partition都是一个轨迹)
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")
# step-4: learn graph
cds <- learn_graph(cds)
plot_cells(cds,
           color_cells_by = "cell.type",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE)
# step-5: 细胞排序
# 手动找根节点
cds = order_cells(cds)
plot_cells(cds,
           color_cells_by = "embryo.time.bin",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5)
# 自动化找根节点（首先根据轨迹图中的节点将与它们最相近的细胞分成组，
# 然后计算来自最早时间点的细胞所占组分，最后挑出包含早期细胞数量最多的节点，认为它是就是根节点）
if(F){
  get_earliest_principal_node <- function(cds, time_bin="130-170"){
    cell_ids <- which(colData(cds)[, "embryo.time.bin"] == time_bin)
    
    closest_vertex <-
      cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
      igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
  }
  cds = order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
}
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


# 3D trajectories
cds_3d = reduce_dimension(cds, max_components = 3)
cds_3d = cluster_cells(cds_3d)
cds_3d = learn_graph(cds_3d)
cds_3d = order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))

cds_3d_plot_obj = plot_cells_3d(cds_3d, color_cells_by="partition")


#############################
# 差异分析之Regression analysis
#############################
# 这里构建一个小数据集（实际上可以分析上千差异基因）
ciliated_genes = c("che-1",
                   "hlh-17",
                   "nhr-6",
                   "dmd-6",
                   "ceh-36",
                   "ham-1")
cds_subset = cds[rowData(cds)$gene_short_name %in% ciliated_genes,]

# step-0:
# 其中model_formula_str就是要比较的分组对象，如果相获得不同的cluster或者partition的差异基因，
# 就用model_formula_str = "~cluster"或者model_formula_str = "~partition"；
# 另外还支持添加多个变量，比如考虑到批次效应 model_formula_str = "~embryo.time + batch"
gene_fits = fit_models(cds_subset, model_formula_str = "~embryo.time")

# step-1:
fit_coefs = coefficient_table(gene_fits)
# 挑出时间相关的组分
emb_time_terms = fit_coefs %>% filter(term == "embryo.time")
# coefficient_table()默认使用 Benjamini and Hochberg（BH）方法进行了p值的校正，得到了q值
emb_time_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)

# （可选）小提琴图
plot_genes_violin(cds_subset, group_cells_by="embryo.time.bin", ncol=2) +
  theme(axis.text.x=element_text(angle=45, hjust=1))

# （考虑批次效应）
gene_fits = fit_models(cds_subset, model_formula_str = "~embryo.time + batch")
fit_coefs = coefficient_table(gene_fits)
fit_coefs %>% filter(term != "(Intercept)") %>%
  select(gene_short_name, term, q_value, estimate)

# step-2: 评价模型
evaluate_fits(gene_fits)
# 模型中到底要不要考虑批次呢？使用compare_models()比较判断：
time_batch_models = fit_models(cds_subset,
                               model_formula_str = "~embryo.time + batch",
                               expression_family="negbinomial")
time_models = fit_models(cds_subset,
                         model_formula_str = "~embryo.time",
                         expression_family="negbinomial")
compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)

#############################
# 差异分析之Graph-autocorrelation analysis
#############################
# 挑一部分细胞
neurons_cds = cds[,colData(cds)$assigned_cell_type == "Neurons"]
plot_cells(neurons_cds, color_cells_by="partition")
# 这里才开始分析（graph_test 得到莫兰指数）
pr_graph_test_res = graph_test(neurons_cds, neighbor_graph="knn", cores=8)
pr_deg_ids = row.names(subset(pr_graph_test_res, q_value < 0.05))

#############################
# # 将共同作用的基因合并成模块
#############################
gene_module_df = find_gene_modules(neurons_cds[pr_deg_ids,], resolution=1e-2)
# 第一种可视化
cell_group_df = tibble::tibble(cell=row.names(colData(neurons_cds)), cell_group=partitions(cds)[colnames(neurons_cds)])
agg_mat = aggregate_gene_expression(neurons_cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) = stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
# 第二种可视化
plot_cells(neurons_cds,
           genes=gene_module_df %>% filter(module %in% c(16,38,33,42)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)

#############################
# 找到影响发育轨迹的基因
#############################
ciliated_cds_pr_test_res = graph_test(cds, neighbor_graph="principal_graph", cores=4)
# 使用neighbor_graph="principal_graph"来检验轨迹相邻的细胞的表达是否相关
pr_deg_ids = row.names(subset(ciliated_cds_pr_test_res, q_value < 0.05))
# 从pr_deg_ids中挑选几个基因，然后可视化
plot_cells(cds, genes=c("hlh-4", "gcy-8", "dac-1", "oig-8"),
           show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE)

# 可以将这些在轨迹上变化的pr_deg_ids继续分为小模块
gene_module_df = monocle3:::find_gene_modules(cds[pr_deg_ids,], resolution=c(0,10^seq(-6,-1)))
# 
cell_group_df = tibble::tibble(cell=row.names(colData(cds)), cell_group=colData(cds)$cell.type)
agg_mat = aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) = stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(29,20, 11,22)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

## 另外还有一种方法：
# 先画完整的轨迹图：
plot_cells(cds, show_trajectory_graph=FALSE)
# 然后挑某一个细胞亚群对应的一段轨迹：
AFD_genes = c("gcy-8", "dac-1", "oig-8")
AFD_lineage_cds = cds[rowData(cds)$gene_short_name %in% AFD_genes,
                      clusters(cds) %in% c(22, 28, 35)]

plot_genes_in_pseudotime(AFD_lineage_cds,
                         color_cells_by="embryo.time.bin",
                         min_expr=0.5)


