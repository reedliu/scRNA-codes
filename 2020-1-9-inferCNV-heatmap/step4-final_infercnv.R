### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2020-01-09
### Email: jieandze1314@gmail.com
### Title: 改造inferCNV的热图
### ---------------

# 目的：重复final_infercnv的图：https://jieandze1314-1255603621.cos.ap-guangzhou.myqcloud.com/blog/2020-01-09-014050.png

rm(list = ls())
options(stringsAsFactors = F)

##############################
# step0：数据准备
##############################
infercnv_obj <- readRDS('run.final.infercnv_obj')
# 提取：处理后的表达矩阵
expr <- infercnv_obj@expr.data
dim(expr)
expr[1:4,1:4]

# 提取：表达矩阵样本中正常细胞的位置
(normal_all_loc <- infercnv_obj@reference_grouped_cell_indices)
normal_loc <- normal_all_loc$WT
length(normal_loc)
# 提取：表达矩阵样本中肿瘤细胞的位置
# final数据中的tumor又分成了三类：tumor、mut、control组
tumor_all_loc <- infercnv_obj@observation_grouped_cell_indices

tumor_loc <- tumor_all_loc$Tumor
tumor_mut_loc <- tumor_all_loc$MUT
tumor_ctrl_loc <- tumor_all_loc$Control

length(c(tumor_loc,tumor_mut_loc,tumor_ctrl_loc))

##############################
## Step1: 将基因名与染色体位置对应
##############################
gn <- rownames(expr)
length(gn)
head(gn)
# 加载基因位置信息文件(这个存储的是排序后的)
geneFile <- read.table('geneFile.txt')
geneFile[1:4,1:4]
length(geneFile$V1);length(gn) #下面需要根据gn对geneFile取小

sub_geneFile <-  geneFile[geneFile$V1%in%gn,]
dim(sub_geneFile)
length(sub_geneFile$V1);length(gn) #发现取小后sub_geneFile有500个基因不在gn中

# 那就先按sub_geneFile中包含的基因来做后续分析
sub_expr <- expr[sub_geneFile$V1,]
sub_expr <- as.data.frame(sub_expr)
identical(rownames(sub_expr),sub_geneFile$V1)

##############################
## Step2: 拆分矩阵
##############################
# 整体分成两部分：上面的热图是正常细胞，下面是肿瘤细胞，并且我们知道了各自的位置，就能先把各自的小表达矩阵提取出来

norm_expr <- sub_expr[,normal_loc]
norm_expr$chr <- as.factor(sub_geneFile$V2) #最后加上一列：对应的chr信息
table(norm_expr$chr) #这个信息就与横坐标的间隔对应，chr间隔越大表示其中包含的基因越多

# 原图中是从上而下按照：Control、MUT、Tumor顺序来的，这里也按这个顺序取小
tumor_expr <- sub_expr[,c(tumor_ctrl_loc,tumor_mut_loc,tumor_loc)]
tumor_expr$chr <- as.factor(sub_geneFile$V2)
dim(norm_expr);dim(tumor_expr) # 原来580个样本，现在normal、tumor各自又增加一列chr

##############################
## Step3-4: 画图-Tumor下半部分
##############################
library(ComplexHeatmap)
# 1 首先做一个样本的注释信息，一个注释是图中最左侧注释条：全是Tumor；
#  另一个是右侧的注释条：从上到下分成Control、MUT、Tumor

tumor_all_name <- colnames(sub_expr)[c(tumor_ctrl_loc,tumor_mut_loc,tumor_loc)]
meta <- data.frame(sample=tumor_all_name,
                   all = "Tumor",
                    type=c(rep('Control',length(tumor_ctrl_loc)),
                          rep('MUT',length(tumor_mut_loc)),
                          rep('Tumor',length(tumor_loc))))
head(meta)

## -----表达矩阵部分-----
tmp <- tumor_expr[,-ncol(tumor_expr)]
tmp <- tmp[,tumor_all_name]

n <- t(tmp)
dim(n)
identical(rownames(n),meta$sample)

n[n>1.6]=1.6
n[n<0.4]=0.4
n[0.8<n & n<1.2]=1

## ------设置行注释信息【sample】--------
# https://www.biostars.org/p/317349/
ann <- data.frame(all=meta$all,type=meta$type)
colours <- list("all"=c("Tumor"="#8FD3C7"),
                "type"=c("Control"="#FCEC78",
                         "MUT"="#D8C868",
                         "Tumor"="#8FD3C7"))

rowAnn <- HeatmapAnnotation(df=ann, which="row", col=colours, 
                            annotation_width=unit(c(1, 4), "cm"), 
                            gap=unit(1, "mm"))

##---- 设置列的注释信息【chr】------
color <- c('#8FD3C7','#CFEBBD','#F4F3BC','#CFCCCF','#D0A789','#F2877E','#BF989F',
           '#87B2CC','#CDB18D','#ECBB6A','#C3D46E','#F7CDDE','#E9D3DE','#D5D0D6',
            '#C49DC4','#BE9DBE','#C9DAC3','#E1EAA4','#FCEC72')
# 设定不同chr出现的顺序
new_cluster <- tumor_expr$chr
names(color) <- levels(new_cluster)

top_color <- HeatmapAnnotation(
  cluster = anno_block(gp = gpar(fill = color), # 设置填充色
                       labels = levels(new_cluster), 
                       labels_gp = gpar(cex = 0.8, col = "black"))) 


# 如果要指定图例的颜色（说明最小值0.4对应blue；中间1对应white）
library("circlize")
col_fun = colorRamp2(c(0.4, 1, 1.6), c("blue", "white", "red"))

##------ 如何在热图上画线----------
# 参考：https://jokergoo.github.io/ComplexHeatmap-reference/book/heatmap-decoration.html
# 下面👇为示例
# set.seed(123)
# Heatmap(matrix(rnorm(100), 10), name = "mat")
# decorate_heatmap_body("mat", {
#   #在x轴上画线
#   grid.lines(c(0.5, 0.5), c(0, 1), gp = gpar(lty = 2, lwd = 2))
#   #在y轴上画线
#   grid.lines(c(0, 1), c(0.1, 0.1), gp = gpar(lty = 2, lwd = 2))
# })


##------ 最后开始作图 ----------
pdf("final-infercnv-tumor-heatmap-10.pdf",width = 15,height = 10)
if(T){
  ht_tumor = Heatmap(as.matrix(n),
                     name="ht_tumor",
                     col = col_fun,
                     cluster_rows = F,
                     #clustering_method_rows = 'ward.D2',
                     cluster_columns = F,
                     show_column_names = F,
                     show_row_names = F,
                     column_split = new_cluster,
                     column_gap = unit(0, "mm"),
                     heatmap_legend_param = list(
                       title = "Modified Expression",
                       title_position = "leftcenter-rot", # 图例标题位置
                       at=c(0.4,0.8,1,1.2,1.6), #图例范围
                       legend_height = unit(3, "cm") #图例长度
                     ),
                     top_annotation = top_color,
                     row_title = "Observations (Cells)",
                     row_title_side = c("right"),
                     row_title_rot = 90,
                     column_title = "Genomic Region",
                     column_title_side = c("bottom"),
                     left_annotation=rowAnn)
  draw(ht_tumor, heatmap_legend_side = "left") # 图例位置
  
}

## -----如何画出【x轴】各个chr的区分线？-------
# 要画x轴上的分割线，就要得到每个区间右侧的坐标
# 因为之前把各个chr进行了column_split，所以画图会认为每个分开的chr为一个单独的区间
# decorate_heatmap_body会以第一个chr1为基准，它的范围是[0,1]。要得到之后各个chr的右侧坐标，就要先累加，然后除以第一个chr1的长度
# 而且前提是每个单独的区间之间不存在间隔，因此前面把column_gap设成了0
#
if(T){
  chr_info <- as.numeric(table(sub_geneFile$V2))
  x=cumsum(chr_info)/chr_info[1]
  # 然后自定义一个函数，传递给sapply，给每个chr画图
  draw_x_line=function(x,ht,1color){
    decorate_heatmap_body(ht1, {
      grid.lines(c(x, x), c(0, 1), gp = gpar(lty = 2, lwd = 2, col=color))
    })
  }
  sapply(x, function(x)draw_x_line(x,ht1="ht_tumor",color="black"))
}

## -----如何画出【y轴】各个sample type的区分线？-------
if(T){
  # 这个就比较简单了，y轴和x轴不同，它是作为一整个区间，所以不需要像x轴一样去累加，从上到下就是[0,1]
  type_info=as.numeric(table(ann$type))
  # 又因为y轴是从下往上画，所以顺序是tumor=》mut=》control。只需要得到两条线的位置即可（也就是下面的y）
  type_info=rev(type_info)
  y_dat=cumsum(type_info/sum(type_info))[1:2]
  y_dat
  # 另一个坑就是：grid.lines中的第一个指定x不能是c(0, 1)，这样只会画出chr1中一条短短的线，要得到所有chr的总长度
  # 也就是存储在之前cumsum计算的x的最后一个值
  draw_y_line=function(y,ht2,color){
    decorate_heatmap_body(ht2, {
      grid.lines(c(0, x[length(x)]), c(y, y), gp = gpar(lty = 2, lwd = 2, col=color))
    })
  }
  sapply(y_dat, function(y)draw_y_line(y,ht2="ht_tumor",color="black"))
}

dev.off()


##############################
## Step3-5: 画图-整合
##############################
## -----整合画图-------
m <- t(norm_expr[,-ncol(norm_expr)])
dim(m)
ann2 <- data.frame(type=rep("WT",nrow(m)))
colours2 <- list("type"=c("WT"="#8FD3C7"))

colAnn2 <- HeatmapAnnotation(df=ann2, which="row", col=colours2, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))

m[m>1.6]=1.6
m[m<0.4]=0.4
m[0.8<m & m<1.2]=1


if(T){
  
  ht_normal = Heatmap(as.matrix(m),
                      name = "ht_normal",
                      col = col_fun,
                      cluster_rows = F,
                      cluster_columns = F,
                      show_column_names = F,
                      show_row_names = F,
                      column_split = new_cluster,
                      column_gap = unit(0, "mm"),
                      row_title = "References (Cells)",
                      row_title_side = c("right"),
                      row_title_rot = 90,
                      # row_title_gp = gpar(fontsize = 25),
                      column_title = NULL, 
                      heatmap_legend_param = list(
                        title = "Modified Expression",
                        title_position = "leftcenter-rot", # 图例标题位置
                        # title_gp = gpar(fontsize = 20),# 图例标题大小
                        at=c(0.4,0.8,1,1.2,1.6), #图例范围
                        legend_height = unit(3, "cm")),#图例长度
                      left_annotation=colAnn2,
                      width = 20, height = 5) 
  
  ht_tumor = Heatmap(as.matrix(n),
                     name="ht_tumor",
                     col = col_fun,
                     cluster_rows = F,
                     #clustering_method_rows = 'ward.D2',
                     cluster_columns = F,
                     show_column_names = F,
                     show_row_names = F,
                     column_split = new_cluster,
                     column_gap = unit(0, "mm"),
                     show_heatmap_legend=F,
                     top_annotation = top_color,
                     row_title = "Observations (Cells)",
                     row_title_side = c("right"),
                     row_title_rot = 90,
                     column_title = "Genomic Region",
                     column_title_side = c("bottom"),
                     left_annotation=rowAnn,
                     width = 20, height = 10,
                     heatmap_height = 15)

  
}
pdf("final-infercnv-tumor-heatmap-11.pdf",width = 15,height = 10)
draw(ht_normal %v% ht_tumor)
# 不要最后一条线
dat=x[-length(x)]
sapply(dat, function(x)draw_x_line(x,
                                 ht1="ht_normal",color="bluack))
sapply(dat, function(x)draw_x_line(x,
                                 ht1="ht_tumor",color="bluack))

sapply(y, function(y)draw_y_line(y,
                                 ht2="ht_tumor",color="black"))

dev.off()







