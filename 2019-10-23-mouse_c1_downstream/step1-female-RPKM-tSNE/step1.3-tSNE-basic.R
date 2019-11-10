### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-23
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-更简单的basic函数-tSNE分群
### ---------------
rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
load('../female_rpkm.Rdata')

females[1:4,1:4]
dim(females)

## 下面代码参考之前的：19-10-15-GSE117988-10X免疫：step5-R-pkgs-or-not.R
## 拿到分群结果
cluster <- read.csv('step1.1-D-female_clustering.csv')
table(cluster[,2])
color <- rainbow(4)[as.factor(cluster[,2])]
table(color)

if(T){
  choosed_count <- females
  # 表达矩阵过滤
  choosed_count <- choosed_count[apply(choosed_count, 1, sd)>0,]
  choosed_count <- choosed_count[names(head(sort(apply(choosed_count, 1, sd),decreasing = T),1000)),]
  # 然后进行PCA分析
  pca_out <- prcomp(t(choosed_count),scale. = T)
  library(ggfortify)
  autoplot(pca_out, col=color) +theme_classic()+ggtitle('PCA plot')
  
  str(pca_out)
  pca_out$x[1:3,1:3]
  library(Rtsne)
  tsne_out <- Rtsne(pca_out$x[,1:9], perplexity = 10,
                    pca = F, max_iter = 2000,
                    verbose = T)
  tsnes_cord <- tsne_out$Y
  colnames(tsnes_cord) <- c('tSNE1','tSNE2')
  ggplot(tsnes_cord, aes(x=tSNE1, y = tSNE2)) + geom_point(col=color) + theme_classic()+ggtitle('tSNE plot')
}

# 除了使用HCPC分群，还可以用DBSCAN、kmeans
if(T){
  # 这个运行会非常慢！
  library(Rtsne)
  N_tsne <- 50
  tsne_out <- list(length = N_tsne)
  KL <- vector(length = N_tsne)
  set.seed(1234)
  for(k in 1:N_tsne)
  {
    tsne_out[[k]]<-Rtsne(t(log2(females+1)),initial_dims=30,verbose=FALSE,check_duplicates=FALSE,
                         perplexity=27, dims=2,max_iter=5000)
    KL[k]<-tail(tsne_out[[k]]$itercosts,1)
    print(paste0("FINISHED ",k," TSNE ITERATION"))
  }
  names(KL) <- c(1:N_tsne)
  opt_tsne <- tsne_out[[as.numeric(names(KL)[KL==min(KL)])]]$Y
}
# DBSCAN
library(dbscan)
plot(opt_tsne,  col=dbscan(opt_tsne,eps=3.5)$cluster, 
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")
# kmeans
plot(opt_tsne,  col=kmeans(opt_tsne,centers = 4)$clust, 
     pch=19, xlab="tSNE dim 1", ylab="tSNE dim 2")

table(kmeans(opt_tsne,centers = 4)$clust,dbscan(opt_tsne,eps=3.5)$cluster)

# 还可以和之前的cluster结果比较
# HCPC结果
cluster1 <- read.csv('step1.1-D-female_clustering.csv')
# seurat结果
load('step1.2-B-seurat3-female-tsne.Rdata')
cluster2 <- as.data.frame(Idents(sce_female_tsne))


