### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-谱系推断-Slingshot
### ---------------

rm(list = ls()) 
options(warn=-1) 
options(stringsAsFactors = F)
source("../analysis_functions.R")
load('../female_rpkm.Rdata')
dim(females)

#############################
# 运行slingshot
#############################
load('step4.1-diffusionMap_output.Rdata')
# 使用包装的函数
female_lineage <- get_lineage(
  dm=female_dm, 
  dim=c(1:4), 
  condition=factor(female_clustering),
  start="C1",
  end=c("C2", "C4"),
  shrink.method="cosine"
)

# 如果不使用包装的函数
if(F){
  dm=female_dm
  dim=c(1:4)
  condition=factor(female_clustering)
  data <- data.frame(
    dm@eigenvectors[,dim]
  )
  crv <- slingshot(
    data, 
    condition, 
    start.clus = "C1", 
    end.clus=c("C2", "C4"),
    maxit=100000,
    shrink.method="cosine"
    # shrink.method="tricube"
  )
}


#############################
# 画图
#############################
##  如果要画3D图
# 探索6个发育时间
female_stagePalette=c(
  E10.5="#2754b5", 
  E11.5="#8a00b0", 
  E12.5="#d20e0f", 
  E13.5="#f77f05", 
  E16.5="#f9db21",
  P6="#43f14b"
)
plot_dm_3D(
  dm=female_dm, 
  dc=c(1:3),
  condition=female_stages, 
  colour=female_stagePalette
)

#############################
# 运行结果
#############################
female_pseudotime <- get_pseudotime(female_lineage, wthres=0.9)
rownames(female_pseudotime) <- colnames(females)
save(female_pseudotime,file = 'step4.2-female_pseudotime.Rdata')
write.csv(female_pseudotime, file="step4.2-female_pseudotime.csv")




