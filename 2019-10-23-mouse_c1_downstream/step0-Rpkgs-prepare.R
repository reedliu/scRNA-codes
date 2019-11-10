### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-11-10
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-R包准备工作（补充）
### ---------------

##############################
## install cran pkgs
##############################


if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
cran_packages <- c('taRifx',
                   'matrixStats',
                   'ggplot2',
                   'Rtsne',
                   'fpc',
                   'factoextra',
                   'viridis',
                   'gplots',
                   'RColorBrewer',
                   'rgl',
                   'scatterplot3d',
                   'pheatmap',
                   'matrixStats',
                   'statmod',
                   'FactoMineR',
                   'jackstraw',
                   'arulesViz',
                   'ggpubr'
                   )
for (pkg in cran_packages){
  if (! require(pkg,character.only=T) ) {
    install.packages(pkg,ask = F,update = F)
    suppressMessages(require(pkg,character.only=T)) 
  }
}

##############################
## install bioconductor pkgs
##############################

# first prepare BioManager on CRAN
if(length(getOption("CRAN"))==0) options(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/")

if(!require("BiocManager")) install.packages("BiocManager",update = F,ask = F)

if(length(getOption("BioC_mirror"))==0) options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")


# use BiocManager to install
Biocductor_packages <- c('monocle',
                         'destiny',
                         'slingshot',
                         'made4',
                         'ReactomePA',
                         'org.Mm.eg.db',
                         'clusterProfiler',
                         'GOSemSim',
                         'lfa')
for (pkg in Biocductor_packages){
  if (! require(pkg,character.only=T) ) {
    BiocManager::install(pkg,ask = F,update = F)
    suppressMessages(require(pkg,character.only=T)) 
  }
}





