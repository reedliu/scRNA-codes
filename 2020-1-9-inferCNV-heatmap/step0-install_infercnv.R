### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2020-01-03
### Email: jieandze1314@gmail.com
### Title: 改造inferCNV的热图
### ---------------

# 安装说明：https://github.com/broadinstitute/inferCNV/wiki/Installing-infercnv

#--- first install the JAGS package ------
# for mac: brew install jags
# for windows: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/JAGS-4.3.0.exe
# for linux: conda install

#--- second install the infercnv package ------
BiocManager::install("infercnv")
