### ---------------
### Creator: Yunze Liu (Reed Liu)
### Date: 2019-10-24
### Email: jieandze1314@gmail.com
### Title: Smartseq2-C1-小鼠性腺发育-差异分析-比较
### ---------------
rm(list=ls())
options(stringsAsFactors = F)

#############################
# 准备两种差异分析结果
#############################
load(file = 'step3.1-DEG-monocle_summary.Rdata')
load(file = 'step3.2-DEG-ROTS_summary.Rdata')

# monocle差异基因
mnc_DEG <- subset(de_clusters, qval<0.05)

# ROTS差异基因
head(summary_pop1)
head(summary_pop2)
head(summary_pop3)
head(summary_pop4)

# 对第一群比较
g1=rownames(summary_pop1[summary_pop1$FDR>0.05,])
g2=rownames(mnc_DEG[mnc_DEG$cluster=='C1',])
length(intersect(g1,g2));length(g1);length(g2)

# 对第二群比较
g1=rownames(summary_pop2[summary_pop2$FDR>0.05,])
g2=rownames(mnc_DEG[mnc_DEG$cluster=='C2',])
length(intersect(g1,g2));length(g1);length(g2)

# 对第三群比较
g1=rownames(summary_pop3[summary_pop3$FDR>0.05,])
g2=rownames(mnc_DEG[mnc_DEG$cluster=='C3',])
length(intersect(g1,g2));length(g1);length(g2)

## 对第四群比较
g1=rownames(summary_pop2[summary_pop2$FDR>0.05,])
g2=rownames(mnc_DEG[mnc_DEG$cluster=='C2',])
length(intersect(g1,g2));length(g1);length(g2)





