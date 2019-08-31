#############################################################
############# Single Cell RNA-Seq for MCC ###################
# Author: Reed Liu 
# Mail: jieandze1314@gmail.com
# version-1: 2019-5-14  Firstly Build
# version-2: 2019-7-17  Add aggr from 2-1
#############################################################


###############################################
## 1-1: mkfastq (from bcl to fastq) 
###############################################
# 如果自己有了fq文件，就不需要这一步
# 要是从BCL初始测序文件开始，先新建一个目录
wkd=/home/single-cell/MCC
mkdir $wkd/bcl && cd "$_"

# 标准模式（将bcl文件放在bcl目录下）
cellranger mkfastq   --id=new-bcl \
                     --run=$wkd/bcl \
                     --csv=cellranger-tiny-bcl-simple-1.2.0.csv

# 当从一个大的flowcell中拆分少量样本时，可以去掉一些无关的fastq(undetermined)
cellranger mkfastq   --id=new2-bcl \
                     --run=$wkd/bcl \
                     --csv=cellranger-tiny-bcl-simple-1.2.0.csv \
                     --delete-undetermined
# --id这里定义输出的目录名，默认按照--run的名称 

# 可能的报错：
# 1.[error] No bcl2fastq found on path. demux requires bcl2fastq v2.17 or greater for RTA version: 1.18.66.4
# 原因：cellranger依赖bcl2fastq，但是没有安装
# 解决：conda 下载一个：conda install -c dranew bcl2fastq

# 输出结果在$wkd/bcl/outs/fastq_path/H35KCBCXY/test_sample


###############################################
## 2-1: count 定量
###############################################
# 小练习
# make file
for i in $(seq 1 4);do echo SRR811169$i;done >id.txt
# download sra(每个文件3G左右)
cat id.txt | while read i;do prefetch $i -O ./; echo "** ${i}.sra done **";done
# sra2fq
cat id.txt | while read i;do \
time fastq-dump --gzip --split-files -A $i ${i}.sra && echo "** ${i}.sra to fastq done **";\
done


# 文章数据
cd $wkd/count
cat $wkd/raw/P2586-4/SRR_Acc_List-2586-4.txt| \
while read i;do \
(cellranger count --id=$i \
                   --transcriptome=$wkd/ref/hg38_and_mcv \
                   --fastqs=$wkd/raw \
                   --sample=$i \
                   --nosecondary \
                   --localcores=44 \
                   --localmem=200); \
done

# 其中nosecondary表示只获得矩阵，不进行后续的降维、聚类和可视化分析

# 使用localcores之前，应该先检查ulimit -u，看看服务器最大支持多少用户同时在线，因为cellranger使用一个核就会产生64个用户队列，不能超过这个限定值。例如检查ulimit -u为4096，那么最多设置64个核心，也就是64 * 64 = 4096 ，才不会因为这个报错

###############################################
## 2-2: aggr 整合
###############################################
# 当分析多个生物学样本或者一个样本的多个文库/技术重复时，最好先单独跑count，再aggr合并。并且需要将tumor、PBMC分别整合
# 首先创建一个Aggregation CSV：需要指定两列并包含表头，library_id、molecule_h5
# （示例）
# library_id,molecule_h5
# LV123,/opt/runs/LV123/outs/molecule_info.h5
# LB456,/opt/runs/LB456/outs/molecule_info.h5
# LP789,/opt/runs/LP789/outs/molecule_info.h5

# 然后跑流程
cd $wkd/count
cellranger aggr --id=2586-4-tumor-2 \
                  --csv=2586-4-tumor_libraries.csv \
                  --normalize=raw

# 注意这里有个参数：mapped。这个是默认参数，当组合不同的GEMs时，考虑到测序深度不同导致的批次效应，

# 另外，如果要整合不同试剂版本的数据，还需要增加一列——batch（cellranger 3.0版本才能使用）
# library_id,molecule_h5,batch
# LV123,/opt/runs/LV123/outs/molecule_info.h5,v2_lib
# LB456,/opt/runs/LB456/outs/molecule_info.h5,v3_lib
# LP789,/opt/runs/LP789/outs/molecule_info.h5,v3_lib





















