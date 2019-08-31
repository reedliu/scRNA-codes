#############################################################
############# Single Cell RNA-Seq for MCC ###################
# Author: Reed Liu 
# Mail: jieandze1314@gmail.com
# version-1: 2019-5-3  Firstly Build
# version-2: 2019-5-15 Rearrange the pipeline 
#############################################################

###############################################
## 1-1: Download SRA raw data
###############################################
wkd=/home/single-cell/MCC
mkdir -p $wkd/{raw,ref,qc,count,biosoft}

cd $wkd/raw
mkdir P2586-4 P9245-3

# for patient P2586-4
# https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA483959
# 37(Pre)、38(AR) Tumor组；39-42PBMC组(pre,early,resp,ar)
cd P2586-4
cat >SRR_Acc_List-2586-4.txt
SRR7722937 
SRR7722938
SRR7722939
SRR7722940
SRR7722941
SRR7722942

cat SRR_Acc_List-2586-4.txt |while read i
do prefetch $i -O `pwd` && echo "** ${i}.sra done **"
done

# for patient P9245-3
# 88、89 Tumor组；86、87 PBMC组
cd P9245-3
cat >SRR_Acc_List-9245-3.txt
SRR7692286
SRR7692287
SRR7692288
SRR7692289


cat SRR_Acc_List-9245-3.txt |while read i
do prefetch $i -O `pwd` && echo "** ${i}.sra done **"
done



# 附：对于NCBI没有的，或者下载失败的，可以去EBI搜索
# https://www.ebi.ac.uk/ena/data/view/SRR7722939
ascp -QT -l 300m -P33001 -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh \
	era-fasp@fasp.sra.ebi.ac.uk:vol1/srr/SRR772/009/SRR7722939 ./



###############################################
## 1-2: Raw data convert
###############################################

# sra to fastq
cat $wkd/raw/P2586-4/SRR_Acc_List-2586-4.txt |while read i ; \
do \
time fastq-dump --gzip --split-files -A $i ${i}.sra && echo "** ${i}.sra to fastq done **" ;\
done

cat $wkd/raw/P9245-3/SRR_Acc_List-9245-3.txt |while read i ; \
do \
time fastq-dump --gzip --split-files -A $i ${i}.sra && echo "** ${i}.sra to fastq done **" ;\
done
# --gzip将生成的fastq文件压缩
# --split-3表示双端测序
# -A指定输出的文件名

# rename fq to cellranger format in batch 
# e.g. 
# /Sample1
    #- Sample1_S1_L001_I1_001.fastq.gz
    #- Sample1_S1_L001_R1_001.fastq.gz
    #- Sample1_S1_L001_R2_001.fastq.gz

cat  SRR_Acc_List-9245-3.txt | while read i ;do \
	(mv ${i}_1*.gz ${i}_S1_L001_I1_001.fastq.gz;\
	mv ${i}_2*.gz ${i}_S1_L001_R1_001.fastq.gz;\
	mv ${i}_3*.gz ${i}_S1_L001_R2_001.fastq.gz);done


###############################################
## 1-3: Raw data QC
###############################################
cd $wkd/qc
# for patient P2586-4
find $wkd/raw/P2586-4 -name '*R2*.gz'>P2586-4-id-2.txt

cat P2586-4-id-2.txt | xargs fastqc -t 10 -o ./

# for patient P9245-3
find $wkd/raw/P9245-3 -name '*R2*.gz'>P9245-3-id-2.txt

cat P9245-3-id-2.txt| xargs fastqc -t 20 -o ./

###############################################
## 2-1: Software Prepare
###############################################
cd $wkd/biosoft
# 2.0版本下载(732M)
# https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/2.0
tar zxvf cellranger-2.0.2.tar.gz 

# 添加到环境变量
export PATH=$PATH:$wkd/biosoft/cellranger-2.0.2
source ~/.bashrc

# 软件自检
cellranger testrun --id=tiny

###############################################
## 2-2: Reference Prepare
###############################################
cd $wkd/ref
# 下载参考序列
curl -O http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-1.2.0.tar.gz
# 然后解压
tar -xzvf refdata-cellranger-GRCh38-1.2.0.tar.gz

# 软件构建注释(中间a--ttribute是可选的)
# mkgtf <input_gtf> <output_gtf> [--attribute=KEY:VALUE...]
cellranger mkgtf Homo_sapiens.GRCh38.84.gtf Homo_sapiens.GRCh38.84.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene

# 软件利用构建好的注释，去构建需要的基因组
cellranger mkref --genome=GRCh38 \
                 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh38.84.filtered.gtf \
                 --ref-version=1.2.0

# 如果对于多个物种组合(本文的数据其实就应该这样组合起来)
cellranger mkref --genome=hg19 --fasta=hg19.fa --genes=hg19-filtered-ensembl.gtf \
                 --genome=mm10 --fasta=mm10.fa --genes=mm10-filtered-ensembl.gtf

cellranger mkref --genome=hg38 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    --genes=Homo_sapiens.GRCh38.84.filtered.gtf \
    --genome=mcv --fasta=mcv.fasta --genes=mcv_filter.gtf
# 结果保存在hg38_and_mcv目录中

##############################################
## 我看到这里写了这么多gene_biotype（也就是基因的生物类型）的键值对，不禁好奇，GTF中存在多少种基因类型，分别有多少种？
cat Homo_sapiens.GRCh38.84.filtered.gtf  |grep -v "#" |awk -v FS='gene_biotype ' 'NF>1{print $2}'|awk -F ";" '{print $1}'|sort | uniq -c

    213 "IG_C_gene"
     33 "IG_C_pseudogene"
    152 "IG_D_gene"
     76 "IG_J_gene"
      9 "IG_J_pseudogene"
   1209 "IG_V_gene"
    646 "IG_V_pseudogene"
    125 "TR_C_gene"
     16 "TR_D_gene"
    316 "TR_J_gene"
     12 "TR_J_pseudogene"
    848 "TR_V_gene"
    110 "TR_V_pseudogene"
  45662 "antisense"
  58181 "lincRNA"
2337766 "protein_coding"
##############################################

##############################################
## 或者自己尝试构建参考数据

# 下载基因组(这里以人类为例)
wget ftp://ftp.ensembl.org/pub/release-84/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

# (如果可以下载到注释的话)下载注释
wget ftp://ftp.ensembl.org/pub/release-84/gtf/homo_sapiens/Homo_sapiens.GRCh38.84.gtf.gz
gunzip Homo_sapiens.GRCh38.84.gtf.gz
# 注意：count只识别gtf，如果下载了gff，这里推荐使用genometools进行转换
# http://genometools.org/pub/binary_distributions/gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz
#  gt gff3_to_gtf xxx.gff >xxx.gtf 2>/dev/null
# 其他方法如：cufflinks下面的gffread，可能导致生成的gtf有的record没有gene_id
# 简单说一下使用方法：gffread xxx.gff -T -o xxx.gtf

# (如果不能下载到注释的话)自己构建注释【文章数据--mcv病毒的5个基因】
ADE45414.1_1	 Gnomon  exon    465    1190    .       -       .       gene_id "1"; transcript_id "1.1";
ADE45415.1_2     Gnomon  exon    1156   2427    .       -       .       gene_id "1"; transcript_id "1.1";
ADE45416.1_3     Gnomon  exon    2503   4722    .       -       .       gene_id "1"; transcript_id "1.1";
ADE45416.1_3     Gnomon  exon    5154   5387    .       -       .       gene_id "1"; transcript_id "1.1";
ADE45417.1_4     Gnomon  exon    4827   5387    .       -       .       gene_id "1"; transcript_id "1.1";
# 注意：自己构建要保证GTF文件分隔符是tab，如何检查？
awk -F '\t' '{print NF}' mcv.gtf #返回值为9正确
# 如果中间出现了多个空格而非一个tab，需要替换掉
sed 's/ \+ /\t/g' mcv.gtf > new_mcv.gtf

##############################################
















