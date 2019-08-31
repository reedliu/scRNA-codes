# 单细胞实战(六) 试着自己去探索cellranger count过程

> 刘小泽写于19.5.14
> 想要知根知底，就不能满足于跑一个简单的流程，需要从头到尾把每一步走一遍，踩踩坑。因为原文数据太大，所以这里选择了一个小的数据作为测试
> <https://isugenomics.github.io/bioinformatics-workbook/dataAnalysis/RNA-Seq/Single_Cell_RNAseq/Chromium_Cell_Ranger.html>

我们需要用到4个sra文件：SRR8111691到SRR8111694，为了保证代码的重复性，虽然只有四个文件，但还是使用脚本构建名称

### 数据下载

![image-20190514172821530](单细胞实战(六) 试着自己去探索cellranger count过程.assets/image-20190514172821530.png)

```shell
wkd=/YOUR_PATH/
mkdir -p $wkd/{raw,ref,qc,count}
cd raw
# 制作一个sra下载id号
for i in $(seq 1 4);do echo SRR811169$i;done >id.txt
# 下载 sra(每个文件3G左右)，需要提前装好sratools和aspera
cat id.txt | while read i;do prefetch $i -O ./; echo "** ${i}.sra done **";done
# sra转为faastq(全部大约有16G)
cat id.txt | while read i;do \
time fastq-dump --gzip --split-files -A $i ${i}.sra && echo "** ${i}.sra to fastq done **";\
done
# 将数据名称转换
mv  SRR8111691_1.fastq.gz SRR8111691_S1_L001_R1_001.fastq.gz
mv  SRR8111691_2.fastq.gz SRR8111691_S1_L001_R2_001.fastq.gz
mv  SRR8111692_1.fastq.gz SRR8111692_S1_L002_R1_001.fastq.gz
mv  SRR8111692_2.fastq.gz SRR8111692_S1_L002_R2_001.fastq.gz
mv  SRR8111693_1.fastq.gz SRR8111693_S1_L003_R1_001.fastq.gz
mv  SRR8111693_2.fastq.gz SRR8111693_S1_L003_R2_001.fastq.gz
mv  SRR8111694_1.fastq.gz SRR8111694_S1_L004_R1_001.fastq.gz
mv  SRR8111694_2.fastq.gz SRR8111694_S1_L004_R2_001.fastq.gz
```

### 下载参考信息

这里选取了有点难度的，因为从官网下载是很方便的，就不做过多解释。我们这里是要自己下载fasta和gff(注意：还不是平常使用的gtf哦)

```shell
cd ref
# 基因组fasta(112M)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/224/145/GCF_000224145.3_KH/GCF_000224145.3_KH_genomic.fna.gz
gunzip GCF_000224145.3_KH_genomic.fna.gz
# 注释文件gff(154M)
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/224/145/GCF_000224145.3_KH/GCF_000224145.3_KH_genomic.gff.gz
gunzip  GCF_000224145.3_KH_genomic.gff.gz
```

**CellRanger的要求比较严格，它需要gtf格式**，那么究竟为什么只要gtf不要gff格式呢？就要稍微了解一下这二者的区别

#### GFF与GTF的概念

- **GFF：General Feature Format**，目前常用3版本，即GFF3。主要用来注释基因组，包含tab分隔的9列：

  ```shell
  1. seq_id:序列编号(可以是chr或者scaffold)
  2. source：注释来源(一般是数据库或机构)，未知用.表示
  3. type: 注释类型，如Gene、cDNA、mRNA、CDS等
  4. start：基因/转录本在参考序列起始位置(1-based)
  5. end：（同上）终止位置
  6. score：（数字表示）得分，如序列相似性的E值或者基因预测的P值，未知用.表示
  7. strand：基因/转录本在参考序列的+/-链
  8. phase：当第3列type为编码蛋白的CDS时，这一列才有意义。表示下一个密码子开始的位置，或者说到下一个密码子需要跳过的碱基数，用0、1、2表示
  # 其中，前8列在gff的三个版本(gff1、gff2、gff3)中记录的信息都是一样的，只是名称有差异：比如第一列虽然都记录序列编号，但gff1叫seqname，gff2叫reference  sequence；type在gff1、gff2中叫feature；phase在gff1、gff2中叫frame
  9. attributes：键值对表示的属性，格式是：key=value，不同的键值对用分号分隔；一个键内的多个值用逗号隔开。其中有一些已经定义好的键名：如
  	ID表示注释类型(也就是type的名称)；
  	Name：注释名称，可以重复；
  	Parent：type属于的上一级名称，如exon上一级是transcript，而transcript的上一级是gene
  ```

- **GTF：General Transfer Format**，主要对基因进行注释，目前主要使用gtf2。与GFF有两个主要的地方存在不同之处：

  1. 第三列的feature中一定包含CDS、start_condon、stop_codon
  2. 第九列的attribute中必须以gene_id和transcript_id开头，并且键值对之间用空格隔开，而不是用等号；另外每个键值对之间用分号分开，并且最后一个键值对的末尾也要有分号

  ```shell
  NC_020166.2	Gnomon	exon	2398	3338	.	-	.	gene_id "1"; transcript_id "1.1";
  NC_020166.2	Gnomon	exon	3781	3884	.	-	.	gene_id "1"; transcript_id "1.1";
  ```

#### GTF与GFF的转换

如果以后做项目时发现，自己研究的物种没有gtf，只有gff，不要管，只管下载下来，然后进行转换即可，而转换工具也有不少

##### 工具一：gffread

```shell
# 它是由cufflinks开发的，使用如下
gffread -T xxx.gff  -o xxx.gtf
# 生成的结果
NC_000001.11    BestRefSeq      exon    11874   12227   .       +       .       transcript_id "rna0"; gene_id "gene0"; gene_name "DDX11L1";
```

生成的GTF文件中只有`exon`与`CDS`信息，第九列只包含`transcript_id`、`gene_id`、`gene_name` 

##### 工具二：genometools中的gff3_to_gtf 

```shell
# 工具下载安装
http://genometools.org/pub/binary_distributions/gt-1.5.10-Linux_x86_64-64bit-complete.tar.gz
# 它会将一些非exon、CDS的gene_type作为warning信息，可以直接不输出这些无用信息
gt gff3_to_gtf xxx.gff >xxx.gtf 2>/dev/null
```

#### 参考信息使用

> 这个教程是**基于CellRanger 3.0.2版本**，如果之前安装过2.X版本，只需要将它的路径加在bashrc中的2.X版本下方，就可以成功调用版本3了

```shell
# 根据gene_biotype来过滤gtf(这里没有选择具体的gene_biotype，是因为它没有gene_biotype，只有gene_id和transcript_id)
cellranger mkgtf GCF_000224145.3_KH_genomic.gtf GCF_000224145.3_KH_genomicfiltered.gtf
# 参考基因组
cellranger mkref --genome C_robusta --fasta=GCF_000224145.3_KH_genomic.fna --genes=GCF_000224145.3_KH_genomic.gtf --ref-version=3.0.0
```

这里会同时尝试另一种方式，直接使用fasta加上转换后的gtf文件

---

> 以上都是准备工作，下面从质控开始

### 原始数据质控

```shell
cd qc
ls *.gz | xargs fastqc -t 10 -o ./
multiqc ./
```

结果导出看一看，测序reads质量还不错，都在Q30以上，图中红色和橙色的reads是R1

![](https://upload-images.jianshu.io/upload_images/9376801-0ba7e0afa814e4d2.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)







































