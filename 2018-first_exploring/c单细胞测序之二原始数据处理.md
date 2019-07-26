# 单细胞测序之二原始数据处理

> 这是单细胞测序第二部分，开始通过真实数据进行数据准备

本次的参考文献是2015年Molecular Cell的文献

```
http://dx.doi.org/10.1016/j.molcel.2015.04.005
```

### 原始数据下载

```shell
mkdir raw && cd raw
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR522/ERR522959/ERR522959_1.fastq.gz
ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR522/ERR522959/ERR522959_2.fastq.gz
```

### 数据质控

```shell
mkdir fastqc && cd fastqc
fastqc -t 10 -o /fastqc ../raw/ERR522959_1.fastq ../raw/ERR522959_2.fastq
```

结果显示，两个read测试文件用的是Nextera Transposase Sequence的接头，并且都存在污染![](https://upload-images.jianshu.io/upload_images/9376801-8828c676ff1c4772.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

---

### 数据过滤

使用trim-galore软件，它可以用来切掉存在接头污染或者低质量序列的read末尾，另外过滤后要再用fastqc质控一次

> 注意：安装时的名称是trim-galore，使用的时候要用trim_galore

```shell
mkdir trimmed_fastqc && cd trimmed_fastqc
trim_galore --nextera --fastqc -o ./ ../../raw/ERR522959_1.fastq ../../raw/ERR522959_2.fastq
## --nextera指定污染序列接头类型
## --fastqc在trim完后，设置再次进行fastqc
```

![](https://upload-images.jianshu.io/upload_images/9376801-9135b05a719da30d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

这样的reads才是接下来比对需要的

---

### 比对须知

> 比对才能获得测序片段在基因组或转录组上的位点，才能继续分析问题；
>
> 要比对reads，就需要有参考基因组和注释文件gtf/gff

- **对于模式物种**，参考基因组和注释文件都可以在三大数据库中下载到：Ensembl、NCBI、UCSC，但不同的是：
  - Ensemble是最容易使用且是最大的注释文件集
  - NCBI只有高质量的基因注释，也就是说，如果研究的物种比较小众，那么NCBI不会入选
  - UCSC也包括大量的基因集、注释集，但标准有点杂

目前做单细胞测序最多的也应该属人和小鼠了，因此参考文件很齐全

- 简单说明一下**GTF文件格式**：

  （1）seqname:染色体或scaffold编号；

  （2）source：注释来源【比如人类一般采用的都是HAVANA团队注释的】；

  （3）feature：这段区域代表什么基因结构（基因、转录本、外显子？）；

  （4）start；（5）end：feature的起始与终止坐标;

  （6）score：feature存在性和它坐标的可信度（可有可无）；

  （7）strand：+正链-负链；

  （8）frame：0表示从一个完整地密码子开始，起始于密码子5’端的第一个碱基；1表示在第一个完整密码子前有一个碱基；2表示前面还有2个碱基；而反向链，5’端的第一个碱基在end坐标处；

  （9）attribute：分号分隔的额外信息【可选】包括基因id、转录本id、基因生物类型等等内容

  【内容为空，用`.`标记】 

- **比对的目的**是对基因表达量进行定量或者找出样本间表达差异的基因

- **使用的软件一：STAR** 

  它的**原理**是找出比对到参考基因组中一条或多条序列的**最长的测序read** ；
  **示意图**如下：蓝色的是一条测序read，它横跨了两个外显子；黄色的是可变剪切区域；            
  STAR软件先找到read的第一部分，和参考序列的第一个外显子匹配，另外又发现read第二部分和第二个外显子匹配。它就是用这种方式发现可变剪切事件或者染色体重排的。
  **缺点**：当参考基因组很大的时候，比如人和小鼠，它是非常占内存的                       

  ![c3.png](https://upload-images.jianshu.io/upload_images/9376801-f706dcd4da79610d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

  第一步：构建索引

  ```shell
  mkdir index && cd index 
  ## vi create_index.sh
  fasta=/DIR to hg19/hg19.fa
  STAR --runThreadN 20 \ ##这里我用了72线程，总共使用了23分钟
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles $fasta
  ```

  ![](https://upload-images.jianshu.io/upload_images/9376801-3311613e3b65638d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

  第二部：比对

  ```shell
  mkdir star_align && cd star_align
  ## vi star_align.sh
  STAR --runThreadN 20 \
  --readFilesIn ../trimmed_fastqc/ERR522959_1_trimmed.fq ../trimmed_fastqc/ERR522959_1_trimmed.fq \
  --genomeDir ./ \
  --outFileNamePrefix ./STAR ## output files name prefix (including full or relative path)
  ```

- **使用的软件二：Kallisto** 

  > 上面用的STAR是将reads直接比对到参考序列，因此比对速度比较慢；而Kallisto是利用k-mer比对参考序列，在转录组测序中也可以用它来分析，它不需要拼接转录本，速度很快，也被称作"pseudo-aligner"

  **什么是k-mer？**

  k-mer是来自于一条read的长度为k的序列，例如有一条序列是ATCCCGGGTTAT，要得到7-mers可以从头开始数，每7个就是一个7-mers

  ![](https://upload-images.jianshu.io/upload_images/9376801-6a37687c25cb9e20.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

  **为何k-mers比直接比对reads速度要快？** 
  
1. 得益于假比对的算法，先准备不同基因的k-mers的索引，通过将read的k-mers和k-mers索引比较，从而对基因进行计数
  2. 处理错误碱基效率高，例如一条序列的开头比对时就出错了，如果是传统的直接比对，可能要反应一会进行纠错，但是k-mer直接跳过第一个7-mer，继续下面的几个7-mer比对

  **关于Kallisto**  

  该软件2015年被发表在NATURE BIOTECHNOLOGY上，“Near-optimal probabilistic RNA-seq quantification”，它处理一个3千万read的模拟数据只需要单核7.5分钟，另外单个模拟样本只需要3.2GB的内存。https://pachterlab.github.io/kallisto/singlecell.html

  假比对不像STAR可以比对参考基因组或者参考转录组，它的比对模式只是参考转录组，这也就是说，它将reads比对到了剪切转录本而不是基因，这种模式对于单细胞测序还是很有有一定问题的：
  
  1. 单细胞测序比大量细胞RNA测序（bulk-RNA-seq）覆盖度更低，我们能得到的有效信息也更少；
2. 许多单细胞测序手段都存在3‘端的偏差（bias），因此基本不能通过3’端辅助判断。如果两条转录本只在5‘端存在差异，那么就很难找到reads从那一条转录本而来；
  3. 一些单细胞测序数据reads长度比较短，不好判断reads来自哪个转录本

  Kallisto与一般的“假比对（pseudo-alignment）”的算法有一些差异，它不是直接比对的转录本，而是比对到了“equivalence classes”，大概意思就是：如果一条reads比对到好几个转录本，Kallisto将这个reads的比对记录存储在equivalence classes这个空间中，再后来进行下游分析估计表达量时，Kallisto不使用基因或转录本表达量，而使用这个equivalence classes counts。![](https://upload-images.jianshu.io/upload_images/9376801-387063099a917468.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

  
  
  