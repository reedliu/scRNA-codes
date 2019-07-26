# 单个细胞的测序？

> 刘小泽写于18.9.27 手痒又想接触新技术，看看前沿研究，充实自己

### 为何需要研究单个细胞？

> 世界上没有两片相同的树叶🍂

基因组、转录组测序是从大量细胞中获得一定量的DNA/RNA然后进行测序（可以理解为Bulk sequencing），结果当然也是反应细胞群体的特征。

> 比如有三组细胞，每组10个，其中第一组有一个致病能力很强悍的恶性细胞，第二组有5个致病能力中等的恶性细胞，第三组全部是致病能力较差的恶性细胞。这三组用一般的Western blot做出来的结果可能一样，因为他们整体情况差不多，但是就个体而言，我们可能对第一组那个很强悍的细胞感兴趣，为什么它能用一己之力，将整个群体带坏？

这就说明细胞是有**异质性**的（即：**相同表型的细胞遗传信息**由于生物过程的随机性和环境扰动的原因**可能存在显著差异**），有许多低丰度的信息会在分析整体中丢失。为此，在单个细胞水平上，对基因组、转录组、表观组进行高通量测序分析的单细胞测序方法诞生（**Single cell sequencing**）。这样就可以解释这样的现象：不同肿瘤细胞在生长速度、免疫特性、侵染能力等表型的差异，以及它们对不同肿瘤药物的敏感性不同。又或者可以研究一个组织的不同细胞表达和功能，例如血液中免疫细胞存在T细胞、B细胞、自然杀伤细胞、巨噬细胞、嗜酸粒细胞等，而且T细胞还能继续细分；胚胎、大脑的不同脑区的细胞转变更为复杂，因此需要研究某一种细胞。

单细胞测序开始应用于2010年，从2013年MALBAC技术应用以后开始迅速发展并被Nature Methods评为年度技术，2017年10月16日，与“人类基因组计划”相媲美的“人类细胞图谱计划”，其中就需要关键的“单细胞测序”技术；对于肿瘤高异质性、单核微生物细胞、大脑神经细胞、免疫细胞检测、重新定义细胞亚型、寻找罕见的细胞类别等研究都很有帮助。

> 来自Nature Method评论：“Every cell is unique—it occupies an exclusive position in space, carries distinct errors in its copied genome and is subject to programmed and induced changes in gene expression. Yet most DNA and RNA sequencing is performed on tissue samples or cell populations, in which biological differences between cells can be obscured by averaging or mistaken for technical noise.”

### 怎么得到单个细胞的信息？

#### 1. 获得单个细胞- 细胞分离

https://www.jianshu.com/p/9fc521bf82ba

首先需要单细胞悬液（如组织解离液或者细胞培养悬液）做检测样品

**（传统）有限稀释法：**不需要高级设备，通过稀释悬液（稀释完理论上只有一个细胞），移液器吸取一定体积悬液。效率仅20%左右

**显微操作法**：直接观察到细胞然后吸取，不利于高通量，并且实验技术要求高

**（较常用）流式细胞法**：FACS（Fluorescence activated cell sorting）利用流式细胞仪，根据细胞表面标记进行筛选，可以选出单个细胞或者特定细胞亚群。当然它需要较大的细胞量

**微流控技术：**细胞分选效果较好，但需要固定的芯片

#### 2. 实现高通量

早期（2010-2014）的Fluidigm C1平台（基于Smart-Seq + Fluidigm），最早被应用于单细胞领域，几个小时只能捕获96个左右的细胞；
后来2015年左右发展出了MARS-Seq、CytoSeq、Drop-Seq、inDrop技术；
2017-18年，10X Genomics、BD Rhapsody兴起【基于UMI（Unique Molecular Identifier）+ CL（Cell Label），UMI的引入可以使定量结果不受PCR bias的影响】。它们实现了测定1-6K不等的细胞量，覆盖大部分组织的单细胞群体类型

**10X Genomics**起源自Drop-Seq技术， 横向孔道逐个导入凝胶微珠Gel beads，第一个纵向道输入细胞。当凝胶微珠和细胞碰撞会被吸附在微珠上，然后通过微流控技术运送到第二个纵向通道（“油管”）。这时就会形成一个个的油滴GEMs（一个油滴就是一个凝胶微珠，也就是一个单细胞），然后收集在EP管中。每一个凝胶微珠都布满了不同的Barcode和UMI连接的序列，然后再加上PolyT就形成了像“刺”一样的捕获抓手，随后细胞裂解，利用3'端 poly(A) 碱基互补特定抓取mRNA构建转录文库。据说可以7分钟内完成100~80,000个细胞的捕获

![](https://upload-images.jianshu.io/upload_images/9376801-c5aea74fb3032a4a.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

比10X晚一年发布的**BD Rhapsody**平台源于Cytoseq蜂窝板技术，采用20W+的微孔（远超过细胞数量，目的就是为了一孔一个细胞），因此可以叫做“微孔捕获”，效率会相对更高一些。捕获后也是裂解细胞，再进行mRNA抓取。

今年还诞生了一种SPLiT-seq方法，可以兼容与甲醛固定的细胞或细胞核，并且不需要专用设备，省去用定制化微流体或微孔分选单细胞的过程，号称“1美分单细胞测序”。

#### 3. 怎么测序

一个细胞的DNA/RNA微乎其微，在皮克量级上（万亿分之一克），这么细微的量怎么能满足测序仪要求？因此，测序前需要扩增，使用全基因组扩增技术（WGA），主流三种方法：

**PCR扩增法：** 对特殊位点可能存在偏好，导致覆盖度偏低；

**等温法：**MDA（Multiple displacement amplification）为代表，覆盖度较好，但一致性较差；

**MALBAC：** 多重退火和成环循环扩增技术，覆盖度、一致性保持中等水平

### 单细胞分析内容

Marker基因（区分细胞用的标记基因）分析、细胞分群、细胞类型鉴定、拟时序分析（分析细胞分化轨迹）、基因相关性分析（基因调控网络）
详细一点，比如：外周血单个核细胞细胞亚型的分析（为发现细胞亚型新Marker基因）；研究干细胞不同亚型的特征，构建谱系，预测分化；寻找异常增殖的细胞，结合病理进行疾病分型；判断细胞嵌合状态，尤其是骨髓干细胞移植前后，来判断手术效果；新细胞发现

> 目前10X Gemonics提供了大量的测试数据，可以试一试。关于转录组的数据集放在这里：
> https://support.10xgenomics.com/single-cell-gene-expression/datasets
>
> 并且提供了处理流程Cell Ranger 2.1.0
>
> ![](https://upload-images.jianshu.io/upload_images/9376801-2c9a46d6d6c07e51.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

#### 关于10X处理软件的使用：

https://support.10xgenomics.com/single-cell-gene-expression/software/overview/welcome
处理+可视化三件套：Pipelines、Loupe™ Cell Browser、R Kit（下载需要简单的邮箱注册）2018-7-31刚更新了2.2.0版本

> **最低系统要求：**（比一般的转录组苛刻多了）
> 最低8核Intel/AMD（16核以上比较好）
> 内存64G起步（128以上比较好）
> 硬盘1T
> 64位Linux



![](https://upload-images.jianshu.io/upload_images/9376801-fb3413b65db9565b.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

流程大体上还是比对、定量、表达矩阵

```shell
#下载 Cell Ranger - 2.2.0（663M）
wget -O cellranger-2.2.0.tar.gz "http://cf.10xgenomics.com/releases/cell-exp/cellranger-2.2.0.tar.gz?Expires=1538100646&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cDovL2NmLjEweGdlbm9taWNzLmNvbS9yZWxlYXNlcy9jZWxsLWV4cC9jZWxscmFuZ2VyLTIuMi4wLnRhci5neiIsIkNvbmRpdGlvbiI6eyJEYXRlTGVzc1RoYW4iOnsiQVdTOkVwb2NoVGltZSI6MTUzODEwMDY0Nn19fV19&Signature=jQAzPkqHZ8fwjucWuT0viw1LkZxmOkpNbdMYkVLI9PBhNaQrAR9S1UlZgtG2mGXQGYOy3J3H3I-tprCaoZZkjUk-xn3cn-KV-MOF5dqn2BU6wLdX33fbw-jtreZqP160C6f0qryBZytMRXwDPCvRZGuYvWzdkD3ZROER2smCM9WeHsNLEaCGKfrnIqR5h5owWpbvYyi6OsLZYGudNk6soFp0G8IMDyDFrPUZMdHmN9hBTrwWWylFavwWL555Xtkdo3IIR33u5nfOsiOyz86GcPs1SXJdaTb4TzQ7qODMlVADWqcAkE2eZg5Yg0lxMiWM5L3BwDCsUGWClQlrwXD4wg__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"

# 下载人类参考基因组 GRCh38（11G）
nohup wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-GRCh38-1.2.0.tar.gz &
#人类 hg19（11G）
nohup wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-1.2.0.tar.gz &
# 小鼠（9.6G）
nohup wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-mm10-1.2.0.tar.gz &
# 人和小鼠（9.6G）
nohup wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-and-mm10-1.2.0.tar.gz &
# ERCC
wget http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-ercc92-1.2.0.tar.gz
# Chromium Single Cell 3' v2
wget http://cf.10xgenomics.com/supp/cell-exp/chromium-shared-sample-indexes-plate.csv
# Chromium Single Cell 3' v1
http://cf.10xgenomics.com/supp/cell-exp/chromium-single-cell-sample-indexes-plate-v1.csv
# Gemcode Single Cell 3'
http://cf.10xgenomics.com/supp/cell-exp/gemcode-single-cell-sample-indexes-plate.csv
```

除了这个配套的流程，还有一些R包可以进行下游分析，例如**Seurat包**
https://satijalab.org/seurat/pbmc3k_tutorial.html（教程更新到18-3-22）

可以进行单细胞测序数据的质控、过滤、标准化、减弱批次效应、PCA、tNSE亚群聚类分析、差异基因分析等【目前就学了聚类分析🧐路还很长哦】

```shell
#提供了处理好的表达矩阵( 10X得到的外周血单个核细胞 Peripheral Blood Mononuclear Cells，共有2700个细胞，用Illumina NextSeq 500进行测序)
wget https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz
```

![](https://upload-images.jianshu.io/upload_images/9376801-0e5ae9b5ca69ed97.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

---

> 刘小泽写于18.9.29
> 上一次是理论知识和准备工作，这次开始软件安装和测试

### 软件安装和检测

cellranger这个软件内容十分丰富，整合了大量的第三方工具，因此解压需要一段时间，解压完成后导入环境变量，按照官方要求，还要进行安装检测，看一下安装是否完整；另外把下载的数据库文件也解压一下

```shell
cd /db/10X
tar -xzvf refdata-cellranger-ercc92-1.2.0.tar.gz
tar -xzvf refdata-cellranger-hg19-1.2.0.tar.gz
tar -xzvf refdata-cellranger-hg19-and-mm10-1.2.0.tar.gz
cd /opt
tar -xzvf cellranger-2.2.0.tar.gz
export PATH=/opt/cellranger-2.2.0:$PATH
cellranger testrun --id=tiny # 32核检测大约8分钟,检查结束如下图，会生成tiny/tiny.mri.tgz这样的文件
```

![](https://upload-images.jianshu.io/upload_images/9376801-ffe430d65d5a08d4.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

cellranger主要包括：

```shell
# Usage:
cellranger mkfastq #将Illumina得到的原始BCL文件转为FASTQ

cellranger count # 比对、过滤、条形码和UMI计数
cellranger aggr # 针对多个样本的情况，把count合并而且标准化成相同的测序深度之后，再计算gene-barcode矩阵
cellranger reanalyze #将count或者aggr得到的gene-barcode 矩阵进行降维、聚类

# 10X Genomics的专属算法和RNA测序比对软件STAR结合，可以得到BAM、MEX、CSV、HDF5、HTML的标准格式的结果
```

![](https://upload-images.jianshu.io/upload_images/9376801-98debe36911d47d8.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### 下载测序数据

> cellranger要求fastq格式的数据，可以通过cellranger mkfastq转换、illumina的bcl2fastq转换、已发布数据集、cellranger bamtofastq转换得到

下载已有的数据集：https://support.10xgenomics.com/single-cell-gene-expression/datasets，选择小鼠1k Brain Cells from an E18 Mouse数据集，来自E18小鼠皮层、海马区和脑室下区，结果检测到了931个细胞

```shell
nohup wget http://cf.10xgenomics.com/samples/cell-exp/2.1.0/neurons_900/neurons_900_fastqs.tar &

# total 5.4G
37M Aug 25  2017 neurons_900_S1_L001_I1_001.fastq.gz
643M Aug 25  2017 neurons_900_S1_L001_R1_001.fastq.gz
1.8G Aug 25  2017 neurons_900_S1_L001_R2_001.fastq.gz
239M Aug 25  2017 neurons_900_S1_L002_I1_001.fastq.gz
646M Aug 25  2017 neurons_900_S1_L002_R1_001.fastq.gz
1.8G Aug 25  2017 neurons_900_S1_L002_R2_001.fastq.gz
```

文件的命名规则：`[Sample Name]` _S1_L00 `[Lane Number]`_`[Read Type]`_001.fastq.gz。
比如这里sample name是neurons_900，lane有两个1和2，
Read type有三种：`I1`Sample index read也就是cell-barcode；`R1`read1((UMI) reads)；`R2`read2

> 与普通fastq文件相比，单细胞RNASeq fastq文件包含条形码和唯一分子标识符（UMI）的额外信息。从文件大小也能看出来，只有read2是转录本序列

```shell
cellranger count --id= mm_neurons \ #生成的文件都放在这个名字的目录下（必选）
--fastqs=/project/scRNA-seq/10X/raw/neurons_900_fastqs \ #（必选）
--transcriptome=/db/10X/refdata-cellranger-mm10-1.2.0 \ #（必选）
--expect-cells=900 #（可选）期望得到的细胞数
--localcores 10 \ # CPU
```

如果数据包括许多sample，可以指定`--sample=SMAPLENAME`，另外还可以指定lane的编号，如`--lanes=1`

运行成功会提示：

![](https://upload-images.jianshu.io/upload_images/9376801-cd5b3713bd7c59b3.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

|              目录               |                             描述                             |
| :-----------------------------: | :----------------------------------------------------------: |
|            analysis             |            降维PCA、聚类、差异分析（全是CSV矩阵）            |
|          cloupe.cloupe          |              Loupe Cell Browser可视化及分析文件              |
|    filtered_gene_bc_matrices    |          过滤后的gene-barcode矩阵（只包含MEX格式）           |
| filtered_gene_bc_matrices_h5.h5 |             过滤后的gene-barcode矩阵（HDF5格式）             |
|        molecule_info.h5         | 使用cellranger aggr产生的信息，作用是把样本组合成更大的数据集 |
|    possorted_genome_bam.bam     |          reads比对到带有barcode注释的基因组和转录组          |
|  possorted_genome_bam.bam.bai   |                        bam的index信息                        |
|      raw_gene_bc_matrices       |                   未过滤的gene-barcode矩阵                   |
|        web_summary.html         |                      网页版总结（下图）                      |



![](https://upload-images.jianshu.io/upload_images/9376801-488116c07a066416.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

### Cellranger的一些知识

#### 比对流程

- 基因组比对：使用STAR将reads比对到基因组的过程是考虑剪切位点的，然后cellranger将转录组注释信息GTF分解成外显子、内含子以及基因间的区域，并给出比对类型的显著性。如果比对的位置与外显子有超过50%的交叉，那么就认为它比对到了外显子；如果不是外显子并且和内含子有交叉，就认为是内含子，否则就是基因间区域
- MAPQ调整：`MapQ = -10 log10(P)`，比如结果为30，那就是1/1000的概率会出现这个比对结果。对于比对到一个外显子位点但同时还比对到一或多个的非外显子位点，优先考虑比对到外显子，MAPQ 255值为255时果断认为read比对上了外显子
- 比对转录组：比对上外显子的read继续与有注释的转录本比对，寻找兼容性。与转录组的外显子匹配并且比对到同一条链上，就可以被认为比对到了转录组；如果只匹配一个基因的注释，那么它的比对是唯一的并且可信度高。只有比对到转录组可信度高的reads才能用于UMI计数

#### 了解下分子条形码/标签

> 分子条形码又称分子标签（MolecularBarcode, 又称UID Unique identifiers, **UMI Unique molecularidentifiers**）是对原始样本基因组打断后的每一个片段都加上一段特有的标签序列，来**区分同一样本中成千上万的不同的片段**，在后续的数据分析中可以通过这些标签序列来排除 DNA 聚合酶、扩增以及测序过程中所引入的错误

一般UMI由大约10nt的随机序列（如：NNNNNNNNN）或者简并碱基（根据密码子的兼并性,常用一个符号代替某两个或者更多碱基，如NNNRNYN）。它和样本标签（sample barcode）不同，**UMI是针对一个样本**的不同片段，而**样本标签**是为区分**不同样本** 加上的标签序列。

> 一个样本只能有一个相同的样品标签，但可以有成千上万的分子条形码

![](https://upload-images.jianshu.io/upload_images/9376801-0d3d7e55289c32a7.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

- 同一个样本的 DNA 片段，每一个片段都接上一个特定的标签序列；
- 随目标序列一起经过文库构建、PCR 扩增，然后被一同测序；
- 最终测序结果中，带有不同UMI的序列，代表它们来自不同的原始 DNA 片段分子；带有相同UMI的序列，表示它们是从同一条原始的 DNA 片段扩增而

**设置UMI目的**：PCR 和测序过程中的错误是随机发生的，根据UMI可以去除冗余，降低低频突变的假阳性率

