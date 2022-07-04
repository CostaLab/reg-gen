# Introduction

<img src="../_static/hint/logo.png" width="175" align="right">

## Method
HINT (**H**mm-based **I**de**N**tification of **T**ranscription factor footprints) is a framework that uses open chromatin data to identify the active transcription factor binding sites. This method is originally proposed to model the active binding sites by simultaneous analysis of DNase-seq and the ChIP-seq profiles of histone modifications on a genome-wide level (paper). 

The HMM has as input a normalized and a slope signal of DNase-seq and one of the histone marks. It can, therefore, detect the increase, top, and decrease regions of either histone modification and DNase signals. And we next modified HINT to allow only DNase-seq data by removing the three histone-level states and the use of bias-corrected DNase-seq signal before normalization steps. 

Recently, we extended HINT to ATAC-seq, a new assay to identify accessible DNA regions, taking the protocol-specificity into consideration.

## Basic Usage
We describe here how to detect footprints using HINT for ATAC-seq, DNase-seq, and histone modifications data. To perform footprinting, you need at least two files, one with the aligned reads of your chromatin data and another describing the regions to detect footprints. You can use a peak caller, such as [MACS2](https://github.com/macs3-project/MACS),  to define these regions of interest.


### Footprinting for ATAC-seq data
Download [here](https://costalab.ukaachen.de/open_data/hint/tutorial/HINT_ATACTest.tar.gz) the example data for ATAC-seq based on chromosome 1 of the GM12878 cell. Execute the following commands to extract the data from the download file:

```shell
tar xvfz HINT_ATACTest.tar.gz
cd HINT_ATACTest
```
and the below command to perform footprinting:

```shell
rgt-hint footprinting --atac-seq ATAC.bam ATACPeaks.bed
```
For simplicity, we use only the first 1000 peaks from chromosome 1. The above commands will output a BED file containing the footprints in your **current folder** with **footprints** as the prefix. Moreover, You can set the below arguments

```shell
--output-location=your_directory  --output-prefix=your_prefix
```

to tell HINT your preferred output directory and name. Each footprint, i.e. each line of the BED file, will contain information regarding the tag-count score (number of reads) of each footprint. This score can be used as a footprint quality assessment (the higher values indicates better candidates). In addition, a file including the details of reads and footprints will also be written in the same folder of BED file.

If your data is paired-end, you may want to try another model which is optimized for paired-end sequencing data:

```shell
rgt-hint footprinting --atac-seq --paired-end --output-prefix=fp_paired ATAC.bam ATACPeaks.bed
```

**Note**: HINT performs bias correction for ATAC-seq by default, so you must download the genomes following these instructions and correctly specify the genome references with the following command before footprinting:

```shell
--organism=genome_version
```
Currently, the default setting is **hg19**. Find here for more information.

### Footprinting for DNase-seq
You can find [here](http://134.130.18.8/open_data/hint/tutorial/HINT_DNaseTest.tar.gz) example DNase-seq data. Execute the following commands to extract the data from a compressed file:
```shell
tar xvfz HINT_DNaseTest.tar.gz
cd HINT_DNaseTest
```

and the following command to call the footprints:
```shell
rgt-hint footprinting --dnase-seq DNase.bam DNasePeaks.bed
```

We recommend you to use cleavage bias correction. This can be done by using the following command:

```shell
rgt-hint footprinting --dnase-seq --bias-correction DNase.bam DNasePeaks.bed
```

Don’t forget to define the proper genome references using :
```shell
--organism=genome_version
```
Currently, the default setting is **hg19**.

### Footprinting for histone modification data
Download [here](http://134.130.18.8/open_data/hint/tutorial/HINT_HistoneTest.tar.gz) the example data for histone modification. Execute the following commands to extract data:

```shell
tar xvfz HINT_HistoneTest.tar.gz
cd HINT_HistoneTest 
```

and call footprints
```shell
rgt-hint footprinting --histone histone.bam histonePeaks.bed
```


## Citation

If you use HINT with ATAC-seq should cite the following publication:

* Li, Z., Schulz, M. H., Look, T., Begemann, M., Zenke, M., & Costa, I. G. (2019). [Identification of transcription factor binding sites using ATAC-seq](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1642-2). Genome Biology, 20(1), 45.

HINT with DNase with bias correction should cite:

* Gusmao EG, Allhoff M, Zenke M and Costa IG. [Analysis of computational footprinting methods for DNase sequencing experiments](https://www.nature.com/articles/nmeth.3772). Nature Methods, 13(4):303-309, 2016.

HINT with DNase-seq or histones cite the following publication:

* Gusmao EG, Dieterich C, Zenke M and Costa IG. [Detection of active transcription factor binding sites with the combination of DNase hypersensitivity and histone modifications](https://academic.oup.com/bioinformatics/article/30/22/3143/2390674). Bioinformatics, 30(22):3143-3151, 2014.



