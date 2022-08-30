# Tool Usage

Command:

```shell
rgt-THOR [option] CONFIG
```

## Required Input
| **Option Name** | **Type** | **Description** |
| :-------------- |:-------- | :--------------- |
| CONFIG          | File     | Configuration file for the experiment |

### Options

| **Option Name** | **Type** | **Default** | **Description** |
| :--------------- | :-------- | :----------- | :--------------- |
| \-\-name | string | none | Experiment’s name and prefix for all files that are created. |
| \-\-merge | boolean | False |  	Merge peaks which have a distance less than the estimated mean fragment size (recommended for histone data). | 
| \-\-no-merge-bin | boolean | True | Merge the overlapping bin before filtering by p-value. |
| \-\-housekeeping-genes | file (bed format) | none | Define housekeeping genes (BED format) used for normalizing. |
| \-\-output-dir | string | none | Store files in output directory. |
| \-\-report | boolean | False | Generate HTML report about experiment. |
| \-\-deadzones | file (bed format) | none | Define blacklisted genomic regions avoided  for analysis (BED format). |
| \-\-no-correction | boolean | False | Do not use multipe test correction for p-values  (Benjamini/Hochberg). |
| \-\-pvalue | float | 0.1 | P-value cutoff for peak detection. Call only peaks  with p-value lower than cutoff. | 
| \-\-exts | interger1,interger2,... | estimate | Read’s extension size for BAM files (comma separated list  for each BAM file in config file). If option is not chosen, estimate extension sizes from reads. |
| \-\-factors-inputs | float1,float2, ... | estimate | Normalization factors for input-DNA (comma separated list  for each BAM file in config file). If option is not chosen, estimate factors. |
| \-\-scaling-factors | float1,float2, ... | estimate | Scaling factor for each BAM file (not control input-DNA)  as comma separated list for each BAM file in config file. If option is not chosen,  follow normalization strategy (TMM or HK approach) |
| \-\-save-input | boolean | False | Save input DNA bigwig (if input was provided). |

### Advanced options

| **Option Name** | **Type** | **Default** | **Description** |
| :--------------- | :-------- | :----------- | :--------------- |
| \-\-regions | file (bed format) | none | Define regions (BED format) to restrict the analysis, that is, where to train  the HMM and search for DPs. It is faster, but less precise. |
| \-\-binsize | integer | 100 | Size of bins for creating the signal. | 
| \-\-step | integer | 50 | Stepsize with which the window consecutively slides across the genome  to create the signal. |
| \-\-debug | boolean | False | Output debug information. Warning: space consuming! |
| \-\-no-gc-content | boolean | False | Do not normalize towards GC content. | 
| \-\-norm-regions | file (bed format) | none | Restrict normalization to particular regions (BED format). |
| \-\-threshold | float | 95 | Minimum signal support for differential peaks to define training set  as percentage (t<sub>2</sub>, see [paper](http://nar.oxfordjournals.org/content/early/2016/08/01/nar.gkw680.abstract)). |
| \-\-size | integer | 10000 | Number of bins the HMM's training set constists of. |
| \-\-par | integer | 1 | Percentile for p-value postprocessing filter. |
| \-\-poisson | boolean | False | Use binomial distribution as emmission. |
| \-\-single-strand | boolean | False | Allow single strand BAM file as input. | 
| \-\-m\_threshold | integer | 80 | Define the M threshold of percentile for training TMM. |
| \-\-a\_threshold | integer | 95 | Define the A threshold of percentile for training TMM. | 
| \-\-rmdup | boolean | False | Remove the duplicate reads. |

### Config File 

The config file contains all necessary information for the experiment.

The file contains several groups, each group is indicated by a header symbol (#) followed by the group’s name. The information of each group is then given in the following lines.

The headers are:

*   #rep1
*   #rep2
*   #chrom\_sizes
*   #genome (optional)
*   #inputs1 (optional)
*   #inputs2 (optional)

The headers #inputs1, #inputs2 and #genome are optional. All other headers are required. An example for a config file is given [here](https://reg-gen.readthedocs.io/en/latest/thor/introduction.html).

<b>\#rep</b>

Headers #rep1 and #rep2 give all BAM files to be analysed. All files are listed line by line for the first (#rep1) and second (#rep2) biological condition.

<b>\#chrom\_sizes</b>

The chromosome sizes is a tab limited file with the chromosome name and the size of the chromosome. To download the chromosome size file of an organism, follow these [instructions](http://wiki.bits.vib.be/index.php/FetchChromSizes).

<b>\#inputs</b>

For each BAM file, we can also provide control input-DNA BAM files by using headers #inputs1 and #inputs2. Input-DNA helps to handle bias in ChIP-seq profiles and can therefore improve the differential peak estimation.

<b>\#genome</b>

Input-DNA and the genome (in [fasta](http://blast.ncbi.nlm.nih.gov/blastcgihelp.shtml) format) is necessary to correct for GC-content. GC-content correction can lead to more precise DP estimates. For instance, to download the human genome hg19.fa, run

```shell
wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*'
gunzip *gz
cat *fa > hg19.fa
rm md5sum.txt chr*fa README.txt
```

## Output

By default, THOR creates four types of files. The file(s)

*   \*.bw give the postprocessed ChIP-seq signal (in bigWig format),
*   \*-setup.info gives information about the setting,
*   \*-diffpeaks.bed describes the estimated differential peaks in a proprietary BED format, and
*   \*-diffpeaks.narrowPeak describes the estimated differential peaks in narrowPeak format.

### setup.info file

The setup file contains information about the experiment and the configuration of THOR used to predict DPs.

This includes the extension sizes used for the reads, the normalization scaling factors for each single ChIP-seq signal and the overdispersion factor c describing the quadratic relationship between the mean and the variance in the data.

Moreover, the setup.info file gives information about the HMM configuration, that is, the transition matrix as well as the emission distribution. The emission distribution is a Negative Binomial distribution described by mu and alpha.

### BigWig files

THOR creates a BigWig file for  each ChIP-seq signal (\*s[1/2]rep[0/1/2/3/…].bw), that is, for each biological condition (s[1/2]) and each replicate (rep[0/1/2/3/…]).

### narrowPeak and BED files

As the narrowPeak format does not give the possibility to store the counts of condition 1 and condition  2 for each peak, we output an BED format with additional information in the 11th column.

The 11th column in the BED file gives a semicolon separated list for each differential peak. The first (second) element of the semicolon separated list contains a comma separated list of the counts of each replicate of the first (second) biological conditions. The third element of the list gives the calculated p-value.

For downstream analysis of this BED file, we provide two [tools](https://costalab.ukaachen.de/open_data/RGT/THOR/THOR-tools.tar.gz). The first tool separates the BED file by differential peaks that gain peaks in condition 1 and that gain peaks in condition 2. The second tool filters the BED file by p-value.

**Please note:** In both files, \*-diffpeaks.bed and \*-diffpeaks.narrowPeak, the strand column (6th column) indicates whether the peak gains condition 1 (positive strand) or gains condition 2 (negative strand).

To output valid narrowPeak and BED files, we add several dummy columns. Only then it is possible to use them in other tools, like e.g. the IGV browser.

For the BED file, columns 7 and 8, give the same information as columns 2 and 3; column 9 gives a colour code for the peaks (red for a differential peak in signal 1, and green for a differential peak in signal 2). Column 10 is always 0.

For the narrowPeak file, column 5, 7, 9 and 10 do not give any information.

## Advanced Example 

In this example, we study the Dendritic Cell development. Dendritic cells (DC) are professional antigen presenting cells that develop from hematopoietic stem cells in bone marrow. This system allows to differentiate ex-vivo mul-
tipotent progenitors (MPP) to DC progenitors (CDP). We want to call differential peaks in our provided ChIP-seq data for histone mark H3K27ac of MPP and CDP (see [here](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE73143) the GEO entry for the fastq files).

We create the config file DC\_THOR.config:

```shell
#rep1
MPP_WT_H3K27ac_1.bam
MPP_WT_H3K27ac_2.bam
#rep2
CDP_WT_H3K27ac_1.bam
CDP_WT_H3K27ac_2.bam
#genome
mm9.fa
#chrom_sizes
mm9.chrom.sizes
#inputs1
Input_MPP_H3K27ac.bam
Input_MPP_H3K27ac.bam
#inputs2
Input_CDP_H3K27ac.bam
Input_CDP_H3K27ac.bam
```

Here, we do not give only the required (#rep1, #rep2 and #chrom\_sizes), but also the optional header (#genome, #inputs1, #inputs2) to take advantage of input-DNA and GC-content normalization.

The following command finds differential peaks in the two conditions MPP and CDP.

```shell
rgt-THOR THOR_DC.config --report --housekeeping-genes hk_genes_new_promotor_mm9.bed --no-correction --output-dir ~/my_MPP_CDP_exp -n THOR_DC
```

THOR\_DC.config contains all information necessary to run THOR. With *\-\-output-dir* we determine the working directory <em>~/my\_MPP\_CDP\_exp</em>. THOR stops, if the folder already exists. All files created in <em>~/my\_MPP\_CDP\_exp</em> have the prefix *THOR\_DC* (specific by the option *-n*). The p-value is not corrected (option *\-\-no-correction*). Instead of using a TMM based approach to normalize the ChIP-seq profiles (default), we use a housekeeping genes approach. For that, please [download](https://costalab.ukaachen.de/open_data/RGT/THOR/hk_THOR.tar.gz) the housekeeping genes and use the option *\-\-house-keeping-genes*.

Moreover, we create a HTML report of our experiment (option *-\-report*). The HTML report contains useful information about our experiment, such as, the the experimental configuration (also described by _\*-setup.info_), the function THOR is using to model the mean-variance relationship, the function used to estimate the fragmentation size, and the quality check of the house-keeping-gene normalization approach.
