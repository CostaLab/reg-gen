# Configuration of Genomic Data

When RGT is installed, it will automatically create a folder that stores additional data (default: ~/rgtdata). This data includes chromosome sizes, position frequency matrices (describing transcription factor motifs), HTML scripts, etc. Some tools require data too big to fit in the installation procedure, such as genomes and genomic annotations. In this section we will describe how to obtain these data.

## Automatic Data Setup
The easiest way to obtain all data sets required by RGT is to run the setupGenomicData.py python script inside the installed data directory. This will download the files from public servers and will take a few minutes. If you use MAC OS, make sure the command “wget” is available (further instructions here).

The following command will install all the necessary human genome (**hg19**) data sets:

```shell
cd ~/rgtdata
python setupGenomicData.py --hg19
```

The following command will install all the necessary mouse genome (**mm9**) data sets:
```shell
cd ~/rgtdata
python setupGenomicData.py --mm9
```

The following command will install all available data sets: (hg19, hg38, mm9, mm10, zv9, and zv10)
```shell
cd ~/rgtdata
python setupGenomicData.py --all
```

This script has further options that can be viewed with:
```shell
python setupGenomicData.py -h
```

## Customize RGT Data Folder

### The data.config File
The data.config file contains the default data set names (inside RGT data path) used by RGT tools. It is divided into sections (with labels in brackets), such as **GenomeData** and **MotifData**.

### GenomeData
For each genome assembly, there are five fields targeting to the relevant files. You can customize these paths by yourself. Below is the example for hg19:

| Field Name      | Default Value | Description |
| ----------- | ----------- | ----------- |
|genome       |genome_hg19.fa|  Sequence of assembly hg19 in FASTA format. This data set is not available upon installation. See instructions above on how to obtain this data set.|
|chromosome_sizes| chrom.sizes.hg19 |Chromosome sizes file of assembly hg19.|
|gene_regions|genes_hg19.bed|Gene locations in BED format (from Gencode annotation file in GTF format).|
|annotation|gencode.v19.annotation.gtf	|Gene annotation from Gencode version 19 for human in GTF format. This data set is not available upon installation. See instructions above on how to obtain this data set.|
|gene_alias|alias_human.txt|Alias file which allows for translation between multiple different gene IDs.|


**You should never modify the data.config file!** This is due to the fact that every RGT installation will overwrite it. You can however customize the data.config.user file, by copying a similar section from the data.config file and modifying it to your wishes. For example, to use data from the organisms the user is interested in studying you simply create a section with the genome name and define all the relevant paths.

For example, here is a customized genome for Arabis thaliana (TAIR10):
```
[tair10]
genome: path/to/genome_tair10.fa
chromosome_sizes: path/to/chrom.sizes.tair10
gene_regions: path/to/genes_tair10.bed
annotation: path/to/tair10.annotation.gtf
gene_alias: path/to/alias_tair10.txt
```

The files that should be defined include:

* Genome fasta file: These files must contain one sequence for each chromosome. Each sequence header must be the chromosome symbol (such as “chr1” for chromosome 1). It can be obtained from several resources, including the [UCSC Downloads Website](http://hgdownload.soe.ucsc.edu/downloads.html).

* Gene annotation file: It is a [BED](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file containing the genomic coordinates of each gene for the selected organism. It can be downloaded, among other places, in the [UCSC Table Browser](http://genome.ucsc.edu/cgi-bin/hgTables).

* Chromosome sizes: It is a tab-separated plain text file with two columns. The first must contain the chromosome alias and the second must contain the length of the chromosome in base pairs. It can be fetched for some organisms using the [fetchChromSizes](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/fetchChromSizes) script available at the [UCSC Utilities Website](http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/).

* GTF annotation file: A [GTF](http://www.ensembl.org/info/website/upload/gff.html) file.

* Gene ID/Symbol aliases: A tab delimited file with three columns. Each row describes a gene and its aliases. The first column contains the gene’s ENSEMBL ID. The second column contains the gene’s official symbol (or user’s symbol of preference). The third column contains an ampersand(&)-separated list of aliases.

### MotifData
The following table describes the data.config path fields:
| Field Name      | Default Value | Description |
| ----------- | ----------- | ----------- |
| pwm_dataset      | motifs       | Contains the path to the motif position  weight matrices (PWM) repositories.       |
| logo_dataset   | logos        | Contains the path to the logo graphs(graphical depiction of PWMs). This data set is not available upon installation. For more information on how to create this data set click here. |
| repositories | jaspar_vertebrates,uniprobe_primary | The PWM repositories that will be used in the analyses. It is a comma-separated list of folders inside <pwm_dataset> (see this option above) folder. For information on how to add additional repositories, click here. |

### RGT Data Folder Structure
After installation, the RGT data folder will contain the following data.

* **Organism Folders**: Currently, we provide data for Homo sapiens (**hg19**, **hg38**) , Mus musculus (**mm9**, **mm10**), and Danio rerio (**zv9**, **zv10**). Inside these folders, you can find information regarding gene annotation and chromosome sizes.

* **fig**: Default figures and HTML style files.

* **fp_hmms**: Contains default hidden Markov model files for HINT tool.

* **motifs**: Contain position weight matrices (PWMs) for vertebrates obtained in many different repositories.

Additional folders may be created regarding tool-specific data.