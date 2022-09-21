# Tutorial of regions versus signals

In this tutorial, we will demonstrate how we can use RGT-Viz to visualize signlas in different regions.

## Download the example files
We will use the epigenetic data from dendritic cell development study as example. There, we have ChIP-Seq data from the transcription factor **PU.1** and **IRF8**, and histone modifications **H3K4me1**, **H3K4me3**, **H3K9me3**, **H3K27me3**, and **H3K27ac** on four cellular states: multipotent progenitors (MPP), dendritic cell progenitors (CDP), common dendritic cells (cDC) and plamatocyte dendritic cells (pDC). The functional annotation of these histone markers are as follows:

* H3K4me1 is enriched at active and primed enhancers;
* H3K4me3 is highly enriched at active promoters near Transcription start site (TSS);
* H3K9me3 is a marker of heterochromatin which has pivotal role during lineage commitement;
* H3K27me3 is associated with the downregulation of nearby genes via the formation of heterochromatic regions;
* H3K27ac is accociated with the higher activation of transcription and defined as an active enhancer marker.

For simplicity, here we only look at data from the differentiation transition CDP to cDC.

Please follow the following steps to download the necessary example files: genomic signals in bigwig (bw) format for histone modifications and PU.1, as well as peaks from PU.1 (bed files).

Download the folder “rgt_viz_example” from [here](https://costalab.ukaachen.de/open_data/RGT/rgt_viz_example.zip).


## Creating Line Plots with RGT-Viz

Here we demonstrate how we can use RGT-Viz for drawing a lineplots. This allows for example to inspect the ChIP-Seq signals around particular genomic regions, as PU.1. peaks. Before you proceed, please [install RGT-Viz](https://reg-gen.readthedocs.io/en/latest/rgt/installation.html).

### Download the example files

Please follow the following steps to download the necessary example files: genomic signals in bigwig (bw) format for histone modifications and PU.1, as well as peaks from PU.1 (bed files).

Download the folder “rgt_viz_example” from [here](https://costalab.ukaachen.de/open_data/RGT/rgt_viz_example.zip).

```shell
cd rgt_viz_example
```

Now you have the files as described below:

```shell
rgt_viz_example
    ├── Matrix_CDP.txt
    ├── Matrix_CDP_cDC.txt
    ├── Matrix_H3K4me3.txt
    ├── Matrix_PU1.txt
    ├── data
    │   ├── CDP_H3K27me3.bw
    │   ├── CDP_H3K4me1.bw
    │   ├── CDP_H3K4me3.bw
    │   ├── CDP_PU1.bw
    │   ├── CDP_PU1_peaks.bed
    │   ├── CDP_WT_H3K27me3.bw
    │   ├── CDP_WT_H3K4me1.bw
    │   ├── CDP_WT_H3K4me3.bw
    │   ├── H3K4me3
    │   │   ├── CDP_WT_peaks.bed
    │   │   ├── MPP_WT_peaks.bed
    │   │   ├── cDC_WT_peaks.bed
    │   │   └── pDC_WT_peaks.bed
    │   ├── cDC_H3K27me3.bw
    │   ├── cDC_H3K4me1.bw
    │   ├── cDC_H3K4me3.bw
    │   ├── cDC_PU1.bw
    │   ├── cDC_PU1_peaks.bed
    │   ├── cDC_WT_H3K27me3.bw
    │   ├── cDC_WT_H3K4me1.bw
    │   └── cDC_WT_H3K4me3.bw
```

These files include the genomic signals of histone modifications (files with a .bw ending) and the genomic regions of PU.1 peaks (files with .bed endings).

### Understand experimental matrix

Before we use the RGT-Viz, you must define an experimental matrix. This tab separated file includes information necessary for RGT to understand your data, i.e. file paths, protein measured in the ChIP-Seq experiment, type of file and so on.

For example “<em>Matrix_CDP.txt</em>”  includes the files, which we need for finding the association of genomic signals on the genomic peaks of PU.1 transcription factor.

| name       	 | type	   | file			| factor   | cell |
| :------------- | :------ | :------------------------- | :------- | :--- |
| CDP\_PU1\_peaks| regions | ./data/CDP\_PU1\_peaks.bed | PU1      | CDP  |
| CDP\_PU1       | reads   | ./data/CDP\_PU1.bw		| PU.1     | CDP  |
| CDP\_H3K4me1   | reads   | ./data/CDP\_H3K4me1.bw	| H3K4me1  | CDP  |
| CDP\_H3K4me3   | reads   | ./data/CDP\_H3K4me3.bw	| H3K4me3  | CDP  |
| CDP\_H3K27me3  | reads   | ./data/CDP\_H3K27me3.bw 	| H3K27me3 | CDP  |

The first column (name) is a unique name for labeling the data; the second column indicate the type of experiment. Here we have either  “regions” (genomic regions in bed format) or “reads” (genomic signals in bigwig or bam format). The third column is the file path to the data. You can include additional columns to annotate your data.  In our example, the 4th column (factor) indicates the protein measured by the ChIP-Seq and the 5th collumn indicates the cell, where experiments were performed. You can add any more columns and the column names identify the feature.

### Run RGT-Viz

After defining the experiment matrix, now you can simply run RGT-Viz under “<em>rgt\_viz\_example</em>” directory by:

```shell
rgt-viz lineplot Matrix_CDP.txt -o results -t lineplot_CDP
```

- Matrix\_CDP.txt is the experimental matrix which contains the design of the data;
- -o indicates the output directory;
- -t defines the title of this experiment.

This command will generate a directory “<em>results</em>” with figures and html pages.

### Line plot

You can check the result by opening <em>results/index.html</em>

<p align="center">
<img src="../_static/rgt-viz/lineplot_CDP.png" width="450" height="350" align="center">
</p>

This lineplot shows the genomic signals of different histone modifications on the PU.1 genomic regions. The histone modifications are shown in different colors, and the window is centered by the midpoint of each genomic regions from PU.1 peaks.

### Add one more cell type

Lineplot is designed to compare more categories of data. Here we show another example to include one more cell type, cDC.

```shell
rgt-viz lineplot Matrix_CDP_cDC.txt -o results -t lineplot_CDP_cDC -col cell -row regions -srow
```

- <em>Matrix\_CDP\_cDC.txt</em> is the experimental matrix which contains the design of the data;
- -col defines the way to group data in columns, here we use “cell”, which is one of the headers in <em>Matrix\_CDP\_cDC.txt</em>;
- -row defines the way to group data in rows, here we use “regions”;
- -sx shares the y-axis for the plots in the same row.

<p align="center">
<img src="../_static/rgt-viz/lineplot_CDP_cDC.png" width="900" height="350" align="center">
</p>

This lineplot shows the difference of histone signatures on the PU.1 peaks among two cells. This plot indicates an increase in PU.1 and H3K4me3 levels on cDC cells compared to CDP cells.

For better comparison of each genomic signal, we can also plot them in different way, such as:

```shell
rgt-viz lineplot Matrix_CDP_cDC.txt -o results -t lineplot_CDP_cDC_2 -c cell -row reads -col regions -sx
```


- -c defines the way to color the lines, here we use “cell” as the tag to show different cells in different colors;
- -row defines the way to group data in rows, here we use “reads”;
- -col defines the way to group data in columns, here we use “regions”.

<p align="center">
<img src="../_static/rgt-viz/lineplot_CDP_cDC_2.png" width="450" height="540" align="center">
</p>

This design offer better comparison between cells by separating different histone modification and show cells in different colors.

Therefore, by changing the experimental matrix or the way to present, you can generate more complicated lineplot for comparison of your data across cell types, treatments, histone modification, or any other designs. RGT-Viz allows several other plots variants.

## References

1. Lin Q, Chauvistre H, Costa IG, Mitzka S, Gusmao EG, Haenzelmann S, Baying B, Hennuy B, Smeets H, Hoffmann K, Benes V, Sere K, Zenke M, Epigenetic and Transcriptional Architecture of Dendritic Cell Development, Nucleic Acids Research, 43:9680-9693, [[paper]](http://nar.oxfordjournals.org/content/early/2015/10/15/nar.gkv1056.full)[[data]](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64767)[[genome tracks]](http://www.molcell.rwth-aachen.de/dc/)







