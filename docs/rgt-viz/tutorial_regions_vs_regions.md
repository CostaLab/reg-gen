# Tutorial of regions versus regions

In this tutorial, we will demonstrate how we can use RGT-Viz to visualize association among different region sets.

## Download the example files

This tutorial will use the same dataset as another tutorial. Please check the explanation and download it according to [Tutorial of regions versus signals](https://reg-gen.readthedocs.io/en/latest/rgt-viz/tutorial_regions_vs_signals.html).

Below are the related files for this tutorial:
```shell
rgt_viz_example
├── data
│   ├── bw
│   └── peaks
│       ├── H3K4me3_cDC_WT_peaks.bed
│       ├── H3K4me3_CDP_WT_peaks.bed
│       ├── H3K4me3_MPP_WT_peaks.bed
│       ├── H3K4me3_pDC_WT_peaks.bed
│       ├── PU1_cDC_peaks.narrowPeak
│       ├── PU1_CDP_peaks.narrowPeak
│       ├── PU1_MPP_peaks.narrowPeak
│       └── PU1_pDC_peaks.narrowPeak
├── Matrix_H3K4me3.txt
└── Matrix_PU1.txt
```

The idea is to investigate the association of PU1 peaks and H3K4me3 peaks in each cell type quantitively.

## Installing RGT

Here we demonstrate how we can use RGT-Viz for investigate the association among region sets. Before you proceed, please [install RGT-Viz](https://reg-gen.readthedocs.io/en/latest/rgt/installation.html).

### Understand experimental matrix

Before we use the RGT-Viz, you must define an experimental matrix. This tab separated file includes information necessary for RGT to understand your data, i.e. file paths, protein measured in the ChIP-Seq experiment, type of file and so on.

For example “<em>Matrix\_PU1.txt</em>”  includes the BED files of the genomic peaks of PU.1 transcription factor.

| name            | type    | file                                    | factor | cell |
| :-------------- | :------ | :---------------------------------------| :----- | :--- |
| MPP\_PU1\_peaks | regions	| ./data/peaks/PU1\_MPP\_peaks.narrowPeak |  PU1   |  MPP |
| CDP\_PU1\_peaks | regions	| ./data/peaks/PU1\_CDP\_peaks.narrowPeak |  PU1   |  CDP |
| cDC\_PU1\_peaks | regions	| ./data/peaks/PU1\_cDC\_peaks.narrowPeak |  PU1   |  cDC |
| pDC\_PU1\_peaks | regions	| ./data/peaks/PU1\_pDC\_peaks.narrowPeak |  PU1   |  pDC |


See [Tutorial of regions versus signals](https://reg-gen.readthedocs.io/en/latest/rgt-viz/tutorial_regions_vs_signals.html) for more details.

After defining the experiment matrix, now you can simply run RGT-Viz under “<em>rgt\_viz\_example</em>” directory with the following examples.

## Projection test

Projection test evaluates the association level by comparing to the random binomial model.

```shell
rgt-viz projection -r Matrix_H3K4me3.txt -q Matrix_cDC_pDC.txt -o results -t projection -c factor -organism mm9 -g cell
```

- -r is reference region set as the base for statistics;
- -q is query region set for testing its association with the reference regions;
- -o indicates the output directory;
- -t defines the title of this experiment;
- -c defines the color tag for cloring the test;
- -organism defines the genome assembly used here;
- -g defines the group tag for grouping the test.

This command will generate a directory “<em>results/projection</em>” with figures and html pages.

<p align="center">
<img src="../_static/rgt-viz/projection_test.png" height=500 align="center">
</p>

## Jaccard test

Jaccard test evaluates the association level by comparing with jaccard index from repeating randomization.

```shell
rgt-viz jaccard -r Matrix_H3K4me3_cDC_pDC.txt -q Matrix_cDC_pDC.txt -o results -t jaccard -organism mm9 -g cell
```

This command will generate a directory “<em>results/jaccard</em>” with figures and html pages.

<p align="center">
<img src="../_static/rgt-viz/jaccard_test1.png" width=500 align="center">
</p>

<p align="center">
<img src="../_static/rgt-viz/jaccard_test2.png" width=500 align="center">
</p>

## Intersection test

Intersection test provides various modes of intersection to test the association between references and queries.

Firstly, you can run intersection test without statistical test. This mode is faster to get a descriptive result:

```shell
rgt-viz intersect -r Matrix_H3K4me3_cDC_pDC.txt -q Matrix_cDC_pDC.txt -o results -t intersection -organism mm9
```

<p align="center">
<img src="../_static/rgt-viz/intersection_bar.png" width=500 align="center">
</p>


However, you can also run it with statistical test by randomization of query regions. This take more compuational resources.

```shell
rgt-viz intersect -r Matrix_H3K4me3_cDC_pDC.txt -q Matrix_cDC_pDC.txt -o results -t intersection -organism mm9 -stest 100
```

- -stest defines the repitition times of random subregion test between reference and query. The more repitition times are, the more reliable the result is. However, it take time to run.