# Tool Usage

To run the motif analysis tool you have to specify the analysis type (matching or enrichment). Furthermore, in order to run the motif enrichment tool you should run the motif match in advance and then the motif enrichment.

## Motif Matching

Run this command to see all available options:
```shell
rgt-motifanalysis matching --help
```

### Input
Inputs are genomic regions in BED or BigBed format that you wish to perform motif matching on. If an experimental matrix is provided, the inputs are ignored. If a Gene List is provided, then both the inputs and the experimental matrix are ignored, since the regions will be extracted from the promoters of the organism’s genes.

### Output
This analysis will populate the specified matching folder with the following files:

* **X_mpbs.[bed or bb]**: BED or BigBed file containing the matched motif instances for each input region X (or for the selected genes’ promoter regions, if **--gene-list** is provided). MPBS stands for Motif-Predicted Binding Sites. If or **--make-background** are set, then either a **random_regions_mpbs.bed** or a **background_regions_mpbs.bed** file will also be present.

* **random_regions.[bed or bb]**: BED or BigBed file containing the random background regions (this file is generated only if random regions are requested with **--rand-proportion**), plus the corresponding _mpbs file.

* **background_regions.[bed or bb]**: BED or BigBed file containing the genomic regions corresponding to all genes of the target organism (this file is generated only if such a background is requested with **--make-background**), plus the corresponding _mpbs file.

## Motif Enrichment
Run this command to see all available options:

```shell
rgt-motifanalysis enrichment --help
```

### Input
A required input is the file with the regions to be used as a “background” in the enrichment statistics. The file must be in either BED or BigBed format, and it must have a corresponding MPBS file as produced from motif matching. If the background file is named “background.bed”, then the MPBS file must be named “background_mpbs.bed” and be located either in the matching location, or in the same directory as the background.

Input files can be specified just like for matching, in either BED or BigBed format. If an experimental matrix is provided, input files are ignored.

Alternatively, you can specify an experimental matrix that defines the regions in which the motif matching will be applied to. Learn more about the Experimental Matrix format here.

In case gene sets are being used, the experimental matrix may have one additional column. This column will group regions by gene sets. See the example below:
```
name type file genegroup
# Regions
H1-hESC regions ./Input/regions_H1hesc.bed geneset1
HeLa-S3 regions ./Input/regions_HeLaS3.bed geneset1
HepG2 regions ./Input/regions_HepG2.bed geneset2
K562 regions ./Input/regions_K562.bed geneset2
# Genes
UP_REG genes ./Input/up_regulated_genes.txt geneset1
DW_REG genes ./Input/down_regulated_genes.txt geneset2
```

In this example, H1-hESC and HeLaS3 cell types will be associated to the up-regulated (UP_REG) genes and HepG2 and K562 will be associated to down-regulated (DW_REG) genes. Commentary lines (starting with #) were added only for clarity but not required.

### Output
This analysis will populate the folder specified with **--output-location**. This folder will contain subfolders with names depending on the analysis type. If gene sets are provided, each folder will be named X__Y, where X is the region file name and Y is the gene set file name associated to that region. Else, if no gene sets are provided, each folder will be named after each region name only. Each of these folders will contain the following files:

* **coord_association.[bed or bb]**: File containing the association between gene sets and regions. If no gene sets was provided, this file will associate all regions with all genes for the organism being used.

* **mpbs_ev.[bed or bb]**: Contain all Motif-Predicted Binding Sites (MPBS) that occurred in regions that were matched to genes, in case gene sets were given; or in all regions in case no gene sets was provided.

* **mpbs_nev.[bed or bb]**: Contain all Motif-Predicted Binding Sites (MPBS) that occured in regions that were NOT matched to genes, in case gene sets were given. In case gene sets were not given, such output is not available.

* **fulltest_statistics.[html and txt]**: Contains the enrichment results in html and tab-separated text format. The results include, in order, the motif name, enrichment p-value, corrected enrichment p-value, Fisher’s exact test A value, Fisher’s exact test B value, Fisher’s exact test C value, Fisher’s exact test D value, foreground frequency (A/(A+B)), background frequency (C/(C+D)) and a list of genes that were matched to the regions in which that motif occurred separated by comma.

* **genetest_statistics.[html and txt]**: Same as the above files, but regarding the gene association test. Available only if gene sets are specified in the experimental matrix.