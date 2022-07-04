# Tutorial

## Full-site Test
[Download here](http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_FullSiteTest.tar.gz). Execute the following commands in order to perform a motif matching followed by motif enrichment:

```shell
cd RGT_MotifAnalysis_FullSiteTest
rgt-motifanalysis matching --input-files input/regions_K562.bed input/background.bed 
rgt-motifanalysis enrichment input/background.bed input/regions_K562.bed
```

In the enrichment step, the order of the bed files matters. The background must always come first.

You can also reduce the amount of motifs used via the **–filter** parameter:

```shell
rgt-motifanalysis matching --filter "species:sapiens;name:EG" [other arguments..]
rgt-motifanalysis enrichment --filter "species:sapiens;name:EG" [other arguments..]
```

This will restrict the search to only the motifs with the string “EG” in their name (by default, the search is inexact; use **–filter-type** to explore different modes) and whose “species” metadata contains the “sapiens” string.

See the help file for all the available keys:

```
rgt-motifanalysis matching --help
```

You should also look into the metadata files – [read here](https://reg-gen.readthedocs.io/en/latest/motif_analysis/additional_motif_data.html) for more information.

## Promoter Test

[Download here](http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_PromoterTest.tar.gz). Execute the following commands in order to perform a motif matching followed by motif enrichment:
```shell
cd RGT_MotifAnalysis_PromoterTest
rgt-motifanalysis matching --target-genes input/genes.txt --input-files input/background.bed
rgt-motifanalysis enrichment input/background.bed match/target_regions.bed
```

In the enrichment step, the order of the bed files matters. The background must always come first.

### Gene association Test
[Download here](http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_GeneAssocTest.tar.gz). Execute the following commands in order to perform a motif matching followed by motif enrichment:

```
cd RGT_MotifAnalysis_GeneAssocTest
rgt-motifanalysis matching --input-matrix input_matrix.txt --rand-proportion 10
rgt-motifanalysis enrichment --input-matrix input_matrix.txt match/random_regions.bed
```

The matching command will read the experimental matrix, which specifies the PATH to the genomic regions and the genes to make an association test on. It also creates a background made of random regions of size 10 times the biggest genomic region in input. It might take between 10 and 30 minutes to run. It if is taking too long, use **--rand-proportion 1**.

The enrichment command will calculate the enrichment statistics for the all the input regions over the random background. This step should take about 5 minutes to complete.

Further usage instructions are found [here](https://reg-gen.readthedocs.io/en/latest/motif_analysis/tool_usage.html).