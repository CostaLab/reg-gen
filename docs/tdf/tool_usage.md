# Tool Usage

TDF includes the following submodes:

* **promotertest**: Promoter test evaluates the association between the given lncRNA to the target promoters.
* **regiontest**: Genomic region test evaluates the association between the given lncRNA to the target regions by randomization.
* **get_dbss**: Get TTSs in BED format from the single BED file
* **integrate**: Integrate the project’s links and generate project-level statistics.

Here we introduce the common parameters for two main tests (promoter test and genomic region test) first and then describe their test-specific parameters. The last three scripts as the tools are introduced afterward.

## Common Inputs for both tests
TDF can be executed with the following command:
```shell
rgt-TDF {promotertest,regiontest} [required inputs] [options]
```

Where:

* **{promotertest,regiontest}**: Define the applying test, either promoter test, or genomic region test.
* **required inputs**: Required inputs files and paths.
* **options**: Additional input parameters or output options.
There are some inputs common for both tests shown below:

### Required Input for both tests
| **Option Name** | **Type** | **Description** |
| :-------------- |:-------- | :--------------- |
| -h, –help       |      | Show the help message and exit |
| -r      |   PATH   | Input file name for RNA sequence (in fasta format)|
| -rn|String|Define the RNA name|
|-o|PATH|Output directory name for all the results and temporary files|
|-organism|String|Define the organism (hg19, hg38, mm9, mm10… etc)|

### Options
| **Option Name** | **Type** | **Default** | **Description** |
| :-------------- |:-------- | :--------------- | :--------------- |
|-t	| String	| RNA name	|Define the title name for the results under the Output name.
| -a	| Float	| 0.05	|Define significance level for rejection null hypothesis
| -ccf	 |Integer	| 20|	Define the cut off value for promoter counts
|-rt	| Boolean	| False	|Remove temporary files (fa, txp…etc)
| -log	 |Boolean	| False	|Set the plots in log scale
|-ac	 |PATH	| None	|Input file for RNA accecibility
|-accf	 |Integer	| 500	|Define the cut off value for RNA accecibility
|-obed	 |Boolean	| False	|Output the BED files for DNA binding sites.
| -showpa	| Boolean|	 False|	Show parallel and antiparallel bindings in the plot separately
|-filter_havana	| Boolean|	 False|	Apply filtering to remove HAVANA entries.
|-protein_coding|	 Boolean|	 False|	Apply filtering to get only protein coding genes.
|-known_only	| Boolean	| False|	Apply filtering to get only known genes.
|-nofile	| Boolean	| False|	Don’t save any files in the output folder, except the statistics.

### Options for TRIPLEXES
The arguments of the TRIPLEXES can be adjusted by the options below.
| **Option Name** | **Type** | **Default** | **Description** |
| :-------------- |:-------- | :--------------- | :--------------- |
| -l	| Integer	| 20|	Define the minimum length of triplex
| -e	| Integer	| 20|	Set the maximal error-rate in % tolerated
| -c	 |Integer	| 2	|Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as were found to greatly destabilize triplexes in vitro.
|-fr	| String	| off	|Activates the filtering of low complexity regions and repeats in the sequence data
|-fm	 |Integer	| 0|	Method to quickly discard non-hits (Default 0).’0′ = greedy approach; ‘1’ = q-gram filtering.
|-of	| Integer	| 1	|Define output formats of Triplexator
| -mf	| Boolean	| False|	Merge overlapping features into a cluster and report the spanning region.
|-rm	| Boolean	| False|	Set the multiprocessing
|-par	| String	| False	|Define other parameters for TRIPLEXES. Please ignore the first “-” and replace space with underline. For example, when you want to add “-G 80 -g 20”, please do “-par G_80_-g_20”.

## Particular Inputs for promoter test

### Required Input for promoter test
The target promoters can be defined in two ways:

* A gene list, which contains gene symbols or Ensembl IDs, one gene per line in plain text format. The argument, -de, should be used;
* Two BED files containing the regions of target promoters and non-target promoters (background). Two arguments, -bed and -bg, should be used together.

| **Option Name** | **Type** | **Description** |
| :-------------- |:--------| :--------------- |
|-de	| PATH|	Input file for gene list (gene symbols or Ensembl ID)
| -bed	| PATH	|Input BED file of the promoter regions of genes
|-bg	| PATH	|Input BED file of the promoter regions of background genes

### Options for promoter test

| **Option Name** | **Type** | **Default** | **Description** |
| :-------------- |:-------- | :--------------- | :--------------- |
| -pl	| Integer	| 1000	|Define the promotor length
|-score	| Boolean	| False|	Load score column from input gene list of BED file for analysis.
|-scoreh	| Boolean	| False|	Use the header of scores from the given gene list or BED file.

## Particular Inputs for region set test

### Required Input for region set test
| **Option Name** | **Type** | **Description** |
| :-------------- |:--------| :--------------- |
| -bed	 |PATH	|Input BED file for interesting regions on DNA

### Options for region set test
-mp Integer 0Define the number of threads for multiprocessing.

| **Option Name** | **Type** | **Default** | **Description** |
| :-------------- |:-------- | :--------------- | :--------------- |
| -n	 |Integer	| 10000|	Number iterations (randomization)
| -f	 |PATH	| None|	Input BED file as mask for randomization
|-score	 |Boolean	| False|	Load score column from input BED file

### get_ttss
Get TTSs of the given RNA sequence with the single BED file.
```shell
rgt-TDF get_ttss [options]
```
| **Option Name** | **Type** | **Default** | **Description** |
| :-------------- |:-------- | :--------------- | :--------------- |
|-h,| –help	||		show this help message and exit
|-i	 |PATH	||	Input BED file of the target regions
| -tts|	 PATH	||	Output BED file of the TTSs
|-tfo|	 PATH	||	Output BED file of the TFOs
|-tfo|	 PATH	||	Output BED file of the TFOs
|-r	 |PATH	||	Input FASTA file of the RNA
| -organism|	 PATH	||	Define the organism
| -l	| Integer	| 20|	**Triplexes** Define the minimum length of triplex
|-e	| Integer	| 20|	**Triplexes** Set the maximal error-rate in % tolerated
| -c	| Integer	| 2	|**Triplexes** Sets the tolerated number of consecutive errors with respect to the canonical triplex rules as such were found to greatly destabilize triplexes in vitro
| -fr	| on/off	| off	|**Triplexes** Activates the filtering of low complexity regions and repeats in the sequence data
| -fm	| Integer|	 0|	**Triplexes** Method to quickly discard non-hits (default: 0).’0′ = greedy approach; ‘1’ = q-gram filtering.
| -of	| Integer	| 1	|**Triplexes** Define output formats of Triplexes
| -mf	| Boolean	| False	|**Triplexes** Merge overlapping features into a cluster and report the spanning region.
| -rm	 |Integer	| 1|	**Triplexes** Set the multiprocessing

### integrate

Integrate the project’s links and generate project-level statistics.
```shell
rgt-TDF integrate [options]
````
| **Option Name** | **Type** | **Default** | **Description** |
| :-------------- |:-------- | :--------------- | :--------------- |
| -h, –help		|||	show this help message and exit
|-path	| PATH	||	Define the path of the project.
|-exp	| PATH	||	Include expression score for ranking.




