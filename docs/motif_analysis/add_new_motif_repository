# Add a new Motif Repository
There are two ways to add a motif repository: an easy one, and a more complex one. We will cover both, with their pros and cons, in this section.

In short, the easy way requires only adding the PWM files into a directory. The harder way requires obtaining and formatting the motif annotation.

## Easy way: no annotation

Go inside your **rgtdata/motifs** folder and create a directory where you’ll put the PWM files:
```shell
├── motifs
│   ├── hocomoco
│   ├── jaspar_vertebrates
│   ├── my_custom_repository
│   ├── uniprobe_primary
│   └── uniprobe_secondary
```

It is important that the PWM files be named as they are in the default repositories, with the “.pwm” extension, eg:

```shell
TF1.pwm
```

The format we use is the old JASPAR format (before 2016):
```
13 13 3 1 54 1 1 1 0 3 2 5
13 39 5 53 0 1 50 1 0 37 0 17
17 2 37 0 0 52 3 0 53 8 37 12
11 0 9 0 0 0 0 52 1 6 15 20
```

It must have exactly four rows (A, C, G, T from top to bottom – although we are adding support for two more letters, m and 1, for methylation).
The columns are the counts for each base of the motif, left to right.

**NB: if you want to add one of the JASPAR repositories, you can follow these steps:**

1. Go to http://jaspar.genereg.net/downloads and download the PFMs you want (one of the files marked as **Single batch file (txt)**)

2. Run the createPwm.py script we provide inside your RGTDATA folder (**motifs** subfolder), like this:

```shell
python createPwm.py -f jaspar-2016 -i jaspar_file.txt -o my_custom_repository
```

You can now create the logos as explained in the previous section:
```shell
python setupLogoData.py my_custom_repository
```

To use it in your experiments, you can pass it as an argument to the **matching** and **enrichment** commands:
```
rgt-motifanalysis matching --motif-dbs ~/rgtdata/motifs/my_custom_repository [other arguments..]
rgt-motifanalysis enrichment --motif-dbs ~/rgtdata/motifs/my_custom_repository [other arguments..]
```

There are two main issues with using this approach:

* the repositories set in your **data.config** and **data.config.user** are completely ignored

* the **–filter** argument utility is much reduced, because it relies on metadata (you can still filter on motif names, though)

A more organic approach is to add the repository and also add its corresponding **annotation**.

## Hard way: with annotation

The initial steps are the same as before, but we need to write an annotation file: **my_custom_repository.mtf**, still inside the motifs directory.

But what is an MTF file? It’s a format we use to unique all the scattered annotation formats each TF database has.

It is strictly a **tab separated file**, so careful with whitespaces. Let’s look at the first line in the **hocomoco.mtf** file:

```shell
AHR	AHR_HUMAN.H11MO.0.B	0.B	AHR	PAS domain factors	P35869	Integrative	vertebrates	Homo sapiens	2.8475,6.6065,8.544,11.1115,11.6185,12.7235
```

* The first field simply the clean “name” of the motif.
* The second field is the full, unique, original name of the motif.
* The third field is the version of this motif (it changes slightly across different repositories)
* The fourth field is the TF gene name, the so-called “gene symbol”
* The fifth field is the so-called “motif family”, a description that varies a lot across repositories.
* The sixth field is one or more Uniprot identifiers
* The seventh field is the data source (Chip-Seq, Selex, Integrative)
* The eight field is the taxonomic group
* The ninth field is the species (usually the full name, eg “Home sapiens”, not “hs”)
* The tenth field is a list of precomputed thresholds for several FPR values: 0.005, 0.001, 0.0005, 0.0001, 0.00005, 0.00001

The good news is that you don’t need to set them all. Any field you do not want to set you can replace with a single dot.

This is an acceptable mtf entry:

```shell
AHR	AHR_HUMAN.H11MO.0.B	0.B	AHR	.	.	.	.	.	.
```

After you have added the metadata for each of the PWMs inside **my_custom_repository**, you can set the repository inside your **data.config.user** file:
```shell
[MotifData]
repositories: my_custom_repository
```

This will make the next runs of rgt-motifanalysis use only and exclusively your custom repository. You can also add more than one repository, eg:
```shell
[MotifData]
repositories: my_custom_repository, jaspar_vertebrates, hocomoco
```

and then use the **--filter** parameter as appropriate to select one or more of those, as your experiment requires (this only needed during the **matching** stage):

```shell
rgt-motifanalysis matching --filter "database:my_custom_repository,hocomoco" [other arguments..]
```

This command would only use your custom repository and hocomoco, ignoring jaspar vertebrates without the need to modify the data.config.user file.


