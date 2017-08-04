# Python Libraries
from __future__ import print_function
import os
import re
import sys
import math
import glob
import pysam
import numpy
import argparse
import natsort
import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
from os.path import expanduser
home = expanduser("~")

# Local Libraries
from rgt.GeneSet import GeneSet
from rgt.CoverageSet import CoverageSet
from rgt.GenomicRegion import GenomicRegion
from rgt.AnnotationSet import AnnotationSet
from rgt.Util import OverlapType, GenomeData
from rgt.GenomicRegionSet import GenomicRegionSet, GRSFileIO
tag = "RGT-tools"


def get_sequence(sequence, ch, ss, es, reverse=False, complement=False, rna=False, ex=0, strand=None):
    import pysam
    sequence = pysam.Fastafile(sequence)
    if not ch:
        seq = sequence.fetch(max(0, ss - ex), es + ex)
    else:
        seq = sequence.fetch(ch, max(0, ss-ex), es + ex)
    seq = seq.upper()

    if strand == "-" and not reverse and not complement:
        reverse = True
        complement = True

    if complement:
        t = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N" }
        ns = ""
        for s in seq: ns += t[s]
        seq = ns

    if rna:
        t = {"T":"U", "A":"A", "C":"C", "G":"G"}
        ns = ""
        for s in seq: ns += t[s]
        seq = ns

    if reverse:
        seq = seq[::-1]
        return seq
    else:
        return seq

if __name__ == "__main__":
    ##########################################################################
    ##### PARAMETERS #########################################################
    ##########################################################################
    
    parser = argparse.ArgumentParser(description='RGT-convertor is for converting various data format \
                                                  in bioinformatic research\
                                                  Author: Chao-Chung Kuo\
                                                  \nVersion: 0.1.1', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    subparsers = parser.add_subparsers(help='sub-command help', dest='mode')

    ############### GTF add transcripts ######################################
    parser_gtfat = subparsers.add_parser('gtf_add_transcripts', 
                                         help="[GTF] Add transcripts from the existed exons")
    parser_gtfat.add_argument('-i', metavar='input', type=str, help="Input GTF file")
    parser_gtfat.add_argument('-o', metavar='output', type=str, help="Output GTF file")


    ############### GTF to BED ###############################################
    # python rgt-convertor.py gtf_to_bed -i -o
    parser_gtf2bed = subparsers.add_parser('gtf_to_bed',
                                           help="[GTF] Convert GTF file to BED by the given biotype")
    parser_gtf2bed.add_argument('-i', metavar='input', type=str, help="Input GTF file")
    parser_gtf2bed.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_gtf2bed.add_argument('-s', metavar='source', type=str, default="All", help="Define the source {ENSEMBL,HAVANA,All}")
    parser_gtf2bed.add_argument('-f', metavar='feature', type=str, default="gene",
                                help="Define the feature {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine,All}")
    parser_gtf2bed.add_argument('-t', metavar='type', type=str, default="All",
                                help="Define gene type e.g. 'protein_coding' more: http://www.gencodegenes.org/gencode_biotypes.html")
    parser_gtf2bed.add_argument('-st', metavar='status', type=str, default="All",
                                help="Define gene status {KNOWN, NOVEL, PUTATIVE,All}")
    parser_gtf2bed.add_argument('-g', metavar='gene', type=str, default=None,
                                help="Define the gene list for filtering, default is None.")
    parser_gtf2bed.add_argument('-id', action="store_true",
                                help="Use gene ID as region name, instead of gene symbol.")
    parser_gtf2bed.add_argument('-b', action="store_true",
                                help="Save exons into entries with block in BED")

    ############### GTF to FASTA #############################################
    # python rgt-convertor.py
    parser_gtf2fasta = subparsers.add_parser('gtf_to_fasta', 
                                             help="[GTF] Convert GTF file to FASTA (exons) by the given gene name")
    parser_gtf2fasta.add_argument('-i', metavar='input', type=str, help="Input GTF file")
    parser_gtf2fasta.add_argument('-o', metavar='output', type=str, help="Output FASTA file")
    parser_gtf2fasta.add_argument('-t', metavar='transcript', type=str, help="Define the target transcript")
    parser_gtf2fasta.add_argument('-g', metavar='gene', type=str, help="Define the target gene")
    parser_gtf2fasta.add_argument('-genome', metavar='   ', type=str, help="Define the FASTA file of the genome")
    
    ############### GTF add chr on each entry #################################
    # python rgt-convertor.py
    parser_gtfachr = subparsers.add_parser('gtf_add_chr', 
                                             help="[GTF] Add 'chr' to each line in GTF for proper chromosome name")
    parser_gtfachr.add_argument('-i', metavar='input', type=str, help="Input GTF file")

    ############### GTF get intergenic regions in BED #################################
    # python rgt-convertor.py
    parser_gtfintergenic = subparsers.add_parser('gtf_intergenic',
                                           help="[GTF] Generate BED files for exon, intron, and intergenic regions")
    parser_gtfintergenic.add_argument('-i', metavar='input', type=str, help="Input GTF file")
    parser_gtfintergenic.add_argument('-o', metavar='output', type=str, help="Output directory for BED file")
    parser_gtfintergenic.add_argument('-organism', metavar='  ', type=str, help="Define the organism")

    ############### BED add score ############################################
    # python rgt-convertor.py
    parser_bedac = subparsers.add_parser('bed_add_score', help="[BED] Add score column")
    parser_bedac.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedac.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedac.add_argument('-v', metavar='value', type=str, help="Define value to add")

    ############### BED merge  ############################################
    # python rgt-convertor.py
    parser_bedmerge = subparsers.add_parser('bed_merge', help="[BED] Merge regions by name")
    parser_bedmerge.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedmerge.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedmerge.add_argument('-s', action="store_true", help="Strand specific")
    parser_bedmerge.add_argument('-b', action="store_true", help="BED12 format")


    ############### BED merge by name ############################################
    # python rgt-convertor.py
    parser_bedmn = subparsers.add_parser('bed_merge_by_name', help="[BED] Merge regions by name")
    parser_bedmn.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedmn.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedmn.add_argument('-b', action="store_true", help="BED12 format")

    ############### BED rename ###############################################
    # python rgt-tools.py bed_rename -i -o -s -d -organism
    parser_bedrename = subparsers.add_parser('bed_rename', help="[BED] Rename regions by associated genes")
    parser_bedrename.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedrename.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedrename.add_argument('-s', action="store_true", help="Strand specific")
    parser_bedrename.add_argument('-d', action="store_true", help="Show the distance")
    parser_bedrename.add_argument('-organism', metavar='  ', type=str, help="Define the organism")
    parser_bedrename.add_argument('-l', metavar='length', type=int, default=1000,
                                  help="Define the length of promoter region (default:1000 bp)")
    parser_bedrename.add_argument('-t', metavar='threshold', type=int, default=50000,
                                  help="Define the threshold of distance (default:50000bp")
    parser_bedrename.add_argument('-target', metavar='  ', default=False, type=str,
                                  help="Target BED file")
    parser_bedrename.add_argument('-genes', metavar='  ', default=False, type=str,
                                  help="Target gene list")

    ############### BED change strand ###############################################
    # python rgt-convertor.py
    parser_bedchstrand = subparsers.add_parser('bed_change_strand', help="[BED] Change strand of regions by the target BED file")
    parser_bedchstrand.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedchstrand.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedchstrand.add_argument('-d', metavar='distance', type=int, default=0,
                                    help="Define the threshold of distance (default:0 bp")
    parser_bedchstrand.add_argument('-t', metavar='target', type=str, default=None, help="Target BED file")
    parser_bedchstrand.add_argument('-r', action="store_true", help="Reverse the strand")
    parser_bedchstrand.add_argument('-a', metavar='all', type=str, default=None,
                                    help="Define the stand for all regions")

    ############### BED extend ###############################################
    # python rgt-convertor.py
    parser_bedex = subparsers.add_parser('bed_extend', help="[BED] Extend the regions")
    parser_bedex.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedex.add_argument('-o', metavar='output', type=str, help="Output BED name.")
    parser_bedex.add_argument('-oz', "--onlyzero", action="store_true", default=False,
                              help="Extend only the zero-length regions")
    parser_bedex.add_argument('-len', metavar='length', type=int, help="Define the length to extend.")
    parser_bedex.add_argument('-c', action="store_true", default=False,
                              help="Extend from the center to both directions.")
    parser_bedex.add_argument('-l', action="store_true", default=False,
                              help="Extend from the left ends.")
    parser_bedex.add_argument('-r', action="store_true", default=False,
                              help="Extend from the right ends.")
    parser_bedex.add_argument('-up', action="store_true", default=False,
                              help="Extend from the upstream ends.")
    parser_bedex.add_argument('-down', action="store_true", default=False,
                              help="Extend from the downstream ends.")

    ############### BED subtract ###############################################
    # python rgt-convertor.py
    parser_bedsub = subparsers.add_parser('bed_subtract', help="[BED] Subtract the regions")
    parser_bedsub.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedsub.add_argument('-o', metavar='output', type=str, help="Output BED name.")
    parser_bedsub.add_argument('-t', metavar="target", type=str,
                               help="Define the target BED file to subtract.")
    parser_bedsub.add_argument('-all', action="store_true", default=False,
                               help="Subtract the whole region when it overlaps.")
    parser_bedsub.add_argument('-blocki', action="store_true", default=False,
                               help="Read the blocks in input.")
    parser_bedsub.add_argument('-blockt', action="store_true", default=False,
                               help="Read the blocks in target.")

    ############### BED cut ###############################################
    # python rgt-convertor.py
    parser_bedcut = subparsers.add_parser('bed_cut', help="[BED] Cut the regions")
    parser_bedcut.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedcut.add_argument('-o', metavar='output', type=str, help="Output BED name.")
    parser_bedcut.add_argument('-t', metavar="target", type=str,
                               help="Define the target BED file for cutting.")
    parser_bedcut.add_argument('-s', action="store_true", default=False,
                               help="Strand-specific.")

    ############### BED get promoters ########################################
    # python rgt-convertor.py bed_get_promoters -i -o -organism
    parser_bedgp = subparsers.add_parser('bed_get_promoters', 
                       help="[BED] Get promoters from the given genes")
    parser_bedgp.add_argument('-i', metavar='input', type=str, help="Input file (BED or gene list)")
    parser_bedgp.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedgp.add_argument('-organism', type=str, help="Define the organism (necessary if input is a gene list)")
    parser_bedgp.add_argument('-l', metavar='length', type=int, default=1000,
                              help="Define length of promoters (default:1000bp)")
    
    ############### BED get upstream regions #################################
    # python rgt-convertor.py bed_upstream -i -o
    parser_bedupstream = subparsers.add_parser('bed_upstream', 
                       help="[BED] Get regions upstream from the given BED file")
    parser_bedupstream.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedupstream.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedupstream.add_argument('-l', metavar='length', type=int, default=100, help="Define length (default:100bp)")
    parser_bedupstream.add_argument('-d', metavar='distance', type=int, default=100, help="Define distance (default:100bp)")
    parser_bedupstream.add_argument('-min', metavar='minimum', type=int, default=0,
                                    help="Define minimum length of gene to filter out the small genes (default:0)")
    parser_bedupstream.add_argument('-r', '--reverse', action="store_true", default=False,
                                    help="Reverse the strand.")

    ############### BED to FASTA #############################################
    parser_bed2fasta = subparsers.add_parser('bed_to_fasta', 
                       help="[BED] Export the sequences in FASTA according to the given BED file")
    parser_bed2fasta.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bed2fasta.add_argument('-o', metavar='output', type=str, help="Output directory for FASTA files")
    parser_bed2fasta.add_argument('-genome', metavar='   ', type=str, help="Define the FASTA file of the genome sequence")
    parser_bed2fasta.add_argument('-loci', action="store_true", default=False, help="Make genomic loci as sequence name")
    parser_bed2fasta.add_argument('-order', action="store_true", default=False, help="Make ranking number as sequence name")
    parser_bed2fasta.add_argument('-block', action="store_true", default=False,
                                  help="Read blocks")
    parser_bed2fasta.add_argument('-score', action="store_true", default=False,
                                  help="Load the score column in BED into FASTA")

    ############### BED filtered by gene name ################################
    parser_bed_filter = subparsers.add_parser('bed_filter',
                       help="[BED] Filter by the given gene list or minimal size or maximal size.")
    parser_bed_filter.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bed_filter.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bed_filter.add_argument('-gene', type=str, default=False, help="Define file for the gene list")
    parser_bed_filter.add_argument('-min', type=int, default=False, help="Define minimal length")
    parser_bed_filter.add_argument('-max', type=int, default=False, help="Define maximal length")
    parser_bed_filter.add_argument('-score', action="store_true", default=False, help="Add the score from gene list to BED file")
    parser_bed_filter.add_argument('-background', action="store_true", default=False,
                                   help="Get the genes not in the given gene list.")

    ############### BED remove if overlap ################################
    parser_bedro = subparsers.add_parser('bed_remove_if_overlap', 
                       help="[BED] Remove the regions if they overlap with the target regions")
    parser_bedro.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedro.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedro.add_argument('-t', metavar='target', type=str, help="Define BED file for target regions")
    parser_bedro.add_argument('-k', "--keep", action="store_true", default=False, help="Keep the overlapped regions, and remove the non-overlapped ones.")
    parser_bedro.add_argument('-b', "--block", action="store_true", default=False, help="Read and write BED12 format.")

    ############### BED add columns ################################
    parser_bedaddcol = subparsers.add_parser('bed_add_columns', 
                       help="[BED] Add extra columns to the BED file by gene name")
    parser_bedaddcol.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedaddcol.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedaddcol.add_argument('-ref', metavar='refer', type=str, help="Define file for referring the extra columns ")
    parser_bedaddcol.add_argument('-f', '--field', type=int, help="Which field of the reference file is compared for names.")

    ############### BED average size ################################
    parser_bedaversize = subparsers.add_parser('bed_size',
                                             help="[BED] Calculate the average size.")
    parser_bedaversize.add_argument('-i', metavar='input', type=str, help="Input BED file")

    ############### BED complementary regions from Genome ############################
    parser_bedcomplement = subparsers.add_parser('bed_complement',
                                               help="[BED] Get the complementary regions from genome.")
    parser_bedcomplement.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedcomplement.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bedcomplement.add_argument('-organism', type=str, help="Define the organism (necessary if input is a gene list)")

    ############### BED detect polyA reads within the regions ################################
    parser_bedpolya = subparsers.add_parser('bed_polya',
                                               help="[BED] Detect polyA reads within the regions.")
    parser_bedpolya.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bedpolya.add_argument('-b', metavar='bam', type=str, help="Input BAM file")
    parser_bedpolya.add_argument('-o', metavar='output', type=str, help="Output file")

    ############### BED to GTF ################################
    parser_bed2gtf = subparsers.add_parser('bed_to_gtf',
                                            help="[BED] Convert BED file to GTF format.")
    parser_bed2gtf.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bed2gtf.add_argument('-o', metavar='output', type=str, help="Output GTF file")

    ############### BED overlaps ################################
    parser_bedoverlap = subparsers.add_parser('bed_overlap', help="[BED] Output statistics of overlaps among the BED files.")
    parser_bedoverlap.add_argument('-i', metavar='input', type=str, help="Input BED files or directory")
    parser_bedoverlap.add_argument('-o', metavar='output', type=str, help="Output text file")

    ############### BED distance ###############################################
    # python rgt-convertor.py
    parser_beddis = subparsers.add_parser('bed_distance', help="[BED] Show the distance between two region sets")
    parser_beddis.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_beddis.add_argument('-o', metavar='output', type=str, help="Output table.")
    parser_beddis.add_argument('-t', metavar="target", type=str,
                               help="Define the target BED file to define the distance.")

    ############### BED standardize chromosome ################################
    parser_bedschrom = subparsers.add_parser('bed_standard_chrom',
                                            help="[BED] Standardize the chromosomes.")
    parser_bedschrom.add_argument('-i', metavar='input', type=str, help="Input BED file")

    ############### BED add overlapping region name ################################
    parser_adddata = subparsers.add_parser('bed_add_data',
                                             help="[BED] Add overlapping region name")
    parser_adddata.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_adddata.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_adddata.add_argument('-t', metavar='target', type=str, help="Target BED file")

    ############### BED sampling regions randomly ################################
    parser_sampling = subparsers.add_parser('bed_sampling',
                                           help="[BED] Sampling the regions in the given BED file randomly")
    parser_sampling.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_sampling.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_sampling.add_argument('-s', metavar='size', type=int, help="Number of the output regions")

    ############### BED bed12tobed6 ################################
    parser_bed12tobed6 = subparsers.add_parser('bed12tobed6',
                                            help="[BED] Convert BED12 to BED6")
    parser_bed12tobed6.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_bed12tobed6.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_bed12tobed6.add_argument('-e', action="store_true", help="Add exon number or not")

    ############### Divide regions in BED by expression #######################
    # python rgt-convertor.py divideBED -bed -t -o1 -o1 -c -m
    parser_divideBED = subparsers.add_parser('bed_divide', 
                       help="[BED] Divide the BED files by the expression.")
    parser_divideBED.add_argument('-bed', metavar='   ', type=str, help="Input BED file")
    parser_divideBED.add_argument('-t', metavar='table', type=str, help="Input expression table (Gene name should match the region name.")
    parser_divideBED.add_argument('-o1', metavar='output1', type=str, help="Output first BED file")
    parser_divideBED.add_argument('-o2', metavar='output2', type=str, help="Output second BED file")
    parser_divideBED.add_argument('-c', metavar='cutoff', type=int, help="Define the cutoff")
    parser_divideBED.add_argument('-m', metavar='mode', type=str, help="Define the mode, such as mean, max, or min.")

    ############### BAM Filter reads by BED file #######################
    parser_filterBAM = subparsers.add_parser('bam_filter',
                                             help="[BAM] Filter BAM file by the given regions in BED.")
    parser_filterBAM.add_argument('-i', metavar='input', type=str, help="Input BAM file")
    parser_filterBAM.add_argument('-bed', metavar='   ', type=str, help="Input BED file for the regions for filtering")
    parser_filterBAM.add_argument('-o', metavar='output', type=str, help="Output prefix for BAM file")

    ############### THOR MAplot ################################
    parser_thorma = subparsers.add_parser('thor_ma',
                       help="[THOR] Create the MA plot for understanding the effect of normalization.")
    parser_thorma.add_argument('-i', metavar='input', type=str, help="Input data config.")
    parser_thorma.add_argument('-thor', metavar='THOR result', type=str, default=".", help="Output directory of THOR.")
    parser_thorma.add_argument('-o', metavar='output', default=None, type=str, help="Output directory")
    parser_thorma.add_argument('-e', metavar='extension', default=None, type=str, help="Define the extension size.")
    parser_thorma.add_argument('-b', metavar='bin size', default=None, type=str, help="Define the bin size.")

    ############### THOR split and filter ################################
    parser_thorsf = subparsers.add_parser('thor_split', 
                       help="[THOR] Split and filter the differential peaks from rgt-THOR")
    parser_thorsf.add_argument('-i', metavar='input', type=str, help="Input BED file")
    parser_thorsf.add_argument('-o', metavar='output', default=None, type=str, help="Output directory.")
    parser_thorsf.add_argument('-p', metavar='p-value', type=int, help="Define the cut-off of p-value (-log10) for filtering.")
    parser_thorsf.add_argument('-fc', metavar='fold-change', type=int,default=0,
                               help="Define the cut-off of foldchange for filtering.")
    parser_thorsf.add_argument('-rn', '--rename', action="store_true",
                               help="Rename the peak names by associated genes.")
    parser_thorsf.add_argument('-g', metavar='genome', type=str, help="Define the genome")
    parser_thorsf.add_argument('-b', metavar='bin', type=int, help="Define the bin size")
    parser_thorsf.add_argument('-s', metavar='step', type=int, help="Define the step size")


    ############### GENOME get sequence ####################################################
    # python /projects/reg-gen/tools/rgt-tools.py getseq -d /data/rgt
    parser_getseq = subparsers.add_parser('getseq', 
                       help="[FASTA] Get sequence from genome FASTA")
    parser_getseq.add_argument('-b', metavar='bed', type=str, default=None, help="Input BED file")
    parser_getseq.add_argument('-o', metavar='output', type=str, default=None, help="Output FASTA file")
    parser_getseq.add_argument('-d', metavar='dna', type=str, help="DNA sequence in FASTA format")
    parser_getseq.add_argument('-p', metavar='pos', type=str, help="position")
    parser_getseq.add_argument('-s', metavar='strand', type=str, default="both", help="strand (+, -, or both)")
    parser_getseq.add_argument('-ch', metavar='chr', type=str, help="chromosome")
    parser_getseq.add_argument('-ss', metavar='start', type=int, help="start site")
    parser_getseq.add_argument('-es', metavar='end', type=int, help="end site")
    parser_getseq.add_argument('-ex', metavar='extention', default=0, type=int, help="extention")
    parser_getseq.add_argument('-re', '--reverse', action="store_true", help="Reverse the sequence")
    parser_getseq.add_argument('-c', '--complement', action="store_true", help="Get the complement of the sequence")
    parser_getseq.add_argument('-r', '--rna', default=False, action="store_true", help="Convert T to U")

    ############### WIG trim ends by chromosome #############################################
    parser_wig_trim = subparsers.add_parser('wig_trim_end', 
                       help="[WIG] Trim the WIG file according to the given chromosome size")
    parser_wig_trim.add_argument('-i', metavar='input', type=str, help="Input WIG file")
    parser_wig_trim.add_argument('-o', metavar='output', type=str, help="Output WIG file")
    parser_wig_trim.add_argument('-chrosize', metavar='  ', type=str, help="Define path to the chromosome size file")

    ############### GENE list convertion #############################################
    parser_ensembl2symbol = subparsers.add_parser('ensembl2symbol', 
                       help="[GENE] Convert the gene list from ensembl ID to gene symbol")
    parser_ensembl2symbol.add_argument('-i', metavar='input', type=str, help="Input gene list")
    parser_ensembl2symbol.add_argument('-o', metavar='output', type=str, help="Output gene list")
    parser_ensembl2symbol.add_argument('-organism', metavar='  ', type=str, help="Define the organism")

    parser_sumbol2ensembl = subparsers.add_parser('symbol2ensembl',
                                                  help="[GENE] Convert the gene list from gene symbol to ensembl ID")
    parser_sumbol2ensembl.add_argument('-i', metavar='input', type=str, help="Input gene list")
    parser_sumbol2ensembl.add_argument('-o', metavar='output', type=str, help="Output gene list")
    parser_sumbol2ensembl.add_argument('-organism', metavar='  ', type=str, help="Define the organism")

    ############### Get length in bp from FASTA files ################################
    parser_fasta2bp = subparsers.add_parser('fasta2bp',
                                            help="[FASTA] Get the length in bp of FASTA files")
    parser_fasta2bp.add_argument('-i', metavar='input', type=str, help="Input FASTA file or directory")
    parser_fasta2bp.add_argument('-o', metavar='output', type=str, default="", help="Output file with a table")

    ############### STAR junction to BED #############################################
    parser_circRNA = subparsers.add_parser('circRNA', 
                       help="[junction] Convert the Chimeric junction from STAR to BED file")
    parser_circRNA.add_argument('-i', metavar='input', type=str, help="Input chimeric junction file from STAR")
    parser_circRNA.add_argument('-t', metavar='tcons', type=str, help="Input BED file of tcons")
    parser_circRNA.add_argument('-o', metavar='output', type=str, help="Output BED file")
    parser_circRNA.add_argument('-c', metavar='circ', type=str, help="Output BED file of circular RNA")

    ############### FASTA slicing #############################################
    # python rgt-convertor.py sliceFASTA -i -o -l -p
    parser_sliceFASTA = subparsers.add_parser('sliceFASTA', 
                       help="[FASTA] Slice the sequence by given position and length")
    parser_sliceFASTA.add_argument('-i', metavar='input', type=str, help="Input FASTA file")
    parser_sliceFASTA.add_argument('-l', metavar='length', type=int, help="Length of the slice sequence")
    parser_sliceFASTA.add_argument('-o', metavar='output', type=str, help="Output FASTA file")
    parser_sliceFASTA.add_argument('-p', metavar='position', type=int, help="The start position")
    parser_sliceFASTA.add_argument('-r', '--reverse', default=False, action="store_true", help="Reverse the sequence")

    ############### TXP to BED #############################################
    # python rgt-convertor.py txp2bed -i -o
    parser_txp2bed = subparsers.add_parser('txp2bed',
                                           help="[BED] Convert TXP file into BED format")
    parser_txp2bed.add_argument('-i', metavar='input', type=str, help="Input TXP file")
    parser_txp2bed.add_argument('-o', metavar='output', type=str, help="Output BED file")

    ##########################################################################
    ##########################################################################
    ##########################################################################
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    elif len(sys.argv) == 2:  
        # retrieve subparsers from parser
        subparsers_actions = [action for action in parser._actions if isinstance(action, argparse._SubParsersAction)]
        # there will probably only be one subparser_action,but better save than sorry
        for subparsers_action in subparsers_actions:
            # get all subparsers and print help
            for choice, subparser in subparsers_action.choices.items():
                if choice == sys.argv[1]:
                    print("\nYou need more arguments.")
                    print("\nSubparser '{}'".format(choice))        
                    subparser.print_help()
        sys.exit(1)
    else:   
        args = parser.parse_args()

        try:
            if not os.path.exists(args.o.rpartition("/")[0]):
                os.makedirs(args.o.rpartition("/")[0])
            if "~" in args.i: args.i = args.i.replace("~", home)
            if "~" in args.o: args.o = args.o.replace("~", home)
            if "~" in args.genome: args.genome = args.genome.replace("~", home)
        except:
            pass
    ##########################################################################


    ############### GTF add transcripts ######################################
    if args.mode == "gtf_add_transcripts":
        print(tag+": [GTF] Add transcripts from the existed exons")
        print("input:\t" + args.i)
        #NT_166433   protein_coding  start_codon 12026   12028   .+0  gene_id "ENSMUSG00000000702"; transcript_id "ENSMUST00000105216"; exon_number "1"; gene_name "AC007307.1"; gene_biotype "protein_coding"; transcript_name "AC007307.1-201";
        #NT_166433   protein_coding  exon    16677   16841   .   +. gene_id "ENSMUSG00000000702"; transcript_id "ENSMUST00000105216"; exon_number "2"; gene_name "AC007307.1"; gene_biotype "protein_coding"; transcript_name "AC007307.1-201";
        #NT_166433   protein_coding  CDS 16677   16841   .   +0 gene_id "ENSMUSG00000000702"; transcript_id "ENSMUST00000105216"; exon_number "2"; gene_name "AC007307.1"; gene_biotype "protein_coding"; transcript_name "AC007307.1-201"; protein_id "ENSMUSP00000100851";

        entries = []
        with open(args.i) as gtf:
            for line in gtf:
                line = line.strip().split("\t")
                if line[2] == "exon":
                    adddata = line[8].split(";")
                    #print([ line[:7] + adddata])
                    entries.append( line[:8] + adddata )
        # sort
        entries.sort(key=lambda x: x[9])

        # iterate
        transcriptid = ""
        transcript = []

        for exon in entries:
            #print(exon[9])
            if transcriptid == exon[9]:
                left = min(left, int(exon[3]))
                right = max(right, int(exon[4]))
            else:
                if transcriptid != "":
                    transcript.append( [ pre_t[0], pre_t[1], "transcript", str(left), str(right), 
                                         pre_t[5], pre_t[6], pre_t[7], ";".join( pre_t[8:10] + pre_t[11:] ) ] )
                transcriptid = exon[9]
                left = int(exon[3])
                right = int(exon[4])
                pre_t = exon
                
        transcript.append( [ exon[0], exon[1], "transcript", str(left), str(right), 
                             exon[5], exon[6], exon[7], ";".join( exon[8:10] + exon[11:] ) ] )

        #print(transcript[0])
        print("\tNumber of entries: "+ str(len(entries)))
        print("\tNumber of transcripts: "+str(len(transcript)))
        print("output:\t" + args.o)
        with open(args.o, "w") as out:
            for t in transcript:
                print("\t".join(t), file=out)

    ############### GTF to BED ######################################
    elif args.mode == "gtf_to_bed":
        print(tag + ": [GTF] Convert GTF to BED")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("source:\t\t" + args.s)
        print("feature:\t" + args.f)
        print("type:\t\t" + args.t)
        print("status:\t\t" + args.st)
        if args.g: print("gene list:\t" + args.g)
        # Sorting all parameters

        if args.g:
            select_genes = GeneSet("genes")
            select_genes.read(args.g)

        regions = GenomicRegionSet("output")
        with open(args.i, "r") as f,open(args.o, "w") as g:
            for line in f:
                if line[0] == "#": continue
                line = line.strip().split("\t")
                if len(line) < 5: continue
                if args.s == "All": pass
                elif args.s != line[1]: continue
                if args.f == "All": pass
                elif args.f != line[2]: continue
                # Further details
                info = line[8].split("; ")

                gi = [s for s in info if "gene_id" in s][0].partition("\"")[2][:-1].partition(".")[0]
                try:
                    gs = [s for s in info if "gene_name" in s][0].partition("\"")[2][:-1].partition(" (")[0]
                except:
                    gs = gi

                if args.id: gn = gi
                else: gn = gs

                if args.t == "All": pass
                elif args.t != [s for s in info if "gene_type" in s][0].partition("\"")[2][:-1].partition(" (")[0]: continue

                if args.st == "All": pass
                elif args.st != [s for s in info if "gene_status" in s][0].partition("\"")[2][:-1].partition(" (")[0]: continue

                if line[0].isdigit(): ch = "chr" + line[0]
                else: ch = line[0]
                # seq = "\t".join([ch, line[3], line[4], gn, ".", line[6]])
                if not args.g:
                    # print(seq, file=g)
                    regions.add(GenomicRegion(chrom=ch, initial=int(line[3]), final=int(line[4]),
                                              name=gn, orientation=line[6]))
                elif select_genes.check(gs) or select_genes.check(gi):
                    # print(seq, file=g)
                    regions.add(GenomicRegion(chrom=ch, initial=int(line[3]), final=int(line[4]),
                                              name=gn, orientation=line[6]))

        regions.write(filename=args.o, io=GRSFileIO.Bed12 if args.b else GRSFileIO.Bed)
        print("Number:\t\t" + str(len(regions)))



    ############### GTF to FASTA #############################################
    elif args.mode == "gtf_to_fasta":
        print(tag+": [GTF] Export certain gene or transcripts into FASTA sequence")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("genome:\t" + args.genome)
        print("gene list:\t" + args.g)

        ann = AnnotationSet(gene_source=args.genome,
                            filter_havana=False, protein_coding=False, known_only=False)
        geneset = GeneSet("target")
        geneset.read(args.g)

        genes = ann.get_genes(gene_set = geneset)


    ############### GTF add chr ##############################################
    elif args.mode == "gtf_add_chr":
        print(tag+": [GTF] Add 'chr' to each entry")

        with open("temp.gtf", "w") as f:
            with open(args.i) as g:
                for line in g:
                    if line.startswith("#"):
                        print(line, file=f)
                    else:
                        print("chr"+line, file=f)
        # rewrite gtf file
        #os.move("temp.gtf")


    ############### GTF get intergenic regions in BED ########################
    elif args.mode == "gtf_intergenic":
        print(tag + ": [GTF] Generate BED files for exon, intron, and intergenic regions")
        ann = AnnotationSet(gene_source=args.i,filter_havana=False, protein_coding=False, known_only=False)
        genome = GenomicRegionSet(args.organism)
        genome.get_genome_data(organism=args.organism)
        exons = ann.get_exons()
        genes = ann.get_genes()
        introns = genes.subtract(exons)
        interg = genome.subtract(genes)
        if not os.path.exists(args.o):
            os.makedirs(args.o)
        exons.write(os.path.join(args.o, "exons.bed"))
        introns.write(os.path.join(args.o, "introns.bed"))
        interg.write(os.path.join(args.o, "intergenic.bed"))


    ############### BED add score ############################################
    elif args.mode == "bed_add_score":
        print(tag+": [BED] Add scores into BED format")

        with open(args.i) as f, open(args.o, "w") as g:
            for line in f:
                line = line.strip() 
                print(line+"\t"+args.v, file=g)


    ############### BED merge  ########################################
    elif args.mode == "bed_merge":
        print(tag + ": [BED] Merge regions")
        bed1 = GenomicRegionSet("input")
        bed1.read(args.i, io=GRSFileIO.Bed12 if args.b else GRSFileIO.Bed)
        bed1.merge(strand_specific=args.s)
        bed1.write(args.o, io=GRSFileIO.Bed12 if args.b else GRSFileIO.Bed)

    ############### BED merge by name ########################################
    elif args.mode == "bed_merge_by_name":
        print(tag+": [BED] Merge regions by name")

        bed1 = GenomicRegionSet("input")
        bed1.read(args.i, io=GRSFileIO.Bed12 if args.b else GRSFileIO.Bed)
        bed2 = bed1.mergebyname()
        bed2.write(args.o, io=GRSFileIO.Bed12 if args.b else GRSFileIO.Bed)

    ############### BED rename regions #######################################
    elif args.mode == "bed_rename":
        print(tag+": [BED] Rename regions by associated genes")

        if args.target:
            print("target:\t" + args.target)
        else:
            print("organism:\t" + args.organism)

        bed = GenomicRegionSet(args.i)
        bed.read(args.i)
        if args.target:
            target = GenomicRegionSet(args.target)
            target.read(args.target)
            bed.replace_region_name(regions=target)
            bed.write(args.o)
        else:
            if not args.genes:
                renamebed = bed.gene_association(gene_set=None, organism=args.organism,
                                                 promoterLength=args.l, strand_specific=args.s,
                                                 threshDist=args.t, show_dis=args.d)
            else:
                genes = GeneSet("genes")
                genes.read(args.genes)
                renamebed = bed.gene_association(gene_set=genes, organism=args.organism,
                                                 promoterLength=args.l, strand_specific=args.s,
                                                 threshDist=args.t, show_dis=args.d)

            renamebed.write(args.o)


    ############### BED change strands #######################################
    elif args.mode == "bed_change_strand":
        print(tag + ": [BED] Change strands by target BED file")

        bed = GenomicRegionSet(args.i)
        bed.read(args.i)
        # print(len(bed))
        if args.t:
            print("target:\t" + args.t)
            target = GenomicRegionSet(args.t)
            target.read(args.t)
            if args.d != "0":
                target.extend(left=int(args.d), right=int(args.d))
            bed.replace_region_strand(regions=target, reverse=args.r)
        elif args.a:
            bed.replace_region_strand(all=args.a)
        else:
            bed.replace_region_strand(regions=None, reverse=args.r)
        bed.write(args.o)


    ############### BED extend ###############################################
    elif args.mode == "bed_extend":
        print(tag + ": [BED] Entend the regions")
        bed = GenomicRegionSet("bed")
        bed.read(args.i)
        if args.l:
            bed.extend(left=args.len, right=0)
        if args.r:
            bed.extend(left=0, right=args.len)
        if args.up:
            bed.extend_upstream(length=args.len)
        if args.down:
            bed.extend_downstream(length= args.len)
        if args.c:
            bed = bed.relocate_regions(center="midpoint", left_length=args.len, right_length=args.len)

        bed.write(args.o)


    ############### BED subtract ###############################################
    elif args.mode == "bed_subtract":
        print(tag + ": [BED] Subtract the regions")
        bed = GenomicRegionSet("bed")
        bed.read(args.i)
        if args.blocki:
            bed.extract_blocks()
        target = GenomicRegionSet("target")
        target.read(args.t)
        if args.blockt:
            target.extract_blocks()
        out = bed.subtract(y=target, whole_region=args.all)
        out.write(args.o)


    ############### BED cut ###############################################
    elif args.mode == "bed_cut":
        print(tag + ": [BED] Cut the regions and neglect the down stream ones")
        bed = GenomicRegionSet("bed")
        bed.read(args.i)
        target = GenomicRegionSet("target")
        target.read(args.t)
        out = bed.cut_regions(target)
        out.write(args.o)

    ############### BED get promoters #########################################
    elif args.mode == "bed_get_promoters":
        print(tag + ": [BED] Get promoter regions from the genes")
        gene = GenomicRegionSet("genes")
        ### Input BED file
        if args.i.endswith(".bed"):
            gene.read(args.i)
            promoter = GenomicRegionSet("promoter")
            promoterLength = int(args.l)
            for s in gene:
                if s.orientation == "+": 
                    s.initial, s.final = max(s.initial-promoterLength, 0), s.initial

                else: s.initial, s.final = s.final, s.final+promoterLength
                promoter.add(s)
        ### Input gene list
        else:
            ann = AnnotationSet(gene_source=args.organism, alias_source=args.organism,
                                filter_havana=False, protein_coding=False, known_only=False)
            de_gene = GeneSet("de genes")
            de_gene.read(args.i)
            print(len(de_gene))
            promoter = ann.get_promoters(promoterLength=args.l, gene_set=de_gene, unmaplist=False)
            #print(len(de_prom))

        
        #print(len(promoter))
        promoter.write(args.o)



    ############### BED get upstream regions ####################################
    elif args.mode == "bed_upstream":
        print(tag + ": [BED] Get upstream regions from the given BED file")

        gene = GenomicRegionSet("genes")
        ### Input BED file
        
        gene.read(args.i)
        print(len(gene))
        target = GenomicRegionSet("target")
        # if args.min == 0: cut = float("inf")
        # elif args.min > 0: cut = args.min

        for s in gene:
            if s.orientation == "+" and len(s) > args.min: 
                s.initial, s.final = max(s.initial-args.d-args.l, 0), max(s.initial-args.d, 0)
                if s.initial > s.final: s.initial, s.final = s.final, s.initial
                if args.reverse: s.orientation = "-"
                target.add(s)

            elif s.orientation == "-" and len(s) > args.min: 
                s.initial, s.final = s.final+args.d, s.final+args.d+args.l
                if s.initial > s.final: s.initial, s.final = s.final, s.initial
                if args.reverse: s.orientation = "+"
                target.add(s)
        
        print(len(target))
        target.write(args.o)


    ############### BED to FASTA #############################################
    elif args.mode == "bed_to_fasta":
        print(tag+": [BED] BED to FASTA")

        if ".fa" not in args.o and not os.path.exists(args.o):
            os.makedirs(args.o)
        elif ".fa" in args.o:
            fasta = open(args.o, "w")
        regions = GenomicRegionSet("regions")
        regions.read(args.i)

        for region in regions:
            if "/" in region.name:
                region.name = region.name.replace("/", "_")
        if args.block:
            regions.extract_blocks()
        if args.order:
            ranking = []
            with open(args.i) as f:
                for line in f:
                    if line.startswith("#"): continue
                    else:
                        l = line.strip().split()
                        ranking.append([l[0],int(l[1]),int(l[2])])

        for i, r in enumerate(regions):
            if args.order:
                for j, reg in enumerate(ranking):
                    if reg[0] == r.chrom and reg[1] == r.initial and reg[2] == r.final:
                        name = "peak_"+str(j+1)
            elif args.loci:
                name = r.toString(underline=True, strand=True)
            else: name = r.name

            if r.data and len(r.data.split()) == 7:
                target = r.extract_blocks()
                #print("*** using block information in BED file")
                writelines = []
                for exon in target:
                    if not args.score:
                        writelines.append("> "+" ".join([exon.name, exon.toString(), exon.orientation]))
                    else:
                        writelines.append("> " + " ".join([exon.name, exon.toString(), exon.orientation]) +
                                          " score="+str(exon.data.split()[0]))
                    if exon.orientation == "+":
                        s = get_sequence(sequence=args.genome, ch=exon.chrom, 
                                         ss=exon.initial, es=exon.final, 
                                         strand=exon.orientation)
                    else:
                        s = get_sequence(sequence=args.genome, ch=exon.chrom, 
                                         ss=exon.initial, es=exon.final, 
                                         strand=exon.orientation, 
                                         reverse=True, complement=True)
                    ss = [s[i:i+70] for i in range(0, len(s), 70)]
                    writelines += ss

                with open(os.path.join(args.o, name + ".fa"), "w") as f:
                    for line in writelines:
                        print(line, file=f)
            else:

                s = get_sequence(sequence=args.genome, ch=r.chrom, ss=r.initial, es=r.final, 
                                 strand=r.orientation)
                ss = [s[i:i+70] for i in range(0, len(s), 70)]

                if ".fa" not in args.o:
                    with open(os.path.join(args.o, name + ".fa"), "w") as f:
                        if not r.orientation: r.orientation = "."
                        if not args.score:
                            print("> " + name + " " + r.toString() + " " + r.orientation, file=f)
                        else:
                            print("> " + name + " " + r.toString() + " " + r.orientation +
                                  " score=" + str(r.data.split()[0]), file=f)
                        for seq in ss: print(seq, file=f)
                elif ".fa" in args.o:
                    if not r.orientation: r.orientation = "."
                    if not args.score:
                        print("> " + name + " " + r.toString() + " " + r.orientation, file=f)
                    else:
                        print("> " + name + " " + r.toString() + " " + r.orientation +
                              " score=" + str(r.data.split()[0]), file=f)
                    for seq in ss: print(seq, file=fasta)
        if ".fa" in args.o:
            fasta.close()


    ############### BED filter gene #############################################
    elif args.mode == "bed_filter":
        print(tag + ": [BED] Filter genes")
        bed = GenomicRegionSet(args.i)
        bed.read(args.i)

        if args.gene:
            if os.path.isfile(args.gene):
                gg = GeneSet("genes")
                if args.score:
                    gg.read_expression(args.gene)
                else:
                    gg.read(args.gene)
                # print(len(gg))
                if not args.background:
                    bed = bed.by_names(gg, load_score=args.score)
                else:
                    bed = bed.by_names(gg, load_score=args.score, background=True)
                # print(len(bed))
            else:
                genes = [ x.upper() for x in args.gene.split(",") ]
                if not args.background:
                    bed = bed.by_names(genes, load_score=args.score)
                else:
                    bed = bed.by_names(genes, load_score=args.score, background=True)

        if isinstance(args.min, int ) and not args.max:
            bed = bed.filter_by_size(minimum=args.min)
        elif not args.min and isinstance(args.max, int ):
            bed = bed.filter_by_size(maximum=args.max)
        elif isinstance(args.min, int ) and isinstance(args.max, int ):
            bed = bed.filter_by_size(minimum=args.min, maximum=args.max)


        bed.write(args.o)
                    
        print("complete.")


    ############### BED remove if overlap ########################################
    elif args.mode == "bed_remove_if_overlap":
        print(tag + ": [BED] Remove the overlapping regions")
        
        if not args.t:
            print("Please define the file for target regions.")
            sys.exit(1)
        else:
            print("target:\t" + args.t)

        # with open(args.target) as f:
        t = GenomicRegionSet("targets")
        t.read(args.t, io=GRSFileIO.Bed12 if args.block else GRSFileIO.Bed)

        # with open(args.i) as fi, open(args.o, "w") as fo:
        input_regions = GenomicRegionSet("input")
        input_regions.read(args.i, io=GRSFileIO.Bed12 if args.block else GRSFileIO.Bed)
        if args.keep:
            output_regions = input_regions.intersect(t, mode=OverlapType.ORIGINAL)
        else:
            output_regions = input_regions.subtract(t, whole_region=True)

        output_regions.write(args.o, io=GRSFileIO.Bed12 if args.block else GRSFileIO.Bed)
        print("input regions:\t"+str(len(input_regions)))
        print("target regions:\t" + str(len(t)))
        print("output regions:\t" + str(len(output_regions)))
        print("complete.")

    ############### BED add columns #############################################
    elif args.mode == "bed_add_columns":
        print(tag + ": [BED] Add column")
        
        if not args.ref:
            print("Please define the file for reference.")
            sys.exit(1)

        with open(args.ref) as f:
            genes = {}
            for line in f:
                line = line.strip().split()
                # print(line[args.f-1].upper())
                # print(line[args.f:])
                try:
                    genes[line[args.f-1].upper()] = line[args.f:]
                except:
                    print("Error: indexing error. Please check -f argument.")
                    sys.exit(1)
        print(len(genes.keys()))

        with open(args.i) as fi, open(args.o, "w") as fo:
            c_add = 0
            c_miss = 0
            for line in fi:
                line = line.strip().split()
                try:
                    print("\t".join(line+genes[line[3].upper()]), file=fo)
                    c_add += 1
                except:
                    print("\t".join(line), file=fo)
                    c_miss += 1
        print("Modified genes:\t"+str(c_add))            
        print("Missed genes:\t"+str(c_miss))
        print("complete.")


    ############### BED divide by expression ###########################
    elif args.mode == "bed_divide":
        print(tag + ": [BED] Divide by expression data")
        
        gene1 = []
        gene2 = []
        
        with open(args.t) as t:
            for line in t:
                l = line.strip().split()
                g = l[0]
                try:
                    exp = [ float(x) for x in l[1:] ]
                except:
                    continue
                if args.m == "max":
                    if max(exp) > args.c: gene1.append(g)
                    else: gene2.append(g)
                elif args.m == "min":
                    if min(exp) > args.c: gene1.append(g)
                    else: gene2.append(g)
                elif args.m == "mean":
                    mean = sum(exp)/len(exp)
                    if mean > args.c: gene1.append(g)
                    else: gene2.append(g)

        bed = GenomicRegionSet(args.bed)
        bed.read(args.bed)
        o1 = GenomicRegionSet(args.o1)
        o2 = GenomicRegionSet(args.o2)
        for r in bed:
            if r.name in gene1: o1.add(r)
            elif r.name in gene2: o2.add(r)
        o1.write(args.o1)
        o2.write(args.o2)


    ############### BED average size ###########################
    elif args.mode == "bed_size":
        print(tag + ": [BED] Average size")
        bed = GenomicRegionSet("bed")
        bed.read(args.i)
        print("Average size:\t"+str(bed.average_size()))
        print("Size variance:\t" + str(bed.size_variance()))
        print()


    ############### BED complement ###########################
    elif args.mode == "bed_complement":
        print(tag + ": [BED] Get complementary regions from genome")
        bed = GenomicRegionSet("bed")
        bed.read(args.i)
        genome = GenomicRegionSet(args.organism)
        genome.get_genome_data(organism=args.organism, chrom_X=True, chrom_Y=False, chrom_M=False)
        res = genome.subtract(bed)
        res.write(args.o)
        print()

    ############### BED Detect polyA reads ###########################
    elif args.mode == "bed_polya":

        non_available = 0
        print(tag + ": [BED] Detect the reads with poly-A tail on the regions")
        def count_polyA_on_bam(bed, bam):
            pattern = "AAAAA"
            samfile = pysam.AlignmentFile(bam, "rb")
            win_width = int(numpy.mean([ read.qlen for read in samfile]))
            res = []
            for r in bed:
                count_polyA = 0
                all_read = 0
                if r.orientation == "-":
                    start = r.initial
                    end = r.initial + win_width
                else:
                    start = r.final - win_width
                    end = r.final
                for pileupcolumn in samfile.pileup(r.chrom, start, end):
                    for pileupread in pileupcolumn.pileups:
                        if not pileupread.is_del and not pileupread.is_refskip:
                            # query position is None if is_del or is_refskip is set.
                            all_read += 1
                            if pileupread.alignment.query_sequence.endswith(pattern):
                                # print(pileupread.alignment.query_sequence)
                                count_polyA += 1


                # all transcript
                rr = GenomicRegionSet(r.name)
                rr.add(r)
                rr.extract_blocks()
                all_a = 0
                all_r = 0
                for exon in rr:
                    for pileupcolumn in samfile.pileup(exon.chrom, exon.initial, exon.final):
                        for pileupread in pileupcolumn.pileups:
                            if not pileupread.is_del and not pileupread.is_refskip:
                                # query position is None if is_del or is_refskip is set.
                                all_r += 1
                                if pileupread.alignment.query_sequence.endswith(pattern):
                                    # print(pileupread.alignment.query_sequence)
                                    all_a += 1

                if all_read == 0:
                    if all_r == 0:
                        res.append([r.name, count_polyA, all_read, non_available, all_a, all_r, non_available])
                    else:
                        res.append([r.name, count_polyA, all_read, non_available, all_a, all_r, float(all_a)/all_r])
                else:
                    if all_r == 0:
                        res.append([r.name, count_polyA, all_read, float(count_polyA)/all_read, all_a, all_r, non_available])
                    else:
                        res.append([r.name, count_polyA, all_read, float(count_polyA)/all_read, all_a, all_r, float(all_a)/all_r])
            samfile.close()
            return res

        bed = GenomicRegionSet("bed")
        bed.read(args.i)

        if os.path.isfile(args.b):
            res = count_polyA_on_bam(bed=bed, bam=args.b)
            with open(args.o, "w") as f:
                print("\t".join(["name", "polyA_reads_in_window", "all_reads_in_window","proportion_in_window",
                                 "polyA_reads_on_transcript", "all_reads_on_transcript","proportion_on_transcript"]), file=f)
                for l in res:
                    print("\t".join([ str(x) for x in l ]), file=f)
        elif os.path.isdir(args.b):
            col_res = {}
            # bams = []
            for r in bed:
                col_res[r.name] = []

            for root, dirs, files in os.walk(args.b):
                for f in files:
                    if f.endswith(".bam"):
                        # bams.append(f.rpartition(".")[0])
                        res = count_polyA_on_bam(bed=bed, bam=os.path.join(root, f))
                        for line in res:
                            col_res[line[0]].append(line[1:])

            for r in bed:
                ar = numpy.array(col_res[r.name])
                # print(ar)
                # print(numpy.mean(ar, axis=1))
                col_res[r.name] = numpy.mean(ar, axis=1).tolist()
                # print(col_res[r.name])
                # sys.exit(1)

            with open(args.o, "w") as f:
                print("\t".join(["name", "polyA_reads_in_window_ave", "all_reads_in_window_ave",
                                 "polyA_reads_on_transcript_ave", "all_reads_on_transcript_ave"]), file=f)
                for gene, l in col_res.items():
                    print("\t".join([ gene ] + [ str(x) for x in l ]), file=f)
        print()


    ############### BED to GTF ###########################
    #
    elif args.mode == "bed_to_gtf":
        print(tag + ": [BED] Convert BED to GTF")
        inf = open(args.i, 'r')
        outf = open(args.o, 'w')

        for linea in inf:
            linea_split = linea.split()
            chrom = linea_split[0]
            ini_pos = int(linea_split[1])
            fin_pos = int(linea_split[2])
            peak = linea_split[3].split("_")[0]
            print(peak)

            outf.write(chrom + "\tjoseph\texon\t" + str(ini_pos) + "\t" + str(fin_pos) + '\t.\t+\t.\t'+
                       'gene_name "' + peak + '";' +
                       'gene_id "' + peak + '";' +
                       'transcript_name "' + peak + '";' +
                       'transcript_id "' + peak + '";' +
                       '\n')
        inf.close()
        outf.close()


    ############### BED overlaps ###########################
    #
    elif args.mode == "bed_overlap":
        beds = []
        bednames = []
        if os.path.isdir(args.i):
            for dirpath, dnames, fnames in os.walk(args.i):
                for f in fnames:
                    if f.endswith(".bed"):
                        name = os.path.basename(f).replace(".bed", "")
                        bed = GenomicRegionSet(name)
                        bed.read(os.path.join(dirpath, f))
                        bed.sort()
                        beds.append(bed)
                        bednames.append(name)

            index = natsort.index_natsorted(bednames)
            beds = natsort.order_by_index(beds, index)
            bednames = natsort.order_by_index(bednames, index)

        else:
            line = args.i.split(",")
            for b in line:
                if b.endswith(".bed"):
                    name = os.path.basename(b).replace(".bed", "")
                    bed = GenomicRegionSet(name)
                    bed.read(b)
                    bed.sort()
                    beds = [bed]
                    bednames = [name]

        count_dic = {}
        for bed1 in beds:
            count_dic[bed1.name] = {}
            for bed2 in beds:
                if bed1.name == bed2.name:
                    count_dic[bed1.name][bed2.name] = len(bed1)
                else:
                    inter = bed1.intersect(bed2, mode=OverlapType.ORIGINAL)
                    count_dic[bed1.name][bed2.name] = len(inter)

        out_table = [["Overlaps"]+bednames]
        for b1 in bednames:
            out_table.append([b1]+ [str(count_dic[b1][b2]) for b2 in bednames])

        if not args.o:
            for l in out_table:
                print("\t".join(l))
        else:
            with open(args.o, "w") as f:
                for l in out_table:
                    print("\t".join(l), file=f)




    ############### BED distance ###########################
    #
    elif args.mode == "bed_distance":
        print(tag + ": [BED] Calculate the distances between two region sets")
        bed = GenomicRegionSet(args.i)
        bed.read(args.i)
        target = GenomicRegionSet(args.t)
        target.read(args.t)

        res_dis = []
        for s in bed:
            dis = []
            for a in target:
                d = s.distance(a)
                if d != None:
                    dis.append(s.distance(a))
            res_dis.append([ s.name, str(min(dis)) ])

        # res_dis = bed.get_distance(target, ignore_overlap=False)
        with open(args.o, "w") as f:
            for line in res_dis:
                print("\t".join(line), file=f)


    ############### BED standardize chromsomes ###########################
    #
    elif args.mode == "bed_standard_chrom":
        print(tag + ": [BED] Standardize chomosomes")
        bed = GenomicRegionSet(args.i)
        bed.read(args.i)
        nbed = bed.standard_chrom()
        nbed.write(args.i)


    ############### BED add overlapping region name ###########################
    #
    elif args.mode == "bed_add_data":
        print(tag + ": [BED] Add overlapping region name")
        bed = GenomicRegionSet(args.i)
        bed.read(args.i)
        target = GenomicRegionSet(args.t)
        target.read(args.t)
        overlap_regions = target.intersect(bed, mode=OverlapType.ORIGINAL)
        with open(args.i) as fin:
            with open(args.o, "w") as fout:
                for line in fin:
                    if line.startswith("chr"):
                        line = line.strip()
                        l = line.split()
                        overlapping = overlap_regions.covered_by_aregion(GenomicRegion(chrom=l[0], initial=int(l[1]), final=int(l[2])))
                        if len(overlapping) > 0:
                            print(line + "\t" + ",".join([g.name for g in overlapping]), file=fout)
                        else:
                            print(line + "\t.", file=fout)


    ############### BED Sampling regions randomly ###########################
    #
    elif args.mode == "bed_sampling":
        print(tag + ": [BED] Sampling the regions randomly")
        bed = GenomicRegionSet(args.i)
        bed.read(args.i)
        beds = bed.random_subregions(size=args.s)
        beds.write(args.o)


    ############### BED bed12tobed6 ###########################
    #
    elif args.mode == "bed12tobed6":
        print(tag + ": [BED] Convert BED12 to BED6")
        bed = GenomicRegionSet(args.i)
        bed.read(args.i, io=GRSFileIO.Bed12)
        bed.write(args.o)


    ############### BAM filtering by BED ###########################
    #
    elif args.mode == "bam_filter":
        print(tag + ": [BED] Filtering BAM file by the regions in BED file")

        bed = GenomicRegionSet("bed")
        bed.read(args.bed)
        bam = pysam.AlignmentFile(args.i, "rb")
        outbam = pysam.AlignmentFile(args.o+".sam", "wh", template=bam)
        # in_reads = []
        for region in bed:
            if "_" not in region.chrom:
                for read in bam.fetch(region.chrom, region.initial, region.final):
                    if read.is_proper_pair:
                        if read.reference_start > region.initial and read.reference_end < region.final:
                            m = bam.mate(read)
                            if m.reference_start > region.initial and m.reference_end < region.final:
                                outbam.write(read)

        bam.close()
        outbam.close()
        os.system("samtools view -Sb "+args.o+".sam"+" > "+args.o+"_temp.bam")
        os.system("samtools sort " + args.o+"_temp.bam " + args.o)
        os.system("samtools index " + args.o + ".bam")
        pysam.index(args.o+".bam")
        os.remove(args.o+"_temp.bam")
        os.remove(args.o + ".sam")

    ############### THOR MAplot #############################################
    elif args.mode == "thor_ma":
        print(tag + ": [THOR] Generate MA plot")
        tag = os.path.basename(args.i).split(".")[0].split("/")[-1]
        if not os.path.exists(args.o):
            os.makedirs(args.o)

        chr = "chr10"
        # config_dir = os.path.dirname(args.i)
        # Parse data config
        data = {}
        with open(args.i) as f:
            t = None
            for line in f:
                l = line.strip()
                if l.startswith("#"):
                    data[l[1:]] = []
                    t = l[1:]
                elif t and l:
                    data[t].append(l)
        organism = data["chrom_sizes"][0].split(".")[0]
        # Get one chromosome
        chrom_size = None
        with open(data["chrom_sizes"][0]) as c:
            for line in c:
                l = line.strip().split()
                if l[0] == chr: chrom_size = int(l[1])
            if not chrom_size:
                print("Chromosome is not found in chrom.size file")
                sys.exit(1)
        gr_chrom = GenomicRegionSet(chr)
        gr_chrom.add(GenomicRegion(chrom=chr, initial=0,final=chrom_size))
        # Load factors
        info_file = glob.glob(os.path.join(args.thor, '*-setup.info'))
        with open(info_file[0]) as f:
            c = False
            ce = False
            for line in f:
                if line.startswith("#Scaling factors"):
                    c = True
                elif c:
                    factors = [ float(x) for x in line.strip().split() ]
                    c = False
                if line.startswith("#Extension size"):
                    ce = True
                elif ce:
                    extension = [ int(x) for x in line.strip().split() ]
                    ce = False
        # Calculate coverage rep1
        # rep = data["rep1"][0]
        bin_cov1 = []
        for i, rep in enumerate(data["rep1"]+data["rep2"]):
            cov = CoverageSet("rep", gr_chrom)
            cov.coverage_from_bam(bam_file=rep,extension_size=extension[i], binsize=int(args.b), stepsize=int(args.e))
            c = cov.coverage[0] + 1
            bin_cov1.append(c)

        # bin_cov1 = [ i.tolist() for i in bin_cov1 ]

        bin_cov2 = []
        for i, rep in enumerate(data["rep1"] + data["rep2"]):
            cov = CoverageSet("rep", gr_chrom)
            cov.coverage_from_bam(bam_file=rep, extension_size=extension[i], binsize=int(args.b), stepsize=int(args.e))
            # bin_cov2.append(cov.coverage[0]*factors[i])
            c = cov.coverage[0]*factors[i] + 1
            bin_cov2.append(c)
        # bin_cov1 = [i.tolist() for i in bin_cov1]


        bin_cov1 = [ numpy.log2(x) for x in bin_cov1]
        bin_cov2 = [ numpy.log2(x) for x in bin_cov2]
        # print(bin_cov1[0].max())
        s1_1_M = 0.5 * (bin_cov1[2] + bin_cov1[0])
        s1_1_A = bin_cov1[2] - bin_cov1[0]
        s1_2_M = 0.5 * (bin_cov2[2] + bin_cov2[0])
        s1_2_A = bin_cov2[2] - bin_cov2[0]
        s2_1_M = 0.5 * (bin_cov1[3] + bin_cov1[1])
        s2_1_A = bin_cov1[3] - bin_cov1[1]
        s2_2_M = 0.5 * (bin_cov2[3] + bin_cov2[1])
        s2_2_A = bin_cov2[3] - bin_cov2[1]
        max_x = max([s1_1_M.max(), s1_2_M.max(), s2_1_M.max(), s2_2_M.max()])
        max_y = max([numpy.absolute(s1_1_A).max(), numpy.absolute(s1_2_A).max(),
                     numpy.absolute(s2_1_A).max(), numpy.absolute(s2_2_A).max()])
        # print([max_x, max_y])

        plt.figure(1)
        plt.subplot(221)
        plt.title('Sample 1 before norm.')
        plt.scatter(s1_1_M.tolist(), s1_1_A.tolist(), s=1, alpha=0.5)
        plt.xlim([0, max_x])
        plt.ylim([-max_y, max_y])
        plt.ylabel('M (log ratio)')
        # plt.xlabel('A (mean average)')
        plt.subplot(222)
        plt.title('Sample 1 after norm.')
        plt.scatter(s1_2_M, s1_2_A, s=1, alpha=0.5)
        plt.xlim([0, max_x])
        plt.ylim([-max_y, max_y])
        plt.ylabel('M (log ratio)')
        # plt.xlabel('A (mean average)')
        plt.subplot(223)
        plt.title('Sample 2 before norm.')
        plt.scatter(s2_1_M, s2_1_A, s=1, alpha=0.5)
        plt.xlim([0, max_x])
        plt.ylim([-max_y, max_y])
        # plt.ylabel('M (log ratio)')
        plt.xlabel('A (mean average)')
        plt.subplot(224)
        plt.title('Sample 2 after norm.')
        plt.scatter(s2_2_M, s2_2_A, s=1, alpha=0.5)
        plt.xlim([0, max_x])
        plt.ylim([-max_y, max_y])
        # plt.ylabel('M (log ratio)')
        plt.xlabel('A (mean average)')
        # plt.savefig(pp, format='pdf')
        # pp.close()
        plt.savefig(os.path.join(args.o,tag+'_MAplot.png'), bbox_inches='tight')

        print("finish")

    ############### THOR split #############################################
    elif args.mode == "thor_split":
        print(tag + ": [THOR] Split the differential peaks")
        if not args.o:
            args.o = os.path.dirname(args.i)

        name = os.path.basename(args.i).split(".")[0]
        if args.fc == 0: tag = "_p" + str(args.p)
        else: tag = "_p"+str(args.p)+"_fc"+str(args.fc)

        bed = GenomicRegionSet("input")
        bed.read(args.i)
        print("Number of input peaks:\t"+str(len(bed)))

        if args.rename and args.g:
            bed2 = bed.gene_association(organism=args.g, strand_specific=True)
        else:
            bed2 = bed

        for region in bed2:
            data = region.data.split()
            stat = data[4].split(";")
            s1 = [float(x) + 1 for x in stat[0].split(":")]
            s2 = [float(x) + 1 for x in stat[1].split(":")]
            fc = math.log((sum(s2) / len(s2)) / (sum(s1) / len(s1)), 2)
            region.data = "\t".join([str(fc)] + data[1:])

        gain_peaks = GenomicRegionSet("gain_peaks")
        lose_peaks = GenomicRegionSet("lose_peaks")
        gain_table = GenomicRegionSet("gain_table")
        lose_table = GenomicRegionSet("lose_table")

        for region in bed2:
            l = region.data.split()
            s = l[4].split(";")
            if abs(float(l[0])) > args.fc and float(s[2]) > args.p:

                s1 = sum([int(x) for x in s[0].split(":")]) / len(s[0].split(":"))
                s2 = sum([int(x) for x in s[1].split(":")]) / len(s[1].split(":"))

                # print([len(region), args.s])
                if args.s:
                    nbins = int(len(region)/args.s)
                    ns1 = float(s1) / nbins
                    ns2 = float(s2) / nbins
                    data = "\t".join([l[0], str(s1), str(s2), str(len(region)),
                                      str(ns1), str(ns2), str(abs(ns1 + ns2)), str(abs(ns1 - ns2)), s[2]])
                else:
                    data = "\t".join([l[0], str(s1), str(s2), str(len(region)), s[2]])

                # Chromosome	Start	End	Name	FC	Strand	Ave. Count 1	Ave. Count 2
                # Length	Norm count 1	Norm count 2	Sum norm count	Diff norm count	P-value

                if float(l[0]) > 0:
                    gain_table.add(GenomicRegion(chrom=region.chrom, initial=region.initial, final=region.final,
                                                 orientation=region.orientation, data=data, name=region.name))
                    gain_peaks.add(region)
                elif float(l[0]) < 0:
                    lose_table.add(GenomicRegion(chrom=region.chrom, initial=region.initial, final=region.final,
                                                 orientation=region.orientation, data=data, name=region.name))
                    lose_peaks.add(region)

        gain_peaks.write(os.path.join(args.o, name + tag + "_gain.bed"))
        lose_peaks.write(os.path.join(args.o, name + tag + "_lose.bed"))
        gain_table.write(os.path.join(args.o, name + tag + "_gain.table"))
        lose_table.write(os.path.join(args.o, name + tag + "_lose.table"))

        print("Number of gain peaks:\t" + str(len(gain_peaks)))
        print("Number of lose peaks:\t" + str(len(lose_peaks)))
        
        
    ############### getseq #############################################
    elif args.mode == "getseq":
        print(tag + ": [FASTA] Get sequence from the given regions")
        
        if args.b:
            regions = GenomicRegionSet("regions")
            regions.read(args.b)
            for r in regions:
                print(r.name)
                s = get_sequence(ch=r.chrom, ss=r.initial, es=r.final, strand=r.orientation, 
                                 rna=args.r, ex=args.ex)
                print(s[:20])
                ss = [s[i:i+70] for i in range(0, len(s), 70)]

                with open(r.name + ".fa", "w") as f:
                    print("> "+ r.name + " "+r.toString()+ " "+r.orientation, file=f)
                    for seq in ss:
                        print(seq, file=f)

        else:
            if args.p:
                # print(args.p)
                if ":" not in args.p:
                    args.ch = None
                    args.ss = int(args.p.split("-")[0])
                    args.es = int(args.p.split("-")[1])
                else:
                # print(args.p.partition(":")[2])
                    args.ch = args.p.partition(":")[0]
                    args.ss = int(args.p.partition(":")[2].split("-")[0])
                    args.es = int(args.p.partition(":")[2].split("-")[1])
                    print([args.ch, args.ss, args.es])
            if args.s == "both":
                seq1 = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand="+",
                                    reverse=args.reverse, complement=False, rna=args.rna, ex=args.ex)
                seq2 = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand="-",
                                    reverse=args.reverse, complement=True, rna=args.rna, ex=args.ex)
                print("5'- " + seq1 + " -3'")
                print("3'- " + seq2 + " -5'")
            elif args.s == "+" and not args.reverse:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.reverse, complement=args.complement, rna=args.rna, ex=args.ex)
                print("5'- " + seq + " -3'")
            elif args.s == "+" and args.reverse:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.reverse, complement=args.complement, rna=args.rna, ex=args.ex)
                print("3'- " + seq + " -5'")
            elif args.s == "-" and not args.reverse:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.reverse, complement=args.complement, rna=args.rna, ex=args.ex)
                print("3'- "+seq+" -5'")
            elif args.s == "-" and args.reverse:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.reverse, complement=args.complement, rna=args.rna, ex=args.ex)
                print("5'- " + seq + " -3'")
            
        print()
    ############### WIG trim end #############################################
    elif args.mode == "wig_trim_end":
        
        chromD = {}
        with open(args.chrosize) as c:
            for line in c:
                line = line.strip().split()
                if len(line) == 2:
                    chromD[line[0]] = int(line[1])

        with open(args.i) as f:
            with open(args.o, "w") as g:
                for line in f:
                    if "chrom=" in line and "span=" in line:
                        line = line.strip()
                        l = line.split()
                        chrom = l[1].partition("=")[2]
                        step = int(l[2].partition("=")[2])
                        print(line, file=g)
                    elif line.startswith("track"):
                        print(line.strip(), file=g)
                    else:
                        line = line.strip()
                        l = line.split()
                        if int(l[0]) + step < chromD[chrom]:
                            print(line, file=g)
                        else:
                            pass


    ############### GENE  ensembl2symbol #############################################
    elif args.mode == "ensembl2symbol":

        g = GeneSet("ensembl_id")
        g.read(args.i)

        ann = AnnotationSet(gene_source=args.organism, tf_source=None, alias_source=args.organism, 
                            filter_havana=False, protein_coding=False, known_only=False)
        mapped_list, unmapped_list = ann.get_official_symbol(gene_name_source=g.genes)
        print("\t"+str(len(g.genes))+"\tgenes are loaded.")
        print("\t"+str(len(mapped_list))+"\tgenes are mapped.")
        print("\t"+str(len(unmapped_list))+"\tgenes are not mapped.")
        g.genes = mapped_list
        g.save(args.o)


    elif args.mode == "symbol2ensembl":
        g = GeneSet("symbol")
        g.read(args.i)
        print(len(g.genes))

        ann = AnnotationSet(gene_source=args.organism, tf_source=None, alias_source=args.organism,
                            filter_havana=False, protein_coding=False, known_only=False)
        mapped_list, unmapped_list = ann.fix_gene_names(g)
        print("\t"+str(len(g.genes))+"\tgenes are loaded.")
        print("\t"+str(len(mapped_list))+"\tgenes are mapped.")
        print("\t"+str(len(unmapped_list))+"\tgenes are not mapped.")
        g.genes = mapped_list
        print(len(g.genes))
        g.save(args.o)


    elif args.mode == "fasta2bp":

        def fasta2bp(filename):
            s = ""
            exon = 0
            with open(filename) as f:
                for line in f:
                    l = line.strip()
                    if l.startswith(">"):
                        exon += 1
                    elif not l: continue
                    else: s += l
            return [str(len(s)), str(exon)]

        if os.path.isfile(args.i):
            l = fasta2bp(args.i)
            print("Length:\t\t" + l[0]+" bp; Number of exons:\t\t"+ l[1])

        elif os.path.isdir(args.i) and args.o:
            list_bp = []
            for root, dirs, files in os.walk(args.i):
                for f in files:
                    if f.endswith(".fa") or f.endswith(".fasta"):
                        # print(f.partition(".")[0])
                        # print(parser_fasta2bp(filename=os.path.join(root,f)))
                        list_bp.append([f.partition(".")[0]] + fasta2bp(os.path.join(root,f)))
        if args.o:
            with open(args.o, "w") as g:
                print("\t".join(["name", "sequence_length","number_of_exons"]), file=g)
                for l in list_bp:
                    print("\t".join(l), file=g)

    ############### STAR junction to BED #######################################
    elif args.mode == "circRNA":
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        
        with open(args.i) as f:
            with open(args.o, "w") as g:
                for line in f:
                    line = line.strip().split()
                    #  same chromosome        same strand
                    if line[0] == line[3] and line[2] == line[5] and int(line[6]) >= 0:
                        donor = [int(line[1]), int(line[10])]
                        acceptor = [int(line[4]), int(line[12])]
                        # Limit the distance
                        if abs(int(donor[0]) - int(acceptor[0])) < 20000:
                            
                            if line[2] == "+" and min(donor) > min(acceptor):
                                #acceptor_end = str(min([int(acceptor[0]),int(acceptor[1])]))
                                #donor_end = str(max([int(donor[0]),int(donor[1])]))
                                print("\t".join([ line[0],
                                                  str(min(acceptor)), str(max(donor)), line[9],
                                                  line[6], line[2], str(min(acceptor)), str(max(donor)),
                                                  "255,0,0", "2", 
                                                  ",".join([ str(abs(acceptor[1] - acceptor[0])),
                                                             str(abs(donor[1] - donor[0])) ]),
                                                  "0,"+str(abs(min(donor)-min(acceptor))) ]), file=g)

                            elif line[2] == "-" and min(acceptor) > min(donor):
                                #acceptor_end = str(max([int(acceptor[0]),int(acceptor[1])]))
                                #donor_end = str(min([int(donor[0]),int(donor[1])]))
                                print("\t".join([ line[0],
                                                  str(min(donor)), str(max(acceptor)), line[9],
                                                  line[6], line[2], str(min(donor)), str(max(acceptor)), 
                                                  "255,0,0", "2", 
                                                  ",".join([ str(abs(donor[1] - donor[0])),
                                                             str(abs(acceptor[1] - acceptor[0])) ]),
                                                  "0,"+str(abs(min(donor)-min(acceptor))) ]), file=g)
                        else:
                            pass
                            #print(line)
                            #sys.exit()

        print("tcons:\t" + args.t)
        tcons = GenomicRegionSet("tcons")
        tcons.read(args.t)
        circrna = GenomicRegionSet("circRNA")
        circrna.read(args.o)
        circ_inTCON = circrna.intersect(y=tcons, mode = OverlapType.COMP_INCL)
        circ_TCONs = tcons.intersect(y=circ_inTCON, mode = OverlapType.ORIGINAL)
        #print(len(circ_TCONs))
        circ_TCONs.write(args.c)


#  0        1       2     3         4       5   6   7   8
# chr1  39449029    +   chr1    39448068    +   0   0   0
#                     9                          10          11      12            13   
# 97ZZTR1:411:C4VC3ACXX:5:1102:16097:34171    39448994    35M15S  39448069    35S15M868p50M

    ############### FASTA slicing #######################################
    elif args.mode == "sliceFASTA":
        print(os.path.basename(args.i) + " -start "+str(args.p)+" -end "+str(args.p+args.l))
        from rgt.SequenceSet import SequenceSet
        seq = SequenceSet(name=args.i, seq_type="RNA")
        seq.read_fasta(fasta_file=args.i)
        start = int(args.p)
        end = int(start + args.l)
        if args.r:
            print("3' - "+ seq.sequences[0].seq[end:start:-1]+ " - 5'")
        else:
            print("5' - "+ seq.sequences[0].seq[start:end]+ " - 3'")



    ############### FASTA slicing #######################################
    elif args.mode == "txp2bed":
        from rgt.tdf.RNADNABindingSet import RNADNABindingSet
        txp = RNADNABindingSet("txp")
        txp.read_txp(filename=args.i, dna_fine_posi=True)
        tmp = os.path.join(os.path.dirname(args.o), "temp.bed")
        txp.write(filename=tmp)
        os.system("sort -k1,1V -k2,2n " + tmp + " > " + args.o)
        os.remove(tmp)
