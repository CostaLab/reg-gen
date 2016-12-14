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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from os.path import expanduser
home = expanduser("~")

# Local Libraries
from rgt.GeneSet import GeneSet
from rgt.CoverageSet import CoverageSet
from rgt.GenomicRegion import GenomicRegion
from rgt.AnnotationSet import AnnotationSet
from rgt.Util import OverlapType, GenomeData
from rgt.GenomicRegionSet import GenomicRegionSet
tag = "RGT-tools"

#print(os.getcwd())

def get_sequence(sequence, ch, ss, es, reverse=False, complement=False, rna=False, ex=0, strand=None):
    import pysam
    sequence = pysam.Fastafile(sequence)
    seq = sequence.fetch(ch, max(0, ss-ex), es+ex )
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
    
    subparsers = parser.add_subparsers(help='sub-command help',dest='mode')

    ############### GTF add transcripts ######################################
    parser_gtfat = subparsers.add_parser('gtf_add_transcripts', 
                                         help="[GTF] Add transcripts from the existed exons")
    parser_gtfat.add_argument('-i', metavar='  ', type=str, help="Input GTF file")
    parser_gtfat.add_argument('-o', metavar='  ', type=str, help="Output GTF file")


    ############### GTF to BED ###############################################
    # python rgt-convertor.py gtf_to_bed -i -o
    parser_gtf2bed = subparsers.add_parser('gtf_to_bed',
                                           help="[GTF] Convert GTF file to BED by the given biotype")
    parser_gtf2bed.add_argument('-i', metavar='  ', type=str, help="Input GTF file")
    parser_gtf2bed.add_argument('-o', metavar='  ', type=str, help="Output BED file")
    parser_gtf2bed.add_argument('-s', '-source', type=str, default="All", help="Define the source {ENSEMBL,HAVANA,All}")
    parser_gtf2bed.add_argument('-f', '-feature', type=str, default="gene",
                                help="Define the feature {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine,All}")
    parser_gtf2bed.add_argument('-t', '-type', type=str, default="All",
                                help="Define gene type e.g. 'protein_coding' more: http://www.gencodegenes.org/gencode_biotypes.html")
    parser_gtf2bed.add_argument('-st', '-status', type=str, default="All",
                                help="Define gene status {KNOWN, NOVEL, PUTATIVE,All}")
    parser_gtf2bed.add_argument('-g', '-gene', type=str, default=None,
                                help="Define the gene list for filtering, default is None.")
    parser_gtf2bed.add_argument('-id', action="store_true",
                                help="Use gene ID as region name, instead of gene symbol.")
    parser_gtf2bed.add_argument('-b', action="store_true",
                                help="Save exons into entries with block in BED")

    ############### GTF to FASTA #############################################
    # python rgt-convertor.py
    parser_gtf2fasta = subparsers.add_parser('gtf_to_fasta', 
                                             help="[GTF] Convert GTF file to FASTA (exons) by the given gene name")
    parser_gtf2fasta.add_argument('-i', metavar='  ', type=str, help="Input GTF file")
    parser_gtf2fasta.add_argument('-o', metavar='  ', type=str, help="Output FASTA file")
    parser_gtf2fasta.add_argument('-t', metavar='  ', type=str, help="Define the target transcript")
    parser_gtf2fasta.add_argument('-g', metavar='  ', type=str, help="Define the target gene")
    parser_gtf2fasta.add_argument('-genome', type=str, help="Define the FASTA file of the genome")
    
    ############### GTF add chr on each entry #################################
    # python rgt-convertor.py
    parser_gtfachr = subparsers.add_parser('gtf_add_chr', 
                                             help="[GTF] Add 'chr' to each line in GTF for proper chromosome name")
    parser_gtfachr.add_argument('-i', metavar='  ', type=str, help="Input GTF file")
    
    ############### BED add score ############################################
    # python rgt-convertor.py
    parser_bedac = subparsers.add_parser('bed_add_score', help="[BED] Add score column")
    parser_bedac.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedac.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedac.add_argument('-v', type=str, help="Define value to add")

    ############### BED merge by name ############################################
    # python rgt-convertor.py
    parser_bedmn = subparsers.add_parser('bed_merge_by_name', help="[BED] Merge regions by name")
    parser_bedmn.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedmn.add_argument('-o', '-output', type=str, help="Output BED file")

    ############### BED rename ###############################################
    # python rgt-convertor.py
    parser_bedrename = subparsers.add_parser('bed_rename', help="[BED] Rename regions by associated genes")
    parser_bedrename.add_argument('-i', metavar='  ', type=str, help="Input BED file")
    parser_bedrename.add_argument('-o', metavar='  ', type=str, help="Output BED file")
    parser_bedrename.add_argument('-d', action="store_true", help="Show the distance")
    parser_bedrename.add_argument('-organism',metavar='  ', type=str, help="Define the organism")
    parser_bedrename.add_argument('-l', metavar='  ', type=int, default=1000, 
                                  help="Define the length of promoter region (default:1000 bp)")
    parser_bedrename.add_argument('-t', metavar='  ', type=int, default=50000, 
                                  help="Define the threshold of distance (default:50000bp")
    parser_bedrename.add_argument('-target', metavar='  ', default=False, type=str, help="Target BED file")
    ############### BED change strand ###############################################
    # python rgt-convertor.py
    parser_bedchstrand = subparsers.add_parser('bed_change_strand', help="[BED] Change strand of regions by the target BED file")
    parser_bedchstrand.add_argument('-i', metavar='  ', type=str, help="Input BED file")
    parser_bedchstrand.add_argument('-o', metavar='  ', type=str, help="Output BED file")
    parser_bedchstrand.add_argument('-d', metavar='  ', type=int, default=0,
                                    help="Define the threshold of distance (default:0 bp")
    parser_bedchstrand.add_argument('-t', metavar='  ', type=str, help="Target BED file")
    parser_bedchstrand.add_argument('-r', action="store_true", help="Reverse the strand")

    ############### BED extend ###############################################
    # python rgt-convertor.py
    parser_bedex = subparsers.add_parser('bed_extend', help="[BED] Extend the regions")
    parser_bedex.add_argument('-i', type=str, help="Input BED file")
    parser_bedex.add_argument('-o', type=str, help="Output BED name.")
    parser_bedex.add_argument('-oz', "-onlyzero", action="store_true", default=False, 
                              help="Extend only the zero-length regions")
    parser_bedex.add_argument('-l', type=int, help="Define the length to extend.")
    parser_bedex.add_argument('-both',action="store_true", default=False, 
                              help="Extend from the both ends.")

    ############### BED subtract ###############################################
    # python rgt-convertor.py
    parser_bedsub = subparsers.add_parser('bed_subtract', help="[BED] Subtract the regions")
    parser_bedsub.add_argument('-i', type=str, help="Input BED file")
    parser_bedsub.add_argument('-o', type=str, help="Output BED name.")
    parser_bedsub.add_argument('-t', "-target", type=str,
                              help="Define the target BED file to subtract.")

    ############### BED cut ###############################################
    # python rgt-convertor.py
    parser_bedcut = subparsers.add_parser('bed_cut', help="[BED] Cut the regions")
    parser_bedcut.add_argument('-i', type=str, help="Input BED file")
    parser_bedcut.add_argument('-o', type=str, help="Output BED name.")
    parser_bedcut.add_argument('-t', "-target", type=str,
                               help="Define the target BED file for cutting.")

    ############### BED get promoters ########################################
    # python rgt-convertor.py bed_get_promoters -i -o -organism
    parser_bedgp = subparsers.add_parser('bed_get_promoters', 
                       help="[BED] Get promoters from the given genes")
    parser_bedgp.add_argument('-i', '-input', type=str, help="Input file (BED or gene list)")
    parser_bedgp.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedgp.add_argument('-organism', type=str, help="Define the organism (necessary if input is a gene list)")
    parser_bedgp.add_argument('-l', type=int, default=1000, 
                              help="Define length of promoters (default:1000bp)")
    
    ############### BED get upstream regions #################################
    # python rgt-convertor.py bed_upstream -i -o
    parser_bedupstream = subparsers.add_parser('bed_upstream', 
                       help="[BED] Get regions upstream from the given BED file")
    parser_bedupstream.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedupstream.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedupstream.add_argument('-l', type=int, default=100, help="Define length (default:100bp)")
    parser_bedupstream.add_argument('-d', type=int, default=100, help="Define distance (default:100bp)")
    parser_bedupstream.add_argument('-min', type=int, default=0, 
                                    help="Define minimum length of gene to filter out the small genes (default:0)")
    parser_bedupstream.add_argument('-r', action="store_true", default=False, help="Reverse the strand.")

    ############### BED to FASTA #############################################
    parser_bed2fasta = subparsers.add_parser('bed_to_fasta', 
                       help="[BED] Export the sequences in FASTA according to the given BED file")
    parser_bed2fasta.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bed2fasta.add_argument('-o', '-output', type=str, help="Output directory for FASTA files")
    parser_bed2fasta.add_argument('-genome', type=str, help="Define the FASTA file of the genome sequence")

    ############### BED filtered by gene name ################################
    parser_bed2fasta = subparsers.add_parser('bed_filter_gene', 
                       help="[BED] Filter by the given gene list")
    parser_bed2fasta.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bed2fasta.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bed2fasta.add_argument('-gene', type=str, help="Define file for the gene list")

    ############### BED remove if overlap ################################
    parser_bedro = subparsers.add_parser('bed_remove_if_overlap', 
                       help="[BED] Remove the regions if they overlap with the target regions")
    parser_bedro.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedro.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedro.add_argument('-t', '-target', type=str, help="Define BED file for target regions")
    parser_bedro.add_argument('-k', action="store_true", default=False, help="Keep the overlapped regions, and remove the non-overlapped ones.")

    ############### BED add columns ################################
    parser_bedaddcol = subparsers.add_parser('bed_add_columns', 
                       help="[BED] Add extra columns to the BED file by gene name")
    parser_bedaddcol.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedaddcol.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedaddcol.add_argument('-ref', type=str, help="Define file for referring the extra columns ")
    parser_bedaddcol.add_argument('-f', '-field', type=int, help="Which field of the reference file is compared for names.")

    ############### Divide regions in BED by expression #######################
    # python rgt-convertor.py divideBED -bed -t -o1 -o1 -c -m
    parser_divideBED = subparsers.add_parser('bed_divide', 
                       help="[BED] Divide the BEd files by the expression.")
    parser_divideBED.add_argument('-bed', type=str, help="Input BED file")
    parser_divideBED.add_argument('-t','-table', type=str, help="Input expression table (Gene name should match the region name.")
    parser_divideBED.add_argument('-o1', '-output1', type=str, help="Output first BED file")
    parser_divideBED.add_argument('-o2', '-output2', type=str, help="Output second BED file")
    parser_divideBED.add_argument('-c', '-cutoff', type=int, help="Define the cutoff")
    parser_divideBED.add_argument('-m', type=str, help="Define the mode, such as mean, max, or min.")

    ############### BAM Filter reads by BED file #######################
    parser_filterBAM = subparsers.add_parser('bam_filter',
                                             help="[BAM] Filter BAM file by the given regions in BED.")
    parser_filterBAM.add_argument('-i', type=str, help="Input BAM file")
    parser_filterBAM.add_argument('-bed', type=str, help="Input BED file for the regions for filtering")
    parser_filterBAM.add_argument('-o', type=str, help="Output prefix for BAM file")

    ############### THOR MAplot ################################
    parser_thorma = subparsers.add_parser('thor_ma',
                       help="[THOR] Create the MA plot for understanding the effect of normalization.")
    parser_thorma.add_argument('-i', '-input', type=str, help="Input data config.")
    parser_thorma.add_argument('-thor', '-thor--result', type=str, default=".", help="Output directory of THOR.")
    parser_thorma.add_argument('-o', '-output', default=None, type=str, help="Output directory")
    parser_thorma.add_argument('-e', '-ext', default=None, type=str, help="Define the extension size.")
    parser_thorma.add_argument('-b', '-bin', default=None, type=str, help="Define the bin size.")

    ############### THOR split and filter ################################
    parser_thorsf = subparsers.add_parser('thor_split', 
                       help="[THOR] Split and filter the differential peaks from rgt-THOR")
    parser_thorsf.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_thorsf.add_argument('-o', '-output', default=None, type=str, help="Output directory.")
    parser_thorsf.add_argument('-p', '-p--value', type=int, help="Define the cut-off of p-value (-log10) for filtering.")
    parser_thorsf.add_argument('-fc', '-fold--change', type=int,default=0,
                               help="Define the cut-off of foldchange for filtering.")
    parser_thorsf.add_argument('-rn', '-rename', action="store_true",
                               help="Rename the peak names by associated genes.")
    parser_thorsf.add_argument('-g', '-genome', type=str, help="Define the genome")

    ############### GENOME get sequence ####################################################
    parser_getseq = subparsers.add_parser('getseq', 
                       help="[FASTA] Get sequence from genome FASTA")
    parser_getseq.add_argument('-o', type=str, help="Output FASTA file")
    parser_getseq.add_argument('-d', '-dna', type=str, help="DNA sequence in FASTA format")
    parser_getseq.add_argument('-b', '-bed', type=str, help="Input BED file")
    parser_getseq.add_argument('-p', '-pos', type=str, help="position")
    parser_getseq.add_argument('-s', '-strand', type=str, default="both", help="strand (+, -, or both)")
    parser_getseq.add_argument('-ch', '-chr', type=str, help="chromosome")
    parser_getseq.add_argument('-ss', '-start', type=int, help="start site")
    parser_getseq.add_argument('-es', '-end', type=int, help="end site")
    parser_getseq.add_argument('-ex', '-extention', default=0, type=int, help="extention")
    parser_getseq.add_argument('-re', '-reverse', action="store_true", help="Reverse the sequence")
    parser_getseq.add_argument('-c', '-complement', action="store_true", help="Get the complement of the sequence")
    parser_getseq.add_argument('-r', '-rna', default=False, action="store_true", help="Convert T to U")

    ############### WIG trim ends by chromosome #############################################
    parser_wig_trim = subparsers.add_parser('wig_trim_end', 
                       help="[WIG] Trim the WIG file according to the given chromosome size")
    parser_wig_trim.add_argument('-i', '-input', type=str, help="Input WIG file")
    parser_wig_trim.add_argument('-o', '-output', type=str, help="Output WIG file")
    parser_wig_trim.add_argument('-chrosize', type=str, help="Define path to the chromosome size file")

    ############### GENE list convertion #############################################
    parser_ensembl2symbol = subparsers.add_parser('ensembl2symbol', 
                       help="[GENE] Convert the gene list from ensembl ID to gene symbol")
    parser_ensembl2symbol.add_argument('-i', '-input', type=str, help="Input gene list")
    parser_ensembl2symbol.add_argument('-o', '-output', type=str, help="Output gene list")
    parser_ensembl2symbol.add_argument('-organism', type=str, help="Define the organism")

    parser_sumbol2ensembl = subparsers.add_parser('symbol2ensembl',
                                                  help="[GENE] Convert the gene list from gene symbol to ensembl ID")
    parser_sumbol2ensembl.add_argument('-i', '-input', type=str, help="Input gene list")
    parser_sumbol2ensembl.add_argument('-o', '-output', type=str, help="Output gene list")
    parser_sumbol2ensembl.add_argument('-organism', type=str, help="Define the organism")

    ############### Get length in bp from FASTA files ################################
    parser_fasta2bp = subparsers.add_parser('fasta2bp',
                                            help="[FASTA] Get the length in bp of FASTA files")
    parser_fasta2bp.add_argument('-i', '-input', type=str, help="Input FASTA file or directory")
    parser_fasta2bp.add_argument('-o', '-output', type=str, default="", help="Output file with a table")

    ############### STAR junction to BED #############################################
    parser_circRNA = subparsers.add_parser('circRNA', 
                       help="[junction] Convert the Chimeric junction from STAR to BED file")
    parser_circRNA.add_argument('-i', '-input', type=str, help="Input chimeric junction file from STAR")
    parser_circRNA.add_argument('-t', '-tcons', type=str, help="Input BED file of tcons")
    parser_circRNA.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_circRNA.add_argument('-c', '-circ', type=str, help="Output BED file of circular RNA")

    ############### FASTA slicing #############################################
    # python rgt-convertor.py sliceFASTA -i -o -l -p
    parser_sliceFASTA = subparsers.add_parser('sliceFASTA', 
                       help="[FASTA] Slice the sequence by given position and length")
    parser_sliceFASTA.add_argument('-i', '-input', type=str, help="Input FASTA file")
    parser_sliceFASTA.add_argument('-l', type=int, help="Length of the slice sequence")
    parser_sliceFASTA.add_argument('-o', '-output', type=str, help="Output FASTA file")
    parser_sliceFASTA.add_argument('-p', type=int, help="The start position")
    parser_sliceFASTA.add_argument('-r', default=False, action="store_true", help="Reverse the sequence")




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
        #print(entries[0][8])
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
                gs = [s for s in info if "gene_name" in s][0].partition("\"")[2][:-1].partition(" (")[0]

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

        if args.b:
            # exons = GenomicRegionSet("output")
            # exons.read_bed(args.o)
            # exons.write_bed_blocks(args.o)
            regions.write_bed_blocks(filename=args.o)
        else:
            regions.write_bed(filename=args.o)
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
        print("input:\t" + args.i)
        with open("temp.gtf", "w") as f:
            with open(args.i) as g:
                for line in g:
                    if line.startswith("#"):
                        print(line, file=f)
                    else:
                        print("chr"+line, file=f)
        # rewrite gtf file
        #os.move("temp.gtf")
        

    ############### BED add score ############################################
    elif args.mode == "bed_add_score":
        print(tag+": [BED] Add scores")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        with open(args.i) as f, open(args.o, "w") as g:
            for line in f:
                line = line.strip() 
                print(line+"\t"+args.v, file=g)


    ############### BED merge by name ########################################
    elif args.mode == "bed_merge_by_name":
        print(tag+": [BED] Merge regions by name")
        print("input:\t" + args.i)
        print("output:\t" + args.o)

        bed1 = GenomicRegionSet("input")
        bed1.read_bed(args.i)
        bed2 = bed1.mergebyname()
        bed2.write_bed(args.o)

    ############### BED rename regions #######################################
    elif args.mode == "bed_rename":
        print(tag+": [BED] Rename regions by associated genes")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        if args.target:
            print("target:\t" + args.target)
        else:
            print("organism:\t" + args.organism)

        bed = GenomicRegionSet(args.i)
        bed.read_bed(args.i)
        if args.target:
            target = GenomicRegionSet(args.target)
            target.read_bed(args.target)
            bed.replace_region_name(regions=target)
            bed.write_bed(args.o)
        else:
            renamebed = bed.gene_association(gene_set=None, organism=args.organism,
                                             promoterLength=args.l,
                                             threshDist=args.t, show_dis=args.d)
            renamebed.write_bed(args.o)


    ############### BED change strands #######################################
    elif args.mode == "bed_change_strand":
        print(tag + ": [BED] Change strands by target BED file")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("target:\t" + args.t)

        bed = GenomicRegionSet(args.i)
        bed.read_bed(args.i)
        target = GenomicRegionSet(args.t)
        target.read_bed(args.t)
        if args.d != "0":
            target.extend(left=int(args.d), right=int(args.d))
        bed.replace_region_strand(regions=target, reverse=args.r)
        bed.write_bed(args.o)


    ############### BED extend ###############################################
    elif args.mode == "bed_extend":
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        bed = GenomicRegionSet("bed")
        bed.read_bed(args.i)
        for region in bed:
            if args.oz:
                if region.initial == region.final:
                    region.final += args.l
            else:
                if args.both:
                    region.initial -= args.l
                else: pass
                region.final += args.l

        bed.write_bed(args.o)


    ############### BED subtract ###############################################
    elif args.mode == "bed_subtract":
        print("input:\t" + args.i)
        print("target:\t" + args.t)
        print("output:\t" + args.o)
        bed = GenomicRegionSet("bed")
        bed.read_bed(args.i)
        target = GenomicRegionSet("target")
        target.read_bed(args.t)
        out = bed.subtract(y=target)
        out.write_bed(args.o)


    ############### BED subtract ###############################################
    elif args.mode == "bed_cut":
        print("input:\t" + args.i)
        print("target:\t" + args.t)
        print("output:\t" + args.o)
        bed = GenomicRegionSet("bed")
        bed.read_bed(args.i)
        target = GenomicRegionSet("target")
        target.read_bed(args.t)
        out = GenomicRegionSet("output")
        for b in bed:
            z = GenomicRegionSet("temp")
            z.add(b)
            y = z.subtract(y=target)
            y.sort()
            if len(y) > 0:
                if b.orientation == "+":
                    out.add(y[-1])
                else:
                    out.add(y[0])
        out.write_bed(args.o)

    ############### BED get promoters #########################################
    elif args.mode == "bed_get_promoters":
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("organism:\t" + args.organism)
        gene = GenomicRegionSet("genes")
        ### Input BED file
        if args.i.endswith(".bed"):
            gene.read_bed(args.i)
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
        promoter.write_bed(args.o)



    ############### BED get upstream regions ####################################
    elif args.mode == "bed_upstream":
        print("input:\t" + args.i)
        print("output:\t" + args.o)

        gene = GenomicRegionSet("genes")
        ### Input BED file
        
        gene.read_bed(args.i)
        print(len(gene))
        target = GenomicRegionSet("target")
        # if args.min == 0: cut = float("inf")
        # elif args.min > 0: cut = args.min

        for s in gene:
            if s.orientation == "+" and len(s) > args.min: 
                s.initial, s.final = max(s.initial-args.d-args.l, 0), max(s.initial-args.d, 0)
                if s.initial > s.final: s.initial, s.final = s.final, s.initial
                if args.r: s.orientation = "-"
                target.add(s)

            elif s.orientation == "-" and len(s) > args.min: 
                s.initial, s.final = s.final+args.d, s.final+args.d+args.l
                if s.initial > s.final: s.initial, s.final = s.final, s.initial
                if args.r: s.orientation = "+"
                target.add(s)
        
        print(len(target))
        target.write_bed(args.o)


    ############### BED to FASTA #############################################
    elif args.mode == "bed_to_fasta":
        print(tag+": [BED] BED to FASTA")
        print("input:\t\t" + args.i)
        print("output directory:\t" + args.o)
        if not os.path.exists(args.o): os.makedirs(args.o)
        regions = GenomicRegionSet("regions")
        regions.read_bed(args.i)
        for r in regions:
            if len(r.data.split()) == 7:
                target = r.extract_blocks()
                #print("*** using block information in BED file")
                writelines = []
                for exon in target:
                    writelines.append("> "+" ".join([exon.name, exon.toString(), exon.orientation]))
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

                with open(os.path.join(args.o, r.name + ".fa"), "w") as f:
                    for line in writelines:
                        print(line, file=f)
            else:

                s = get_sequence(sequence=args.genome, ch=r.chrom, ss=r.initial, es=r.final, 
                                 strand=r.orientation)
                ss = [s[i:i+70] for i in range(0, len(s), 70)]
                with open(os.path.join(args.o, r.name + ".fa"), "w") as f:
                    print("> "+ r.name + " "+r.toString()+ " "+r.orientation, file=f)
                    for seq in ss:
                        print(seq, file=f)

    ############### BED filter gene #############################################
    elif args.mode == "bed_filter_gene":
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        
        if not args.gene:
            print("Please define the file for gene list.")
            sys.exit(1)

        with open(args.gene) as f:
            genes = f.read().splitlines()
            genes = map(lambda x: x.split("\t")[0].upper(), genes)
            print(str(len(genes))+" genes are loaded.")

        with open(args.i) as fi, open(args.o, "w") as fo:
            for line in fi:
                line = line.strip().split()
                if line[3].upper() in genes:
                    print("\t".join(line), file=fo)
                    
        print("complete.")


    ############### BED remove if overlap ########################################
    elif args.mode == "bed_remove_if_overlap":

        print("input:\t" + args.i)
        print("output:\t" + args.o)
        
        if not args.t:
            print("Please define the file for target regions.")
            sys.exit(1)
        else:
            print("target:\t" + args.t)

        # with open(args.target) as f:
        t = GenomicRegionSet("targets")
        t.read_bed(args.t)

        # with open(args.i) as fi, open(args.o, "w") as fo:
        input_regions = GenomicRegionSet("input")
        input_regions.read_bed(args.i)
        if args.k:
            output_regions = input_regions.intersect(t, mode=OverlapType.ORIGINAL)
        else:
            output_regions = input_regions.subtract(t, whole_region=True)
        output_regions.write_bed(args.o)
        print("complete.")

    ############### BED add columns #############################################
    elif args.mode == "bed_add_columns":
        print("input:\t" + args.i)
        print("reference:\t" + args.ref)
        print("output:\t" + args.o)
        
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


    ############### BED divide by erxpression ###########################
    elif args.mode == "bed_divide":
        print("input:\t" + args.bed)
        print("table:\t" + args.t)
        
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
        bed.read_bed(args.bed)
        o1 = GenomicRegionSet(args.o1)
        o2 = GenomicRegionSet(args.o2)
        for r in bed:
            if r.name in gene1: o1.add(r)
            elif r.name in gene2: o2.add(r)
        o1.write_bed(args.o1)
        o2.write_bed(args.o2)


    ############### BAM filtering by BED ###########################
    #
    elif args.mode == "bam_filter":
        print("input:\t" + args.i)
        print("regions:\t" + args.bed)
        print("output prefix:\t" + args.o)

        bed = GenomicRegionSet("bed")
        bed.read_bed(args.bed)
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
        print("input:\t"+args.i)
        print("result from THOR:\t"+args.thor)
        print("output:\t"+args.o)
        if not os.path.exists(args.o):
            os.makedirs(args.o)
        print("extension:\t" + args.e)
        print("bin size:\t" + args.b)

        chr = "chr1"
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

        # print(max(bin_cov1))
        # print(max(bin_cov2))
        #
        # print(bin_cov1[0:10])
        # print(bin_cov2[0:10])

        bin_cov1 = [ numpy.log2(x) for x in bin_cov1]
        bin_cov2 = [ numpy.log2(x) for x in bin_cov2]

        # M =
        # A =
        #
        # print(max(M))
        # print(max(A))

        # from matplotlib.backends.backend_pdf import PdfPages
        # pp = PdfPages('MAplot.pdf')
        # print(bin_cov1)

        # plt.scatter(bin_cov1[0], bin_cov2[0], s=1, alpha=0.5)
        plt.figure(1)
        plt.subplot(221)
        plt.title('Sample 1 before norm.')
        plt.scatter([ 0.5*(bin_cov1[2] + bin_cov1[0]) ], [ bin_cov1[2] - bin_cov1[0] ], s=1, alpha=0.5)
        plt.xlim([0, 6])
        plt.ylim([-6, 6])
        plt.ylabel('M (log ratio)')
        # plt.xlabel('A (mean average)')
        plt.subplot(222)
        plt.title('Sample 1 after norm.')
        plt.scatter([ 0.5*(bin_cov2[2] + bin_cov2[0]) ], [ bin_cov2[2] - bin_cov2[0] ], s=1, alpha=0.5)
        plt.xlim([0, 6])
        plt.ylim([-6, 6])
        plt.ylabel('M (log ratio)')
        # plt.xlabel('A (mean average)')
        plt.subplot(223)
        plt.title('Sample 2 before norm.')
        plt.scatter([0.5 * (bin_cov1[3] + bin_cov1[1])], [bin_cov1[3] - bin_cov1[1]], s=1, alpha=0.5)
        plt.xlim([0, 6])
        plt.ylim([-6, 6])
        # plt.ylabel('M (log ratio)')
        plt.xlabel('A (mean average)')
        plt.subplot(224)
        plt.title('Sample 2 after norm.')
        plt.scatter([0.5 * (bin_cov2[3] + bin_cov2[1])], [bin_cov2[3] - bin_cov2[1]], s=1, alpha=0.5)
        plt.xlim([0, 6])
        plt.ylim([-6, 6])
        # plt.ylabel('M (log ratio)')
        plt.xlabel('A (mean average)')
        # plt.savefig(pp, format='pdf')
        # pp.close()
        plt.savefig(os.path.join(args.o,'MAplot.png'), bbox_inches='tight')


        print("finish")







    ############### THOR split #############################################
    elif args.mode == "thor_split":
        print("input:\t" + args.i)
        if not args.o:
            args.o = os.path.dirname(args.i)
        print("output:\t" + args.o)

        name = os.path.basename(args.i).split(".")[0]
        if args.fc == 0: tag = "_p" + str(args.p)
        else: tag = "_p"+str(args.p)+"_fc"+str(args.fc)

        bed = GenomicRegionSet("input")
        bed.read_bed(args.i)
        print("Number of input peaks:\t"+str(len(bed)))

        if args.rn and args.g:
            bed2 = bed.gene_association(organism=args.g)
        else:
            bed2 = bed

        for region in bed2:
            data = region.data.split()
            stat = data[5].split(";")
            s1 = [float(x) + 1 for x in stat[0].split(":")]
            s2 = [float(x) + 1 for x in stat[1].split(":")]
            fc = math.log((sum(s2) / len(s2)) / (sum(s1) / len(s1)), 2)
            region.data = "\t".join([str(fc)] + data[1:])
        # bed2.write_bed(args.o)

        # gain_f = open(os.path.join(args.o,name+"_"+str(args.p)+"_gain.bed"), "w")
        # lose_f = open(os.path.join(args.o,name+"_"+str(args.p)+"_lose.bed"), "w")
        gain_peaks = GenomicRegionSet("gain_peaks")
        lose_peaks = GenomicRegionSet("lose_peaks")
        mix = GenomicRegionSet("mix")

        for region in bed2:
            l = region.data.split()
            s = l[5].split(";")
            if abs(float(l[0])) > args.fc and float(s[2]) > args.p:
                if float(l[0]) > 0:
                    gain_peaks.add(region)
                elif float(l[0]) < 0:
                    lose_peaks.add(region)
        gain_peaks.write_bed(os.path.join(args.o, name + tag + "_gain.bed"))
        lose_peaks.write_bed(os.path.join(args.o, name + tag + "_lose.bed"))

        for region in bed2:
            l = region.data.split()
            s = l[5].split(";")
            if abs(float(l[0])) > args.fc and float(s[2]) > args.p:
                s1 = sum([int(x) for x in s[0].split(":")]) / len(s[0].split(":"))
                s2 = sum([int(x) for x in s[1].split(":")]) / len(s[1].split(":"))
                length = abs(region.final - region.initial)
                ns1 = float(s1) / length
                ns2 = float(s2) / length
                region.data = "\t".join([l[0], str(s1), str(s2), str(length),
                                         str(ns1), str(ns2), str(ns1 + ns2), str(ns1 - ns2), s[2]])
                mix.add(region)

        mix.write_bed(os.path.join(args.o, name+tag+".table"))
        print("Number of gain peaks:\t" + str(len(gain_peaks)))
        print("Number of lose peaks:\t" + str(len(lose_peaks)))
        
        
    ############### getseq #############################################
    elif args.mode == "getseq":
        #print("input:\t" + args.i)
        #print("output:\t" + args.o)
        
        if args.b:
            regions = GenomicRegionSet("regions")
            regions.read_bed(args.b)
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
                print(args.p)
                args.ch = args.p.partition(":")[0]
                args.ss = int(args.p.partition(":")[2].partition("-")[0])
                args.es = int(args.p.partition(":")[2].partition("-")[2])
            if args.s == "both":
                seq1 = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand="+",
                                    reverse=False, complement=False, rna=args.r, ex=args.ex)
                seq2 = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand="-",
                                    reverse=False, complement=True, rna=args.r, ex=args.ex)
                print("5'- " + seq1 + " -3'")
                print("3'- " + seq2 + " -5'")
            elif args.s == "+" and not args.re:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.re, complement=args.c, rna=args.r, ex=args.ex)
                print("5'- " + seq + " -3'")
            elif args.s == "+" and args.re:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.re, complement=args.c, rna=args.r, ex=args.ex)
                print("3'- " + seq + " -5'")
            elif args.s == "-" and not args.re:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.re, complement=args.c, rna=args.r, ex=args.ex)
                print("3'- "+seq+" -5'")
            elif args.s == "-" and args.re:
                seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s,
                                   reverse=args.re, complement=args.c, rna=args.r, ex=args.ex)
                print("5'- " + seq + " -3'")
            
        print()
    ############### WIG trim end #############################################
    elif args.mode == "wig_trim_end":
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("chromosome size:\t" + args.chrosize)
        
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
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("organism:\t" + args.organism)
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
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("organism:\t" + args.organism)
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
        print("input:\t" + args.i)
        print("output:\t" + args.o)

        def fasta2bp(filename):
            s = ""
            with open(filename) as f:
                for line in f:
                    l = line.strip()
                    if l.startswith(">"): continue
                    elif not l: continue
                    else: s += l
            return len(s)

        if os.path.isfile(args.i):
            l = fasta2bp(args.i)
            print("Length:\t\t" + str(l)+" bp")

        elif os.path.isdir(args.i) and args.o:
            list_bp = []
            for root, dirs, files in os.walk(args.i):
                for f in files:
                    if f.endswith(".fa") or f.endswith(".fasta"):
                        # print(f.partition(".")[0])
                        # print(parser_fasta2bp(filename=os.path.join(root,f)))
                        list_bp.append([f.partition(".")[0],
                                        str(fasta2bp(os.path.join(root,f)))])
            with open(args.o, "w") as g:
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
        tcons.read_bed(args.t)
        circrna = GenomicRegionSet("circRNA")
        circrna.read_bed(args.o)
        circ_inTCON = circrna.intersect(y=tcons, mode = OverlapType.COMP_INCL)
        circ_TCONs = tcons.intersect(y=circ_inTCON, mode = OverlapType.ORIGINAL)
        #print(len(circ_TCONs))
        circ_TCONs.write_bed(args.c)


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

