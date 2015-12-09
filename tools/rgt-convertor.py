# Python Libraries
from __future__ import print_function
import sys
import argparse
import os

# Local Libraries
from rgt.GeneSet import GeneSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.AnnotationSet import AnnotationSet
from rgt.Util import OverlapType
tag = "RGT-convertor"

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
        #print("3'- "+seq+" -5'")
    else:
        return seq        
        #print("5'- "+seq+" -3'")

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
    parser_gtfat.add_argument('-i', '-input', type=str, help="Input GTF file")
    parser_gtfat.add_argument('-o', '-output', type=str, help="Output GTF file")

    ############### GTF to BED ###############################################
    parser_gtf2bed = subparsers.add_parser('gtf_to_bed', 
                                           help="[GTF] Convert GTF file to BED by the given biotype")
    parser_gtf2bed.add_argument('-i', '-input', type=str, help="Input GTF file")
    parser_gtf2bed.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_gtf2bed.add_argument('-t', '-target', type=str, help="Define the target Biotype\
                                                                 (e.g. gene or exon or transcript...)")
    parser_gtf2bed.add_argument('-b', '-block', action="store_true", 
                                help="Save exons into entries with block in BED")

    ############### GTF to FASTA #############################################
    parser_gtf2fasta = subparsers.add_parser('gtf_to_fasta', 
                                             help="[GTF] Convert GTF file to FASTA (exons) by the given gene name")
    parser_gtf2fasta.add_argument('-i', '-input', type=str, help="Input GTF file")
    parser_gtf2fasta.add_argument('-o', '-output', type=str, help="Output FASTA file")
    parser_gtf2fasta.add_argument('-t', '-transcript', type=str, help="Define the target transcript")
    parser_gtf2fasta.add_argument('-g', '-gene', type=str, help="Define the target gene")
    parser_gtf2fasta.add_argument('-genome', type=str, help="Define the FASTA file of the genome")
    
    ############### BED add score ############################################
    parser_bedac = subparsers.add_parser('bed_add_score', help="[BED] Add score column")
    parser_bedac.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedac.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedac.add_argument('-v', type=str, help="Define value to add")

    ############### BED extend ###############################################
    parser_bedex = subparsers.add_parser('bed_extend', help="[BED] Extend the regions")
    parser_bedex.add_argument('-i', type=str, help="Input BED file")
    parser_bedex.add_argument('-o', type=str, help="Output BED name.")
    parser_bedex.add_argument('-oz', "-onlyzero", action="store_true", default=False, 
                              help="Extend only the zero-length regions")
    parser_bedex.add_argument('-l', type=int, help="Define the length to extend.")
    parser_bedex.add_argument('-both',action="store_true", default=False, 
                              help="Extend from the both ends.")

    ############### BED get promoters ########################################
    parser_bedgp = subparsers.add_parser('bed_get_promoters', 
                       help="[BED] Get promoters from the given genes")
    parser_bedgp.add_argument('-i', '-input', type=str, help="Input file (BED or gene list)")
    parser_bedgp.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedgp.add_argument('-organism', type=str, help="Define the organism (necessary if input is a gene list)")
    parser_bedgp.add_argument('-l', type=int, default=1000, 
                              help="Define length of promoters (default:1000bp)")
    
    ############### BED to FASTA #############################################
    parser_bed2fasta = subparsers.add_parser('bed_to_fasta', 
                       help="[BED] Export the sequences in FASTA according to the given BED file")
    parser_bed2fasta.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bed2fasta.add_argument('-o', '-output', type=str, help="Output directory for FASTA files")
    parser_bed2fasta.add_argument('-genome', type=str, help="Define the FASTA file of the genome sequence")

    ############### WIG trim ends by chromosome #############################################
    parser_wig_trim = subparsers.add_parser('wig_trim_end', 
                       help="[WIG] Trim the WIG file according to the given chromosome size")
    parser_wig_trim.add_argument('-i', '-input', type=str, help="Input WIG file")
    parser_wig_trim.add_argument('-o', '-output', type=str, help="Output WIG file")
    parser_wig_trim.add_argument('-chrosize', type=str, help="Define path to the chromosome size file")


    ############### STAR junction to BED #############################################
    parser_circRNA = subparsers.add_parser('circRNA', 
                       help="[junction] Convert the Chimeric junction from STAR to BED file")
    parser_circRNA.add_argument('-i', '-input', type=str, help="Input chimeric junction file from STAR")
    parser_circRNA.add_argument('-t', '-tcons', type=str, help="Input BED file of tcons")
    parser_circRNA.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_circRNA.add_argument('-c', '-circ', type=str, help="Output BED file of circular RNA")

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

        if "/" in args.o:
            if not os.path.exists(args.o.rpartition("/")[0]):
                os.makedirs(args.o.rpartition("/")[0])
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
        print(tag+": [GTF] Convert GTF to BED")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("target:\t" +args.t)

        if args.t == "gene" or args.t == "transcript":
            with open(args.i) as f, open(args.o, "w") as g:
                for line in f:
                    if line[0] == "#": continue
                    line = line.split()
                    #print(line)
                    if line[2] == args.t:
                        #print(line[9].split('"')[1])
                        print("\t".join([line[0], line[3], line[4], line[17][1:-2], ".", line[6]]), file=g)

        elif args.t == "exon":
            exons = GenomicRegionSet("exons")
            with open(args.i) as f:
                for line in f:
                    if line[0] == "#": continue
                    line = line.split()
                    #print(line)
                    if line[2] == args.t:
                        #print(line[9].split('"')[1])
                        exons.add(GenomicRegion(chrom=line[0], initial=int(line[3]), final=int(line[4]), 
                                                name=line[17][1:-2], orientation=line[6]))         
            exons.write_bed_blocks(args.o)
        

    ############### GTF to FASTA #############################################
    elif args.mode == "gtf_to_fasta":
        print(tag+": [GTF] Export certain gene or transcripts into FASTA sequence")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        



    ############### BED add score ############################################
    elif args.mode == "bed_add_score":
        print(tag+": [BED] Add scores")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        with open(args.i) as f, open(args.o, "w") as g:
            for line in f:
                line = line.strip() 
                print(line+"\t"+args.v, file=g)

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

    ############### BED get promoters #########################################
    elif args.mode == "bed_get_promoters":
        print("input:\t" + args.i)
        print("output:\t" + args.o)
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
            promoter = ann.get_promoters(promoterLength=args.l, gene_set=de_gene, unmaplist=False)
            #print(len(de_prom))

        
        #print(len(promoter))
        promoter.write_bed(args.o)


    ############### BED to FASTA #############################################
    elif args.mode == "bed_to_fasta":
        print("input:\t" + args.i)
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



    ############### BED to FASTA #############################################
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