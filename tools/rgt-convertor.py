# Python Libraries
from __future__ import print_function
import sys
import argparse
import os

# Local Libraries
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
tag = "RGT-convertor"

def main():
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
    parser_gtfat = subparsers.add_parser('gtf_add_transcripts', help="[GTF] Add transcripts from the existed exons")
    parser_gtfat.add_argument('-i', '-input', type=str, help="Input GTF file")
	parser_gtfat.add_argument('-o', '-output', type=str, help="Output GTF file")

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
    parser_bedgp.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedgp.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedgp.add_argument('-l', type=int, default=1000, 
                              help="Define length of promoters (default:1000bp)")

    ##########################################################################
	args = parser.parse_args()
    ##########################################################################


    ############### GTF add transcripts ######################################
    if args.mode == "gtf_add_transcripts":
        print(tag+": [GTF] Add transcripts from the existed exons")
        print("input: " + args.i)
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
        print("Number of entries: "+ str(len(entries)))
        print("Number of transcripts: "+str(len(transcript)))
        print("output: " + args.o)
        with open(args.o, "w") as out:
            for t in transcript:
                print("\t".join(t), file=out)

    ############### BED add score ############################################
    elif args.mode == "bed_add_score":
        print(tag+": [BED] Add scores")
        print("input: " + args.i)
        with open(args.i) as f, open(args.o, "w") as g:
            for line in f:
                line = line.strip() 
                print(line+"\t"+args.v, file=g)

    ############### BED extend ###############################################
    elif args.mode == "bed_extend":
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

        gene = GenomicRegionSet("genes")
        gene.read_bed(args.i)

        promoter = GenomicRegionSet("promoter")

        promoterLength = int(args.l)

        for s in gene:
            if s.orientation == "+": 
                s.initial, s.final = max(s.initial-promoterLength, 0), s.initial

            else: s.initial, s.final = s.final, s.final+promoterLength
            promoter.add(s)

        #print(len(promoter))
        promoter.write_bed(args.o)


