# Python Libraries
from __future__ import print_function
import sys
import argparse
import os
from os.path import expanduser
home = expanduser("~")

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
                                help="Define the feature {gene,transcript,exon,CDS,UTR,start_codon,stop_codon,Selenocysteine}")
    parser_gtf2bed.add_argument('-t', '-type', type=str, default="All", 
                                help="Define gene type e.g. 'protein_coding' more: http://www.gencodegenes.org/gencode_biotypes.html")
    parser_gtf2bed.add_argument('-st', '-status', type=str, default="All", 
                                help="Define gene status {KNOWN, NOVEL, PUTATIVE,All}")
    parser_gtf2bed.add_argument('-g', '-gene', type=str, default=None, 
                                help="Define the gene list for filtering, default is None.")
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

    ############### BED add columns ################################
    parser_bedaddcol = subparsers.add_parser('bed_add_columns', 
                       help="[BED] Add extra columns to the BED file by gene name")
    parser_bedaddcol.add_argument('-i', '-input', type=str, help="Input BED file")
    parser_bedaddcol.add_argument('-o', '-output', type=str, help="Output BED file")
    parser_bedaddcol.add_argument('-ref', type=str, help="Define file for referring the extra columns ")
    parser_bedaddcol.add_argument('-f', '-field', type=int, help="Which field of the reference file is compared for names.")

    ############### GENOME get sequence ####################################################
    parser_getseq = subparsers.add_parser('getseq', 
                       help="[FASTA] Get sequence from genome FASTA")
    parser_getseq.add_argument('-o', type=str, help="Output FASTA file")
    parser_getseq.add_argument('-d', '-dna', type=str, help="DNA sequence in FASTA format")
    parser_getseq.add_argument('-b', '-bed', type=str, help="Input BED file")
    parser_getseq.add_argument('-p', '-pos', type=str, help="position")
    parser_getseq.add_argument('-s', '-strand', type=str, default="+", help="strand (+ or -)")
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
        print(tag+": [GTF] Convert GTF to BED")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("source:\t\t" + args.s)
        print("feature:\t" + args.f)
        print("type:\t\t" + args.t)
        print("status:\t\t" + args.st)
        if args.g:
            print("gene list:\t" + args.g)
        
        if args.s == "ENSEMBL": tag_s = "ENSEMBL"
        elif args.s == "HAVANA": tag_s = "HAVANA"
        elif args.s == "All": tag_s = None
        else: 
            print("Please redefine the argument -s.")
            sys.exit(1)

        if args.f == "gene": tag_f = "gene"
        elif args.f == "transcript": tag_f = "transcript"
        elif args.f == "exon": tag_f = "exon"
        elif args.f == "CDS": tag_f = "CDS"
        elif args.f == "UTR": tag_f = "UTR"
        elif args.f == "start_codon": tag_f = "start_codon"
        elif args.f == "stop_codon": tag_f = "stop_codon"
        elif args.f == "Selenocysteine": tag_f = "Selenocysteine"
        elif args.f == "All": tag_f = None
        else: 
            print("Please redefine the argument -f.")
            sys.exit(1)

        if args.t == "protein_coding": tag_t = ' \"protein_coding\"'
        elif args.t == "non_coding": tag_t = ' \"non_coding\"'
        elif args.t == "known_ncrna": tag_t = ' \"known_ncrna\"'
        elif args.t == "pseudogene": tag_t = ' \"pseudogene\"'
        elif args.t == "All": tag_t = None
        else: 
            print("Please redefine the argument -t.")
            sys.exit(1)

        if args.st == "KNOWN": tag_st = 'gene_status \"KNOWN\"'
        elif args.st == "NOVEL": tag_st = 'gene_status \"NOVEL\"'
        elif args.st == "PUTATIVE": tag_st = 'gene_status \"PUTATIVE\"'
        elif args.st == "All": tag_st = None
        else: 
            print("Please redefine the argument -st.")
            sys.exit(1)

        if args.g:
            select_genes = GeneSet("genes")
            select_genes.read(args.g)

        ### Indexing
        with open(args.i, "r") as f,open(args.o, "w") as g:
            for line in f:
                if line[0] == "#": continue
                line = line.strip().split("\t")

                if len(line) < 5: continue
                if tag_s and tag_s != line[1]: continue
                if tag_f and tag_f != line[2]: continue
                if tag_t and tag_t in line[8]: pass
                elif not tag_t: pass
                else: continue              
                if tag_st and tag_st in line[8]: pass
                elif not tag_st: pass
                else: continue
                # print(line)
                # print(line[8].split("; "))
                info = line[8].split("; ")
                
                gn = [s for s in info if "gene_name" in s][0].partition("\"")[2][:-1]
                gn = gn.partition(" (")[0]
                gi = [s for s in info if "gene_id" in s][0].partition("\"")[2][:-1]
                # print(gn)
                # print("\t".join([line[0], line[3], line[4], line[ind+1][1:-2], ".", line[6]]))
                if int(line[3]) < int(line[4]):
                    if line[0].isdigit():
                        ch = "chr" + line[0]
                        seq = "\t".join([ch, line[3], line[4], gn, ".", line[6]])
                    elif line[0].startswith("chr"):
                        ch = line[0]
                        seq = "\t".join([ch, line[3], line[4], gn, ".", line[6]])
                    else:
                        continue
                    
                else:
                    if line[0].isdigit():
                        ch = "chr" + line[0]
                        seq = "\t".join([ch, line[4], line[3], gn, ".", line[6]])
                    elif line[0].startswith("chr"):
                        ch = line[0]
                        seq = "\t".join([ch, line[4], line[3], gn, ".", line[6]])
                    else:
                        continue
                # print(seq)

                if not args.g:
                    print(seq, file=g)
                elif select_genes.check(gn) or select_genes.check(gi):
                    
                    print(seq, file=g)
                else:
                    continue

        if args.b:
            exons = GenomicRegionSet("output")
            exons.read_bed(args.o)
            exons.write_bed_blocks(args.o)

        # sys.exit(1)

        # if args.g:
        #     select_genes = GeneSet("genes")
        #     select_genes.read(args.g)

        # # if args.t == "gene" or args.t == "transcript":
        # with open(args.i, "r") as f,open(args.o, "w") as g:
        #     find_ind = False
        #     for line in f:
        #         if line[0] == "#": 
        #             continue
        #         elif args.known_only:
        #             if "gene_status \"KNOWN\"" in line:
        #                 pass
        #             else:
        #                 continue

        #         line = line.split()
        #         #print(line)
        #         if not find_ind:
        #             for i, l in enumerate(line):
        #                 if "gene_name" in l: 
        #                     ind = i
        #                     find_ind = True

        #         if line[2] == args.t:
                    
        #             if not args.g:
        #                 print(line[ind+1][1:-2])
        #                 print("\t".join([line[0], line[3], line[4], line[ind+1][1:-2], ".", line[6]]), file=g)
        #             elif select_genes.check(line[ind+1][1:-2]):
        #                 print(line[ind+1][1:-2])
        #                 print("\t".join([line[0], line[3], line[4], line[ind+1][1:-2], ".", line[6]]), file=g)
        #             else:
        #                 continue


        # elif args.t == "exon":
        #     exons = GenomicRegionSet("exons")
        #     with open(args.i) as f:
        #         for line in f:
        #             if line[0] == "#": continue
        #             line = line.split()
        #             #print(line)
        #             if line[2] == args.t:
        #                 #print(line[9].split('"')[1])
        #                 exons.add(GenomicRegion(chrom=line[0], initial=int(line[3]), final=int(line[4]), 
        #                                         name=line[15][1:-2], orientation=line[6]))         
        #     exons.write_bed_blocks(args.o)
        

    ############### GTF to FASTA #############################################
    elif args.mode == "gtf_to_fasta":
        print(tag+": [GTF] Export certain gene or transcripts into FASTA sequence")
        print("input:\t" + args.i)
        print("output:\t" + args.o)
        print("genome:\t" + args.genome)
        print("gene list:\t" + args.g)

        ann = AnnotationSet(gene_source=genome,filter_havana=False, protein_coding=False, known_only=False)
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
        print("organism:\t" + args.organism)

        bed = GenomicRegionSet(args.i)
        bed.read_bed(args.i)
        renamebed = bed.gene_association(gene_set=None, organism=args.organism, 
                                         promoterLength=args.l, 
                                         threshDist=args.t, show_dis=args.d)
        renamebed.write_bed(args.o)

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
        
        #print(len(promoter))
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
        output_regions = input_regions.subtract(t,whole_region=True)
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

    ############### BED add columns #############################################
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
            seq = get_sequence(sequence=args.d, ch=args.ch, ss=args.ss, es=args.es, strand=args.s, 
                               reverse=args.re, complement=args.c, rna=args.r, ex=args.ex)

            print(seq)
            
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
        print("5' - "+ seq.sequences[0].seq[start:end]+ " - 3'")

