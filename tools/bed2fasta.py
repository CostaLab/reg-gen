from __future__ import print_function
from rgt.GenomicRegionSet import *
from rgt.Util import GenomeData
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse 
import os


##################################################################################
def load_exon_sequence(bed, directory, genome_path):
    """Load the exon sequence from the the transcripts. 
    Input BED format should contain:
        blockCount - The number of blocks (exons) in the BED line.
        blockSizes - A comma-separated list of the block sizes.
        blockStarts - A comma-separated list of block starts. 
        see details: http://genome.ucsc.edu/FAQ/FAQformat#format1

    Output:
        Each FASTA file represants a transcript and contains all the exons within the file.

    """
    regionset = GenomicRegionSet("bed")
    regionset.read(bed)
    regionset.sort()

    
    genome = pysam.Fastafile(genome_path)
    
    try:
        if len(regionset.sequences[0].data.split("\t")) == 7: 
            blockinfor = True
            no_exon = False
    except:
        blockinfor = False
        regionset.sequences.sort(key=lambda g: g.name)
        no_exon = True

    if blockinfor:
        
        for gr in regionset:
            if not gr.name:
                print("Error: For fetching exon sequences, please define the transcript name.")
                sys.exit()
            else:
                if not os.path.exists(directory):
                    os.makedirs(directory)
                f = open(os.path.join(directory, gr.name+".fa"), "w")
                data = gr.data.split("\t")
                #print(len(data))
                if len(data) == 7:
                    #print(data)
                    n = int(data[4])
                    
                    blocks = [ int(b) for b in filter(None, data[5].split(",")) ]
                    starts = [ int(s) for s in filter(None, data[6].split(",")) ]
                    printstr = []

                    for i in range(n):
                        start = gr.initial + starts[i]
                        end = start + blocks[i]
                        if no_exon and i == 0:
                            ex = ""
                        elif gr.orientation == "-":
                            ex = "exon:"+str(n-i)
                        else:
                            ex = "exon:"+str(i+1)

                        if gr.orientation == "-":
                            seq = Seq(genome.fetch(gr.chrom, start-1, end-1), IUPAC.unambiguous_dna)
                            seq = seq.reverse_complement()
                            p = [ ">"+ " ".join([ gr.name, 
                                                  ex, 
                                                  "_".join(["REGION",gr.chrom,
                                                            str(start),str(end), 
                                                            gr.orientation]) ]),
                                  seq ]
                            
                            printstr.append(p)
                            

                        else:
                            p = [ ">"+ " ".join([gr.name, ex, 
                                  "_".join(["REGION",gr.chrom,str(start),str(end), gr.orientation]) ]),
                                  genome.fetch(gr.chrom, start-1, end-1)
                                ]
                            printstr.append(p)
                            

                    if gr.orientation == "-": printstr = printstr[::-1]
                    for i in range(n):
                        print(printstr[i][0], file=f)
                        print(printstr[i][1], file=f)
                        

                else:
                    print("Warning: The given regions have no block information, please try write_bed_blocks")
                f.close()
    else:
        pre_id = ""
        for gr in regionset:
            if not gr.name: 
                gr.name = gr.toString()

            if pre_id == "": 
                pre_id = gr.name
                z = GenomicRegionSet(gr.name)
                z.add(gr)
            elif gr.name == pre_id:
                z.add(gr)
            else:
                f = open(os.path.join(directory, pre_id+".fa"), "w")
                for i, g in enumerate(z):
                    try: regiontag = "_".join(["REGION", g.chrom, str(g.initial),str(g.final), gr.orientation])
                    except: regiontag = "_".join(["REGION", g.chrom, str(g.initial),str(g.final)] )

                    print( ">"+ " ".join([g.name,  
                                          regiontag ]), file=f)
                    print(genome.fetch(g.chrom, g.initial, g.final), file=f)
                f.close()

                pre_id = gr.name
                z = GenomicRegionSet(gr.name)
                z.add(gr)

        # Last TX
        f = open(os.path.join(directory, pre_id+".fa"), "w")
        for i, g in enumerate(z):
            try: regiontag = "_".join(["REGION", g.chrom, str(g.initial),str(g.final), gr.orientation])
            except: regiontag = "_".join(["REGION", g.chrom, str(g.initial),str(g.final)] )
            print( ">"+ " ".join([g.name, 
                                  regiontag ]), file=f)
            print(genome.fetch(g.chrom, g.initial, g.final), file=f)
        f.close()
        



##################################################################################
parser = argparse.ArgumentParser(description='Convert BED files into FASTAs', 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-bed', type=str, help="BED file or a directory containing BED files")
parser.add_argument('-output', type=str, help="Define the output directory")
parser.add_argument('-organism', type=str, help="Define the organism")
args = parser.parse_args()

if not os.path.exists(args.output):
    os.makedirs(args.output)

genome = GenomeData(args.organism)

if os.path.isfile(args.bed):
    load_exon_sequence(bed=args.bed, directory=args.output, genome_path=genome.get_genome())

elif os.path.isdir(args.bed):

    for root, dirnames, filenames in os.walk(args.bed):
            
        for filename in filenames:
            if ".bed" in filename:
                print(filename)
                fn = os.path.basename(filename)
                fn = fn.partition(".bed")[0]
                try:
                    load_exon_sequence(bed=os.path.join(args.bed,filename), 
                                       directory=os.path.join(args.output,fn), 
                                       genome_path=genome.get_genome())
                except:
                    pass
