from __future__ import print_function
from rgt.GenomicRegionSet import *
from rgt.Util import GenomeData
import argparse 
import os        

##################################################################################
parser = argparse.ArgumentParser(description='Replace TCONs in BED file by assoicated gene names', 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-bed', type=str, help="BED file or a directory containing BED files")
parser.add_argument('-output', type=str, help="Define the output directory")
parser.add_argument('-organism', type=str, help="Define the organism")
args = parser.parse_args()




genome = GenomeData(args.organism)

if os.path.isfile(args.bed):
    regionset = GenomicRegionSet("bed")
    regionset.read(args.bed)
    gr = regionset.gene_association(organism=args.organism, promoterLength=1000, 
                                    threshDist=500000, show_dis=True)
    regionset.replace_region_name(gr,combine=True)
    
    regionset.write(args.output)

elif os.path.isdir(args.bed):
    if not os.path.exists(args.output):
        os.makedirs(args.output)
    for root, dirnames, filenames in os.walk(args.bed):
            
        for filename in filenames:
            if ".bed" in filename:
                print(filename)
                fnn = os.path.basename(filename)
                fn = fnn.partition(".bed")[0]
                try:
                    regionset = GenomicRegionSet("bed")
                    regionset.read(os.path.join(args.bed,fnn))
                    gr = regionset.gene_association(organism=args.organism, promoterLength=1000, 
                                                    threshDist=500000, show_dis=True)
                    regionset.replace_region_name(gr,combine=True)
                    regionset.write(os.path.join(args.output, fn+"_associated.bed"))
                    
                except:
                    pass
