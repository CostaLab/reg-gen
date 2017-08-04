from __future__ import print_function
from __future__ import division
import os, sys
import argparse
import random
import copy
import os
import time
from scipy import stats
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
import os
from rgt.Util import GenomeData
from rgt.Util import OverlapType

present_path = os.getcwd()
bed_dir = present_path+"/bed/"
fasta_dir = present_path+"/fasta/"
genome_path = 'genome.fa'
chromosome_size_path = '/media/931ef578-eebe-4ee8-ac0b-0f4690f126ed/data/rgtdata/hg19/chrom.sizes'

"""
Updated on 22 May 2014 by Joseph
"""

#############################  Parameters    ##############################
parser = argparse.ArgumentParser(description='Return the random sequences according to the given parameters.')
parser.add_argument('-o','-organism', default= "hg19", help='Define the organism. Default: hg19')
parser.add_argument('-l','-length', type=int, help='Define the length of each sequence.')
parser.add_argument('-n','-number', type=int, help='Define the number of random regions.')
parser.add_argument('-f','-filter', default=None, help='Given the path to the BED file as the filter.')

args = parser.parse_args()

# Setup the entries
region = GenomicRegion("sample", initial=0, final=args.l)
template = GenomicRegionSet("tamplate")
template.add(region)
    
if not os.path.exists(bed_dir):
    os.makedirs(bed_dir)
    
# Random region
result = template.random_regions(organism= "hg19", total_size=args.n, multiply_factor=0, overlap_result=True, overlap_input=True, chrom_X=False, chrom_M=False, filter_path=args.f)
result.write(os.path.join(bed_dir, "00total.bed"))
chrom = GenomicRegionSet("chrom")
chrom.get_genome_data(organism=args.o, chrom_X=False, chrom_M=False)
            
chrom_list = []
for r in chrom.sequences:
    chrom_list.append(r.chrom)
        
print("Settings:\n\tAllowing overlapping within random regions.")
print("Randomizing method:\n\tChoose chromosome with the same possibility (1/23 for each chromosome)\n" +
      "\tChoose random position with the given length without the regard of previous chosen regions.")
print("The distribution of the results are:")
for c in chrom_list:
    print("\t"+str(len(result.any_chrom(c)))+" random regions on " + c)


# Write the bed file indivdually

for i, r in enumerate(result.sequences):
    f = open(bed_dir+str(i+1)+".bed",'w')
    f.write(r.__str__())
    f.close()
    
# FASTA

if not os.path.exists(fasta_dir):
    os.makedirs(fasta_dir)

for i in range(1, args.n+1):
    os.system("fastaFromBed -fi " + genome_path + " -bed " + bed_dir + str(i) + ".bed -fo " + fasta_dir + str(i) + ".fa")
