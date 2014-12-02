#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%filterVCF [options] <CONFIG> 

Filter VCF files by pipeline described in README.
<CONFIG> file must be a space separated list of:

<sample's name> <path to VCF>

@author: Manuel Allhoff

"""
from __future__ import print_function
#import pydevd;pydevd.settrace()
from optparse import OptionParser
import sys
from rgt.GenomicVariantSet import GenomicVariantSet
from rgt.GenomicRegionSet import GenomicRegionSet
from max_density import AlgGoldwasser
from itertools import combinations
from copy import copy

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def load_data(path_sample_vcf):
    """Return list of GenomicVarientSets described by config file"""
    data = []
    with open(path_sample_vcf) as file:
        for line in file:
            tmp = line.split(" ")
            name, path = tmp[0], tmp[1]
            path = path.strip()
            data.append(GenomicVariantSet(path, name=name))
            
    return data

def print_length(l, info):
    """Print length of each GenomicVariant set in list l"""
    names = map(lambda x: x.name, l)
    print(info, file=sys.stderr)
    for i, name in enumerate(names):
        print(name, len(l[i]), file=sys.stderr, sep='\t')

def get_max_density(GenomicVariantSets, lowerBound=20, upperBound=50, max_it=5):
    """Compute max. density of homozygous SNPs and write results to stderr"""
    for GenomicVariantSet in GenomicVariantSets:
        max_pos = 0
        for record in GenomicVariantSet:
            if 'Mask' not in record.filter and not record.samples[0].is_het:
                if record.pos > max_pos:
                    max_pos = record.pos

    density_seq = [ [0, 1] for _ in range(max_pos + 2 * upperBound) ] 
    
    i=0
    for GenomicVariantSet in GenomicVariantSets:
        for record in GenomicVariantSet:
            if 'Mask' not in record.filter and not record.samples[0].is_het:
                density_seq[record.pos][0] += 1
                i += 1
        
    density_seq = map(lambda x: tuple(x), density_seq)
    
    print('Max. Density Iterations, lower bound: %s, upper bound: %s' %(lowerBound, upperBound), file=sys.stderr)
    print('(Regions to consider: %s)' %i, file=sys.stderr)
    
    for i in range(max_it):
        den, coord = AlgGoldwasser(density_seq, lowerBound, upperBound)
	s = sum(map(lambda x: x[0], density_seq[coord[0]:coord[1]]))
	print("%s. It.: density %s at %s - %s with %s SNPs" %(i+1, round(den, 5), coord[0], coord[1], s), file=sys.stderr)
        density_seq[coord[0]:coord[1]] = [(0,1)] * (coord[1] - coord[0])
        

def input():
    """Parse options"""
    parser = HelpfulOptionParser(usage=__doc__)
     
    parser.add_option("--t-mq", dest="t_mq", default=20, help="Threshold for mapping quality (MQ) [default: %default]")
    parser.add_option("--t-dp", dest="t_dp", default=20, help="Threshold for combined depth (DP) [default: %default]")
    parser.add_option("--dbSNP", dest="c_dbSNP", default=None, help="Check for dbSNP [default: %default]")
    parser.add_option("--list-WT", dest="list_wt", default=None, help="List of WildTypes [default: %default]")
    parser.add_option("--bed", dest="list_bed", default=None, help="Filter against BED file (e.g. TFBS) [default: %default]")
    parser.add_option("--max-density", dest="max_density", default=False, action="store_true", help="Perform max. density search [default: %default]")
    parser.add_option("--lowerBound", dest="lower_bound", default=15000, help="lower window bound for max. density search [default: %default]")
    parser.add_option("--upperBound", dest="upper_bound", default=30000, help="upper window bound for max. density search [default: %default]")
    
    (options, args) = parser.parse_args()
    
    vcf_list = args[0]
     
    i = 1
    if len(args) != i:
        parser.error("Exactly %s parameters are needed" %i)
        
    return options, vcf_list

def pipeline(sample_data, options):
    """Filter sample data by MQ, DP and dbSNP"""
    print_length(sample_data, "#unfiltered variants")
    
    #filter MQ
    for sample in sample_data:
        sample.filter('MQ', '>=', options.t_mq)
    
    print_length(sample_data, "#variants with MQ >= %s" %options.t_mq)
    
    #filter DP
    for sample in sample_data:
        sample.filter('DP', '>=', options.t_dp)
    
    print_length(sample_data, "#variants with DP >= %s" %options.t_dp)
    
    if options.c_dbSNP:
        for sample in sample_data:
            sample.filter_dbSNP()
    
        print_length(sample_data, "#variants after filtering by dbSNP")
    else:
        print("#Do not filter by dbSNP", file=sys.stderr)

def output_intersections(sample_data):
    """Compute and write VCF of all possible sample data subsets' intersection."""
    for i in range(2, len(sample_data)+1):
        for tuple in combinations(sample_data, i): 
        #tuples have size i and are all subsets of sample_data with size i
            for j, sample in enumerate(list(tuple)):
                if j == 0:
                    intersect_variants = GenomicVariantSet(name = 'intersection')
                    intersect_variants.reader = copy(sample.reader)
                    intersect_variants.sequences = copy(sample.sequences)
                else:
                    intersect_variants.intersect(sample)
            
            name = "-".join(map(lambda x: x.name, list(tuple)))
            print("Intersection:", name, len(intersect_variants), sep='\t', file=sys.stderr)  
            intersect_variants.write_vcf('intersect-%s.vcf' %name)

def main():
    options, vcf_list = input()
    
    #thres_mq = 20
    #thres_dp = 20 
    #filter_dbSNP = True
    #tfbs_motifs_path = '/home/manuel/workspace/cluster_p/human_genetics/exp/exp01_motifsearch_sox2/humangenetics_motifs/Match/chr11_mpbs.bed'
    
    sample_data = load_data(vcf_list)
    print("##Filter variants of samples", file=sys.stderr)
    pipeline(sample_data, options)
    
    if options.list_wt:
        wt_data = load_data(options.list_wt)
        print("##Filter variants of wildtypes", file=sys.stderr)
        pipeline(wt_data, options)
        union_wt = GenomicVariantSet(name = "union_wt")
        for wt in wt_data:
            union_wt.sequences += wt.sequences 
        
        print("#wildtype variants:", file=sys.stderr)
        print("union WT", len(union_wt), file=sys.stderr, sep="\t")
        
        #delete Wildtype
        for sample in sample_data:
            sample.subtract(union_wt)
        
        print_length(sample_data, "#variants after subtracting wildtypes")
    else:
        print("#Do not filter by wildtype", file=sys.stderr)
    
    if options.max_density:
        get_max_density(GenomicVariantSets=sample_data, lowerBound=options.lower_bound, upperBound=options.upper_bound)
    else:
        print("#Do not perform max. density search", file=sys.stderr)
        
    if options.list_bed:
        tfbs_motifs = GenomicRegionSet('tfbs_motifs')   
        tfbs_motifs.read_bed(options.list_bed)
        
        for sample in sample_data:
            sample.intersect(tfbs_motifs)
    
        print_length(sample_data, "#variants after filtering by BED file")
    else:
        print("#Do not filter by BED file", file=sys.stderr)
    
    print("#Compute intersection of sample's subsets (give intersection's name and size)")
    output_intersections(sample_data)
    
    print("#Write filtered sample files")
    for sample in sample_data:
        sample.write_vcf("%s-filtered.vcf" %sample.name)
    
if __name__ == '__main__':
    main()   
