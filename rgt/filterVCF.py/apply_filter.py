#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog

tool describtion

@author: Manuel Allhoff

"""
from __future__ import print_function
#import pydevd;pydevd.settrace()
from optparse import OptionParser
import sys
from rgt.GenomicVariantSet import GenomicVariantSet
from rgt.GenomicRegionSet import GenomicRegionSet
from max_density import AlgGoldwasser


class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def load_data(cluster = False):
    if cluster:
        pass
    else:
        p01 = GenomicVariantSet('/home/manuel/data/humangenetics/01_S1_L001_R1_001.filtered.vcf', name='p01')
        p11 = GenomicVariantSet('/home/manuel/data/humangenetics/11_S2_L001_R1_001.filtered.vcf', name='p11')
        p12 = GenomicVariantSet('/home/manuel/data/humangenetics/12_S3_L001_R1_001.filtered.vcf', name='p12')
        p18 = GenomicVariantSet('/home/manuel/data/humangenetics/18_S4_L001_R1_001.filtered.vcf', name='p18')
        p25 = GenomicVariantSet('/home/manuel/data/humangenetics/25_S5_L001_R1_001.filtered.vcf', name='p25')
        pWT = GenomicVariantSet('/home/manuel/data/humangenetics/K28_S8_L001_R1_001.filtered.vcf', name='pWT')
        
    return p01, p11, p12, p18, p25, pWT

def print_length(l, info):
    names = map(lambda x: x.name, l)
    for i, name in enumerate(names):
        print('Lengths ', info, name, len(l[i]), file=sys.stderr, sep='\t')
    print("", file=sys.stderr)

def get_max_density(GenomicVariantSets, lowerBound=20, upperBound=50, max_it=5):

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
        

if __name__ == '__main__':
#     parser = HelpfulOptionParser(usage=__doc__)
#     
#     parser.add_option("--arg1", dest="arg1", default=None, help="help msg")
#     (options, args) = parser.parse_args()
#     variable1 = args[0]
#     
#     i = 5
#     if len(args) != i:
#         parser.error("Exactly %s parameters are needed" %i)
    
    
    thres_mq = 20
    thres_dp = 20 
    filter_dbSNP = True
    tfbs_motifs_path = '/home/manuel/workspace/cluster_p/human_genetics/exp/exp01_motifsearch_sox2/humangenetics_motifs/Match/chr11_mpbs.bed'
    
    #Load data
    p01, p11, p12, p18, p25, pWT = load_data()
    
    print_length([p01, p11, p12, p18, p25, pWT], "original")

    
    #filter MQ
    p01.filter('MQ', '>=', thres_mq)
    p11.filter('MQ', '>=', thres_mq)
    p12.filter('MQ', '>=', thres_mq)
    p18.filter('MQ', '>=', thres_mq)
    p25.filter('MQ', '>=', thres_mq)
    pWT.filter('MQ', '>=', thres_mq)
 
    print_length([p01, p11, p12, p18, p25, pWT], "after keeping MQ >= %s" %thres_mq)
     
    #filter DP
    p01.filter('DP', '>=', thres_dp)
    p11.filter('DP', '>=', thres_dp)
    p12.filter('DP', '>=', thres_dp)
    p18.filter('DP', '>=', thres_dp)
    p25.filter('DP', '>=', thres_dp)
    pWT.filter('DP', '>=', thres_dp)
 
    print_length([p01, p11, p12, p18, p25, pWT], "after keeping DP >= %s" %thres_dp)

#     #filter dbSNP
#     if filter_dbSNP:
#         p01.filter_dbSNP()
#         p11.filter_dbSNP()
#         p12.filter_dbSNP()
#         p18.filter_dbSNP()
#         p25.filter_dbSNP()
#         pWT.filter_dbSNP()
# 
#     print_length([p01, p11, p12, p18, p25, pWT], "after filtering with dbSNP >= %s" %thres_dp)
    
    
    #get_max_density(GenomicVariantSets=[p01, p11, p12, p18, p25], lowerBound=15000, upperBound=35000)
    
    #delete Wildtype
    p01.subtract(pWT)
    p11.subtract(pWT)
    p12.subtract(pWT)
    p18.subtract(pWT)
    p25.subtract(pWT)
 
    print_length([p01, p11, p12, p18, p25, pWT], "after subtracting Wildtype")
    
    #TFBS sides
    tfbs_motifs = GenomicRegionSet('tfbs_motifs')   
    tfbs_motifs.read_bed(tfbs_motifs_path)
    
    p01.intersect(tfbs_motifs)
    p11.intersect(tfbs_motifs)
    p12.intersect(tfbs_motifs)
    p18.intersect(tfbs_motifs)
    p25.intersect(tfbs_motifs)
    
    print_length([p01, p11, p12, p18, p25, pWT], "after considering TBFS motifs")
    
    p01.write_vcf('p01.final.vcf')
    p11.write_vcf('p11.final.vcf')
    p12.write_vcf('p12.final.vcf')
    p18.write_vcf('p18.final.vcf')
    p25.write_vcf('p25.final.vcf')
    pWT.write_vcf('pWT.final.vcf')
    
    
