#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog <BAM> <BAM> <FASTA>

Find differential peaks between two <BAM> files in <FASTA> genome.

Author: Manuel Allhoff (allhoff@aices.rwth-aachen.de)

"""

from __future__ import print_function
from optparse import OptionParser
from .. CoverageSet import CoverageSet
from .. GenomicRegion import GenomicRegion
from .. GenomicRegionSet import GenomicRegionSet
from .. Util import GenomeData

from HTSeq import FastaReader
from get_extension_size import get_extension_size
import os.path
import sys
from DualCoverageSet import DualCoverageSet
from get_fast_gen_pvalue import get_log_pvalue_new

SIGNAL_CUTOFF = 10000

def dump_posteriors_and_viterbi(name, posteriors, DCS, states):
    indices_of_interest = DCS.indices_of_interest
    first_overall_coverage = DCS.first_overall_coverage
    second_overall_coverage = DCS.second_overall_coverage
    
    c1 = list(first_overall_coverage)
    c2 = list(second_overall_coverage)
    
    #print("Computing info...", file=sys.stderr)
    f = open(name + '-posts.bed', 'w')
    g = open(name + '-states-viterbi.bed', 'w')
    
    for i in range(len(indices_of_interest)):
        cov1 = c1[indices_of_interest[i]]
        cov2 = c2[indices_of_interest[i]]
        
        p1, p2, p3 = posteriors[i][0], posteriors[i][1], posteriors[i][2]
        chrom, start, end = DCS._index2coordinates(indices_of_interest[i])
        
        print(chrom, start, end, states[i], cov1, cov2, sep = '\t', file=g)
        print(chrom, start, end, max(p3, max(p1,p2)), p1, p2, p3, cov1, cov2, sep = '\t', file=f)

    f.close()
    g.close()


def _compute_pvalue((x, y, side, distr)):
    if x == 'NA':
        return sys.maxint
    else:
        return -get_log_pvalue_new(x, y, side, distr)
    
def get_peaks(name, DCS, states, distr):
    indices_of_interest = DCS.indices_of_interest
    first_overall_coverage = DCS.first_overall_coverage
    second_overall_coverage = DCS.second_overall_coverage
    
    c1 = list(first_overall_coverage)
    c2 = list(second_overall_coverage)
    
    
    tmp_peaks = []
    
    for i in range(len(indices_of_interest)):
        if states[i] not in [1,2]:
            continue #ignore background states
        
        strand = '+' if states[i] == 1 else '-'
        
        cov1 = c1[indices_of_interest[i]]
        cov2 = c2[indices_of_interest[i]]
        chrom, start, end = DCS._index2coordinates(indices_of_interest[i])
        
        tmp_peaks.append((chrom, start, end, cov1, cov2, strand))

    i, j, r = 0, 0, 0
    
    peaks = []
    pvalues = []
    
    while i < len(tmp_peaks):
        j+=1
        c, s, e, c1, c2, strand = tmp_peaks[i]
        v1 = [c1]
        v2 = [c2]
        
        #merge bins
        while i+1 < len(tmp_peaks) and e == tmp_peaks[i+1][1] and strand == tmp_peaks[i+1][5]:
            e = tmp_peaks[i+1][2]
            v1.append(tmp_peaks[i+1][3])
            v2.append(tmp_peaks[i+1][4])
            i += 1
        
        s1 = sum(v1)
        s2 = sum(v2)

        if s1 + s2 > SIGNAL_CUTOFF:
            pvalues.append(('NA', 'NA', 'NA', 'NA'))
        else:
            if strand == '+':
                pvalues.append((s1, s2, 'l', distr))
            else:
                pvalues.append((s1, s2, 'r', distr))

        peaks.append((c, s, e, s1, s2, strand))
        i += 1
    
    print('Number of Peaks where p-value is not calculated: ', pvalues.count(('NA', 'NA', 'NA', 'NA')), file=sys.stderr)
    
    #pool = multiprocessing.Pool(processes=2)#multiprocessing.cpu_count() * 3/2)
    pvalues = map(_compute_pvalue, pvalues)
    
    assert len(pvalues) == len(peaks)
    
    f = open(name + '-diffpeaks.bed', 'w')
    
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000
    
    for i in range(len(pvalues)):
        c, s, e, c1, c2, strand = peaks[i]
        color = colors[strand]

        print(c, s, e, 'Peak' + str(i), bedscore, strand, s, e, \
              color, 0, str(c1) + ',' + str(c2) + ',' + str(pvalues[i]), sep='\t', file=f)

    f.close()


def initialize(name, genome_path, regions, stepsize, binsize, bam_file_1, bam_file_2, ext_1, ext_2, \
               input_1, input_factor_1, ext_input_1, input_2, input_factor_2, ext_input_2, chrom_sizes, verbose, norm_strategy, no_gc_content, deadzones):
    
    regionset = GenomicRegionSet(name)
    
    #if regions option is set, take the values, otherwise the whole set of 
    #chromosomes as region to search for DPs
    if regions is not None:
        with open(regions) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                c, s, e = line[0], int(line[1]), int(line[2])
                regionset.add(GenomicRegion(chrom=c, initial=s, final=e))
    else:
        chromosomes = []
        for s in FastaReader(genome_path):
            chromosomes.append(s.name)
        chromosomes.sort()
        print(chromosomes, file=sys.stderr)
        with open(chrom_sizes) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                chr, end = line[0], int(line[1])
                if chr in chromosomes:
                    regionset.add(GenomicRegion(chrom=chr, initial=0, final=end))
    
    regionset.sequences.sort()
    
    start = 0
    end = 600
    ext_stepsize = 5
    #TODO: maybe for-loops?
    #compute extension size
    if [ext_1, ext_2, ext_input_1, ext_input_2].count(None) > 0:
        print("Computing read extension sizes...", file=sys.stderr)
    if ext_1 is None:
        ext_1, values_1 = get_extension_size(bam_file_1, start=start, end=end, stepsize=ext_stepsize)
    print("Read extension for first file: %s" %ext_1, file=sys.stderr)
    
    if ext_2 is None:
        ext_2, values_2 = get_extension_size(bam_file_2, start=start, end=end, stepsize=ext_stepsize)
    print("Read extension for second file: %s" %ext_2, file=sys.stderr)

    if input_1 is not None and ext_input_1 is None:
        ext_input_1, values_input_1 = get_extension_size(input_1, start=start, end=end, stepsize=ext_stepsize)
    print("Read extension for first input file: %s" %ext_input_1, file=sys.stderr)
    
    if input_1 is not None and input_2 is not None and input_1 == input_2 and 'ext_input_1' in locals() and 'values_input_1' in locals():
        ext_input_2, values_input_2 = ext_input_1, values_input_1
    elif input_2 is not None and ext_input_2 is None:
        ext_input_2, values_input_2 = get_extension_size(input_2, start=start, end=end, stepsize=ext_stepsize)
    print("Read extension for second input file: %s" %ext_input_2, file=sys.stderr)
    
    if verbose:
        if 'values_1' in locals() and values_1 is not None:
            with open(name + '-read-ext-1', 'w') as f:
                for v, i in values_1:
                    print(i, v, sep='\t', file=f)
        
        if 'values_2' in locals() and values_2 is not None:
            with open(name + '-read-ext-2', 'w') as f:
                for v, i in values_2:
                    print(i, v, sep='\t', file=f)
        
        if 'values_input_1' in locals() and values_input_1 is not None:
            with open(name + '-read-ext-input-1', 'w') as f:
                for v, i in values_input_1:
                    print(i, v, sep='\t', file=f)
        
        if 'values_input_2' in locals() and values_input_2 is not None:
            with open(name + '-read-ext-input-2', 'w') as f:
                for v, i in values_input_2:
                    print(i, v, sep='\t', file=f)

    cov_cdp_mpp = DualCoverageSet(name=name, region=regionset, genome_path=genome_path, binsize=binsize, stepsize=stepsize,rmdup=True,\
                                  file_1=bam_file_1, ext_1=ext_1,\
                                  file_2=bam_file_2, ext_2=ext_2, \
                                  input_1=input_1, ext_input_1=ext_input_1, input_factor_1=input_factor_1, \
                                  input_2=input_2, ext_input_2=ext_input_2, input_factor_2=input_factor_2, \
                                  chrom_sizes=chrom_sizes, verbose=verbose, norm_strategy=norm_strategy, no_gc_content=no_gc_content, deadzones=deadzones)
    
    return cov_cdp_mpp


def _get_chrom_list(file):
    l = []
    with open(file) as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            if line[0] not in l:
                l.append(line[0])
    return l

# def _check_order(deadzones, regions, parser):
#     chrom_dz = _get_chrom_list(deadzones)
#     chrom_regions= _get_chrom_list(regions)
#     #chrom_regions may be subset of chrom_dz, but in same ordering
#     pos_old = -1
#     tmp_dz = []
#     #the list should have the same element
#     for r_chrom in chrom_dz:
#         if r_chrom in chrom_regions:
#             tmp_dz.append(r_chrom)
#     #they should be in the same order
#     if not tmp_dz == chrom_regions:
#         parser.error("Please make sure the deadzone file has the same order as the region file.")
    

class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def input(laptop):
    parser = HelpfulOptionParser(usage=__doc__)
    
    parser.add_option("--input-1", dest="input_1", default=None, \
                      help="Input control file for first parameter [default: %default]")
    parser.add_option("--input-2", dest="input_2", default=None, \
                      help="Input control file for second parameter [default: %default]")
    parser.add_option("-p", "--pvalue", dest="pcutoff", default=0.01, type="float",\
                      help="Call only peaks with p-value lower than cutoff. [default: %default]")
    parser.add_option("-b", "--binsize", dest="binsize", default=100, type="int",\
                      help="Size of underlying bins for creating the signal.  [default: %default]")
    parser.add_option("-s", "--step", dest="stepsize", default=50, type="int",\
                      help="Stepsize with which the window consecutively slides across the genome to create the HMM signal.")
    parser.add_option("-n", "--name", default=None, dest="name", type="string",\
                      help="Experiment's name and prefix for all files that are created.")
    parser.add_option("--ext-1", default=None, dest="ext_1", type="int",\
                      help="Read's extension size for first BAM file. If None, estimate extension size. [default: %default]")
    parser.add_option("--ext-2", default=None, dest="ext_2", type="int",\
                      help="Read's extension size for second BAM file. If None, estimate extension size [default: %default]")
    parser.add_option("--ext-input-1", default=None, dest="ext_input_1", type="int",\
                      help="Read's extension size for first input file. If None, estimate extension size. [default: %default]")
    parser.add_option("--ext-input-2", default=None, dest="ext_input_2", type="int",\
                      help="Read's extension size for second input file. If None, estimate extension size. [default: %default]")
    parser.add_option("--factor-input-1", default=None, dest="input_factor_1", type="float",\
                      help="Normalization factor for first input. If None, estimate factor. [default: %default]")
    parser.add_option("--factor-input-2", default=None, dest="input_factor_2", type="float",\
                      help="Normalization factor for first input. If None, estimate factor. [default: %default]")
    parser.add_option("-v", "--verbose", default=False, dest="verbose", action="store_true", \
                      help="Output among others initial state distribution, putative differential peaks and genomic signal. [default: %default]")
    parser.add_option("--version", dest="version", default=False, action="store_true", help="Show script's version.")
    #parser.add_option("--norm-strategy", dest="norm_strategy", default=5, type="int", help="1: naive; 2: Diaz; 3: own; 4: Diaz and own; 5: diaz and naive")
    parser.add_option("--no-gc-content", dest="no_gc_content", default=False, action="store_true", \
                      help="Do not compute GC content model  [default: %default]")
    parser.add_option("-r", "--regions", default=None, dest="regions", \
                     help="Regions where to search for differential peaks.  [default: %default]")
    parser.add_option("-c", "--chrom-sizes", default=None, dest="chrom_sizes", \
                      help="Sizes of chromosomes.  [default: %default]")
#     parser.add_option("--deadzones", dest="deadzones", default=None, \
#                       help="Deadzones (BED) [default: %default]")
    
    (options, args) = parser.parse_args()
    options.norm_strategy = 5
    
    if options.version:
        version = "version \"0.01alpha\""
        print("")
        print(version)
        sys.exit()
    
    if len(args) != 3:
        parser.error("Exactly three parameters are needed")
        
    bamfile_1 = args[0]
    bamfile_2 = args[1]
    genome = args[2]
    
    if options.chrom_sizes is not None and not os.path.isfile(options.chrom_sizes):
        parser.error("Please give a proper file for chromosome files.")
    
    if options.chrom_sizes is None:
        options.chrom_sizes = GenomeData.CHROMOSOME_SIZES
    
    if options.regions is not None and not os.path.isfile(options.regions):
        parser.error("Please give a proper region file.")
    
    if not os.path.isfile(bamfile_1) or not os.path.isfile(bamfile_2):
        parser.error("At least one input parameter is not a file")
    
    if options.name is None:
        prefix = os.path.splitext(os.path.basename(bamfile_1))[0]
        suffix = os.path.splitext(os.path.basename(bamfile_2))[0]
        options.name = "-".join(['exp', prefix, suffix])
    
    if (options.input_1 is None and options.ext_input_1 is not None) \
        or (options.input_2 is None and options.ext_input_2 is not None):
        parser.error("Read extension size without input file (use -i)")
    
    if options.input_1 is None and options.input_2 is None:
        print("GC content is not calculated as there is no input file.", file=sys.stderr)
        
    if options.norm_strategy in [2, 4, 5] and (options.input_1 is None or options.input_2 is None):
        parser.error("Please define input files for this normalization strategy!")
        
    if options.norm_strategy is not None and (options.input_factor_1 is not None or options.input_factor_1 is not None ):
        parser.error("Input factors are not allowed for this normalization strategy!")
    
    if not options.no_gc_content and (options.input_1 is None or options.input_2 is None):
        parser.error("GC content can only be computed with both input files.")
    
#     if options.deadzones is not None:
#         #check the ordering of deadzones and region
#         _check_order(options.deadzones, regions, parser)
    
    return options, bamfile_1, bamfile_2, genome

