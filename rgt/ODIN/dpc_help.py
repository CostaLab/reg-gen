#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
%prog <BAM> <BAM> <FASTA> <CHROM SIZES>

Find differential peaks between two <BAM> files in <FASTA> genome.

Author: Manuel Allhoff (allhoff@aices.rwth-aachen.de)

"""

from __future__ import print_function
from optparse import OptionParser, OptionGroup
from rgt.CoverageSet import CoverageSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet

from HTSeq import FastaReader
from get_extension_size import get_extension_size
import os.path
import sys
from DualCoverageSet import DualCoverageSet
from get_fast_gen_pvalue import get_log_pvalue_new
from postprocessing import merge_delete
from math import log10
from rgt.motifanalysis.Statistics import multiple_test_correction
import os

SIGNAL_CUTOFF = 30000

def get_bibtex_entry():
    print("@article{Allhoff2014ODIN,", file=sys.stderr)
    print("    author = {Allhoff, Manuel and Seré, Kristin and Chauvistré,", file=sys.stderr)
    print("              Heike and Lin, Qiong and Zenke, Martin and Costa, Ivan G.},", file=sys.stderr)
    print("    title = {Detecting differential peaks in ChIP-seq signals with ODIN},", file=sys.stderr)
    print("    journal = {Bioinformatics},", file=sys.stderr)
    print("    year = {2014},", file=sys.stderr)
    print("    volume = {30}, ", file=sys.stderr)
    print("    number = {24}, ", file=sys.stderr)
    print("    pages = {3467-3475}, ", file=sys.stderr)
    print("    doi = {10.1093/bioinformatics/btu722}, ", file=sys.stderr)
    print("    URL = {http://bioinformatics.oxfordjournals.org/content/30/24/3467}", file=sys.stderr)
    print("}", file=sys.stderr)

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

def _output_BED(name, pvalues, peaks, pv_pass):
    f = open(name + '-diffpeaks.bed', 'w')
     
    colors = {'+': '255,0,0', '-': '0,255,0'}
    bedscore = 1000
    
    for i in range(len(pvalues)):
        c, s, e, c1, c2, strand = peaks[i]
        color = colors[strand]
            
        if pv_pass[i]:
            print(c, s, e, 'Peak' + str(i), bedscore, strand, s, e, \
                  color, 0, str(c1) + ',' + str(c2) + ',' + str(pvalues[i]), sep='\t', file=f)
    
    f.close()

def _output_narrowPeak(name, pvalues, peaks, pv_pass):
    """Output in narrowPeak format,
    see http://genome.ucsc.edu/FAQ/FAQformat.html#format12"""
    f = open(name + '-diffpeaks.narrowPeak', 'w')
    for i in range(len(pvalues)):
        c, s, e, _, _, strand = peaks[i]
        if pv_pass[i]:
            print(c, s, e, 'Peak' + str(i), 0, strand, 0, pvalues[i], 0, -1, sep='\t', file=f)
    f.close()

def get_peaks(name, DCS, states, ext_size, merge, distr, pcutoff, no_correction):
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

    i, j = 0, 0
    
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
    
    merge_delete(ext_size, merge, peaks, pvalues, name)
    
    #peaks = [(c, s, e, s1, s2, strand)]
    if not no_correction:
        pvalues = map(lambda x: 10**-x, pvalues)
        pv_pass, pvalues = multiple_test_correction(pvalues, alpha=pcutoff)
        pvalues = map(_get_log10pvalue, pvalues)
    else:
        pv_pass = [True] * len(pvalues)
    
    _output_BED(name, pvalues, peaks, pv_pass)
    _output_narrowPeak(name, pvalues, peaks, pv_pass)

def _get_log10pvalue(x):
    """Return -log10(p)"""
    if x < 1e-320:
        return sys.maxint
    else:
        return -log10(x)

def initialize(name, genome_path, regions, stepsize, binsize, bam_file_1, bam_file_2, ext_1, ext_2, \
               input_1, input_factor_1, ext_input_1, input_2, input_factor_2, ext_input_2, chrom_sizes, verbose, norm_strategy, no_gc_content, deadzones,\
               factor_input_1, factor_input_2, debug, tracker):
    
    regionset = GenomicRegionSet(name)
    chrom_sizes_dict = {}
    #if regions option is set, take the values, otherwise the whole set of 
    #chromosomes as region to search for DPs
    if regions is not None:
        with open(regions) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                c, s, e = line[0], int(line[1]), int(line[2])
                regionset.add(GenomicRegion(chrom=c, initial=s, final=e))
                chrom_sizes_dict[c] = e
    else:
        with open(chrom_sizes) as f:
            for line in f:
                line = line.strip()
                line = line.split('\t')
                chrom, end = line[0], int(line[1])
                regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))
                chrom_sizes_dict[chrom] = end
    
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
    
    tracker.write(text=str(ext_1) + "," + str(ext_2), header="Extension size IP1, IP2")
    if input_1 is not None and input_2 is not None:
        tracker.write(text=str(ext_input_1) + "," + str(ext_input_2), header="Extension size Control1, Control2")
    
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
                                  chrom_sizes=chrom_sizes, verbose=verbose, norm_strategy=norm_strategy, no_gc_content=no_gc_content, deadzones=deadzones,\
                                  factor_input_1=factor_input_1, factor_input_2=factor_input_2, chrom_sizes_dict=chrom_sizes_dict, debug=debug, tracker=tracker)
    
    return cov_cdp_mpp, [ext_1, ext_2]


class HelpfulOptionParser(OptionParser):
    """An OptionParser that prints full help on errors."""
    def error(self, msg):
        self.print_help(sys.stderr)
        self.exit(2, "\n%s: error: %s\n" % (self.get_prog_name(), msg))

def which(program):
    """Return path of program or None, see
    http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python"""
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)
    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None

def input(test):
    parser = HelpfulOptionParser(usage=__doc__)
    
    if test:
        (options, args) = parser.parse_args()
#        bamfile_1 = '/home/manuel/data/project_chipseq_norm/data/PU1_CDP_1000000.bam'
        #bamfile_1 = '/home/manuel/data/project_chipseq_norm/data/PU1_CDP_10k.bam'
#        bamfile_1 = '/home/manuel/data/project_chipseq_norm/data/PU1_MPP_chr1-2.bam'
#        bamfile_1 = '/home/manuel/data/project_chipseq_norm/data/PU1_CDP_1m.bam'
        bamfile_1 = '/home/manuel/data/project_chipseq_norm/data/pu1_mpp_1-12-18.bam'
#        bamfile_2 = '/home/manuel/data/project_chipseq_norm/data/PU1_MPP_100000.bam'
       # bamfile_2 = '/home/manuel/data/project_chipseq_norm/data/PU1_MPP_10k.bam'
#        bamfile_2 = '/home/manuel/data/project_chipseq_norm/data/PU1_CDP_chr1-2.bam'
        bamfile_2 = '/home/manuel/data/project_chipseq_norm/data/pu1_cdp_1-12-18.bam'

        #options.regions = '/home/manuel/data/mm9_features/mm9.extract.sizes'
        options.merge=True
        options.mag = 3
        options.regions = None
        genome = '/home/manuel/data/mm/mm9/mm9.fa'
        options.ext_1 = None
        options.ext_2 = None
        options.ext_input_1 = None
        options.ext_input_2 = None
        options.input_2 = None #'/home/manuel/data/project_chipseq_norm/data/PU1_Input_10k.bam'
        options.input_1 = None #'/home/manuel/data/project_chipseq_norm/data/PU1_Input_10k.bam'
        options.confidence_threshold=0.7
        options.foldchange=1.05
        options.pcutoff = 0.1
        options.name='test'
        options.distr='binom'
        options.constchrom = None #'chr1'
        options.stepsize=50
        options.binsize=100
        options.input_factor_1= None #0.7
        options.input_factor_2= None #0.7
        options.norm_strategy = 5
        options.verbose=False
        options.debug=False
        #chrom_sizes='/home/manuel/data/mm/mm9/mm9.chrom.sizes'
        chrom_sizes = '/home/manuel/workspace/cluster_p/genomes/hg/hg19.sizes'
        options.no_gc_content = False
        options.deadzones = None #"/home/manuel/dz.bed"
        options.version=False
        options.factor_input_1=None #for BAM
        options.factor_input_2=None
        options.no_correction = False
        
        if options.debug:
            options.verbose = True
        print(options, file=sys.stderr)
    else:
        parser.add_option("--input-1", dest="input_1", default=None, \
                          help="Input control for first parameter [default: %default]")
        parser.add_option("--input-2", dest="input_2", default=None, \
                          help="Input control for second parameter [default: %default]")
        parser.add_option("-p", "--pvalue", dest="pcutoff", default=0.1, type="float",\
                          help="p-value cutoff: call only peaks with p-value lower than cutoff [default: %default]")
        parser.add_option("--no-correction", default=False, dest="no_correction", action="store_true", \
                          help="No Benjamini/Hochberg p-value multiple testing correction [default: %default]")
        parser.add_option("-m", "--merge", default=False, dest="merge", action="store_true", \
                          help="Merge peaks which have a distance less than the estimated fragment size (recommended for histone data). [default: %default]")
        parser.add_option("-n", "--name", default=None, dest="name", type="string",\
                          help="Experiment's name and prefix for all files that are created.")
        parser.add_option("--ext-1", default=None, dest="ext_1", type="int",\
                          help="Read's extension size for first BAM file. If option is not chosen, estimate extension size. [default: %default]")
        parser.add_option("--ext-2", default=None, dest="ext_2", type="int",\
                          help="Read's extension size for second BAM file. If option is not chosen, estimate extension size [default: %default]")
        parser.add_option("--ext-input-1", default=None, dest="ext_input_1", type="int",\
                          help="Read's extension size for first input file. If option is not chosen, estimate extension size. [default: %default]")
        parser.add_option("--ext-input-2", default=None, dest="ext_input_2", type="int",\
                          help="Read's extension size for second input file. If option is not chosen, estimate extension size. [default: %default]")
        #parser.add_option("--factor-1", default=None, dest="factor_input_1", type="float",\
        #                  help="Factor for first BAM. [default: %default]")
        #parser.add_option("--factor-2", default=None, dest="factor_input_2", type="float",\
        #                  help="Factor for second BAM. [default: %default]")
        
        parser.add_option("--distr", default="binom", dest="distr", type="string", \
                          help="HMM's emission distribution. [Binomial (binom), (constraint) mixture of poisson (poisson-c)] [default: %default]")
        parser.add_option("--mag", default=3, dest="mag", type="int", help="Magnitude of Poisson mixture model. [default: %default]")
        
        parser.add_option("-v", "--verbose", default=False, dest="verbose", action="store_true", \
                          help="output further information of DP-Calling progress [default: %default]")
        parser.add_option("--version", dest="version", default=False, action="store_true", help="show version [default: %default]")
        parser.add_option("--bibtex", dest="bibtex", default=False, action="store_true", help="show BibTeX entry [default: %default]")
        #parser.add_option("--norm-strategy", dest="norm_strategy", default=5, type="int", help="1: naive; 2: Diaz; 3: own; 4: Diaz and own; 5: diaz and naive")
        parser.add_option("--output-dir", dest="outputdir", default=None, type="string", \
                          help="All files are stored in output directory which is created if necessary.")
        
        group = OptionGroup(parser, "Advanced options")
        group.add_option("--regions", dest="regions", default=None, help="regions (BED) where to search for DPs [default: entire genome]")
        group.add_option("--deadzones", dest="deadzones", default=None, help="Deadzones (BED) [default: %default]")
        group.add_option("--no-gc-content", dest="no_gc_content", default=False, action="store_true", \
                          help="turn of GC-content calculation (faster, but less accurate) [default: %default]")
        group.add_option("--const-chrom", default=None, dest="constchrom", type="string",\
                          help="Constrain HMM learning process to chromosome. [default: %default]")
        group.add_option("--factor-input-1", default=None, dest="input_factor_1", type="float",\
                          help="Normalization factor for first input. If option is not chosen, estimate factor. [default: %default]")
        group.add_option("--factor-input-2", default=None, dest="input_factor_2", type="float",\
                          help="Normalization factor for second input. If option is not chosen, estimate factor. [default: %default]")
        group.add_option("-c", "--confidence_threshold", dest="confidence_threshold", default=0.7, type="float",\
                          help="Threshold that each observation's posterior probability must exceed to be considered as a differential peak. [default: %default]")
        group.add_option("-f", "--foldchange", default=1.05, dest="foldchange", type="float",\
                          help="Minimum fold change which a potential differential peak must exhibit. [default: %default]")
        group.add_option("-b", "--binsize", dest="binsize", default=100, type="int",\
                          help="Size of underlying bins for creating the signal [default: %default]")
        group.add_option("-s", "--step", dest="stepsize", default=50, type="int",\
                          help="Stepsize with which the window consecutively slides across the genome to create the signal. [default: %default]")
        group.add_option("--debug", default=False, dest="debug", action="store_true", \
                          help="Output debug information. Warning: space consuming! [default: %default]")

        parser.add_option_group(group)
        
        (options, args) = parser.parse_args()
        
        options.norm_strategy = 5 #get rid of other options, this is an ugly but efficient solution
        options.factor_input_1 = None #supposed to vary the BAM input with a predefined factor
        options.factor_input_2 = None
        
        if options.debug:
            options.verbose = True
        
        if options.version:
            version = "version \"0.1alpha\""
            print("")
            print(version)
            sys.exit()
        
        if options.bibtex:
            get_bibtex_entry()
            sys.exit(0)
        
        if len(args) != 4:
            parser.error("Exactly four parameters are needed")
            
        bamfile_1 = args[0]
        bamfile_2 = args[1]
        genome = args[2]
        chrom_sizes = args[3]

    
    if not os.path.isfile(bamfile_1) or not os.path.isfile(bamfile_2) \
        or (options.regions is not None and not os.path.isfile(options.regions)) or not os.path.isfile(genome):
        parser.error("At least one input parameter is not a file")
    
    if (options.input_1 is not None and not os.path.isfile(options.input_1)) or \
        (options.input_2 is not None and not os.path.isfile(options.input_2)):
        parser.error("At least one control input is not a file")
    
    if options.name is None:
        prefix = os.path.splitext(os.path.basename(bamfile_1))[0]
        suffix = os.path.splitext(os.path.basename(bamfile_2))[0]
        options.name = "-".join(['exp', prefix, suffix])
    
    if (options.input_1 is None and options.ext_input_1 is not None) \
        or (options.input_2 is None and options.ext_input_2 is not None):
        parser.error("Read extension size without input file (use -i)")
    
    if options.input_1 is None and options.input_2 is None:
        print("GC content is not calculated as there is no input file.", file=sys.stderr)
        options.no_gc_content = True
    
    if options.input_1 is None and options.input_2 is None:
        options.norm_strategy = 1
    
    if options.norm_strategy in [2, 4, 5] and (options.input_1 is None or options.input_2 is None):
        parser.error("Please define input files for this normalization strategy!")
        
    #if options.norm_strategy is not None and (options.input_factor_1 is not None or options.input_factor_1 is not None ):
    #    parser.error("Input factors are not allowed for this normalization strategy!")
    
    if not options.no_gc_content and (options.input_1 is None or options.input_2 is None):
        parser.error("GC content can only be computed with both input files.")
    
#     if options.deadzones is not None:
#         #check the ordering of deadzones and region
#         _check_order(options.deadzones, regions, parser)
    
    if not which("wigToBigWig"):
        print("Warning: wigToBigWig programm not found! Signal will not be stored!", file=sys.stderr)
    
    
    if options.outputdir:
        options.outputdir = os.path.expanduser(options.outputdir) #replace ~ with home path
        if os.path.isdir(options.outputdir) and sum(map(lambda x: x.startswith(options.name), os.listdir(options.outputdir))) > 0:
            parser.error("Output directory exists and contains files with names starting with your chosen experiment name! Do nothing to prevend file overwriting!")
        if not os.path.exists(options.outputdir):
            os.mkdir(options.outputdir)
    else:
        options.outputdir = os.getcwd() 
    
    options.name = os.path.join(options.outputdir, options.name)
    
    return options, bamfile_1, bamfile_2, genome, chrom_sizes
