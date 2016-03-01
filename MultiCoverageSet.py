"""
MultiCoverageSet
===================
MultiCoverageSet represents a set of CoverageSet.

"""

from __future__ import print_function
from rgt.CoverageSet import CoverageSet
import numpy as np
from random import sample
from rgt.ODIN.gc_content import get_gc_context
import sys
from rgt.ODIN.normalize import get_normalization_factor
from math import fabs
from rgt.ODIN.DualCoverageSet import DualCoverageSet
from rgt.GenomicRegionSet import GenomicRegionSet
from copy import deepcopy
import gc
from math import fabs
from rgt.THOR.norm_genelevel import norm_gene_level

EPSILON = 1**-320
ROUND_PRECISION = 3
DEBUG = None
VERBOSE = None

class MultiCoverageSet:
    """*Keyword arguments:*

        - name -- name.
        - dims -- tuple of dimension of CoverageSet of first and second condition
	- regions -- RegionSet the signal is defined on
	- genome_path -- path to genome
	- binsize -- size of genomic bins
	- stepsize -- size of steps
	- chrom_sizes -- file with chromosome sizes
	- norm_regionset -- RegionSet the normalization approach is restricted to
	- verbose -- boolean
	- debug -- boolean
	- no_gc_content -- boolean, do not normalize against GC content
	- rmdup - boolean, remove duplicates in CoverageSet
	- path_bamfiles - list of bamfiles
	- exts - extension sizes of bam files
	- path_input - list of input-dna bam files
	- exts_inputs - extension sizes of input-dna bam files
	- factors_input - predefined normalization factors for input-dna files
	- chrom_sizes_dict - dict of chromosome sizes
	- scaling_factors_ip - predefined normalization factors for IP channel
	- save_wig - boolean
	- strand_cov - save strand information
	- housekeeping_genes - normalize ip channel with control regions
	- tracker - tracker object
	- end - deprecated
	- counter - deprecated
	- gc_content_cov - predefined values for GC content normalization
	- avg_gc_content - predefined values for GC content normalization
	- gc_hist - predefined values for GC content normalization
	- output_bw - output bigwig files after initialization
	- folder - output folder
	- report - html report
    """

    def _help_init(self, path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, dim, regions, norm_regionset, strand_cov):
        """Return self.covs and self.inputs as CoverageSet"""
        self.exts = exts
        self.covs = [CoverageSet('file' + str(i), regions) for i in range(dim)]
        for i, c in enumerate(self.covs):
            c.coverage_from_bam(bam_file=path_bamfiles[i], read_size=exts[i], rmdup=rmdup, binsize=binsize,\
                                stepsize=stepsize, get_strand_info = strand_cov)
        self.covs_avg = [CoverageSet('cov_avg'  + str(i) , regions) for i in range(2)]
        if path_inputs:
            self.inputs = [CoverageSet('input' + str(i), regions) for i in range(len(path_inputs))]
            for i, c in enumerate(self.inputs):
                c.coverage_from_bam(bam_file=path_inputs[i], read_size=exts_inputs[i], rmdup=rmdup, binsize=binsize,\
                                stepsize=stepsize, get_strand_info = strand_cov)
            self.input_avg = [CoverageSet('input_avg'  + str(i), regions) for i in range(2)]
        else:
            self.inputs = []
            
        if norm_regionset:
            self.norm_regions = [CoverageSet('norm_region' + str(i), norm_regionset) for i in range(dim)]
            for i, c in enumerate(self.norm_regions):
                c.coverage_from_bam(bam_file=path_bamfiles[i], read_size=exts[i], rmdup=rmdup, binsize=binsize,\
                                    stepsize=stepsize, get_strand_info = strand_cov)
            self.input_avg = [CoverageSet('input_avg'  + str(i), regions) for i in range(2)]
        else:
            self.norm_regions = None
    
    def _compute_gc_content(self, no_gc_content, path_inputs, stepsize, binsize, genome_path, name, chrom_sizes, chrom_sizes_dict):
        """Compute GC-content"""
        if not no_gc_content and path_inputs and self.gc_content_cov is None:
            print("Compute GC-content", file=sys.stderr)
            for i, cov in enumerate(self.covs):
                inputfile = self.inputs[i] #1 to 1 mapping between input and cov
                rep = i if i < self.dim_1 else i-self.dim_1
                sig = 1 if i < self.dim_1 else 2
                self.gc_content_cov, self.avg_gc_content, self.gc_hist = get_gc_context(stepsize, binsize, genome_path, inputfile.coverage, chrom_sizes_dict)
                self._norm_gc_content(cov.coverage, self.gc_content_cov, self.avg_gc_content)
                self._norm_gc_content(inputfile.coverage, self.gc_content_cov, self.avg_gc_content)
            
    def output_bw(self, name, chrom_sizes, save_wig):
        """Output bigwig file.
        
        *Keyword arguments:*
        
        - name -- filename
        - chrom_size -- file with chromosome sizes
        - save_wig -- BOOLEAN
        
        .. note::
            save_wig may produce large filess.
        """
        
        for i in range(len(self.covs)):
            rep = i if i < self.dim_1 else i-self.dim_1
            sig = 1 if i < self.dim_1 else 2
            self.covs[i].write_bigwig(name + '-' + str(self.counter) + '-s%s-rep%s.bw' %(sig, rep), chrom_sizes, save_wig=save_wig, end=self.end)
        
        self.covs_avg = None
        self.input_avg = None
        if self.inputs:
            for i in range(len(self.covs)):
                self.inputs[i] = None #last time that we need this information, delete it
        gc.collect()
        
    
    def _help_get_data(self, i, type):
        if type != 'normregion':
            for j in range(len(self.covs[i].genomicRegions)):
                if type == 'cov':
                    yield self.covs[i].coverage[j]
                elif type == 'strand':
                    yield self.covs[i].cov_strand_all[j]
        elif type == 'normregion':
            for j in range(len(self.norm_regions[i].genomicRegions)):
                yield self.norm_regions[i].coverage[j]
    
    def _help_init_overall_coverage(self, cov_strand=True):
        """Convert coverage data (and optionally strand data) to matrix list"""
        
        tmp = [[], []]
        tmp2 = [[[], []], [[], []]]
        
        for k in range(2):
            it = range(self.dim_1) if k == 0 else range(self.dim_1, self.dim_1 + self.dim_2)
            for i in it:
                if cov_strand:
                    tmp_el = reduce(lambda x,y: np.concatenate((x,y)), self._help_get_data(i, 'cov'))
                    tmp[k].append(tmp_el)
                 
                    tmp_el = map(lambda x: (x[0], x[1]), reduce(lambda x,y: np.concatenate((x,y)), self._help_get_data(i, 'strand')))
                    tmp2[k][0].append(map(lambda x: x[0], tmp_el))
                    tmp2[k][1].append(map(lambda x: x[1], tmp_el))
                else:
                    tmp_el = reduce(lambda x,y: np.concatenate((x,y)), self._help_get_data(i, 'normregion'))
                    tmp[k].append(tmp_el)

        if cov_strand:
            #1. or 2. signal -> pos/neg strand -> matrix with rep x bins
            overall_coverage_strand = [[np.matrix(tmp2[0][0]), np.matrix(tmp2[0][1])], [np.matrix(tmp2[1][0]), np.matrix(tmp2[0][1])]]
            #list of matrices: #replicates (row) x #bins (columns)
            overall_coverage = [np.matrix(tmp[0]), np.matrix(tmp[1])]
         
            return overall_coverage, overall_coverage_strand
        else:
            return [np.matrix(tmp[0]), np.matrix(tmp[1])]
    
    def count_positive_signal(self):
        return np.sum([self.covs[i].coverage for i in range(self.dim_1 + self.dim_2)])
    
    def __init__(self, name, dims, regions, genome_path, binsize, stepsize, chrom_sizes, norm_regionset, \
                 verbose, debug, no_gc_content, rmdup, path_bamfiles, exts, path_inputs, exts_inputs, \
                 factors_inputs, chrom_sizes_dict, scaling_factors_ip, save_wig, strand_cov, housekeeping_genes,\
                 tracker, end, counter, gc_content_cov=None, avg_gc_content=None, gc_hist=None, output_bw=True,\
                 folder_report=None, report=None):
        """Compute CoverageSets, GC-content and normalize input-DNA and IP-channel"""
        self.genomicRegions = regions
        self.binsize = binsize
        self.stepsize = stepsize
        self.name = name
        self.chrom_sizes_dict = chrom_sizes_dict
        self.dim_1, self.dim_2 = dims
        self.exts = exts
        self.exts_inputs = exts_inputs
        self.gc_content_cov = gc_content_cov
        self.avg_gc_content = avg_gc_content
        self.gc_hist = gc_hist
        self.scaling_factors_ip = scaling_factors_ip
        self.factors_inputs = factors_inputs
        self.end = end
        self.counter = counter
        self.no_data = False
        self.FOLDER_REPORT = folder_report
        
        global DEBUG, VERBOSE
        DEBUG = debug
        VERBOSE = verbose
        
        #make data nice
        self._help_init(path_bamfiles, exts, rmdup, binsize, stepsize, path_inputs, exts_inputs, sum(dims), regions, norm_regionset, strand_cov = strand_cov)
        if self.count_positive_signal() < 1:
            self.no_data = True
            return None
        self._compute_gc_content(no_gc_content, path_inputs, stepsize, binsize, genome_path, name, chrom_sizes, chrom_sizes_dict)
        self._normalization_by_input(path_bamfiles, path_inputs, name, factors_inputs)
        
        self.overall_coverage, self.overall_coverage_strand = self._help_init_overall_coverage(cov_strand=True)
        
        self._normalization_by_signal(name, scaling_factors_ip, path_bamfiles, housekeeping_genes, tracker, norm_regionset, report)
        
        self.scores = np.zeros(len(self.overall_coverage[0]))
        self.indices_of_interest = []
    
    def _normalization_by_input(self, path_bamfiles, path_inputs, name, factors_inputs):
        """Normalize input-DNA. Use predefined factors or follow Diaz et al, 2012"""
        
        if VERBOSE:
            print("Normalize input-DNA", file=sys.stderr)
        
        if factors_inputs:
            if VERBOSE:
                print("Use with predefined factors", file=sys.stderr)
            for i in range(len(path_bamfiles)):
                self.inputs[i].scale(factors_inputs[i])
                self.covs[i].subtract(self.inputs[i])
        elif path_inputs:
            factors_inputs = []
            print("Compute factors", file=sys.stderr)
            for i in range(len(path_bamfiles)):
                rep = i if i < self.dim_1 else i-self.dim_1
                sig = 0 if i < self.dim_1 else 1
                j = 0 if i < self.dim_1 else 1
                _, n = get_normalization_factor(path_bamfiles[i], path_inputs[i], step_width=1000, zero_counts=0, \
                                                filename=name + '-norm' + str(i), debug=DEBUG, chrom_sizes_dict=self.chrom_sizes_dict, two_sample=False, stop=True)
                if n is not None:
                    print("Normalize input of Signal %s, Rep %s with factor %s"\
                           %(sig, rep, round(n, ROUND_PRECISION)) , file=sys.stderr)
                    self.inputs[i].scale(n)
                    self.covs[i].subtract(self.inputs[i])
                    factors_inputs.append(n)
        
        self.factors_inputs = factors_inputs
        
                    
    def _trim4TMM(self, m_values, a_values, m_threshold=80, a_threshold=95):
        """q=20 or q=5"""
        assert len(m_values) == len(a_values)
        
        mask = np.asarray([not x for x in np.isinf(m_values) + np.isinf(a_values)])
        
        m_values = m_values[mask]
        a_values = a_values[mask]
        
        perc_m_l = np.percentile(m_values, 100-m_threshold)
        perc_m_h = np.percentile(m_values, m_threshold)
        perc_a_l = np.percentile(a_values, 100-a_threshold)
        perc_a_h = np.percentile(a_values, a_threshold)
        
        try:
            res = filter(lambda x: not(x[0]>perc_m_h or x[0]<perc_m_l),\
                     filter(lambda x: not(x[1]>perc_a_h or x[1]<perc_a_l), zip(list(m_values.squeeze()),list(a_values.squeeze()))))
        except:
            print('something wrong %s %s' %(len(m_values), len(a_values)), file=sys.stderr)
            return np.asarray(m_values), np.asarray(a_values)
        
        if res:
            return np.asarray(map(lambda x: x[0], res)), np.asarray(map(lambda x: x[1], res))
        else:
            print('TMM normalization: nothing trimmed...', file=sys.stderr)
            return np.asarray(m_values), np.asarray(a_values)
    
    def _norm_TMM(self, overall_coverage):
        """Normalize with TMM approach, based on PePr"""
        scaling_factors_ip = []
        for j, cond_max in enumerate([self.dim_1, self.dim_2]):
            for i in range(cond_max): #normalize all replicates
                ref = np.asarray(np.sum(overall_coverage[0], axis=0) + np.sum(overall_coverage[1], axis=0), dtype='float')/ (self.dim_1 + self.dim_2)
                
                mask_ref = ref > 0
                ref = ref[mask_ref]
                data_rep = np.asarray(overall_coverage[j][i,:])[mask_ref]
                tmp = zip(data_rep, ref, data_rep + ref)
                tmp.sort(key = lambda x: x[2], reverse=True)
                tmp = tmp[:min(len(tmp), 10000)]
                
                data_rep = np.asarray(map(lambda x: x[0], tmp))
                ref = np.asarray(map(lambda x: x[1], tmp))
                data_rep = data_rep[data_rep > 0]
                ref = ref[data_rep > 0]
                
                m_values = np.log(ref / data_rep)
                a_values = 0.5 * np.log(data_rep * ref)
                try:
                    m_values, a_values = self._trim4TMM(m_values, a_values)
                    f = 2 ** (np.sum(m_values * a_values) / np.sum(a_values))
                    scaling_factors_ip.append(f)
                except:
                    print('TMM not sucessfully', file=sys.stderr)
                    scaling_factors_ip.append(1)
                
        return scaling_factors_ip
    
    def _normalization_by_signal(self, name, scaling_factors_ip, bamfiles, housekeeping_genes, tracker, norm_regionset, report):
        """Normalize signal"""
        
        if VERBOSE:
            print('Normalize ChIP-seq profiles', file=sys.stderr)
        
        if not scaling_factors_ip and housekeeping_genes:
            print('Use housekeeping gene approach', file=sys.stderr)
            scaling_factors_ip, _ = norm_gene_level(bamfiles, housekeeping_genes, name, verbose=True, folder = self.FOLDER_REPORT, report=report)
        elif not scaling_factors_ip:
            if norm_regionset:
                print('Use TMM approach based on peaks', file=sys.stderr)
                norm_regionset_coverage = self._help_init_overall_coverage(cov_strand=False) #TMM approach based on peaks
                scaling_factors_ip = self._norm_TMM(norm_regionset_coverage)
            else:
                print('Use global TMM approach ', file=sys.stderr)
                scaling_factors_ip = self._norm_TMM(self.overall_coverage) #TMM approach
        
        for i in range(len(scaling_factors_ip)):
            self.covs[i].scale(scaling_factors_ip[i]) 
        
        if scaling_factors_ip:
            for j, cond in enumerate([self.dim_1, self.dim_2]):
                for i in range(cond): #normalize all replicates
                    k = i if j == 0 else i+self.dim_1
                    self.overall_coverage[j][i,:] *= scaling_factors_ip[k]
                    if DEBUG:
                        print('Use scaling factor %s' %round(scaling_factors_ip[k], ROUND_PRECISION), file=sys.stderr)
        
        self.scaling_factors_ip = scaling_factors_ip
                              
    def __len__(self):
        """Return number of observations."""
        return len(self.indices_of_interest)
    
    def get_observation(self, mask=np.array([])):
	"""Return array of indices of observations.
        
        *Keyword arguments:*
        
        - mask -- array like, ignore masked entries in outupt

        """

        """Return indices of observations. Do not consider indices contained in <mask> array"""
        mask = np.asarray(mask)
        if not mask.size:
            mask = np.array([True]*self._get_bin_number())
        return np.asarray(np.concatenate((self.overall_coverage[0][:,mask].T, self.overall_coverage[1][:,mask].T), axis=1))
    
    
    def _get_bin_number(self):
        """Return number of bins"""
        return self.overall_coverage[0].shape[1]
    

