"""This file would be used to get statistics about bam files; This information could be used in
 <1> compute extension_size
 <2> simply judge which chromosomes we have and decide if we prpceed the process:
    <2.1 > get training data according to different distributions
    <2.2 > then peak_calling in different areas
 """
from __future__ import print_function
import pysam
import sys
import os
import numpy as np
import time
from scipy.signal import correlate

from rgt.THOR.RegionGiver import RegionGiver


def get_read_statistics(fname, chrom_fname):
    """ return how many chromosomes are in each files and how many reads are correspondes to each file
    it should return a dictionary, mapped information
    use pysam 0.12. there is function,
    get_index_statistics(self): self, means the BAM file,
    return statistics about mapped/unmapped reads per chromosome as they are stored in the index.
    Returns: list  mapped, unmapped and total.
    Return type: a list of records for each chromosome. Each record has the attributes contig,
    """
    # in pysam 0.9, it uses idxstats...
    # f = pysam.Samfile(fname, "rb")
    # print(fname)
    # print(os.path.isfile(fname))
    stats_string = pysam.idxstats(fname)  # return in a string form and each info (chrom_name, total_length, mapped, unmapped) separated by '\n'
    stats_tmp = stats_string.split('\n')  # separate each info
    # using stats, we initialize the distribution of sample data, all we use is read_mapped data
    chroms = RegionGiver(chrom_fname, None).get_chrom_dict().keys()
    # change stats into dictionary list with only chromosome in chroms
    stats_total = [] # whole information in file including 0 reads
    stats_data = []  # with none zero data statistic

    isspatial = False
    isspatial_sum = 0
    isspatial_max = 0

    for chrom in chroms:
        for each_stats in stats_tmp:
            tmp = each_stats.split()  # transform string into separate list
            tmp[1:] = map(int, tmp[1:]) # change num from string to int
            if chrom in tmp:
                if tmp[2] > 0: # read_mapped data >0 will be recorded in stats_data
                    stats_total.append(tmp)
                    isspatial_sum += tmp[2]
                    if tmp[2] > isspatial_max:
                        isspatial_max = tmp[2]
                break

    if isspatial_max < isspatial_sum * 0.15:
         isspatial = True
         # data is spatial, and we use biggest 5 stats_data as new stats_data
         new_stats_nums = [tmp[2] for tmp in stats_total]
         idxs = np.argsort(new_stats_nums)[-1:]
         stats_data = [stats_total[idx] for idx in idxs]
    else:
         stats_data = stats_total
    return stats_total, stats_data, isspatial


def get_sample_dis(stats_data, sample_size=None, data_spatial=False):
    """according to stats_data, we get sample distribution;
       get small total size from 0.7* data_len and sample size.
    """
    read_sum = float(sum(zip(*stats_data)[2]))
    read_dis = [x / read_sum for x in zip(*stats_data)[2]]
    if sample_size:
        num = min(read_sum * 0.7, sample_size)
    else:
        num = read_sum * 0.7
    print(num)
    sample_dis = [int(x * num) for x in read_dis]
    return sample_dis


def get_read_size(fname, stats_data):
    """to get  read size for each files by using statistical data for each file """
    # sample len w.r.t distribution of reads in file
    # sample_dis = get_sample_dis(stats_data, 500)
    # sample_dis = get_sample_dis(stats_data, 1000)
    # sample_dis = get_sample_dis(stats_data, 2500)
    print('get read_size')
    sample_dis = get_sample_dis(stats_data, 10000)
    # sample_dis = get_sample_dis(stats_data, 10000)
    # sample_dis = get_sample_dis(stats_data, 20000)
    f = pysam.Samfile(fname, "rb")
    s = []

    for idx in range(len(stats_data)):
        i = 0
        for read in f.fetch(stats_data[idx][0]):
            i += 1
            if i == sample_dis[idx]:
                break
            if not read.is_unmapped:
                s.append(read.rlen)
    return sum(s) / len(s)


def init(bam_filename, stats_data):
    """:return cov_f and cov_r in dictionary
     cov_f forward coverage, left-most position
     cov_r reverse coverage, right-most position
    """
    cov_f, cov_r = {},{} # store forward and reverse strand data
    sample_dis = get_sample_dis(stats_data, 200000)

    f = pysam.Samfile(bam_filename, "rb")

    for chrom_idx in range(len(stats_data)):
        i = 0
        for read in f.fetch(stats_data[chrom_idx][0]):
            # print(stats_data[chrom_idx][0])
            if i >= sample_dis[chrom_idx]:
                break
            if not read.is_unmapped:
                if read.is_reverse: # choose right-most pos for reverse strand
                    cov_r[read.pos] = 1
                else:
                    cov_f[read.pos] = 1
            i += 1
        # print(i)
    if not cov_f and not cov_r:
        return None, None
    else:
        return cov_f, cov_r


def get_hvalue(cov,  pos):

    return cov[pos] if cov.has_key(pos) else 0


def get_hvalue2(cov_f, cov_r, pos_f, pos_r):
    if cov_f.has_key(pos_f) and cov_r.has_key(pos_r):
        return 2
    elif not cov_f.has_key(pos_f) and not cov_r.has_key(pos_r):
        return 0
    else:
        return 1


def get_hvalue3(cov_f, cov_r, pos):
    """define function for cross correlation"""
    if cov_f.has_key(pos) and cov_r.has_key(pos):
        return 2
    elif not cov_f.has_key(pos) and not cov_r.has_key(pos):
        return 0
    else:
        return 1


def ccf(cov_f, cov_r, k, small_step=1):
    """Return value of cross-correlation function"""
    sums = np.zeros(small_step*2 +1)
    forward_keys = set(cov_f.keys())

    for idx in range(-small_step, small_step+1, 1):
        new_k = k + idx
        reverse_keys = set(map(lambda x: x - new_k, cov_r.keys()))
        keys = forward_keys & reverse_keys  # union of positions
        # I think we need at least 1000 data here to get proper result
        # print(len(keys))
        #for p in keys:
        #    sums[idx + small_step] += get_hvalue(cov_f, p) & get_hvalue(cov_r, p + new_k)
            # sums[idx + small_step] += get_hvalue2(cov_f, p, cov_r, p+new_k)
        # assert len(keys) == sums[idx + small_step], 'equal ??'
        sums[idx + small_step] = len(keys)

    return np.mean(sums), k


def ccf2(cov_f, cov_r, k):
    """Return value of cross-correlation function, return is then for fragment size"""
    sums = 0
    keys = set(cov_f.keys()) | set(cov_r.keys())  # union of positions
    # print(len(keys))

    for p in keys:
        # sums[idx + small_step] += get_hvalue(cov_f, p) & get_hvalue(cov_r, p + new_k)
        tmp = get_hvalue2(cov_f, cov_r, p, p+k)
        sums += tmp
    return sums, k


def ccf3(cov_f, cov_r, k):
    """Return value of cross-correlation function"""
    #sums = 0
    forward_keys = set(cov_f.keys())
    reverse_keys = set(map(lambda x: x - k, cov_r.keys()))
    keys = forward_keys & reverse_keys  # union of positions

    return len(keys), k


def get_extension_size(fname, stats_data, start=0, end=600, stepsize=3):
    """return extension size for each bam files.. But it should be independent

    Before its name is get_fragment_size but later changed into get_extension_szie
    """
    cov_f, cov_r = init(fname, stats_data)
    if cov_f and cov_r:
        read_size = get_read_size(fname, stats_data)
        # start = max(read_size, start - read_size)
        start = max(0, start - read_size)
        # start -= read_size
        # r = [ccf(cov_f,cov_r,k) for k in range(start, end, stepsize)]
        # r = [ccf2(cov_f,cov_r,k) for k in range(start, end, stepsize)]
        r = [ccf3(cov_f, cov_r, k) for k in range(start, end, stepsize)]
        cov_f.clear()
        cov_r.clear()
        return read_size, max(r)[1], sorted(r, reverse=True)
    else:
        cov_r.clear()
        cov_f.clear()
        return None, None, None


def compute_extension_sizes(signal_files, signal_extension_sizes, inputs_files, inputs_extension_sizes, report):
    """Compute Extension sizes for bamfiles and input files"""
    start = 0
    end = 600
    ext_stepsize = 5

    ext_data_list = []
    read_sizes, inputs_read_sizes = [],[]

    # chrom_fname = '/home/kefang/programs/THOR_example_data/bug/debugData/GRCh38_no_alt_analysis_set.sizes.genome.txt'
    chrom_fname = '/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes'
    # compute extension size for all signal files
    if not signal_extension_sizes:
        print("Computing read extension sizes for ChIP-seq profiles", file=sys.stderr)
        for bamfile in signal_files:
            start_time = time.time()
            stats_total, stats_data, isspatial = get_read_statistics(bamfile, chrom_fname)
            print(stats_data)
            read_size, e, ext_data = get_extension_size(bamfile, stats_data,  start=start, end=end, stepsize=ext_stepsize)
            read_sizes.append(read_size)
            signal_extension_sizes.append(e)
            ext_data_list.append(ext_data)
            end_time = time.time() - start_time
            print('used time: ', end_time)
            print(ext_data)

    if report and ext_data_list:
        print(ext_data_list, signal_files)
    print('for input file')
    if inputs_files and not inputs_extension_sizes:
        for bamfile in inputs_files:
            start_time = time.time()
            stats_total, stats_data, isspatial = get_read_statistics(bamfile, chrom_fname)

            print(stats_data)

            read_size, ie, scores = get_extension_size(bamfile, stats_data, start=start, end=end, stepsize=ext_stepsize)
            inputs_extension_sizes.append(ie)
            inputs_read_sizes.append(read_size)
            end_time = time.time() - start_time
            print('used time: ', end_time)
            print(scores)

    return signal_extension_sizes, read_sizes, inputs_extension_sizes,inputs_read_sizes


if __name__ == "__main__":
    """
    wname = 'ATAC.bam'
    generate_test_files(wname)
    print('end')
    
    
    fname = '/home/kefang/programs/THOR_example_data/bug/evaluate_extension_size/cDC_WT_H3K27ac_1.bam'
    chrom_fname = '/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes'
    stats_total, stats_data = get_read_statistics(fname, chrom_fname)

    read_size = get_read_size(fname, stats_data)
    print(read_size)

    cov_f, cov_r = init(fname, stats_data)


    # hold on, one question, f = extension_size + read_size, what we get is actually fragment_size, not extension size
    # extension_size, fragments = get_extension_size(fname,stats_data)
    # print(extension_size)
    # print(fragments)

    # extension_size = fragment_size - read_size
    # print(extension_size)

    """
    os.chdir("/home/kefang/programs/THOR_example_data/bug/evaluate_extension_size")
    # path = "/home/kefang/programs/THOR_example_data/bug/test/"

    # test compute_extension_sizes for signal files and inputs files
    # signal_files = ["rep1_A.bam","rep1_B.bam","rep2_A.bam","rep2_B.bam"]
    signal_files = ['cDC_WT_H3K27ac_1.bam','cDC_WT_H3K27ac_2.bam' ,'CDP_WT_H3K27ac_1.bam','CDP_WT_H3K27ac_2.bam']
    # inputs_files = ["inpt1_A.bam","inpt1_B.bam"]  # ,"inpt2_A.bam","inpt2_B.bam"
    # signal_files = ["FL5_H3K27ac.100k.bam", "FL8_H3K27ac.100k.bam", "CC4_H3K27ac.100k.bam", "CC5_H3K27ac.100k.bam"]
    # signal_files = [path + tmp for tmp in signal_files]
    # signal_files = []
    inputs_files = ['Input_cDC_H3K27ac.bam', 'Input_CDP_H3K27ac.bam']
    # inputs_files = []
    signal_extension_sizes = []
    inputs_extension_sizes = []
    report = True
    
    signal_extension_sizes, read_sizes, inputs_extension_sizes, inputs_read_sizes = compute_extension_sizes(signal_files, signal_extension_sizes, inputs_files, inputs_extension_sizes, report)
    print(signal_extension_sizes)
    print(inputs_extension_sizes)
    print(read_sizes)
    print(inputs_read_sizes)
    # """


