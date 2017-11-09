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
    stats_string = pysam.idxstats(fname)  # return in a string form and each info (chrom_name, total_length, mapped, unmapped) separated by '\n'
    stats_tmp = stats_string.split('\n')  # separate each info
    # using stats, we initialize the distribution of sample data, all we use is read_mapped data
    chroms = RegionGiver(chrom_fname).get_chrom_dict().keys()
    # change stats into dictionary list with only chromosome in chroms
    stats_total = [] # whole information in file including 0 reads

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
         idxs = np.argsort(new_stats_nums)[-3:]
         stats_data = [stats_total[idx] for idx in idxs]
    else:
         stats_data = stats_total
    return stats_total, stats_data, isspatial


def get_file_statistics(fnames, chrom_fname):
    """ get read statistical data for a file list, return a dictionary for each file"""


    if isinstance(fnames[0], list):
        file_dimension = (len(fnames), len(fnames[0]))
    else:
        file_dimension = len(fnames)

    statics = []
    for i in range(file_dimension[0]):
        statics.append([])
        for j in range(file_dimension[1]):
            stats_total, stats_data, isspatial = get_read_statistics(fnames[i][j], chrom_fname)
            statics[i].append({'fname':fnames[i][j], 'stats_total':stats_total,'stats_data':stats_data, 'isspatial':isspatial})
    # we could get different stats_data and we need to unify them...
    return {'data':statics, 'dim':file_dimension}


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
    # print(num)
    sample_dis = [int(x * num) for x in read_dis]
    return sample_dis


def get_read_size(fname, stats_data):
    """to get  read size for each files by using statistical data for each file """
    # sample len w.r.t distribution of reads in file

    # print('get read_size')
    sample_dis = get_sample_dis(stats_data, 10000)
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


def ccf(cov_f, cov_r, k):
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
        r = [ccf(cov_f, cov_r, k) for k in range(start, end, stepsize)]
        cov_f.clear()
        cov_r.clear()
        return read_size, max(r)[1], sorted(r, reverse=True)
    else:
        cov_r.clear()
        cov_f.clear()
        return None, None, None


def compute_extension_sizes(signal_statics, inputs_statics, report):
    """Compute Extension sizes for bamfiles and input files
    Argument: signal_files are in a list format are: [[sample1_file1, sample1_file2], [sample2_file1, sample2_file2]]
          inputs_files are in the same format
    Return:
        signal_extension_sizes, read_sizes, inputs_extension_sizes, inputs_read_sizes

    signal_statics includes stats, fname data, and dimesions, the same as inputs_statics
    """
    start = 0
    end = 600
    ext_stepsize = 3
    ext_data_list = []
    file_dimension = 0
    signal_extension_sizes, read_sizes = None, None
    inputs_extension_sizes, inputs_read_sizes = None, None

    # compute extension size for all signal files
    if signal_statics:
        print("Computing read extension sizes for ChIP-seq profiles", file=sys.stderr)
        file_dimension = signal_statics['dim']
        signal_extension_sizes = np.ones(file_dimension, int) * -1
        read_sizes = np.ones(file_dimension, int) * -1

        for i in range(file_dimension[0]):
            for j in range(file_dimension[1]):

                read_size, e, ext_data = get_extension_size(signal_statics[i][j]['data']['fname'], signal_statics[i][j]['data']['stats_data'],  start=start, end=end, stepsize=ext_stepsize)
                read_sizes[i][j] = read_size
                signal_extension_sizes[i][j] = e
                ext_data_list.append(ext_data)
                signal_statics['data']['read_size'] = read_size
                signal_statics['data']['extension_size'] = e

    if inputs_statics:

        inputs_extension_sizes = np.ones(file_dimension, int) * (-1)
        inputs_read_sizes = np.ones(file_dimension, int) * -1

        for i in range(file_dimension[0]):
            for j in range(file_dimension[1]):

                read_size, ie, scores = get_extension_size(inputs_statics[i][j]['data']['fname'], inputs_statics[i][j]['data']['fname'], start=start, end=end, stepsize=ext_stepsize)
                inputs_extension_sizes[i][j] = ie
                inputs_read_sizes[i][j] = read_size

                if inputs_extension_sizes[i][j] + inputs_read_sizes[i][j] < (signal_extension_sizes[i][j] + read_sizes[i][j]) * 0.3:
                    # we need to make sure how many files for one sample, it's on dims and divide it
                    inputs_extension_sizes[i][j] = np.mean(signal_extension_sizes[i] + read_sizes[i]) - \
                                                   inputs_read_sizes[i][j]
                    print('adjust input extension size')
                inputs_statics['data']['read_size'] = read_size
                inputs_statics['data']['extension_size'] = inputs_extension_sizes[i][j]
        """
        for i in range(file_dimension[0]):
            for j in range(file_dimension[1]): 
                if inputs_extension_sizes[i][j] + inputs_read_sizes[i][j] < (signal_extension_sizes[i][j] + read_sizes[i][j]) * 0.3:
                    # we need to make sure how many files for one sample, it's on dims and divide it
                    inputs_extension_sizes[i][j] = np.mean(signal_extension_sizes[i] + read_sizes[i]) - \
                                                   inputs_read_sizes[i][j]
        """
    if report and ext_data_list:
        print(ext_data_list, signal_files)

    # when data is in same format, we choose to use np.array to save it
    return signal_extension_sizes, read_sizes, inputs_extension_sizes, inputs_read_sizes


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
    chrom_fname = '/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes'
    # path = "/home/kefang/programs/THOR_example_data/bug/test/"

    # test compute_extension_sizes for signal files and inputs files
    # signal_files = ["rep1_A.bam","rep1_B.bam","rep2_A.bam","rep2_B.bam"]
    signal_files = [['cDC_WT_H3K27ac_1.bam','cDC_WT_H3K27ac_2.bam'] ,['CDP_WT_H3K27ac_1.bam','CDP_WT_H3K27ac_2.bam']]

    stats = get_file_statistics(signal_files,chrom_fname)
    print(stats)
    # inputs_files = ["inpt1_A.bam","inpt1_B.bam"]  # ,"inpt2_A.bam","inpt2_B.bam"
    # signal_files = ["FL5_H3K27ac.100k.bam", "FL8_H3K27ac.100k.bam", "CC4_H3K27ac.100k.bam", "CC5_H3K27ac.100k.bam"]
    # signal_files = [path + tmp for tmp in signal_files]
    # signal_files = []
    """
    inputs_files = [['Input_cDC_H3K27ac.bam','Input_cDC_H3K27ac.bam'], ['Input_CDP_H3K27ac.bam','Input_CDP_H3K27ac.bam']]
    inputs_files = []
    signal_extension_sizes = []
    inputs_extension_sizes = []
    report = True
    
    signal_extension_sizes, read_sizes, inputs_extension_sizes, inputs_read_sizes = compute_extension_sizes(signal_files, signal_extension_sizes, inputs_files, inputs_extension_sizes, report)
    print(signal_extension_sizes)
    print(inputs_extension_sizes)
    print(read_sizes)
    print(inputs_read_sizes)

    dims = (2,2)
    for i in range(dims[0]):
        for j in range(dims[1]):
            if inputs_extension_sizes[i][j] + inputs_read_sizes[i][j] < (signal_extension_sizes[i][j] + read_sizes[i][j])*0.3:
                # we need to make sure how many files for one sample, it's on dims and divide it
                inputs_extension_sizes[i][j] = np.mean(signal_extension_sizes[i] + read_sizes[i]) - inputs_read_sizes[i][j]

    print(signal_extension_sizes)
    print(inputs_extension_sizes)
    print(read_sizes)
    print(inputs_read_sizes)

    # """


