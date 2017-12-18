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


def get_read_statistics(fname, chroms):
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
    # change stats into dictionary list with only chromosome in chroms
    stats_total = [] # whole information in file including 0 reads

    isspatial = False
    isspatial_sum = 0
    isspatial_max = 0

    for chrom in chroms:
        for each_stats in stats_tmp:
            tmp = each_stats.split()  # transform string into separate list
            # tmp in format [[chrom],[gene_len],[mapped_read_num],[unmapped_read_nums]]
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


def get_file_statistics(fnames, region_giver):
    """ get read statistical data for a file list, return a dictionary for each file"""
    if isinstance(fnames[0], list):
        file_dimension = (len(fnames), len(fnames[0]))
    else:
        file_dimension = len(fnames)

    statics = []
    chroms = region_giver.chrom_sizes_dict.keys()
    for i in range(file_dimension[0]):
        statics.append([])
        for j in range(file_dimension[1]):
            stats_total, stats_data, isspatial = get_read_statistics(fnames[i][j], chroms)
            read_size = get_read_size(fnames[i][j], stats_data)
            statics[i].append({'fname':fnames[i][j], 'stats_total':stats_total,'stats_data':stats_data, 'isspatial':isspatial, 'read_size':read_size})
    # we could get different stats_data and we need to unify them...
    return {'data':statics, 'dim':file_dimension}


def is_chrom_valid(stats_total,chrom):
    """judge if it's valid to use this chromsome to get data for one file;
    how about many files??? If one is not valid, will it affect the other tests??
    If four are 0 reads we ignore it, else we still initialize it
    """
    for stat_data in stats_total:
        if chrom == stat_data[0] and stat_data[2] > stat_data[1] * 0.00001:
            return True
    return False


def is_stats_valid(statics, chrom):
    valid = False
    if statics:
        for i in range(statics['dim'][0]):
            for j in range(statics['dim'][1]):
                valid |= is_chrom_valid(statics['data'][i][j]['stats_total'], chrom)
    return valid


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
     One problem here is that we don't need to dip into many chroms
     Since there will be sth repetive; So we just read one data
    """
    cov_f, cov_r = {},{} # store forward and reverse strand data

    f = pysam.Samfile(bam_filename, "rb")
    i =0
    for read in f.fetch(stats_data[-1][0]):
        if i >= 200000:
            break
        if not read.is_unmapped:
            if read.is_reverse: # not choose right-most pos for reverse strand
                cov_r[read.pos] = 1
            else:
                cov_f[read.pos] = 1
        i += 1
    if not cov_f and not cov_r:
        return None, None
    else:
        return cov_f, cov_r


def h(forward_keys, reverse_keys, pos):
    """:return h value for pos in fforward keys an reverese_keys"""
    if pos in forward_keys and pos in reverse_keys:
        return 2
    else: # no zeros, since h value we uses are based on one value
        return 1


def ccf(forward_keys, reverse_keys, k):
    """Return value of cross-correlation function, but here we need a better method to estimate fragment sizes
    w.r.t thesis from M, it uses h() function for one position in Union of forward and reverse keys
    """
    # sums = 0
    tmp_reverse_keys = set([x-k for x in reverse_keys])
    keys = forward_keys & tmp_reverse_keys  # union of positions
    return len(keys), k


def get_extension_size(fname, stats_data, read_size, start=0, end=600, stepsize=3):
    """return extension size for each bam files.. But it should be independent

    Before its name is get_fragment_size but later changed into get_extension_szie
    """
    cov_f, cov_r = init(fname, stats_data)
    if cov_f and cov_r:
        # start = max(read_size, start - read_size)
        start = max(0, start - read_size)
        forward_keys = set(cov_f.keys()) # left-most position for forward strand
        reverse_keys = set(cov_r.keys())  # not right-most position for reverse strand
        r = [ccf(forward_keys, reverse_keys, k) for k in range(start, end, stepsize)]
        cov_f.clear()
        cov_r.clear()
        r = sorted(r, reverse=True)
        return r[0][1], r
    else:
        cov_r.clear()
        cov_f.clear()
        return None, None


def compute_extension_sizes(signal_statics, report=False):
    """Compute Extension sizes for bamfiles and input files
    Argument: signal_files are in a list format are: [[sample1_file1, sample1_file2], [sample2_file1, sample2_file2]]
          inputs_files are in the same format
    Return:
        signal_extension_sizes, read_sizes, inputs_extension_sizes, inputs_read_sizes
    Bad Aspect: we include read_size and extension_sizes together,
     so only if we use compute_extension_size we are going to use it;;
    signal_statics includes stats, fname data, and dimesions, the same as inputs_statics
    """
    start = 0
    end = 600
    ext_stepsize = 3
    ext_data_list = []

    if signal_statics:
        print("Computing read extension sizes", file=sys.stderr)
        file_dimension = signal_statics['dim']
        signal_extension_sizes = []

        for i in range(file_dimension[0]):
            signal_extension_sizes.append([])
            for j in range(file_dimension[1]):
                ext_size, ext_data = get_extension_size(signal_statics['data'][i][j]['fname'], signal_statics['data'][i][j]['stats_data'], \
                                                 read_size=signal_statics['data'][i][j]['read_size'], start=start, end=end, stepsize=ext_stepsize)
                signal_extension_sizes[i].append(ext_size)
                ext_data_list.append(ext_data)
                signal_statics['data'][i][j]['extension_size'] = ext_size
        print('end of compute extension size for signal files ', file=sys.stderr)

        if report and ext_data_list:
            print(ext_data_list, signal_extension_sizes)
        return signal_extension_sizes


def adjust_extension_sizes(signal_statics, inputs_statics):
    """we test if the inputs data are good for process, if not we need to adjust values for it"""
    file_dimension = signal_statics['dim']
    exts_inputs = []
    for i in range(file_dimension[0]):
        exts_inputs.append([])
        for j in range(file_dimension[1]):
            if inputs_statics['data'][i][j]['extension_size'] + inputs_statics['data'][i][j]['read_size'] < (
                        signal_statics['data'][i][j]['extension_size'] + signal_statics['data'][i][j]['read_size']) * 0.3:
                # we need to make sure how many files for one sample, it's on dims and divide it
                inputs_statics['data'][i][j]['extension_size'] = np.mean([ signal_statics['data'][i][k]['extension_size'] +
                                                                           signal_statics['data'][i][k]['read_size'] for k in range(file_dimension[1]) ]) \
                                                                 - inputs_statics['data'][i][j]['read_size']
            exts_inputs[i].append(inputs_statics['data'][i][j]['extension_size'])

    return exts_inputs


def update_statics_extension_sizes(statics, exts_list):
    """update extension sizes for statics data from given 1-D list"""
    if statics and exts_list:
        for i in range(statics['dim'][0]):
            for j in range(statics['dim'][1]):
                statics['data'][i][j]['extension_size'] = exts_list[i][j]


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


