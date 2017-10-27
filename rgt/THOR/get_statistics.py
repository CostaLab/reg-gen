"""This file would be used to get statistics about bam files; This information could be used in
 <1> compute extension_size
 <2> simply judge which chromosomes we have and decide if we prpceed the process:
    <2.1 > get training data according to different distributions
    <2.2 > then peak_calling in different areas
 """
from __future__ import print_function
import pysam
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
    chroms = RegionGiver(chrom_fname, None).get_chrom_dict().keys()
    # change stats into dictionary list with only chromosome in chroms
    stats_total = [] # whole information in file including 0 reads
    stats_data = []  # with none zero data statistic
    for chrom in chroms:
        for each_stats in stats_tmp:
            tmp = each_stats.split()  # transform string into separate list
            tmp[1:] = map(int, tmp[1:]) # change num from string to int
            if chrom in tmp:
                stats_total.append(tmp)
                if tmp[2] > 0: # read_mapped data >0 will be recorded in stats_data
                    stats_data.append(tmp)
                break
    return stats_total, stats_data


def get_sample_dis(stats_data, sample_size=None):
    """according to stats_data, we get sample distribution; then we get sample into sample_size"""
    read_sum = float(sum(zip(*stats_data)[2]))
    read_dis = [x / read_sum for x in zip(*stats_data)[2]]
    if sample_size is not None:
        num = min(read_sum * 0.7, sample_size)
    else:
        num = min(read_sum * 0.7)
    sample_dis = [int(x * num) for x in read_dis]
    return sample_dis


def get_read_size(fname, stats_data):
    """to get  read size for each files by using statistical data for each file """
    # sample len w.r.t distribution of reads in file
    # sample_dis = get_sample_dis(stats_data, 500)
    # sample_dis = get_sample_dis(stats_data, 1000)
    # sample_dis = get_sample_dis(stats_data, 2500)
    sample_dis = get_sample_dis(stats_data, 5000)
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
    # sample_dis = get_sample_dis(stats_data, 1000)
    # sample_dis = get_sample_dis(stats_data, 500)
    # sample_dis = get_sample_dis(stats_data, 2500)
    sample_dis = get_sample_dis(stats_data, 5000)
    # sample_dis = get_sample_dis(stats_data, 10000)
    # sample_dis = get_sample_dis(stats_data, 20000)

    print(sum(sample_dis))

    f = pysam.Samfile(bam_filename, "rb")

    for chrom_idx in range(len(stats_data)):
        i = 0
        for read in f.fetch(stats_data[chrom_idx][0]):
            # print(stats_data[chrom_idx][0])
            i += 1
            if i >= sample_dis[chrom_idx]:
                break
            if not read.is_unmapped:
                if read.is_reverse: # choose right-most pos for reverse strand
                    # pos = read.pos + read.rlen - 1 # for right_most position
                    cov_r[read.pos] = 1
                else:
                    cov_f[read.pos] = 1
        # print(i)
    if not cov_f and not cov_r:
        return None, None
    else:
        return cov_f, cov_r


def get_hvalue(cov,  pos):

    return cov[pos] if cov.has_key(pos) else 0


def get_hvalue2(cov_f, pos_f, cov_r, pos_r):
    if cov_f.has_key(pos_f) and cov_r.has_key(pos_r):
        return 2
    elif not cov_f.has_key(pos_f) and not cov_r.has_key(pos_r):
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

        for p in keys:
            # sums[idx + small_step] += get_hvalue(cov_f, p) & get_hvalue(cov_r, p + new_k)
            sums[idx + small_step] += get_hvalue2(cov_f, p, cov_r, p+new_k)

    return np.mean(sums), k


def ccf2(cov_f, cov_r, k):
    """Return value of cross-correlation function"""
    sums = 0
    forward_keys = set(cov_f.keys())
    reverse_keys = set(map(lambda x: x - k, cov_r.keys()))
    keys = forward_keys & reverse_keys  # union of positions

    for p in keys:
        # sums[idx + small_step] += get_hvalue(cov_f, p) & get_hvalue(cov_r, p + new_k)
        sums += get_hvalue2(cov_f, p, cov_r, p+k)

    return np.mean(sums), k



def ccf3(cov_f, cov_r, k):
    """Return value of cross-correlation function"""
    sums = 0
    forward_keys = set(cov_f.keys())
    reverse_keys = set(map(lambda x: x - k, cov_r.keys()))
    keys = forward_keys & reverse_keys  # union of positions

    for p in keys:
        # sums[idx + small_step] += get_hvalue(cov_f, p) & get_hvalue(cov_r, p + new_k)
        sums += get_hvalue(cov_f, p) & get_hvalue(cov_r, p + k)

    return sums, k



def get_fragment_size(fname, stats_data, start=0, end=600, stepsize=5):
    """return extension size for each bam files.. But it should be independent"""
    cov_f, cov_r = init(fname, stats_data)
    if cov_f and cov_r:
        read_size = get_read_size(fname, stats_data)
        # start = max(read_size, start - read_size)
        start = max(0, start - read_size)
        # start -= read_size
        r = [ccf(cov_f,cov_r,k) for k in range(start, end, stepsize)]
        # r = [ccf2(cov_f,cov_r,k) for k in range(start, end, stepsize)]
        # r = [ccf3(cov_f, cov_r, k) for k in range(start, end, stepsize)]
        cov_f.clear()
        cov_r.clear()
        return read_size, max(r)[1], sorted(r)
    else:
        cov_r.clear()
        cov_f.clear()
        return None, None, None


if __name__ == "__main__":
    """
    wname = 'ATAC.bam'
    generate_test_files(wname)
    print('end')
    """
    
    fname = '/home/kefang/programs/THOR_example_data/bug/extension_size/test_172.bam'
    chrom_fname = '/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes'
    stats_total, stats_data = get_read_statistics(fname, chrom_fname)

    read_size = get_read_size(fname, stats_data)
    print(read_size)

    # hold on, one question, f = extension_size + read_size, what we get is actually fragment_size, not extension size
    fragment_size, fragments = get_fragment_size(fname,stats_data)
    print(fragment_size)
    # print(fragments)

    extension_size = fragment_size - read_size
    print(extension_size)
    # """
