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
    stats_string = pysam.idxstats(fname)  # in a list form (chrom_name, total_length, mapped, unmapped)
    stats_tmp = stats_string.split('\n')
    # using stats, we initialize the distribution of sample data, all we use is read_mapped data
    chroms = RegionGiver(chrom_fname, None).get_chrom_dict().keys()
    # change stats into dictionary list with only chromosome in chroms
    stats_total = []
    stats_data = []  # with none zero data statistic
    for chrom in chroms:
        for each_stats in stats_tmp:
            tmp = each_stats.split()  # transform string into separate list
            tmp[1:] = map(int, tmp[1:])
            if chrom in tmp:
                stats_total.append(tmp)  # change string to int
                if tmp[2] > 0:
                    stats_data.append(tmp)
                break
    return stats_total, stats_data


def get_sample_dis(stats_data, sample_size):
    """according to stats_data, we get sample distribution"""
    read_sum = float(sum(zip(*stats_data)[2]))
    read_dis = [x / read_sum for x in zip(*stats_data)[2]]
    if sample_size is not None:
        num = min(read_sum * 0.7, sample_size)
    else:
        num = min(read_sum * 0.7)
    sample_dis = [int(x * num) for x in read_dis]
    return sample_dis


def get_read_size(fname, stats_data):
    """to get  read size for each files, and maybe extension size, if we don't consider the """
    # sample len w.r.t distribution of reads in file
    sample_dis = get_sample_dis(stats_data, 1000)
    f = pysam.Samfile(fname, "rb")
    s = []

    for idx in range(len(stats_data)):
        i =0
        for read in f.fetch(stats_data[idx][0]):
            i += 1
            if i == sample_dis[idx]:
                break
            if not read.is_unmapped:
                s.append(read.rlen)
    return sum(s) / len(s)


def init(bam_filename, stats_data):
    """:return boolean for reading cov_f and cov_r for each bam files; True, if cov_f reading right, else False
        and the read length for each bamfiles
       Now there are some differences between old and new init.. I would use the old firstly and see how it's.
    """
    cov_f, cov_r = {},{} # store forward and reverse strand data
    sample_dis = get_sample_dis(stats_data, 10000)

    f = pysam.Samfile(bam_filename, "rb")

    for chrom_idx in range(len(stats_data)):
        i = 0
        for read in f.fetch(stats_data[chrom_idx][0]):
            i += 1
            if i == sample_dis[chrom_idx]:
                break
            # to find right way to do it !!! get position, why do we use h??
            if not read.is_unmapped:
                if read.is_reverse: # choose right-most pos for reverse strand
                    pos = read.pos + read_size - 1  # not sure to use read_size or read.rlen..
                    cov_r[pos] = 1
                else:
                    cov_f[read.pos] = 1
        # print(i)
    if not cov_f and not cov_r:
        return None, None
    else:
        return cov_f, cov_r


def get_hvalue(cov, pos):
    """
    if cov_f.has_key(pos) and cov_r.has_key(pos):
        return 2
    elif not cov_f.has_key(pos) and not cov_r.has_key(pos):
        return 0
    else:
        return 1
    """
    return cov[pos] if cov.has_key(pos) else 0


def ccf(cov_f, cov_r, k, small_step=1):
    """Return value of cross-correlation function"""
    sums = np.zeros(small_step*2 +1)
    forward_keys = set(cov_f.keys())
    reverse_keys = set(map(lambda x: x - k, cov_r.keys()))
    keys = forward_keys | reverse_keys # union of positions

    for p in keys:
        for idx in range(-small_step, small_step+1, 1):
            new_k = k + idx
            sums[idx + small_step] += get_hvalue(cov_f, p) * get_hvalue(cov_r, p + new_k)
    return np.mean(sums), k


def get_extension_size(fname, stats_data, start=0, end=600, stepsize=5):
    """return extension size for each bam files.. But it should be independent"""

    cov_f, cov_r = init(fname, stats_data)
    if cov_f and cov_r:
        read_size = get_read_size(fname, stats_data)
        start = max(0, start - read_size)
        # start -= read_size
        r = [ccf(cov_f,cov_r,k) for k in range(start, end, stepsize)]
        return max(r)[1], sorted(r)

    return None, None


if __name__ == "__main__":
    # unit test of codes.
    fname = '/home/kefang/programs/THOR_example_data/bug/test/CC4_H3K27ac.100k.bam'
    chrom_fname = '/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes'
    stats_total, stats_data = get_read_statistics(fname, chrom_fname)

    read_size = get_read_size(fname, stats_data)
    print(read_size)

    extension_size = get_extension_size(fname,stats_data)
    print(extension_size)

