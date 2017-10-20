"""This file would be used to get statistics about bam files; This information could be used in
 <1> compute extension_size
 <2> simply judge which chromosomes we have and decide if we prpceed the process:
    <2.1 > get training data according to different distributions
    <2.2 > then peak_calling in different areas
 """
from __future__ import print_function
import pysam
import matplotlib.pyplot as plt
import numpy as np
import os

from rgt.THOR.RegionGiver import RegionGiver
from get_extension_size import get_extension_size

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
    sample_dis = get_sample_dis(stats_data, 1000)
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
    sample_dis = get_sample_dis(stats_data, 10000)

    f = pysam.Samfile(bam_filename, "rb")

    for chrom_idx in range(len(stats_data)):
        i = 0
        for read in f.fetch(stats_data[chrom_idx][0]):
            i += 1
            if i == sample_dis[chrom_idx]:
                break
            if not read.is_unmapped:
                if read.is_reverse: # choose right-most pos for reverse strand
                    pos = read.pos + read.rlen - 1
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

    for idx in range(-small_step, small_step+1, 1):
        new_k = k + idx
        reverse_keys = set(map(lambda x: x - new_k, cov_r.keys()))
        keys = forward_keys | reverse_keys  # union of positions

        for p in keys:
            sums[idx + small_step] += get_hvalue(cov_f, p) * get_hvalue(cov_r, p + new_k)

    return np.mean(sums), k


def get_fragment_size(fname, stats_data, start=0, end=600, stepsize=5):
    """return extension size for each bam files.. But it should be independent"""

    cov_f, cov_r = init(fname, stats_data)
    if cov_f and cov_r:
        read_size = get_read_size(fname, stats_data)
        start = max(0, start - read_size)
        # start -= read_size
        r = [ccf(cov_f,cov_r,k) for k in range(start, end, stepsize)]
        return max(r)[1], sorted(r)

    return None, None


def generate_test_files(whole_bamfile):
    """we need to get test files according to specific fragment sizes from whole bam file
    fragment_sizes is a list of fragment size
    """
    f = pysam.Samfile(whole_bamfile, "rb")
    print('read while right')
    # we need to create from 40 to 500 bam files and open it..
    files_seq = [pysam.Samfile('test_'+str(idx)+'.bam', 'wb',template=f) for idx in range(40,501,1)]
    for read in f.fetch():
        fs = abs(read.tlen)
        if fs >39 and fs< 501:
            files_seq[fs-40].write(read)

    for idx, file in enumerate(files_seq):
        # print('index output file %d'%(idx+40))
        file.close()
        pysam.index('test_'+str(idx+40)+'.bam')
    f.close()


def generate_test_result_plot(result_file):
    """generate test result plot by given results text file with three columns"""
    true_fsize =  []
    m_fsize, m_rsize =[],[]
    d_fsize , d_rsize= [],[]

    with open(result_file) as f:
        for line in f:
            ll = line.split()
            true_fsize.append(int(ll[0]))
            d_fsize.append(int(ll[1]))
            m_fsize.append(int(ll[4]))
            d_rsize.append(int(ll[3]))
            m_rsize.append(int(ll[2]))

    print('the most varied data from true fragment size ')
    dif = np.asarray(true_fsize) - np.asarray(d_fsize)
    zidx = np.argsort(abs(dif))[-10:]
    for item in zip(zidx, dif[zidx],np.asarray(true_fsize)[zidx],np.asarray(d_fsize)[zidx]):
        print(item)

    x= range(40, 40 + len(true_fsize))

    fig = plt.figure(2)
    fig1 = fig.add_subplot(211)

    plt.gca().set_color_cycle(['black','red', 'green', 'blue', 'yellow'])
    fig1.plot(x,true_fsize)
    fig1.plot(x, d_fsize)
    fig1.plot(x, m_fsize)
    fig1.legend(['true_fsize', 'd_fsize', 'm_fsize'], loc='upper left')
    ## for read_size, we could print out in another figure
    fig2 = fig.add_subplot(212)
    plt.plot(x, d_rsize)
    plt.plot(x, m_rsize)

    fig2.legend([ 'd_rsize','m_rsize'], loc='upper left')

    plt.show()


def test():
    """generate result and plot it"""
    # open all data
    chrom_fname = '/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes'

    result = open('extension_result2.txt','a',1)
    #result.write('#true fragment size\tD_fsize\tM_fsize\tD_rsize\tM_rsize\n')

    for idx in range(40, 501):
        fname = 'test_' + str(idx) + '.bam'
        stats_total, stats_data = get_read_statistics(fname, chrom_fname)
        read_size = get_read_size(fname, stats_data)
        fragment_size, fragments = get_fragment_size(fname, stats_data)
        read_size2, extension_size, _ = get_extension_size(fname)

        # print('%d\t%d\t%d\t%d\t%d\n'%(idx,fragment_size, extension_size + read_size2, read_size, read_size2))
        result.write('%d\t%d\t%d\t%d\t%d\n'%(idx,fragment_size, extension_size + read_size2, read_size,read_size2))

    result.close()
    # generate_test_result_plot('extension_result.txt')


if __name__ == "__main__":
    # unit test of codes.
    # get all bam names in current directory
    os.chdir('/home/kefang/programs/THOR_example_data/bug/extension_size/')
    test()
    # fname = '/home/kefang/programs/THOR_example_data/bug/extension_size/test_172.bam'
    # read_size2, extension_size, _ = get_extension_size(fname)

    # generate_test_result_plot('extension_result.txt')

    """
    wname = 'ATAC.bam'
    generate_test_files(wname)
    print('end')

    
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
    """
