"""A file to test different extension size methods"""

import matplotlib.pyplot as plt
import numpy as np
import pysam
import os
import time
from get_statistics import get_extension_size, get_read_statistics


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
    m_fsize, m_rsize , m_time=[],[],[]
    d_fsize , d_rsize, d_time= [],[],[]

    with open(result_file) as f:
        for line in f:
            ll = line.split()
            true_fsize.append(int(ll[0]))
            d_fsize.append(int(ll[1]))
            d_rsize.append(int(ll[2]))
            d_time.append(float(ll[3]))

            m_fsize.append(int(ll[4]))
            m_rsize.append(int(ll[5]))
            m_time.append(float(ll[6]))

    print('the most varied data from true fragment size by D-method ')
    dif = np.asarray(true_fsize) - np.asarray(d_fsize)
    zidx = np.argsort(abs(dif))[-10:]
    for item in zip(zidx, dif[zidx],np.asarray(true_fsize)[zidx],np.asarray(d_fsize)[zidx]):
        print(item)

    print('the most varied data from true fragment size by M-method ')
    dif = np.asarray(true_fsize) - np.asarray(m_fsize)
    zidx = np.argsort(abs(dif))[-10:]
    for item in zip(zidx, dif[zidx],np.asarray(true_fsize)[zidx],np.asarray(m_fsize)[zidx]):
        print(item)

    x= range(45, 45 + len(true_fsize))

    fig = plt.figure(3)
    fig1 = fig.add_subplot(311)

    plt.gca().set_color_cycle(['red', 'green', 'blue', 'yellow','red'])
    fig1.plot(x,true_fsize)
    fig1.plot(x, d_fsize)
    fig1.plot(x, m_fsize)
    fig1.legend(['true_fsize', 'd_fsize', 'm_fsize'], loc='upper left')
    ## for read_size, we could print out in another figure
    fig2 = fig.add_subplot(312)
    fig2.plot(x, d_rsize)
    fig2.plot(x, m_rsize)
    fig2.legend([ 'd_rsize','m_rsize'], loc='upper left')

    fig3 = fig.add_subplot(313)
    fig3.plot(x, d_time)
    fig3.plot(x, m_time)
    fig3.legend(['d_time', 'm_time'], loc='upper left')

    """
    fig3 = fig.add_subplot(413)
    fig3.plot(x, d_fsize)
    fig3.plot(x, d_time)
    fig3.legend(['d_fsize', 'd_fsize2'], loc='upper left')

    fig4 = fig.add_subplot(414)
    fig4.plot(x, true_fsize)
    fig4.plot(x, m_fsize)
    fig4.plot(x, m_time)
    fig4.legend(['true_fsize','m_fsize', 'm_fsize2'], loc='upper left')
    """

    plt.show()


def simulate_cross_correlation(cov_f, cov_r):
    """generate simulation data and then find the extension size by using cross correlation"""




def test():
    """generate result and plot it"""
    # open all data
    chrom_fname = '/home/kefang/programs/THOR_example_data/bug/test/hg19.chrom.sizes'

    result = open('extension_result_merged_len3.txt','a',1)
    #result.write('#true fragment size\tD_fsize\tM_fsize\tD_rsize\tM_rsize\n')

    for idx in range(45, 121,5):

        fname = 'test_merged_' + str(idx) + '.bam'
        start_time_d = time.time()
        stats_total, stats_data, isspatial = get_read_statistics(fname, chrom_fname)
        # read_size = get_read_size(fname, stats_data)
        read_size, de_size,  frag_data = get_extension_size(fname, stats_data)
        end_time_d = time.time() - start_time_d
        #start_time_m = time.time()
        #read_size2, extension_size, _ = get_extension_size(fname)
        #end_time_m = time.time() - start_time_m

        # print('%d\t%d\t%d\t%d\t%d\n'%(idx,fragment_size, extension_size + read_size2, read_size, read_size2))
        #result.write('%d\t%d\t%d\t%f\t%d\t%d\t%f \n'%(idx,de_size + read_size, read_size, end_time_d, extension_size + read_size2,  read_size2, end_time_m))
        #result.write('%d\t%d\n'%(idx, extension_size + read_size2))
        result.write('%d\t%d\t%d\t%f\n' % (idx, de_size + read_size, read_size, end_time_d))

    result.close()
    # generate_test_result_plot('extension_result.txt')


if __name__ == '__main__':
    os.chdir("/home/kefang/programs/THOR_example_data/bug/extension_size/merged")
    test()

    #test_file = "/home/kefang/programs/THOR_example_data/bug/extension_size/merged/extension_result_merged_1.txt"
    # generate_test_result_plot(test_file)