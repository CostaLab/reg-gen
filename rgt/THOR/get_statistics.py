"""This file would be used to get statistics about bam files; This information could be used in
 <1> compute extension_size
 <2> simply judge which chromosomes we have and decide if we prpceed the process:
    <2.1 > get training data according to different distributions
    <2.2 > then peak_calling in different areas
 """
from __future__ import print_function
import pysam

def get_read_size(filename):
    """to get  read size for each files, and maybe extension size, if we don't consider the """
    f = pysam.Samfile(filename, "rb")

    s = []
    i = 0
    for read in f.fetch():
        i += 1
        if i == 1000:
            break
        if not read.is_unmapped:
            s.append(read.rlen)
    return sum(s) / len(s)

def get_extension_size(filename):
    """return extension size for each bam files.. But it should be independent from get_read_size"""
    pass


def get_read_statistics(filename):
    """ return how many chromosomes are in each files and how many reads are correspondes to each file
    it should return a dictionary, mapped information
    """
    pass



