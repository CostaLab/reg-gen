"""
helper codes
===================
some other extra codes
"""

from .GenomicRegion import GenomicRegion
from .GenomicRegionSet import GenomicRegionSet


def get_chrom_sizes_as_genomicregionset(chrom_size_path):
    regionset = GenomicRegionSet('')
    with open(chrom_size_path) as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            chrom, end = line[0], int(line[1])
            regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))

    return regionset


def pretty(d, indent=0):
    for key, value in list(d.items()):
        print(('\t' * indent + str(key)))
        if isinstance(value, dict):
            pretty(value, indent + 1)
        else:
            for t in value:
                print(('\t' * (indent + 1) + str(t)))
