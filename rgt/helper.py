from GenomicRegionSet import GenomicRegionSet
from GenomicRegion import GenomicRegion

def get_chrom_sizes_as_genomicregionset(chrom_size_path):
    regionset = GenomicRegionSet('')
    with open(chrom_size_path) as f:
        for line in f:
            line = line.strip()
            line = line.split('\t')
            chrom, end = line[0], int(line[1])
            regionset.add(GenomicRegion(chrom=chrom, initial=0, final=end))
    
    return regionset