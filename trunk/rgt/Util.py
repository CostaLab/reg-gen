import os

class GenomePath:
    p = os.path.dirname(__file__) + '/packagePathFile.txt'
    if os.path.exists(p):
        package_path_file = open(p,'r')
        rootDir = package_path_file.readline()
        package_path_file.close()
        MM9 = rootDir + '/trunk/data/mm9/chrom.sizes'
        HG19 = rootDir + '/trunk/data/hg19/chrom.sizes'
    
class OverlapType:
    OVERLAP = 0 
    ORIGINAL = 1
    COMP_INCL = 2