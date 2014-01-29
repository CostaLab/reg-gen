import os

class GenomePath:
    package_path_file = open(os.path.dirname(__file__) + '/packagePathFile.txt','r')
    
    rootDir = package_path_file.readline()
    package_path_file.close()
    MM9 = rootDir + '/trunk/data/mm9/chrom.sizes'
    HG19 = rootDir + '/trunk/data/hg19/chrom.sizes'
    
class OverlapType:
    OVERLAP = 0 
    ORIGINAL = 1
    COMP_INCL = 2