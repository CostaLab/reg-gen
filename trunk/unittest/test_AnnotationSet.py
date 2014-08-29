from rgt.AnnotationSet import *
from Util import GenomeData

annot = AnnotationSet("hg19")

gd = GenomeData(organism="hg19")
print(gd.get_gencode_annotation())

