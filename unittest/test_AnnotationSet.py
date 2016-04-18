from rgt.AnnotationSet import *
from rgt.Util import GenomeData

genome = "hg19"
print("Cheching " + genome)
annot = AnnotationSet(genome,filter_havana=True,protein_coding=True,known_only=True)
print("\tloading AnnotationSet... succeeds")
gd = GenomeData(organism=genome)
print("\t"+gd.get_annotation())
print("\tloading GenomeData... succeeds")

genome = "hg38"
print("Cheching " + genome)
annot = AnnotationSet(genome,filter_havana=True,protein_coding=True,known_only=True)
print("\tloading AnnotationSet... succeeds")
gd = GenomeData(organism=genome)
print("\t"+gd.get_annotation())
print("\tloading GenomeData... succeeds")

genome = "mm9"
print("Cheching " + genome)
annot = AnnotationSet(genome,filter_havana=True,protein_coding=True,known_only=True)
print("\tloading AnnotationSet... succeeds")
gd = GenomeData(organism=genome)
print("\t"+gd.get_annotation())
print("\tloading GenomeData... succeeds")

genome = "zv9"
print("Cheching " + genome)
annot = AnnotationSet(genome,filter_havana=True,protein_coding=True,known_only=True)
print("\tloading AnnotationSet... succeeds")
gd = GenomeData(organism=genome)
print("\t"+gd.get_annotation())
print("\tloading GenomeData... succeeds")


