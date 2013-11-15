possible_types=["genes","regions","reads"]
from rgt.GenomicRegionSet import *
from rgt.GeneSet import *
import numpy

class ExperimentalMatrix:

    def __init__(self):
        self.fields=[]
        self.fieldsDict={}
        self.names=[]
        self.files={}
        self.types=[]
        self.objectsDict={}
    
    def read(self,file):
        f=open(file)
        header=f.readline()
        header=header.strip("\n")
        header=header.split("\t")
        assert(header[0]=="name")
        assert(header[1]=="type")
        assert(header[2]=="file")
        self.fields=header
        for fi in range(3,len(self.fields)):
            self.fieldsDict[header[fi]]={}
        for line in f:
            line=line.strip("\n")
            line=line.split("\t")
            self.names.append(line[0])
            self.files[line[0]]=line[2]
            self.types.append(line[1])
            for fi in range(3,len(self.fields)):
                dict=self.fieldsDict[header[fi]]
                try:
                    dict[line[fi]].append(line[0])
                except:
                    dict[line[fi]]=[line[0]]
	self.types=numpy.array(self.types)
	self.names=numpy.array(self.names)
        self.load_objects()

    def get_genesets(self):
	return [self.objectsDict[i] for i in self.names[self.types=="genes"]]

    def get_regionsets(self):
        return [self.objectsDict[i] for i in self.names[self.types=="regions"]]


    def load_objects(self):
        for i,t in enumerate(self.types):
            if t == "regions":
                bed = GenomicRegionSet(self.names[i])
                bed.read_bed(self.files[self.names[i]])
                self.objectsDict[self.names[i]]=bed
            if t == "genes":
                genes= GeneSet(self.names[i])
                genes.read(self.files[self.names[i]])
                self.objectsDict[self.names[i]]=genes

            
