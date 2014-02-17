from __future__ import print_function
from rgt.GenomicRegionSet import *
from rgt.GeneSet import *
import numpy


possible_types=["genes","regions","reads"]

class ExperimentalMatrix:

    def __init__(self):
        self.names=[] # the name of experiment (filename)
        self.types=[] # the type of data
        self.files={} # the path of the related file with its filename
        self.fields=[] # list types of informations including names, types, files and others
        self.fieldsDict={} # its keys are just self.fields, and the values are extra informations        
        self.objectsDict={} # key is the names; value is GenomicRegionSet or GeneSet
    
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
            print(line)
            if len(line) <= 1:
		continue
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
        self.load_objects(dir_path=os.path.dirname(f.name))

    def get_genesets(self):
	    return [self.objectsDict[i] for i in self.names[self.types=="genes"]]

    def get_regionsets(self):
        return [self.objectsDict[i] for i in self.names[self.types=="regions"]]

    def get_readsfiles(self):
        return [self.files[i] for i in self.names[self.types=="reads"]]

    def get_readsnames(self):
        return [i for i in self.names[self.types=="reads"]]

    def load_objects(self,dir_path):
        for i,t in enumerate(self.types):
            print(self.files[self.names[i]])
            if t == "regions":
                bed = GenomicRegionSet(self.names[i])
                bed.read_bed(os.path.join(dir_path,self.files[self.names[i]]))
                self.objectsDict[self.names[i]]=bed
            if t == "genes":
                genes= GeneSet(self.names[i])
                genes.read(self.files[self.names[i]])
                self.objectsDict[self.names[i]]=genes
