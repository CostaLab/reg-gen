from __future__ import print_function
from rgt.GenomicRegionSet import *
from rgt.GeneSet import *
import numpy
import sys

"""
Describe an experiment.

Authors: Ivan G. Costa, Manuel Allhoff

"""

possible_types=["genes","regions","reads"]

class ExperimentalMatrix:

    def __init__(self):
        self.names = [] # the name of experiment (filename)
        self.types = [] # the type of data
        self.files = {} # the path of the related file with its filename as keys
        self.fields = [] # list types of informations including names, types, files and others
        self.fieldsDict = {} # its keys are just self.fields, and the values are extra informations        
        self.objectsDict = {} # key is the names; value is GenomicRegionSet or GeneSet
    
    def read(self, file_path):
        """Read Experimental matrix file, which looks like:
        name    type    file    further1
        MPP_PU1    regions    fil21.bed    addidional_info1
        CDP_PU1    regions    file2.bed    addidional_info2
        [ ... ]
        """
        f = open(file_path)
        
        #read and check header
        header = f.readline()
        header = header.strip("\n")
        header = header.split("\t")
        
        assert(header[0] == "name")
        assert(header[1] == "type")
        assert(header[2] == "file")
        self.fields = header
        
        #initialize further header files        
        for fi in range(3,len(self.fields)):
            self.fieldsDict[ header[fi] ] = {}
        
        for line in f:
            line = line.strip("\n")
            line = line.split("\t")
            self.names.append(line[0])
            self.files[line[0]] = line[2] #dict: filename -> filepath
            self.types.append(line[1])
            
            for fi in range(3, len(self.fields)): #read further information
                d = self.fieldsDict[ header[fi] ]
                try:
                    d[line[fi]].append(line[0])
                except:
                    d[line[fi]] = [line[0]]
        
        self.types = numpy.array(self.types)
        self.names = numpy.array(self.names)
        self.load_objects()

    def get_genesets(self):
        """Return GeneSets"""
        return [self.objectsDict[i] for i in self.names[self.types=="genes"]]

    def get_regionsets(self):
        """Return RegionSets"""
        return [self.objectsDict[i] for i in self.names[self.types=="regions"]]

    def get_readsfiles(self):
        return [self.files[i] for i in self.names[self.types=="reads"]]

    def get_readsnames(self):
        return [i for i in self.names[self.types=="reads"]]

    def load_objects(self):
        """Load files and initialize object"""
        for i, t in enumerate(self.types):
            print("Loading file ", self.files[self.names[i]], file = sys.stderr)
            
            if t == "regions":
                regions = GenomicRegionSet(self.names[i])
                regions.read_bed(self.files[self.names[i]])
                self.objectsDict[self.names[i]] = regions
            
            if t == "genes":
                genes = GeneSet(self.names[i])
                genes.read(self.files[self.names[i]])
                self.objectsDict[self.names[i]] = genes
                

