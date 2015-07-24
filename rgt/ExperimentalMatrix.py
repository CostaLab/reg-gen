from __future__ import print_function
from rgt.GenomicRegionSet import *
from rgt.GeneSet import *
import numpy
import sys
from collections import *

"""
Describe an experiment.

Authors: Ivan G. Costa, Manuel Allhoff

"""

possible_types=["genes","regions","reads"]

class ExperimentalMatrix:

    def __init__(self):
        self.names = [] # the unique name of experiment (filename)
        self.types = [] # the type of data
        self.files = {} # the path of the related file with its filename as keys
        self.fields = [] # list types of informations including names, types, files and others
        self.fieldsDict = {} # its keys are just self.fields, and the values are extra informations        
        self.objectsDict = {} # key is the names; value is GenomicRegionSet or GeneSet
        
    def read(self, file_path, is_bedgraph=False, verbose=False):
        """Read Experimental matrix file, which looks like:
        name    type    file    further1
        MPP_PU1    regions    fil21.bed    addidional_info1
        CDP_PU1    regions    file2.bed    addidional_info2
        [ ... ]
        """
        f = open(file_path,'rU')
        
        #read and check header
        #header = f.readline()
        #header = header.strip("\n")
        #header = header.split("\t")
        
        #assert(header[0] == "name")
        #assert(header[1] == "type")
        #assert(header[2] == "file")
        #self.fields = header
        
        
        
        for line in f:
            # Neglect comment lines
            if line[0] == "#": continue
            
            # Read header
            elif line[:4] == "name":
                header = line.strip("\n")
                header = line.strip(" ")
                header = line.split()
                
                assert(header[0] == "name")
                assert(header[1] == "type")
                assert(header[2] == "file")
                self.fields = header
                
                #initialize further header files        
                for fi in range(3,len(header)):
                    self.fieldsDict[ header[fi] ] = OrderedDict()
                
            # Read further information    
            else:
                line = line.strip("\n")
                line = line.strip(" ")
                line = line.split()
                
                if len(line) < 3:  # Skip the row which has insufficient information
                    #print("Ignore line, as tab-separated number of fields < 3s: %s" %line, file=sys.stderr)
                    continue
                if verbose: print("Reading: ", line, file=sys.stderr)
                
                self.names.append(line[0])
                self.files[line[0]] = line[2] #dict: filename -> filepath
                self.types.append(line[1])
                
                for fi in range(3, len(self.fields)): #read further information
                    d = self.fieldsDict[ self.fields[fi] ]
                    try:
                        d[line[fi]].append(line[0])
                    except:
                        try:
                            d[line[fi]] = [line[0]]
                        except:
                            continue
        self.types = numpy.array(self.types)
        self.names = numpy.array(self.names)
        self.load_objects(is_bedgraph, verbose=verbose)
        
    def get_genesets(self):
        """Return GeneSets"""
        return [self.objectsDict[i] for i in self.names[self.types=="genes"]]

    def get_regionsets(self):
        """Return RegionSets"""
        return [self.objectsDict[i] for i in self.names[self.types=="regions"]]
    
    def get_regionsnames(self):
        return [i for i in self.names[self.types=="regions"]]
    
    def get_readsfiles(self):
        return [self.files[i] for i in self.names[self.types=="reads"]]

    def get_readsnames(self):
        return [i for i in self.names[self.types=="reads"]]

    def load_objects(self, is_bedgraph, verbose=False):
        """Load files and initialize object"""
        for i, t in enumerate(self.types):
            if verbose: print("Loading file ", self.files[self.names[i]], file = sys.stderr)
            
            if t not in ["regions", "genes"] and verbose:
                print("Cannot load objects", file=sys.stderr)
            
            if t == "regions":
                regions = GenomicRegionSet(self.names[i])
                if is_bedgraph:
                    regions.read_bedgraph(os.path.abspath(self.files[self.names[i]]))
                    
                else:
                    regions.read_bed(os.path.abspath(self.files[self.names[i]]))  # Here change the relative path into absolute path
                self.objectsDict[self.names[i]] = regions
            
            elif t == "genes":
                genes = GeneSet(self.names[i])
                genes.read(os.path.abspath(self.files[self.names[i]]))  # Here change the relative path into absolute path
                self.objectsDict[self.names[i]] = genes
            
    def get_type(self,name,field):
        """ Return the type according to the given name and field. """
        for f in self.fieldsDict.keys():
            if f == field:
                for t in self.fieldsDict[f].keys(): 
                    if name in self.fieldsDict[f][t]:
                        return t
                    
    def get_types(self,name):
        """ Fetch all extra informations as a list according to the given name """
        result = []
        for c in self.fieldsDict.keys():
            for t in self.fieldsDict[c].keys():
                if name in self.fieldsDict[c][t]:
                    result.append(t)
        return result

    def remove_name(self, name):
        i = self.name.index(name)
        del self.types[i]
        del self.names[i]
        self.files.pop(name, None)
        for f in self.fieldsDict.keys():
            for t in self.fieldsDict[f].keys(): 
                try: self.fieldsDict[f][t].remove(name)
                except: continue
        try: self.objectsDict.pop(name, None)
        except: pass
                    

    def match_ms_tags(self,field):
        """Add more entries to match the missing tags of the given field. For example, there are tags 
        for cell like 'cell_A' and 'cell_B' for reads, but no these tag for regions. 
        Then the regions are repeated for each tags from reads to match all reads."""

        # check regions or reads have empty tag
        beds = self.get_regionsnames()
        ms_bed = False
        for bed in beds:
            if self.get_type(bed,field) == "": ms_bed = True

        bams = self.get_readsnames()
        ms_bam = False
        for bam in bams:
            if self.get_type(bam,field) == "": ms_bam = True
        
        # generate new entries
        if ms_bed:
            for bed in beds:
                for t in self.fieldDict[field].keys():
                    n = bed+"_"+t
                    self.names.append(n)
                    self.types.append("regions")
                    self.files[n] = self.files[bed]
                    self.fieldsDict[field][t].append(n)        
                    g = GenomicRegionSet(n)
                    g.read_bed(self.files[bed])
                    self.objectsDict[n] = g
                self.remove_name(bed)

        elif ms_bam:
            for bam in bams:
                for t in self.fieldDict[field].keys():
                    n = bam+"_"+t
                    self.names.append(n)
                    self.types.append("reads")
                    self.files[n] = self.files[bam]
                    self.fieldsDict[field][t].append(n)        
                    #self.objectsDict
                self.remove_name(bam)