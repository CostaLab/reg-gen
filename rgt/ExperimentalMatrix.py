"""
ExperimentalMatrix
===================
ExperimentalMatrix describes an experiment.

"""

# Python
from __future__ import print_function
import os
import sys
from collections import *

# Internal
from rgt.GenomicRegionSet import *
from rgt.GeneSet import *
# from rgt.helper import pretty

# External

possible_types=["genes","regions","reads"]

class ExperimentalMatrix:
    """ Describes an experimental matrix.

    *Variables:*

        - names -- The unique name of experiment (filename).
        - types -- The type of data.
        - files -- The path of the related file with its filename as keys.
        - fields -- List types of informations including names, types, files and others.
        - fieldsDict -- Its keys are just self.fields, and the values are extra informations.
        - objectsDict -- Key is the names; value is GenomicRegionSet or GeneSet.
        - trash -- List of names being deleted
    """

    def __init__(self):
        self.names = []
        self.types = []
        self.files = {}
        self.fields = []
        self.fieldsDict = {}
        self.objectsDict = {}
        self.trash = []
        
    def read(self, file_path, is_bedgraph=False, verbose=False, test=False):
        """Read Experimental matrix file.

        *Keyword arguments:*

            - file_path -- Experimental matrix file path + name.
            - is_bedgraph -- Whether regions are in bedgraph format (default = False).
            - verbose -- Verbose output (default = False).
            - test -- Fetch only 10 regions form each BED files for test.

        *Example of experimental matrix file:*

            ======== ======== ========= ==================
              name      type     file     further1
            ======== ======== ========= ==================
            MPP_PU1  regions  file1.bed  addidional_info1
            CDP_PU1  regions  file2.bed  addidional_info2
            [ ... ]
            ======== ======== ========= ==================
        """
        f = open(file_path,'rU')

        base_dir = ""
        for line in f:
            # Neglect comment lines
            line = line.strip()
            if not line: continue
            elif line[0] == "#": continue
            
            # Read header
            elif line[:4] == "name":
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
                line = line.split()
                if line[0].startswith("BASE_DIR"):
                    base_dir = line[1]
                if len(line) < 3:  # Skip the row which has insufficient information
                    #print("Ignore line, as tab-separated number of fields < 3s: %s" %line, file=sys.stderr)
                    continue
                if verbose: print("Reading: ", line, file=sys.stderr)
                
                self.names.append(line[0])
                self.files[line[0]] = os.path.join(base_dir,line[2]) #dict: filename -> filepath
                self.types.append(line[1])
                
                curr_id = None
                for fi in range(3, len(self.fields)): #read further information
                    d = self.fieldsDict[ self.fields[fi] ]
                    # print(line[fi])
                    if "," in line[fi] and "(" not in line[fi]:
                        for t in line[fi].split(","):
                            try: d[t].append(line[0]+t)
                            except: d[t] = [line[0]+t]
                            self.names.append(line[0]+t)
                            self.files[line[0]+t] = line[2]
                            self.types.append(line[1])
                            for f in range(3, len(self.fields)):
                                if f != fi:
                                    try: self.fieldsDict[ self.fields[f] ][line[f]].append(line[0]+t)
                                    except: self.fieldsDict[ self.fields[f] ][line[f]] = [line[0]+t]
                            if line[0] not in self.trash: 
                                self.trash.append(line[0])
                            
                    else:
                        if curr_id:
                            try: d[line[fi]] += curr_id
                            except: d[line[fi]] = curr_id
                        else:
                            try: d[line[fi]].append(line[0])
                            except: d[line[fi]] = [line[0]]
                            
        # self.types = numpy.array(self.types)
        # self.names = numpy.array(self.names)
        self.remove_name()
        self.load_bed_url(".")
        self.load_objects(is_bedgraph, verbose=verbose, test=test)

        
    def get_genesets(self):
        """Returns the GeneSets."""
        return [self.objectsDict[n] for i,n in enumerate(self.names) if self.types[i] =="genes"]

    def get_regionsets(self):
        """Returns the RegionSets."""
        return [self.objectsDict[n] for i,n in enumerate(self.names) if self.types[i] =="regions"]
    
    def get_regionsnames(self):
        """Returns the region names."""
        return [n for i,n in enumerate(self.names) if self.types[i] =="regions"]
    
    def get_readsfiles(self):
        """Returns the 'read' type files."""
        return [self.files[n] for i,n in enumerate(self.names) if self.types[i]=="reads"]

    def get_readsnames(self):
        """Returns the 'read' type names."""
        return [n for i,n in enumerate(self.names) if self.types[i]=="reads"]

    def load_objects(self, is_bedgraph, verbose=False, test=False):
        """Load files and initialize object.

        *Keyword arguments:*

            - is_bedgraph -- Whether regions are in bedgraph format (default = False).
            - verbose -- Verbose output (default = False).
            - test -- Fetch only 10 regions form each BED files for test.
        """
        for i, t in enumerate(self.types):
            if verbose: print("Loading file ", self.files[self.names[i]], file = sys.stderr)
            
            if t not in ["regions", "genes"] and verbose:
                print("Cannot load objects", file=sys.stderr)
            
            if t == "regions":
                regions = GenomicRegionSet(self.names[i])
                if is_bedgraph:
                    regions.read_bedgraph(os.path.abspath(self.files[self.names[i]]))
                else:
                    regions.read_bed(os.path.abspath(self.files[self.names[i]]))
                    regions.sort()
                    if test: regions.sequences = regions.sequences[0:11]
                self.objectsDict[self.names[i]] = regions
            
            elif t == "genes":
                genes = GeneSet(self.names[i])
                genes.read(os.path.abspath(self.files[self.names[i]]))  # Here change the relative path into absolute path
                self.objectsDict[self.names[i]] = genes
            
    def get_type(self,name,field):
        """Return the type according to the given name and field.

        *Keyword arguments:*

            - name -- Name to return.
            - field -- Field to return.
        """
        
        if field == "regions" or field == "reads":
            field = "factor"
        
        for t in self.fieldsDict[field].keys():
            if name in self.fieldsDict[field][t]:
                return t
        
    def get_types(self,name,skip_all=False):
        """Fetch all extra informations as a list according to the given name.

        *Keyword arguments:*

            - name -- Name to return.
        """
        result = []
        for f in self.fieldsDict.keys():
            for t in self.fieldsDict[f].keys():
                if skip_all and t == "ALL": continue
                elif name in self.fieldsDict[f][t]:
                    result.append(t)
        return result

    def remove_name(self):
        """Removes experiments by name.

        *Keyword arguments:*

            - name -- Name to remove.
        """
        self.trash = list(set(self.trash))
        # for i, name in enumerate(self.names):
        for name in self.trash:
            i = self.names.index(name)
            del self.types[i]
            del self.names[i]
            self.files.pop(name, None)

            for f in self.fieldsDict.keys():
                for t in self.fieldsDict[f].keys(): 
                    # try:
                    if name in self.fieldsDict[f][t]:
                        self.fieldsDict[f][t].remove(name)
                    if not self.fieldsDict[f][t]:
                        self.fieldsDict[f].pop(t, None)
            try: self.objectsDict.pop(name, None)
            except: pass
        self.trash = []

    def match_ms_tags(self,field, test=False):
        """Add more entries to match the missing tags of the given field. For example, there are tags for cell like 'cell_A' and 'cell_B' for reads, but no these tag for regions. Then the regions are repeated for each tags from reads to match all reads.

        *Keyword arguments:*

            - field -- Field to add extra entries.
        """

        # check regions or reads have empty tag
        altypes = self.fieldsDict[field].keys()
        if "ALL" in altypes:
            altypes.remove("ALL")
            for name in self.fieldsDict[field]["ALL"]:
                i = self.names.index(name)
                for t in altypes:
                    # print("\t"+t)
                    n = name+"_"+t
                    # print("\t\t"+n)
                    self.names.append(n)
                    self.types.append(self.types[i])
                    self.files[n] = self.files[name]
                    # types = self.get_types(name,skip_all=True)
                    # print("************")
                    # print(types)

                    for f in self.fields[3:]:
                        if f == field: 
                            try: self.fieldsDict[f][t].append(n)
                            except: self.fieldsDict[f][t] = [n]
                        else:
                            try: self.fieldsDict[f][self.get_type(name=name,field=f)].append(n)
                            except: self.fieldsDict[f][self.get_type(name=name,field=f)] = [n]
                    # for f in self.fieldsDict.keys():
                    #     for ty in types:
                    #         try: self.fieldsDict[f][ty].append(n)
                    #         except: pass
                    if self.types[i] == "regions":
                        g = GenomicRegionSet(n)
                        g.read_bed(self.files[name])
                        if test: g.sequences = g.sequences[0:11]
                        self.objectsDict[n] = g
                    self.trash.append(name)

    def remove_empty_regionset(self):
        """Remove the entry with zero regions."""
        for r in self.get_regionsnames():
            if len(self.objectsDict[r]) == 0:
                self.trash.append(r)
                print("***Warning: "+r+" has zero regions and is ignored.")
        self.remove_name()


    def load_bed_url(self, temp_dir):
        """Load the BED files which contains url as file path to temporary directory."""
        import urllib
        for key, value in self.files.iteritems():
            if self.types[self.names.index(key)] == "regions" and "http" in value:
                tmpfile = os.path.join(temp_dir, value.split('/')[-1])
                dest_name = tmpfile.partition(".gz")[0]
                if not os.path.isfile(dest_name):
                    if not os.path.isfile(tmpfile):
                        urllib.urlretrieve(value, tmpfile)
                    if value.endswith(".gz"):
                        import gzip
                        with gzip.open(tmpfile, 'rb') as infile:
                            with open(dest_name, 'w') as outfile:
                                for line in infile:
                                    outfile.write(line)
                        os.remove(tmpfile)
                    else:
                        dest_name = tmpfile

                self.files[key] = dest_name

    def add_factor_col(self):
        """Add factor column by the entry name"""
        self.fields.append("factor")
        self.fieldsDict["factor"] = {}
        for n in self.names:
            try:
                self.fieldsDict["factor"][n].append(n)
            except:
                self.fieldsDict["factor"][n] = [n]