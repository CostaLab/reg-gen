# Python Libraries

import sys
import argparse
import os
from os.path import expanduser
home = expanduser("~")

# Local Libraries
from rgt.GeneSet import GeneSet
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.AnnotationSet import AnnotationSet
from rgt.Util import OverlapType
tag = "emtool"

##### FUNCTIONS #########################################################
class EM:
    def __init__(self, path):
        self.path = path
        if os.path.isfile(path): 
            print("\tLoad file: "+path)
            self.load()
        else:
            print("\tCreate file: "+path)
            self.create()
    def load(self):
        self.matrix = []
        with open(self.path) as f:
            for i, line in enumerate(f):
                l = line.strip().split()
                if i == 0:
                    self.header = l
                elif not line: continue
                else:
                    self.matrix.append(l)
    def create(self, col=4, row=1):
        if col == 3:
            self.header = ["name","type","file","factor"]
        elif col > 3:
            self.header = ["name","type","file","factor"] + ["header"]*(col - 4)
        else:
            print("Error: The number of columns must be at least 3.")
            sys.exit(1)

        self.matrix = []
        for i in range(row):
            self.matrix.append([""]*col)
    def save(self):
        print("\tSave the matrix in: "+self.path)
        with open(self.path, "w") as f:
            print("\t".join(self.header), file=f)
            for line in self.matrix:
                print("\t".join(line), file=f)
    def resize(self, col):
        print("\tResize the matrix to "+str(col)+" columns.")
        if len(self.header) < col: 
            self.header += ["header"] * (col-len(self.header))
        for row in self.matrix:
            if len(row) < col: 
                row += [""] * (col-len(row))
    def sortfile(self):
        print("\tSort the file paths and add the types and names.")
        self.resize(col=len(self.header))
        for i in range(len(self.matrix)):            
            for j, c in enumerate(self.matrix[i]):
                # if os.path.isfile(c):
                if ".bam" in c or ".bw" in c or ".bigWig" in c:
                    self.matrix[i][2] = c
                    self.matrix[i][1] = "reads"
                    self.matrix[i][0] = os.path.basename(c).partition(".")[0]
                elif ".bed" in c:
                    self.matrix[i][2] = c
                    self.matrix[i][1] = "regions"
                    self.matrix[i][0] = os.path.basename(c).partition(".")[0]
    def replace(self, old, new, col=None):
        print("\tReplace the string: "+old+" into: "+new)
        for i in range(len(self.matrix)):            
            for j, c in enumerate(self.matrix[i]):
                if col and col == (j+1) and old in c:
                    self.matrix[i][j] = self.matrix[i][j].replace(old, new)
                elif not col and old in c:
                    self.matrix[i][j] = self.matrix[i][j].replace(old, new)
    def relative_path(self):
        print("\tConvert the paths into relative path")
        base = os.path.dirname(self.path)
        # base = os.path.join(os.getcwd(),base)
        base = os.path.abspath(base)
        # print(base)
        for i in range(len(self.matrix)):            
            for j, c in enumerate(self.matrix[i]):
                if j == 2:
                    commonp = os.path.commonprefix([base, os.path.abspath(self.matrix[i][j])])
                    # print(commonp)
                    self.matrix[i][j] = os.path.relpath(os.path.abspath(self.matrix[i][j]), commonp)




if __name__ == "__main__":
    ##########################################################################
    ##### PARAMETERS #########################################################
    ##########################################################################
    
    parser = argparse.ArgumentParser(description='emtool is for manipulating experimental matrix\
                                                  Author: Chao-Chung Kuo\
                                                  \nVersion: 0.0.1', 
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('filename', metavar='  ', type=str, help="Define the filename of experimental matrix.")
    parser.add_argument('-ncol', metavar='  ', type=int, default=False, help="Define the number of columns.")
    parser.add_argument('-sortfile', action="store_true", default=False, help="Move the file paths to the file column and add the types.")
    parser.add_argument('-repath', action="store_true", default=False, help="Convert the absolute paths into relative paths.")
    parser.add_argument('-replace', metavar='  ', type=str, default=None, help="Replace the string. For example, 'control>treat' will replace all string 'control' into 'treat'; 'control>treat:4' will performs the same function only in column 4.")

    # parser.add_argument('-pcol', metavar='  ', type=int, default=False, help="Locate the position for changing, such as 2,3 or 2:5,4")
    # parser.add_argument('-row', metavar='  ', type=int, default=False, help="Locate the position of row.")
    
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    else:
        args = parser.parse_args()

    em = EM(path=args.filename)
    if args.ncol: em.resize(col=args.col)
    if args.sortfile: em.sortfile()
    if args.replace:
        old = args.replace.partition(">")[0]
        new = args.replace.partition(">")[1].partition(":")[0]
        col = args.replace.partition(">")[1].partition(":")[1]
        em.replace(old, new, col)
    if args.repath: em.relative_path()

    em.save()
