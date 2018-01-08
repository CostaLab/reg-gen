
###################################################################################################
# Libraries
###################################################################################################

# Python 3 compatibility
from __future__ import print_function
from __future__ import division

# Python
import os

# Internal

# External

###################################################################################################
# Functions
###################################################################################################


def is_bed(filename):
    _, ext = os.path.splitext(filename)

    if ext.lower() == ".bed":
        return True

    return False


def is_bb(filename):
    _, ext = os.path.splitext(filename)

    if ext.lower() == ".bb":
        return True

    return False


def bb_to_bed(filename):
    path, ext = os.path.splitext(filename)

    if ext.lower() == ".bb":
        # convert BB to BED
        bed_filename = os.path.join(path + ".bed")
        os.system(" ".join(["bigBedToBed", filename, bed_filename]))

        return bed_filename
    elif ext.lower() == ".bed":
        return filename
    else:
        raise ValueError("{} is neither a BED nor a BB".format(filename))


def bed_to_bb(filename, chrom_sizes_filename):
    path, ext = os.path.splitext(filename)

    if ext.lower() == ".bed":
        # convert BED to BB
        bb_filename = os.path.join(path + ".bb")
        os.system(" ".join(["bedToBigBed", filename, chrom_sizes_filename, bb_filename, "-verbose=0"]))

        return bb_filename
    elif ext.lower() == ".bb":
        return filename
    else:
        raise ValueError("{} is neither a BED nor a BB".format(filename))


def write_bed_color(region_set, filename, color):
    with open(filename, 'w') as f:
        for s in region_set:
            append = "\t".join([str(s.initial), str(s.final), color])
            print(str(s) + append, file=f)

###################################################################################################
# Classes
###################################################################################################


class Input:
    def __init__(self, gene_set, region_list):
        self.gene_set = gene_set
        self.region_list = region_list


class Result:
    def __init__(self):
        self.name = ""  # String
        self.p_value = 0.0  # Float
        self.corr_p_value = 0.0  # Float
        self.a = 0  # Integer
        self.b = 0  # Integer
        self.c = 0  # Integer
        self.d = 0  # Integer
        self.percent = 0.0  # Float
        self.back_percent = 0.0  # Float
        self.genes = None  # GeneSet

    def __str__(self):
        return "\t".join([str(e) for e in [self.name, self.p_value, self.corr_p_value, self.a, self.b, self.c, self.d,
                                           self.percent, self.back_percent, ",".join(self.genes.genes)]])
