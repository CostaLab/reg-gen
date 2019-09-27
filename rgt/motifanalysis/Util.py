
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


def parse_filter(pattern):
    # Converting given filter lists to a dictionary that can be used by the filter function
    # dictionaries might contain invalid keys which will raise in error when applying the filter function
    valid_keys = ["name", "gene_names", "family", "uniprot_ids", "data_source", "tax_group", "species", "database",
                  "name_file", "gene_names_file"]
    filter_values = {}
    print("correct version")
    if pattern:
        items = pattern.strip().split(";")
        for i, item in enumerate(items):
            if len(item.strip().split(":")) == 1:
                raise ValueError("Could not process given filter. Please use this format: "
                                 "\"species:human,mus;data_source=selex\". Separate key and possible values by \":\", "
                                 "the values by \",\" and put a \";\" in front of a new key.")
            items[i] = item.strip().split(":")
            if not items[i][0] in valid_keys:
                raise ValueError(items[i][0] + " is not a valid key for the filter function")

        names = []
        gene_names = []

        # iterate over keys passed to filter option
        for item in items:
            key = item[0]
            values = item[1].strip().split(",")
            # key is a string, values is either a string or a list of strings

            # process name_file and gene_names_file differently
            if key == "name_file":
                if not os.path.exists(values):
                    print("invalid name_file passed to filter")
                else:
                    with open(values, "r") as f:
                        # read TF names specified in file
                        content = f.readline()
                        for line in content:
                            names.append(line.strip())

            elif key == "gene_names_file":
                if not os.path.exists(values):
                    print("invalid gene_names_file passed to filter")
                else:
                    with open(values, "r") as f:
                        # read gene names specified in file
                        content = f.readline()
                        for line in content:
                            gene_names.append(line.strip())

            else:
                filter_values[key] = values

        # ensure that filter_values["name"] and filter_values["gene_names"] are correct
        # should now contain the intersection of passed (gene-) names and content of respective file (if both is passed)
        if "name" in filter_values and names:
            names = list(set(filter_values["name"]) & set(names))
            filter_values["name"] = names
        elif names:
            filter_values["name"] = names

        if "gene_names" in filter_values and gene_names:
            gene_names = list(set(filter_values["gene_names"]) & set(gene_names))
            filter_values["gene_names"] = gene_names
        elif gene_names:
            filter_values["gene_names"] = gene_names

    return filter_values

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
        self.genes = []  # GeneSet

    def __str__(self):
        return "\t".join([str(e) for e in [self.name, self.p_value, self.corr_p_value, self.a, self.b, self.c, self.d,
                                           self.percent, self.back_percent, ",".join(self.genes)]])
