
###################################################################################################
# Libraries
###################################################################################################

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
    filter_values = {}
    if pattern:
        items = pattern.strip().split(";")

        names = []
        gene_names = []

        # iterate over keys passed to filter option
        for i in range(0, len(items)):

            cur_item = items[i].strip().split(":")  # cur_item=[key,list of values]
            key = cur_item[0].strip()

            # process name_file and gene_names_file differently
            if key == "name_file":
                file_name = cur_item[1].strip()
                if not os.path.exists(file_name):
                    print("invalid name_file passed to filter")
                else:
                    with open(file_name, "r") as f:
                        # read TF names specified in file
                        content = f.readline()
                        for line in content:
                            names.append(line.strip())

            elif key == "gene_names_file":
                file_name = cur_item[1].strip()
                if not os.path.exists(file_name):
                    print("invalid gene_names_file passed to filter")
                else:
                    with open(file_name, "r") as f:
                        # read gene names specified in file
                        content = f.readline()
                        for line in content:
                            gene_names.append(line.strip())

            else:
                filter_values[key] = cur_item[1].strip().split(",")

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
