###################################################################################################
# Libraries
###################################################################################################
from rgt.GenomicRegionSet import GenomicRegionSet

class BiasTable:
    """
    Represent a bias table.

    Authors: Eduardo G. Gusmao.
    """

    def __init__(self, reads_file=None, regions_file=None, genome_file=None, k_nb=None,
                 forward_shift=None, reverse_shift=None, bias_type=None, output_location=None):
        """ 
        Initializes BiasTable.
        """
        self.regions = GenomicRegionSet("Bias Regions")
        self.reads_file = reads_file
        if regions_file is not None:
            self.regions.read(regions_file)
        self.genome_file = genome_file
        self.k_nb = k_nb
        self.forward_shift = forward_shift
        self.reverse_shift = reverse_shift
        self.bias_type = bias_type
        self.output_location = output_location

    def load_table(self, table_file_name_F, table_file_name_R):
        """ 
        Creates a bias table from a tab separated file with a k-mer and bias estimate in each line.

        Keyword arguments:
        table_file_name -- Table file name.
        
        Return:
        bias_table_F, bias_table_R -- Bias tables.
        """
        bias_table_F = dict()
        table_file_F = open(table_file_name_F, "r")
        for line in table_file_F:
            ll = line.strip().split("\t")
            bias_table_F[ll[0]] = float(ll[1])
        table_file_F.close()
        bias_table_R = dict()
        table_file_R = open(table_file_name_R, "r")
        for line in table_file_R:
            ll = line.strip().split("\t")
            bias_table_R[ll[0]] = float(ll[1])
        table_file_R.close()
        return [bias_table_F, bias_table_R]
