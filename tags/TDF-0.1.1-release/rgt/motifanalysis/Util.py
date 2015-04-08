
###################################################################################################
# Libraries
###################################################################################################

# Python

# Internal

# External

###################################################################################################
# Classes
###################################################################################################

class Input:
    def __init__(self,gene_set,region_list):
        self.gene_set = gene_set
        self.region_list = region_list

class Result:
    def __init__(self):
        self.name = "" # String
        self.p_value = 0.0 # Float
        self.corr_p_value = 0.0 # Float
        self.a = 0 # Integer
        self.b = 0 # Integer
        self.c = 0 # Integer
        self.d = 0 # Integer
        self.percent = 0.0 # Float
        self.back_percent = 0.0 # Float
        self.genes = None # GeneSet
    def __str__(self):
        return "\t".join([str(e) for e in [self.name,self.p_value,self.corr_p_value,self.a,self.b,self.c,self.d,
                                               self.percent,self.back_percent,",".join(self.genes.genes)]])


