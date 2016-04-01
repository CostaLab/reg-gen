
###################################################################################################
# Libraries
###################################################################################################

# Python
import warnings
warnings.filterwarnings("ignore")
from itertools import product

# Internal
from .. Util import ErrorHandler
from .. Util import AuxiliaryFunctions

# External
from pysam import __version__ as ps_version
from pysam import Samfile
from pysam import Fastafile

###################################################################################################
# Classes
###################################################################################################

class BiasTable:
    """
    Represent a bias table.

    Authors: Eduardo G. Gusmao.
    """

    def __init__(self,regions=None,dnase_file_name=None,genome_file_name=None,table_file_F=None,table_file_R=None):
        """ 
        Initializes BiasTable.
        """

        self.regions = regions
        self.dnase_file_name = dnase_file_name
        self.genome_file_name = genome_file_name
        self.table_file_F = table_file_F
        self.table_file_R = table_file_R

        if(self.table_file_F): self.table = self.load_table(self.table_file_F,self.table_file_R)
        else: self.table = self.estimate_table(self.regions,self.dnase_file_name,self.genome_file_name)
        
    def load_table(self, table_file_name_F, table_file_name_R):
        """ 
        Creates a bias table from a tab separated file with a k-mer and bias estimate in each line.

        Keyword arguments:
        table_file_name -- Table file name.
        
        Return:
        bias_table_F, bias_table_R -- Bias tables.
        """
        bias_table_F = dict()
        table_file_F = open(table_file_name_F,"r")
        for line in table_file_F:
            ll = line.strip().split("\t")
            bias_table_F[ll[0]] = float(ll[1])
        table_file_F.close()
        bias_table_R = dict()
        table_file_R = open(table_file_name_R,"r")
        for line in table_file_R:
            ll = line.strip().split("\t")
            bias_table_R[ll[0]] = float(ll[1])
        table_file_R.close()
        return [bias_table_F, bias_table_R]

    def estimate_table(self, regions, dnase_file_name, genome_file_name):
        """ 
        Estimates bias based on HS regions, DNase-seq signal and genomic sequences.

        Keyword arguments:
        regions -- DNase-seq HS regions.
        dnase_file_name -- DNase-seq file name.
        genome_file_name -- Genome to fetch genomic sequences from.
        
        Return:
        bias_table_F, bias_table_R -- Bias tables.
        """

        # Parameters
        maxDuplicates = 100
        k_nb = 6
        pseudocount = 1.0

        # Initializing bam and fasta
        if(dnase_file_name.split(".")[-1].upper() != "BAM"): return None # TODO ERROR
        bamFile = Samfile(dnase_file_name, "rb")
        fastaFile = Fastafile(genome_file_name)

        # Initializing dictionaries
        obsDictF = dict(); obsDictR = dict()
        expDictF = dict(); expDictR = dict()

        # Iterating on HS regions
        for region in regions:

            # Initialization
            prevPos = -1
            trueCounter = 0

            # Evaluating observed frequencies ####################################

            # Fetching reads
            for r in bamFile.fetch(region.chrom, region.initial, region.final):

                # Calculating positions
                if(not r.is_reverse): p1 = r.pos - (k_nb/2) - 1
                else: p1 = r.aend - (k_nb/2) + 1
                p2 = p1 + k_nb

                # Verifying PCR artifacts
                if(p1 == prevPos): trueCounter += 1
                else:
                    prevPos = p1
                    trueCounter = 0
                if(trueCounter > maxDuplicates): continue

                # Fetching k-mer
                try: currStr = str(fastaFile.fetch(region.chrom, p1, p2)).upper()
                except Exception: continue
                if(r.is_reverse): currStr = AuxiliaryFunctions.revcomp(currStr)

                # Counting k-mer in dictionary
                if(not r.is_reverse):
                    try: obsDictF[currStr] += 1
                    except Exception: obsDictF[currStr] = 1
                else:
                    try: obsDictR[currStr] += 1
                    except Exception: obsDictR[currStr] = 1

            # Evaluating expected frequencies ####################################

            # Fetching whole sequence
            try: currStr = str(fastaFile.fetch(region.chrom, region.initial, region.final)).upper()
            except Exception: continue
            currRevComp = AuxiliaryFunctions.revcomp(currStr)

            # Iterating on each sequence position
            for i in range(0,len(currStr)-k_nb):

                # Counting k-mer in dictionary
                s = currStr[i:i+k_nb]
                try: expDictF[s] += 1
                except Exception: expDictF[s] = 1

                # Counting k-mer in dictionary for reverse complement
                s = currRevComp[i:i+k_nb]
                try: expDictR[s] += 1
                except Exception: expDictR[s] = 1

        # Closing files
        bamFile.close()
        fastaFile.close()

        # Creating bias dictionary
        alphabet = ["A","C","G","T"]
        kmerComb = ["".join(e) for e in product(alphabet, repeat=k_nb)]
        bias_table_F = dict([(e,0.0) for e in kmerComb]) 
        bias_table_R = dict([(e,0.0) for e in kmerComb]) 
        for kmer in kmerComb:
            try: obsF = obsDictF[kmer] + pseudocount
            except Exception: obsF = pseudocount
            try: expF = expDictF[kmer] + pseudocount
            except Exception: expF = pseudocount
            bias_table_F[kmer] = round(float(obsF)/float(expF),6)
            try: obsR = obsDictR[kmer] + pseudocount
            except Exception: obsR = pseudocount
            try: expR = expDictR[kmer] + pseudocount
            except Exception: expR = pseudocount
            bias_table_R[kmer] = round(float(obsR)/float(expR),6)

        # Return
        return [bias_table_F, bias_table_R]


