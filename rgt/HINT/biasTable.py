###################################################################################################
# Libraries
###################################################################################################

# Python
import os
import warnings

warnings.filterwarnings("ignore")
from itertools import product
from math import floor

# Internal
from rgt.Util import ErrorHandler
from rgt.Util import AuxiliaryFunctions
from rgt.GenomicRegionSet import GenomicRegionSet

# External
from pysam import __version__ as ps_version
from pysam import Samfile
from pysam import Fastafile

from Bio import motifs
from Bio.Seq import Seq


###################################################################################################
# Classes
###################################################################################################

class BiasTable:
    """
    Represent a bias table.

    Authors: Eduardo G. Gusmao.
    """

    def __init__(self, original_regions=None, dnase_file_name=None, genome_file_name=None,
                 k_nb=None, forward_shift=None, reverse_shift=None, estimate_bias_type=None, output_loc=None):
        """ 
        Initializes BiasTable.
        """
        self.regions = GenomicRegionSet("Bias Regions")
        if original_regions != None:
            if original_regions.split(".")[-1] == "bed":
                self.regions.read_bed(original_regions)
            if original_regions.split(".")[-1] == "fa":
                self.regions.read_sequence(original_regions)
        self.dnase_file_name = dnase_file_name
        self.genome_file_name = genome_file_name
        self.k_nb = k_nb
        self.forward_shift = forward_shift
        self.reverse_shift = reverse_shift
        self.estimate_bias_type = estimate_bias_type
        self.output_loc = output_loc

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

    def write_tables(self, file_name, table):
        f = open(file_name + "_F.txt", "w")
        for t in table[0].keys():
            f.write(t + "\t" + str(table[0][t]) + "\n")
        f.close()
        f = open(file_name + "_R.txt", "w")
        for t in table[1].keys():
            f.write(t + "\t" + str(table[1][t]) + "\n")
        f.close()

    def estimate_table(self):
        bias_table = None
        if self.estimate_bias_type == "FRE":
            bias_table = self.estimate_table_fre()
        elif self.estimate_bias_type == "PWM":
            bias_table = self.estimate_table_pwm()
        return bias_table

    def estimate_table_fre(self):
        """ 
        Estimates bias based on HS regions or whole genome, DNase-seq signal and genomic sequences.

        Keyword arguments:
        regions -- background regions.
        dnase_file_name -- DNase-seq file name.
        genome_file_name -- Genome to fetch genomic sequences from.
        
        Return:
        bias_table_F, bias_table_R -- Bias tables.
        """

        # Parameters
        maxDuplicates = 100
        pseudocount = 1.0

        # Initializing bam and fasta
        if (self.dnase_file_name.split(".")[-1].upper() != "BAM"): return None  # TODO ERROR
        bamFile = Samfile(self.dnase_file_name, "rb")
        fastaFile = Fastafile(self.genome_file_name)

        # Initializing dictionaries
        obsDictF = dict()
        obsDictR = dict()
        expDictF = dict()
        expDictR = dict()

        ct_reads_r = 0
        ct_reads_f = 0
        ct_kmers = 0

        # Iterating on HS regions
        for region in self.regions:

            # Initialization
            prevPos = -1
            trueCounter = 0

            # Evaluating observed frequencies ####################################
            # Fetching reads
            for r in bamFile.fetch(region.chrom, region.initial, region.final):

                # Calculating positions
                if (not r.is_reverse):
                    cut_site = r.pos + self.forward_shift - 1
                    p1 = cut_site - int(floor(self.k_nb / 2))
                else:
                    cut_site = r.aend + self.reverse_shift + 1
                    p1 = cut_site - int(floor(self.k_nb / 2))
                p2 = p1 + self.k_nb

                # Verifying PCR artifacts
                if (p1 == prevPos):
                    trueCounter += 1
                else:
                    prevPos = p1
                    trueCounter = 0
                if (trueCounter > maxDuplicates): continue

                # Fetching k-mer
                try:
                    currStr = str(fastaFile.fetch(region.chrom, p1, p2)).upper()
                except Exception:
                    continue
                if (r.is_reverse): currStr = AuxiliaryFunctions.revcomp(currStr)

                # Counting k-mer in dictionary
                if (not r.is_reverse):
                    ct_reads_f += 1
                    try:
                        obsDictF[currStr] += 1
                    except Exception:
                        obsDictF[currStr] = 1
                else:
                    ct_reads_r += 1
                    try:
                        obsDictR[currStr] += 1
                    except Exception:
                        obsDictR[currStr] = 1

            # Evaluating expected frequencies ####################################
            # Fetching whole sequence
            try:
                currStr = str(fastaFile.fetch(region.chrom, region.initial, region.final)).upper()
            except Exception:
                continue
            currRevComp = AuxiliaryFunctions.revcomp(currStr)

            # Iterating on each sequence position
            for i in range(0, len(currStr) - self.k_nb):
                ct_kmers += 1
                # Counting k-mer in dictionary
                s = currStr[i:i + self.k_nb]
                try:
                    expDictF[s] += 1
                except Exception:
                    expDictF[s] = 1

                # Counting k-mer in dictionary for reverse complement
                s = currRevComp[i:i + self.k_nb]
                try:
                    expDictR[s] += 1
                except Exception:
                    expDictR[s] = 1

        # Closing files
        bamFile.close()
        fastaFile.close()

        # Creating bias dictionary
        alphabet = ["A", "C", "G", "T"]
        kmerComb = ["".join(e) for e in product(alphabet, repeat=self.k_nb)]
        bias_table_F = dict([(e, 0.0) for e in kmerComb])
        bias_table_R = dict([(e, 0.0) for e in kmerComb])
        for kmer in kmerComb:
            try:
                obsF = obsDictF[kmer] + pseudocount
            except Exception:
                obsF = pseudocount
            try:
                expF = expDictF[kmer] + pseudocount
            except Exception:
                expF = pseudocount
            if ct_reads_f == 0:
                bias_table_F[kmer] = 1
            else:
                bias_table_F[kmer] = round(float(obsF / ct_reads_f) / float(expF / ct_kmers), 6)
            try:
                obsR = obsDictR[kmer] + pseudocount
            except Exception:
                obsR = pseudocount
            try:
                expR = expDictR[kmer] + pseudocount
            except Exception:
                expR = pseudocount
            if ct_reads_r == 0:
                bias_table_R[kmer] = 1
            else:
                bias_table_R[kmer] = round(float(obsR / ct_reads_r) / float(expR / ct_kmers), 6)

        # Return
        return [bias_table_F, bias_table_R]

    def get_pwm_score(self, sequence, pwm, k_nb):
        score = 1.0
        for position in range(k_nb):
            letter = sequence[position]
            score *= pwm[letter][position]
        return score

    def estimate_table_pwm(self):
        """
        Estimates bias based on HS regions or whole genome, DNase-seq signal and genomic sequences.

        Keyword arguments:
        regions -- background regions.
        atac_file_name -- DNase-seq file name.
        genome_file_name -- Genome to fetch genomic sequences from.

        Return:
        bias_table_F, bias_table_R -- Bias tables.
        """

        # Initializing bam and fasta
        if (self.dnase_file_name.split(".")[-1].upper() != "BAM"): return None  # TODO ERROR
        bamFile = Samfile(self.dnase_file_name, "rb")
        fastaFile = Fastafile(self.genome_file_name)

        obs_f_pwm_dict = dict([("A", [0.0] * self.k_nb), ("C", [0.0] * self.k_nb),
                        ("G", [0.0] * self.k_nb), ("T", [0.0] * self.k_nb), ("N", [0.0] * self.k_nb)])
        exp_f_pwm_dict = dict([("A", [0.0] * self.k_nb), ("C", [0.0] * self.k_nb),
                        ("G", [0.0] * self.k_nb), ("T", [0.0] * self.k_nb), ("N", [0.0] * self.k_nb)])
        obs_r_pwm_dict = dict([("A", [0.0] * self.k_nb), ("C", [0.0] * self.k_nb),
                        ("G", [0.0] * self.k_nb), ("T", [0.0] * self.k_nb), ("N", [0.0] * self.k_nb)])
        exp_r_pwm_dict = dict([("A", [0.0] * self.k_nb), ("C", [0.0] * self.k_nb),
                        ("G", [0.0] * self.k_nb), ("T", [0.0] * self.k_nb), ("N", [0.0] * self.k_nb)])

        # Iterating on HS regions
        for region in self.regions:
            # Evaluating observed frequencies
            # Fetching reads
            for r in bamFile.fetch(region.chrom, region.initial, region.final):
                # Calculating positions
                # if(not r.is_reverse): p1 = r.pos - (k_nb/2) - 1 + shift
                # else: p1 = r.aend - (k_nb/2) + 1 - shift
                if (not r.is_reverse):
                    cut_site = r.pos + self.forward_shift - 1
                    p1 = cut_site - int(floor(self.k_nb / 2))
                else:
                    cut_site = r.aend + self.reverse_shift + 1
                    p1 = cut_site - int(floor(self.k_nb / 2))
                p2 = p1 + self.k_nb

                # Fetching k-mer
                try:
                    currStr = str(fastaFile.fetch(region.chrom, p1, p2)).upper()
                except Exception:
                    continue
                if (r.is_reverse): currStr = AuxiliaryFunctions.revcomp(currStr)

                # Counting k-mer in dictionary
                if (not r.is_reverse):
                    for i in range(0, len(currStr)):
                        obs_f_pwm_dict[currStr[i]][i] += 1
                else:
                    for i in range(0, len(currStr)):
                        obs_r_pwm_dict[currStr[i]][i] += 1

            # Evaluating expected frequencies
            # Fetching whole sequence
            try:
                currStr = str(fastaFile.fetch(region.chrom, region.initial, region.final)).upper()
            except Exception:
                continue

            # Iterating on each sequence position
            s = None
            for i in range(0, len(currStr) - self.k_nb):
                # Counting k-mer in dictionary
                s = currStr[i:i + self.k_nb]
                for i in range(0, len(s)):
                    exp_f_pwm_dict[s[i]][i] += 1

                # Counting k-mer in dictionary for reverse complement
                    s = AuxiliaryFunctions.revcomp(s)
                for i in range(0, len(s)):
                    exp_r_pwm_dict[s[i]][i] += 1

        # Closing files
        bamFile.close()
        fastaFile.close()

        # Output pwms
        pwm_dict_list = [obs_f_pwm_dict, obs_r_pwm_dict, exp_f_pwm_dict, exp_r_pwm_dict]
        pwm_file_list = []
        pwm_obs_f = os.path.join(self.output_loc, "Bias", "pwm",
                                       "obs_{}_{}_f.pwm".format(str(self.k_nb), str(self.forward_shift)))
        pwm_obs_r = os.path.join(self.output_loc, "Bias", "pwm",
                                       "obs_{}_{}_r.pwm".format(str(self.k_nb), str(self.forward_shift)))
        pwm_exp_f = os.path.join(self.output_loc, "Bias", "pwm",
                                       "exp_{}_{}_f.pwm".format(str(self.k_nb), str(self.forward_shift)))
        pwm_exp_r = os.path.join(self.output_loc, "Bias", "pwm",
                                       "exp_{}_{}_r.pwm".format(str(self.k_nb), str(self.forward_shift)))

        pwm_file_list.append(pwm_obs_f)
        pwm_file_list.append(pwm_obs_r)
        pwm_file_list.append(pwm_exp_f)
        pwm_file_list.append(pwm_exp_r)

        for i in range(len(pwm_dict_list)):
            with open(pwm_file_list[i], "w") as pwm_file:
                for e in ["A", "C", "G", "T"]:
                    pwm_file.write(" ".join([str(int(f)) for f in pwm_dict_list[i][e]]) + "\n")

        motif_obs_f = motifs.read(open(pwm_obs_f), "pfm")
        motif_obs_r = motifs.read(open(pwm_obs_r), "pfm")
        motif_exp_f = motifs.read(open(pwm_exp_f), "pfm")
        motif_exp_r = motifs.read(open(pwm_exp_r), "pfm")

        # Output logos
        logo_obs_f = os.path.join(self.output_loc, "Bias", "logo",
                                       "obs_{}_{}_f.pdf".format(str(self.k_nb), str(self.forward_shift)))
        logo_obs_r = os.path.join(self.output_loc, "Bias", "logo",
                                       "obs_{}_{}_r.pdf".format(str(self.k_nb), str(self.forward_shift)))
        logo_exp_f = os.path.join(self.output_loc, "Bias", "logo",
                                       "exp_{}_{}_f.pdf".format(str(self.k_nb), str(self.forward_shift)))
        logo_exp_r = os.path.join(self.output_loc, "Bias", "logo",
                                       "exp_{}_{}_r.pdf".format(str(self.k_nb), str(self.forward_shift)))
        motif_obs_f.weblogo(logo_obs_f, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.2, yaxis_tic_interval=0.1)
        motif_obs_r.weblogo(logo_obs_r, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.2, yaxis_tic_interval=0.1)
        motif_exp_f.weblogo(logo_exp_f, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.02, yaxis_tic_interval=0.01)
        motif_exp_r.weblogo(logo_exp_r, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.02, yaxis_tic_interval=0.01)



        # Creating bias dictionary
        print(motif_exp_r.pwm)
        alphabet = ["A", "C", "G", "T"]
        k_mer_comb = ["".join(e) for e in product(alphabet, repeat=self.k_nb)]
        bias_table_F = dict([(e, 0.0) for e in k_mer_comb])
        bias_table_R = dict([(e, 0.0) for e in k_mer_comb])
        for k_mer in k_mer_comb:
            obsF = self.get_pwm_score(k_mer, motif_obs_f.pwm, self.k_nb)
            expF = self.get_pwm_score(k_mer, motif_exp_f.pwm, self.k_nb)
            bias_table_F[k_mer] = round(obsF / expF, 6)
            obsR = self.get_pwm_score(k_mer, motif_obs_r.pwm, self.k_nb)
            expR = self.get_pwm_score(k_mer, motif_exp_r.pwm, self.k_nb)
            bias_table_R[k_mer] = round(obsR / expR, 6)

        # Return
        return [bias_table_F, bias_table_R]
