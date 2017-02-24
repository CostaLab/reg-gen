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

    def __init__(self, output_loc=None):
        """ 
        Initializes BiasTable.
        """
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

    def estimate_table(self, regions, dnase_file_name, genome_file_name, k_nb, forward_shift, reverse_shift):
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
        pseudocount = 1.0

        # Initializing bam and fasta
        if (dnase_file_name.split(".")[-1].upper() != "BAM"): return None  # TODO ERROR
        bamFile = Samfile(dnase_file_name, "rb")
        fastaFile = Fastafile(genome_file_name)

        # Initializing dictionaries
        obsDictF = dict()
        obsDictR = dict()
        expDictF = dict()
        expDictR = dict()

        ct_reads_r = 0
        ct_reads_f = 0
        ct_kmers = 0

        # Iterating on HS regions
        for region in regions:

            # Initialization
            prevPos = -1
            trueCounter = 0

            # Evaluating observed frequencies ####################################
            # Fetching reads
            for r in bamFile.fetch(region.chrom, region.initial, region.final):

                # Calculating positions
                if (not r.is_reverse):
                    cut_site = r.pos + forward_shift - 1
                    p1 = cut_site - int(floor(k_nb / 2))
                else:
                    cut_site = r.aend + reverse_shift + 1
                    p1 = cut_site - int(floor(k_nb / 2))
                p2 = p1 + k_nb

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
            for i in range(0, len(currStr) - k_nb):
                ct_kmers += 1
                # Counting k-mer in dictionary
                s = currStr[i:i + k_nb]
                try:
                    expDictF[s] += 1
                except Exception:
                    expDictF[s] = 1

                # Counting k-mer in dictionary for reverse complement
                s = currRevComp[i:i + k_nb]
                try:
                    expDictR[s] += 1
                except Exception:
                    expDictR[s] = 1

        # Closing files
        bamFile.close()
        fastaFile.close()

        # Creating bias dictionary
        alphabet = ["A", "C", "G", "T"]
        kmerComb = ["".join(e) for e in product(alphabet, repeat=k_nb)]
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

    def estimate_table_pwm(self, regions, dnase_file_name, genome_file_name, k_nb, forward_shift, reverse_shift):
        """
        Estimates bias based on HS regions, DNase-seq signal and genomic sequences.

        Keyword arguments:
        regions -- DNase-seq HS regions.
        atac_file_name -- DNase-seq file name.
        genome_file_name -- Genome to fetch genomic sequences from.

        Return:
        bias_table_F, bias_table_R -- Bias tables.
        """

        # Initializing bam and fasta
        if (dnase_file_name.split(".")[-1].upper() != "BAM"): return None  # TODO ERROR
        bamFile = Samfile(dnase_file_name, "rb")
        fastaFile = Fastafile(genome_file_name)

        obsSeqsF = []
        obsSeqsR = []
        expSeqsF = []
        expSeqsR = []

        # Iterating on HS regions
        for region in regions:
            # Evaluating observed frequencies
            # Fetching reads
            for r in bamFile.fetch(region.chrom, region.initial, region.final):
                # Calculating positions
                # if(not r.is_reverse): p1 = r.pos - (k_nb/2) - 1 + shift
                # else: p1 = r.aend - (k_nb/2) + 1 - shift
                if (not r.is_reverse):
                    cut_site = r.pos + forward_shift - 1
                    p1 = cut_site - int(floor(k_nb / 2))
                else:
                    cut_site = r.aend + reverse_shift + 1
                    p1 = cut_site - int(floor(k_nb / 2))
                p2 = p1 + k_nb

                # Fetching k-mer
                try:
                    currStr = str(fastaFile.fetch(region.chrom, p1, p2)).upper()
                except Exception:
                    continue
                if (r.is_reverse): currStr = AuxiliaryFunctions.revcomp(currStr)

                # Counting k-mer in dictionary
                if 'N' not in currStr:
                    if (not r.is_reverse):
                        obsSeqsF.append(Seq(currStr))
                    else:
                        obsSeqsR.append(Seq(currStr))

            # Evaluating expected frequencies
            # Fetching whole sequence
            try:
                currStr = str(fastaFile.fetch(region.chrom, region.initial, region.final)).upper()
            except Exception:
                continue
            currRevComp = AuxiliaryFunctions.revcomp(currStr)

            # Iterating on each sequence position
            for i in range(0, len(currStr) - k_nb):
                s = currStr[i:i + k_nb]
                if 'N' not in currStr:
                    # Counting k-mer in dictionary
                    expSeqsF.append(Seq(s))

                    # Counting k-mer in dictionary for reverse complement
                    s = currRevComp[i:i + k_nb]
                    expSeqsR.append(Seq(s))

        # Closing files
        bamFile.close()
        fastaFile.close()

        obsMotifsF = motifs.create(obsSeqsF)
        obsMotifsR = motifs.create(obsSeqsR)
        expMotifsF = motifs.create(expSeqsF)
        expMotifsR = motifs.create(expSeqsR)

        obsPwmF = obsMotifsF.pwm
        obsPwmR = obsMotifsR.pwm
        expPwmF = expMotifsF.pwm
        expPwmR = expMotifsR.pwm

        # Output logos
        logo_obs_f = os.path.join(self.output_loc, "Bias", "logo",
                                       "obs_{}_{}_f.pdf".format(str(k_nb), str(forward_shift)))
        logo_obs_r = os.path.join(self.output_loc, "Bias", "logo",
                                       "obs_{}_{}_r.pdf".format(str(k_nb), str(forward_shift)))
        logo_exp_f = os.path.join(self.output_loc, "Bias", "logo",
                                       "exp_{}_{}_f.pdf".format(str(k_nb), str(forward_shift)))
        logo_exp_r = os.path.join(self.output_loc, "Bias", "logo",
                                       "exp_{}_{}_r.pdf".format(str(k_nb), str(forward_shift)))
        obsMotifsF.weblogo(logo_obs_f, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.2, yaxis_tic_interval=0.1)
        obsMotifsR.weblogo(logo_obs_r, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.2, yaxis_tic_interval=0.1)
        expMotifsF.weblogo(logo_exp_f, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.02, yaxis_tic_interval=0.01)
        expMotifsR.weblogo(logo_exp_r, format="pdf", stack_width="large", color_scheme="color_classic",
                           yaxis_scale=0.02, yaxis_tic_interval=0.01)

        # Output pwms
        pwm_data_list = [obsPwmF, obsPwmR, expPwmF, expPwmR]
        pwm_file_list = []
        pwm_obs_f = os.path.join(self.output_loc, "Bias", "pwm",
                                       "obs_{}_{}_f.pwm".format(str(k_nb), str(forward_shift)))
        pwm_obs_r = os.path.join(self.output_loc, "Bias", "pwm",
                                       "obs_{}_{}_r.pwm".format(str(k_nb), str(forward_shift)))
        pwm_exp_f = os.path.join(self.output_loc, "Bias", "pwm",
                                       "exp_{}_{}_f.pwm".format(str(k_nb), str(forward_shift)))
        pwm_exp_r = os.path.join(self.output_loc, "Bias", "pwm",
                                       "exp_{}_{}_r.pwm".format(str(k_nb), str(forward_shift)))

        pwm_file_list.append(pwm_obs_f)
        pwm_file_list.append(pwm_obs_r)
        pwm_file_list.append(pwm_exp_f)
        pwm_file_list.append(pwm_exp_r)

        for i in range(len(pwm_data_list)):
            with open(pwm_file_list[i], "w") as f:
                f.write(str(pwm_data_list[i]))

        # Creating bias dictionary
        alphabet = ["A", "C", "G", "T"]
        k_mer_comb = ["".join(e) for e in product(alphabet, repeat=k_nb)]
        bias_table_F = dict([(e, 0.0) for e in k_mer_comb])
        bias_table_R = dict([(e, 0.0) for e in k_mer_comb])
        for k_mer in k_mer_comb:
            obsF = self.get_pwm_score(k_mer, obsPwmF, k_nb)
            expF = self.get_pwm_score(k_mer, expPwmF, k_nb)
            bias_table_F[k_mer] = round(obsF / expF, 6)
            obsR = self.get_pwm_score(k_mer, obsPwmR, k_nb)
            expR = self.get_pwm_score(k_mer, expPwmR, k_nb)
            bias_table_R[k_mer] = round(obsR / expR, 6)

        # Return
        return [bias_table_F, bias_table_R]

# if __name__ == "__main__":
#  import sys
#  from rgt.GenomicRegionSet import *
#  bam_file=sys.argv[1]
#  fasta_file=sys.argv[2]
#  bed_file=sys.argv[3]
#  kmer=int(sys.argv[4])
#  shift=int(sys.argv[5])
#  out=sys.argv[6]
#  regions=GenomicRegionSet("regions")
#  regions.read_bed(bed_file)
#  table=BiasTable(regions=regions,dnase_file_name=bam_file,genome_file_name=fasta_file,k_nb=kmer,shift=shift)
#  table.write_tables(out)
