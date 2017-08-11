import os
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
from scipy import stats
from collections import OrderedDict
from rgt.tdf.RNADNABindingSet import RNADNABindingSet
from rgt.motifanalysis.Statistics import multiple_test_correction
from rgt.tdf.triplexTools import dbd_regions
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.Util import OverlapType
from rgt.tdf.triplexTools import rna_associated_gene

class Statistics(object):
    def __init__(self, pars):
        self.pars = pars
        self.sig_DBD = []
        self.stat = { "name":self.pars.rn, "genome":self.pars.organism,
                      "exons":1, "seq_length": None,
                      "target_regions": 0, "background_regions": 0,
                      "DBD_all": 0, "DBD_sig": 0,
                      "DBSs_target_all": 0, "DBSs_target_DBD_sig": 0,
                      "DBSs_background_all": 0, "DBSs_background_DBD_sig": 0, "p_value": "-",
                      "associated_gene": ".", "expression": "n.a.", "loci": "-", "autobinding": 0,
                      "MA_G": 0, "MA_T": 0, "MP_G": 0, "MP_T": 0, "RA_A": 0, "RA_G": 0, "YP_C": 0, "YP_T": 0,
                      "uniq_MA_G": 0, "uniq_MA_T": 0, "uniq_MP_G": 0, "uniq_MP_T": 0,
                      "uniq_RA_A": 0, "uniq_RA_G": 0, "uniq_YP_C": 0, "uniq_YP_T": 0 }

    def count_frequency_promoters(self, target_regions, background, file_tpx_de, file_tpx_nde):
        # count_frequency(self, temp, remove_temp, cutoff, l, obedp=False):
        """Count the frequency between DE genes and non-DE genes with the given BindingSiteSet"""

        len_de = len(target_regions)
        len_nde = len(background)
        # print([len_de, len_nde])

        self.frequency = {}
        self.frequency["promoters"] = {"de": OrderedDict(), "nde": OrderedDict()}

        ########################################################
        # Count the number of targeted promoters on each merged DNA Binding Domain
        print("\tCounting frequency of promoters on DBD...")
        self.tpx_de = RNADNABindingSet("DE")

        self.tpx_de.read_tpx(file_tpx_de, dna_fine_posi=False)
        self.tpx_de.remove_duplicates()
        self.tpx_de.merge_rbs(rm_duplicate=True, cutoff=self.pars.ccf,
                              region_set=target_regions, name_replace=target_regions)
        self.rbss = self.tpx_de.merged_dict.keys()

        for rbs in self.rbss:
            # DE
            self.tpx_de.merged_dict[rbs].remove_duplicates()
            l1 = len(self.tpx_de.merged_dict[rbs])
            self.frequency["promoters"]["de"][rbs] = [l1, len_de - l1]

        self.tpx_nde = RNADNABindingSet("non-DE")
        self.tpx_nde.read_tpx(file_tpx_nde, dna_fine_posi=False)

        self.tpx_nde.remove_duplicates()
        self.tpx_nde.merge_rbs(rbss=self.rbss, region_set=background, rm_duplicate=True)

        for rbs in self.rbss:
            l2 = len(self.tpx_nde.merged_dict[rbs])
            self.frequency["promoters"]["nde"][rbs] = [l2, len_nde - l2]

        print("\t\t" + str(len(self.tpx_de)) + "\tBinding target promoters")
        print("\t\t" + str(len(self.tpx_nde)) + "\tBinding non-target promoters")
        # self.stat["DBSs_target_all"] = str(len(self.txp_de))
        # self.stat["DBSs_background_all"] = str(len(self.txp_nde))

        ########################################################
        # Count the number of hits on the promoters from each merged DBD
        print("\tCounting frequency of binding sites on DBD...")
        ####################
        # DE
        self.tpx_def = RNADNABindingSet("DE")
        self.tpx_def.read_tpx(file_tpx_de, dna_fine_posi=True, seq=True)
        self.tpx_def.merge_rbs(rbss=self.rbss, rm_duplicate=True, name_replace=target_regions)


        # Promoter profiling
        self.promoter = {"de": {}, "nde": {}}
        self.promoter["de"]["rd"] = self.tpx_def.sort_rd_by_regions(regionset=target_regions)
        self.promoter["de"]["dbs"] = {}
        self.promoter["de"]["dbs_coverage"] = {}

        for promoter in target_regions:
            dbs = self.promoter["de"]["rd"][promoter.toString()].get_dbs()
            m_dbs = dbs.merge(w_return=True)
            self.promoter["de"]["dbs"][promoter.toString()] = len(dbs)
            # self.promoter["de"]["merged_dbs"][promoter.toString()] = len(m_dbs)
            self.promoter["de"]["dbs_coverage"][promoter.toString()] = float(m_dbs.total_coverage()) / len(promoter)
        print("\t\t" + str(len(self.tpx_def)) + "\tBinding sites on de promoters")

        ######################
        # nDE
        self.promoter["nde"]["dbs"] = {}
        self.promoter["nde"]["dbs_coverage"] = {}

        self.tpx_nde = RNADNABindingSet("non-DE")
        self.tpx_nde.read_tpx(file_tpx_nde, dna_fine_posi=False)
        self.tpx_nde.merge_rbs(rbss=self.rbss, rm_duplicate=True)
        ndef_dbs = self.tpx_nde.get_dbs()

        counts = background.counts_per_region(regionset=ndef_dbs)
        coverage = background.coverage_per_region(regionset=ndef_dbs)
        for i, p in enumerate(target_regions):
            self.promoter["nde"]["dbs"][p.toString()] = counts[i]
            self.promoter["nde"]["dbs_coverage"][p.toString()] = coverage[i]

        if self.pars.showdbs:
            self.frequency["hits"] = {"de": OrderedDict(), "nde": OrderedDict()}
            numdbs_def = len(self.tpx_def.get_dbs(rm_duplicate=True))
            numdbs_ndef = len(self.tpx_ndef.get_dbs(rm_duplicate=True))
            for rbs in self.rbss:
                # DE
                l1 = len(self.tpx_def.merged_dict[rbs])
                self.frequency["hits"]["de"][rbs] = [l1, numdbs_def - l1]
                # non-DE
                l2 = len(self.tpx_ndef.merged_dict[rbs])
                self.frequency["hits"]["nde"][rbs] = [l2, numdbs_ndef - l2]

        if self.pars.rm:
            os.remove(file_tpx_de)
        os.remove(file_tpx_nde)

        target_regions.write(filename=os.path.join(self.pars.o, self.pars.rn + "_target_promoters.bed"))
        self.tpx_de.write_bed(filename=os.path.join(self.pars.o, self.pars.rn + "_target_promoters_dbs.bed"),
                              associated=self.pars.organism)
        self.tpx_def.write_bed(filename=os.path.join(self.pars.o, self.pars.rn + "_dbss.bed"),
                               remove_duplicates=False, associated=self.pars.organism)
    def fisher_exact_de(self,):
        """Return oddsratio and pvalue"""
        self.oddsratio = {}
        self.pvalue = {}
        pvalues = []

        for rbs in self.frequency["promoters"]["de"]:
            table = numpy.array([self.frequency["promoters"]["de"][rbs], self.frequency["promoters"]["nde"][rbs]])
            # print(table)
            self.oddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
            pvalues.append(p)

        # correction
        if len(self.frequency["promoters"]["de"]) > 1:
            reject, pvals_corrected = multiple_test_correction(pvalues, alpha=self.pars.a, method='indep')
        else:
            pvals_corrected = pvalues
        for i, rbs in enumerate(self.frequency["promoters"]["de"]):
            self.pvalue[rbs] = pvals_corrected[i]
            if pvals_corrected[i] < self.pars.a:
                self.sig_DBD.append(rbs)

        try:
            self.min_p_value = str(min(pvals_corrected))
        except:
            self.min_p_value = "1"

        if self.pars.showdbs:
            self.hoddsratio = {}
            self.hpvalue = {}
            pvalues = []
            for rbs in self.frequency["hits"]["de"]:
                table = numpy.array([self.frequency["hits"]["de"][rbs], self.frequency["hits"]["nde"][rbs]])
                self.hoddsratio[rbs], p = stats.fisher_exact(table, alternative="greater")
                pvalues.append(p)
            # correction
            if len(self.frequency["hits"]["de"]) > 1:
                reject, pvals_corrected = multiple_test_correction(pvalues, alpha=self.pars.a, method='indep')
            else:
                pvals_corrected = pvalues
            for i, rbs in enumerate(self.frequency["hits"]["de"]):
                self.hpvalue[rbs] = pvals_corrected[i]
                if pvals_corrected[i] < self.pars.a:
                    self.sig_DBD.append(rbs)

    def dbd_regions(self, rna_exons):
        dbd_regions(exons=rna_exons, sig_region=self.sig_DBD,
                    rna_name=self.pars.rn, output=self.pars.o)


    def summary_stat(self, input):
        self.stat["target_regions"] = str(len(input.dna.target_regions))
        self.stat["background_regions"] = str(len(input.dna.nontarget_regions))
        if input.rna.rna_regions:
            r_genes = rna_associated_gene(rna_regions=input.rna.rna_regions,
                                          name=self.pars.rn, organism=self.pars.organism)
            self.stat["associated_gene"] = r_genes
            self.stat["loci"] = input.rna.rna_regions[0][0] + ":" + str(input.rna.rna_regions[0][1]) + "-" + \
                                str(input.rna.rna_regions[-1][2]) + "_" + input.rna.rna_regions[0][3]
        else:
            self.stat["associated_gene"] = "."
            self.stat["loci"] = "-"
        self.stat["expression"] = str(input.rna.rna_expression)
        self.stat["exons"] = input.rna.num_exons
        self.stat["seq_length"] = input.rna.seq_length
        self.stat["DBSs_target_all"] = str(len(self.tpx_de))
        self.stat["DBSs_background_all"] = str(len(self.tpx_nde))
        self.stat["p_value"] = self.min_p_value
        self.stat["DBD_all"] = str(len(self.rbss))
        self.stat["DBD_sig"] = str(len(self.sig_DBD))
        sigDBD = GenomicRegionSet("DBD_sig")
        sigDBD.sequences = self.sig_DBD
        rbss = self.tpx_de.get_rbs()
        overlaps = rbss.intersect(y=sigDBD, mode=OverlapType.ORIGINAL)
        self.stat["DBSs_target_DBD_sig"] = str(len(overlaps))
        rbss = self.tpx_nde.get_rbs()
        overlaps = rbss.intersect(y=sigDBD, mode=OverlapType.ORIGINAL)
        self.stat["DBSs_background_DBD_sig"] = str(len(overlaps))

        self.stat["autobinding"] = len(self.autobinding)