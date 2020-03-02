# Python Libraries


import os
import numpy
numpy.seterr(divide='ignore', invalid='ignore')
from scipy import stats
import multiprocessing
from collections import OrderedDict
from ..tdf.RNADNABindingSet import RNADNABindingSet
from ..motifanalysis.Statistics import multiple_test_correction
from ..tdf.triplexTools import dbd_regions
from ..GenomicRegionSet import GenomicRegionSet
from ..Util import OverlapType
from ..tdf.triplexTools import rna_associated_gene
from ..tdf.Triplexes import random_each




class Statistics(object):
    def __init__(self, pars):
        self.pars = pars
        self.sig_DBD = []
        self.stat = { "name": self.pars.rn,
                      "genome": self.pars.organism,
                      "exons": 1,
                      "seq_length": None,
                      "target_regions": 0,
                      "background_regions": 0,
                      "DBD_all": 0,
                      "DBD_sig": 0,
                      "DBSs_target_all": 0,
                      "DBSs_target_DBD_sig": 0,
                      "DBSs_background_all": 0,
                      "DBSs_background_DBD_sig": 0,
                      "p_value": "1",
                      "associated_gene": ".",
                      "expression": "n.a.",
                      "loci": "-",
                      "autobinding": 0,
                      "MA_G": 0, "MA_T": 0, "MP_G": 0, "MP_T": 0,
                      "RA_A": 0, "RA_G": 0, "YP_C": 0, "YP_T": 0,
                      "uniq_MA_G": 0, "uniq_MA_T": 0, "uniq_MP_G": 0, "uniq_MP_T": 0,
                      "uniq_RA_A": 0, "uniq_RA_G": 0, "uniq_YP_C": 0, "uniq_YP_T": 0,
                      "target_in_trans": 0, "target_in_cis": 0, "target_local": 0,
                      "background_in_trans": 0, "background_in_cis": 0, "background_local": 0}

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
        print("\t\t" + str(len(self.tpx_de)) + "\tBinding target promoters")
        self.tpx_de.remove_duplicates()
        self.tpx_de.merge_rbs(rm_duplicate=True, cutoff=self.pars.ccf,
                              region_set=target_regions, name_replace=target_regions)
        self.rbss = list(self.tpx_de.merged_dict.keys())

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

        if self.pars.rt:
            os.remove(file_tpx_de)
        os.remove(file_tpx_nde)

        # target_regions.write(filename=os.path.join(self.pars.o, self.pars.rn + "_target_promoters.bed"))
        self.tpx_de.write_bed(filename=os.path.join(self.pars.o, self.pars.rn + "_target_promoters_dbs.bed"),
                              associated=self.pars.organism)
        # self.tpx_def.write_bed(filename=os.path.join(self.pars.o, self.pars.rn + "_dbss.bed"),
        #                        remove_duplicates=False, associated=self.pars.organism)
    def fisher_exact_de(self):
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


    def summary_stat(self, input, triplexes, mode, no_binding=False):
        """Summerize the statistics according to its mode: promotertest or regiontest."""
        self.stat["target_regions"] = str(len(input.dna.target_regions))

        if mode == "promotertest":
            self.stat["background_regions"] = str(len(input.dna.nontarget_regions))
            self.stat["DBSs_target_all"] = str(len(self.tpx_de))

            if not no_binding:
                self.stat["DBSs_background_all"] = str(len(self.tpx_nde))
                self.stat["p_value"] = self.min_p_value

                self.stat["DBD_sig"] = str(len(self.sig_DBD))
                sigDBD = GenomicRegionSet("DBD_sig")
                sigDBD.sequences = self.sig_DBD
                rbss = self.tpx_de.get_rbs()
                overlaps = rbss.intersect(y=sigDBD, mode=OverlapType.ORIGINAL)
                self.stat["DBSs_target_DBD_sig"] = str(len(overlaps))
                rbss = self.tpx_nde.get_rbs()
                overlaps = rbss.intersect(y=sigDBD, mode=OverlapType.ORIGINAL)
                self.stat["DBSs_background_DBD_sig"] = str(len(overlaps))
        elif mode == "regiontest":
            self.stat["background_regions"] = str(len(input.dna.target_regions))
            if not no_binding:
                self.stat["DBSs_target_all"] = str(len(self.tpxf))
                self.stat["DBSs_background_all"] = str(numpy.mean(self.counts_dbss))
                self.stat["p_value"] = str(min(self.data["region"]["p"]))

                self.stat["DBD_sig"] = str(len(self.data["region"]["sig_region"]))
                sigDBD = GenomicRegionSet("DBD_sig")
                sigDBD.sequences = self.data["region"]["sig_region"]
                rbss = self.tpx.get_rbs()
                overlaps = rbss.intersect(y=sigDBD, mode=OverlapType.ORIGINAL)
                self.stat["DBSs_target_DBD_sig"] = str(len(overlaps))



        if input.rna.regions:
            r_genes = rna_associated_gene(rna_regions=input.rna.regions,
                                          name=self.pars.rn, organism=self.pars.organism)
            self.stat["associated_gene"] = r_genes
            self.stat["loci"] = input.rna.regions[0][0] + ":" + str(input.rna.regions[0][1]) + "-" + \
                                str(input.rna.regions[-1][2]) + "_" + input.rna.regions[0][3]
        # else:
        #     self.stat["associated_gene"] = "."
        #     self.stat["loci"] = "-"

        self.stat["expression"] = str(input.rna.expression)
        self.stat["exons"] = input.rna.num_exons
        self.stat["seq_length"] = input.rna.seq_length

        if not no_binding:

            self.stat["DBD_all"] = str(len(self.rbss))
            self.stat["DBD_sig"] = str(len(self.sig_DBD))
        self.stat["autobinding"] = len(triplexes.autobinding)

    def output_bed(self, input, tpx):
        # print(input.rna.regions)
        if len(input.rna.regions) > 0:
            # print(self.rna_regions)
            rna_regionsets = GenomicRegionSet(name=self.pars.rn)
            rna_regionsets.load_from_list(input.rna.regions)

            autobinding_loci = tpx.get_overlapping_regions(regionset=rna_regionsets)
            autobinding_loci.write(filename=os.path.join(self.pars.o, self.pars.rn+"_autobinding.bed"))

        tpx.write_bed(filename=os.path.join(self.pars.o, self.pars.rn + "_dbss.bed"),
                      remove_duplicates=False, associated=self.pars.organism)


    def write_stat(self, filename):
        """Write the statistics into file"""
        order_stat = ["title", "name", "genome",
                      "exons", "seq_length",
                      "target_regions", "background_regions",
                      "DBD_all", "DBD_sig",
                      "DBSs_target_all", "DBSs_target_DBD_sig",
                      "DBSs_background_all", "DBSs_background_DBD_sig", "p_value",
                      "Norm_DBD", "Norm_DBS", "Norm_DBS_sig",
                      "associated_gene", "expression", "loci", "autobinding",
                      "MA_G", "MA_T", "MP_G", "MP_T", "RA_A", "RA_G", "YP_C", "YP_T",
                      "uniq_MA_G", "uniq_MA_T", "uniq_MP_G", "uniq_MP_T",
                      "uniq_RA_A", "uniq_RA_G", "uniq_YP_C", "uniq_YP_T",
                      "target_in_trans", "target_in_cis", "target_local",
                      "background_in_trans", "background_in_cis", "background_local"]

        with open(filename, "w") as f:
            for k in order_stat:
                try:
                    print("\t".join([k, str(self.stat[k])]), file=f)
                except:
                    continue




    def dbs_motif(self, tpx):
        tpx.motif_statistics()
        for i, mode in enumerate(tpx.motifs):
            for con in tpx.motifs[mode]:
                self.stat[mode+"_"+con] = str(tpx.motifs[mode][con])

    def uniq_motif(self, tpx, rnalen):
        tpx.uniq_motif_statistics(rnalen=rnalen)
        for k, v in tpx.uniq_motifs.items():
            self.stat[k] = sum(v)
            self.stat["uniq_" + k] = sum([1 for x in v if x > 0])

    def random_test(self, repeats, target_regions, filter_bed, mp, genome_fasta):
        """Perform randomization for the given times"""
        self.repeats = repeats
        marks = numpy.round(numpy.linspace(0, repeats - 1, num=41)).tolist()
        # Prepare the input lists for multiprocessing
        mp_input = []
        for i in range(repeats):
            mp_input.append([str(i), os.path.join(self.pars.o, "rna_temp.fa"), target_regions,
                             self.pars.o, self.pars.organism, self.rbss, str(marks.count(i)),
                             str(self.pars.l), str(self.pars.e), str(self.pars.c), str(self.pars.fr),
                             str(self.pars.fm), str(self.pars.of), str(self.pars.mf), 0,
                             filter_bed, genome_fasta, self.pars.par])
        # Multiprocessing
        print("\t\t|0%                  |                100%|")
        print("\t\t[", end="")
        pool = multiprocessing.Pool(processes=mp)
        mp_output = pool.map(random_each, mp_input)
        # print(mp_output)
        pool.close()
        pool.join()
        print("]")

        # Processing the result
        self.region_matrix = []
        self.dbss_matrix = []
        self.data = {"region": {"ave": [],
                                "sd": [],
                                "p": [],
                                "sig_region": [],
                                "sig_boolean": []},
                     "dbs": {"ave": [],
                             "sd": [],
                             "p": [],
                             "sig_region": [],
                             "sig_boolean": []}}

        region_counts = [v[0] for v in mp_output]
        dbss_counts = [v[1] for v in mp_output]

        for i, rbs in enumerate(self.rbss):

            counts_regions = [v[i] for v in region_counts]

            self.data["region"]["ave"].append(numpy.mean(counts_regions))
            self.data["region"]["sd"].append(numpy.std(counts_regions))
            num_sig = len([h for h in counts_regions if h > self.counts_tr[rbs][0]])
            p_region = float(num_sig) / repeats
            self.data["region"]["p"].append(p_region)
            self.region_matrix.append(counts_regions)

            if p_region < self.pars.a:
                self.sig_DBD.append(rbs)
                self.data["region"]["sig_region"].append(rbs)
                self.data["region"]["sig_boolean"].append(True)
            else:
                self.data["region"]["sig_boolean"].append(False)

            try:
                if p_region < self.topDBD[1]:
                    self.topDBD = [rbs.str_rna(pa=False), p_region]
            except:
                self.topDBD = [rbs.str_rna(pa=False), p_region]

            # Analysis based on DBSs
            if self.pars.showdbs:
                counts_dbss = [v[i] for v in dbss_counts]

                self.data["dbs"]["ave"].append(numpy.mean(counts_dbss))
                self.data["dbs"]["sd"].append(numpy.std(counts_dbss))
                num_sig = len([h for h in counts_dbss if h > self.counts_dbs[rbs]])
                p_dbs = float(num_sig) / repeats
                self.data["dbs"]["p"].append(p_dbs)
                self.dbss_matrix.append(counts_dbss)
                if p_dbs < self.pars.a:
                    self.data["dbs"]["sig_region"].append(rbs)
                    self.data["dbs"]["sig_boolean"].append(True)
                else:
                    self.data["dbs"]["sig_boolean"].append(False)

        self.region_matrix = numpy.array(self.region_matrix)

        # if self.pars.showdbs:
        # self.dbss_matrix = numpy.array(self.dbss_matrix)

        self.counts_dbss = [v[i] for v in dbss_counts]
        # self.stat["DBSs_background_all"] = numpy.mean(counts_dbss)
        # try: self.stat["p_value"] = str(min(self.data["region"]["p"]))
        # except: self.stat["p_value"] = "1"

        with open(os.path.join(self.pars.o, "counts_random_matrix.txt"), "w") as f:
            for l in self.region_matrix:
                print("\t".join([str(x) for x in l]), file=f)
        with open(os.path.join(self.pars.o, "counts_dbs.txt"), "w") as f:
            print("\t".join([str(x) for x in list(self.counts_dbs.values())]), file=f)

        # Integarte distances
        self.stat["background_in_trans"] = float(sum([v[2]["in_trans"] for v in mp_output])) / len(mp_output)
        self.stat["background_in_cis"] = float(sum([v[2]["in_cis"] for v in mp_output])) / len(mp_output)
        self.stat["background_local"] = float(sum([v[2]["local"] for v in mp_output])) / len(mp_output)

    def target_stat(self, target_regions, tpx, tpxf):
        # self.stat["DBSs_target_all"] = str(len(self.txpf))
        tpx.merge_rbs(rm_duplicate=True, region_set=target_regions,
                      asgene_organism=self.pars.organism, cutoff=self.pars.ccf)
        self.rbss = list(tpx.merged_dict.keys())
        self.counts_tr = OrderedDict()
        self.counts_dbs = OrderedDict()

        self.tpxf = tpxf
        tpxf.merge_rbs(rm_duplicate=True, region_set=target_regions, rbss=self.rbss,
                       asgene_organism=self.pars.organism)

        for rbs in self.rbss:
            tr = len(tpx.merged_dict[rbs])
            self.counts_tr[rbs] = [tr, len(target_regions) - tr]
            self.counts_dbs[rbs] = len(tpx.merged_dict[rbs])
        self.region_dbd = tpxf.sort_rbs_by_regions(target_regions)
        self.region_dbs = tpxf.sort_rd_by_regions(regionset=target_regions)
        self.region_normdbs = {}
        for t in target_regions:
            self.region_normdbs[t.toString()] = len(self.region_dbs[t.toString()]) * 1000 / len(t)
        self.region_dbsm = {}
        self.region_coverage = {}

        for region in target_regions:
            self.region_dbsm[region.toString()] = self.region_dbs[region.toString()].get_dbs().merge(w_return=True)
            self.region_coverage[region.toString()] = float(self.region_dbsm[region.toString()].total_coverage()) / len(region)
        # self.stat["target_regions"] = str(len(target_regions))

    def distance_distribution(self, tpx):
        dis_count = tpx.distance_distribution()
        self.stat["target_in_trans"] = dis_count["in_trans"]
        self.stat["target_in_cis"] = dis_count["in_cis"]
        self.stat["target_local"] = dis_count["local"]