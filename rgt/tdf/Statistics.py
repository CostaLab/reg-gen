
class Statistics(object):
    def __init__(self, pars):
        self.pars = pars
    def count_frequency(self):
        # count_frequency(self, temp, remove_temp, cutoff, l, obedp=False):
        """Count the frequency between DE genes and non-DE genes with the given BindingSiteSet"""

        len_de = len(self.de_regions)
        len_nde = len(self.nde_regions)
        # print([len_de, len_nde])

        self.frequency = {}
        self.frequency["promoters"] = {"de": OrderedDict(), "nde": OrderedDict()}

        ########################################################
        # Count the number of targeted promoters on each merged DNA Binding Domain
        print("\tCounting frequency of promoters on DBD...")
        self.txp_de = RNADNABindingSet("DE")

        self.txp_de.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=False)
        print("\t\t" + str(len(self.txp_de)) + "\tBinding de promoters")
        self.txp_de.remove_duplicates()
        self.txp_de.merge_rbs(rm_duplicate=True, cutoff=cutoff,
                              region_set=self.de_regions, name_replace=self.de_regions)

        self.rbss = self.txp_de.merged_dict.keys()

        for rbs in self.rbss:
            # DE
            self.txp_de.merged_dict[rbs].remove_duplicates()
            l1 = len(self.txp_de.merged_dict[rbs])
            self.frequency["promoters"]["de"][rbs] = [l1, len_de - l1]

        self.txp_nde = RNADNABindingSet("non-DE")
        self.txp_nde.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=False)
        print("\t\t" + str(len(self.txp_nde)) + "\tBinding nde promoters")
        self.txp_nde.remove_duplicates()
        self.txp_nde.merge_rbs(rbss=self.rbss, region_set=self.nde_regions,
                               rm_duplicate=True)  # , asgene_organism=self.organism)

        for rbs in self.rbss:
            l2 = len(self.txp_nde.merged_dict[rbs])
            self.frequency["promoters"]["nde"][rbs] = [l2, len_nde - l2]

        self.stat["DBSs_target_all"] = str(len(self.txp_de))
        self.stat["DBSs_background_all"] = str(len(self.txp_nde))

        ########################################################
        # Count the number of hits on the promoters from each merged DBD

        print("\tCounting frequency of binding sites on DBD...")
        ####################
        # DE
        self.txp_def = RNADNABindingSet("DE")
        self.txp_def.read_txp(os.path.join(temp, "de.txp"), dna_fine_posi=True, seq=True)
        self.txp_def.merge_rbs(rbss=self.rbss, rm_duplicate=True,
                               name_replace=self.de_regions)  # asgene_organism=self.organism
        print("\t\t" + str(len(self.txp_def)) + "\tBinding sites on de promoters")

        # Promoter profiling
        self.promoter = {"de": {}, "nde": {}}

        self.promoter["de"]["rd"] = self.txp_def.sort_rd_by_regions(regionset=self.de_regions)
        # self.promoter["de"]["merged_dbs"] = {}
        self.promoter["de"]["dbs"] = {}
        self.promoter["de"]["dbs_coverage"] = {}

        for promoter in self.de_regions:
            dbs = self.promoter["de"]["rd"][promoter.toString()].get_dbs()
            m_dbs = dbs.merge(w_return=True)
            self.promoter["de"]["dbs"][promoter.toString()] = len(dbs)
            # self.promoter["de"]["merged_dbs"][promoter.toString()] = len(m_dbs)
            self.promoter["de"]["dbs_coverage"][promoter.toString()] = float(m_dbs.total_coverage()) / len(promoter)

        ######################
        # nDE
        self.promoter["nde"]["dbs"] = {}
        self.promoter["nde"]["dbs_coverage"] = {}

        # ndef_dbs = GenomicRegionSet("nde_dbs")

        self.txp_nde = RNADNABindingSet("non-DE")
        self.txp_nde.read_txp(os.path.join(temp, "nde.txp"), dna_fine_posi=False)
        self.txp_nde.merge_rbs(rbss=self.rbss, rm_duplicate=True)  # , asgene_organism=self.organism)

        ndef_dbs = self.txp_nde.get_dbs()

        counts = self.nde_regions.counts_per_region(regionset=ndef_dbs)

        coverage = self.nde_regions.coverage_per_region(regionset=ndef_dbs)
        for i, p in enumerate(self.de_regions):
            self.promoter["nde"]["dbs"][p.toString()] = counts[i]
            # self.promoter["nde"]["merged_dbs"][p.toString()] = mcounts[i]
            self.promoter["nde"]["dbs_coverage"][p.toString()] = coverage[i]

        if self.showdbs:
            self.frequency["hits"] = {"de": OrderedDict(), "nde": OrderedDict()}
            numdbs_def = len(self.txp_def.get_dbs(rm_duplicate=True))
            numdbs_ndef = len(self.txp_ndef.get_dbs(rm_duplicate=True))
            for rbs in self.rbss:
                # DE
                l1 = len(self.txp_def.merged_dict[rbs])
                self.frequency["hits"]["de"][rbs] = [l1, numdbs_def - l1]
                # non-DE
                l2 = len(self.txp_ndef.merged_dict[rbs])
                self.frequency["hits"]["nde"][rbs] = [l2, numdbs_ndef - l2]

        if remove_temp:
            os.remove(os.path.join(temp, "de.txp"))

        os.remove(os.path.join(temp, "nde.txp"))

        if obedp:
            try:
                output = self.de_regions.change_name_by_dict(convert_dict=self.ensembl2symbol)
            except:
                output = self.de_regions
            output.write(filename=os.path.join(temp, obedp + "_target_promoters.bed"))

            self.txp_de.write_bed(filename=os.path.join(temp, obedp + "_target_promoters_dbs.bed"),
                                  associated=self.organism)
            self.txp_def.write_bed(filename=os.path.join(temp, obedp + "_dbss.bed"),
                                   remove_duplicates=False, associated=self.organism)
