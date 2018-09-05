# Python Libraries
from __future__ import print_function
import os
import numpy
import os
# Local Libraries

# Distal Libraries
from ..GenomicRegionSet import *
from ..Util import OverlapType
from .BindingSiteSet import BindingSite, BindingSiteSet

class RNADNABinding:
    """Describe a binding event of RNA and DNA to form a triplex 
    
    Attributes:
        self.rna        A BindingSite object which includes RNA name, \
                        region (start position and end position), \
                        score, error and its Sequence
        self.dna        A BindingSite object which includes DNA name, \
                        region (start position and end position), \
                        score, error and its Sequence
        self.score      The score from comparing RNA region and DNA region 
        self.err_rate   The ration of error base pair between RNA and DNA
        self.err        *************
        self.motif      The motif of triplex formation:
                        R - the purine motif that permit guanines (G) and adenines (A)
                        Y - the pyrimidine motif that permit cytosines (C) and thymines (T)
                        M - the mixed motif, purine-pyrimdine, that permit guanines (G) and thymines (T)
                        P - parallel binding, i.e. motifs facilitating Hoogsten bonds; 
                            this covers the pyrimidine motif and the purine-pyrimidine motif in parallel configuration
                        A - anti-parallel binding, i.e. motifs facilitating reverse Hoogsten bonds; 
                            this covers the purine motif and the purine-pyrimidine motif in anti-parallel configuration
        self.strand     The strand of DNA
        self.orient     The orientation of RNA regarding to DNA
        self.guan_rate  *************
        self.rna_seq    A string of RNA extended by gaps to match the position of DNA
        self.dna_seq    A string of DNA extended by gaps to match the position of RNA
        self.match      A string with '|' for perfect match, '*' for mismatch
    
    """

    __slots__ = ['rna', 'dna', 'motif', 'orient', 'score', 'err_rate', 'err', 'guan_rate', 
                 'match' ]

    def __init__(self, rna, dna, score, err_rate, err, guan_rate, match=None, rna_gap=None, dna_gap=None):
        """Initiation"""
        self.rna = rna  # BindingSite
        self.dna = dna  # GenomicRegion
        self.motif = rna.motif
        self.strand = dna.orientation
        self.orient = rna.orientation
        self.score = score
        self.err_rate = err_rate
        self.err = err
        self.guan_rate = guan_rate
        self.match = match
        self.rna_gap = rna_gap
        self.dna_gap = dna_gap
    
    def __str__(self):
        return "\t".join([self.rna.chrom, str(self.rna.initial), str(self.rna.final), self.dna.toString(), "0",
                            str(abs(self.dna.final - self.dna.initial)), str(self.score), str(self.err_rate), str(self.err), str(self.motif),
                            str(self.strand), str(self.orient), str(self.guan_rate), str(self.rna_gap), str(self.dna_gap)])

    def __eq__(self, other):
        return (self.dna, self.rna) == (other.dna, other.rna)

    def motif_statistics(self, cm=False):
        if not self.match:
            return None
        else:
            res = {"Mix_Antiparallel": {"G": 0, "T": 0},
                    "Mix_Parallel": {"G": 0, "T": 0},
                    "Purine_Antiparallel": {"A": 0, "G": 0},
                    "Pyrimidine_Parallel": {"C": 0, "T": 0}}
            if cm:
                res_c = {"Mix_Purine_AA": {"G": 0, "T": 0, "A":0},
                         "Purine_Mix_AA": {"G": 0, "T": 0, "A": 0},
                         "Pyrimidine_Mix_PP": {"G": 0, "T": 0, "C": 0},
                         "Mix_Pyrimidine_PP": {"G": 0, "T": 0, "C": 0},
                         "Mix_Pyrimidine_PA": {"G": 0, "T": 0, "C": 0},
                         "Pyrimidine_Mix_PA": {"G": 0, "T": 0, "C": 0},
                         "Purine_Mix_PA": {"A": 0, "G": 0, "T": 0},
                         "Mix_Purine_PA": {"A": 0, "G": 0, "T": 0},
                         "Purine_Pyrimidine_PA": {"A": 0, "G": 0, "C": 0, "T": 0},
                         "Pyrimidine_Purine_PA": {"A": 0, "G": 0, "C": 0, "T": 0},
                         "Mix_Mix_PA": {"G": 0, "T": 0}
                         }

                res.update(res_c)
            if self.strand == "+":
                rna = self.match[0].split()[1]
                linking = self.match[1].strip()
                # dnap = self.match[2].split()[1]
                # dnan = self.match[3].split()[1]
            elif self.strand == "-":
                rna = self.match[3].split()[1]
                linking = self.match[2].strip()
                # dnap = self.match[0].split()[1]
                # dnan = self.match[1].split()[1]
            elif self.strand == "+_-":
                rna = self.match[0].split()[1] + self.match[3].split()[3]
                linking = self.match[2].strip().split()[3] + self.match[1].strip().split()[0]
            elif self.strand == "-_+":
                rna = self.match[3].split()[1] + self.match[0].split()[3]
                linking = self.match[1].strip().split()[3] + self.match[2].strip().split()[0]
            elif self.strand == "+_+":
                rna = self.match[0].split()[1] + self.match[0].split()[3]
                linking = self.match[1].strip().split()[0] + self.match[1].strip().split()[1]
            elif self.strand == "-_-":
                rna = self.match[3].split()[1] + self.match[3].split()[3]
                linking = self.match[2].strip().split()[0] + self.match[2].strip().split()[1]

            if self.motif == "M":
                motif = "Mix"
            elif self.motif == "R":
                motif = "Purine"
            elif self.motif == "Y":
                motif = "Pyrimidine"

            if self.orient == "P":
                orient = "Parallel"
            elif self.orient == "A":
                orient = "Antiparallel"
            if cm:
                if self.motif == "M_R":
                    motif = "Mix_Purine"
                elif self.motif == "R_M":
                    motif = "Purine_Mix"
                elif self.motif == "M_Y":
                    motif = "Mix_Pyrimidine"
                elif self.motif == "Y_M":
                    motif = "Pyrimidine_Mix"
                elif self.motif == "M_M":
                    motif = "Mix_Mix"
                elif self.motif == "Y_R":
                    motif = "Pyrimidine_Purine"
                elif self.motif == "R_Y":
                    motif = "Purine_Pyrimidine"

            if cm:
                if self.orient == "P_P":
                    orient = "PP"
                elif self.orient == "A_A":
                    orient = "AA"
                elif self.orient == "A_P" or self.orient == "P_A":
                    orient = "PA"
            #########################################333
            # Overall counts
            for i, bp in enumerate(linking):
                if bp == "|":
                    if rna[i].upper() == "A":
                        res[motif + "_" + orient]["A"] += 1
                    elif rna[i].upper() == "G":
                        res[motif + "_" + orient]["G"] += 1
                    elif rna[i].upper() == "T":
                        res[motif + "_" + orient]["T"] += 1
                    elif rna[i].upper() == "C":
                        res[motif + "_" + orient]["C"] += 1
            return res


class RNADNABindingSet:
    """A framework to map DNA binding sites to corresponding RNA binding sites 
    
    Attributes:
        self.rna    A BindingSequence-ID   
        [1] TFO start   
        [2] TFO end 
        [3] Duplex-ID   
        [4] TTS start   
        [5] TTS end 
        [6] Score   
        [7] Error-rate  
        [8] Errors  
        [9] Motif   
        [10] Strand  
        [11] Orientation 
        [12] Guanine-rate
    
    """
    __slots__ = ['name', 'sequences', 'sorted_dna', 'sorted_rna' ]

    def __init__(self, name):
        """Initialize"""
        self.name = name       # RNA name
        self.sequences = []    # A list to contain all RNA/DNA interactions and extra information 
        self.sorted_dna = False
        self.sorted_rna = False
 

    def __len__(self):
        """Return the number of interactions between DNA and RNA """
        return len(self.sequences)

    def __iter__(self):
        return iter(self.sequences)

    def __getitem__(self, key):
        return self.sequences[key]

    def add(self, rnadnabinding):
        """Add one pair of BindingSite on RNA and GenomicRegion on DNA"""
        self.sequences.append(rnadnabinding)
        
    def get_rbs(self, sort=False, orientation=None):
        """Return a BindingSiteSet which contains all RNA binding sites"""
        rna_set = BindingSiteSet(name=self.name)
        for rd in self.sequences:
            if not orientation:
                rna_set.add(rd.rna)
            elif orientation == rd.orient:
                rna_set.add(rd.rna)
            else: pass

        if sort: rna_set.sort()
        return rna_set

    def get_dbs(self, sort=False, orientation=None, rm_duplicate=False, dbd_tag=False):
        """Return GenomicRegionSet which contains all DNA binding sites"""
        dna_set = GenomicRegionSet(name="DNA_binding_sites")
        if len(self) == 0: return dna_set
        for rd in self.sequences:
            if dbd_tag:
                dbs = GenomicRegion(chrom=rd.dna.chrom, initial=rd.dna.initial, final=rd.dna.final,
                                    name=rd.rna.str_rna(), orientation=rd.dna.orientation, 
                                    data=rd.score)
            else:
                dbs = GenomicRegion(chrom=rd.dna.chrom, initial=rd.dna.initial, final=rd.dna.final,
                                    name=rd.dna.name, orientation=rd.dna.orientation, 
                                    data=rd.score)

            if not orientation:
                dna_set.add(dbs)
            else:
                if orientation == rd.orient:
                    dna_set.add(dbs)
                else: pass
        if sort: dna_set.sort()
        if rm_duplicate: dna_set.remove_duplicates()
        return dna_set

    def sort_rbs(self):
        """Sort the dictionary by RNA"""
        self.sequences = sorted(self.sequences, key=lambda x: x.rna, cmp=GenomicRegion.__cmp__)
        self.sorted_rna = True
        
    def sort_dbs(self):
        """Sort the dictionary by DNA"""
        self.sequences = sorted(self.sequences, key=lambda x: x.dna, cmp=GenomicRegion.__cmp__)
        self.sorted_dna = True
    
    def sort_dbs_by_regions(self, regionset):
        """Sort the DBS by given GenomicRegionSet"""
        dbss = self.get_dbs(sort=True)

        result = {}

        if not regionset.sorted: regionset.sort()
        
        iter_dbs = iter(dbss)
        dbs = iter_dbs.next()

        last_j = len(regionset)-1
        j = 0
        cont_loop = True
        pre_inter = 0
        cont_overlap = False

        while cont_loop:
            # When the regions overlap
            if dbs.overlap(regionset[j]):
                result[regionset[j].toString()].add(dbs)

                if cont_overlap == False: pre_inter = j
                if j == last_j: 
                    try: dbs = iter_dbs.next()
                    except: cont_loop = False 
                else: 
                    j = j + 1
                    result[regionset[j].toString()] = GenomicRegionSet("RBS_"+regionset[j].toString())
                cont_overlap = True
            
            elif dbs < regionset[j]:
                try: 
                    dbs = iter_dbs.next()
                    j = pre_inter
                    cont_overlap = False
                except: cont_loop = False 
            
            elif dbs > regionset[j]:
                if j == last_j:
                    cont_loop = False
                else:
                    j = j + 1
                    result[regionset[j].toString()] = GenomicRegionSet("RBS_"+regionset[j].toString())
                    cont_overlap = False
        return result

    def sort_rbs_by_regions(self, regionset, merge_rbs=True):
        """Sort the DBS by given GenomicRegionSet"""
        
        result = OrderedDict()

        #txp_copy = copy.deepcopy(self)
        #txp_copy.sort_dbs()
        self.sort_dbs()

        if not regionset.sorted: regionset.sort()

        for region in regionset:
            result[region.toString()] = BindingSiteSet("binding site:"+region.toString())
        
        if len(self) == 0: 
            return result
        #iter_rd = iter(txp_copy)
        iter_rd = iter(self)
        rd = iter_rd.next()

        last_j = len(regionset)-1
        j = 0
        cont_loop = True
        pre_inter = 0
        cont_overlap = False

        while cont_loop:
            # When the regions overlap
            if rd.dna.overlap(regionset[j]):
                result[regionset[j].toString()].add(rd.rna)

                if cont_overlap == False: pre_inter = j
                if j == last_j: 
                    try: rd = iter_rd.next()
                    except: cont_loop = False 
                else: 
                    j = j + 1
                    
                cont_overlap = True
            
            elif rd.dna < regionset[j]:
                try: 
                    rd = iter_rd.next()
                    j = pre_inter
                    cont_overlap = False
                except: cont_loop = False 
            
            elif rd.dna > regionset[j]:
                if j == last_j:
                    cont_loop = False
                else:
                    j = j + 1
                    cont_overlap = False
         

        return result


    def sort_rd_by_regions(self, regionset):
        """Sort RNADNA binding information by a given GenomicRegionSet"""
        """Sort the DBS by given GenomicRegionSet"""

        result = OrderedDict()

        #txp_copy = copy.deepcopy(self)
        #txp_copy.sort_dbs()
        self.sort_dbs()

        if not regionset.sorted: regionset.sort()

        for region in regionset:
            result[region.toString()] = RNADNABindingSet("RNADNA_interaction:"+region.toString())

        if len(self) == 0: 
            return result
        #iter_rd = iter(txp_copy)
        iter_rd = iter(self)
        
        rd = iter_rd.next()

        last_j = len(regionset)-1
        j = 0
        cont_loop = True
        pre_inter = 0
        cont_overlap = False
        

        while cont_loop:
            # When the regions overlap
            if rd.dna.overlap(regionset[j]):
                result[regionset[j].toString()].add(rd)

                if cont_overlap == False: pre_inter = j
                if j == last_j: 
                    try: rd = iter_rd.next()
                    except: cont_loop = False 
                else: 
                    j = j + 1
                cont_overlap = True
            
            elif rd.dna < regionset[j]:
                try: 
                    rd = iter_rd.next()
                    j = pre_inter
                    cont_overlap = False
                except: cont_loop = False 
            
            elif rd.dna > regionset[j]:
                if j == last_j:
                    cont_loop = False
                else:
                    j = j + 1
                    cont_overlap = False
        return result

    def sort_rd_by_rbs(self, rbss):
        """Sort RNADNA binding information by a given BindingSiteSet"""
        
        result = OrderedDict()

        #txp_copy = copy.deepcopy(self)
        #txp_copy.sort_rbs()

        
        self.sort_rbs()

        if not rbss.sorted: rbss.sort()

        #iter_rd = iter(txp_copy)
        iter_rd = iter(self)
        rd = iter_rd.next()

        last_j = len(rbss)-1
        j = 0
        cont_loop = True
        pre_inter = 0
        cont_overlap = False

        while cont_loop:
            # When the regions overlap
            if rd.rna.overlap(rbss[j]):
                result[rbss[j].str_rna()].add(rd)

                if cont_overlap == False: pre_inter = j
                if j == last_j: 
                    try: rd = iter_rd.next()
                    except: cont_loop = False 
                else: 
                    j = j + 1
                    result[rbss[j].str_rna()] = RNADNABindingSet("RNADNA_interaction:"+rbss[j].str_rna())
                cont_overlap = True
            
            elif rd.dna < rbss[j]:
                try: 
                    rd = iter_rd.next()
                    j = pre_inter
                    cont_overlap = False
                except: cont_loop = False 
            
            elif rd.dna > rbss[j]:
                if j == last_j:
                    cont_loop = False
                else:
                    j = j + 1
                    result[rbss[j].str_rna()] = RNADNABindingSet("RNADNA_interaction:"+rbss[j].str_rna())
                    cont_overlap = False
        return result


    def merge_rbs(self, rbss=None, rm_duplicate=False, asgene_organism=None, region_set=None, 
                  name_replace=None, cutoff=0):
        """Merge the RNA binding regions which have overlap to each other and 
           combine their corresponding DNA binding regions.
        
        extend -> Define the extending length in basepair of each RNA binding regions
        perfect_match -> Merge only the exactly same RNA binding regions
        """
        # A dict: RBS as key, and GenomicRegionSet as its value
        self.merged_dict = OrderedDict()
        # reg = copy.deepcopy(region_set)
        # res = copy.deepcopy(name_replace)
        # reg = region_set
        # res = name_replace

        if not rbss:
            # Merge RBS
            rna_merged = self.get_rbs()
            rna_merged.merge()
        else:
            rna_merged = rbss

        for r in rna_merged:
            self.merged_dict[r] = GenomicRegionSet(r.toString())

        rbsm = iter(rna_merged)
        try: r = rbsm.next()
        except: return
                
        if self.sorted_rna: pass
        else: self.sort_rbs()

        con = iter(self)
        try: rd = con.next()
        except: return

        con_loop = True
        while con_loop:
            #print(".", end="")
            if r.overlap(rd.rna):
                self.merged_dict[r].add(rd.dna)
                try: rd = con.next()
                except: 
                    try:
                        r = rbsm.next()
                        # self.merged_dict[r] = GenomicRegionSet(r.toString())
                    except:
                        if rm_duplicate: self.merged_dict[r].remove_duplicates()
                        con_loop = False
            elif rd.rna < r:
                try: rd = con.next()
                except: 
                    try:
                        r = rbsm.next()
                        # self.merged_dict[r] = GenomicRegionSet(r.toString())
                    except: 
                        if rm_duplicate: self.merged_dict[r].remove_duplicates()
                        con_loop = False
            elif rd.rna > r:
                if rm_duplicate: self.merged_dict[r].remove_duplicates()
                    #print(self.merged_dict[r].sequences[0].name)
                try:
                    r = rbsm.next()
                    # self.merged_dict[r] = GenomicRegionSet(r.toString())
                except: 
                    try: rd = con.next()
                    except: con_loop = False

        if region_set:
            for r in self.merged_dict:
                s = region_set.intersect(self.merged_dict[r],
                                         mode=OverlapType.ORIGINAL,
                                         rm_duplicates=rm_duplicate)
                self.merged_dict[r] = s

        if not region_set and rm_duplicate:
            for r in self.merged_dict:
                self.merged_dict[r].remove_duplicates()

        if cutoff:
            for r in self.merged_dict:
                if len(self.merged_dict[r]) < cutoff:
                    n = self.merged_dict.pop(r, None)
        if name_replace:
            for r in self.merged_dict:
                self.merged_dict[r].replace_region_name(regions=name_replace)
        #self.merged_dict = new_dict
        elif asgene_organism:
            for r in self.merged_dict:
                self.merged_dict[r] = self.merged_dict[r].gene_association(organism=asgene_organism)
                # try: self.merged_dict[r] = self.merged_dict[r].gene_association(organism=asgene_organism)
                # except: pass

    def read_tpx(self, filename, dna_fine_posi=False, shift=None, seq=False):
        """Read txp file to load all interactions. """
        with open(filename) as f:
            for line in f:
                # print(line)
                line = line.strip()
                # skip the comment line
                if line.startswith("#"):
                    continue
                l = line.split()
                # Load binding site
                if len(l) == 12:
                    l.insert(8, "_")
                if len(l) == 13:
                    if "\tchrM:" in line: continue # skip chromosome Mitocondria
                    # Format of each line in txp
                    #     [0] Sequence-ID
                    #     [1] TFO start
                    #     [2] TFO end
                    #     [3] Duplex-ID
                    #     [4] TTS start
                    #     [5] TTS end
                    #     [6] Score
                    #     [7] Error-rate
                    #     [8] Errors
                    #     [9] Motif
                    #     [10] Strand
                    #     [11] Orientation
                    #     [12] Guanine-rate

                    # RNA binding site
                    if not self.name: self.name = l[0]

                    if shift:
                        rna_start, rna_end = int(l[1])+shift, int(l[2])+shift
                    else:
                        rna_start, rna_end = int(l[1]), int(l[2])

                    if rna_start > rna_end: rna_start, rna_end =  rna_end, rna_start

                    rna = BindingSite(chrom=l[0], initial=rna_start, final=rna_end, score=l[6],
                                      errors_bp=l[8], motif=l[9], orientation=l[11],
                                      guanine_rate=l[12])
                    # DNA binding site
                    if ":" in l[3]:
                        rg = l[3].split(":")[1].split("-")
                        if dna_fine_posi:
                            dna_start = int(rg[0]) + int(l[4])
                            dna_end = int(rg[0]) + int(l[5])
                            dna = GenomicRegion(chrom=l[3].split(":")[0], initial=dna_start, final=dna_end,
                                                name=l[3], orientation=l[10],
                                                data=[l[6], l[9], l[11]])  # score, motif, orientation
                        else:
                            dna = GenomicRegion(chrom=l[3].split(":")[0], initial=int(rg[0]), final=int(rg[1]),
                                                name=l[3], orientation=l[10],
                                                data=[l[6], l[9], l[11]])
                    else:
                        dna = GenomicRegion(chrom=l[3], initial=int(l[4]), final=int(l[5]),
                                            name=l[3], orientation=l[10],
                                            data=[l[6], l[9], l[11]])
                    if seq:
                        binding = []
                        rd = [l[6],l[7],l[8],l[12]]
                    else:
                        self.add(RNADNABinding(rna=rna, dna=dna, score=l[6], err_rate=l[7], err=l[8], guan_rate=l[12]))
                elif 0 < len(l) < 10 and seq:
                    if "TFO: " in line:
                        b = line.replace("TFO: ", "")
                    elif "TTS: " in line:
                        b = line.replace("TTS: ", "")
                    elif "|" in line:
                        b = "    " + line
                    else:
                        b = line
                    binding.append(b)
                elif len(l) == 0 and seq:

                    self.add(RNADNABinding(rna=rna, dna=dna, score=rd[0], err_rate=rd[1], err=rd[2],
                                           guan_rate=rd[3], match=binding))


    def map_promoter_name(self, promoters):
        """Give each DNA region the corresponding name from the given promoter"""
        self.sort_dbs()
        if not promoters.sorted: promoters.sort()
        con_p = iter(promoters)
        p = con_p.next()

        for rd in self:
            if p.overlap(rd.dna):
                rd.dna.name = p.name
            elif p < rd.dna:
                try: p = con_p.next()
                except: break

    def write_txp(self, filename, match=False):
        """Write RNADNABindingSet into the file"""
        d = os.path.dirname(filename)
        if d:
            try: os.stat(d)
            except: os.mkdir(d)
        with open(filename, "w") as f:
            print("# RNA-ID\tRBS-start\tRBS-end\tDNA-ID\tDBS-start\tDBS-end\tScore\tError-rate\tErrors\tMotif\tStrand\tOrientation\tGuanine-rate", file=f)
            for rd in self:
                print(str(rd), file=f)
                if match:
                    for i in range(4):
                        print(rd.match[i], file=f)
                print("", file=f)

    def write_bed(self, filename, dbd_tag=True, remove_duplicates=False, convert_dict=None, associated=False):
        """Write BED file for all the DNA Binding sites
        filename: define the output filename
        remove_duplicates: remove all exact duplicates
        convert_dict: given a dictionary to change the region name
        """
        dbss = self.get_dbs(dbd_tag=dbd_tag)
        if remove_duplicates:
            dbss.remove_duplicates()
        if convert_dict:
            dbss = dbss.change_name_by_dict(convert_dict=convert_dict)
        if associated:
            dbss.add_associated_gene_data(organism=associated)
        dbss.write(filename)

    def get_overlapping_regions(self, regionset):
        """Return a GenomicRegionSet which overlapping the given regions"""
        dbss = self.get_dbs(dbd_tag=True, sort=True)
        overlaps = dbss.intersect(regionset, mode=OverlapType.ORIGINAL)
        return overlaps

    def count_dbs_on_dbd(self, rna_len):
        """Return the frequency list of the appearance of DBS on RNA profile"""
        count = [0]*rna_len

        for rd in self.sequences:
            for i in range(rd.rna.initial, rd.dna.final):
                count[i] += 1

        return count

    def remove_duplicates(self):
        """Remove the duplicate RNA-DNA interactions and remain the unique ones. (No return)"""
        if self.sorted_dna == False: self.sort_dbs()
        
        for i in range(len(self.sequences) - 1):
            loop = True
            while loop:
                try:
                    if self.sequences[i] == self.sequences[i+1]:
                        del self.sequences[i+1]
                        loop = True
                    else:
                        loop = False
                except:
                    loop = False
                    continue

    def remove_duplicates_by_dbs(self):
        """Remove the duplicate RNA-DNA interactions and remain the unique ones. (No return)"""
        if self.sorted_dna == False: self.sort_dbs()
        
        for i in range(len(self.sequences) - 1):
            loop = True
            while loop:
                try:
                    if self.sequences[i].dna.toString() == self.sequences[i+1].dna.toString():
                        del self.sequences[i+1]
                        loop = True
                    else:
                        loop = False
                except:
                    loop = False
                    continue

    def generate_jaspar_matrix(self, genome_path, rbss=None):
        
        nucl = {"A":0 , "C":1 , "G":2 , "T":3 ,
                "a":0 , "c":1 , "g":2 , "t":3 }
        genome = pysam.Fastafile(genome_path)
        # A dict: RBS as key, and matrix as its value
        self.pwm_dict = OrderedDict()
        #print("merge_rbs")

        if not rbss:
            # Merge RBS
            rna_merged = self.get_rbs()
            rna_merged.merge()
        else:
            rna_merged = rbss

        for r in rna_merged:
            # Two same size matrix: P & A
            self.pwm_dict[r] = [ numpy.zeros((4, len(r))), numpy.zeros((4, len(r))) ]

        rbsm = iter(rna_merged)
        try: r = rbsm.next()
        except: return
                
        if self.sorted_rna: pass
        else: self.sort_rbs()

        con = iter(self)
        try: rd = con.next()
        except: return

        con_loop = True
        while con_loop:
            #print(".", end="")
            if r.overlap(rd.rna):
                seq = genome.fetch(rd.dna.chrom, max(0, rd.dna.initial), rd.dna.final)
                if rd.orient == "A": 
                    pa = 1
                    seq = seq[::-1]
                else: 
                    pa = 0
                print(seq)
                ind = rd.rna.initial - r.initial
                for i, a in enumerate(seq):
                    self.pwm_dict[r][pa][ nucl[a], ind + i ] += 1

                try: rd = con.next()
                except: 
                    try:
                        r = rbsm.next()
                    except:
                        con_loop = False
            elif rd.rna < r:
                try: rd = con.next()
                except: 
                    try:
                        r = rbsm.next()
                    except: 
                        con_loop = False
            elif rd.rna > r:
                try:
                    r = rbsm.next()
                except: 
                    try: rd = con.next()
                    except: con_loop = False

    def overlap_rbss(self, rbss):
        z = RNADNABindingSet(self.name)
        for rd in self.sequences:
            for rbs in rbss:
                if rd.rna.overlap(rbs):
                    z.add(rd)
        return z

    def motif_statistics(self, cm=False):
        self.motifs = {"Mix_Antiparallel": {"G": 0, "T": 0},
                       "Mix_Parallel": {"G": 0, "T": 0},
                       "Purine_Antiparallel": {"A": 0, "G": 0},
                       "Pyrimidine_Parallel": {"C": 0, "T": 0}}
        if cm:
            res_c = {"Mix_Purine_AA": {"G": 0, "T": 0, "A":0},
                     "Purine_Mix_AA": {"G": 0, "T": 0, "A": 0},
                     "Pyrimidine_Mix_PP": {"G": 0, "T": 0, "C": 0},
                     "Mix_Pyrimidine_PP": {"G": 0, "T": 0, "C": 0},
                     "Mix_Pyrimidine_PA": {"G": 0, "T": 0, "C": 0},
                     "Pyrimidine_Mix_PA": {"G": 0, "T": 0, "C": 0},
                     "Purine_Mix_PA": {"A": 0, "G": 0, "T": 0},
                     "Mix_Purine_PA": {"A": 0, "G": 0, "T": 0},
                     "Purine_Pyrimidine_PA": {"A": 0, "G": 0, "C": 0, "T": 0},
                     "Pyrimidine_Purine_PA": {"A": 0, "G": 0, "C": 0, "T": 0},
                     "Mix_Mix_PA": {"G": 0, "T": 0}
                     }
            self.motifs.update(res_c)
        if len(self) > 0:
            for s in self:
                if cm:
                    m = s.motif_statistics(cm=True)
                else:
                    m = s.motif_statistics()
                for mode in m:
                    for com in m[mode]:
                        self.motifs[mode][com] += m[mode][com]
        # print(self.motifs)

    def uniq_motif_statistics(self, rnalen, cm=False):
        self.uniq_motifs = {"MA_G": [0] * rnalen,
                            "MA_T": [0] * rnalen,
                            "MP_G": [0] * rnalen,
                            "MP_T": [0] * rnalen,
                            "RA_A": [0] * rnalen,
                            "RA_G": [0] * rnalen,
                            "YP_C": [0] * rnalen,
                            "YP_T": [0] * rnalen}
        if cm:
            res_c = {
                "M_RA_A_T": [0] * rnalen,
                "M_RA_A_A": [0] * rnalen,
                "M_RA_A_G": [0] * rnalen,
                "R_MA_A_T": [0] * rnalen,
                "R_MA_A_A": [0] * rnalen,
                "R_MA_A_G": [0] * rnalen,
                "M_RP_A_T": [0] * rnalen,
                "M_RP_A_A": [0] * rnalen,
                "M_RP_A_G": [0] * rnalen,
                "R_MA_P_T": [0] * rnalen,
                "R_MA_P_A": [0] * rnalen,
                "R_MA_P_G": [0] * rnalen,
                "M_YP_P_T": [0] * rnalen,
                "M_YP_P_C": [0] * rnalen,
                "M_YP_P_G": [0] * rnalen,
                "Y_MP_P_T": [0] * rnalen,
                "Y_MP_P_C": [0] * rnalen,
                "Y_MP_P_G": [0] * rnalen,
                "M_YA_P_T": [0] * rnalen,
                "M_YA_P_C": [0] * rnalen,
                "M_YA_P_G": [0] * rnalen,
                "Y_MP_A_T": [0] * rnalen,
                "Y_MP_A_C": [0] * rnalen,
                "Y_MP_A_G": [0] * rnalen,
                "R_YA_P_T": [0] * rnalen,
                "R_YA_P_C": [0] * rnalen,
                "R_YA_P_G": [0] * rnalen,
                "R_YA_P_A": [0] * rnalen,
                "Y_RP_A_C": [0] * rnalen,
                "Y_RP_A_G": [0] * rnalen,
                "Y_RP_A_T": [0] * rnalen,
                "Y_RP_A_A": [0] * rnalen,
                "M_MP_P_T": [0] * rnalen,
                "M_MP_P_G": [0] * rnalen,
                "M_MP_A_G": [0] * rnalen,
                "M_MP_A_T": [0] * rnalen,
                "M_MA_P_G": [0] * rnalen,
                "M_MA_P_T": [0] * rnalen,
                "M_MA_A_G": [0] * rnalen,
                "M_MA_A_T": [0] * rnalen,
            }
            self.uniq_motifs.update(res_c)

        if len(self) > 0:
            for rd in self:
                if rd.strand == "+":
                    rna = rd.match[0].split()[1]
                    linking = rd.match[1].strip()
                elif rd.strand == "-":
                    rna = rd.match[3].split()[1]
                    linking = rd.match[2].strip()
                if cm:
                    if rd.strand == "+_-":
                        rna = rd.match[0].split()[1] + rd.match[3].split()[3]
                        linking = rd.match[2].strip().split()[3] + rd.match[1].strip().split()[0]
                    elif rd.strand == "-_+":
                        rna = rd.match[3].split()[1] + rd.match[0].split()[3]
                        linking = rd.match[1].strip().split()[3] + rd.match[2].strip().split()[0]
                    elif rd.strand == "+_+":
                        rna = rd.match[0].split()[1] + rd.match[0].split()[3]
                        linking = rd.match[1].strip().split()[0] + rd.match[1].strip().split()[1]
                    elif rd.strand == "-_-":
                        rna = rd.match[3].split()[1] + rd.match[3].split()[3]
                        linking = rd.match[2].strip().split()[0] + rd.match[2].strip().split()[1]
                for i, l in enumerate(linking):
                    if l == "|":
                        self.uniq_motifs[rd.motif+rd.orient+"_"+rna[i].upper()][rd.rna.initial+i] += 1

    def rna_track(self, rnalen):
        self.rna_track = [0] * rnalen
        for rd in self:
            for i in range(rd.rna.initial, rd.rna.final):
                self.rna_track[i] += 1

    def get_RNA_DNA_counts(self, DNA_regions, filename):
        """Return a list of the number of binding sites between RNA (row) and DNAs (column)."""

        dbss = self.get_dbs()
        cov = DNA_regions.counts_per_region(regionset=dbss)

        with open(filename, "w") as f:
            print("\t".join([x.toString(underline=True) for x in DNA_regions]), file=f)
            print("\t".join([str(x) for x in cov]), file=f)

    def distance_distribution(self):
        dis_count = {"in_trans": 0, "in_cis": 0, "local": 0}
        # chr_1_954955_955150__REV
        for rd in self:
            if rd.rna.chrom.startswith("chr_"):
                if rd.rna.chrom.split("_")[3] == "FWD":
                    strand = "+"
                else:
                    strand = "-"
                r = GenomicRegion(chrom="chr" + rd.rna.chrom.split("_")[1],
                                  initial=int(rd.rna.chrom.split("_")[2])+rd.rna.initial,
                                  final=int(rd.rna.chrom.split("_")[3])+rd.rna.final,
                                  orientation=strand)
                d = r.distance(rd.dna)
                if not d: dis_count["in_trans"] += 1
                elif d > 10000: dis_count["in_cis"] += 1
                else: dis_count["local"] += 1
        return dis_count

    def region_distance(self, rbs1, rbs2):
        distance = 0
        if rbs1.chrom == rbs2.chrom:
            if rbs2.initial > rbs1.final:
                distance = rbs2.initial - rbs1.final
            elif rbs1.initial > rbs2.final:
                distance = rbs1.initial - rbs2.final
            return distance
        else:
            return None

    def combine_motif(self, name="rdb"):
        num_YM, num_MY, num_RM, num_MR = 0, 0, 0, 0
        s_PA_pos, s_PA_neg, s_AP_pos, s_AP_neg = 0, 0, 0, 0
        merged_tpx = RNADNABindingSet("combine_motif")
        for num_of_rbs in range(len(self.sequences) - 2):
            next_of_rbs = num_of_rbs + 1
            cur_rbs = self.sequences[num_of_rbs]
            next_rbs = self.sequences[next_of_rbs]
            while cur_rbs.rna.distance(next_rbs.rna) < 5 and next_of_rbs < len(self.sequences) - 2:
                # different RNA binding motif and same parallel
                if cur_rbs.rna.motif != next_rbs.rna.motif or cur_rbs.rna.orientation != next_rbs.rna.orientation:
                    # no RNA overlap
                    if self.region_distance(cur_rbs.dna, next_rbs.dna) > 0 and self.region_distance(cur_rbs.rna, next_rbs.rna) > 0:
                        # DNA gap smaller than 10bp and no overlap
                        if 0 < self.region_distance(cur_rbs.dna, next_rbs.dna) <= self.region_distance(cur_rbs.rna, next_rbs.rna): #5) and (cur_rbs.dna.chrom == next_rbs.dna.chrom):
                            # same strand on the DNA
                            if cur_rbs.dna.orientation == next_rbs.dna.orientation:
                                # combine motif - Parallel
                                if cur_rbs.rna.orientation == next_rbs.rna.orientation and next_rbs.rna.orientation == "P":
                                    # pos strand
                                    if cur_rbs.dna.orientation == "+":
                                        if cur_rbs.dna.final < next_rbs.dna.initial:
                                            merged_rbs = self.merge_combined_motif(cur_rbs, next_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)
                                            if cur_rbs.rna.motif == "Y":
                                                num_YM += 1
                                            else:
                                                num_MY += 1
                                    # neg strand
                                    if cur_rbs.dna.orientation == "-":
                                        if next_rbs.dna.final < cur_rbs.dna.initial:
                                            merged_rbs = self.merge_combined_motif(next_rbs, cur_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)
                                            if cur_rbs.rna.motif == "Y":
                                                num_YM += 1
                                            else:
                                                num_MY += 1

                                # combine motif - Anti-parallel
                                if cur_rbs.rna.orientation == next_rbs.rna.orientation and next_rbs.rna.orientation == "A":
                                    # pos strand
                                    if cur_rbs.dna.orientation == "+":
                                        if next_rbs.dna.final < cur_rbs.dna.initial:
                                            merged_rbs = self.merge_combined_motif(next_rbs, cur_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)
                                            if cur_rbs.rna.motif == "M":
                                                num_MR += 1
                                            else:
                                                num_RM += 1
                                    # neg strand
                                    if cur_rbs.dna.orientation == "-":
                                        if cur_rbs.dna.final < next_rbs.dna.initial:
                                            merged_rbs = self.merge_combined_motif(cur_rbs, next_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)
                                            if cur_rbs.rna.motif == "M":
                                                num_MR += 1
                                            else:
                                                num_RM += 1


                            # strand switch - different RNA orientation -different DNA strand
                            else:
                                if cur_rbs.rna.orientation != next_rbs.rna.orientation:
                                    # first rna ori == "P"
                                    if cur_rbs.rna.orientation == "P":
                                        if cur_rbs.dna.orientation == "+" and cur_rbs.dna.final < next_rbs.dna.initial:
                                            merged_rbs = self.merge_combined_motif(cur_rbs, next_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)
                                            s_PA_pos += 1
                                        elif cur_rbs.dna.orientation == "-" and next_rbs.dna.final < cur_rbs.dna.initial:
                                            s_PA_neg += 1
                                            merged_rbs = self.merge_combined_motif(next_rbs, cur_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)

                                    # first rna ori == "A"
                                    else:
                                        if cur_rbs.dna.orientation == "+" and next_rbs.dna.final < cur_rbs.dna.initial:
                                            s_AP_pos += 1
                                            merged_rbs = self.merge_combined_motif(next_rbs, cur_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)
                                        elif cur_rbs.dna.orientation == "-" and cur_rbs.dna.final < next_rbs.dna.initial:
                                            s_AP_neg += 1
                                            merged_rbs = self.merge_combined_motif(cur_rbs, next_rbs, rname=name)
                                            merged_tpx.add(merged_rbs)

                next_of_rbs += 1
                next_rbs = self.sequences[next_of_rbs]

        # print ("------  "+name + "\t\t[Done]  ---------------------------------")
        print([name, num_YM, num_MY, num_RM, num_MR], [s_PA_pos, s_PA_neg, s_AP_pos, s_AP_neg])
        # return [num_YM,num_MY,num_RM,num_MR,s_PA_pos,s_PA_neg,s_AP_pos,s_AP_neg],merged_tpx
        return merged_tpx

    def merge_match(self, rbs1, rbs2, dna_gap=0, rna_gap=0):
        if rbs1.strand == rbs2.strand and rbs1.strand == "+":
            match1 = rbs1.match
            match2 = rbs2.match
            match_a = ['', '', '', '']
            match_a[0] = match1[0][:len(match1[0]) - 4] + ' ' + str(rna_gap) + ' ' + match2[0][4:]
            match_a[1] = match1[1] + '   ' + match2[1][4:]
            match_a[2] = match1[2][:len(match1[2]) - 4] + ' ' + str(dna_gap) + ' ' + match2[2][4:]
            match_a[3] = match1[3][:len(match1[3]) - 4] + ' ' + str(dna_gap) + ' ' + match2[3][4:]
        elif rbs1.strand == rbs2.strand and rbs1.strand == "-":
            match1 = rbs1.match
            match2 = rbs2.match
            match_a = ['', '', '', '']
            match_a[0] = match1[0][:len(match1[0]) - 4] + ' ' + str(dna_gap) + ' ' + match2[0][4:]
            match_a[1] = match1[1][:len(match1[1]) - 4] + ' ' + str(dna_gap) + ' ' + match2[1][4:]
            match_a[2] = match1[2] + '   ' + match2[2][4:]
            match_a[3] = match1[3][:len(match1[3]) - 4] + ' ' + str(rna_gap) + ' ' + match2[3][4:]
        else:
            match1 = rbs1.match
            match2 = rbs2.match
            match_a = ['', '', '', '', '', '']
            if rbs1.rna.orientation == "+" and rbs2.rna.orientation == "-":
                match_a[0] = match1[0][:len(match1[0]) - 4]
                match_a[1] = match1[1]
                match_a[2] = match1[2][:len(match1[2]) - 4] + ' ' + str(dna_gap) + ' ' + match2[0][4:]
                match_a[3] = match1[3][:len(match1[3]) - 4] + ' ' + str(dna_gap) + ' ' + match2[1][4:]
                match_a[4] = ' '*len(match1[3]) + ' ' + match2[2][4:]
                match_a[5] = ' '*len(match1[3]) + ' ' + match2[3][4:] + 'r' + str(rna_gap)
                # match_a[0] = match1[0] + ' ' + match2[0] + "-r" + str(rna_gap)
                # match_a[1] = match1[1] + '     ' + match2[1]
                # match_a[2] = match1[2] + ' ' + match2[2]
                # match_a[3] = match1[3] + ' ' + match2[3] + "-d" + str(dna_gap)
            elif rbs1.rna.orientation == "-" and rbs2.rna.orientation == "+":
                match_a[0] = ' '*(len(match1[0]) - 4) + ' ' + match2[0][4:] + 'r'+str(rna_gap)
                match_a[1] = ' '*(len(match1[0]) - 4) + ' ' + match2[1][4:]
                match_a[2] = match1[0][:len(match1[0]) - 4] + ' ' + str(dna_gap) + ' ' + match2[2][4:]
                match_a[3] = match1[1][:len(match1[1]) - 4] + ' ' + str(dna_gap) + ' ' + match2[3][4:]
                match_a[4] = match2[2]
                match_a[5] = match2[3][:len(match1[3]) - 4]
                # match_a[0] = match1[0] + ' ' + match2[0] + "-r" + str(rna_gap)
                # match_a[1] = match1[1] + ' ' + match2[1]
                # match_a[2] = match1[2] + '     ' + match2[2]
                # match_a[3] = match1[3] + ' ' + match2[3] + "-d" + str(dna_gap)

        return match_a

    def merge_combined_motif(self, rbs1, rbs2, rname="rdb"):
        rna_score = str((int(rbs1.rna.score) + int(rbs2.rna.score)) / 2)
        err_rate = str((float(rbs1.err_rate) + float(rbs2.err_rate)) / 2)
        err = rbs1.err + "_" + rbs2.err
        guan_rate = str((float(rbs1.guan_rate) + float(rbs2.guan_rate)) / 2)
        rna_gap = self.region_distance(rbs1.rna, rbs2.rna)  # rbs2.rna.initial - rbs1.rna.final
        dna_gap = self.region_distance(rbs1.dna, rbs2.dna)
        match = self.merge_match(rbs1, rbs2, dna_gap, rna_gap)
        # dna_data =
        cmotif = rbs1.rna.motif + "_" + rbs2.rna.motif
        cr_orientation = rbs1.rna.orientation + "_" + rbs2.rna.orientation
        cd_orientation = rbs1.dna.orientation + "_" + rbs2.dna.orientation
        if rbs1.dna.final < rbs2.dna.initial:
            dna_ini = rbs1.dna.initial
            dna_fin = rbs2.dna.final
        else:
            dna_ini = rbs2.dna.initial
            dna_fin = rbs1.dna.final

        if rbs1.rna.final < rbs2.rna.initial:
            rna_ini = rbs1.rna.initial
            rna_fin = rbs2.rna.final
        else:
            rna_ini = rbs2.rna.initial
            rna_fin = rbs1.rna.final

        rna_m = BindingSite(rname, rna_ini, rna_fin, name=None, score=rna_score,
                        errors_bp=None, motif=cmotif, orientation=cr_orientation,
                        guanine_rate=None, seq=None)
        dna_m = GenomicRegion(rbs1.dna.chrom, dna_ini, dna_fin, name="Merge DNA region",
                          orientation=cd_orientation, data=rbs1.dna.data, proximity=None)

        rbs_m = RNADNABinding(rna_m, dna_m, rna_score, err_rate, err, guan_rate, match=match, rna_gap = rna_gap, dna_gap = dna_gap)
        return rbs_m
