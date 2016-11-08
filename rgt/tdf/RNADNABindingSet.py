# Python Libraries
from __future__ import print_function
from collections import *
import numpy
import pysam
# Local Libraries

# Distal Libraries
from rgt.GenomicRegion import *
from rgt.GenomicRegionSet import *
from rgt.Util import OverlapType
from BindingSiteSet import BindingSite, BindingSiteSet

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
                 'rna_seq', 'dna_seq', 'match' ]

    def __init__(self, rna, dna, score, err_rate, err, guan_rate, rna_seq=None, dna_seq=None, match=None):
        """Initiation"""
        self.rna = rna  # BindingSite
        self.dna = dna  # GenomicRegion
        self.motif = rna.motif
        #self.strand = dna.orientation
        self.orient = rna.orientation
        self.score = score
        self.err_rate = err_rate
        self.err = err
        self.guan_rate = guan_rate
        self.rna_seq = rna_seq
        self.dna_seq = dna_seq
        self.match = match
    
    def __str__(self):
        return "\t".join( [ self.rna.chrom, str(self.rna.initial), str(self.rna.final), self.dna.toString(), "0", 
                            str(len(self.dna)), self.score, self.err_rate, self.err, self.motif,
                            self.strand, self.orient, self.guan_rate ] )

    def __eq__(self, other):
        return (self.dna, self.rna) == (other.dna, other.rna)
        
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
        res = name_replace

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
            for r in self.merged_dict.keys():
                s = region_set.intersect(self.merged_dict[r],
                                         mode=OverlapType.ORIGINAL,
                                         rm_duplicates=rm_duplicate)
                self.merged_dict[r] = s

        if not region_set and rm_duplicate:
            for r in self.merged_dict.keys():
                self.merged_dict[r].remove_duplicates()

        if cutoff:
            for r in self.merged_dict.keys():
                if len(self.merged_dict[r]) < cutoff:
                    n = self.merged_dict.pop(r, None)
        if res:
            for r in self.merged_dict.keys():
                self.merged_dict[r].replace_region_name(regions=res)
        #self.merged_dict = new_dict
        elif asgene_organism:
            for r in self.merged_dict.keys():
                try: self.merged_dict[r] = self.merged_dict[r].gene_association(organism=asgene_organism)
                except: pass

    def read_txp(self, filename, dna_fine_posi=False, shift=None):
        """Read txp file to load all interactions. """
        
        with open(filename) as f:
            for line in f:
                # skip the comment line
                if line.startswith("TFO:") or line.startswith("TTS:") or line.startswith("    ") or line.startswith("#"):
                    continue
                if "\tchrM:" in line: continue # skip chromosome Mitocondria
                line = line.strip("\n")
                line = line.split("\t")
                
                if len(line) < 10: 
                    #print(line)
                    continue # skip the unimportant lines in txp

                if len(line) == 12:
                    line.insert(8,"_")
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
                if not self.name: self.name = line[0]
                
                if shift:
                    rna_start, rna_end = int(line[1])+shift, int(line[2])+shift
                else:
                    rna_start, rna_end = int(line[1]), int(line[2])
                
                if rna_start > rna_end: rna_start, rna_end =  rna_end, rna_start

                rna = BindingSite(chrom=line[0], initial=rna_start, final=rna_end, score=line[6], 
                                  errors_bp=line[8], motif=line[9], orientation=line[11], 
                                  guanine_rate=line[12])
                # DNA binding site
                rg = line[3].split(":")[1].split("-")
                try: rg.remove("")
                except: pass
                
                if dna_fine_posi:
                    dna_start = int(rg[0]) + int(line[4])
                    dna_end = int(rg[0]) + int(line[5])
                    dna = GenomicRegion(chrom=line[3].split(":")[0], initial=dna_start, final=dna_end, 
                                        name=line[3], orientation=line[10],
                                        data=[line[6], line[9], line[11]]) # score, motif, orientation

                else:
                    try:
                        dna = GenomicRegion(chrom=line[3].split(":")[0], initial=int(rg[0]), final=int(rg[1]), 
                                            name=line[3], orientation=line[10],
                                            data=[line[6], line[9], line[11]])
                    except:
                        print(line)
                # Map RNA binding site to DNA binding site
                self.add(RNADNABinding(rna=rna, dna=dna, score=line[6], err_rate=line[7], err=line[8], guan_rate=line[12]))

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

    def write_txp(self, filename):
        """Write RNADNABindingSet into the file"""
        try: os.stat(os.path.dirname(filename))
        except: os.mkdir(os.path.dirname(filename))   
        with open(filename, "w") as f:
            print("# RNA-ID\tRBS-start\tRBS-end\tDNA-ID\tDBS-start\tDBS-end\tScore\tError-rate\tErrors\
                   \tMotif\tStrand\tOrientation\tGuanine-rate", file=f)
            for rd in self:
                print(str(rd), file=f) 

    def write_bed(self, filename, remove_duplicates=False, convert_dict=None, associated=False):
        """Write BED file for all the DNA Binding sites
        filename: define the output filename
        remove_duplicates: remove all exact duplicates
        convert_dict: given a dictionary to change the region name
        """
        dbss = self.get_dbs(dbd_tag=True)
        if remove_duplicates:
            dbss.remove_duplicates()
        if convert_dict:
            dbss = dbss.change_name_by_dict(convert_dict=convert_dict)
        if associated:
            dbss.add_associated_gene_data(organism=associated)
        dbss.write_bed(filename)


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
