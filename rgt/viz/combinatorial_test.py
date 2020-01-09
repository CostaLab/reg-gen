import sys
import copy
import itertools
from collections import OrderedDict
from ..Util import OverlapType
from ..GenomicRegionSet import GenomicRegionSet
from ..ExperimentalMatrix import ExperimentalMatrix


###########################################################################################
#                    Combinatorial test
###########################################################################################

def posi2region(regions, p):
    all = list(range(len(regions)))
    new_r = GenomicRegionSet(name="")
    for r in p:
        new_r.combine(regions[r])
    return new_r


class Combinatorial:
    def __init__(self, em_path, query):
        self.EM = ExperimentalMatrix()
        self.EM.read(em_path)
        self.regions = self.EM.get_regionsets()
        self.regionnames = self.EM.get_regionsnames()

    def group_refque(self, groupby):
        self.grouped_regions = OrderedDict()  # Store all bed names according to their types
        if groupby:
            for rs in self.regions:
                ty = self.EM.get_type(rs.name, groupby)
                try:
                    self.grouped_regions[ty].append(rs)
                except:
                    self.grouped_regions[ty] = [rs]
        else:
            print("** Please define grouping column '-g'")
            sys.exit(1)

    def posi2set(self, regions, p):
        all_len = list(range(len(regions)))
        inter_r = copy.deepcopy(regions[p[0]])

        for i in all_len:
            # print("inter_r: "+inter_r.name)
            if i in p[1:]:
                inter_r = inter_r.intersect(regions[i], mode=OverlapType.OVERLAP)
            elif i == p[0]:
                pass
            else:
                inter_r = inter_r.subtract(regions[i], whole_region=False)
        # print("inter_r: "+inter_r.name)
        return inter_r

    def combinatorial(self, background=None):
        def p2sign(plist, length):
            output = ["-"] * length
            for p in plist:
                output[p] = "+"
            return output

        new_refsp = OrderedDict()
        new_refs = OrderedDict()
        ref_names = []
        self.comb_ref_infor = {}

        for ty in list(self.groupedreference.keys()):
            n = len(self.groupedreference[ty])
            new_refs[ty] = []
            new_refsp[ty] = []

            for i in range(1, n):
                new_refsp[ty].append(itertools.combinations(list(range(n)), i))
            for posi in new_refsp[ty]:
                posi = [list(i) for i in posi]

                for p in posi:
                    # print("   " + str(p))
                    pr = self.posi2set(self.groupedreference[ty], p)
                    new_refs[ty].append(pr)
                    ref_names.append(pr.name)
                    self.comb_ref_infor[pr.name] = p2sign(p, n)
            all_int = self.posi2set(self.groupedreference[ty], list(range(n)))
            new_refs[ty].append(all_int)
            ref_names.append(all_int.name)
            self.comb_ref_infor[all_int.name] = p2sign(list(range(n)), n)
            """
            # Background
            unions = GenomicRegionSet(name="")
            for r in self.groupedreference[ty]:
                unions.combine(r)
            unions.name = " + ".join([r.name for r in self.groupedreference[ty]])

            nonset = self.backgroung.subtract(unions)
            nonset.name = "!("+"+".join([r.name for r in self.groupedreference[ty]]) + ")"
            new_refs[ty].append(nonset)
            ref_names.append(nonset.name)
            """
        # self.comb_reference = new_refs
        self.groupedreference = copy.deepcopy(new_refs)
        self.orig_refs = copy.deepcopy(self.referencenames)
        self.referencenames = list(set(ref_names))

    def combine_regions(self, background=None):
        new_refsp = OrderedDict()
        new_refs = OrderedDict()
        ref_names = []

        for ty in list(self.groupedreference.keys()):
            n = len(self.groupedreference[ty])
            new_refs[ty] = []
            new_refsp[ty] = []
            for i in range(1, n):
                new_refsp[ty].append(itertools.combinations(list(range(n)), i))
            for posi in new_refsp[ty]:
                posi = [list(i) for i in posi]
                for p in posi:
                    print(("   " + str(p)))
                    pr = posi2region(self.groupedreference[ty], p)
                    new_refs[ty].append(pr)
                    ref_names.append(pr.name)
