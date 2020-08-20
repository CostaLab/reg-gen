"""
GenomicRegionSet
===================
GenomicRegionSet represent a list of GenomicRegions.

"""

###############################################################################
# Libraries
###############################################################################

# Python


import sys
import random
from ctypes import *
from scipy import stats
from copy import deepcopy
from collections import OrderedDict
import numpy

# Internal
from .SequenceSet import *
from .GeneSet import GeneSet
from .GenomicRegion import GenomicRegion
from .Util import GenomeData, OverlapType, LibraryPath

# External
from pysam import Fastafile

random.seed(42)


###############################################################################
# Class
###############################################################################


class GRSFileIO:
    class Bed:
        """
        Each row maps to a GenomicRegion.

        Note: Chrom (1), start (2), end (2), name (4) and orientation (6) is used for GenomicRegion.
              All other columns (5, 7, 8, ...) are put to the data attribute of the GenomicRegion.
              The numbers in parentheses are the columns of the BED format.
        """

        @staticmethod
        def read_to_grs(grs, filename):
            with open(filename) as f:
                error_line = 0  # Count error line
                for line in f:
                    line = line.strip("\n")
                    line = line.split()
                    try:
                        name, orientation, data = None, None, None
                        size = len(line)
                        chrom = line[0]
                        start, end = int(line[1]), int(line[2])

                        if start > end:
                            start, end = end, start
                        if size > 3:
                            name = line[3]

                        if size > 5:
                            orientation = line[5]
                            data = "\t".join([line[4]] + line[6:])
                        if size == 5:
                            data = line[4]

                        if start == end:
                            raise Exception(
                                "zero-length region: " + grs.chrom + "," + str(grs.initial) + "," + str(grs.final))
                        g = GenomicRegion(chrom, start, end, name, orientation, data)

                        grs.add(g)
                    except:
                        if not line:
                            continue
                        else:
                            error_line += 1
                            # if error_line > 2:
                            #     # Skip the first error line which contains the track information
                            #     print("Error at line", line, filename)
                grs.sort()

            return grs

        @staticmethod
        def write_from_grs(grs, filename, mode="w"):
            with open(filename, mode) as f:
                for gr in grs:
                    print(gr, file=f)

    class Bed12:
        """
        Bed file with "block information", eg exons.
        """

        @staticmethod
        def read_to_grs(grs, filename):
            with open(filename) as f:
                error_line = 0  # Count error line
                for line in f:
                    line = line.strip("\n")
                    line = line.split()
                    try:
                        name, orientation, data = None, None, None
                        size = len(line)
                        chrom = line[0]
                        start, end = int(line[1]), int(line[2])

                        if start > end:
                            start, end = end, start
                        if size > 3:
                            name = line[3]

                        if size > 5:
                            orientation = line[5]
                            data = "\t".join([line[4]] + line[6:])
                        if size == 5:
                            data = line[4]

                        if start == end:
                            raise Exception(
                                "zero-length region: " + grs.chrom + "," + str(grs.initial) + "," + str(grs.final))
                        g = GenomicRegion(chrom, start, end, name, orientation, data)

                        if size == 12 and int(line[6]) and int(line[7]) and int(line[9]):
                            gs = g.extract_blocks()
                            for gg in gs:
                                grs.add(gg)
                        else:
                            grs.add(g)
                    except:
                        if not line:
                            continue
                        else:
                            error_line += 1
                            if error_line > 2:
                                # Skip the first error line which contains the track information
                                print("Error at line", line, filename)
                grs.sort()

            return grs

        @staticmethod
        def write_from_grs(grs, filename, mode="w"):
            with open(filename, mode) as f:

                blocks = {}
                for gr in grs:
                    try:
                        blocks[gr.name].add(gr)
                    except:
                        blocks[gr.name] = GenomicRegionSet(gr.name)
                        blocks[gr.name].add(gr)

                for name in list(blocks.keys()):
                    blocks[name].merge()
                    start = min([g.initial for g in blocks[name]])
                    end = max([g.final for g in blocks[name]])

                    block_width = []
                    block_start = []
                    for g in blocks[name]:
                        block_width.append(len(g))
                        block_start.append(g.initial - start)
                    block_count = len(blocks[name])

                    print("\t".join([blocks[name][0].chrom,
                                     str(start),
                                     str(end),
                                     name,
                                     "0",
                                     blocks[name][0].orientation,
                                     str(start),
                                     str(end),
                                     "0",
                                     str(block_count),
                                     ",".join([str(w) for w in block_width]),
                                     ",".join([str(s) for s in block_start])]), file=f)

    class BedGraph:
        """
        BedGraph format, read-only.
        """

        @staticmethod
        def read_to_grs(grs, filename):
            with open(filename) as f:
                for line in f:
                    try:
                        line = line.strip("\n")
                        line = line.split("\t")
                        assert len(line) == 4

                        chrom, start, end, data = line[0], int(line[1]), int(line[2]), str(line[3])

                        grs.add(GenomicRegion(chrom=chrom, initial=start, final=end, data=data))
                    except:
                        print("Error at line", line, filename)

                grs.sort()

            return grs

        @staticmethod
        def write_from_grs(grs, filename, mode="w"):
            raise NotImplementedError

    class Fasta:
        # FIXME: this one is a bit strange. It's based on the assumption that the GRS has already been
        # populated with GenomicRegions, and then it adds the corresponding DNA sequence to each of those.

        @staticmethod
        def read_to_grs(grs, filename):
            fasta = Fastafile(filename)

            for r in grs:
                try:
                    seq = str(fasta.fetch(r.chrom, r.initial, r.final))
                    strand = r.orientation if r.orientation else "+"
                    r.sequence = Sequence(seq=seq, name=str(r), strand=strand)
                except:
                    print("Warning: region '%s' is skipped" % r)

        @staticmethod
        def write_from_grs(grs, filename, mode="w"):
            with open(filename, mode) as f:
                for r in grs:
                    try:
                        f.write(">" + r.chrom + ":" + str(r.initial) + "-" + str(r.final) + "-" +
                                str(r.orientation) + "\n" + r.sequence.seq + "\n")
                    except:
                        print("Warning: region '%s' is skipped" % r)


class GenomicRegionSet:
    """*Keyword arguments:*

        - name -- Name of the GenomicRegionSet
    """

    def __init__(self, name):
        self.name = name
        self.sequences = []
        self.sorted = False

    def read(self, filename, io=GRSFileIO.Bed):
        io.read_to_grs(self, filename)

    def write(self, filename, mode="w", io=GRSFileIO.Bed):
        io.write_from_grs(self, filename, mode)

    def get_chrom(self):
        """Return all chromosomes."""
        return [r.chrom for r in self]

    def get_names(self):
        """Return a list of all region names. If the name is None, it return the region string."""
        names = []
        for r in self:
            if r.name:
                names.append(r.name)
            else:
                names.append(r.toString())
        return names

    def add(self, region):
        """Add GenomicRegion.

        *Keyword arguments:*
            - region -- The GenomicRegion to be added.
        """
        self.sequences.append(region)
        self.sorted = False

    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        return iter(self.sequences)

    def __getitem__(self, key):
        return self.sequences[key]

    def extend(self, left, right, percentage=False, w_return=False):
        """Perform extend step for every element.

        *Keyword arguments:*

            - percentage -- input value of left and right can be any positive value or negative value larger than -50 %
        """
        z = GenomicRegionSet(name=self.name)

        if percentage:
            if percentage > -50:
                for s in self:
                    if w_return:
                        z.add(s.extend(int(len(s) * left / 100), int(len(s) * right / 100), w_return=True))
                    else:
                        s.extend(int(len(s) * left / 100), int(len(s) * right / 100))
            else:
                print("Percentage for extension must be larger than 50%%.")
                sys.exit(0)
        else:
            for s in self:
                if w_return:
                    z.add(s.extend(left, right, w_return=True))
                else:
                    s.extend(left, right)

        if w_return:
            return z
        else:
            return

    def extend_upstream(self, length=1000, w_return=False):
        """Perform extend step upstream for every element.

        *Keyword arguments:*

            - length -- Extending length
        """
        z = GenomicRegionSet(name=self.name)
        for s in self:
            if w_return:
                if s.orientation == "+":
                    z.add(s.extend(left=length, right=0, w_return=True))
                else:
                    z.add(s.extend(left=0, right=length, w_return=True))
            else:
                if s.orientation == "+":
                    s.extend(left=length, right=0)
                else:
                    s.extend(left=0, right=length)
        if w_return:
            return z
        else:
            return

    def extend_downstream(self, length=1000, w_return=False):
        """Perform extend step downstream for every element.

        *Keyword arguments:*

            - length -- Extending length
        """
        z = GenomicRegionSet(name=self.name)
        for s in self:
            if w_return:
                if s.orientation == "+":
                    z.add(s.extend(left=0, right=length, w_return=True))
                else:
                    z.add(s.extend(left=length, right=0, w_return=True))
            else:
                if s.orientation == "+":
                    s.extend(left=0, right=length)
                else:
                    s.extend(left=length, right=0)
        if w_return:
            return z
        else:
            return

    def sort(self, key=None, reverse=False):
        """Sort Elements by criteria defined by a GenomicRegion.

        *Keyword arguments:*

            - key -- given the key for comparison.
            - reverse -- reverse the sorting result.
        """
        if key:
            self.sequences.sort(key=key, reverse=reverse)
        else:
            self.sequences.sort()
            self.sorted = True

    def get_sequences(self, genome_fasta, ex=0):
        """Return SequenceSet according to the given regions and genome"""
        seq = SequenceSet(name=self.name, seq_type="DNA")
        seq.read_regions(regionset=self, genome_fasta=genome_fasta, ex=ex)
        return seq

    def random_subregions(self, size, name=None):
        """Return a subsampling of the genomic region set with a specific number of regions.

        *Keyword arguments:*

            - size -- define number of the subsampling regions.
        """
        if not name: name = self.name + '_random'
        z = GenomicRegionSet(name)
        samp = random.sample(list(range(len(self))), size)
        for i in samp:
            z.add(self.sequences[i])
        return z

    def random_split(self, size):
        """Return two exclusive GenomicRegionSets from self randomly.

        *Keyword arguments:*

            - size -- define number of the spliting regions.
        """
        a, b = GenomicRegionSet('random_split1'), GenomicRegionSet('random_split2')
        samp = random.sample(list(range(len(self))), size)
        for i in range(len(self)):
            if i in samp:
                a.add(self.sequences[i])
            else:
                b.add(self.sequences[i])
        return a, b

    def gene_association(self, organism, gene_set=None, promoter_length=1000,
                         thresh_dist=100000, show_dis=False, strand_specific=False):
        """Associates coordinates to genes given the following rules:

            1. If the peak is inside gene (promoter+coding) then this peak is associated with that gene.
            2. If a peak is inside overlapping genes, then the peak is annotated with both genes.
            3. If peak is between two genes (not overlapping neither), then both genes are annotated.
            4. If the distance between peak and gene is greater than a threshold distance, then it is not annotated.

        *Keyword arguments:*

            - gene_set -- List of gene names as a GeneSet object. If None, then consider all genes to be enriched. (default None)
            - organism -- Organism in order to fetch genomic data. (default hg19)
            - promoter_length -- Length of the promoter region. (default 1000)
            - thresh_dist -- Threshold maximum distance for a coordinate to be considered associated with a gene. (default 50000)
            - show_dis -- Show distance to the closest genes in parentheses.

        *Return:*

            - result_grs -- GenomicRegionSet exactly as self, but with the following additional information:

                1. name: String of genes associated with that coordinate separated by ':'
                2. data: String of proximity information (if the coordinate matched to the corresponding gene in the previous list in a proximal position (PROX) or distal position (DIST)) separated by ':'

                The gene will contain a '.' in the beginning of its name if it is not in the gene_set given.
        """

        z = GenomicRegionSet(self.name)
        if len(self) == 0:
            return self

        else:
            # If there is overlap within self or y, they should be merged first.
            if not self.sorted: self.sort()

            genome = GenomeData(organism)

            genes = GenomicRegionSet("genes")
            genes.read(genome.get_gene_regions())

            if gene_set:
                new_genes = GenomicRegionSet("genes")
                name_list = [g.upper() for g in gene_set.genes]
                for g in genes:
                    if g.name.upper() in name_list: new_genes.add(g)
                genes = new_genes
            genes.extend_upstream(length=promoter_length)

            if not genes.sorted: genes.sort()

            last_j = len(genes) - 1
            j = 0
            if strand_specific:
                pre_inter = [0, 0]
            else:
                pre_inter = 0

            for s in self:
                cont_loop = True
                cont_overlap = False
                asso_names = {"overlap": [], "close_l": [], "close_r": []}

                while cont_loop:
                    if strand_specific and s.orientation != genes[j].orientation:
                        if j == last_j:
                            cont_loop = False
                        else:
                            j += 1
                    else:
                        d = s.distance(genes[j])
                        if d == 0:
                            asso_names["overlap"].append(genes[j].name)
                            if not cont_overlap:
                                if strand_specific and s.orientation == "+":
                                    pre_inter[0] = j
                                elif strand_specific and s.orientation == "-":
                                    pre_inter[1] = j
                                elif not strand_specific:
                                    pre_inter = j
                            if j == last_j:
                                cont_loop = False
                            else:
                                j += 1
                                cont_overlap = True

                        elif asso_names["overlap"]:
                            if strand_specific:
                                if pre_inter[0] > 0 and pre_inter[1] > 0:
                                    j = min(pre_inter)
                                elif pre_inter[0] == 0 and pre_inter[1] > 0:
                                    j = pre_inter[1]
                                elif pre_inter[0] > 0 and pre_inter[1] == 0:
                                    j = pre_inter[0]
                            elif s.chrom == genes[j].chrom and pre_inter > 0:
                                j = pre_inter
                            cont_loop = False

                        elif d is not None and 0 < d < thresh_dist:
                            if show_dis:
                                if s > genes[j]:
                                    if asso_names["close_l"] and d < asso_names["close_l"][0]:
                                        asso_names["close_l"] = [d, genes[j].name + "(-" + str(d) + ")"]
                                    else:
                                        asso_names["close_l"] = [d, genes[j].name + "(-" + str(d) + ")"]

                                elif s < genes[j]:
                                    asso_names["close_r"] = [d, genes[j].name + "(+" + str(d) + ")"]
                                    cont_loop = False
                            else:
                                if s > genes[j]:

                                    if asso_names["close_l"] and d < asso_names["close_l"][0]:
                                        asso_names["close_l"] = [d, genes[j].name + "(-)"]
                                    else:
                                        asso_names["close_l"] = [d, genes[j].name + "(-)"]
                                elif s < genes[j]:
                                    asso_names["close_r"] = [d, genes[j].name + "(+)"]
                                    cont_loop = False
                            if j == last_j:
                                cont_loop = False
                            else:
                                j += 1

                        elif s < genes[j]:
                            if strand_specific and s.orientation == "+":
                                if s.chrom == genes[j].chrom and pre_inter[0] > 0:
                                    j = pre_inter[0]
                            elif strand_specific and s.orientation == "-":
                                if s.chrom == genes[j].chrom and pre_inter[1] > 0:
                                    j = pre_inter[1]
                            elif s.chrom == genes[j].chrom and pre_inter > 0:
                                j = pre_inter
                            cont_loop = False

                        elif s > genes[j]:
                            if j == last_j:
                                cont_loop = False
                            else:
                                j += 1

                if asso_names["overlap"]:
                    z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final,
                                        orientation=s.orientation, name=":".join(asso_names["overlap"]),
                                        data=s.data, proximity=0))
                elif asso_names["close_l"] and asso_names["close_r"]:
                    ss = [asso_names["close_l"][1],
                          asso_names["close_r"][1]]
                    dd = [str(asso_names["close_l"][0]),
                          str(asso_names["close_r"][0])]
                    z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final,
                                        orientation=s.orientation, name=":".join(ss),
                                        data=s.data, proximity=":".join(dd)))
                elif asso_names["close_l"] and not asso_names["close_r"]:
                    z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final,
                                        orientation=s.orientation, name=asso_names["close_l"][1],
                                        data=s.data, proximity=str(asso_names["close_l"][0])))
                elif not asso_names["close_l"] and asso_names["close_r"]:
                    z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final,
                                        orientation=s.orientation, name=asso_names["close_r"][1],
                                        data=s.data, proximity=str(asso_names["close_r"][0])))
                    # ss = []
                    # try: ss.append(asso_names["close_l"][1])
                    # except: pass
                    # try: ss.append(asso_names["close_r"][1])
                    # except: pass
                    # z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final,
                    #                     orientation=s.orientation, name=":".join(ss),
                    #                     data=s.data, proximity=s.proximity))
                else:
                    z.add(GenomicRegion(chrom=s.chrom, initial=s.initial, final=s.final,
                                        orientation=s.orientation, name=".",
                                        data=s.data, proximity=s.proximity))

            return z

    def filter_by_gene_association(self, gene_set=None, organism="hg19", promoter_length=1000, thresh_dist=50000):
        """Updates self in order to keep only the coordinates associated to genes which are in gene_set.
        
        It also returns information regarding the mapped genes and proximity information (if mapped in proximal (PROX)
        or distal (DIST) regions). If a region was associated with two genes, one in the gene_set and the other not in
        gene_set, then the updated self still keeps information about both genes. However, the latter won't be reported
        as mapped gene.

        *Keyword arguments:*

            - gene_set -- List of gene names as a GeneSet object. If None, then consider all genes to be enriched.
              (default None)
            - organism -- Organism in order to fetch genomic data. (default hg19)
            - promoter_length -- Length of the promoter region. (default 1000)
            - thresh_dist -- Threshold maximum distance for a coordinate to be considered associated with a gene.
              (default 50000)

       *Return:*

            - None -- Updates self in order to keep only the coordinates associated to genes which are in gene_set
            - all_genes = GeneSet that contains all genes associated with the coordinates
            - mapped_genes = GeneSet that contains the genes associated with the coordinates which are in gene_set
            - all_proxs = List that contains all the proximity information of genes associated with the coordinates
            - mapped_proxs = List that contains all the proximity information of genes associated with the coordinates
              which are in gene_set
        """

        # Making association with gene_set
        assoc_grs = self.gene_association(organism=organism, gene_set=gene_set, promoter_length=promoter_length,
                                          thresh_dist=thresh_dist)

        # Initializing updated genomic region set
        updated_grs = GenomicRegionSet(
            assoc_grs.name)  # New GRS will contain only coordinates associated to genes which are in gene_set

        # Initializing resulting gene sets
        all_genes = GeneSet("all_genes")  # Contains all genes associated with the coordinates
        mapped_genes = GeneSet(
            "mapped_genes")  # Contains the genes associated with the coordinates which are in gene_set
        all_proxs = []  # Contains all the proximity information of genes associated with the coordinates
        mapped_proxs = []  # Contains all the proximity information of genes associated with the coordinates which are in gene_set

        # Iterating on associated genomic regions
        for gr in assoc_grs:

            # If the coordinate wasn't associated with any gene, it is not considered
            if gr.name == ".": continue

            # Fetching list of genes and proximity information
            curr_genes = gr.name.split(":")
            curr_proxs = gr.proximity.split(":")

            # Updating all/mapped - genes/proxs
            flag_assoc = False
            for i in range(0, len(curr_genes)):
                g = curr_genes[i]
                p = curr_proxs[i]
                if g[0] != ".":
                    mapped_genes.genes.append(g)
                    mapped_proxs.append(p)
                    flag_assoc = True
                else:
                    g = g[1:]
                all_genes.genes.append(g)
                all_proxs.append(p)

            # Adding current genomic region to updated GRS if there was a gene associated with it in gene_set
            if flag_assoc: updated_grs.add(deepcopy(gr))

        # Removing duplicates
        all_remdup = list(set([all_genes.genes[i] + "_" + all_proxs[i] for i in range(0, len(all_genes.genes))]))
        mapped_remdup = list(
            set([mapped_genes.genes[i] + "_" + mapped_proxs[i] for i in range(0, len(mapped_genes.genes))]))
        all_genes = GeneSet("all_genes")
        mapped_genes = GeneSet("mapped_genes")
        all_proxs = []
        mapped_proxs = []
        for e in all_remdup:
            all_genes.genes.append(e.split("_")[0])
            all_proxs.append(e.split("_")[1])
        for e in mapped_remdup:
            mapped_genes.genes.append(e.split("_")[0])
            mapped_proxs.append(e.split("_")[1])

        # Updating self
        self.name = updated_grs.name
        self.sequences = updated_grs.sequences
        self.sorted = False

        return all_genes, mapped_genes, all_proxs, mapped_proxs

    def intersect(self, y, mode=OverlapType.OVERLAP, rm_duplicates=False):
        """Return the overlapping regions with three different modes.

        *Keyword arguments:*

            - y -- the GenomicRegionSet which to compare with.
            - mode -- OverlapType.OVERLAP, OverlapType.ORIGINAL or OverlapType.COMP_INCL.
            - rm_duplicates -- remove duplicates within the output GenomicRegionSet

        *Return:*
        
            - A GenomicRegionSet according to the given overlapping mode.

        *mode = OverlapType.OVERLAP*
        
            Return new GenomicRegionSet including only the overlapping regions with y.

            .. note:: it will merge the regions.
        
            ::

                self       ----------              ------
                y                 ----------                    ----
                Result            ---

        *mode = OverlapType.ORIGINAL*
        
            Return the regions of original GenomicRegionSet which have any intersections with y.

            ::

                self       ----------              ------
                y              ----------                    ----
                Result     ----------

        *mode = OverlapType.COMP_INCL*
        
            Return region(s) of the GenomicRegionSet which are 'completely' included by y.

            ::

                self        -------------             ------
                y              ----------      ---------------              ----
                Result                                ------
        """

        return self.intersect_c(y, mode, rm_duplicates)

    def intersect_python(self, y, mode=OverlapType.OVERLAP, rm_duplicates=False):
        z = GenomicRegionSet(self.name)
        if len(self) == 0 or len(y) == 0:
            return z

        else:
            # If there is overlap within self or y, they should be merged first.
            a = copy.deepcopy(self)
            b = copy.deepcopy(y)
            if not a.sorted: a.sort()
            if not b.sorted: b.sort()
            if mode == OverlapType.OVERLAP:
                a.merge()
                b.merge()

            iter_a = iter(a)
            s = next(iter_a)
            last_j = len(b) - 1
            j = 0
            cont_loop = True
            pre_inter = 0
            cont_overlap = False
            ####################### OverlapType.OVERLAP ###############################
            if mode == OverlapType.OVERLAP:
                while cont_loop:
                    # When the regions overlap
                    if s.overlap(b[j]):
                        z.add(GenomicRegion(chrom=s.chrom,
                                            initial=max(s.initial, b[j].initial),
                                            final=min(s.final, b[j].final),
                                            name=s.name,
                                            orientation=s.orientation,
                                            data=s.data,
                                            proximity=s.proximity))

                        if not cont_overlap:
                            pre_inter = j
                        if j == last_j:
                            try:
                                s = next(iter_a)
                            except:
                                cont_loop = False
                        else:
                            j += 1
                        cont_overlap = True

                    elif s < b[j]:
                        try:
                            s = next(iter_a)
                            if s.chrom == b[j].chrom and pre_inter > 0:
                                j = pre_inter
                            cont_overlap = False
                        except:
                            cont_loop = False

                    elif s > b[j]:
                        if j == last_j:
                            cont_loop = False
                        else:
                            j += 1
                            cont_overlap = False
                    else:
                        try:
                            s = next(iter_a)
                        except:
                            cont_loop = False

            ####################### OverlapType.ORIGINAL ###############################
            if mode == OverlapType.ORIGINAL:
                while cont_loop:
                    # When the regions overlap
                    if s.overlap(b[j]):
                        z.add(s)
                        try:
                            s = next(iter_a)
                        except:
                            cont_loop = False
                    elif s < b[j]:
                        try:
                            s = next(iter_a)
                        except:
                            cont_loop = False
                    elif s > b[j]:
                        if j == last_j:
                            cont_loop = False
                        else:
                            j += 1
                    else:
                        try:
                            s = next(iter_a)
                        except:
                            cont_loop = False
            ####################### OverlapType.COMP_INCL ###############################
            if mode == OverlapType.COMP_INCL:
                while cont_loop:
                    # When the regions overlap
                    if s.overlap(b[j]):
                        if s.initial >= b[j].initial and s.final <= b[j].final:
                            z.add(s)
                        if not cont_overlap: pre_inter = j
                        if j == last_j:
                            try:
                                s = next(iter_a)
                            except:
                                cont_loop = False
                        else:
                            j += 1
                        cont_overlap = True

                    elif s < b[j]:
                        try:
                            s = next(iter_a)
                            if s.chrom == b[j].chrom and pre_inter > 0:
                                j = pre_inter
                            cont_overlap = False
                        except:
                            cont_loop = False

                    elif s > b[j]:
                        if j == last_j:
                            cont_loop = False
                        else:
                            j += 1
                            cont_overlap = False
                    else:
                        try:
                            s = next(iter_a)
                        except:
                            cont_loop = False

            if rm_duplicates: z.remove_duplicates()
            # z.sort()
            return z

    def intersect_c(self, y, mode=OverlapType.OVERLAP, rm_duplicates=False):
        # Determine path of shared library
        Lib = LibraryPath()
        lib = cdll.LoadLibrary(Lib.get_c_rgt())
        # C-Binding of intersect overlap function
        intersect_overlap_c = lib.intersectGenomicRegionSetsOverlap
        intersect_overlap_c.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_char_p),
                                        POINTER(c_int), POINTER(c_int), c_int, POINTER(POINTER(c_int)),
                                        POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int)]
        intersect_overlap_c.restype = None

        # C-Binding of intersect original function
        intersect_original_c = lib.intersectGenomicRegionSetsOriginal
        intersect_original_c.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_char_p),
                                         POINTER(c_int), POINTER(c_int), c_int, POINTER(POINTER(c_int)),
                                         POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int)]
        intersect_original_c.restype = None

        # C-Binding of intersect completely function
        intersect_completely_included_c = lib.intersectGenomicRegionSetsCompletelyIncluded
        intersect_completely_included_c.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int,
                                                    POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int,
                                                    POINTER(POINTER(c_int)), POINTER(POINTER(c_int)),
                                                    POINTER(POINTER(c_int)), POINTER(c_int)]
        intersect_completely_included_c.restype = None

        # If one of the sets is empty, the intersection is trivially empty as well
        result = GenomicRegionSet(self.name)
        if len(self) == 0 or len(y) == 0:
            return result

        # Otherwise
        else:
            a = self
            b = y
            # Sort sets if necessary
            if not a.sorted:
                a.sort()
            if not b.sorted:
                b.sort()
            assert a.sorted
            assert b.sorted

            # If there is overlap within a or b, they should be merged first.
            if mode == OverlapType.OVERLAP:
                a = a.merge(w_return=True)
                b = b.merge(w_return=True)

            # Convert to ctypes
            len_self = len(a)
            len_y = len(b)
            max_len_result = len_self + len_y

            chromosomes_self_python = [gr.chrom.encode("utf-8") for gr in a.sequences]
            chromosomes_self_c = (c_char_p * len_self)(*chromosomes_self_python)

            chromosomes_y_python = [gr.chrom.encode("utf-8") for gr in b.sequences]
            chromosomes_y_c = (c_char_p * len_y)(*chromosomes_y_python)

            initials_self_python = [gr.initial for gr in a.sequences]
            initials_self_c = (c_int * len_self)(*initials_self_python)

            initials_y_python = [gr.initial for gr in b.sequences]
            initials_y_c = (c_int * len_y)(*initials_y_python)

            finals_self_python = [gr.final for gr in a.sequences]
            finals_self_c = (c_int * len_self)(*finals_self_python)

            finals_y_python = [gr.final for gr in b.sequences]
            finals_y_c = (c_int * len_y)(*finals_y_python)

            indices_c = POINTER(c_int)((c_int * max_len_result)())
            initials_result_c = POINTER(c_int)((c_int * max_len_result)())
            finals_result_c = POINTER(c_int)((c_int * max_len_result)())
            size_result_c = c_int()

            # Call C-function
            if mode == OverlapType.OVERLAP:
                intersect_overlap_c(chromosomes_self_c, initials_self_c, finals_self_c, len_self, chromosomes_y_c,
                                    initials_y_c, finals_y_c, len_y, pointer(indices_c), pointer(initials_result_c),
                                    pointer(finals_result_c), byref(size_result_c))
            elif mode == OverlapType.ORIGINAL:
                intersect_original_c(chromosomes_self_c, initials_self_c, finals_self_c, len_self, chromosomes_y_c,
                                     initials_y_c, finals_y_c, len_y, pointer(indices_c), pointer(initials_result_c),
                                     pointer(finals_result_c), byref(size_result_c))
            elif mode == OverlapType.COMP_INCL:
                intersect_completely_included_c(chromosomes_self_c, initials_self_c, finals_self_c, len_self,
                                                chromosomes_y_c, initials_y_c, finals_y_c, len_y, pointer(indices_c),
                                                pointer(initials_result_c), pointer(finals_result_c),
                                                byref(size_result_c))

            # Construct result set
            for i in range(size_result_c.value):
                ci = indices_c[i]
                result.add(GenomicRegion(str(chromosomes_self_c[ci], "ascii"), initials_result_c[i],
                                         finals_result_c[i], name=a.sequences[ci].name,
                                         orientation=a.sequences[ci].orientation, data=a.sequences[ci].data,
                                         proximity=a.sequences[ci].proximity))
            if rm_duplicates:
                result.remove_duplicates()
            return result

    def intersect_count(self, regionset, mode_count="count", threshold=False):
        """Return the number of regions in regionset A&B in following order: (A-B, B-A, intersection)

        *Keyword arguments:*

            - regionset -- the GenomicRegionSet which to compare with.
            - mode_count -- count the number of regions or to measure the length of intersection.
            - threshold -- Define the cutoff of the proportion of the intersecting region (0~50%)

        *Return:*

            - A tupple of numbers: (A-B, B-A, intersection)
        """

        if len(self) == 0:
            return 0, len(regionset), 0
        elif len(regionset) == 0:
            return len(self), 0, 0

        else:
            a = deepcopy(self)
            b = deepcopy(regionset)
            a.merge()
            b.merge()
            if mode_count == "count":

                inter = a.intersect(b, mode=OverlapType.ORIGINAL)
                inter2 = b.intersect(a, mode=OverlapType.ORIGINAL)
                c_a = len(a) - len(inter)
                c_b = len(b) - len(inter2)
                c_ab = len(inter)
                return c_a, c_b, c_ab

            elif mode_count == "bp":
                intersect_r = a.intersect(b, mode=OverlapType.OVERLAP)
                len_inter = intersect_r.total_coverage()
                allbed1 = a.total_coverage()
                allbed2 = b.total_coverage()
                len_12 = allbed1 - len_inter
                len_21 = allbed2 - len_inter
                return len_12, len_21, len_inter

    def closest(self, y, max_dis=10000, return_list=False, top_N=None):
        """Return a new GenomicRegionSet including the region(s) of y which is closest to any self region.
        If there are intersection, return False.

        *Keyword arguments:*

            - y -- the GenomicRegionSet which to compare with
            - max_dis -- maximum distance (default=10000 bp)
            - return_list -- return a list of the distances
            - top_N -- return a dictionary with region names as keys and the GenomicRegionSet containing N clostest regions as values.

        *Return:*

            - A GenomicRegionSet which contains the nearest regions to the self
        """
        if not self.sorted: self.sort()
        if not y.sorted: y.sort()

        targets = self.window(y=y, adding_length=max_dis)

        chroms = self.get_chrom()
        uni_chrom = list(set(chroms))
        target_dict = {}
        for ch in uni_chrom:
            target_dict[ch] = targets.any_chrom(chrom=ch)

        if not top_N:
            if return_list: z_list = []

            z = GenomicRegionSet(self.name)

            for region in self:
                distances = []
                if len(target_dict[region.chrom]) == 0:
                    continue
                elif len(target_dict[region.chrom]) == 1:
                    d = region.distance(target_dict[region.chrom][0])
                    z.add(target_dict[region.chrom][0])
                    if return_list: z_list.append(d)
                else:

                    for r in target_dict[region.chrom]:
                        d = region.distance(r)
                        distances.append(d)
                    min_ind = distances.index(min(distances))
                    z.add(target_dict[region.chrom][min_ind])
                    if return_list: z_list.append(distances[min_ind])

            if return_list:
                return z, z_list
            else:
                return z

        else:
            res_dict = OrderedDict()
            if return_list: res_dist = OrderedDict()

            for region in self:
                if region.name:
                    tag = region.name
                else:
                    tag = region.toString()

                distances = []
                if len(target_dict[region.chrom]) == 0:
                    continue
                elif len(target_dict[region.chrom]) <= top_N:
                    res_dict[tag] = GenomicRegionSet("closest regions to: " + region.name)
                    for t in target_dict[region.chrom]:
                        res_dict[tag].add(t)
                    if return_list: res_dist[tag] = [region.distance(r) for r in target_dict[region.chrom]]
                else:
                    res_dict[tag] = GenomicRegionSet("closest regions to: " + region.name)

                    for r in target_dict[region.chrom]:
                        d = region.distance(r)
                        distances.append(d)

                    sel_ind = sorted(list(range(len(distances))), key=lambda i: distances[i], reverse=False)[:top_N]
                    for s in sel_ind:
                        g = GenomicRegion(chrom=region.chrom, initial=target_dict[region.chrom][s].initial,
                                          final=target_dict[region.chrom][s].final,
                                          name=target_dict[region.chrom][s].name,
                                          orientation=target_dict[region.chrom][s].orientation, data=str(distances[s]))
                        res_dict[tag].add(g)
                    if return_list: res_dist[tag] = distances[sel_ind]

            if return_list:
                return res_dict, res_dist
            else:
                return res_dict

    def remove_duplicates(self, sort=True):
        """
        Remove any duplicate regions, and also returns the sequence list (sorted, by default).
        """
        self.sequences = list(set(self.sequences))

        if sort:
            self.sort()

        return self.sequences

    def window(self, y, adding_length=1000):
        """Return the overlapping regions of self and y with adding a specified number (1000, by default) of base pairs
           upstream and downstream of each region in self. In effect, this allows regions in y that are near regions
           in self to be detected.

        *Keyword arguments:*

            - y -- the GenomicRegionSet which to compare with
            - adding_length -- the length of base pairs added to upstream and downstream of self (default 1000)

        *Return:*

            - A GenomicRegionSet including the regions of overlapping between extended self and original y.
        """
        if len(self) == 0 or len(y) == 0:
            return GenomicRegionSet('None region')
        # Establish an extended GenomicRegionSet
        extended_self = deepcopy(self)
        extended_self.extend(adding_length, adding_length)
        # Find their intersections
        return extended_self.intersect(y)

    def subtract(self, y, whole_region=False, merge=True, exact=False):
        """Return a copy of this GenomicRegionSet minus all the regions overlapping with y.

        *Keyword arguments:*

            - y -- the GenomicRegionSet which to subtract by
            - whole_region -- subtract the whole region, not partially
            - merge -- before subtracting, it merges any overlapping sequence within self or y
            - exact --  only regions which match exactly with a region within y are subtracted
                        if set, whole_region and merge are completely ignored
                        if set, the returned GRS is sorted and does not contain duplicates

        *Return:*

            - A GenomicRegionSet which contains the remaining regions of self after subtraction

        ::

            self     ----------              ------
            y               ----------                    ----
            Result   -------                 ------
        """
        if exact:
            if len(self.sequences) > 1:
                self.sort()
                y.sort()
                small_self = GenomicRegionSet("small_self")
                count_target_regions = 0
                finished = 0
                last_region = self.sequences[-1]

                for region in self.sequences:
                    if not finished:
                        if region != last_region:  # handling of duplicates within self
                            if region.chrom == y.sequences[count_target_regions].chrom:
                                if region.initial == y.sequences[count_target_regions].initial:

                                    if region.final == y.sequences[count_target_regions].final:
                                        # cur y  |-----|
                                        # region |-----|
                                        count_target_regions = count_target_regions + 1
                                        if count_target_regions == len(y.sequences):
                                            finished = 1

                                    elif region.final < y.sequences[count_target_regions].final:
                                        # cur y    |------|
                                        # region   |---|
                                        small_self.add(region)
                                    else:
                                        # cur y   |----|
                                        # region  |------|
                                        loop = 1
                                        while y.sequences[count_target_regions].final < region.final and loop:
                                            if count_target_regions < len(y.sequences) - 1:
                                                count_target_regions = count_target_regions + 1
                                                if y.sequences[count_target_regions].chrom > region.chrom:
                                                    small_self.add(region)
                                                    loop = 0
                                                else:
                                                    # still the same chromosome as current region from self
                                                    if y.sequences[count_target_regions].final == region.final:
                                                        # cur y  |-----|
                                                        # region |-----|
                                                        count_target_regions = count_target_regions + 1
                                                        if count_target_regions == len(y.sequences):
                                                            loop = 0
                                                            finished = 1
                                                    elif y.sequences[count_target_regions].final > region.final:
                                                        # cur y   |----|
                                                        # region  |---|
                                                        small_self.add(region)
                                            else:
                                                # reached the last region in y and this is not exactly mapping
                                                small_self.add(region)
                                                loop = 0

                                elif region.initial < y.sequences[count_target_regions].initial:
                                    # cur y       |---| |        |----| |     |---|     |    |---|
                                    # region    |---|   |  |---|        |   |-----|     |   |------|
                                    small_self.add(region)

                                else:
                                    # cur y     |---|        |  |---|      | |-----| |   |----|
                                    # region          |---|  |     |---|   |   |--|  |     |--|
                                    loop = 1
                                    while y.sequences[count_target_regions].initial < region.initial and loop:
                                        if count_target_regions < len(y.sequences) - 1:
                                            count_target_regions = count_target_regions + 1
                                            if y.sequences[count_target_regions].chrom > region.chrom:
                                                small_self.add(region)
                                                loop = 0
                                            else:
                                                if y.sequences[count_target_regions].final == region.final:
                                                    # cur y    |-----|
                                                    # region   |-----|
                                                    count_target_regions = count_target_regions + 1
                                                    if count_target_regions == len(y.sequences):
                                                        loop = 0
                                                        finished = 1
                                                elif y.sequences[count_target_regions].initial > region.initial:
                                                    # cur y      |---| |        |----| |     |---|     |    |---|
                                                    # region   |---|   |  |---|        |   |-----|     |   |------|
                                                    small_self.add(region)
                                        else:
                                            # reached the last region in y and this is not exactly mapping
                                            small_self.add(region)
                                            loop = 0
                                            finished = 1

                            elif region.chrom > y.sequences[count_target_regions].chrom:
                                count_target_regions = count_target_regions + 1
                                if count_target_regions == len(y.sequences):
                                    finished = 1
                            else:
                                # region.chrom < target.sequences[count_target_regions].chrom:
                                small_self.add(region)
                    else:
                        # finished -> reached end of y, simply add all remaining regions of self that are not equal
                        # to the last subtracted or a region that has already been added
                        if region != last_region:
                            small_self.add(region)
                    last_region = region
                return small_self

            elif len(self.sequences) == 1:
                # GRS only contains 1 region, only check if this matches exactly with any region within y
                for target_region in y.sequences:
                    if target_region == self.sequences[0]:
                        return GenomicRegionSet("small_self")  # return empty GRS
                return self
            else:
                # self is empty GRS
                return self
        else:
            # exact = False

            z = GenomicRegionSet(self.name + ' - ' + y.name)
            if len(self) == 0 or len(y) == 0:
                return self

            if not self.sorted:
                self.sort()

            # if there is overlap within self or y, and the `merge` option is set,
            # we merge any overlapping sequence and create two different GenomicRegionSets
            if merge:
                a = self.merge(w_return=True)
                b = y.merge(w_return=True)
            else:
                a = self
                b = y

            iter_a = iter(a)
            s = next(iter_a)
            last_j = len(b) - 1
            j = 0
            cont_loop = True
            pre_inter = 0
            cont_overlap = False

            while cont_loop:
                # print("Compare: "+s.__repr__()+"\t"+b[j].__repr__())

                # ----------------------
                # -----  --    ----    -----  ----

                # When the regions overlap
                if s.overlap(b[j]):
                    if not cont_overlap: pre_inter = j
                    if whole_region:  # Simply jump to next region
                        try:
                            s = next(iter_a)
                            j = pre_inter
                            cont_overlap = False
                            continue
                        except:
                            cont_loop = False

                    # ------        -----      -------       ---        ------    -----  --
                    #   ------        --        ---       --------  ------       -----  -----
                    if s.initial < b[j].initial:

                        # ------        -----      -------
                        #   ------        --        ---
                        if s.final > b[j].final:
                            s1 = GenomicRegion(chrom=s.chrom, initial=s.initial, final=b[j].initial,
                                               name=s.name, orientation=s.orientation, data=s.data,
                                               proximity=s.proximity)
                            s2 = GenomicRegion(chrom=s.chrom, initial=b[j].final, final=s.final,
                                               name=s.name, orientation=s.orientation, data=s.data,
                                               proximity=s.proximity)
                            z.add(s1)
                            s = s2
                            if j < last_j: j = j + 1
                            cont_overlap = True
                            continue
                        else:
                            s1 = GenomicRegion(chrom=s.chrom, initial=s.initial, final=b[j].initial,
                                               name=s.name, orientation=s.orientation, data=s.data,
                                               proximity=s.proximity)
                            z.add(s1)
                            try:
                                s = next(iter_a)
                                j = pre_inter
                                cont_overlap = False
                                continue
                            except:
                                cont_loop = False

                    elif s.final > b[j].final:

                        #     ------
                        # ------
                        s2 = GenomicRegion(chrom=s.chrom, initial=b[j].final, final=s.final,
                                           name=s.name, orientation=s.orientation, data=s.data, proximity=s.proximity)
                        s = s2
                        if j < last_j: j = j + 1
                        cont_overlap = True
                        continue
                    else:

                        #     ---       -----  --
                        #   --------    -----  -----
                        try:
                            s = next(iter_a)
                            j = pre_inter
                        except:
                            cont_loop = False

                    if j == last_j:
                        try:
                            s = next(iter_a)
                            j = pre_inter
                        except:
                            cont_loop = False
                    else:
                        j = j + 1
                    cont_overlap = True

                elif s < b[j]:
                    z.add(s)
                    try:
                        s = next(iter_a)
                        j = pre_inter
                        cont_overlap = False
                    except:
                        cont_loop = False

                elif s > b[j]:
                    if j == last_j:
                        z.add(s)
                        try:
                            s = next(iter_a)
                        except:
                            cont_loop = False
                    else:
                        j = j + 1
                        cont_overlap = False

            return z

    def subtract_aregion(self, y):
        """Return a GenomicRegionSet excluded the overlapping regions with y.

        *Keyword arguments:*

            - y -- the GenomicRegion which to subtract by

        *Return:*

            - the remaining regions of self after subtraction

        ::

            self     ----------              ------
            y               ----------
            Result   -------                 ------
        """
        if len(self) == 0:
            return GenomicRegionSet('None region')

        else:
            z = GenomicRegionSet('Subtracted RegionSet')
            z.add(y)
            result = self.subtract(z)
            return result

    def mergebyname(self):
        """Merge the regions regardless the intersection by names"""
        names = self.get_names()

        dict_re = {}
        for name in names:
            dict_re[name] = None

        for r in self:
            if dict_re[r.name]:
                if dict_re[r.name].chrom == r.chrom and dict_re[r.name].initial > r.initial:
                    dict_re[r.name].initial = r.initial
                if dict_re[r.name].chrom == r.chrom and dict_re[r.name].final < r.final:
                    dict_re[r.name].final = r.final
            else:
                dict_re[r.name] = r

        z = GenomicRegionSet(self.name)
        for r in list(dict_re.values()):
            z.add(r)

        return z

    def merge(self, w_return=False, namedistinct=False, strand_specific=False):
        """Merge the regions within the GenomicRegionSet

        *Keyword arguments:*

            - w_return -- If TRUE, it returns a GenomicRegionSet; if FALSE, it merges the regions in place.
            - namedistinct -- Merge the regions which have the same names only.
        """
        if not self.sorted:
            self.sort()

        if len(self.sequences) in [0, 1]:
            if w_return:
                return self
            else:
                pass
        else:
            z = GenomicRegionSet(name=self.name)
            prev_region = self.sequences[0]

            if not namedistinct and not strand_specific:
                for cur_region in self.sequences[1:]:
                    if prev_region.overlap(cur_region):
                        prev_region.initial = min(prev_region.initial, cur_region.initial)
                        prev_region.final = max(prev_region.final, cur_region.final)
                    else:
                        z.add(prev_region)
                        prev_region = cur_region
                z.add(prev_region)

            elif namedistinct and not strand_specific:
                for cur_region in self.sequences[1:]:
                    if prev_region.overlap(cur_region) and prev_region.name == cur_region.name:
                        prev_region.initial = min(prev_region.initial, cur_region.initial)
                        prev_region.final = max(prev_region.final, cur_region.final)
                    else:
                        z.add(prev_region)
                        prev_region = cur_region
                z.add(prev_region)

            elif not namedistinct and strand_specific:
                for cur_region in self.sequences[1:]:
                    if prev_region.overlap(cur_region) and prev_region.orientation == cur_region.orientation:
                        prev_region.initial = min(prev_region.initial, cur_region.initial)
                        prev_region.final = max(prev_region.final, cur_region.final)
                    else:
                        z.add(prev_region)
                        prev_region = cur_region
                z.add(prev_region)

            elif namedistinct and strand_specific:
                for cur_region in self.sequences[1:]:
                    if prev_region.overlap(
                            cur_region) and prev_region.name == cur_region.name and prev_region.orientation == cur_region.orientation:
                        prev_region.initial = min(prev_region.initial, cur_region.initial)
                        prev_region.final = max(prev_region.final, cur_region.final)
                    else:
                        z.add(prev_region)
                        prev_region = cur_region
                z.add(prev_region)

            if w_return:
                return z
            else:
                self.sequences = z.sequences

    def combine(self, region_set, change_name=True, output=False):
        """Adding another GenomicRegionSet without merging the overlapping regions.

        *Keyword arguments:*

            - region_set -- the GenomicRegion which to combine with
            - change_name -- Combine the names as a new name for the combined regions
            - output -- If TRUE, it returns a GenomicRegionSet; if FASLSE, it merge the regions in place.
        """
        if output:
            a = GenomicRegionSet(name="")
            for s in self.sequences:
                a.add(s)
            for s in region_set.sequences:
                a.add(s)
            if change_name:
                if a.name == "":
                    a.name = region_set.name
                else:
                    a.name = a.name + " + " + region_set.name
            a.sorted = False
            return a
        else:
            self.sequences.extend(region_set.sequences)
            if change_name:
                if self.name == "":
                    self.name = region_set.name
                else:
                    self.name = self.name + " + " + region_set.name
            self.sorted = False

    def cluster(self, max_distance):
        """Cluster the regions with a certain distance and return the result as a new GenomicRegionSet.

        *Keyword arguments:*

            - max_distance -- the maximum distance between regions within the same cluster

        *Return:*

            - A GenomicRegionSet including clusters

        ::

            self           ----           ----            ----
                              ----             ----                 ----
            Result(d=1)    -------        ---------       ----      ----
            Result(d=10)   ---------------------------------------------
        """

        if not self.sorted: self.sort()

        if len(self) == 0:
            return GenomicRegionSet('None region')
        elif len(self) == 1:
            return self
        else:
            z = GenomicRegionSet('Clustered region set')
            previous = self.sequences[0]
            for s in self.sequences[1:]:
                s_ext = deepcopy(s)
                s_ext.extend(max_distance, max_distance)
                if s_ext.overlap(previous):
                    previous.initial = min(previous.initial, s.initial)
                    previous.final = max(previous.final, s.final)
                else:
                    z.add(previous)
                    previous = s
            z.add(previous)
            return z

    def flank(self, size):
        """Return two flanking intervals with given size from both ends of each region.

        *Keyword arguments:*

            - size -- the length of flanking intervals (default = SAME length as the region)

        *Return:*

            - z -- A GenomicRegionSet including all flanking intervals

        ::

            self        -----           --            ---
            Result -----     -----    --  --       ---   ---
        """
        if len(self) == 0:
            return GenomicRegionSet("Empty")
        else:
            z = GenomicRegionSet("Flanking intervals")
            for s in self:
                s1 = GenomicRegion(name='upstream', chrom=s.chrom,
                                   initial=max(0, s.initial - size),
                                   final=s.initial, data=s.data)
                s2 = GenomicRegion(name='downstream', chrom=s.chrom,
                                   initial=max(0, s.final),
                                   final=s.final + size, data=s.data)  # Adding the limit of chromosome length
                z.add(s1)
                z.add(s2)
            return z

    def jaccard(self, query):
        """Return jaccard index, a value of similarity of these two GenomicRegionSet.

        *Keyword arguments:*

            - query -- the GenomicRegionSet which to compare with.

        *Return:*

            - similarity -- (Total length of overlapping regions)/(Total length of original regions)

        ::

            self              --8--      ---10---    -4-
            query        ---10---             ---10---
            intersect         -5-             -4-    2
            similarity: (5+4+2)/[(8+10+4)+(10+10)-(5+4+2)] = 11/31
        """

        # if sys.platform == "darwin":
        #     return self.jaccard_python(query)
        # else:
        return self.jaccard_c(query)

    def jaccard_python(self, query):
        a = copy.deepcopy(self)
        b = copy.deepcopy(query)
        if a.total_coverage() == 0 and len(a) > 0:
            print(" ** Warning: \t" + a.name + " has zero length.")
            return a.name
        if b.total_coverage() == 0 and len(b) > 0:
            print(" ** Warning: \t" + b.name + " has zero length.")
            return b.name

        intersects = a.intersect(b)
        intersects.merge()
        inter = intersects.total_coverage()

        a.combine(b, change_name=False)
        a.merge()
        uni = a.total_coverage()
        similarity = inter / uni
        return similarity

    def jaccard_c(self, query):
        # Determine path of shared library
        Lib = LibraryPath()
        lib = cdll.LoadLibrary(Lib.get_c_rgt())

        # Bind library
        ctypes_jaccardC = lib.jaccard

        # Specify data types
        ctypes_jaccardC.argtypes = [POINTER(c_char_p), POINTER(c_int), POINTER(c_int), c_int, POINTER(c_char_p),
                                    POINTER(c_int), POINTER(c_int), c_int]
        ctypes_jaccardC.restype = c_double

        if not self.sorted:
            self.sort()
        if not query.sorted:
            query.sort()

        assert self.sorted
        assert query.sorted

        # Convert to ctypes
        chroms_self_python = [gr.chrom.encode("utf-8") for gr in self.sequences]
        chroms_self_c = (c_char_p * len(chroms_self_python))(*chroms_self_python)
        # print('Converted self.chroms to c', str(chroms_self_python[:4]), '...')

        chroms_query_python = [gr.chrom.encode("utf-8") for gr in query.sequences]
        chroms_query_c = (c_char_p * len(chroms_query_python))(*chroms_query_python)
        # print('Converted query.chroms to c', str(chroms_query_python[:4]), '...')

        initials_self_python = [gr.initial for gr in self.sequences]
        initials_self_c = (c_int * len(initials_self_python))(*initials_self_python)
        # print('Converted self.initials to c', str(initials_self_python[:4]), '...')

        initials_query_python = [gr.initial for gr in query.sequences]
        initials_query_c = (c_int * len(initials_query_python))(*initials_query_python)
        # print('Converted query.initials to c', str(initials_query_python[:4]), '...')

        finals_self_python = [gr.final for gr in self.sequences]
        finals_self_c = (c_int * len(finals_self_python))(*finals_self_python)
        # print('Converted self.finals to c', str(finals_self_python[:4]), '...')

        finals_query_python = [gr.final for gr in query.sequences]
        finals_query_c = (c_int * len(finals_query_python))(*finals_query_python)
        # print('Converted query.finals to c', str(finals_query_python[:4]), '...')

        # print('Converted to ctypes')

        # Call C-function
        return ctypes_jaccardC(chroms_self_c, initials_self_c, finals_self_c, len(self), chroms_query_c,
                               initials_query_c, finals_query_c, len(query))

    def within_overlap(self):
        """Check whether there is overlapping within or not."""
        refer_posi = GenomicRegion(name="reference", chrom="chr1", initial=0, final=0)
        if not self.sorted: self.sort()

        for s in self:
            if s.overlap(refer_posi):
                return True
            refer_posi = s
        return False

    def total_coverage(self):
        """Return the sum of all lengths of regions."""
        length = 0
        for s in self:
            try:
                length = length + len(s)
            except:
                print("Warning: cannot get length of " + str(s) + " skipped")
                continue
        return length

    def get_genome_data(self, organism, chrom_X=True, chrom_Y=False, chrom_M=False):
        """Add genome data from database into the GenomicRegionSet.

        *Keyword arguments:*

            - organism -- Define the organism
            - chrom_X -- Include chromosome X
            - chrom_Y -- Include chromosome Y
            - chrom_M -- Include mitochondrial chromosome
        """
        genome = GenomeData(organism)
        chromosome_file = open(genome.get_chromosome_sizes(), 'r')
        for line in chromosome_file:
            if "random" not in line and "_" not in line:
                if not chrom_X and "chrX" in line: continue
                if not chrom_Y and "chrY" in line: continue
                if not chrom_M and "chrM" in line: continue
                chrom_region = GenomicRegion(chrom=line.split("\t")[0],
                                             initial=0,
                                             final=int(line.split("\t")[1]))
                self.add(chrom_region)
                continue
        chromosome_file.close()

    def random_regions(self, organism, total_size=None, multiply_factor=1,
                       overlap_result=True, overlap_input=True,
                       chrom_X=False, chrom_M=False, filter_path=None):
        """Return a GenomicRegionSet which contains the random regions generated by given entries and given number
           on the given organism.
        
        *Keyword arguments:*

            - organism -- Define organism's genome to use. (hg19, mm9)
            - total_size -- Given the number of result random regions.
            - multiply_factor -- This factor multiplies to the number of entries is the number of exporting random
              regions. ** total_size has higher priority than multiply_factor. **
            - overlap_result -- The results whether overlap with each other or not. (True/False)
            - overlap_input -- The results whether overlap with input entries or not. (True/False)
            - chrom_X -- The result covers chromosome X or not. (True/False)
            - chrom_M -- The result covers mitochondria chromosome or not. (True/False)
            - filter_path -- Given the path of filter BED file
        
        *Return:*

            - z -- A GenomicRegionSet which contains the random regions
        """

        def list_max(chrom_list, result_map):
            """Generate a list containing maximum length for each chromosome."""
            map_max = []  # Store the minimum length corresponding to chrom_list
            for ch in chrom_list:
                map_max.append(max([len(s) for s in result_map.any_chrom(ch)]))
            return map_max

        def randoming(length, result_map, chrom_list, map_max, choices):
            """Return a new GenomicRegion as the result of randomization."""
            cont_loop = True
            candidate_chrom = [chrom_list[i] for i, l in enumerate(map_max) if l > length]

            while True:
                ch = weighted_choice(choices)
                if ch in candidate_chrom:
                    break

            while cont_loop:
                # try:
                regions = result_map.any_chrom(ch, len_min=length)
                sample = random.choice(regions)
                if len(sample) < length:
                    print("1", end="")
                    continue
                else:
                    # print("2", end="")
                    random_posi = random.randint(sample.initial, sample.final - length)
                    # print(sample.toString())
                    # print(random_posi)
                    cont_loop = False

                # except:
                # print("3", end="")
                #    continue

            return GenomicRegion(chrom=sample.chrom, initial=random_posi, final=random_posi + length)

            # raise Exception("There is no further space for randomization on the genome.")

        def weighted_choice(choices):
            total = sum(w for c, w in choices)
            r = random.uniform(0, total)
            upto = 0
            for c, w in choices:
                if upto + w > r:
                    return c
                upto += w
            assert False, "Shouldn't get here"

        ############################################
        # input_map = copy.deepcopy(self)

        input_map = self
        input_num = len(self)
        input_list = [len(x) for x in self.sequences]
        result_list = []
        # Total number and lengths of random regions
        if total_size:
            for i in range(int(total_size)):
                result_list.append(input_list[i % input_num])
        elif multiply_factor > 0:
            for i in range(int(multiply_factor * input_num)):
                result_list.append(input_list[i % input_num])

        # Maps
        # Fetching the chromosome length from data
        chrom_map = GenomicRegionSet("chrom_map")
        chrom_map.get_genome_data(organism, chrom_X=chrom_X, chrom_M=chrom_M)
        chrom_list = chrom_map.get_chrom()
        if filter_path:
            filter_map = GenomicRegionSet('filter')
            filter_map.read(filter_path)
            chrom_map = chrom_map.subtract(filter_map)

        if overlap_input:
            result_map = chrom_map
        else:
            result_map = chrom_map.subtract(input_map)

        z = GenomicRegionSet(name="random regions")
        map_max = list_max(chrom_list, result_map)

        # Generate pk which stores the total length of all regions in each chromosome
        choices = []
        for ch in chrom_list:
            cov = result_map.any_chrom(ch)
            pk = sum([len(s) for s in cov])

            choices.append([ch, pk])

        for length in result_list:
            new_region = randoming(length, result_map, chrom_list, map_max, choices)
            z.add(new_region)
            if not overlap_result:
                result_map = result_map.subtract_aregion(new_region)
                map_max[chrom_list.index(new_region.chrom)] = max(
                    [len(s) for s in result_map.any_chrom(new_region.chrom)])
                choices[chrom_list.index(new_region.chrom)][1] -= len(new_region)
        return z

    def trim_by(self, background):
        """Trim a GenomicRegionSet by a given background, another GenomicRegionSet."""
        return self.intersect(background, mode=OverlapType.ORIGINAL)
        # self = self.intersect(background, mode = OverlapType.OVERLAP)
        # self = s

    def projection_test(self, query, organism, extra=None, background=None):
        """"Return the p value of binomial test.

        *Keyword arguments:*

            - query -- A GenomicRegionSet as query
            - organism -- Define the organism
            - extra -- Return the extra statistics
            - background -- Use a GenomicRegionSet as the background
            - return_intersected_query -- Return a GenomicRegionSet containing the intersected regions of query

        *Return:*

            - if extra=True, returns (possibility, ration, p-value, intersected_query)
            - if extra=False, returns p-value
        """
        chrom_map = GenomicRegionSet("Genome")
        chrom_map.get_genome_data(organism=organism)
        if self.total_coverage() == 0 and len(self) > 0:
            print(" ** Warning: \t" + self.name + " has zero length.")
            if extra:
                return 0, 0, "na"
            else:
                return "na"
        if query.total_coverage() == 0 and len(query) > 0:
            query.extend(0, 1)
        # print("coverage of reference: ",self.total_coverage(),"\tcoverage of genome: ",chrom_map.total_coverage())
        if background:  # backgound should be a GenomicRegionSet
            ss = self.intersect(background, OverlapType.OVERLAP)
            possibility = ss.total_coverage() / background.total_coverage()
        else:
            possibility = self.total_coverage() / chrom_map.total_coverage()  # The average likelihood

        nquery = query.relocate_regions(center='midpoint', left_length=0, right_length=0)
        intersect_regions = nquery.intersect(self, mode=OverlapType.ORIGINAL)
        if extra:
            intersect_q = query.intersect(self, mode=OverlapType.ORIGINAL)
        n = len(nquery)
        k = len(intersect_regions)
        try:
            r = k / n
        except:
            r = 0
        # print("intersections: ",k,"\tnumber of query",n,"\tgenetic coverage: ",possibility)
        p = float(stats.binom_test(k, n, possibility))
        if extra:
            return possibility, r, p, intersect_q
        else:
            return p

    def any_chrom(self, chrom, len_min=False, len_max=False,
                  return_list=True, return_regionset=False):
        """Return a list of regions which belongs to given chromosome.

        *Keyword arguments:*

            - chrom -- Define chromosome
            - len_min -- minimum length
            - len_max -- maximum length

        *Return:*

            - A list of regions which belongs to given chromosome.
        """
        if not len_min and not len_max:
            res = [s for s in self if s.chrom == chrom]
        elif len_min > 0 and not len_max:
            res = [s for s in self if s.chrom == chrom and len(s) >= len_min]
        elif len_max > 0 and not len_min:
            res = [s for s in self if s.chrom == chrom and len(s) <= len_max]
        else:
            res = [s for s in self if s.chrom == chrom and len_min <= len(s) <= len_max]

        if return_list:
            return res
        elif return_regionset:
            z = GenomicRegionSet(chrom)
            z.sequences = res
            return z
        else:
            self.sequences = res

    def relocate_regions(self, center='midpoint', left_length=2000, right_length=2000):
        """Return a new GenomicRegionSet which relocates the regions by given center and extend length.
        
        *Keyword arguments:*

            - center -- Define the referring point of each region

                1. midpointra -- locate the new region's center as original region's midpoint
                2. leftend -- locate the new region's center as original region's 5' end (if no orientation information, default is left end)
                3. rightend -- locate the new region's center as original region's 3' end (if no orientation information, default is right end)
                4. bothends -- locate the new region's center as original region's both ends
                5. downstream -- rightend in positive strand and leftend in negative strand
                6. upstream -- leftend in positive strand and rightend in negative strand

            - left_length -- Define the length to extend on the left side
            - right_length -- Define the length to extend on the right side
        """
        new_regions = GenomicRegionSet("relocated_" + self.name)
        for r in self.sequences:
            # Define the position
            if center == 'midpoint':
                mp = int(0.5 * (r.initial + r.final))
                nr = GenomicRegion(chrom=r.chrom, initial=mp, final=mp, name=r.name,
                                   orientation=r.orientation, data=r.data, proximity=r.proximity)

            elif center == 'leftend':
                nr = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial, name=r.name,
                                   orientation=r.orientation, data=r.data, proximity=r.proximity)
            elif center == 'rightend':
                nr = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final, name=r.name,
                                   orientation=r.orientation, data=r.data, proximity=r.proximity)
            elif center == 'downstream':
                if r.orientation == "-":
                    nr = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial, name=r.name,
                                       orientation=r.orientation, data=r.data, proximity=r.proximity)
                else:
                    nr = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final, name=r.name,
                                       orientation=r.orientation, data=r.data, proximity=r.proximity)
            elif center == 'upstream':
                if r.orientation == "-":
                    nr = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final, name=r.name,
                                       orientation=r.orientation, data=r.data, proximity=r.proximity)
                else:
                    nr = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial, name=r.name,
                                       orientation=r.orientation, data=r.data, proximity=r.proximity)
            elif center == 'bothends':
                nrL = GenomicRegion(chrom=r.chrom, initial=r.initial, final=r.initial, name=r.name,
                                    orientation=r.orientation, data=r.data, proximity=r.proximity)  # newbed on left end
                nrR = GenomicRegion(chrom=r.chrom, initial=r.final, final=r.final, name=r.name,
                                    orientation=r.orientation, data=r.data,
                                    proximity=r.proximity)  # newbed on right end

            if center == 'bothends':
                new_regions.add(nrL)
                new_regions.add(nrR)
            else:
                new_regions.add(nr)

        # Extend the region
        new_regions.extend(left_length, right_length)
        return new_regions

    def maximum_length(self):
        """Return the length of the maximum region from the GenomicRegionSet."""
        try:
            maxl = len(self.sequences[0])
            for s in self:
                maxl = max(len(s), maxl)
            return maxl
        except:
            print("There is no region in the given GenomicRegionSet, {}".format(self.name))
            return 0

    def include(self, region):
        """Check whether the given region has intersect with the original regionset.

        *Keyword arguments:*

            - region -- A GenomicRegion to be checked.
        """
        for s in self:
            if s.overlap(region):
                return True
            else:
                continue
        return False

    def complement(self, organism, chrom_X=True, chrom_Y=False, chrom_M=False):
        """Return the complement GenomicRegionSet for the given organism.

        *Keyword arguments:*

            - organism -- Define organism's genome to use. (hg19, mm9)
            - chrom_X -- The result covers chromosome X or not. (True/False)
            - chrom_Y -- The result covers chromosome Y or not. (True/False)
            - chrom_M -- The result covers mitochondrial chromosome or not. (True/False)

        *Return:*

            - z -- A GenomicRegionSet which contains the complement regions
        """
        g = GenomicRegionSet("complement_" + self.name)
        g.get_genome_data(organism, chrom_X, chrom_Y, chrom_M)
        h = g.subtract(self)
        return h

    def count_by_region(self, region):
        """Return the number of intersection regions with the given GenomicRegion.
        
        *Keyword arguments:*

            - region -- A GenomicRegion defining the interval for counting.
        """
        query = GenomicRegionSet("query")
        query.add(region)
        return len(self.intersect(query))

    def count_by_regionset(self, regionset):
        """Return the number of intersection regions with the given GenomicRegionSet.

        *Keyword arguments:*

            - regionset -- A GenomicRegionSet defining the interval for counting.
        """
        return len(self.intersect(regionset, mode=OverlapType.ORIGINAL))

    def counts_per_region(self, regionset):
        """Return a list of counting numbers of the given GenomicRegionSet based on the self.
        
        *Keyword arguments:*

            - regionset -- A GenomicRegionSet defining the interval for counting.

        .. note:: The length of the result list is the same as self GenomicRegionSet
        """
        if len(self) == 0: return None
        if len(regionset) == 0: return [0] * len(self)

        # a = copy.deepcopy(con_self)
        # b = copy.deepcopy(regionset)
        if not self.sorted: self.sort()
        if not regionset.sorted: regionset.sort()
        counts = []

        iter_a = iter(self)
        s = next(iter_a)
        last_j = len(regionset) - 1
        j = 0
        cont_loop = True
        pre_inter = 0
        cont_overlap = False
        c = 0
        while cont_loop:
            # When the regions overlap
            # print(str(s)+"\t\t"+str(b[j]))
            if s.overlap(regionset[j]):
                c += 1
                if not cont_overlap: pre_inter = j
                if j == last_j:
                    try:
                        s = next(iter_a)
                        counts.append(c)
                        c = 0
                        j = pre_inter
                    except:
                        cont_loop = False
                        counts.append(c)
                else:
                    j += 1
                cont_overlap = True

            elif s < regionset[j]:
                try:
                    s = next(iter_a)
                    counts.append(c)
                    c = 0
                    j = pre_inter
                    cont_overlap = False
                except:
                    cont_loop = False
                    counts.append(c)

            elif s > regionset[j]:
                if j == last_j:
                    try:
                        s = next(iter_a)
                        counts.append(c)
                        c = 0

                    except:
                        cont_loop = False
                        counts.append(c)
                else:
                    j = j + 1
                    cont_overlap = False
        return counts

    def covered_by_aregion(self, region):
        """Return a GenomicRegionSet which includes all the regions covered by a given region.

        *Keyword arguments:*

            - region -- A GenomicRegion defining the interval.

        *Return:*

            - A GenomicRegionSet containing the regions within the defined interval.
        """
        region_set = GenomicRegionSet("Query")
        region_set.add(region)
        return self.intersect(region_set, mode=OverlapType.ORIGINAL)

    def replace_region_name(self, regions, combine=False):
        """Replace the region names by the given GenomicRegionSet.

        *Keyword arguments:*

            - regions -- A GenomicRegionSet as the source for the names.
            - combine -- Combine the names from the old and new regions.
        """

        if len(self) == 0 or len(regions) == 0:
            return

        else:
            # res = copy.deepcopy(regions)
            # res = regions
            # If there is overlap within self or y, they should be merged first.

            if not self.sorted: self.sort()
            if not regions.sorted: regions.sort()

            iter_a = iter(self)
            s = next(iter_a)
            last_j = len(regions) - 1
            j = 0
            cont_loop = True

            while cont_loop:
                # print(str(s),"\t",str(b[j]))
                # When the res overlap
                if s.overlap(regions[j]):
                    if combine:
                        s.name = s.name + "_" + regions[j].name
                    else:
                        s.name = regions[j].name
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False

                elif s < regions[j]:
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False
                elif s > regions[j]:
                    if j == last_j:
                        cont_loop = False
                    else:
                        j = j + 1
                else:
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False
            return

    def replace_region_strand(self, regions=None, reverse=False, all=None):
        """Replace the region strand by the given GenomicRegionSet.

        *Keyword arguments:*

            - regions -- A GenomicRegionSet as the source for the strand.
        """
        if not regions and reverse:
            for r in self:
                if r.orientation == "+":
                    r.orientation = "-"
                elif r.orientation == "-":
                    r.orientation = "+"
            return
        elif all == "+" or all == "-":
            for r in self:
                r.orientation = all
            return

        elif len(self) == 0 or len(regions) == 0:
            return

        elif regions:
            if not self.sorted: self.sort()
            if not regions.sorted: regions.sort()

            iter_a = iter(self)
            s = next(iter_a)
            last_j = len(regions) - 1
            j = 0
            cont_loop = True

            while cont_loop:
                # When the res overlap
                if s.overlap(regions[j]):
                    if not reverse:
                        s.orientation = regions[j].orientation
                    elif reverse and regions[j].orientation == "+":
                        s.orientation = "-"
                    elif reverse and regions[j].orientation == "-":
                        s.orientation = "+"
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False

                elif s < regions[j]:
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False
                elif s > regions[j]:
                    if j == last_j:
                        cont_loop = False
                    else:
                        j = j + 1
                else:
                    try:
                        s = next(iter_a)
                    except:
                        cont_loop = False
            return

    def change_name_by_dict(self, convert_dict):
        """Change the names of each region by the given dictionary.

        *Keyword arguments:*

            - convert_dict -- A dictionary having original names as its keys and new names as its values.
        """
        z = GenomicRegionSet(self.name)
        for s in self:
            try:
                n = convert_dict[s.name]
                if not n: n = s.name
            except:
                n = s.name
            z.add(GenomicRegion(s.chrom, s.initial, s.final, n, s.orientation, s.data, s.proximity))
        return z

    def by_names(self, names, load_score=False, background=False):
        """Subset the GenomicRegionSet by the given list of names.

        *Keyword arguments:*

            - names -- A list of names as targets or a GeneSet.

        *Return:*

            - A GenomicRegionSet containing the regions with the target names.
        """
        z = GenomicRegionSet(self.name)
        if isinstance(names, list):
            targets = names
        elif isinstance(names.genes, list):
            targets = names.genes
        targets = [x.upper() for x in targets]
        for gr in self:

            if gr.name.upper() in targets and not background:
                if load_score:
                    d = gr.data.split()
                    gr.data = "\t".join([str(names.values[gr.name.upper()])] + d[1:])
                z.add(gr)
            elif gr.name.upper() not in targets and background:
                z.add(gr)
        return z

    def coverage_per_region(self, regionset):
        """Return a list of coverage of the given GenomicRegionSet based on the self GenomicRegionSet.

        *Keyword arguments:*

            - regionset -- A GenomicRegionSet as the signal for calculate the coverage.

        .. note:: The length of the result list is the same as self GenomicRegionSet.
        """
        if len(self) == 0: return None
        if len(regionset) == 0: return [0] * len(self)

        if not self.sorted: self.sort()
        if not regionset.sorted: regionset.sort()
        coverages = []

        iter_a = iter(self)
        s = next(iter_a)
        last_j = len(regionset) - 1
        j = 0
        cont_loop = True
        pre_inter = 0
        cont_overlap = False
        c = GenomicRegionSet("coverage")

        while cont_loop:
            # When the regions overlap
            # print(str(s)+"\t\t"+str(b[j]))
            if s.overlap(regionset[j]):

                c.add(regionset[j])
                if not cont_overlap: pre_inter = j
                if j == last_j:
                    coverages.append(c.total_coverage() / len(s))
                    try:
                        s = next(iter_a)
                        c = GenomicRegionSet("coverage")
                        j = pre_inter
                    except:
                        cont_loop = False
                else:
                    j += 1
                cont_overlap = True

            elif s < regionset[j]:
                overlapping = False
                coverages.append(c.total_coverage() / len(s))
                try:
                    s = next(iter_a)
                    c = GenomicRegionSet("coverage")
                    j = pre_inter
                    cont_overlap = False
                except:
                    cont_loop = False

            elif s > regionset[j]:
                overlapping = False
                if j == last_j:
                    coverages.append(c.total_coverage() / len(s))
                    try:
                        s = next(iter_a)
                        c = GenomicRegionSet("coverage")

                    except:
                        cont_loop = False
                else:
                    j = j + 1
                    cont_overlap = False
        return coverages

    def extract_blocks(self, keep_name=False):
        """Extract the exon information from self.data and add them into the self GenomicRegionSet."""
        regions = []
        for rg in self:
            try:
                z = rg.extract_blocks(keep_name)
            except:
                z = [rg]
            regions = regions + z
        self.sequences = regions

    def sort_score(self):
        """Sort the regions by their scores."""
        self.sort(key=lambda x: float(x.data.split("\t")[0]), reverse=True)

    def filter_strand(self, strand="+"):
        """Return the defined strands"""
        z = GenomicRegionSet(self.name)
        for region in self:
            if region.orientation == strand:
                z.add(region)
        return z

    def add_associated_gene_data(self, organism):
        """Add the associated gene symbol to data"""
        a = self.gene_association(organism=organism, show_dis=True)
        for i, r in enumerate(self):
            r.data = r.data + "\t" + a[i].name

    def longest_region(self, return_set=False):
        """Return the longest region(s)"""
        length_list = [len(region) for region in self]
        max_len = max(length_list)
        longest_ind = [i for i, j in enumerate(length_list) if j == max_len]
        if not return_set:
            return self[longest_ind[0]]
        else:
            z = GenomicRegionSet("longest regions")
            for i in longest_ind:
                z.add(self[i])

    def sample_close_region(self, background, min_dis=200, len_each=100):
        """
        Return a GenomicRegionSet which includes the regions close to the self and within the background
        under the limit of distance.
        """

        def random_choose(col_regionset):
            options = []
            for i, region in enumerate(col_regionset):
                if len(region) > len_each:
                    options.append(i)
            if options:
                ind = random.sample(options, 1)[0]
                ssite = random.sample(list(range(col_regionset[ind].initial, col_regionset[ind].final - len_each)), 1)[
                    0]
                res = GenomicRegion(chrom=col_regionset[ind].chrom,
                                    initial=ssite,
                                    final=ssite + len_each,
                                    orientation=col_regionset[ind].orientation)
            else:
                print("Find not enough space in " + col_regionset[0].chrom +
                      ":" + str(col_regionset[0].initial) + "-" + str(
                    col_regionset[-1].final) + "\tSample from the largest region.")
                longest = col_regionset.longest_region()
                mid = int(0.5 * (longest.final + longest.initial))
                res = GenomicRegion(chrom=longest.chrom,
                                    initial=mid - int(0.5 * len_each),
                                    final=mid + int(0.5 * len_each),
                                    orientation=longest.orientation)
            return res

        a = background.intersect(self, mode=OverlapType.ORIGINAL)
        b = self.extend(left=min_dis, right=min_dis, w_return=True)
        b.merge()
        c = a.subtract(b)
        # Iteration
        iter_a = iter(a)
        sa = next(iter_a)
        iter_c = iter(c)
        sc = next(iter_c)
        # Loop
        z = GenomicRegionSet("sample")
        q_coll = GenomicRegionSet(sa.toString())
        cont_loop = True
        while cont_loop:
            # print("\t".join([sa.toString(), sc.toString()]))
            if sa.overlap(sc):
                q_coll.add(sc)
                try:
                    sc = next(iter_c)
                except:
                    if len(q_coll):
                        z.add(random_choose(col_regionset=q_coll))
                    try:
                        sa = next(iter_a)
                    except:
                        cont_loop = False
            elif sa < sc:
                if len(q_coll):
                    z.add(random_choose(col_regionset=q_coll))
                q_coll = GenomicRegionSet(sa.toString())
                try:
                    sa = next(iter_a)
                except:
                    cont_loop = False

            else:
                if len(q_coll):
                    z.add(random_choose(col_regionset=q_coll))
                q_coll = GenomicRegionSet(sa.toString())
                try:
                    sc = next(iter_c)
                except:
                    cont_loop = False

        return z

    def split_by_chromosome(self):
        """Return a list of many GenomicRegionSets, each one for one chromosome """
        chrom_list = set(self.get_chrom())
        z = []
        for c in chrom_list:
            g = self.any_chrom(chrom=c, return_list=False, return_regionset=True)
            z.append(g)
        return z

    def average_size(self):
        """Return the average size of the regions"""
        size = [len(r) for r in self.sequences]
        return sum(size) / len(size)

    def median_size(self):
        """Return the average size of the regions"""
        size = [len(r) for r in self.sequences]
        return numpy.median(size)

    def max_size(self):
        """Return the maximum size of the regions"""
        size = [len(r) for r in self.sequences]
        return max(size)

    def min_size(self):
        """Return the minimum size of the regions"""
        size = [len(r) for r in self.sequences]
        return min(size)

    def size_variance(self):
        """Return the average size of the regions"""
        import numpy as np
        size = [len(r) for r in self.sequences]
        return np.std(size)

    def filter_by_size(self, maximum=None, minimum=1):
        """Return a GenomicRegionSet containing filtered regions by the given limits. """
        z = GenomicRegionSet("filtered")
        for r in self:
            if maximum:
                if minimum < len(r) < maximum:
                    z.add(r)
            else:
                if minimum < len(r):
                    z.add(r)
        return z

    def get_distance(self, y, ignore_overlap=False, strand_specific=False, thresh_dist=50000):
        """Return a list of distances between the closest regions from two region sets."""
        if not self.sorted: self.sort()
        if not y.sorted: y.sort()

        last_j = len(y) - 1
        j = 0
        if strand_specific:
            pre_inter = [0, 0]
        else:
            pre_inter = 0
        res = []
        for s in self:
            cont_loop = True
            cont_overlap = False
            asso_names = {"overlap": [], "close_l": [], "close_r": []}
            while cont_loop:
                if strand_specific and s.orientation != y[j].orientation:
                    if j == last_j:
                        cont_loop = False
                    else:
                        j += 1
                else:
                    d = s.distance(y[j])
                    if d == 0:
                        asso_names["overlap"].append([s.name, str(0), y[j].name])
                        if not cont_overlap:
                            if strand_specific and s.orientation == "+":
                                pre_inter[0] = j
                            elif strand_specific and s.orientation == "-":
                                pre_inter[1] = j
                            elif not strand_specific:
                                pre_inter = j
                        if j == last_j:
                            cont_loop = False
                        else:
                            j += 1
                            cont_overlap = True
                    elif asso_names["overlap"] and d != 0:
                        if strand_specific:
                            if pre_inter[0] > 0 and pre_inter[1] > 0:
                                j = min(pre_inter)
                            elif pre_inter[0] == 0 and pre_inter[1] > 0:
                                j = pre_inter[1]
                            elif pre_inter[0] > 0 and pre_inter[1] == 0:
                                j = pre_inter[0]
                        elif s.chrom == y[j].chrom and pre_inter > 0:
                            j = pre_inter
                        cont_loop = False
                    elif 0 < d < thresh_dist:
                        if s > y[j]:
                            asso_names["close_l"] = [[s.name, "-" + str(d), y[j].name]]
                        elif s < y[j]:
                            asso_names["close_r"] = [[s.name, "+" + str(d), y[j].name]]
                            cont_loop = False
                        if j == last_j:
                            cont_loop = False
                        else:
                            j += 1
                    elif s < y[j]:
                        if strand_specific and s.orientation == "+":
                            if s.chrom == y[j].chrom and pre_inter[0] > 0:
                                j = pre_inter[0]
                        elif strand_specific and s.orientation == "-":
                            if s.chrom == y[j].chrom and pre_inter[1] > 0:
                                j = pre_inter[1]
                        elif s.chrom == y[j].chrom and pre_inter > 0:
                            j = pre_inter
                        cont_loop = False
                    elif s > y[j]:
                        if j == last_j:
                            cont_loop = False
                        else:
                            j += 1

            if asso_names["overlap"] and not ignore_overlap:
                last_one = ["x", 0, "x"]
                for line in asso_names["overlap"]:
                    if line[2] == last_one[2]:
                        continue
                    else:
                        res += [line]
                        last_one = line
            if asso_names["close_l"]:
                res += asso_names["close_l"]
            if asso_names["close_r"]:
                res += asso_names["close_r"]

        return res

    def cut_regions(self, y, keep="upstream"):
        z = self.subtract(y)
        genes = {}
        for r in z:
            if r.name in genes:
                if keep == "upstream" and r.orientation == "+" and r.initial < genes[r.name].initial:
                    genes[r.name] = r
                elif keep == "upstream" and r.orientation == "-" and r.final > genes[r.name].final:
                    genes[r.name] = r
                elif keep == "downstream" and r.orientation == "+" and r.initial > genes[r.name].initial:
                    genes[r.name] = r
                elif keep == "downstream" and r.orientation == "-" and r.initial.final < genes[r.name].final:
                    genes[r.name] = r
            else:
                genes[r.name] = r

        zz = GenomicRegionSet("out")
        for k in list(genes.values()):
            zz.add(k)
        return zz

    def get_score_dict(self):
        """Get a dictionary of scores"""
        d = {}
        for r in self:
            #
            if isinstance(r.data, str):
                d[r.toString()] = float(r.data.split("\t")[0])
            elif isinstance(r.data, float) or isinstance(r.data, int):
                d[r.toString()] = r.data

        # except:
        # continue
        return d

    def standard_chrom(self):
        """Remove the random chromosomes and leave only the standard chromosomes, e.g. chr1-22, chrX"""
        z = GenomicRegionSet(self.name)
        for r in self:
            if "_" not in r.chrom:
                z.add(r)
        return z

    def get_promoters(self, length=1000):
        promoters = GenomicRegionSet("promoters")
        for s in self:
            if s.orientation == "+":
                s.initial, s.final = max(s.initial - length, 0), s.initial
            else:
                s.initial, s.final = s.final, s.final + length
            promoters.add(s)
        return promoters

    def get_GeneSet(self):
        genes = GeneSet(self.name)
        for r in self:
            genes.add(gene_name=r.name, value=float(r.data.split("\t")[0]))
        return genes

    def load_from_list(self, loci_list):
        """Load the regions from a list, such as [['chr1', 1000, 1500, '+']]"""
        for l in loci_list:
            self.add(GenomicRegion(chrom=l[0], initial=int(l[1]), final=int(l[2]), orientation=l[3]))

    def map_names(self, target, strand=False, convert_nt=False):
        """Return a list of the target names overlapping the regions in the self in order"""
        if len(target) == 0:
            return ["."] * len(self)
        elif len(self) == 0:
            return None
        else:
            names = []
            convert_dic = {"A": "T", "T": "A", "C": "G", "G": "C"}
            iter_a = iter(self)
            s = next(iter_a)
            last_j = len(target) - 1
            j = 0
            cont_loop = True
            # pre_j =
            # print(target[0].name)
            if target[0].name:
                if convert_nt and ")n" not in target[0].name:
                    convert_nt = False
            elif target[0].name == None:
                convert_nt = False
                # for r in target:
                #     r.name = ""

            while cont_loop:
                # When the regions
                if not target[j].name:
                    target[j].name = "."
                if s.overlap(target[j]):
                    if strand:
                        if s.orientation == target[j].orientation:
                            names.append(target[j].name)
                            try:
                                s = next(iter_a)
                                # j = pre_j
                            except:
                                cont_loop = False
                        else:
                            if j == last_j:
                                names.append(".")
                                cont_loop = False
                            else:
                                j += 1
                    elif not strand:
                        if convert_nt and s.orientation == "-":
                            seq = target[j].name.partition("(")[2].partition(")")[0]
                            nseq = [convert_dic[r] for r in seq]
                            n = "(" + "".join(nseq) + ")n"

                        else:
                            n = target[j].name
                        names.append(n)
                        try:
                            s = next(iter_a)
                            # j = pre_j
                        except:
                            cont_loop = False
                    else:
                        if j == last_j:
                            names.append(".")
                            cont_loop = False
                        else:
                            j += 1
                elif s < target[j]:
                    names.append(".")
                    try:
                        s = next(iter_a)
                        # j = pre_j
                    except:
                        cont_loop = False
                elif s > target[j]:
                    # pre_j = j
                    if j == last_j:
                        names.append(".")
                        cont_loop = False
                    else:
                        j += 1
                else:
                    names.append(".")
                    try:
                        s = next(iter_a)
                        # j = pre_j
                    except:
                        cont_loop = False
            # print([len(self), len(names)])
            while len(names) < len(self):
                # print(".", end="")
                names.append(".")

            return names

    def is_stranded(self):
        """Check the GenomicRegionSet is stranded or not"""
        if self.sequences[0].orientation == "+" or self.sequences[0].orientation == "-":
            return True
        else:
            return False

    def fragmentize(self, size):
        """Fragmentize each region by the given size"""
        z = GenomicRegionSet(self.name)
        for r in self:
            l = len(r)
            if l < size:
                z.add(r)
            else:
                for i in range(int(l / size)):
                    ff = r.initial + (i + 1) * size
                    if ff > r.final:
                        ff = r.final
                    g = GenomicRegion(chrom=r.chrom, initial=r.initial + i * size,
                                      final=ff, name=r.name, orientation=r.orientation)
                    z.add(g)
        return z

    def intersect_merge_pvalue(self, y, strandness=False):
        z = GenomicRegionSet(self.name)
        new_z = GenomicRegionSet(self.name)
        if len(self) == 0 or len(y) == 0:
            return z

        else:
            pull_regions = []
            for i, s in enumerate(self):
                # pull_regions.append([s])
                # print(s.data)
                for ss in y:
                    if s.overlap(ss, strandness=strandness):
                        # pull_regions[i].append(ss)
                        # print(ss.data)
                        pull_regions.append([s, ss])

            # print(len(pull_regions))
            for cc in pull_regions:
                # print(cc)
                # ps = [float(x.data) for x in cc]
                ps = [float(cc[0].data), float(cc[1].data)]
                ps.sort(reverse=True)
                # print(ps)
                # print(len(ps))

                p = - min([numpy.log10(len(ps)) - x - numpy.log10(j + 1) for j, x in enumerate(ps)])
                # print(p)
                z.add(GenomicRegion(chrom=cc[0].chrom,
                                    initial=max(cc[0].initial, cc[0].initial),
                                    final=min(cc[0].final, cc[0].final),
                                    orientation=cc[0].orientation,
                                    name=cc[0].name + ":" + cc[1].name,
                                    data=str(p)))
            z.sort()
            # print(z.sequences[1:10])
            for i, zz in enumerate(z):
                if i == 0:
                    previous_z = zz
                    continue
                else:
                    if zz.overlap(previous_z, strandness=True):
                        ps = [float(previous_z.data), float(zz.data)]
                        ps.sort(reverse=True)
                        p = - min([numpy.log10(len(ps)) - x - numpy.log10(j + 1) for j, x in enumerate(ps)])
                        previous_z = GenomicRegion(chrom=zz.chrom,
                                                   initial=min(previous_z.initial, zz.initial),
                                                   final=max(previous_z.final, zz.final),
                                                   orientation=previous_z.orientation,
                                                   name=previous_z.name,
                                                   data=str(p))

                    else:
                        new_z.add(previous_z)
                        previous_z = zz
            new_z.add(previous_z)
            return new_z
