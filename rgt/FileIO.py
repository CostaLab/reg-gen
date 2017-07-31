
# Python
from __future__ import print_function
from abc import ABCMeta, abstractmethod
import os

# Internal
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.SequenceSet import SequenceSet, Sequence, SequenceType


class FileIO:
    """
    Base (abstract) class defining the contract of the file backend for GenomicRegionSets.
    Should not, obviously, be instantiated.
    """

    __metaclass__ = ABCMeta

    def __init__(self, filename):
        self.filename = filename

    @abstractmethod
    def read(self, data):
        """
        Concrete classes should return the (modified) data object. The input data is considered stale.
        """
        return data

    @abstractmethod
    def write(self, data, mode):
        pass


class BedFile(FileIO):
    """
    Each row maps to a GenomicRegion.

    Note: Chrom (1), start (2), end (2), name (4) and orientation (6) is used for GenomicRegion.
          All other columns (5, 7, 8, ...) are put to the data attribute of the GenomicRegion.
          The numbers in parentheses are the columns of the BED format.
    """

    def __init__(self, filename):
        super(BedFile, self).__init__(filename)

    def read(self, grs):

        if not isinstance(grs, GenomicRegionSet):
            raise TypeError("Must be a GenomicRegionSet")

        with open(self.filename) as f:
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
                        raise Exception("zero-length region: " + grs.chrom + "," + str(grs.initial) + "," + str(grs.final))
                    g = GenomicRegion(chrom, start, end, name, orientation, data)

                    grs.add(g)
                except:
                    if not line:
                        continue
                    else:
                        error_line += 1
                        if error_line > 2:
                            # Skip the first error line which contains the track information
                            print("Error at line", line, self.filename)
            grs.sort()

        return grs

    def write(self, grs, mode="w"):

        with open(self.filename, mode) as f:
            for gr in grs:
                print(gr, file=f)


class Bed12File(FileIO):
    """
    Bed file with "block information", eg exons.
    """

    def __init__(self, filename):
        super(Bed12File, self).__init__(filename)

    def read(self, grs):
        if not isinstance(grs, GenomicRegionSet):
            raise TypeError("Must be a GenomicRegionSet")

        with open(self.filename) as f:
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
                        raise Exception("zero-length region: " + grs.chrom + "," + str(grs.initial) + "," + str(grs.final))
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
                            print("Error at line", line, self.filename)
            grs.sort()

        return grs

    def write(self, grs, mode="w"):

        if not isinstance(grs, GenomicRegionSet):
            raise TypeError("Must be a GenomicRegionSet")

        with open(self.filename, mode) as f:

            blocks = {}
            for gr in grs:
                try:
                    blocks[gr.name].add(gr)
                except:
                    blocks[gr.name] = GenomicRegionSet(gr.name)
                    blocks[gr.name].add(gr)

            for name in blocks.keys():
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


class BedGraphFile(FileIO):
    """
    BedGraph format, read-only.
    """

    def __init__(self, filename):
        super(BedGraphFile, self).__init__(filename)

    def read(self, grs):

        if not isinstance(grs, GenomicRegionSet):
            raise TypeError("Must be a GenomicRegionSet")

        with open(self.filename) as f:
            for line in f:
                try:
                    line = line.strip("\n")
                    line = line.split("\t")
                    assert len(line) == 4

                    chrom, start, end, data = line[0], int(line[1]), int(line[2]), str(line[3])

                    grs.add(GenomicRegion(chrom=chrom, initial=start, final=end, data=data))
                except:
                    print("Error at line", line, self.filename)

            grs.sort()

        return grs

    def write(self, grs, mode):
        raise NotImplementedError


class FastaFile(FileIO):
    """
    FIXME: do we actually need this?
    """

    def __init__(self, filename):
        super(FastaFile, self).__init__(filename)

    def read(self, grs):

        if not isinstance(grs, GenomicRegionSet):
            raise TypeError("Must be a GenomicRegionSet")

        # Parse each chromosome and fetch the defined region in this chromosome
        chroms = list(set(grs.get_chrom()))

        chro_files = [x.split(".")[0] for x in os.listdir(self.filename)]

        for ch in chroms:
            if ch not in chro_files:
                print(" *** There is no genome FASTA file for: "+ch)

            # Read genome in FASTA according to the given chromosome
            ch_seq = SequenceSet(name=ch, seq_type=SequenceType.DNA)
            try:
                ch_seq.read_fasta(os.path.join(self.filename, ch+".fa"))
            except:
                continue

            # Regions in given chromosome
            beds = grs.any_chrom(chrom=ch)

            for s in beds:
                seq = ch_seq[0].seq[s.initial:s.final]

                try:
                    strand = s.strand
                except:
                    strand = "+"

                s.sequence = (Sequence(seq=seq, name=s.__repr__(), strand=strand))

        return grs

    def write(self, grs, mode="w"):

        with open(self.filename, mode) as f:
            for r in grs:
                try:
                    f.write(">" + r.chrom + ":" + str(r.initial) + "-" + str(r.final) + "-" +
                            str(r.orientation) + "\n" + r.sequence.seq + "\n")
                except:
                    pass
