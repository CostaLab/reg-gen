
# Python 3 compatibility


# Python
import unittest
import os

# Internal
from rgt.GenomicRegion import GenomicRegion
from rgt.Util import GenomeData
from rgt.motifanalysis.Match import match_multiple
from rgt.MotifSet import MotifSet

# External
from MOODS import tools, scan
from pysam.libcfaidx import Fastafile


class MatchTest(unittest.TestCase):
    def setUp(self):
        # the genome must be available
        # TODO: we could make this test pure by manually using the sequence corresponding to the input region
        self.genome_data = GenomeData("hg19")
        self.genome_data.genome = os.path.join(os.path.dirname(__file__), "hg19_chr1_710000_715000.fa")
        self.genome_file = Fastafile(self.genome_data.get_genome())

    def test_match_multiple(self):

        ms = MotifSet(preload_motifs="default")
        ms = ms.filter({'database': ["jaspar_vertebrates"], 'name': ["MA0139.1.CTCF"]}, search="inexact")

        self.assertEqual(len(ms), 1)

        motif = ms.get_motif_list(1, 0.0001)[0]

        scanner = scan.Scanner(7)

        pssm_list, thresholds = [], []

        thresholds.append(motif.threshold)
        thresholds.append(motif.threshold)
        pssm_list.append(motif.pssm)
        pssm_list.append(motif.pssm_rc)

        bg = tools.flat_bg(4)
        scanner.set_motifs(pssm_list, bg, thresholds)

        genomic_region = GenomicRegion("chr1", 0, 5022)

        # Reading sequence associated to genomic_region
        sequence = str(self.genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

        grs = match_multiple(scanner, [motif], sequence, genomic_region)

        self.assertSequenceEqual(grs.sequences,
                                 [GenomicRegion("chr1", 4270, 4289, name="MA0139.1.CTCF", orientation="+"),
                                  GenomicRegion("chr1", 4180, 4199, name="MA0139.1.CTCF", orientation="-")])
