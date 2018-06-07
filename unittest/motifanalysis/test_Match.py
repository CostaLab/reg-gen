
# Python 3 compatibility
from __future__ import print_function

# Python
import unittest
import os

# Internal
from rgt.GenomicRegion import GenomicRegion
from rgt.Util import GenomeData
from rgt.motifanalysis.Match import match_multiple
from rgt.motifanalysis.Motif import Motif

# External
from MOODS import tools, scan
from pysam.libcfaidx import Fastafile


class MatchTest(unittest.TestCase):
    def setUp(self):
        # the genome must be available
        # TODO: we could make this test pure by manually using the sequence corresponding to the input region
        self.genome_data = GenomeData("hg19")
        self.genome_file = Fastafile(self.genome_data.get_genome())

    def test_match_multiple(self):
        dirname = os.path.dirname(__file__)
        jasp_dir = "../../data/motifs/jaspar_vertebrates/"

        scanner = scan.Scanner(7)

        pssm_list = []
        thresholds = []

        motif = Motif(os.path.join(dirname, jasp_dir, "MA0139.1.CTCF.pwm"), 1, 0.0001, None)

        thresholds.append(motif.threshold)
        thresholds.append(motif.threshold_rc)
        pssm_list.append(motif.pssm)
        pssm_list.append(motif.pssm_rc)

        bg = tools.flat_bg(4)
        scanner.set_motifs(pssm_list, bg, thresholds)

        genomic_region = GenomicRegion("chr1", 710000, 715000)

        # Reading sequence associated to genomic_region
        sequence = str(self.genome_file.fetch(genomic_region.chrom, genomic_region.initial, genomic_region.final))

        grs = match_multiple(scanner, [motif], sequence, genomic_region)

        self.assertSequenceEqual(grs.sequences,
                                 [GenomicRegion("chr1", 714270, 714289, name="MA0139.1.CTCF", orientation="+"),
                                  GenomicRegion("chr1", 714180, 714199, name="MA0139.1.CTCF", orientation="-")])
