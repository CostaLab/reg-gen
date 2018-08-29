# Python Libraries
from __future__ import print_function
from __future__ import division
import os
from ..Util import GenomeData
from ..tdf.triplexTools import save_sequence, run_triplexator, silentremove
from ..tdf.RNADNABindingSet import RNADNABindingSet
# from rgt.GenomicRegionSet import GenomicRegionSet


def random_each(input):
    """Return the counts of DNA Binding sites with randomization
    For multiprocessing.
    Input contains:
    0       1               2                3     4              5          6
    str(i), self.rna_fasta, self.dna_region, temp, self.organism, self.rbss, str(marks.count(i)),
    number, rna,            region,          temp, organism,      rbss,      number of mark

    7  8  9  10  11  12  13  14  15          16                 17
    l, e, c, fr, fm, of, mf, rm, filter_bed, self.genome_path,  par
    """
    import sys
    # Filter BED file
    if input[15]:
        random = input[2].random_regions(organism=input[4], multiply_factor=1,
                                         overlap_result=True, overlap_input=True,
                                         chrom_X=True, chrom_M=False, filter_path=input[15])
    else:
        random = input[2].random_regions(organism=input[4], multiply_factor=1,
                                         overlap_result=True, overlap_input=True,
                                         chrom_X=True, chrom_M=False)
    # Generate FASTA
    save_sequence(dir=input[3], filename="random_" + input[0] + ".fa",
                  regions=random, genome_path=input[16])
    # Triplexator
    run_triplexator(ss=input[1], ds=os.path.join(input[3], "random_" + input[0] + ".fa"),
                    output=os.path.join(input[3], "random_" + input[0] + ".tpx"),
                    l=int(input[7]), e=int(input[8]), c=input[9], fr=input[10],
                    fm=input[11], of=input[12], mf=input[13], rm=input[14],
                    par=input[17])

    # Read txp
    tpx = RNADNABindingSet("random")
    tpx.read_tpx(os.path.join(input[3], "random_" + input[0] + ".tpx"),
                 dna_fine_posi=False, seq=True)
    tpx.remove_duplicates()


    tpx.merge_rbs(rbss=input[5], rm_duplicate=True)

    tpxf = RNADNABindingSet("random")
    tpxf.read_tpx(os.path.join(input[3], "random_" + input[0] + ".tpx"),
                  dna_fine_posi=True, seq=True)
    tpxf.merge_rbs(rbss=input[5], rm_duplicate=True)
    sys.stdout.flush()
    print("".join(["="] * int(input[6])), end="")
    distances = tpxf.distance_distribution()
    silentremove(os.path.join(input[3], "random_" + input[0] + ".fa"))
    silentremove(os.path.join(input[3], "random_" + input[0] + ".tpx"))

    return [[len(tr) for tr in tpx.merged_dict.values()], [len(dbss) for dbss in tpxf.merged_dict.values()], distances]


class Triplexes(object):

    def __init__(self, organism, pars, l=None, e=None):
        self.genome = GenomeData(organism=organism)
        self.organism = organism
        if not l:
            self.l = pars.l
        else:
            self.l = l
        if not e:
            self.e = pars.e
        else:
            self.e = e
        self.c = pars.c
        self.fr = pars.fr
        self.fm = pars.fm
        self.of = pars.of
        self.mf = pars.mf
        self.rm = pars.rm
        self.pars = pars
        self.outdir = pars.o

    def search_triplex(self, target_regions, prefix, remove_temp=False):
        # print("    \tRunning Triplexator...")
        rna_fasta = os.path.join(self.outdir, "rna_temp.fa")
        dna_fasta = os.path.join(self.outdir, prefix+".fa")
        tpx_file = os.path.join(self.outdir, prefix+".tpx")
        # Target
        save_sequence(dir=self.outdir, filename=dna_fasta,
                      regions=target_regions, genome_path=self.genome.get_genome())

        run_triplexator(ss=rna_fasta, ds=dna_fasta, output=tpx_file,
                        l=self.l, e=self.e, c=self.c, fr=self.fr, fm=self.fm,
                        of=self.of, mf=self.mf, rm=self.rm, par=self.pars.par)

        if remove_temp:
            os.remove(dna_fasta)

        return tpx_file

    def find_autobinding(self, rbss):
        rna_fasta = os.path.join(self.pars.o, "rna_temp.fa")
        run_triplexator(ss=None, ds=None, autobinding=rna_fasta,
                        output=os.path.join(self.pars.o, "autobinding.tpx"),
                        l=self.l, e=self.e, c=self.c, fr=self.fr, fm=self.fm, of=self.of, mf=self.mf,
                        rm=self.rm, par="abo_0")
        self.autobinding = RNADNABindingSet("autobinding")
        self.autobinding.read_tpx(filename=os.path.join(self.pars.o, "autobinding.tpx"), dna_fine_posi=True, seq=True)

        self.autobinding.merge_rbs(rbss=rbss, rm_duplicate=False)

    def get_tpx(self, rna_fasta_file, target_regions, dna_fine_posi, prefix="", remove_temp=False,
                autobinding=False):
        """Given a GenomicRegionSet to run Triplexator and return the RNADNABindingSet"""
        # Generate FASTA
        save_sequence(dir=self.outdir, filename="targets_" + prefix + ".fa",
                      regions=target_regions, genome_path=self.genome.get_genome())
        # Triplexator
        run_triplexator(ss=rna_fasta_file, ds=os.path.join(self.outdir, "targets_" + prefix + ".fa"),
                        output=os.path.join(self.outdir, "targets_" + prefix + ".tpx"),
                        l=self.l, e=self.e, c=self.c, fr=self.fr, fm=self.fm, of=self.of,
                        mf=self.mf, rm=self.rm, par=self.pars.par)
        # Autobinding
        if autobinding:
            run_triplexator(ss=rna_fasta_file, ds=os.path.join(self.outdir, "targets_" + prefix + ".fa"),
                            output=os.path.join(self.outdir, "autobinding_" + prefix + ".txp"),
                            l=self.l, e=self.e, c=self.c, fr=self.fr, fm=self.fm, of=self.of,
                            mf=self.mf, rm=self.rm, par=self.pars.par + "_auto-binding-file")
        # Read txp
        tpx = RNADNABindingSet("targets")
        tpx.read_tpx(os.path.join(self.outdir, "targets_" + prefix + ".tpx"), dna_fine_posi=dna_fine_posi, seq=True)
        tpx.remove_duplicates()
        if remove_temp:
            os.remove(os.path.join(self.outdir, "targets_" + prefix + ".fa"))
            os.remove(os.path.join(self.outdir, "targets_" + prefix + ".tpx"))
        return tpx

    def merging_combined_motifs(self, tpx1_file, tpx2_file, output_file, name= "rdb"):
        """Load two tpx files and merge RBSs and DBSs according to the rules for combined motifs and save into a tpx file"""
        tpx1 = RNADNABindingSet(name="tpx1")
        tpx1.read_tpx(tpx1_file, dna_fine_posi=True, seq=True)
        tpx2 = RNADNABindingSet(name="tpx2")
        tpx2.read_tpx(tpx2_file, dna_fine_posi=True, seq=True)
        tpx2_merged = tpx2.combine_motif(name=name)
        tpx_12_merged = RNADNABindingSet("merged tpx")
        for rdb in tpx1:
            tpx_12_merged.add(rdb)
        for rdb in tpx2_merged:
            tpx_12_merged.add(rdb)
        tpx_12_merged.write_txp(output_file, match=True)
