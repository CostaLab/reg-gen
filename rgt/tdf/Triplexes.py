import os
from rgt.Util import GenomeData
from rgt.tdf.triplexTools import save_sequence, run_triplexator
from rgt.tdf.RNADNABindingSet import RNADNABindingSet
from rgt.GenomicRegionSet import GenomicRegionSet

class Triplexes(object):

    def __init__(self, organism, pars):
        self.genome = GenomeData(organism=organism)
        self.l = pars.l
        self.e = pars.e
        self.c = pars.c
        self.fr = pars.fr
        self.fm = pars.fm
        self.of = pars.of
        self.mf = pars.mf
        self.pars = pars.pars
        self.outdir = pars.o

    def search_triplex(self, target_regions, prefix, remove_temp=False):
        print("    \tRunning Triplexator...")
        rna_fasta = os.path.join(self.outdir, "rna_temp.fa")
        dna_fasta = os.path.join(self.outdir, prefix+".fa")
        tpx_file = os.path.join(self.outdir, prefix+".tpx")
        # Target
        save_sequence(dir=self.outdir, filename=dna_fasta,
                      regions=target_regions, genome_path=self.genome.get_genome())

        run_triplexator(ss=rna_fasta, ds=dna_fasta, output=tpx_file,
                        l=self.l, e=self.e, c=self.c, fr=self.fr, fm=self.fm,
                        of=self.of, mf=self.mf, par=self.pars.par)

        if remove_temp:
            os.remove(dna_fasta)

        return tpx_file

    def autobinding(self, rbss):
        rna_fasta = os.path.join(self.pars.o, "rna_temp.fa")
        run_triplexator(ss=None, ds=None, autobinding=rna_fasta,
                        output=os.path.join(self.pars.o, "autobinding.tpx"),
                        l=self.l, e=self.e, c=self.c, fr=self.fr, fm=self.fm, of=self.of, mf=self.mf,
                        par="abo_0")
        self.autobinding = RNADNABindingSet("autobinding")
        self.autobinding.read_tpx(filename=os.path.join(self.pars.o, "autobinding.tpx"), dna_fine_posi=True, seq=True)

        self.autobinding.merge_rbs(rbss=self.rbss, rm_duplicate=False)
        # self.autobinding.motif_statistics()
        # Saving autobinding dbs in BED
