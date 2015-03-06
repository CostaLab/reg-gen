from __future__ import print_function
import os
from rgt.GenomicRegion import GenomicRegion
from rgt.GenomicRegionSet import GenomicRegionSet

# Test
from rgt.Util import SequenceType
####################################################################################
####################################################################################

class BindingSite(GenomicRegion):
    """Describes a binding region on DNA or RNA including the information regarding to this region.

    Authors: Joseph Kuo
    """

    def __init__(self, chrom, initial, final, name=None, score=None, errors_bp=None, motif=None, 
                 strand=None, orientation=None, guanine_rate=None, seq=None):
        """Initialize
        
        name             The name of this binding site (Default: None)
        seq_type         DNA or RNA
        chrm             Define the chromosome for DNA; for RNA its default is "RNA"
        initial          Binding start position
        final            Binding end position
        score            Score of the binding pattern (Default: None)
        errors_bp        Error base pair in this binding (Default: None)
        motif            The motif for this binding (Default: None)
        strand           The strand of DNA (+ or -) (Default: None)
        orientation      Parallel or antiparallel (Default: None)
        guanine_rate     (Default: None)
        seq              Sequence of this region with ATCG as letters

        """
        GenomicRegion.__init__(self, chrom=chrom, initial=initial, final=final)
        
        self.name = name                      # RNA name
        self.score = score                    # Score for pattern matching
        self.errors_bp = errors_bp                  
        self.motif = motif
        #self.strand = strand
        self.orientation = orientation
        self.seq = seq                        # An object (Sequence) not just a string
        if seq:
            self.guanine_rate = "{0:.2f}".format(float(seq.seq.count("G"))/len(seq))

    def __str__(self):
        """Give informal string representation"""
        infos = [ self.name, self.chrom, self.initial, self.final, self.score, self.errors_bp, 
                  self.motif, self.orientation, self.seq ]
        return '-'.join( [str(x) for x in infos if x] )
        
    def __repr__(self):
        """Return official representation of GenomicRegion"""
        infos = [ self.name, self.chrom, self.initial, self.final, self.score, self.errors_bp, 
                  self.motif, self.orientation, self.seq ]
        return ','.join( [str(x) for x in infos if x] ) 
        
    def __len__(self):
        """Return the length of the binding site """
        return self.final - self.initial

    def __eq__(self, other):
        return (self.initial, self.final) == (other.initial, other.final)

    def __hash__(self):
        return hash(tuple([self.chrom, self.initial, self.final]))
                        
    def str_rna(self, pa=True):
        if pa:
            return "{0}-{1}-{2}".format(self.initial, self.final, self.orientation)
        else:
            return "{0}-{1}".format(self.initial, self.final)
            
            
            
            

####################################################################################
####################################################################################

class BindingSiteSet(GenomicRegionSet):
    """Represent a collection of RNABinding with some functions

    Authors: Joseph Kuo
    """

    def __init__(self, name):
        """Initialize"""
        GenomicRegionSet.__init__(self, name = name)

    def sort(self):
        """Sort Elements by criteria defined by a GenomicRegion"""
        self.sequences.sort(cmp = GenomicRegion.__cmp__)
        self.sorted = True
    
    def get_bs(self, orientation):
        """Get the Binding Sites with the given orientation"""
        output = BindingSiteSet(self.name+":"+orientation)
        for bs in self.sequences:
            if bs.orientation == orientation:
                output.add(bs)
        return output

    def write_rbs(self, filename):
        """Write the information into a file with .rbs """
        d = os.path.dirname(filename)
        if not os.path.exists(d):
            os.makedirs(d)

        with open(filename, "w") as f:
            print("# Sequence-ID\tStart\tEnd\tScore\tMotif\tError-rate\tGuanine-rate\tSequence",file=f)
            # Sequence-ID   Start   End Score   Motif   Error-rate  Errors  
            # Guanine-rate    Duplicates  TFO Duplicate locations

            for bs in self.sequences:
                err_rate = 1 - bs.score/(bs.final - bs.initial)
                print("\t".join([bs.name, str(bs.initial), str(bs.final), str(bs.score), 
                                 bs.motif, "{0:.2f}".format(err_rate), bs.guanine_rate, bs.seq.seq]), file=f)

    def write_dbs(self, filename):
        """Write the information into a file with .dbs """
        d = os.path.dirname(filename)
        if not os.path.exists(d):
            os.makedirs(d)

        with open(filename, "w") as f:
            print("# Chromosome\tStart\tEnd\tScore\tError-rate\tGuanine-rate\tSequence",file=f)
            # Sequence-ID   Start   End Score   Motif   Error-rate  Errors  
            # Guanine-rate    Duplicates  TFO Duplicate locations

            for bs in self.sequences:
                print(bs)
                err_rate = 1 - bs.score/(bs.final - bs.initial)
                print("\t".join([bs.chrom, str(bs.initial), str(bs.final), str(bs.score), 
                                 "{0:.2f}".format(err_rate), bs.guanine_rate, bs.seq.seq]), file=f)

    def concatenate(self, another_BindingSiteSet):
        """Concatenate another RNABindingSet without sorting"""
        self.sequences = self.sequences + another_BindingSiteSet.sequences

    def count_rbs_position(self, bp):
        count = 0
        for r in self.sequences:
            if r.initial <= bp and r.final > bp:
                count += 1
        return count



####################################################################################
#   Test
####################################################################################

if __name__ == '__main__':
    a = BindingSiteSet(name="a")
    a.add(BindingSite(SequenceType.RNA, "a",1,5))
    a.add(BindingSite(SequenceType.RNA, "a",10,15))
    
    b = BindingSiteSet(name="b")
    b.add(BindingSite(SequenceType.RNA, "b",4,8))
    
    print(len(a))
    print(len(b))
    print(a.sequences)
    print(b.sequences)
    
    a.subtract(b)
    a.merge()
    print(a.sequences)
    print(b.sequences)
    