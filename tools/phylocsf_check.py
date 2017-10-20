from __future__ import print_function

from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.Util import GenomeData
import argparse
import subprocess 
import os
import sys
sys.path.insert(1, "/home/joseph/Apps/alignio-maf")
from Bio import AlignIO
from Bio.AlignIO import MultipleSeqAlignment
#import Bio.AlignIO import MafIO
from Bio.AlignIO.MafIO import MafIndex

def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

##################################################################################
parser = argparse.ArgumentParser(description='Check the coding potential by PhyloCSF', 
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', metavar='  ', type=str, help="Input BED file")
parser.add_argument('-o', metavar='  ', type=str, help="Output BED file with the coding-potential score")
parser.add_argument('-organism', metavar='  ', type=str, help="Define the organism")
parser.add_argument('-rmcoding', metavar='  ', type=float, help="Define the cutoff to remove the entries with coding potential")
parser.add_argument('-mafdir', metavar='  ', type=str, help="Define the directory to MAF files")
# python /projects/reg-gen/tools/phylocsf_check.py -i
args = parser.parse_args()

bed = GenomicRegionSet("input")
bed.read(args.i)
num = len(bed)

organisms = { "hg18": "Human",
              "panTro2": "Chimp",
              "rheMac2": "Rhesus",
              "tarSyr1": "Tarsier",
              "micMur1": "Mouse_lemur",
              "otoGar1": "Bushbaby",
              "tupBel1": "Shrew",
              "mm9": "Mouse",
              "rn4": "Rat",
              "dipOrd1": "Kangaroo_Rat",
              "cavPor2": "Guinea_Pig",
              "speTri1": "Squirrel",
              "oryCun1": "Rabbit",
              "ochPri2": "Pika",
              "vicPac1": "Alpaca",
              "turTru1": "Dolphin",
              "bosTau3": "Cow",
              "bosTau4": "Cow",
              "equCab1": "Horse",
              "equCab2": "Horse",
              "felCat3": "Cat",
              "canFam2": "Dog",
              "myoLuc1": "Microbat",
              "pteVam1": "Megabat",
              "eriEur1": "Hedgehog",
              "sorAra1": "Shrew",
              "loxAfr1": "Elephant",
              "loxAfr2": "Elephant",
              "proCap1": "Rock_Hyrax",
              "echTel1": "Tenrec",
              "dasNov1": "Armadillo",
              "dasNov2": "Armadillo",
              "choHof1": "Sloth",
              #"ponAbe2": "Orangutan",
              #"calJac1": "Marmoset",
              #"monDom4": "Opossum",
              #"fr2": "Fugu",
              #"ornAna1": "Platypus",
              #"danRer5": "Zebrafish",
              #"galGal3": "Chicken",
              #"gasAcu1": "Stickleback",
              #"xenTro2": "Frog",
              #"anoCar1": "Lizard",
              #"tetNig1": "Tetraodon",
              #"oryLat1": "Medaka"
               }

for i, rg in enumerate(bed):
    print(str(i+1)+"/"+str(num)+":\t"+rg.name+"\t"+rg.toString())
    try:
        idx = MafIndex(os.path.join(args.mafdir, rg.chrom+".mafindex"),
                       os.path.join(args.mafdir, rg.chrom+".maf"), 
                       args.organism+"."+rg.chrom)
    except:
        print("Generating index file...")
        os.remove(os.path.join(args.mafdir, rg.chrom+".mafindex"))
        idx = MafIndex(os.path.join(args.mafdir, rg.chrom+".mafindex"),
                       os.path.join(args.mafdir, rg.chrom+".maf"), 
                       args.organism+"."+rg.chrom)
    exon_s = []
    exon_e = []
    if rg.data:
        data = rg.data.split("\t")
        if len(data) == 7 and int(data[4]) > 1:
            z = rg.extract_blocks()
            for g in z:
                exon_s.append(g.initial)
                exon_e.append(g.final)
        else:
            exon_s = [rg.initial]
            exon_e = [rg.final]

    else: 
        exon_s = [rg.initial]
        exon_e = [rg.final]

    if rg.orientation == "+":
        strand = "+1"
    else:
        strand = "-1"     
    multiple_alignment = idx.get_spliced(exon_s,
                                         exon_e,
                                         strand = strand)
    
    for j, seqrec in enumerate(multiple_alignment):
        if seqrec.id.startswith(args.organism):
            if j == 0: break
            else:
                multiple_alignment._records[0], multiple_alignment._records[j] = multiple_alignment._records[j], multiple_alignment._records[0]
    
    #print(dir(multiple_alignment))
    #sys.exit(0)

    seqs = []
    for seqrec in multiple_alignment:
        #print(dir(seqrec))
        try:
            name_id = seqrec.id.partition(" ")[0].partition(".")[0]
            #print(name_id)
            seqrec.id = organisms[ name_id ]
            if seqrec.id in organisms.values():
                seqs.append(seqrec)
        except:
            continue

    new_alignment = MultipleSeqAlignment(records=seqs)
    #print(len(new_alignment))
            
    AlignIO.write(new_alignment, "mm9_"+rg.name+".fa", "fasta")

    process = subprocess.Popen(["/home/joseph/Apps/PhyloCSF/PhyloCSF", 
                                "29mammals", 
                                "mm9_"+rg.name+".fa", 
                                "--removeRefGaps",
                                "--strategy=omega",
                                "--orf=StopStop3",
                                "--minCodons=25",
                                "--frames=3"], 
                                stdout=subprocess.PIPE)
    out, err = process.communicate()
    print(out)
    #print(out.split("\t")[2])
    #print(out.split("\t")[3])
    #print(out.split("\t")[4])

    data = rg.data.split("\t")
    score = out.split("\t")[2]
    rg.data = "\t".join([score] + data[1:])

bed.write(args.o)

    
# 29/9/2015
# python /projects/reg-gen/tools/phylocsf_check.py -i /projects/ig440396_dendriticcells/exp/RNASeq/expression/isofroms/deseq/new_bed/all_TCONs.bed -o all_TCONS_phyloCSF.bed -organism mm9 -mafdir /data/genome/mm9/multiz30way/maf/
