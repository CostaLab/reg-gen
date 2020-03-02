###################################################################################################
# Setup Genomic Data Script
#   This script will fetch the necessary genomic data for all tools. Be aware that it will download
#   big data sets and it might be storage consuming. Also, this script will need an Internet
#   connection and may take several minutes to execute depending on the Internet speed.
###################################################################################################

# Import
import gzip
import argparse
from os import system, path, remove, mkdir
from sys import platform
import re


def download(url, prefix, output=None):
    if output:
        system("wget -c --tries=0 --read-timeout=30 " + url + " -O " + path.join(prefix, output))
    else:
        system("wget -c --tries=0 --read-timeout=30 " + url + " -P " + prefix)


# Initialising Options
parser = argparse.ArgumentParser()

# Parameter: RGT Data Location
parser.add_argument("--all", dest="all", action="store_true", default=False, help="Fetch all data sets.")
parser.add_argument("--hg19", dest="hg19", action="store_true", default=False, help="Fetch human genome files.")
parser.add_argument("--hg38", dest="hg38", action="store_true", default=False, help="Fetch human genome files.")
parser.add_argument("--mm9", dest="mm9", action="store_true", default=False, help="Fetch mouse files.")
parser.add_argument("--mm10", dest="mm10", action="store_true", default=False, help="Fetch mouse files.")
parser.add_argument("--zv9", dest="zv9", action="store_true", default=False, help="Fetch zebrafish files.")
parser.add_argument("--zv10", dest="zv10", action="store_true", default=False, help="Fetch zebrafish files.")

parser.add_argument("--hg19-genome-path", type=str, metavar="STRING", dest="hg19_genome_path", default=None,
                    help="Path to an already existing hg19 genome (all chromosomes in the same file).")
parser.add_argument("--hg38-genome-path", type=str, metavar="STRING", dest="hg38_genome_path", default=None,
                    help="Path to an already existing hg38 genome (all chromosomes in the same file).")
parser.add_argument("--mm9-genome-path", type=str, metavar="STRING", dest="mm9_genome_path", default=None,
                    help="Path to an already existing mm9 genome (all chromosomes in the same file).")
parser.add_argument("--mm10-genome-path", type=str, metavar="STRING", dest="mm10_genome_path", default=None,
                    help="Path to an already existing mm9 genome (all chromosomes in the same file).")
parser.add_argument("--zv9-genome-path", type=str, metavar="STRING", dest="zv9_genome_path", default=None,
                    help="Path to an already existing zv9 genome (all chromosomes in the same file).")
parser.add_argument("--zv10-genome-path", type=str, metavar="STRING", dest="zv10_genome_path", default=None,
                    help="Path to an already existing zv10 genome (all chromosomes in the same file).")

parser.add_argument("--hg19-gtf-path", type=str, metavar="STRING", dest="hg19_gtf_path", default=None,
                    help="Path to an already existing hg19 GTF file.")
parser.add_argument("--hg38-gtf-path", type=str, metavar="STRING", dest="hg38_gtf_path", default=None,
                    help="Path to an already existing hg38 GTF file.")
parser.add_argument("--mm9-gtf-path", type=str, metavar="STRING", dest="mm9_gtf_path", default=None,
                    help="Path to an already existing mm9 GTF file.")
parser.add_argument("--mm10-gtf-path", type=str, metavar="STRING", dest="mm10_gtf_path", default=None,
                    help="Path to an already existing mm10 GTF file.")
parser.add_argument("--zv9-gtf-path", type=str, metavar="STRING", dest="zv9_gtf_path", default=None,
                    help="Path to an already existing zv9 GTF file.")
parser.add_argument("--zv10-gtf-path", type=str, metavar="STRING", dest="zv10_gtf_path", default=None,
                    help="Path to an already existing zv10 GTF file.")

# Repeat masker
parser.add_argument("--hg38-rm", dest="hg38_rm", action="store_true", default=False,
                    help="Fetch repeat masker files for hg38 from Broad institute.")
parser.add_argument("--hg19-rm", dest="hg19_rm", action="store_true", default=False,
                    help="Fetch repeat masker files for hg19 from Broad institute.")
parser.add_argument("--mm9-rm", dest="mm9_rm", action="store_true", default=False,
                    help="Fetch repeat masker files for mm9 from Broad institute.")

args = parser.parse_args()

if not any([args.all, args.hg19, args.hg38, args.mm9, args.mm10, args.zv9, args.zv10]):
    parser.print_help()
    exit(1)

if args.all:
    args.hg19 = True
    args.hg38 = True
    args.mm9 = True
    args.mm10 = True
    args.zv9 = True
    args.zv10 = True

###################################################################################################
# Parameters
###################################################################################################

# Current rgt data path
curr_dir = path.dirname(path.realpath(__file__))

# Platform
supported_platforms = ["linux", "linux2", "darwin"]
if platform not in supported_platforms:
    print("ERROR: This package currently supports only unix-based systems (Linux and MAC OS X).")
    exit(0)

###################################################################################################
# Genomic Data HG19
###################################################################################################

gencode_url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/"

if args.hg19:

    output_location = path.join(curr_dir, "hg19")
    if not path.exists(output_location):
        mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location, "genome_hg19.fa")
    if args.hg19_genome_path:
        print("Creating symbolic link to HG19 genome")
        system("ln -s " + args.hg19_genome_path + " " + output_genome_file_name)
        print("OK")
    else:
        gen_root_url = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/"
        chr_list = ["chr" + str(e) for e in list(range(1, 23)) + ["X", "Y", "M"]]
        output_genome_file = open(output_genome_file_name, "w")
        to_remove = []
        for chr_name in chr_list:
            print("Downloading hg19 genome (" + chr_name + ")")
            gz_file_name = path.join(output_location, chr_name + ".fa.gz")
            download(gen_root_url + chr_name + ".fa.gz", output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write(gz_file.read().decode("utf-8"))
            gz_file.close()
            to_remove.append(gz_file_name)
            print("OK")
        output_genome_file.close()
        for gz_file_name in to_remove:
            remove(gz_file_name)

    # Fetching GTF
    gtf_output_file_name = path.join(output_location, "gencode.v19.annotation.gtf")
    if args.hg19_gtf_path:
        print("Creating symbolic link to HG19 GTF")
        system("ln -s " + args.hg19_gtf_path + " " + gtf_output_file_name)
        print("OK")
    else:
        gtf_url = gencode_url + "Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location, "gencode.v19.annotation.gtf.gz")
        if path.isfile(gtf_output_file_name_gz): remove(gtf_output_file_name_gz)
        print("Downloading hg19 GTF (gene annotation) from genode")
        download(gtf_url, output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name, "w")
        gtf_output_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print("OK")

###################################################################################################
# Genomic Data HG38
###################################################################################################

if args.hg38:

    output_location = path.join(curr_dir, "hg38")
    if not path.exists(output_location):
        mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location, "genome_hg38.fa")
    if args.hg38_genome_path:
        print("Creating symbolic link to HG38 genome")
        system("ln -s " + args.hg38_genome_path + " " + output_genome_file_name)
        print("OK")
    else:
        gen_root_url = gencode_url + "Gencode_human/release_24/GRCh38.primary_assembly.genome.fa.gz"
        # chr_list = ["chr"+str(e) for e in range(1,23)+["X","Y","M"]]
        output_genome_file = open(output_genome_file_name, "w")
        # for chr_name in chr_list:
        print("Downloading hg38 genome")
        gz_file_name = path.join(output_location, "GRCh38.primary_assembly.genome.fa.gz")
        # if(path.isfile(gz_file_name)): remove(gz_file_name)
        download(gen_root_url, output_location)
        gz_file = gzip.open(gz_file_name, 'rb')
        output_genome_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        remove(gz_file_name)
        print("OK")
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location, "gencode.v24.annotation.gtf")
    if args.hg38_gtf_path:
        print("Creating symbolic link to HG38 GTF")
        system("ln -s " + args.hg38_gtf_path + " " + gtf_output_file_name)
        print("OK")
    else:
        gtf_url = gencode_url + "Gencode_human/release_24/gencode.v24.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location, "gencode.v24.annotation.gtf.gz")
        if path.isfile(gtf_output_file_name_gz): remove(gtf_output_file_name_gz)
        print("Downloading hg19 GTF (gene annotation) from genode")
        download(gtf_url, output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name, "w")
        gtf_output_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print("OK")

    # Get Repeat Masker

###################################################################################################
# Genomic Data MM9
###################################################################################################

if args.mm9:

    output_location = path.join(curr_dir, "mm9")
    if not path.exists(output_location):
        mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location, "genome_mm9.fa")
    if args.mm9_genome_path:
        print("Creating symbolic link to MM9 genome")
        system("ln -s " + args.mm9_genome_path + " " + output_genome_file_name)
        print("OK")
    else:
        gen_root_url = "http://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/"
        chr_list = ["chr" + str(e) for e in list(range(1, 20)) + ["X", "Y", "M"]]
        output_genome_file = open(output_genome_file_name, "w")
        to_remove = []
        for chr_name in chr_list:
            print("Downloading MM9 genome (" + chr_name + ")")
            gz_file_name = path.join(output_location, chr_name + ".fa.gz")
            download(gen_root_url + chr_name + ".fa.gz", output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write(gz_file.read().decode("utf-8"))
            gz_file.close()
            to_remove.append(gz_file_name)
            print("OK")
        output_genome_file.close()
        for gz_file_name in to_remove:
            remove(gz_file_name)

    # Fetching GTF
    gtf_output_file_name = path.join(output_location, "gencode.vM1.annotation.gtf")
    if args.mm9_gtf_path:
        print("Creating symbolic link to MM9 GTF")
        system("ln -s " + args.mm9_gtf_path + " " + gtf_output_file_name)
        print("OK")
    else:

        gtf_url = gencode_url + "Gencode_mouse/release_M1/gencode.vM1.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location, "gencode.vM1.annotation.gtf.gz")
        if path.isfile(gtf_output_file_name_gz): remove(gtf_output_file_name_gz)
        print("Downloading MM9 GTF (gene annotation)")
        download(gtf_url, output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name, "w")
        gtf_output_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print("OK")

###################################################################################################
# Genomic Data MM10
###################################################################################################

if args.mm10:

    output_location = path.join(curr_dir, "mm10")
    if not path.exists(output_location):
        mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location, "genome_mm10.fa")
    if args.mm10_genome_path:
        print("Creating symbolic link to MM10 genome")
        system("ln -s " + args.mm10_genome_path + " " + output_genome_file_name)
        print("OK")
    else:
        gen_root_url = gencode_url + "Gencode_mouse/release_M11/GRCm38.primary_assembly.genome.fa.gz"
        output_genome_file = open(output_genome_file_name, "w")
        print("Downloading mm10 genome")
        gz_file_name = path.join(output_location, "GRCm38.primary_assembly.genome.fa.gz")
        download(gen_root_url, output_location)
        gz_file = gzip.open(gz_file_name, 'rb')
        output_genome_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        remove(gz_file_name)
        print("OK")
        output_genome_file.close()

    # Fetching GTF
    gtf_output_file_name = path.join(output_location, "gencode.vM20.annotation.gtf")
    if args.mm10_gtf_path:
        print("Creating symbolic link to MM10 GTF")
        system("ln -s " + args.mm10_gtf_path + " " + gtf_output_file_name)
        print("OK")
    else:
        gtf_url = gencode_url + "Gencode_mouse/release_M20/gencode.vM20.annotation.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location, "gencode.vM20.annotation.gtf.gz")
        if path.isfile(gtf_output_file_name_gz): remove(gtf_output_file_name_gz)
        print("Downloading MM10 GTF (gene annotation)")
        download(gtf_url, output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name, "w")
        gtf_output_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print("OK")

###################################################################################################
# Genomic Data ZV9
###################################################################################################

if args.zv9:

    output_location = path.join(curr_dir, "zv9")
    if not path.exists(output_location):
        mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location, "genome_zv9_ensembl_release_79.fa")
    if args.zv9_genome_path:
        print("Creating symbolic link to ZV9 genome")
        system("ln -s " + args.zv9_genome_path + " " + output_genome_file_name)
        print("OK")
    else:
        gen_root_url = "ftp://ftp.ensembl.org/pub/release-79/fasta/danio_rerio/dna/"
        chr_list = [str(e) for e in range(1, 25)] + ["MT"]
        output_genome_file = open(output_genome_file_name, "w")
        to_remove = []
        for chr_name in chr_list:
            print("Downloading ZV9 genome (chromosome " + chr_name + ")")
            gz_file_name = path.join(output_location, "Danio_rerio.Zv9.dna.chromosome." + chr_name + ".fa.gz")
            download(gen_root_url + "Danio_rerio.Zv9.dna.chromosome." + chr_name + ".fa.gz", output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write(re.sub("^>.*", ">chr" + chr_name, gz_file.read().decode("utf-8"), flags=re.MULTILINE))
            gz_file.close()
            to_remove.append(gz_file_name)
            print("OK")
        output_genome_file.close()
        for gz_file_name in to_remove:
            remove(gz_file_name)

    # Fetching GTF
    gtf_output_file_name = path.join(output_location, "Danio_rerio.Zv9.79.gtf")
    if args.zv9_gtf_path:
        print("Creating symbolic link to ZV9 GTF")
        system("ln -s " + args.zv9_gtf_path + " " + gtf_output_file_name)
        print("OK")
    else:
        gtf_url = "ftp://ftp.ensembl.org/pub/release-79/gtf/danio_rerio/Danio_rerio.Zv9.79.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location, "Danio_rerio.Zv9.79.gtf.gz")
        if path.isfile(gtf_output_file_name_gz): remove(gtf_output_file_name_gz)
        print("Downloading ZV9 GTF (gene annotation)")
        download(gtf_url, output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name, "w")
        gtf_output_file.write(re.sub('^([0-9a-zA-Z_\\\.]+\t)', "chr\\1", gz_file.read().decode("utf-8"), flags=re.MULTILINE))
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print("OK")

###################################################################################################
# Genomic Data ZV10
###################################################################################################

if args.zv10:

    output_location = path.join(curr_dir, "zv10")
    if not path.exists(output_location):
        mkdir(output_location)

    # Fetching genome
    output_genome_file_name = path.join(output_location, "genome_zv10_ensembl_release_84.fa")
    if args.zv10_genome_path:
        print("Creating symbolic link to ZV10 genome")
        system("ln -s " + args.zv10_genome_path + " " + output_genome_file_name)
        print("OK")
    else:
        gen_root_url = "ftp://ftp.ensembl.org/pub/release-84/fasta/danio_rerio/dna/"
        chr_list = [str(e) for e in range(1, 25)] + ["MT"]
        output_genome_file = open(output_genome_file_name, "w")
        to_remove = []
        for chr_name in chr_list:
            print("Downloading ZV10 genome (chromosome " + chr_name + ")")
            gz_file_name = path.join(output_location, "Danio_rerio.GRCz10.dna.chromosome." + chr_name + ".fa.gz")
            download(gen_root_url + "Danio_rerio.GRCz10.dna.chromosome." + chr_name + ".fa.gz", output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_genome_file.write(re.sub("^>.*", ">chr" + chr_name, gz_file.read().decode("utf-8"), flags=re.MULTILINE))
            gz_file.close()
            to_remove.append(gz_file_name)
            print("OK")
        output_genome_file.close()
        for gz_file_name in to_remove:
            remove(gz_file_name)

    # Fetching GTF
    gtf_output_file_name = path.join(output_location, "Danio_rerio.GRCz10.84.gtf")
    if args.zv10_gtf_path:
        print("Creating symbolic link to ZV10 GTF")
        system("ln -s " + args.zv10_gtf_path + " " + gtf_output_file_name)
        print("OK")
    else:
        gtf_url = "ftp://ftp.ensembl.org/pub/release-84/gtf/danio_rerio/Danio_rerio.GRCz10.84.gtf.gz"
        gtf_output_file_name_gz = path.join(output_location, "Danio_rerio.GRCz10.84.gtf.gz")
        if path.isfile(gtf_output_file_name_gz): remove(gtf_output_file_name_gz)
        print("Downloading ZV10 GTF (gene annotation)")
        download(gtf_url, output_location)
        gz_file = gzip.open(gtf_output_file_name_gz, 'rb')
        gtf_output_file = open(gtf_output_file_name, "w")
        gtf_output_file.write(re.sub("^([0-9a-zA-Z_\\\.]+\t)", "chr\\1", gz_file.read().decode("utf-8"), flags=re.MULTILINE))
        gz_file.close()
        remove(gtf_output_file_name_gz)
        gtf_output_file.close()
        print("OK")

###################################################################################################
# Repeat Maskers
###################################################################################################

if args.hg38_rm:
    root_url = "https://data.broadinstitute.org/igvdata/annotations/hg38/rmsk/"
    gen_root_url = ["LINE.bed", "LTR.bed", "DNA.bed", "Simple_repeat.bed", "Low_complexity.bed",
                    "Satellite.bed", "RC.bed", "Retroposon.bed", "RNA.bed", "rRNA.bed",
                    "scRNA.bed", "snRNA.bed", "srpRNA.bed", "tRNA.bed"]

    output_location = path.join(curr_dir, "hg38", "repeat_maskers")
    if not path.exists(output_location):
        mkdir(output_location)

    to_remove = []
    for masker in gen_root_url:
        output_file = open(path.join(output_location, masker), "w")
        gz_file_name = path.join(output_location, masker + ".gz")
        download(root_url + masker + ".gz", output_location)
        gz_file = gzip.open(gz_file_name, 'rb')
        output_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        print(path.join(output_location, masker) + "\t.....OK")
        output_file.close()
    for gz_file_name in to_remove:
        remove(gz_file_name)

if args.hg19_rm:
    root_url = "https://s3.amazonaws.com/igv.broadinstitute.org/annotations/hg19/variations/"
    gen_root_url = ["SINE.rmask.gz", "LINE.rmask.gz", "LTR.rmask.gz", "DNA.rmask.gz",
                    "Simple_repeat.rmask.gz", "Low_complexity.rmask.gz", "Satellite.rmask.gz", "RC.rmask.gz",
                    "RNA.rmask.gz", "Other.rmask.gz", "Unknown.rmask.gz"]

    output_location = path.join(curr_dir, "hg19", "repeat_maskers")
    if not path.exists(output_location):
        mkdir(output_location)

    to_remove = []
    for masker in gen_root_url:
        bedname = masker.split(".")[0] + ".bed"
        output_file = open(path.join(output_location, bedname), "w")
        gz_file_name = path.join(output_location, masker)
        download(root_url + masker, output_location)
        gz_file = gzip.open(gz_file_name, 'rb')
        output_file.write(gz_file.read().decode("utf-8"))
        gz_file.close()
        print(path.join(output_location, bedname) + "\t.....OK")
        output_file.close()
        to_remove.append(gz_file_name)
    for gz_file_name in to_remove:
        remove(gz_file_name)

if args.mm9_rm:
    gen_root_url = ["http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_SINE.bed",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_LINE.bed",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_LTR.bed",
                    "http://igv.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_DNA.bed.gz",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_Simple_repeat.bed",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_Low_complexity.bed",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_Satellite.bed",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_RNA.bed",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_Other.bed",
                    "http://igvdata.broadinstitute.org/annotations/mm9/variation/repeat_masker/mm9_repmask_Unknown.bed"]

    output_location = path.join(curr_dir, "mm9", "repeat_maskers")
    if not path.exists(output_location):
        mkdir(output_location)

    to_remove = []
    for masker in gen_root_url:
        bedname = path.basename(masker).replace("mm9_repmask_", "").replace(".gz", "")

        if ".gz" in masker:
            output_file = open(path.join(output_location, bedname), "w")
            gz_file_name = path.join(output_location, path.basename(masker))
            download(masker, output_location)
            gz_file = gzip.open(gz_file_name, 'rb')
            output_file.write(gz_file.read().decode("utf-8"))
            gz_file.close()
            output_file.close()
            to_remove.append(gz_file_name)
        else:
            download(masker, output_location, output=bedname)

        print(path.join(output_location, bedname) + "\t.....OK")
    for gz_file_name in to_remove:
        remove(gz_file_name)
