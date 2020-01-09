import os
import shutil
import subprocess
from itertools import product
from math import floor

from Bio import motifs
# External
from pysam import Samfile, Fastafile

from rgt.GenomicRegionSet import GenomicRegionSet
from rgt.Util import AuxiliaryFunctions, GenomeData, HmmData


def estimation_args(parser):
    # Input Options
    parser.add_argument("--organism", type=str, metavar="STRING", default="hg19",
                        help="Organism considered on the analysis. Must have been setup in the RGTDATA folder. "
                             "Common choices are hg19, hg38. mm9, and mm10. DEFAULT: hg19")
    parser.add_argument("--bias-type", type=str, metavar="STRING", default="VOM",
                        help="The methods that used to estimate the bias table "
                             "Available options are: 'KMER', 'PWM' and 'VOM'. DEFAULT: VOM")
    parser.add_argument("--reads-file", type=str, metavar="FILE", default=None,
                        help="The BAM file containing aligned reads. DEFAULT: None")
    parser.add_argument("--regions-file", type=str, metavar="FILE", default=None,
                        help="The BED file containing regions to estimate the bias. DEFAULT: None")
    parser.add_argument("--downstream-ext", type=int, metavar="INT", default=1)
    parser.add_argument("--upstream-ext", type=int, metavar="INT", default=0)
    parser.add_argument("--forward-shift", type=int, metavar="INT", default=4)
    parser.add_argument("--reverse-shift", type=int, metavar="INT", default=-4)
    parser.add_argument("--k-nb", type=int, metavar="INT", default=8,
                        help="Size of k-mer for bias estimation. DEFAULT: 8")

    # Output Options
    parser.add_argument("--output-location", type=str, metavar="PATH", default=os.getcwd(),
                        help="Path where the output bias table files will be written. DEFAULT: current directory")
    parser.add_argument("--output-prefix", type=str, metavar="STRING", default="Bias",
                        help="The prefix for results files. DEFAULT: Bias")


def estimation_run(args):
    if args.bias_type == "KMER":
        estimate_bias_kmer(args)
    if args.bias_type == "PWM":
        estimate_bias_pwm(args)
    if args.bias_type == "VOM":
        estimate_bias_vom(args)


def estimate_bias_kmer(args):
    # Parameters
    maxDuplicates = 100
    pseudocount = 1.0

    # Initializing bam and fasta
    bamFile = Samfile(args.reads_file, "rb")
    genome_data = GenomeData(args.organism)
    fastaFile = Fastafile(genome_data.get_genome())
    regions = GenomicRegionSet("regions")
    regions.read(args.regions_file)

    # Initializing dictionaries
    obsDictF = dict()
    obsDictR = dict()
    expDictF = dict()
    expDictR = dict()

    ct_reads_r = 0
    ct_reads_f = 0
    ct_kmers = 0

    # Iterating on HS regions
    for region in regions:

        # Initialization
        prevPos = -1
        trueCounter = 0

        # Evaluating observed frequencies ####################################
        # Fetching reads
        for r in bamFile.fetch(region.chrom, region.initial, region.final):

            # Calculating positions
            if not r.is_reverse:
                cut_site = r.pos + args.forward_shift - 1
                p1 = cut_site - int(floor(args.k_nb / 2))
            else:
                cut_site = r.aend + args.reverse_shift + 1
                p1 = cut_site - int(floor(args.k_nb / 2))
            p2 = p1 + args.k_nb

            # Verifying PCR artifacts
            if p1 == prevPos:
                trueCounter += 1
            else:
                prevPos = p1
                trueCounter = 0
            if trueCounter > maxDuplicates: continue

            # Fetching k-mer
            try:
                currStr = str(fastaFile.fetch(region.chrom, p1, p2)).upper()
            except Exception:
                continue
            if r.is_reverse: currStr = AuxiliaryFunctions.revcomp(currStr)

            # Counting k-mer in dictionary
            if not r.is_reverse:
                ct_reads_f += 1
                try:
                    obsDictF[currStr] += 1
                except Exception:
                    obsDictF[currStr] = 1
            else:
                ct_reads_r += 1
                try:
                    obsDictR[currStr] += 1
                except Exception:
                    obsDictR[currStr] = 1

        # Evaluating expected frequencies ####################################
        # Fetching whole sequence
        try:
            currStr = str(fastaFile.fetch(region.chrom, region.initial, region.final)).upper()
        except Exception:
            continue
        currRevComp = AuxiliaryFunctions.revcomp(currStr)

        # Iterating on each sequence position
        for i in range(0, len(currStr) - args.k_nb):
            ct_kmers += 1
            # Counting k-mer in dictionary
            s = currStr[i:i + args.k_nb]
            try:
                expDictF[s] += 1
            except Exception:
                expDictF[s] = 1

            # Counting k-mer in dictionary for reverse complement
            s = currRevComp[i:i + args.k_nb]
            try:
                expDictR[s] += 1
            except Exception:
                expDictR[s] = 1

    # Closing files
    bamFile.close()
    fastaFile.close()

    # Creating bias dictionary
    alphabet = ["A", "C", "G", "T"]
    kmerComb = ["".join(e) for e in product(alphabet, repeat=args.k_nb)]
    bias_table_F = dict([(e, 0.0) for e in kmerComb])
    bias_table_R = dict([(e, 0.0) for e in kmerComb])
    for kmer in kmerComb:
        try:
            obsF = obsDictF[kmer] + pseudocount
        except Exception:
            obsF = pseudocount
        try:
            expF = expDictF[kmer] + pseudocount
        except Exception:
            expF = pseudocount
        if ct_reads_f == 0:
            bias_table_F[kmer] = 1
        else:
            bias_table_F[kmer] = round(float(obsF / ct_reads_f) / float(expF / ct_kmers), 6)
        try:
            obsR = obsDictR[kmer] + pseudocount
        except Exception:
            obsR = pseudocount
        try:
            expR = expDictR[kmer] + pseudocount
        except Exception:
            expR = pseudocount
        if ct_reads_r == 0:
            bias_table_R[kmer] = 1
        else:
            bias_table_R[kmer] = round(float(obsR / ct_reads_r) / float(expR / ct_kmers), 6)

    write_table(args.output_location, args.output_prefix, [bias_table_F, bias_table_R])


def estimate_bias_pwm(args):
    # Parameters
    max_duplicates = 100

    # Initializing bam and fasta
    bamFile = Samfile(args.reads_file, "rb")
    genome_data = GenomeData(args.organism)
    fastaFile = Fastafile(genome_data.get_genome())
    regions = GenomicRegionSet("regions")
    regions.read(args.regions_file)

    obs_f_pwm_dict = dict([("A", [0.0] * args.k_nb), ("C", [0.0] * args.k_nb),
                           ("G", [0.0] * args.k_nb), ("T", [0.0] * args.k_nb), ("N", [0.0] * args.k_nb)])
    exp_f_pwm_dict = dict([("A", [0.0] * args.k_nb), ("C", [0.0] * args.k_nb),
                           ("G", [0.0] * args.k_nb), ("T", [0.0] * args.k_nb), ("N", [0.0] * args.k_nb)])
    obs_r_pwm_dict = dict([("A", [0.0] * args.k_nb), ("C", [0.0] * args.k_nb),
                           ("G", [0.0] * args.k_nb), ("T", [0.0] * args.k_nb), ("N", [0.0] * args.k_nb)])
    exp_r_pwm_dict = dict([("A", [0.0] * args.k_nb), ("C", [0.0] * args.k_nb),
                           ("G", [0.0] * args.k_nb), ("T", [0.0] * args.k_nb), ("N", [0.0] * args.k_nb)])

    # Iterating on HS regions
    for region in regions:
        # Initialization
        prev_pos = -1
        true_counter = 0

        # Evaluating observed frequencies
        # Fetching reads
        for r in bamFile.fetch(region.chrom, region.initial, region.final):
            # Calculating positions
            if not r.is_reverse:
                cut_site = r.pos + args.forward_shift - 1
                p1 = cut_site - int(floor(args.k_nb / 2))
            else:
                cut_site = r.aend + args.reverse_shift + 1
                p1 = cut_site - int(floor(args.k_nb / 2))
            p2 = p1 + args.k_nb

            # Verifying PCR artifacts
            if p1 == prev_pos:
                true_counter += 1
            else:
                prev_pos = p1
                true_counter = 0
            if true_counter > max_duplicates: continue

            # Fetching k-mer
            try:
                currStr = str(fastaFile.fetch(region.chrom, p1, p2)).upper()
            except Exception:
                continue
            if r.is_reverse: currStr = AuxiliaryFunctions.revcomp(currStr)

            # Counting k-mer in dictionary
            if not r.is_reverse:
                for i in range(0, len(currStr)):
                    obs_f_pwm_dict[currStr[i]][i] += 1
            else:
                for i in range(0, len(currStr)):
                    obs_r_pwm_dict[currStr[i]][i] += 1

        # Evaluating expected frequencies
        # Fetching whole sequence
        try:
            currStr = str(fastaFile.fetch(region.chrom, region.initial, region.final)).upper()
        except Exception:
            continue

        # Iterating on each sequence position
        s = None
        for i in range(0, len(currStr) - args.k_nb):
            # Counting k-mer in dictionary
            s = currStr[i:i + args.k_nb]
            for i in range(0, len(s)):
                exp_f_pwm_dict[s[i]][i] += 1

            # Counting k-mer in dictionary for reverse complement
            s = AuxiliaryFunctions.revcomp(s)
            for i in range(0, len(s)):
                exp_r_pwm_dict[s[i]][i] += 1

    # Closing files
    bamFile.close()
    fastaFile.close()

    # Output pwms
    os.system("mkdir -p " + os.path.join(args.output_location, "pfm"))
    pwm_dict_list = [obs_f_pwm_dict, obs_r_pwm_dict, exp_f_pwm_dict, exp_r_pwm_dict]
    pwm_file_list = []
    pwm_obs_f = os.path.join(args.output_location, "pfm", "obs_{}_f.pfm".format(str(args.k_nb)))
    pwm_obs_r = os.path.join(args.output_location, "pfm", "obs_{}_r.pfm".format(str(args.k_nb)))
    pwm_exp_f = os.path.join(args.output_location, "pfm", "exp_{}_f.pfm".format(str(args.k_nb)))
    pwm_exp_r = os.path.join(args.output_location, "pfm", "exp_{}_r.pfm".format(str(args.k_nb)))

    pwm_file_list.append(pwm_obs_f)
    pwm_file_list.append(pwm_obs_r)
    pwm_file_list.append(pwm_exp_f)
    pwm_file_list.append(pwm_exp_r)

    for i in range(len(pwm_dict_list)):
        with open(pwm_file_list[i], "w") as pwm_file:
            for e in ["A", "C", "G", "T"]:
                pwm_file.write(" ".join([str(int(f)) for f in pwm_dict_list[i][e]]) + "\n")

    motif_obs_f = motifs.read(open(pwm_obs_f), "pfm")
    motif_obs_r = motifs.read(open(pwm_obs_r), "pfm")
    motif_exp_f = motifs.read(open(pwm_exp_f), "pfm")
    motif_exp_r = motifs.read(open(pwm_exp_r), "pfm")

    # Output logos
    os.system("mkdir -p " + os.path.join(args.output_location, "logo"))
    logo_obs_f = os.path.join(args.output_location, "logo", "obs_{}_f.pdf".format(str(args.k_nb)))
    logo_obs_r = os.path.join(args.output_location, "logo", "obs_{}_r.pdf".format(str(args.k_nb)))
    logo_exp_f = os.path.join(args.output_location, "logo", "exp_{}_f.pdf".format(str(args.k_nb)))
    logo_exp_r = os.path.join(args.output_location, "logo", "exp_{}_r.pdf".format(str(args.k_nb)))

    motif_obs_f.weblogo(logo_obs_f, format="pdf", stack_width="large", color_scheme="color_classic",
                        yaxis_scale=0.2, yaxis_tic_interval=0.1)
    motif_obs_r.weblogo(logo_obs_r, format="pdf", stack_width="large", color_scheme="color_classic",
                        yaxis_scale=0.2, yaxis_tic_interval=0.1)
    motif_exp_f.weblogo(logo_exp_f, format="pdf", stack_width="large", color_scheme="color_classic",
                        yaxis_scale=0.02, yaxis_tic_interval=0.01)
    motif_exp_r.weblogo(logo_exp_r, format="pdf", stack_width="large", color_scheme="color_classic",
                        yaxis_scale=0.02, yaxis_tic_interval=0.01)

    # Creating bias dictionary
    alphabet = ["A", "C", "G", "T"]
    k_mer_comb = ["".join(e) for e in product(alphabet, repeat=args.k_nb)]
    bias_table_F = dict([(e, 0.0) for e in k_mer_comb])
    bias_table_R = dict([(e, 0.0) for e in k_mer_comb])
    for k_mer in k_mer_comb:
        obs_f = get_ppm_score(k_mer, motif_obs_f.pwm, args.k_nb)
        exp_f = get_ppm_score(k_mer, motif_exp_f.pwm, args.k_nb)
        bias_table_F[k_mer] = round(obs_f / exp_f, 6)
        obs_r = get_ppm_score(k_mer, motif_obs_r.pwm, args.k_nb)
        exp_r = get_ppm_score(k_mer, motif_exp_r.pwm, args.k_nb)
        bias_table_R[k_mer] = round(obs_r / exp_r, 6)

    write_table(args.output_location, args.output_prefix, [bias_table_F, bias_table_R])


def get_ppm_score(sequence, ppm, k_nb):
    score = 1.0
    for position in range(k_nb):
        letter = sequence[position]
        score *= ppm[letter][position]
    return score


def estimate_bias_vom(args):
    regions = GenomicRegionSet("regions")
    regions.read(args.regions_file)
    create_signal(args, regions)

    hmm_data = HmmData()

    learn_dependency_model = hmm_data.get_dependency_model()
    slim_dimont_predictor = hmm_data.get_slim_dimont_predictor()
    test_fa = hmm_data.get_default_test_fa()

    shutil.copy(test_fa, args.output_location)
    os.chdir(args.output_location)
    print((os.getcwd()))

    output_fname_f_obs = os.path.join(args.output_location, "{}_f_obs.fa".format(str(args.k_nb)))
    output_fname_f_exp = os.path.join(args.output_location, "{}_f_exp.fa".format(str(args.k_nb)))
    output_fname_r_obs = os.path.join(args.output_location, "{}_r_obs.fa".format(str(args.k_nb)))
    output_fname_r_exp = os.path.join(args.output_location, "{}_r_exp.fa".format(str(args.k_nb)))

    infix = "{}_f_obs".format(str(args.k_nb))
    create_model(args, output_fname_f_obs, infix, learn_dependency_model, slim_dimont_predictor)

    infix = "{}_f_exp".format(str(args.k_nb))
    create_model(args, output_fname_f_exp, infix, learn_dependency_model, slim_dimont_predictor)

    infix = "{}_r_obs".format(str(args.k_nb))
    create_model(args, output_fname_r_obs, infix, learn_dependency_model, slim_dimont_predictor)

    infix = "{}_r_exp".format(str(args.k_nb))
    create_model(args, output_fname_r_exp, infix, learn_dependency_model, slim_dimont_predictor)

    os.remove(os.path.join(args.output_location, "test.fa"))

    compute_bias(args)


def seq2vector(sequence):
    res = list()
    for char in sequence:
        if char == 'A':
            res.append(0)
        elif char == 'C':
            res.append(1)
        elif char == 'G':
            res.append(2)
        elif char == 'T':
            res.append(3)

    return res


def create_model(args, seq_file, infix, learn_dependency_model, slim_dimont_predictor):
    subprocess.call(["java", "-jar", learn_dependency_model, "i=Tabular", "m=LSlim", "is={}".format(seq_file),
                     "outdir={}".format(args.output_location)])

    output_file = os.path.join(args.output_location, "{}.txt".format(infix))
    slim_dimont_classifier = os.path.join(args.output_location, "SlimDimont_classifier.xml")

    with open(output_file, "w") as f:
        subprocess.call(["java", "-jar", slim_dimont_predictor, "slimdimont={}".format(slim_dimont_classifier),
                         "data=test.fa", "infix={}".format(infix)], stdout=f)

    # os.remove(os.path.join(args.output_location, "Dependency_logo.pdf"))
    os.remove(os.path.join(args.output_location, "Predicted_sequence_orientations_and_scores.tsv"))
    os.remove(os.path.join(args.output_location, "protocol_learn.txt"))
    os.remove(os.path.join(args.output_location, "{}-predictions.txt".format(infix)))
    os.remove(slim_dimont_classifier)
    os.remove(seq_file)


def create_signal(args, regions):
    def revcomp(s):
        rev_dict = dict([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"), ("N", "N")])
        return "".join([rev_dict[e] for e in s[::-1]])

    alphabet = ["A", "C", "G", "T"]
    kmer_comb = ["".join(e) for e in product(alphabet, repeat=args.k_nb)]
    f_obs_dict = dict([(e, 0.0) for e in kmer_comb])
    r_obs_dict = dict([(e, 0.0) for e in kmer_comb])
    f_exp_dict = dict([(e, 0.0) for e in kmer_comb])
    r_exp_dict = dict([(e, 0.0) for e in kmer_comb])

    bam_file = Samfile(args.reads_file, "rb")
    genome_data = GenomeData(args.organism)
    fasta_file = Fastafile(genome_data.get_genome())

    for region in regions:
        # Fetching observed reads
        reads = bam_file.fetch(reference=region.chrom, start=region.initial, end=region.final)
        for read in reads:
            if not read.is_reverse:
                p1 = read.pos - int(floor(args.k_nb / 2)) + args.forward_shift - 1
            else:
                p1 = read.aend - int(floor(args.k_nb / 2)) + args.reverse_shift + 1
            p2 = p1 + args.k_nb
            try:
                dna_sequence_obs = str(fasta_file.fetch(region.chrom, p1, p2)).upper()
            except Exception:
                continue
            if 'N' not in dna_sequence_obs:
                if read.is_reverse:
                    dna_sequence_obs = revcomp(dna_sequence_obs)
                    r_obs_dict[dna_sequence_obs] += 1
                else:
                    f_obs_dict[dna_sequence_obs] += 1

        # Fetching whole sequence
        try:
            dna_sequence_exp = str(fasta_file.fetch(region.chrom, region.initial, region.final)).upper()
        except Exception:
            continue
        dna_sequence_exp_rev = revcomp(dna_sequence_exp)
        for i in range(0, len(dna_sequence_exp) - args.k_nb):
            s = dna_sequence_exp[i:i + args.k_nb]
            if "N" not in s:
                f_exp_dict[s] += 1
            s = dna_sequence_exp_rev[i:i + args.k_nb]
            if "N" not in s:
                r_exp_dict[s] += 1

    output_fname_f_obs = os.path.join(args.output_location, "{}_f_obs.fa".format(str(args.k_nb)))
    output_fname_f_exp = os.path.join(args.output_location, "{}_f_exp.fa".format(str(args.k_nb)))
    output_fname_r_obs = os.path.join(args.output_location, "{}_r_obs.fa".format(str(args.k_nb)))
    output_fname_r_exp = os.path.join(args.output_location, "{}_r_exp.fa".format(str(args.k_nb)))

    output_file_f_obs = open(output_fname_f_obs, "w")
    output_file_f_exp = open(output_fname_f_exp, "w")
    output_file_r_obs = open(output_fname_r_obs, "w")
    output_file_r_exp = open(output_fname_r_exp, "w")

    for kmer in list(r_obs_dict.keys()):
        if f_obs_dict[kmer] > 0:
            output_file_f_obs.write(kmer + "\t" + str(f_obs_dict[kmer]) + "\n")
    for kmer in list(r_obs_dict.keys()):
        if f_exp_dict[kmer] > 0:
            output_file_f_exp.write(kmer + "\t" + str(f_exp_dict[kmer]) + "\n")
    for kmer in list(r_obs_dict.keys()):
        if r_obs_dict[kmer] > 0:
            output_file_r_obs.write(kmer + "\t" + str(r_obs_dict[kmer]) + "\n")
    for kmer in list(r_obs_dict.keys()):
        if r_exp_dict[kmer] > 0:
            output_file_r_exp.write(kmer + "\t" + str(r_exp_dict[kmer]) + "\n")

    output_file_f_obs.close()
    output_file_f_exp.close()
    output_file_r_obs.close()
    output_file_r_exp.close()


def write_table(output_location, output_prefix, table):
    output_fname = os.path.join(output_location, "{}_F.txt".format(output_prefix))
    f = open(output_fname, "w")
    for t in list(table[0].keys()):
        f.write(t + "\t" + str(table[0][t]) + "\n")
    f.close()

    output_fname = os.path.join(output_location, "{}_R.txt".format(output_prefix))
    f = open(output_fname, "w")
    for t in list(table[1].keys()):
        f.write(t + "\t" + str(table[1][t]) + "\n")
    f.close()


def compute_proba(vector, sequence, model, k_nb):
    proba0 = model["P_0(x|c=0)"][vector[0]]

    key1 = "P_1(x|y={},c=1)".format(sequence[0])
    proba1 = model["P_1(c=0)"] * model["P_1(x|c=0)"][vector[1]] + \
             model["P_1(c=1)"] * model[key1][vector[1]]

    if k_nb == 2:
        return proba0 * proba1

    key1 = "P_2(x|y={},c=1)".format(sequence[0])
    key2 = "P_2(x|y={},c=1)".format(sequence[1])
    proba2 = model["P_2(c=0)"] * model["P_2(x|c=0)"][vector[2]] + \
             model["P_2(c=1)"] * (model["P_2(p|c=1)"][0] * model[key1][vector[2]] +
                                  model["P_2(p|c=1)"][1] * model[key2][vector[2]])

    key1 = "P_3(x|y={},c=1)".format(sequence[0])
    key2 = "P_3(x|y={},c=1)".format(sequence[1])
    key3 = "P_3(x|y={},c=1)".format(sequence[2])
    proba3 = model["P_3(c=0)"] * model["P_3(x|c=0)"][vector[3]] + \
             model["P_3(c=1)"] * (model["P_3(p|c=1)"][0] * model[key1][vector[3]] +
                                  model["P_3(p|c=1)"][1] * model[key2][vector[3]] +
                                  model["P_3(p|c=1)"][2] * model[key3][vector[3]])
    if k_nb == 4:
        return proba0 * proba1 * proba2 * proba3

    key1 = "P_4(x|y={},c=1)".format(sequence[0])
    key2 = "P_4(x|y={},c=1)".format(sequence[1])
    key3 = "P_4(x|y={},c=1)".format(sequence[2])
    key4 = "P_4(x|y={},c=1)".format(sequence[3])
    proba4 = model["P_4(c=0)"] * model["P_4(x|c=0)"][vector[4]] + \
             model["P_4(c=1)"] * (model["P_4(p|c=1)"][0] * model[key1][vector[4]] +
                                  model["P_4(p|c=1)"][1] * model[key2][vector[4]] +
                                  model["P_4(p|c=1)"][2] * model[key3][vector[4]] +
                                  model["P_4(p|c=1)"][3] * model[key4][vector[4]])

    key1 = "P_5(x|y={},c=1)".format(sequence[0])
    key2 = "P_5(x|y={},c=1)".format(sequence[1])
    key3 = "P_5(x|y={},c=1)".format(sequence[2])
    key4 = "P_5(x|y={},c=1)".format(sequence[3])
    key5 = "P_5(x|y={},c=1)".format(sequence[4])
    proba5 = model["P_5(c=0)"] * model["P_5(x|c=0)"][vector[5]] + \
             model["P_5(c=1)"] * (model["P_5(p|c=1)"][0] * model[key1][vector[5]] +
                                  model["P_5(p|c=1)"][1] * model[key2][vector[5]] +
                                  model["P_5(p|c=1)"][2] * model[key3][vector[5]] +
                                  model["P_5(p|c=1)"][3] * model[key4][vector[5]] +
                                  model["P_5(p|c=1)"][4] * model[key5][vector[5]])
    if k_nb == 6:
        return proba0 * proba1 * proba2 * proba3 * proba4 * proba5
    else:
        key1 = "P_6(x|y={},c=1)".format(sequence[1])
        key2 = "P_6(x|y={},c=1)".format(sequence[2])
        key3 = "P_6(x|y={},c=1)".format(sequence[3])
        key4 = "P_6(x|y={},c=1)".format(sequence[4])
        key5 = "P_6(x|y={},c=1)".format(sequence[5])
        proba6 = model["P_6(c=0)"] * model["P_6(x|c=0)"][vector[6]] + \
                 model["P_6(c=1)"] * (model["P_6(p|c=1)"][0] * model[key1][vector[6]] +
                                      model["P_6(p|c=1)"][1] * model[key2][vector[6]] +
                                      model["P_6(p|c=1)"][2] * model[key3][vector[6]] +
                                      model["P_6(p|c=1)"][3] * model[key4][vector[6]] +
                                      model["P_6(p|c=1)"][4] * model[key5][vector[6]])

        if k_nb == 7:
            return proba0 * proba1 * proba2 * proba3 * proba4 * proba5 * proba6
        else:
            key1 = "P_7(x|y={},c=1)".format(sequence[2])
            key2 = "P_7(x|y={},c=1)".format(sequence[3])
            key3 = "P_7(x|y={},c=1)".format(sequence[4])
            key4 = "P_7(x|y={},c=1)".format(sequence[5])
            key5 = "P_7(x|y={},c=1)".format(sequence[6])
            proba7 = model["P_7(c=0)"] * model["P_7(x|c=0)"][vector[7]] + \
                     model["P_7(c=1)"] * (model["P_7(p|c=1)"][0] * model[key1][vector[7]] +
                                          model["P_7(p|c=1)"][1] * model[key2][vector[7]] +
                                          model["P_7(p|c=1)"][2] * model[key3][vector[7]] +
                                          model["P_7(p|c=1)"][3] * model[key4][vector[7]] +
                                          model["P_7(p|c=1)"][4] * model[key5][vector[7]])
            if k_nb == 8:
                return proba0 * proba1 * proba2 * proba3 * proba4 * proba5 * proba6 * proba7
            else:
                key1 = "P_8(x|y={},c=1)".format(sequence[3])
                key2 = "P_8(x|y={},c=1)".format(sequence[4])
                key3 = "P_8(x|y={},c=1)".format(sequence[5])
                key4 = "P_8(x|y={},c=1)".format(sequence[6])
                key5 = "P_8(x|y={},c=1)".format(sequence[7])
                proba8 = model["P_8(c=0)"] * model["P_8(x|c=0)"][vector[8]] + \
                         model["P_8(c=1)"] * (model["P_8(p|c=1)"][0] * model[key1][vector[8]] +
                                              model["P_8(p|c=1)"][1] * model[key2][vector[8]] +
                                              model["P_8(p|c=1)"][2] * model[key3][vector[8]] +
                                              model["P_8(p|c=1)"][3] * model[key4][vector[8]] +
                                              model["P_8(p|c=1)"][4] * model[key5][vector[8]])
                if k_nb == 9:
                    return proba0 * proba1 * proba2 * proba3 * proba4 * proba5 * proba6 * proba7 * proba8
                else:

                    key1 = "P_9(x|y={},c=1)".format(sequence[4])
                    key2 = "P_9(x|y={},c=1)".format(sequence[5])
                    key3 = "P_9(x|y={},c=1)".format(sequence[6])
                    key4 = "P_9(x|y={},c=1)".format(sequence[7])
                    key5 = "P_9(x|y={},c=1)".format(sequence[8])
                    proba9 = model["P_9(c=0)"] * model["P_9(x|c=0)"][vector[9]] + \
                             model["P_9(c=1)"] * (model["P_9(p|c=1)"][0] * model[key1][vector[9]] +
                                                  model["P_9(p|c=1)"][1] * model[key2][vector[9]] +
                                                  model["P_9(p|c=1)"][2] * model[key3][vector[9]] +
                                                  model["P_9(p|c=1)"][3] * model[key4][vector[9]] +
                                                  model["P_9(p|c=1)"][4] * model[key5][vector[9]])
                    if k_nb == 10:
                        return proba0 * proba1 * proba2 * proba3 * proba4 * proba5 * proba6 * proba7 * proba8 * proba9
                    else:
                        key1 = "P_10(x|y={},c=1)".format(sequence[5])
                        key2 = "P_10(x|y={},c=1)".format(sequence[6])
                        key3 = "P_10(x|y={},c=1)".format(sequence[7])
                        key4 = "P_10(x|y={},c=1)".format(sequence[8])
                        key5 = "P_10(x|y={},c=1)".format(sequence[9])
                        proba10 = model["P_10(c=0)"] * model["P_10(x|c=0)"][vector[10]] + \
                                  model["P_10(c=1)"] * (model["P_10(p|c=1)"][0] * model[key1][vector[10]] +
                                                        model["P_10(p|c=1)"][1] * model[key2][vector[10]] +
                                                        model["P_10(p|c=1)"][2] * model[key3][vector[10]] +
                                                        model["P_10(p|c=1)"][3] * model[key4][vector[10]] +
                                                        model["P_10(p|c=1)"][4] * model[key5][vector[10]])
                        if k_nb == 11:
                            return proba0 * proba1 * proba2 * proba3 * proba4 * proba5 * proba6 * proba7 * proba8 * \
                                   proba9 * proba10
                        else:
                            key1 = "P_11(x|y={},c=1)".format(sequence[6])
                            key2 = "P_11(x|y={},c=1)".format(sequence[7])
                            key3 = "P_11(x|y={},c=1)".format(sequence[8])
                            key4 = "P_11(x|y={},c=1)".format(sequence[9])
                            key5 = "P_11(x|y={},c=1)".format(sequence[10])
                            proba11 = model["P_11(c=0)"] * model["P_11(x|c=0)"][vector[10]] + \
                                      model["P_11(c=1)"] * (model["P_11(p|c=1)"][0] * model[key1][vector[10]] +
                                                            model["P_11(p|c=1)"][1] * model[key2][vector[10]] +
                                                            model["P_11(p|c=1)"][2] * model[key3][vector[10]] +
                                                            model["P_11(p|c=1)"][3] * model[key4][vector[10]] +
                                                            model["P_11(p|c=1)"][4] * model[key5][vector[10]])
                            return proba0 * proba1 * proba2 * proba3 * proba4 * proba5 * proba6 * proba7 * \
                                   proba8 * proba9 * proba10 * proba11


def read_model(input_fname, k_nb):
    model_list = list()
    with open(input_fname, "r") as input_file:
        lines = input_file.readlines()
        for idx in range(len(lines)):
            if lines[idx].startswith("Motif model"):
                setting_dict = dict()

                ll = lines[idx + 1]  # read motif probability
                setting_dict["motif probability"] = float(ll.split(":")[-1])

                ll = lines[idx + 2]  # read P_0
                setting_dict["P_0(c=0)"] = float(ll.split("=")[-1])
                ll = lines[idx + 3]  # read P_0
                setting_dict["P_0(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                ll = lines[idx + 6]  # read P_1
                setting_dict["P_1(c=0)"] = float(ll.split("=")[-1])
                ll = lines[idx + 7]  # read P_1
                setting_dict["P_1(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 9]  # read P_1
                setting_dict["P_1(c=1)"] = float(ll.split("=")[-1])
                ll = lines[idx + 10]  # read P_1
                setting_dict["P_1(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 11]  # read P_1
                setting_dict["P_1(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 12]  # read P_1
                setting_dict["P_1(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 13]  # read P_1
                setting_dict["P_1(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 14]  # read P_1
                setting_dict["P_1(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                if k_nb == 2:
                    model_list.append(setting_dict)
                    continue
                ll = lines[idx + 17]  # read P_2
                setting_dict["P_2(c=0)"] = float(ll.split("=")[-1])
                ll = lines[idx + 18]  # read P_2
                setting_dict["P_2(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 20]  # read P_2
                setting_dict["P_2(c=1)"] = float(ll.split("=")[-1])
                ll = lines[idx + 21]  # read P_2
                setting_dict["P_2(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 22]  # read P_2
                setting_dict["P_2(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 23]  # read P_2
                setting_dict["P_2(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 24]  # read P_2
                setting_dict["P_2(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 25]  # read P_2
                setting_dict["P_2(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                ll = lines[idx + 28]  # read P_3
                setting_dict["P_3(c=0)"] = float(ll.split("=")[-1])
                ll = lines[idx + 29]  # read P_3
                setting_dict["P_3(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 31]  # read P_3
                setting_dict["P_3(c=1)"] = float(ll.split("=")[-1])
                ll = lines[idx + 32]  # read P_3
                setting_dict["P_3(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 33]  # read P_3
                setting_dict["P_3(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 34]  # read P_3
                setting_dict["P_3(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 35]  # read P_3
                setting_dict["P_3(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 36]  # read P_3
                setting_dict["P_3(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                if k_nb == 4:
                    model_list.append(setting_dict)
                    continue

                ll = lines[idx + 39]  # read P_4
                setting_dict["P_4(c=0)"] = float(ll.split("=")[-1])
                ll = lines[idx + 40]  # read P_4
                setting_dict["P_4(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 42]  # read P_4
                setting_dict["P_4(c=1)"] = float(ll.split("=")[-1])
                ll = lines[idx + 43]  # read P_4
                setting_dict["P_4(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 44]  # read P_4
                setting_dict["P_4(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 45]  # read P_4
                setting_dict["P_4(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 46]  # read P_4
                setting_dict["P_4(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                ll = lines[idx + 47]  # read P_4
                setting_dict["P_4(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                if k_nb == 5:
                    model_list.append(setting_dict)
                    continue
                else:
                    ll = lines[idx + 50]  # read P_5
                    setting_dict["P_5(c=0)"] = float(ll.split("=")[-1])
                    ll = lines[idx + 51]  # read P_5
                    setting_dict["P_5(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                    ll = lines[idx + 53]  # read P_5
                    setting_dict["P_5(c=1)"] = float(ll.split("=")[-1])
                    ll = lines[idx + 54]  # read P_5
                    setting_dict["P_5(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                    ll = lines[idx + 55]  # read P_5
                    setting_dict["P_5(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                    ll = lines[idx + 56]  # read P_5
                    setting_dict["P_5(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                    ll = lines[idx + 57]  # read P_5
                    setting_dict["P_5(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                    ll = lines[idx + 58]  # read P_5
                    setting_dict["P_5(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                    if k_nb == 6:
                        model_list.append(setting_dict)
                        continue
                    else:
                        ll = lines[idx + 61]  # read P_6
                        setting_dict["P_6(c=0)"] = float(ll.split("=")[-1])
                        ll = lines[idx + 62]  # read P_6
                        setting_dict["P_6(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                        ll = lines[idx + 64]  # read P_6
                        setting_dict["P_6(c=1)"] = float(ll.split("=")[-1])
                        ll = lines[idx + 65]  # read P_6
                        setting_dict["P_6(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                        ll = lines[idx + 66]  # read P_6
                        setting_dict["P_6(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                        ll = lines[idx + 67]  # read P_6
                        setting_dict["P_6(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                        ll = lines[idx + 68]  # read P_6
                        setting_dict["P_6(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                        ll = lines[idx + 69]  # read P_6
                        setting_dict["P_6(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                        if k_nb == 7:
                            model_list.append(setting_dict)
                            continue
                        else:
                            ll = lines[idx + 72]  # read P_7
                            setting_dict["P_7(c=0)"] = float(ll.split("=")[-1])
                            ll = lines[idx + 73]  # read P_7
                            setting_dict["P_7(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                            ll = lines[idx + 75]  # read P_7
                            setting_dict["P_7(c=1)"] = float(ll.split("=")[-1])
                            ll = lines[idx + 76]  # read P_7
                            setting_dict["P_7(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                            ll = lines[idx + 77]  # read P_7
                            setting_dict["P_7(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                            ll = lines[idx + 78]  # read P_7
                            setting_dict["P_7(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                            ll = lines[idx + 79]  # read P_7
                            setting_dict["P_7(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                            ll = lines[idx + 80]  # read P_7
                            setting_dict["P_7(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                            if k_nb == 8:
                                model_list.append(setting_dict)
                                continue
                            else:
                                ll = lines[idx + 83]  # read P_8
                                setting_dict["P_8(c=0)"] = float(ll.split("=")[-1])
                                ll = lines[idx + 84]  # read P_8
                                setting_dict["P_8(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                ll = lines[idx + 86]  # read P_8
                                setting_dict["P_8(c=1)"] = float(ll.split("=")[-1])
                                ll = lines[idx + 87]  # read P_8
                                setting_dict["P_8(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                ll = lines[idx + 88]  # read P_8
                                setting_dict["P_8(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                ll = lines[idx + 89]  # read P_8
                                setting_dict["P_8(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                ll = lines[idx + 90]  # read P_8
                                setting_dict["P_8(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                ll = lines[idx + 91]  # read P_8
                                setting_dict["P_8(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                                if k_nb == 9:
                                    model_list.append(setting_dict)
                                    continue
                                else:
                                    ll = lines[idx + 94]  # read P_9
                                    setting_dict["P_9(c=0)"] = float(ll.split("=")[-1])
                                    ll = lines[idx + 95]  # read P_9
                                    setting_dict["P_9(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                    ll = lines[idx + 97]  # read P_9
                                    setting_dict["P_9(c=1)"] = float(ll.split("=")[-1])
                                    ll = lines[idx + 98]  # read P_9
                                    setting_dict["P_9(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                    ll = lines[idx + 99]  # read P_9
                                    setting_dict["P_9(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                    ll = lines[idx + 100]  # read P_9
                                    setting_dict["P_9(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                    ll = lines[idx + 101]  # read P_9
                                    setting_dict["P_9(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                    ll = lines[idx + 102]  # read P_9
                                    setting_dict["P_9(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                                    if k_nb == 10:
                                        model_list.append(setting_dict)
                                        continue
                                    else:
                                        ll = lines[idx + 105]  # read P_10
                                        setting_dict["P_10(c=0)"] = float(ll.split("=")[-1])
                                        ll = lines[idx + 106]  # read P_10
                                        setting_dict["P_10(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                        ll = lines[idx + 108]  # read P_10
                                        setting_dict["P_10(c=1)"] = float(ll.split("=")[-1])
                                        ll = lines[idx + 109]  # read P_10
                                        setting_dict["P_10(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                        ll = lines[idx + 110]  # read P_10
                                        setting_dict["P_10(x|y=A,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                        ll = lines[idx + 111]  # read P_10
                                        setting_dict["P_10(x|y=C,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                        ll = lines[idx + 112]  # read P_10
                                        setting_dict["P_10(x|y=G,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                        ll = lines[idx + 113]  # read P_10
                                        setting_dict["P_10(x|y=T,c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))

                                        if k_nb == 11:
                                            model_list.append(setting_dict)
                                            continue
                                        else:
                                            ll = lines[idx + 116]  # read P_11
                                            setting_dict["P_11(c=0)"] = float(ll.split("=")[-1])
                                            ll = lines[idx + 117]  # read P_11
                                            setting_dict["P_11(x|c=0)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                            ll = lines[idx + 119]  # read P_11
                                            setting_dict["P_11(c=1)"] = float(ll.split("=")[-1])
                                            ll = lines[idx + 120]  # read P_11
                                            setting_dict["P_11(p|c=1)"] = list(map(float, ll.split("=")[-1][1:-2].split(",")))
                                            ll = lines[idx + 121]  # read P_11
                                            setting_dict["P_11(x|y=A,c=1)"] = list(map(float,
                                                                                  ll.split("=")[-1][1:-2].split(",")))
                                            ll = lines[idx + 122]  # read P_11
                                            setting_dict["P_11(x|y=C,c=1)"] = list(map(float,
                                                                                  ll.split("=")[-1][1:-2].split(",")))
                                            ll = lines[idx + 123]  # read P_11
                                            setting_dict["P_11(x|y=G,c=1)"] = list(map(float,
                                                                                  ll.split("=")[-1][1:-2].split(",")))
                                            ll = lines[idx + 124]  # read P_11
                                            setting_dict["P_11(x|y=T,c=1)"] = list(map(float,
                                                                                  ll.split("=")[-1][1:-2].split(",")))

                                            model_list.append(setting_dict)
    return model_list


def compute_bias(args):
    input_f_obs = os.path.join(args.output_location, "{}_f_obs.txt".format(str(args.k_nb)))
    input_f_exp = os.path.join(args.output_location, "{}_f_exp.txt".format(str(args.k_nb)))
    input_r_obs = os.path.join(args.output_location, "{}_r_obs.txt".format(str(args.k_nb)))
    input_r_exp = os.path.join(args.output_location, "{}_r_exp.txt".format(str(args.k_nb)))

    model_list_f_obs = read_model(input_f_obs, args.k_nb)
    model_list_f_exp = read_model(input_f_exp, args.k_nb)
    model_list_r_obs = read_model(input_r_obs, args.k_nb)
    model_list_r_exp = read_model(input_r_exp, args.k_nb)

    os.remove(input_f_obs)
    os.remove(input_f_exp)
    os.remove(input_r_obs)
    os.remove(input_r_exp)

    # Creating bias dictionary
    alphabet = ["A", "C", "G", "T"]
    kmer_comb = ["".join(e) for e in product(alphabet, repeat=args.k_nb)]
    bias_table_f = dict([(e, 0.0) for e in kmer_comb])
    bias_table_r = dict([(e, 0.0) for e in kmer_comb])

    for kmer in kmer_comb:
        vector = seq2vector(kmer)

        # Compute observe probability of forward k-mer
        proba_f_obs = compute_proba(vector, kmer, model_list_f_obs[0], args.k_nb)

        # Compute expected probability of forward k-mer
        proba_f_exp = compute_proba(vector, kmer, model_list_f_exp[0], args.k_nb)

        # Compute observe probability of reverse k-mer
        proba_r_obs = compute_proba(vector, kmer, model_list_r_obs[0], args.k_nb)

        # Compute expected probability of reverse k-mer
        proba_r_exp = compute_proba(vector, kmer, model_list_r_exp[0], args.k_nb)

        # Compute bias
        bias_table_f[kmer] = min(round(proba_f_obs / proba_f_exp, 6), 10)
        bias_table_r[kmer] = min(round(proba_r_obs / proba_r_exp, 6), 10)

    write_table(args.output_location, args.output_prefix, [bias_table_f, bias_table_r])
