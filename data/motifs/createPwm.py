from __future__ import print_function
import os
import re
import argparse

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input-file', type=str, required=True, help='name of the input file')
parser.add_argument('-f', '--input-format', choices=['jaspar-2014', 'jaspar-2016', 'hocomoco-pcm', 'meme'], type=str,
                    required=True, help='format of the input file')
parser.add_argument('-o', '--output-folder', type=str, required=True, help='name of output Folder')

args = parser.parse_args()

# read the input file
with open(args.input_file, "r") as f:
    content = f.readlines()

n_lines = len(content)

###################################################################################################
# JASPAR 2014
###################################################################################################

if args.input_format == "jaspar-2014":
    for i in range(n_lines/5):
        motif_name = content[i * 5 + 0].strip()
        count_a = content[i * 5 + 1].strip()
        count_c = content[i * 5 + 2].strip()
        count_g = content[i * 5 + 3].strip()
        count_t = content[i * 5 + 4].strip()
        count_a = re.sub('\s+', ' ', count_a)
        count_c = re.sub('\s+', ' ', count_c)
        count_g = re.sub('\s+', ' ', count_g)
        count_t = re.sub('\s+', ' ', count_t)

        outputFileName = os.path.join(args.output_folder, "{}.pwm".format(motif_name.replace(">", "")))
        with open(outputFileName, "w") as f:
            f.write(count_a + "\n")
            f.write(count_c + "\n")
            f.write(count_g + "\n")
            f.write(count_t + "\n")

###################################################################################################
# JASPAR 2016
###################################################################################################

elif args.input_format == "jaspar-2016":
    for i in range(n_lines/5):
        motif_name = content[i * 5 + 0].replace(">", "").replace("\t", ".").replace("/", "_").strip()
        count_a = content[i * 5 + 1].translate(None, '[A]').strip()
        count_c = content[i * 5 + 2].translate(None, '[C]').strip()
        count_g = content[i * 5 + 3].translate(None, '[G]').strip()
        count_t = content[i * 5 + 4].translate(None, '[T]').strip()
        count_a = re.sub('\s+', ' ', count_a)
        count_c = re.sub('\s+', ' ', count_c)
        count_g = re.sub('\s+', ' ', count_g)
        count_t = re.sub('\s+', ' ', count_t)

        outputFileName = os.path.join(args.output_folder, "{}.pwm".format(motif_name))
        with open(outputFileName, "w") as f:
            f.write(count_a + "\n")
            f.write(count_c + "\n")
            f.write(count_g + "\n")
            f.write(count_t + "\n")

###################################################################################################
# HOCOMOCO
###################################################################################################

elif args.input_format == "hocomoco-pcm":

    count_a, count_c, count_g, count_t = [], [], [], []
    motif_name = content[0].strip(">")

    for i in range(n_lines):
        if content[i].startswith(">"):
            motif_name = content[i][1:].strip()
            count_a, count_c, count_g, count_t = [], [], [], []
        else:
            line = content[i].split("\t")
            count_a.append(str(int(round(float(line[0].strip())))))
            count_c.append(str(int(round(float(line[1].strip())))))
            count_g.append(str(int(round(float(line[2].strip())))))
            count_t.append(str(int(round(float(line[3].strip())))))

            if i < (n_lines-1):

                if content[i+1].startswith(">"):

                    outputFileName = os.path.join(args.output_folder, "{}.pwm".format(motif_name))
                    with open(outputFileName, "w") as f:
                        f.write(' '.join(count_a) + "\n")
                        f.write(' '.join(count_c) + "\n")
                        f.write(' '.join(count_g) + "\n")
                        f.write(' '.join(count_t) + "\n")

###################################################################################################
# MEME
###################################################################################################

elif args.input_format == "meme":

    # input_file: combined.meme
    # note: 3 Motifs per meme file

    count_a, count_c, count_g, count_t, count_m, count_1 = [], [], [], [], [], []
    motif_base_name = (os.path.dirname(args.input_file)).split("/")[-1]
    motif_name = motif_base_name
    # state = 0: header, state=1: reached motif, state=2: skipped empty line, state=3: process letter probability matrix
    state = 0
    nsites = 1

    for i, line in enumerate(content):

        if not content[i] == line:
            print("dp")

        if state == 0:
            line_content = line.strip().split(" ")
            if line_content[0] == "MOTIF":
                # set motif_name and initialize new lists
                # motif_name = <dir_name>_<motif_nr>, where <dir_name> is name of dir containing combined.meme
                motif_name = "_".join([motif_base_name, line_content[1]])
                count_a, count_c, count_g, count_t, count_m, count_1 = [], [], [], [], [], []
                state = 1

        elif state == 1:
            # skip this line
            state = 2

        elif state == 2:
            line_content = line.strip().split(" ")
            if line_content[6] == 'nsites=':
                # save nsites
                nsites = int(line_content[7])
            else:
                print("combined.meme file " + args.input_file + " has strange format")
            state = 3

        elif state == 3:
            line_content = [x.strip() for x in (line.strip().split("\t"))]

            # process motif
            # print(line_content[0])
            count_a.append(str(int(round(float(line_content[0])*nsites))))
            count_c.append(str(int(round(float(line_content[1])*nsites))))
            count_g.append(str(int(round(float(line_content[2])*nsites))))
            count_t.append(str(int(round(float(line_content[3])*nsites))))
            count_m.append(str(int(round(float(line_content[4])*nsites))))
            count_1.append(str(int(round(float(line_content[5])*nsites))))

            if (len(content) > (i+2) and content[i+2].strip().split(" ")[0] == "MOTIF") or i+2 == n_lines:
                # reached end of motif
                state = 0
                # save pwm
                outputFileName = os.path.join(args.output_folder, "{}.pwm".format(motif_name))
                with open(outputFileName, "w") as f:
                    f.write(' '.join(count_a) + "\n")
                    f.write(' '.join(count_c) + "\n")
                    f.write(' '.join(count_g) + "\n")
                    f.write(' '.join(count_t) + "\n")
                    f.write(' '.join(count_m) + "\n")
                    f.write(' '.join(count_1) + "\n")