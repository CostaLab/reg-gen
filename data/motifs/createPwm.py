import os
import re
import argparse

from rgt.Util import npath

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input-file', type=str, required=True, help='name of the input file')
parser.add_argument('-f', '--input-format', choices=['jaspar-2014', 'jaspar-2016', 'hocomoco-pcm'],
                    type=str, required=True, help='format of the input file')
parser.add_argument('-o', '--output-folder', type=str, required=True, help='name of output Folder')

args = parser.parse_args()

# read the input file
with open(npath(args.input_file), "r") as f:
    content = f.readlines()

n_lines = len(content)

output_folder = npath(args.output_folder)

# make output directory path, if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

###################################################################################################
# JASPAR 2014
###################################################################################################

if args.input_format == "jaspar-2014":
    for i in range(n_lines // 5):
        motif_name = content[i * 5 + 0].strip()
        count_a = content[i * 5 + 1].strip()
        count_c = content[i * 5 + 2].strip()
        count_g = content[i * 5 + 3].strip()
        count_t = content[i * 5 + 4].strip()
        count_a = re.sub('\s+', ' ', count_a)
        count_c = re.sub('\s+', ' ', count_c)
        count_g = re.sub('\s+', ' ', count_g)
        count_t = re.sub('\s+', ' ', count_t)

        outputFileName = os.path.join(output_folder, "{}.pwm".format(motif_name.replace(">", "")))
        with open(outputFileName, "w") as f:
            f.write(count_a + "\n")
            f.write(count_c + "\n")
            f.write(count_g + "\n")
            f.write(count_t + "\n")

###################################################################################################
# JASPAR 2016
###################################################################################################

elif args.input_format == "jaspar-2016":
    for i in range(n_lines // 5):
        motif_name = content[i * 5 + 0].replace(">", "").replace("\t", ".").replace("/", "_").strip()
        count_a = content[i * 5 + 1].translate(str.maketrans('','',"[A]")).strip()
        count_c = content[i * 5 + 2].translate(str.maketrans('','',"[C]")).strip()
        count_g = content[i * 5 + 3].translate(str.maketrans('','',"[G]")).strip()
        count_t = content[i * 5 + 4].translate(str.maketrans('','',"[T]")).strip()
        count_a = re.sub('\s+', ' ', count_a)
        count_c = re.sub('\s+', ' ', count_c)
        count_g = re.sub('\s+', ' ', count_g)
        count_t = re.sub('\s+', ' ', count_t)

        outputFileName = os.path.join(output_folder, "{}.pwm".format(motif_name))
        with open(outputFileName, "w") as f:
            f.write(count_a + "\n")
            f.write(count_c + "\n")
            f.write(count_g + "\n")
            f.write(count_t + "\n")

###################################################################################################
# HOCOMOCO
###################################################################################################

elif args.input_format == "hocomoco-pcm":

    count_a = []
    count_c = []
    count_g = []
    count_t = []
    motif_name = content[0].strip(">")

    for i in range(n_lines):
        if content[i].startswith(">"):
            motif_name = content[i][1:].strip()
            count_a = []
            count_c = []
            count_g = []
            count_t = []
        else:
            line = content[i].split("\t")
            count_a.append(str(int(round(float(line[0].strip())))))
            count_c.append(str(int(round(float(line[1].strip())))))
            count_g.append(str(int(round(float(line[2].strip())))))
            count_t.append(str(int(round(float(line[3].strip())))))

            if i < (n_lines-1):

                if content[i+1].startswith(">"):

                    count_a = ' '.join(count_a)
                    count_c = ' '.join(count_c)
                    count_g = ' '.join(count_g)
                    count_t = ' '.join(count_t)

                    outputFileName = os.path.join(output_folder, f"{motif_name}.pwm")
                    with open(outputFileName, "w") as f:
                        f.write(count_a + "\n")
                        f.write(count_c + "\n")
                        f.write(count_g + "\n")
                        f.write(count_t + "\n")

