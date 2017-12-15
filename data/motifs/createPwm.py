import os
import sys

inputFileName = sys.argv[1]
outputFolder = sys.argv[2]
inputFormat = sys.argv[3]


with open(inputFileName, "r") as f:
    content = f.readlines()

n_lines = len(content)

if inputFormat == "jaspar-2014":
    for i in range(n_lines/5):
        motif_name = content[i * 5 + 0].strip()
        count_a = content[i * 5 + 1].strip()
        count_c = content[i * 5 + 2].strip()
        count_g = content[i * 5 + 3].strip()
        count_t = content[i * 5 + 4].strip()

        outputFileName = os.path.join(outputFolder, "{}.pwm".format(motif_name.replace(">", "")))
        with open(outputFileName, "w") as f:
            f.write(count_a + "\n")
            f.write(count_c + "\n")
            f.write(count_g + "\n")
            f.write(count_t + "\n")

elif inputFormat == "jaspar-2016":
    for i in range(n_lines/5):
        motif_name = content[i * 5 + 0].replace(">","").replace("\t", ".").strip()
        count_a = content[i * 5 + 1].replace("A","").replace("[","").replace("]","").strip()
        count_c = content[i * 5 + 2].replace("C","").replace("[","").replace("]","").strip()
        count_g = content[i * 5 + 3].replace("G","").replace("[","").replace("]","").strip()
        count_t = content[i * 5 + 4].replace("T","").replace("[","").replace("]","").strip()

        outputFileName = os.path.join(outputFolder, "{}.pwm".format(motif_name.replace(">", "")))
        with open(outputFileName, "w") as f:
            f.write(count_a + "\n")
            f.write(count_c + "\n")
            f.write(count_g + "\n")
            f.write(count_t + "\n")

else:
    print("Unknown input format.")
    exit(1)

