#!/usr/bin/env bash

set -e

# RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

VizDIR="./viz"
mkdir -p ${VizDIR}
cd ${VizDIR}

echo "**********************************************"
echo "Testing Viz"

# Download the data
file="./rgt_viz_example/scripts.sh"
if [ -f "$file" ]
then
    echo "Example data are loaded."
else
    echo "Downloading example files for rgt-viz"
    wget --no-check-certificate -qO- -O viz_examples.zip https://costalab.ukaachen.de/open_data/RGT/rgt_viz_example.zip
    unzip -o viz_examples.zip
    rm viz_examples.zip
fi

# Download the BW files
file="./rgt_viz_example/data/cDC_H3K4me1.bw"
if [ -f "$file" ]
then
    echo "Example data exist."
else
    echo "Downloading BW files for rgt-viz"
    cd rgt_viz_example/
    sh download_examples_RGT-viz.sh
    cd ..
fi

# Run test script
cd rgt_viz_example/
# Basic lineplot
rgt-viz lineplot Matrix_CDP.txt -o results -t lineplot_CDP -test
# Add one more cell type
rgt-viz lineplot Matrix_CDP_cDC.txt -o results -t lineplot_CDP_cDC -col cell -row regions -srow -test
# Chagne the layout
rgt-viz lineplot Matrix_CDP_cDC.txt -o results -t lineplot_CDP_cDC_2 -c cell -row reads -col regions -srow -test
# Projection test
rgt-viz projection -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t projection -g cell -organism mm9
# Jaccard test
rgt-viz jaccard -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t jaccard -g cell -organism mm9 -rt 10
# Intersection test
rgt-viz intersect -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t intersection -g cell -organism mm9
# With statistical test by randomization
rgt-viz intersect -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t intersection -g cell -organism mm9 -stest 10

echo "********* viz test completed ****************"
