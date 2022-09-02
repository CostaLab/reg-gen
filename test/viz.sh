#!/usr/bin/env bash

set -e

RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

VizDIR="${RGTTEST}/viz"
mkdir -p $VizDIR
cd ${VizDIR}

echo "**********************************************"
echo "******** Testing RGT-Viz *********************"
echo "**********************************************"

# Download the data
DEST="${VizDIR}/rgt_viz_example"
if [ -d "$DEST" ]
then
    echo "dir found."
else
    echo "$DEST not found."
    echo "Downloading example files for rgt-viz"
    wget --no-check-certificate -qO- -O rgt_viz_example.zip https://costalab.ukaachen.de/open_data/RGT/rgt_viz_example.zip
    unzip -o rgt_viz_example.zip
    rm rgt_viz_example.zip
fi

# Run test script
cd ${VizDIR}/rgt_viz_example/
# Basic lineplot
rgt-viz lineplot Matrix_CDP.txt -o results -t lineplot_CDP
# Add one more cell type
rgt-viz lineplot Matrix_CDP_cDC.txt -o results -t lineplot_CDP_cDC -col cell -row regions -srow
# Chagne the layout
rgt-viz lineplot Matrix_CDP_cDC.txt -o results -t lineplot_CDP_cDC_2 -c cell -row reads -col regions -srow
# Projection test
rgt-viz projection -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t projection -g cell -organism mm9
# Jaccard test
rgt-viz jaccard -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t jaccard -g cell -organism mm9 -rt 10
# Intersection test
rgt-viz intersect -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t intersection -g cell -organism mm9
# With statistical test by randomization
rgt-viz intersect -r Matrix_H3K4me3.txt -q Matrix_PU1.txt -o viz_results -t intersection -g cell -organism mm9 -stest 10

echo "********* viz test completed ****************"
