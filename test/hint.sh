#!/usr/bin/env bash

set -e

RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

DIR="${RGTTEST}/HINT"
mkdir -p $DIR
cd ${DIR}

echo "**********************************************"
echo "Testing HINT"

echo "Running HINT using only DNase-seq data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_DNaseTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_DNaseTest.tar.gz $url && tar xvfz HINT_DNaseTest.tar.gz && rm HINT_DNaseTest.tar.gz
cd ${DIR}/HINT_DNaseTest
rgt-hint footprinting --dnase-seq DNase.bam DNasePeaks.bed
echo "Running HINT-BC using only DNase-seq data.."
rgt-hint footprinting --dnase-seq --bias-correction DNase.bam DNasePeaks.bed

echo "Running HINT using only ATAC-seq data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_ATACTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_ATACTest.tar.gz $url && tar xvfz HINT_ATACTest.tar.gz && rm HINT_ATACTest.tar.gz
cd ${DIR}/HINT_ATACTest
echo "Running HINT using ATAC-seq data.."
rgt-hint footprinting --atac-seq ATAC.bam ATACPeaks.bed
echo "Testing the paired-end model of HINT using ATAC-seq data.."
rgt-hint footprinting --atac-seq --paired-end --output-prefix=fp_paired ATAC.bam ATACPeaks.bed

echo "Running HINT using only histone modification data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_HistoneTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_HistoneTest.tar.gz $url && tar xvfz HINT_HistoneTest.tar.gz && rm HINT_HistoneTest.tar.gz
cd ${DIR}/HINT_HistoneTest
rgt-hint footprinting --histone histone.bam histonePeaks.bed

echo "********* HINT test completed ****************"
