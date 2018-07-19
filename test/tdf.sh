#!/usr/bin/env bash

set -e

RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

DIR="${RGTTEST}/TDF"
mkdir -p $DIR
cd ${DIR}

echo "**********************************************"
echo "Testing TDF"

# Download the data
file="${DIR}/TDF_examples/FENDRR_mm9/FENDRR.fasta"
if [ -f "$file" ]
then
echo "$file found."
else
echo "$file not found."
wget -qO- -O TDF_examples.zip http://costalab.org/files/tdf/TDF_examples.zip && unzip TDF_examples.zip && rm TDF_examples.zip
fi

# Run test script
#cd ${DIR}/TDF_examples/FENDRR_mm9/
#rgt-TDF promotertest -r FENDRR.fasta -de fendrr_gene_list.txt -organism mm9 -rn FENDRR -o promoter_test/
#rgt-TDF promotertest -r FENDRR.fasta -de fendrr_gene_list_fold_change.txt -score -organism mm9 -rn FENDRR -o promoter_test -t FENDRR_FC/

cd ${DIR}/TDF_examples/TERC_hg19/
rgt-TDF regiontest -r terc.fasta -bed terc_peaks.bed -rn TERC -f Nregions_hg19.bed -organism hg19 -l 15 -o genomic_region_test/ -n 10 -mp 5
rgt-TDF integrate -path genomic_region_test
echo "********* TDF test completed ****************"
