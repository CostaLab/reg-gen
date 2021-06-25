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
curl -k -O https://costalab.ukaachen.de/open_data/TDF/TDF_examples.zip && unzip TDF_examples.zip && rm TDF_examples.zip
fi


# Run test script
cd ${DIR}/TDF_examples/FENDRR_mm9/
rgt-TDF promotertest -r FENDRR.fasta -de fendrr_gene_list_fold_change.txt -score -organism mm9 -rn FENDRR -o promoter_test -t FENDRR_FC/ -l 15
rgt-TDF integrate -path promoter_test

cd ${DIR}/TDF_examples/MEG3_hg38/
rgt-TDF regiontest -r MEG3_sequence.fa -bed MEG3_hg38_CHOP.bed -rn MEG3 -o genomic_region_test -n 100 -organism hg38 -l 14 -mp 5 -ccf 100
rgt-TDF integrate -path genomic_region_test


echo "********* TDF test completed ****************"
