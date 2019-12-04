#!/usr/bin/env bash

set -e

RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

DIR="${RGTTEST}/THOR"
mkdir -p $DIR
cd ${DIR}

echo "**********************************************"
echo "Testing THOR"

# Download the data
file="${DIR}/THOR_example_data/THOR.config"
if [ -f "$file" ]
then
echo "$file found."
else
echo "$file not found."
curl http://www.regulatory-genomics.org/wp-content/uploads/2015/07/THOR_example_data.tar.gz | tar xz
fi

# Run test script
cd THOR_example_data/
rm -rf report_* sample*
rgt-THOR -n sample --report THOR.config
rgt-THOR -n sample_nomergebin --report THOR.config --no-merge-bin
rgt-THOR -n sample_nocor_pv1_merge --report THOR.config --no-correction --pvalue 1.0 --merge

echo "********* THOR test completed ****************"
