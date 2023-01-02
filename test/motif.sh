#!/usr/bin/env bash

set -e

RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

DIR="${RGTTEST}/motifanalysis"
mkdir -p $DIR
cd ${DIR}

echo "**********************************************"
echo "******** Testing Motif Analysis tool *********"
echo "**********************************************"

# full-site test
echo
echo "# Full-Site Test"
echo

url="https://costalab.ukaachen.de/open_data/RGT/MotifAnalysis/RGT_MotifAnalysis_FullSiteTest.tar.gz"

# Download the data
DEST="RGT_MotifAnalysis_FullSiteTest"
if [ -d "$DEST" ]
then
    echo "dir found."
else
    echo "$DEST not found."
    wget -qO- -O ${DEST}.tar.gz $url && tar xvfz ${DEST}.tar.gz && rm ${DEST}.tar.gz
fi

# Run test script
cd $DEST
echo "## MATCHING"
echo
echo "### with filter"
rgt-motifanalysis matching --organism hg19 --filter "database:jaspar_vertebrates" --input-files input/regions_K562.bed input/background.bed
mv match/background_mpbs.bed match/background_mpbs_filter.bed
mv match/regions_K562_mpbs.bed match/regions_K562_mpbs_filter.bed
echo
echo "### without filter"
rgt-motifanalysis matching --organism hg19 --input-files input/regions_K562.bed input/background.bed

a=$(sort match/background_mpbs.bed | md5sum | awk '{ print $1 }')
b=$(sort match/background_mpbs_filter.bed | md5sum | awk '{ print $1 }')
c=$(sort match/regions_K562_mpbs.bed | md5sum | awk '{ print $1 }')
d=$(sort match/regions_K562_mpbs_filter.bed | md5sum | awk '{ print $1 }')

echo
if [ ${a} != ${b} ]
then
    echo "#### filter error in background_mpbs:"
    echo ${a}
    echo ${b}
fi
echo
if [ ${c} != ${d} ]
then
    echo "#### filter error in regions_K562_mpbs:"
    echo ${c}
    echo ${d}
fi

echo

echo "## ENRICHMENT"
rgt-motifanalysis enrichment --organism hg19 input/background.bed input/regions_K562.bed
echo

# Promoter test
cd $DIR
echo "# Promoter Test"
echo

url="https://costalab.ukaachen.de/open_data/RGT/MotifAnalysis/RGT_MotifAnalysis_PromoterTest.tar.gz"

# Download the data
DEST="RGT_MotifAnalysis_PromoterTest"
if [ -d "$DEST" ]
then
    echo "dir found."
else
    echo "$DEST not found."
    wget -qO- -O ${DEST}.tar.gz $url && tar xvfz ${DEST}.tar.gz && rm ${DEST}.tar.gz
fi

# Run test script
cd $DEST
echo "## MATCHING"
rgt-motifanalysis matching --organism hg19 --make-background --target-genes input/genes.txt --promoters-only
echo "## ENRICHMENT"
rgt-motifanalysis enrichment --organism hg19 --logo-copy match/background_regions.bed match/target_regions.bed

echo

# Gene-association test
echo "# Gene Association Test:"
echo

url="https://costalab.ukaachen.de/open_data/RGT/MotifAnalysis/RGT_MotifAnalysis_GeneAssocTest.tar.gz"

# Download the data
DEST="RGT_MotifAnalysis_GeneAssocTest"
cd $DIR
if [ -d "$DEST" ]
then
    echo "dir found."
else
    echo "$DEST not found."
    wget -qO- -O ${DEST}.tar.gz $url && tar xvfz ${DEST}.tar.gz && rm ${DEST}.tar.gz
fi

# Run test script
cd $DEST
echo "## MATCHING"
rgt-motifanalysis matching --organism hg19 --input-matrix input_matrix.txt --rand-proportion 10
echo "## ENRICHMENT"
rgt-motifanalysis enrichment --organism hg19 --logo-embed --input-matrix input_matrix.txt match/random_regions.bed

echo

echo "********* Motif Analysis test completed ****************"
