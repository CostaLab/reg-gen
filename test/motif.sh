#!/usr/bin/env bash

set -e

RGTTEST=${RGTTEST:-"$HOME/rgt_test"}

DIR="${RGTTEST}/motifanalysis"
mkdir -p $DIR
cd ${DIR}

echo "**********************************************"
echo "Testing Motif Analysis"
echo

# full-site test
echo "Full-Site Test:"

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_FullSiteTest.tar.gz"

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
echo "Running matching.."
rgt-motifanalysis matching --organism hg19 --input-files input/regions_K562.bed input/background.bed
mv match/background_mpbs.bed match/background_mpbs_no_filter.bed
mv match/regions_K562_mpbs.bed match/regions_K562_mpbs_no_filter.bed
rgt-motifanalysis matching --organism hg19 --filter "database:jaspar_vertebrates" --input-files input/regions_K562.bed input/background.bed
mv match/background_mpbs.bed match/background_mpbs_filter.bed
mv match/regions_K562_mpbs.bed match/regions_K562_mpbs_filter.bed

a=$(md5sum match/background_mpbs_no_filter.bed)
b=$(md5sum match/background_mpbs_filter.bed)
c=$(md5sum match/regions_K562_mpbs_no_filter.bed)
d=$(md5sum match/regions_K562_mpbs_filter.bed)

if [ ${a} != ${b} ]
then
    echo "filter error in background_mpbs"
fi
if [ ${c} != ${d} ]
then
    echo "filter error in regions_K562_mpbs"
fi
echo "Running enrichment.."
rgt-motifanalysis enrichment --organism hg19 input/background.bed input/regions_K562.bed

echo

# Promoter test
cd $DIR
echo "Promoter Test:"

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_PromoterTest.tar.gz"

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
echo "Running matching.."
rgt-motifanalysis matching --organism hg19 --target-genes input/genes.txt --input-files input/background.bed
echo "Running enrichment.."
rgt-motifanalysis enrichment --organism hg19 --logo-copy input/background.bed match/target_regions.bed

echo

# Gene-association test
echo "Gene Association Test:"

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_GeneAssocTest.tar.gz"

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
echo "Running matching.."
rgt-motifanalysis matching --organism hg19 --input-matrix input_matrix.txt --rand-proportion 10
echo "Running enrichment.."
rgt-motifanalysis enrichment --organism hg19 --logo-embed --input-matrix input_matrix.txt match/random_regions.bed

echo

echo "********* Motif Analysis test completed ****************"