#!/bin/bash
################################################################
# Test for RGT tools
################################################################
DIR=${HOME}/rgt_test

mkdir -p $DIR

################################################################
# THOR
#
echo "**********************************************"
echo "Testing THOR"
mkdir -p ${DIR}/THOR
cd ${DIR}/THOR/

# Download the data
file="${DIR}/THOR/THOR_example_data/THOR.config"
if [ -f "$file" ]
then
echo "$file found."
else
echo "$file not found."
curl http://www.regulatory-genomics.org/wp-content/uploads/2015/07/THOR_example_data.tar.gz | tar xz
fi

Run test script
cd THOR_example_data/
rm -rf report_* sample-*
rgt-THOR THOR.config -n sample --report

echo "********* THOR test completed ****************"

################################################################
# TDF
echo "**********************************************"
echo "Testing TDF"
mkdir -p ${DIR}/TDF
cd ${DIR}/TDF/

# Download the data
file="${DIR}/TDF/TDF_examples/FENDRR_mm9/FENDRR.fasta"
if [ -f "$file" ]
then
echo "$file found."
else
echo "$file not found."
wget -qO- -O TDF_examples.zip http://costalab.org/files/tdf/TDF_examples.zip && unzip TDF_examples.zip && rm TDF_examples.zip
fi

# Run test script
cd ${DIR}/TDF/TDF_examples/FENDRR_mm9/
rgt-TDF promotertest -r FENDRR.fasta -de fendrr_gene_list.txt -organism mm9 -rn FENDRR -o promoter_test/

cd ${DIR}/TDF/TDF_examples/TERC_hg19/
rgt-TDF regiontest -r terc.fasta -bed terc_peaks.bed -rn TERC -f Nregions_hg19.bed -organism hg19 -l 15 -o genomic_region_test/ -n 100 -mp 5

echo "********* TDF test completed ****************"

#################################################################
## Viz
echo "**********************************************"
echo "Testing Viz"
mkdir -p ${DIR}/viz

cd ${DIR}/viz/

# Download the data
file="${DIR}/viz/viz_examples/scripts.sh"
if [ -f "$file" ]
then
	echo "Example data are loaded."
else
	echo "Downloading example files for rgt-viz"
    wget --no-check-certificate -qO- -O viz_examples.zip http://www.regulatory-genomics.org/wp-content/uploads/2016/09/rgt_viz_example.zip
    unzip -o viz_examples.zip
    rm viz_examples.zip
    mv zip viz_examples
fi

# Download the BW files
file="${DIR}/viz/viz_examples/data/cDC_H3K4me1.bw"
if [ -f "$file" ]
then
    echo "Example data exist."
else
    echo "Downloading BW files for rgt-viz"
    cd viz_examples/
    sh download_examples_RGT-viz.sh
fi

# Run test script
cd ${DIR}/viz/viz_examples/
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

#################################################################
## Motif Analysis
echo "**********************************************"
echo "Testing Motif Analysis"
mkdir -p ${DIR}/motifanalysis

cd ${DIR}/motifanalysis

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_FullSiteTest.tar.gz"

# Download the data
file="${DIR}/motifanalysis/RGT_MotifAnalysis_Test/input_matrix.txt"
if [ -f "$file" ]
then
    echo "$file found."
else
    echo "$file not found."
    wget -qO- -O RGT_MotifAnalysis_FullSiteTest.tar.gz $url && tar xvfz RGT_MotifAnalysis_FullSiteTest.tar.gz && rm RGT_MotifAnalysis_FullSiteTest.tar.gz
fi

# Run test script
cd RGT_MotifAnalysis_FullSiteTest
echo "Running matching.."
rgt-motifanalysis --matching input/regions_K562.bed input/background.bed
echo "Running enrichment.."
rgt-motifanalysis --enrichment input/background.bed input/regions_K562.bed

echo "********* Motif Analysis test completed ****************"

##################################################################
## HINT
echo "**********************************************"
echo "Testing HINT"
mkdir -p ${DIR}/HINT

cd ${DIR}/HINT

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/04/RGT_HINT_Test.tar.gz"

# Download the data
file="${DIR}/HINT/HINT_Test/InputMatrix_HINT_DNase.txt"
if [ -f "$file" ]
then
    echo "Test data exist."
else
    echo "Downloading test data."
    wget -qO- -O RGT_HINT_Test.tar.gz $url && tar xvfz RGT_HINT_Test.tar.gz && rm RGT_HINT_Test.tar.gz
fi

# Run test script
cd RGT_HINT_Test
echo "Running HINT using DNase-seq and histone modification data.."
rgt-hint --output-location ./Output/ --output-fname DU_K562_HINT ./InputMatrix_HINT_DNase+Histone.txt
echo "Running HINT using only DNase-seq data.."
rgt-hint --output-location ./Output/ --output-fname DU_K562_HINT ./InputMatrix_HINT_DNase.txt
echo "Running HINT-BC using only DNase-seq data.."
rgt-hint --default-bias-correction --output-location ./Output/ --output-fname DU_K562_HINT ./InputMatrix_HINTBC_DNase.txt

echo "********* HINT test completed ****************"
