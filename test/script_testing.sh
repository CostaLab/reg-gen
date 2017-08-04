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

# Run test script
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

echo "Running HINT using DNase-seq and histone modification data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_DNaseHistoneTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_DNaseHistoneTest.tar.gz $url && tar xvfz HINT_DNaseHistoneTest.tar.gz && rm HINT_DNaseHistoneTest.tar.gz
cd HINT_DNaseHistoneTest
rgt-hint --dnase-histone-footprints --output-location=./ --output-prefix=test DNase.bam H3K4me1.bam H3K4me3.bam regions.bed
cd ../

echo "Running HINT using only histone modification data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_HistoneTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_HistoneTest.tar.gz $url && tar xvfz HINT_HistoneTest.tar.gz && rm HINT_HistoneTest.tar.gz
cd HINT_HistoneTest
rgt-hint --histone-footprints --output-location=./ --output-prefix=test histone.bam histonePeaks.bed
cd ../

echo "Running HINT using only DNase-seq data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_DNaseTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_DNaseTest.tar.gz $url && tar xvfz HINT_DNaseTest.tar.gz && rm HINT_DNaseTest.tar.gz
cd HINT_DNaseTest
rgt-hint --dnase-footprints --output-location=./ --output-prefix=test DNase.bam DNasePeaks.bed
echo "Running HINT-BC using only DNase-seq data.."
rgt-hint --dnase-footprints --bias-correction --output-location=./ --output-prefix=test_bc DNase.bam DNasePeaks.bed
cd ../

echo "Running HINT using only ATAC-seq data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_ATAC.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_ATAC.tar.gz $url && tar xvfz HINT_ATAC.tar.gz && rm HINT_ATAC.tar.gz
cd HINT_ATAC
mkdir output
rgt-hint --atac-footprints --organism=mm10 input/B_ATAC_chr1.bam input/B_ATACPeaks_chr1.bed --output-location=output --output-prefix=B_ATAC_chr1_footprints

echo "********* HINT test completed ****************"
