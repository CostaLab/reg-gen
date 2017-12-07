#!/bin/bash
################################################################
# Test for RGT tools
################################################################
DIR=${HOME}/rgt_test

# make any error stop the script
set -e

mkdir -p $DIR

################################################################
# THOR

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
rgt-THOR -n sample --report THOR.config

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
rgt-TDF promotertest -r FENDRR.fasta -de fendrr_gene_list_fold_change.txt -score -organism mm9 -rn FENDRR -o promoter_test -t FENDRR_FC/

cd ${DIR}/TDF/TDF_examples/TERC_hg19/
rgt-TDF regiontest -r terc.fasta -bed terc_peaks.bed -rn TERC -f Nregions_hg19.bed -organism hg19 -l 15 -o genomic_region_test/ -n 100 -mp 5

echo "********* TDF test completed ****************"

################################################################
# Viz
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

################################################################
# Motif Analysis
echo "**********************************************"
echo "Testing Motif Analysis"
mkdir -p ${DIR}/motifanalysis

cd ${DIR}/motifanalysis

# full-site test
echo "Full-Site Test:"

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_FullSiteTest.tar.gz"

# Download the data
dir="${DIR}/motifanalysis/RGT_MotifAnalysis_FullSiteTest/"
if [ -d "$dir" ]
then
    echo "dir found."
else
    echo "$dir not found."
    wget -qO- -O RGT_MotifAnalysis_FullSiteTest.tar.gz $url && tar xvfz RGT_MotifAnalysis_FullSiteTest.tar.gz && rm RGT_MotifAnalysis_FullSiteTest.tar.gz
fi

# Run test script
cd RGT_MotifAnalysis_FullSiteTest
echo "Running matching.."
rgt-motifanalysis matching --organism hg19 --input-files input/regions_K562.bed input/background.bed
echo "Running enrichment.."
rgt-motifanalysis enrichment --organism hg19 input/background.bed input/regions_K562.bed
cd ..

# Promoter test
echo "Promoter Test:"

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_PromoterTest.tar.gz"

# Download the data
dir="${DIR}/motifanalysis/RGT_MotifAnalysis_PromoterTest"
if [ -d "$dir" ]
then
    echo "dir found."
else
    echo "$dir not found."
    wget -qO- -O RGT_MotifAnalysis_PromoterTest.tar.gz $url && tar xvfz RGT_MotifAnalysis_PromoterTest.tar.gz && rm RGT_MotifAnalysis_PromoterTest.tar.gz
fi

# Run test script
cd RGT_MotifAnalysis_PromoterTest
echo "Running matching.."
rgt-motifanalysis matching --organism hg19 --target-genes input/genes.txt --make-background
echo "Running enrichment.."
rgt-motifanalysis enrichment --organism hg19 match/background_regions.bed match/target_regions.bed
cd ..

# Gene-association test
echo "Gene Association Test:"

url="http://www.regulatory-genomics.org/wp-content/uploads/2017/03/RGT_MotifAnalysis_GeneAssocTest.tar.gz"

# Download the data
dir="${DIR}/motifanalysis/RGT_MotifAnalysis_GeneAssocTest"
if [ -d "$dir" ]
then
    echo "dir found."
else
    echo "$dir not found."
    wget -qO- -O RGT_MotifAnalysis_GeneAssocTest.tar.gz $url && tar xvfz RGT_MotifAnalysis_GeneAssocTest.tar.gz && rm RGT_MotifAnalysis_GeneAssocTest.tar.gz
fi

# Run test script
cd RGT_MotifAnalysis_GeneAssocTest
echo "Running matching.."
rgt-motifanalysis matching --organism hg19 --input-matrix input_matrix.txt --rand-proportion 10
echo "Running enrichment.."
rgt-motifanalysis enrichment --organism hg19 --input-matrix input_matrix.txt match/random_regions.bed
cd ..

echo "********* Motif Analysis test completed ****************"

#################################################################
# HINT
echo "**********************************************"
echo "Testing HINT"
mkdir -p ${DIR}/HINT

cd ${DIR}/HINT

echo "Running HINT using only DNase-seq data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_DNaseTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_DNaseTest.tar.gz $url && tar xvfz HINT_DNaseTest.tar.gz && rm HINT_DNaseTest.tar.gz
cd HINT_DNaseTest
rgt-hint footprinting --dnase-seq DNase.bam DNasePeaks.bed
echo "Running HINT-BC using only DNase-seq data.."
rgt-hint footprinting --dnase-seq --bias-correction DNase.bam DNasePeaks.bed
cd ../

echo "Running HINT using only ATAC-seq data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_ATACTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_ATACTest.tar.gz $url && tar xvfz HINT_ATACTest.tar.gz && rm HINT_ATACTest.tar.gz
cd HINT_ATACTest
echo "Running HINT using ATAC-seq data.."
rgt-hint footprinting --atac-seq ATAC.bam ATACPeaks.bed
echo "Testing the paired-end model of HINT using ATAC-seq data.."
rgt-hint footprinting --atac-seq --paired-end --output-prefix=fp_paired ATAC.bam ATACPeaks.bed

echo "Running HINT using only histone modification data.."
url="http://134.130.18.8/open_data/hint/tutorial/HINT_HistoneTest.tar.gz"
echo "Downloading test data."
wget -qO- -O HINT_HistoneTest.tar.gz $url && tar xvfz HINT_HistoneTest.tar.gz && rm HINT_HistoneTest.tar.gz
cd HINT_HistoneTest
rgt-hint footprinting --histone histone.bam histonePeaks.bed
cd ../

echo "********* HINT test completed ****************"

#################################################################
# rgt-tools.py
