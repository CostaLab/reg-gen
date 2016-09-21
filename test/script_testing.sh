################################################################
# Test for RGT tools
################################################################
DIR=${HOME}/rgt_test

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
rgt-THOR THOR.config
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
	curl http://costalab.org/files/tdf/TDF_examples.zip | tar xz
fi

# Run test script
cd TDF_examples/FENDRR_mm9/
rgt-TDF promotertest -r FENDRR.fasta -de fendrr_gene_list.txt -organism mm9 -rn FENDRR -o promoter_test/FENDRR/

cd ../TERC_hg19/
rgt-TDF regiontest -r terc.fasta -bed terc_peaks.bed -rn TERC -f Nregions_hg19.bed -organism hg19 -o genomic_region_test/TERC

echo "********* TDF test completed ****************"

################################################################
# Viz
# echo "**********************************************"
# echo "Testing Viz"
# mkdir -p ${DIR}/viz
# cd ${DIR}/viz/

# # Download the data
# file="${DIR}/viz/TDF_examples/FENDRR_mm9/FENDRR.fasta"
# if [ -f "$file" ]
# then
# 	echo "$file found."
# else
# 	echo "$file not found."
# 	curl http://costalab.org/files/tdf/TDF_examples.zip | tar xz
# fi

# # Run test script
# cd TDF_examples/FENDRR_mm9/
# rgt-TDF promotertest -r FENDRR.fasta -de fendrr_gene_list.txt -organism mm9 -rn FENDRR -o promoter_test/FENDRR/

# cd ../TERC_hg19/
# rgt-TDF regiontest -r terc.fasta -bed terc_peaks.bed -rn TERC -f Nregions_hg19.bed -organism hg19 -o genomic_region_test/TERC -n 100

echo "********* TDF test completed ****************"
