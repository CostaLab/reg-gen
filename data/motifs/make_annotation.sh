#!/bin/bash

echo "### Jaspar Vertebrates 7 (2018) Annotation ###"

echo "Downloading SQL dump.."

wget -c http://jaspar.genereg.net/download/database/JASPAR2018.sql.tar.gz --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"

echo "Un-tarring it (it isn't actually gzipped).."
tar xvf JASPAR2018.sql.tar.gz || { echo "failed extracting archive"; exit 1; }

echo "Converting to SQLite3.."

chmod 755 mysql2sqlite
./mysql2sqlite JASPAR2018.sql | sqlite3 jaspar.db || { echo "failed converting DB to SQLite3"; exit 1; }

echo "Exporting annotation.."

sqlite3 jaspar.db '.read jaspar_anno_query.sql' | sed -e 's/"//g' > jaspar_anno.csv || { echo "failed exporting sql table"; exit 1; }

echo "Cleaning up.."
rm ._JASPAR2018.sql JASPAR2018.sql jaspar.db JASPAR2018.sql.tar.gz

echo "Done"

echo

echo "### Hocomoco Full v11 Annotation ###"

echo "Downloading.."

wget -c http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv

echo "Converting to CSV.."

echo "TfName,GeneName,Family,UniProt,Source" > hocomoco_anno.csv
tail -n+2 HOCOMOCOv11_full_annotation_HUMAN_mono.tsv | awk -F"\t" '{print $1","$2","$14","$19","$8","}' | sed -e 's/{.*}//g' >> hocomoco_anno.csv

echo "Cleaning up.."
rm HOCOMOCOv11_full_annotation_HUMAN_mono.tsv

echo "Done"

echo "Re-creating Mtf file.."
python createMtf.py
