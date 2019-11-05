#!/bin/bash
# param -j for JASPAR, -h for HOCOMOCO

echo $1

# check parameters
# JASPAR: source=1, HOCOMOCO: source=2, both : source=0
source=0
if [ $# -eq 1 ]
    then
        if [ $1 == "-j" ]
            then
                source=1
        elif [ $1 == "-h" ]
            then
                source=2
    fi
fi

if [ ${source} -eq 0 ] || [ ${source} -eq 1 ]
    then
        echo "### Jaspar (2020) Annotation ###"

        echo "Downloading SQL dump.."

        wget -c http://jaspar.genereg.net/download/database/JASPAR2020.sql.tar.gz --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"

        echo "Un-tarring it (it isn't actually gzipped).."
        tar xvf JASPAR2020.sql.tar.gz || { echo "failed extracting archive"; exit 1; }

        echo "Converting to SQLite3.."

        chmod 755 mysql2sqlite
        ./mysql2sqlite JASPAR2020.sql | sqlite3 jaspar.db || { echo "failed converting DB to SQLite3"; exit 1; }

        echo "Exporting annotation.."

        sqlite3 jaspar.db '.read jaspar_anno_query.sql' | sed -e 's/"//g' > jaspar_anno.csv || { echo "failed exporting sql table"; exit 1; }

        echo "Cleaning up.."
        rm JASPAR2020.sql jaspar.db JASPAR2020.sql.tar.gz

        echo "Done"

        echo

        echo "Re-creating Mtf file.."
        python createMtf.py --jv --jp
fi

if [ ${source} -eq 0 ]
    then
        echo
fi

if [ ${source} -eq 0 ] || [ ${source} -eq 2 ]
    then
        echo "### Hocomoco Full v11 Annotation ###"
        echo
        echo "Downloading.."
        echo
        echo "Human"
        wget -c http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/HUMAN/mono/HOCOMOCOv11_full_annotation_HUMAN_mono.tsv

        echo "Mouse"
        wget -c http://hocomoco11.autosome.ru/final_bundle/hocomoco11/full/MOUSE/mono/HOCOMOCOv11_full_annotation_MOUSE_mono.tsv

        echo "Converting to CSV.."

        echo "TfName,GeneName,Family,UniProt,Source,TaxGroup,Species" > hocomoco_anno.csv
        tail -n+2 HOCOMOCOv11_full_annotation_HUMAN_mono.tsv | awk -F"\t" '{print $1","$2","$14","$19","$8",vertebrates,Homo sapiens"}' | sed -e 's/{.*}//g' >> hocomoco_anno.csv
        tail -n+2 HOCOMOCOv11_full_annotation_MOUSE_mono.tsv | awk -F"\t" '{print $1","$2","$14","$19","$8",vertebrates,Mus musculus"}' | sed -e 's/{.*}//g' >> hocomoco_anno.csv

        echo
        echo "Cleaning up.."
        rm HOCOMOCOv11_full_annotation_HUMAN_mono.tsv
        rm HOCOMOCOv11_full_annotation_MOUSE_mono.tsv

        echo
        echo "Done"

        echo

        echo "Re-creating Mtf file.."
        python createMtf.py --hoc
fi
