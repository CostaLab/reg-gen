#!/bin/bash

echo "Cleaning up.."
rm jaspar7.db
rm annotation.csv

echo "Downloading SQL dump.."

wget -c http://jaspar.genereg.net/download/database/JASPAR2018.sql.tar.gz --user-agent="Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.8.1.6) Gecko/20070725 Firefox/2.0.0.6"

echo "Un-tarring it (it isn't actually gzipped).."
tar xvf JASPAR2018.sql.tar.gz || { echo "failed extracting archive"; exit 1; }
rm JASPAR2018.sql.tar.gz

echo "Converting to SQLite3.."

chmod 755 mysql2sqlite
./mysql2sqlite JASPAR2018.sql | sqlite3 jaspar7.db || { echo "failed converting DB to SQLite3"; exit 1; }

echo "Exporting annotation.."

sqlite3 jaspar7.db '.read jaspar_anno_query.sql' > annotation.csv || { echo "failed exporting sql table"; exit 1; }
