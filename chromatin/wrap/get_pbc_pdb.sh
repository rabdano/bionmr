#!/bin/bash

SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/wrap

cd 5_run

echo "Extracting volume information..."
cp $SCRIPTS_PATH/get_summary.VOLUME.sh .
chmod u+x get_summary.VOLUME.sh
./get_summary.VOLUME.sh

echo "Generating XST-table (PBC data)..."
cp $SCRIPTS_PATH/get_XST.py .
python3 get_XST.py

cd ..

echo "Collecting PDB files..."
mkdir -p pdb
cp 5_run/run*.pdb pdb

cd pdb

echo "Checking PBC translations..."
cp $SCRIPTS_PATH/wrap.py .
python3 wrap.py

cd ..
dn=${PWD##*/}
cd pdb_unwrapped
n=`ls *.pdb | wc -l`

echo "Combining $n PDB files to one..."
rm -f ../$dn.pdb
$SCRIPTS_PATH/combine_pdb.sh ../$dn.pdb $n

cd ..
