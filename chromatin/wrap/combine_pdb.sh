#!/bin/bash

new_file=$1
n=1000000
if [[ $2 -ne 0 ]]; then
    n=$2
fi
touch $new_file

counter=1
for f in $(ls *.pdb); do
    echo $f
    ((counter++))
    cat $f >> $new_file
    echo "END" >> $new_file
    if [[ "$counter" -gt "$n" ]]; then
        break;
    fi
done