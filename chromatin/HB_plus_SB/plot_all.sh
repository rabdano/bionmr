#!/bin/bash

export PYTHONPATH="/home/seva/Desktop/tmp/HB_plus_SB/v11:/home/seva/bionmr/home/seva/scripts/git/bionmr/chromatin/N15_relaxation"
mkdir -p Figures

for dn in $(ls -d ?_*); do
    echo $dn

    nf=${dn/5_T/3LZ0_T}
    nf=${nf/4_/1KX5_ext_H4_}
    nf=${nf/3_/1KX5_}

#    cp plot_HB_plus_SB_nonuniform.py $dn
    cd $dn
    python3.8 plot_HB_plus_SB_nonuniform.py > stdout.$(date +%s).txt
    cp H4-1_nonuniform.pdf ../Figures/$nf-H4-1_nonuniform.pdf
    cp H4-2_nonuniform.pdf ../Figures/$nf-H4-2_nonuniform.pdf
    cd ..
done
