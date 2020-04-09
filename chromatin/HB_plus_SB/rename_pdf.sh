#!/bin/bash

for f in $(ls ./Figures/*.pdf); do
    echo $f
    nf=${f/5_T/3LZ0_T}
    nf=${nf/4_/1KX5_ext_H4_}
    nf=${nf/3_/1KX5_}
    mv $f $nf
done