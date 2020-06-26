#!/bin/bash

SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/N15_relaxation
export PYTHONPATH="${SCRIPTS_PATH}:${PYTHONPATH}"

cp ${SCRIPTS_PATH}/get_NH_corfun_pyxmolpp2.py .
cp ${SCRIPTS_PATH}/calculate_R1_R2.py .
cp ${SCRIPTS_PATH}/plot_autocorrelation_fit.py .
cp ${SCRIPTS_PATH}/plot_H4_R1.py .
cp ${SCRIPTS_PATH}/plot_H4_R2.py .
cp ${SCRIPTS_PATH}/plot_H4_R1_full.py .
cp ${SCRIPTS_PATH}/plot_H4_R2_full.py .

python3 get_NH_corfun_pyxmolpp2.py
python3 calculate_R1_R2.py
python3 plot_H4_R1_full.py
python3 plot_H4_R2_full.py
python3 plot_H4_R1.py
python3 plot_H4_R2.py
python3 plot_autocorrelation_fit.py