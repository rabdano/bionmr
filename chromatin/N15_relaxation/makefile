FIRST_DAT_FILE=51
LAST_DAT_FILE=300
STRIDE=1
SCALING=1.005
HN_MASK=H
ALING_DNA=True
DILUTE_STRIDE=100
KEEP_ACF_FRACTION=0.8
FFT_ACF=True
REMOVE_FIRST_POINT=True
REP=3
N_EXP=6
SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/N15_relaxation

.PHONY: all clean plot init

get_NH_corfun_pyxmolpp2.py:
	cp ${SCRIPTS_PATH}/get_NH_corfun_pyxmolpp2.py .
	sed -i 's/___FIRST_DAT_FILE___/${FIRST_DAT_FILE}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___STRIDE___/${STRIDE}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___HN_MASK___/${HN_MASK}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___ALING_DNA___/${ALING_DNA}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___SCALING___/${SCALING}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___KEEP_ACF_FRACTION___/${KEEP_ACF_FRACTION}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___FFT_ACF___/${FFT_ACF}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___REMOVE_FIRST_POINT___/${REMOVE_FIRST_POINT}/g' get_NH_corfun_pyxmolpp2.py
	
	cp ${SCRIPTS_PATH}/calculate_R1_R2.py .
	sed -i 's/___FIRST_DAT_FILE___/${FIRST_DAT_FILE}/g' calculate_R1_R2.py
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' calculate_R1_R2.py
	sed -i 's/___DILUTE_STRIDE___/${DILUTE_STRIDE}/g' calculate_R1_R2.py
	sed -i 's/___N_EXP___/${N_EXP}/g' calculate_R1_R2.py
	sed -i 's/___REP___/${REP}/g' calculate_R1_R2.py

	cp ${SCRIPTS_PATH}/plot_autocorrelation_fit.py .
	sed -i 's/___N_EXP___/${N_EXP}/g' plot_autocorrelation_fit.py
	sed -i 's/___FIRST_DAT_FILE___/${FIRST_DAT_FILE}/g' plot_autocorrelation_fit.py
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' plot_autocorrelation_fit.py

	cp ${SCRIPTS_PATH}/plot_H4_R1.py .
	cp ${SCRIPTS_PATH}/plot_H4_R2.py .

cor_NH_${FIRST_DAT_FILE}-${LAST_DAT_FILE}_diluted: | get_NH_corfun_pyxmolpp2.py
	python3 get_NH_corfun_pyxmolpp2.py

R2_results.txt: | cor_NH_${FIRST_DAT_FILE}-${LAST_DAT_FILE}_diluted
	export PYTHONPATH="${SCRIPTS_PATH}:${PYTHONPATH}" && python3 calculate_R1_R2.py

all: | R2_results.txt

plot:
	python3 plot_H4_R1.py
	python3 plot_H4_R2.py
	python3 plot_autocorrelation_fit.py

init: | get_NH_corfun_pyxmolpp2.py

clean:
	rm -rf cor_NH_*
	rm -rf *_results.txt
	rm -rf *.png
	rm -rf *.pdf
	rm -rf *.py