FIRST_DAT_FILE=51
LAST_DAT_FILE=300
STRIDE=1
SCALING=1.005
HN_MASK=H
DILUTE_STRIDE=100
KEEP_ACF_FRACTION=0.8
REP=3
N_EXP=6
SCRIPTS_PATH="/home/seva/scripts/git/bionmr/chromatin/N15_relaxation"

.PHONY: all clean

clean:
	rm -rf cor_NH_*
	rm -rf *_results.txt
	rm -rf *.png
	rm -rf *.pdf
	rm -rf *.py

init:
	cp ${SCRIPTS_PATH}/get_NH_corfun_pyxmolpp2.py .
	sed -i 's/___FIRST_DAT_FILE___/${FIRST_DAT_FILE}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___STRIDE___/${STRIDE}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___HN_MASK___/${HN_MASK}/g' get_NH_corfun_pyxmolpp2.py

	cp ${SCRIPTS_PATH}/dilute_autocorr.py .
	sed -i 's/___FIRST_DAT_FILE___/${FIRST_DAT_FILE}/g' dilute_autocorr.py
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' dilute_autocorr.py
	sed -i 's/___STRIDE___/${STRIDE}/g' dilute_autocorr.py
	sed -i 's/___DILUTE_STRIDE___/${DILUTE_STRIDE}/g' dilute_autocorr.py
	sed -i 's/___KEEP_ACF_FRACTION___/${KEEP_ACF_FRACTION}/g' dilute_autocorr.py
	sed -i 's/___SCALING___/${SCALING}/g' dilute_autocorr.py
	
	cp ${SCRIPTS_PATH}/calculate_R1_R2.py .
	sed -i 's/___FIRST_DAT_FILE___/${FIRST_DAT_FILE}/g' calculate_R1_R2.py
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' calculate_R1_R2.py
	sed -i 's/___DILUTE_STRIDE___/${DILUTE_STRIDE}/g' calculate_R1_R2.py
	sed -i 's/___N_EXP___/${N_EXP}/g' calculate_R1_R2.py
	sed -i 's/___REP___/${REP}/g' calculate_R1_R2.py

cor_NH_${FIRST_DAT_FILE}-${LAST_DAT_FILE}: init
	python3 get_NH_corfun_pyxmolpp2.py

cor_NH_${FIRST_DAT_FILE}-${LAST_DAT_FILE}_diluted: cor_NH_${FIRST_DAT_FILE}-${LAST_DAT_FILE}
	python3 dilute_autocorr.py

R2_results.txt: cor_NH_${FIRST_DAT_FILE}-${LAST_DAT_FILE}_diluted
	python3 calculate_R1_R2.py

all: R2_results.txt

