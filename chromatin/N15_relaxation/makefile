N_DAT_FILES=300
FIRST_DAT_FILE=1
STRIDE=1
SCALING=1.005
HN_MASK=H
DILUTE_STRIDE=100
KEEP_ACF_FRACTION=0.8
REP=3
N_EXP=6

.PHONY: all clean

clean:
	rm -rf cor_NH_n*
	rm -rf vectors.txt
	rm -rf vectors.pickle
	rm -rf *_results.txt
	rm -rf *.png
	rm -rf *.pdf
	rm -rf *.py

init:
	cp /home/seva/scripts/N15_relaxation/get_NH_corfun_pyxmolpp2.py .
	sed -i 's/___N_DAT_FILES___/${N_DAT_FILES}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___STRIDE___/${STRIDE}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___HN_MASK___/${HN_MASK}/g' get_NH_corfun_pyxmolpp2.py
	sed -i 's/___FIRST_DAT_FILE___/${FIRST_DAT_FILE}/g' get_NH_corfun_pyxmolpp2.py

	cp /home/seva/scripts/N15_relaxation/dilute_autocorr.py .
	sed -i 's/___N_DAT_FILES___/${N_DAT_FILES}/g' dilute_autocorr.py
	sed -i 's/___STRIDE___/${STRIDE}/g' dilute_autocorr.py
	sed -i 's/___DILUTE_STRIDE___/${DILUTE_STRIDE}/g' dilute_autocorr.py
	sed -i 's/___KEEP_ACF_FRACTION___/${KEEP_ACF_FRACTION}/g' dilute_autocorr.py
	sed -i 's/___SCALING___/${SCALING}/g' dilute_autocorr.py
	
	cp /home/seva/scripts/N15_relaxation/calculate_R1_R2.py .
	sed -i 's/___N_DAT_FILES___/${N_DAT_FILES}/g' calculate_R1_R2.py
	sed -i 's/___DILUTE_STRIDE___/${DILUTE_STRIDE}/g' calculate_R1_R2.py
	sed -i 's/___N_EXP___/${N_EXP}/g' calculate_R1_R2.py
	sed -i 's/___REP___/${REP}/g' calculate_R1_R2.py

cor_NH_n${N_DAT_FILES}_s${STRIDE}: init
	python3 get_NH_corfun_pyxmolpp2.py

cor_NH_n${N_DAT_FILES}_s${DILUTE_STRIDE}_diluted: cor_NH_n${N_DAT_FILES}_s${STRIDE}
	python dilute_autocorr.py

R2_results.txt: cor_NH_n${N_DAT_FILES}_s${DILUTE_STRIDE}_diluted
	python calculate_R1_R2.py

all: R2_results.txt

