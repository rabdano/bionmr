FIRST_DAT_FILE=51
LAST_DAT_FILE=1000
STRIDE=100

SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/hbplus


.PHONY: all clean plot init


clean:
	rm -rf *.pdf
	rm -rf *.py
	rm -rf *.dat
	rm -rf *.json
	rm -rf *.txt
	rm -rf *.csv

../../5_run/get_summary.VOLUME.sh:
	cd ../../5_run; chmod u+x get_summary.VOLUME.sh; ./get_summary.VOLUME.sh

../../1_build/ref.pdb:
	cp ${SCRIPTS_PATH}/../vmd/make_ref.tcl ../../1_build/make_ref.tcl
	cd ../../1_build; /home/seva/bin/vmd -e make_ref.tcl

input.json:
	cp ${SCRIPTS_PATH}/hb_input.json .
	cp ${SCRIPTS_PATH}/extract_hb.py .

	sed -i 's/___FIRST_FRAME___/${FIRST_FRAME}/g' input.json
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' input.json
	sed -i 's/___STRIDE___/${STRIDE}/g' input.json

hb_trace.dat: | ../../5_run/get_summary.VOLUME.sh ../../1_build/ref.pdb input.json
	python3 extract_hb.py

all: | hb_trace.dat
