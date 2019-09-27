FIRST_DAT_FILE=51
LAST_DAT_FILE=800
STRIDE=100

FIRST_FRAME=$(shell echo $$(( ($(FIRST_DAT_FILE)-1)*1000 )) )

A=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 2)
B=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 3)
C=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 4)
ALPHA=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 5)
BETA=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 6)
GAMMA=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 7)

SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/trace_SB


.PHONY: all clean plot init


clean:
	rm -rf *.pdf
	rm -rf *.py
	rm -rf *.dat
	rm -rf *.json
	rm -rf *.txt
	rm -rf *.csv
	rm -rf Figures

init: | input.json

plot: | sb_trace.dat
	python3 plot_SB.py

all: | plot result.csv

input.json:
	cp ${SCRIPTS_PATH}/input.json .
	cp ${SCRIPTS_PATH}/plot_SB.py .
	cp ${SCRIPTS_PATH}/correlation_lifetime.py .
	cp ${SCRIPTS_PATH}/../wrap/get_summary.VOLUME.sh ../../5_run/get_summary.VOLUME.sh
	cp ${SCRIPTS_PATH}/../vmd/make_ref.tcl ../../1_build/make_ref.tcl

	sed -i 's/___FIRST_FRAME___/${FIRST_FRAME}/g' input.json
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' input.json
	sed -i 's/___STRIDE___/${STRIDE}/g' input.json
	sed -i 's/___A___/${A}/g' input.json
	sed -i 's/___B___/${B}/g' input.json
	sed -i 's/___C___/${C}/g' input.json
	sed -i 's/___ALPHA___/${ALPHA}/g' input.json
	sed -i 's/___BETA___/${BETA}/g' input.json
	sed -i 's/___GAMMA___/${GAMMA}/g' input.json

sb_trace.dat: | input.json
	cd ../../5_run; bash get_summary.VOLUME.sh
	cd ../../1_build; /home/seva/bin/vmd -e make_ref.tcl
	/home/seva/bin/xmol/release_build/src/xmol_traj trace_salt_bridges --input-config=input.json > sb_trace.dat

result.csv: | sb_trace.dat
	rm -rf Figures
	mkdir -p Figures
	python3 correlation_lifetime.py

