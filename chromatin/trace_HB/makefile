FIRST_DAT_FILE=51
LAST_DAT_FILE=300
STRIDE=10

FIRST_FRAME=$(shell echo $$(( ($(FIRST_DAT_FILE)-1)*1000 )) )

A=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 2)
B=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 3)
C=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 4)
ALPHA=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 5)
BETA=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 6)
GAMMA=$(shell tail -n 1 ../../1_build/box.inpcrd | cut -d " " -f 7)

CUT=3.5
ANGLE=30

SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/trace_HB


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

plot: | hb_trace.dat
	python3 plot_HB.py

all: | plot result.csv


input.json:
	cp ${SCRIPTS_PATH}/input.json .
	cp ${SCRIPTS_PATH}/plot_HB.py .
	cp ${SCRIPTS_PATH}/correlation_lifetime.py .
	cp ${SCRIPTS_PATH}/../wrap/get_summary.VOLUME.sh ../../5_run/get_summary.VOLUME.sh

	sed -i 's/___FIRST_FRAME___/${FIRST_FRAME}/g' input.json
	sed -i 's/___LAST_DAT_FILE___/${LAST_DAT_FILE}/g' input.json
	sed -i 's/___STRIDE___/${STRIDE}/g' input.json
	sed -i 's/___A___/${A}/g' input.json
	sed -i 's/___B___/${B}/g' input.json
	sed -i 's/___C___/${C}/g' input.json
	sed -i 's/___ALPHA___/${ALPHA}/g' input.json
	sed -i 's/___BETA___/${BETA}/g' input.json
	sed -i 's/___GAMMA___/${GAMMA}/g' input.json
	sed -i 's/___CUT___/${CUT}/g' input.json
	sed -i 's/___ANGLE___/${ANGLE}/g' input.json


hb_trace.dat: | input.json
	cd ../../5_run; bash get_summary.VOLUME.sh
	xmol_traj trace_hydrogen_bonds --input-config=input.json > hb_trace.dat

result.csv: | hb_trace.dat
	rm -rf Figures
	mkdir -p Figures
	python3 correlation_lifetime.py
