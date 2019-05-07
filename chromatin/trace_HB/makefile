SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/trace_HB

.PHONY: all clean plot init

clean:
	rm *.pdf
	rm *.py
	rm *.dat

init:
	cp ${SCRIPTS_PATH}/input.json .
	cp ${SCRIPTS_PATH}/plot_HB.py .

plot:
	python3 plot_HB.py

all:
	xmol_traj trace_hydrogen_bonds --input-config=input.json > hb_trace.dat
	python3 plot_HB.py
