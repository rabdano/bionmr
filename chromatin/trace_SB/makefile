SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/trace_SB

.PHONY: all clean plot init

clean:
	rm *.pdf
	rm *.py
	rm *.dat

init:
	cp ${SCRIPTS_PATH}/input.json .
	cp ${SCRIPTS_PATH}/plot_SB.py .

plot:
	python3 plot_SB.py

all:
	xmol_traj trace_salt_bridges --input-config=input.json > sb_trace.dat
	python3 plot_HB.py
