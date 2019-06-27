FIRST_PDB_FILE=1
LAST_PDB_FILE=5

SCRIPTS_PATH=/home/seva/scripts/git/bionmr/chromatin/cs

LAST_PDB_FILE_FMT=$(shell printf "%05d" ${LAST_PDB_FILE})

.PHONY: all clean

all: | shiftx2 shifts
	python3 sum_shiftx2_shifts.py
	python3 plot_cs.py

clean:
	rm -f *.txt
	rm -f *.png
	rm -rf pdb_shiftx2
	rm -rf pdb_shifts
	rm -rf ForkPoolWorker-*

init:
	cp ${SCRIPTS_PATH}/prepare_pdb_shiftx2.py .
	cp ${SCRIPTS_PATH}/calculate_shiftx2.py .
	cp ${SCRIPTS_PATH}/analyze_shiftx2.py .
	cp ${SCRIPTS_PATH}/prepare_pdb_shifts.py .
	cp ${SCRIPTS_PATH}/calculate_shifts.py .
	cp ${SCRIPTS_PATH}/analyze_shifts.py .
	cp ${SCRIPTS_PATH}/analyze_shifts_diff.py .
	cp ${SCRIPTS_PATH}/sum_shiftx2_shifts.py .
	cp ${SCRIPTS_PATH}/plot_cs.py .
	sed -i 's|___SCRIPTS_PATH___|${SCRIPTS_PATH}|g' plot_cs.py

../../pdb_unwrapped/run${LAST_PDB_FILE_FMT}.pdb:
	cd ../..; pwd; cp ${SCRIPTS_PATH}/../wrap/get_pbc_pdb.sh .; ./get_pbc_pdb.sh

pdb_shiftx2: | init ../../pdb_unwrapped/run${LAST_PDB_FILE_FMT}.pdb
	mkdir -p pdb_shiftx2
	for i in `seq -f "%05g" ${FIRST_PDB_FILE} ${LAST_PDB_FILE}`; do \
		echo "cp ../../pdb_unwrapped/run$$i.pdb pdb_shiftx2"; \
		cp ../../pdb_unwrapped/run$$i.pdb pdb_shiftx2; \
	done
	python3 prepare_pdb_shiftx2.py

shiftx2: | pdb_shiftx2
	python3 calculate_shiftx2.py
	python3 analyze_shiftx2.py

pdb_shifts: | init ../../pdb_unwrapped/run${LAST_PDB_FILE_FMT}.pdb
	mkdir -p pdb_shifts
	for i in `seq -f "%05g" ${FIRST_PDB_FILE} ${LAST_PDB_FILE}`; do \
		echo "cp ../../pdb_unwrapped/run$$i.pdb pdb_shifts"; \
		cp ../../pdb_unwrapped/run$$i.pdb pdb_shifts; \
	done
	python3 prepare_pdb_shifts.py

shifts: | pdb_shifts
	python3 calculate_shifts.py
	python3 analyze_shifts.py
	python3 analyze_shifts_diff.py
