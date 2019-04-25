from bionmr_utils.md import *
import numpy as np

# Create PDB records description based on standard records
altered_records = AlteredPdbRecords(StandardPdbRecords.instance())

# Expand ATOM.serial record to columns 7-12
altered_records.alter_record(PdbRecordName("ATOM"), PdbFieldName("serial"), [7, 12])

pdb_filename = '../../input.pdb'
frame = PdbFile(pdb_filename, altered_records).get_frames()[0]

resnames = [r.rName.str for r in frame.asResidues]
resids = [r.rId.serial for r in frame.asResidues]

bins = np.zeros((len(resids)))
print(bins)

fn = 'contacts.txt'
N = 0
with open(fn) as f:
    for line in f:
        N += 1
        line = line.rstrip('\n')
        data = [int(x) for x in line.split()]
        for contact in data[1:]:
            bins[contact-1] += 1
print(N)
with open('contacts-details.txt', 'w') as f:
    f.write('rId,rName,number of contacts,percent of time\n')
    for i, r in enumerate(frame.asResidues):
        s = "%d," % (resids[i]) + resnames[i] + (",%d" % (bins[i]))
        s += ",%.2f\n" % (bins[i]/N*100)
        f.write(s)

