shiftx2 = "shiftx2_report.txt"
shifts = "shifts_diff.txt"

with open(shiftx2, "r") as f:
    lines_shiftx2 = f.readlines()

with open(shifts, "r") as f:
    lines_shifts = f.readlines()

second = False
with open("report.txt", "w") as f:
    for line_shiftx2, line_shifts in zip(lines_shiftx2,lines_shifts):
        if not (('First tail' in line_shiftx2) | ('Second tail' in line_shiftx2) | ('Mean for two tails' in line_shiftx2) | second):
            # parse shiftx2
            # print(len(line_shiftx2))
            shiftx2_resid = int(line_shiftx2.split()[0])
            shiftx2_resname = line_shiftx2.split()[1]
            shiftx2_H = float(line_shiftx2.split()[2])
            shiftx2_N = float(line_shiftx2.split()[3])
            shiftx2_CA = float(line_shiftx2.split()[4])
            shiftx2_CB = float(line_shiftx2.split()[5])
            shiftx2_C = float(line_shiftx2.split()[6])
            shiftx2_HA = float(line_shiftx2.split()[7])
            # parse shifts
            shifts_resid = int(line_shifts.split()[0])
            shifts_resname = line_shifts.split()[1]
            shifts_H = float(line_shifts.split()[2])
            shifts_N = float(line_shifts.split()[3])
            shifts_CA = float(line_shifts.split()[4])
            shifts_CB = float(line_shifts.split()[5])
            shifts_C = float(line_shifts.split()[6])
            shifts_HA = float(line_shifts.split()[7])
            # sum shiftx2 and shifts
            shiftx2_H += shifts_H
            shiftx2_N += shifts_N
            shiftx2_CA += shifts_CA
            shiftx2_CB += shifts_CB
            shiftx2_C += shifts_C
            shiftx2_HA += shifts_HA
            # output
            string = "%3d%6s%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n" % \
                (shifts_resid, shifts_resname, shiftx2_H, shiftx2_N,
                    shiftx2_CA, shiftx2_CB, shiftx2_C, shiftx2_HA)
            f.write(string)
        else:
            f.write(line_shiftx2)
            second = not second

