import numpy as np


def average_data(data, window=100):
    # average data over windows
    # e.g. if in window hydrogen boand is present > 50% then value for window is 1
    n, m = data.shape
    nw = int(n/window)
    result = np.zeros((nw, m))
    for start in range(0, len(data), window):
        result[int(start/window), :] = np.mean(data[start:start+window, :], axis=0)
    return result


def hb_split_header(header):
    # THR~3::N--DC~1133::OP1 THR~3::OG1--DC~1133::OP1
    result = []
    header = header.rstrip('\n').strip()
    hbs = header.split(' ')
    for i, hb in enumerate(hbs):
        donor, acceptor = tuple(hb.split('--'))

        donor_r = donor.split('::')[0]
        donor_aName = donor.split('::')[1]
        donor_rName = donor_r.split('~')[0]
        donor_rId = int(donor_r.split('~')[1])

        acceptor_r = acceptor.split('::')[0]
        acceptor_aName = acceptor.split('::')[1]
        acceptor_rName = acceptor_r.split('~')[0]
        acceptor_rId = int(acceptor_r.split('~')[1])

        result.append([[[donor_rName, donor_rId, donor_aName],
                        [acceptor_rName, acceptor_rId, acceptor_aName],
                        "HB"], i])
    return result


def sb_split_header(header, ref_frame):
    # 0.50.OE1+OE2--0.53.NH1 4.593.OD1+OD2--4.618.NH2
    result = []
    header = header.rstrip('\n').strip()
    sbs = header.split(' ')
    for i, sb in enumerate(sbs):
        donor, acceptor = tuple(sb.split('--'))

        donor_rId = int(donor.split('.')[1])
        donor_rName = three_to_one(resname_by_resid(donor_rId, ref_frame))
        donor_aName = str(donor.split('.')[2])

        acceptor_rId = int(acceptor.split('.')[1])
        acceptor_rName = three_to_one(resname_by_resid(acceptor_rId, ref_frame))
        acceptor_aName = str(acceptor.split('.')[2])

        result.append([[[donor_rName, donor_rId, donor_aName],
                        [acceptor_rName, acceptor_rId, acceptor_aName],
                        "SB"], i])
    return result


# def sb_split_header(header, ref_frame):
#     result = []
#     header = header.rstrip('\n').strip()
#     sbs = header.split(' ')
#     for sb in sbs:
#         donor, acceptor = tuple(sb.split('--'))
#         donor_rName = donor.split('~')[0]
#         donor_rId = int(donor.split('~')[1])
#         acceptor_rName = acceptor.split('~')[0]
#         acceptor_rId = int(acceptor.split('~')[1])
#         donor_aName = ""
#         acceptor_aName = ""
#         result.append([[donor_rName, donor_rId, donor_aName],
#                        [acceptor_rName, acceptor_rId, acceptor_aName],
#                        "SB"])
#     return result


def resname_by_resid(rid, ref_frame):
    resnames = []
    resids = []
    for res in ref_frame.asResidues:
        resnames.append(res.rName.str)
        resids.append(res.rId.serial)
    return resnames[resids.index(rid)]


def correct_numbering(interactions, renumber_map=None):
    if not renumber_map:
        renumber_map = {'A': [rid for rid in range(1, 136)],
                        'B': [rid for rid in range(1, 103)],
                        'C': [rid for rid in range(1, 129)],
                        'D': [rid for rid in range(1, 123)],
                        'E': [rid for rid in range(1, 136)],
                        'F': [rid for rid in range(1, 103)],
                        'G': [rid for rid in range(1, 129)],
                        'H': [rid for rid in range(1, 123)],
                        'I': [rid for rid in range(1, 146)],
                        'J': [rid for rid in range(1, 146)]}

    absolute_map = {}
    renumber_array = []
    last_resid = 0

    for k in sorted(renumber_map):
        v = renumber_map[k]
        absolute_map[k] = [rid for rid in range(last_resid + 1, last_resid + len(v) + 1)]
        renumber_array += v
        last_resid = len(renumber_array)

    for b in interactions:
        i1 = b[0][0][1]
        i2 = b[0][1][1]

        for k, v in absolute_map.items():
            if i1 in v:
                b[0][0].insert(0, k)
            if i2 in v:
                b[0][1].insert(0, k)
            if b[0][0][2] == "Na+":
                b[0][0].insert(0, "X")
            if b[0][1][2] == "Na+":
                b[0][1].insert(0, "X")

        if b[0][0][2] <= len(renumber_array):
            b[0][0][2] = renumber_array[b[0][0][2] - 1]
        if b[0][1][2] <= len(renumber_array):
            b[0][1][2] = renumber_array[b[0][1][2] - 1]

    return interactions


def prune_not_h4_contacts(interactions, data, rid_of_interest):
    to_del = []

    for i, b in enumerate(interactions):
        if not ((b[0][0][1] in rid_of_interest) or (b[0][1][1] in rid_of_interest)):
            to_del.append(i)

    interactions = [[b, i] for i, b in enumerate([interaction[0] for interaction in interactions]) if i not in to_del]
    data = np.delete(data, to_del, axis=1)

    return interactions, data


def prune_intra_tail_contacts(interactions, data, n_tail_resid):
    to_del = []

    for i, b in enumerate(interactions):
        dn_chain = b[0][0][0]
        dn_resid = b[0][0][2]
        ac_chain = b[0][1][0]
        ac_resid = b[0][1][2]

        if ((dn_chain in ['B', 'F']) and (dn_resid in range(1, n_tail_resid + 1)) and
                (ac_chain in ['B', 'F']) and (ac_resid in range(1, n_tail_resid + 1))):
            to_del.append(i)

        if ((ac_chain in ['B', 'F']) and (ac_resid in range(1, n_tail_resid + 1)) and
                (dn_chain in ['B', 'F']) and (dn_resid in range(1, n_tail_resid + 1))):
            to_del.append(i)

    interactions = [[b, i] for i, b in enumerate([interaction[0] for interaction in interactions]) if i not in to_del]
    data = np.delete(data, to_del, axis=1)

    return interactions, data


def prune_empty_traces(interactions, data, threshold=0.2):
    to_del = []

    for i in range(len(interactions)):
        if np.max(data[:, i]) <= threshold:
            to_del.append(i)

    interactions = [[b, i] for i, b in enumerate([interaction[0] for interaction in interactions]) if i not in to_del]
    data = np.delete(data, to_del, axis=1)

    return interactions, data


def joint_sb_data_for_same_residue(sbs, sb_avg_data):
    dn_ac_pairs = set()

    for i, sb in enumerate([s[0] for s in sbs]):
        # (dn_cName, dn_rId, ac_cName, ac_rId)
        dn_ac_pairs.add((sb[0][0], sb[0][2], sb[1][0], sb[1][2]))

    sbs_joined = []
    sb_avg_data_joined = np.zeros((len(sb_avg_data[:, 0]), len(dn_ac_pairs)), dtype=float)

    # # average
    # for i, pair in enumerate(dn_ac_pairs):
    #     dn_cname, dn_rid, ac_cname, ac_rid = pair
    #     counter = 0
    #
    #     for j, sb in enumerate([s[0] for s in sbs]):
    #         if (sb[0][0] == dn_cname) & (sb[0][2] == dn_rid) & \
    #                 (sb[1][0] == ac_cname) & (sb[1][2] == ac_rid):
    #             sb_avg_data_joined[:, i] += sb_avg_data[:, j]
    #             if counter == 0:
    #                 sbs_joined.append([sb, i])
    #             counter += 1
    #
    #     sb_avg_data_joined[:, i] /= counter

    # union
    for i, pair in enumerate(dn_ac_pairs):
        dn_cname, dn_rid, ac_cname, ac_rid = pair
        counter = 0

        for j, sb in enumerate([s[0] for s in sbs]):
            if (sb[0][0] == dn_cname) & (sb[0][2] == dn_rid) & \
                    (sb[1][0] == ac_cname) & (sb[1][2] == ac_rid):
                new_trace = np.zeros_like(sb_avg_data_joined[:, i])
                for k, p1, p2 in zip(range(len(sb_avg_data_joined[:, i])), sb_avg_data_joined[:, i], sb_avg_data[:, j]):
                    new_trace[k] = np.max([p1, p2])
                sb_avg_data_joined[:, i] = new_trace
                if counter == 0:
                    sbs_joined.append([sb, i])
                counter += 1

    return sbs_joined, sb_avg_data_joined


def prune_sodium_contacts(interactions, data):
    to_del = []

    for i, b in enumerate(interactions):
        if (b[0][0][0] == "X") | (b[0][1][0] == "X"):
            to_del.append(i)

    interactions = [[b, i] for i, b in enumerate([interaction[0] for interaction in interactions]) if i not in to_del]
    data = np.delete(data, to_del, axis=1)

    return interactions, data


def build_cmap(hbs, sbs, hb_avg_data, sb_avg_data, cname, n_tail_resid):
    cmap = []
    labels = []
    for resid in range(1, n_tail_resid + 1):
        cur_map = []
        cur_labels = []
        cur_partners = []
        n_hb, n_sb = 0, 0

        # find salt bridges for residue
        for i, sb in enumerate(sbs):
            if (sb[0][0][0] == cname) and (sb[0][0][2] == resid):
                cur_partners.append(sb[0][1][2])
            if (sb[0][1][0] == cname) and (sb[0][1][2] == resid):
                cur_partners.append(sb[0][0][2])

        # append hydrogen bonds to cmap and labels if no salt bridge with partner
        for i, hb in enumerate(hbs):
            if (hb[0][0][0] == cname) & (hb[0][0][2] == resid):
                condition = (hb[0][1][2] in cur_partners) and \
                            (hb[0][0][3] in ["NZ", "NH1", "NH2", "ND1", "NE2", "NE", "N"]) and \
                            (hb[0][1][3] in ["OP1", "OP2", "O3'", "O5'", "OD1", "OD2", "OE1", "OE2"])
                if not condition:
                    cur_labels.append(hb)
                    cur_map.append(hb_avg_data[:, i])
                    n_hb += 1
            elif (hb[0][1][0] == cname) & (hb[0][1][2] == resid):
                condition = (hb[0][0][2] in cur_partners) and \
                            (hb[0][1][3] in ["NZ", "NH1", "NH2", "ND1", "NE2", "NE", "N"]) and \
                            (hb[0][0][3] in ["OP1", "OP2", "O3'", "O5'", "OD1", "OD2", "OE1", "OE2"])
                if not condition:
                    cur_labels.append(hb)
                    cur_map.append(hb_avg_data[:, i])
                    n_hb += 1

        # append salt bridges to cmap and labels
        for i, sb in enumerate(sbs):
            if (sb[0][0][0] == cname) and (sb[0][0][2] == resid):
                cur_labels.append(sb)
                cur_map.append(sb_avg_data[:, i])
                n_sb += 1
            if (sb[0][1][0] == cname) and (sb[0][1][2] == resid):
                cur_labels.append(sb)
                cur_map.append(sb_avg_data[:, i])
                n_sb += 1

        # sort traces by population
        # if len(cur_labels) > 0:
        #     cur_map, cur_labels = (list(t) for t in zip(*sorted(zip(cur_map, cur_labels), key=lambda x: np.sum(x[0]), reverse=True)))

        # sort traces by time of first appearance
        if len(cur_labels) > 0:
            cur_map, cur_labels = (list(t) for t in
                                   zip(*sorted(zip(cur_map, cur_labels), key=lambda x: np.nonzero(x[0])[0][0])))

        # print information about number of HB and SB
        msg_fmt = "{}:\tHB={}\tSB={}"
        print(msg_fmt.format(resid, n_hb, n_sb))

        labels.append(cur_labels)
        cmap.append(cur_map)

    # determine number of lines for each residue
    numbers_of_lines = np.zeros(n_tail_resid, dtype=int)

    for i, res_contacts in enumerate(labels):
        numbers_of_lines[i] = len(res_contacts)

    print("Total number of rows: %d" % np.sum(numbers_of_lines))

    return cmap, labels, numbers_of_lines


def build_pixel_map(cmap, labels, n_frames, numbers_of_lines):
    rows = np.sum(numbers_of_lines)

    hb_pixel_map = np.zeros((n_frames, rows))
    sb_pixel_map = np.zeros((n_frames, rows))
    tcks = [''] * rows

    for i, traces in enumerate(cmap):
        j = 0
        for label, trace in zip(labels[i], traces):
            if label[0][2] == "HB":
                tcks[np.sum(numbers_of_lines[:i]) + j] = label_formatter(d_name=label[0][0][1],
                                                                         d_id=label[0][0][2],
                                                                         d_atom=label[0][0][3],
                                                                         d_chain=label[0][0][0],
                                                                         a_name=label[0][1][1],
                                                                         a_id=label[0][1][2],
                                                                         a_atom=label[0][1][3],
                                                                         a_chain=label[0][1][0])
                hb_pixel_map[:, np.sum(numbers_of_lines[:i]) + j] = trace
            elif label[0][2] == "SB":
                tcks[np.sum(numbers_of_lines[:i]) + j] = label_formatter(d_name=label[0][0][1],
                                                                         d_id=label[0][0][2],
                                                                         d_atom="",
                                                                         d_chain=label[0][0][0],
                                                                         a_name=label[0][1][1],
                                                                         a_id=label[0][1][2],
                                                                         a_atom="",
                                                                         a_chain=label[0][1][0])
                sb_pixel_map[:, np.sum(numbers_of_lines[:i]) + j] = trace
            j += 1

    return hb_pixel_map, sb_pixel_map, tcks


def three_to_one(three_letter_resname):
    map_3_1 = {"ALA": "A",
               "ARG": "R",
               "ASN": "N",
               "ASP": "D",
               "ASH": "D",
               "CYS": "C",
               "GLU": "E",
               "GLH": "E",
               "GLN": "Q",
               "GLY": "G",
               "HIS": "H",
               "HIE": "H",
               "HID": "H",
               "HIP": "H",
               "ILE": "I",
               "LEU": "L",
               "LYS": "K",
               "LYN": "K",
               "MET": "M",
               "PHE": "F",
               "PRO": "P",
               "SER": "S",
               "THR": "T",
               "TRP": "W",
               "TYR": "Y",
               "VAL": "V",
               "DA": "DA",
               "DA3": "DA",
               "DA5": "DA",
               "DT": "DT",
               "DT3": "DT",
               "DT5": "DT",
               "DG": "DG",
               "DG3": "DG",
               "DG5": "DG",
               "DC": "DC",
               "DC3": "DC",
               "DC5": "DC",
               "Na+": "Na+"}
    if len(three_letter_resname) > 1:
        return map_3_1[three_letter_resname]
    else:
        return three_letter_resname


def label_formatter(d_name, d_id, d_atom, d_chain, a_name, a_id, a_atom, a_chain):
    lbl_fmt_don = '{}{:d}{} - {}{:d}{} ({})'
    if ((d_chain == "B") & (a_chain != "B")) | ((d_chain == "F") & (a_chain != "F")):
        return lbl_fmt_don.format(three_to_one(d_name),
                                  d_id,
                                  (" " + d_atom) if (len(d_atom) > 0) else "",
                                  three_to_one(a_name),
                                  a_id,
                                  (" " + a_atom) if (len(a_atom) > 0) else "",
                                  a_chain)
    elif (((d_chain == "B") & (a_chain == "B")) | ((d_chain == "F") & (a_chain == "F"))) & (d_id < 25):
        return lbl_fmt_don.format(three_to_one(d_name),
                                  d_id,
                                  (" " + d_atom) if (len(d_atom) > 0) else "",
                                  three_to_one(a_name),
                                  a_id,
                                  (" " + a_atom) if (len(a_atom) > 0) else "",
                                  a_chain)
    else:
        return lbl_fmt_don.format(three_to_one(a_name),
                                  a_id,
                                  (" " + a_atom) if (len(a_atom) > 0) else "",
                                  three_to_one(d_name),
                                  d_id,
                                  (" " + d_atom) if (len(d_atom) > 0) else "",
                                  d_chain)