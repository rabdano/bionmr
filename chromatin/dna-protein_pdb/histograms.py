import matplotlib.pyplot as plt
import numpy as np

fn_lys_po4 = "lys-po4.csv"
fn_nh3_po4 = "nh3-po4.csv"
fn_arg_po4 = "arg-po4.csv"

dist_start = 0
dist_end = 10
dist_nbin = 100

angle_start = 0
angle_end = 180
angle_nbin = 180

for fn in [fn_lys_po4, fn_nh3_po4, fn_arg_po4]:

    data = np.genfromtxt(fn, delimiter=",", skip_header=1)
    with open(fn, "r") as f:
        header = f.readline().rstrip("\n").lstrip("# ").split(",")
    header = [h.strip()[4:] for h in header]
    print(header)

    print("Pairs in %s: %d" % (fn, len(data)))

    # calculate volumes of bins for normalization
    dist_bin_vol = np.ones((dist_nbin,))
    for i in range(dist_nbin):
        r1 = dist_start + i * (dist_end - dist_start) / dist_nbin
        r2 = dist_start + (i + 1) * (dist_end - dist_start) / dist_nbin
        dist_bin_vol[i] = 4.0/3.0 * np.pi * r2**3 - 4.0/3.0 * np.pi * r1**3

    angle_bin_vol = np.ones((angle_nbin,))
    for i in range(angle_nbin):
        th1 = angle_start + i * (angle_end - angle_start) / angle_nbin
        th2 = angle_start + (i + 1) * (angle_end - angle_start) / angle_nbin
        angle_bin_vol[i] = 2.0 * np.pi * np.sin(0.5 * (th2 + th1) * np.pi / 180.0) * (th2 - th1) * np.pi / 180.0

    dist_angle_bin_vol = np.ones((dist_nbin, angle_nbin))
    for i in range(dist_nbin):
        r1 = dist_start + i * (dist_end - dist_start) / dist_nbin
        r2 = dist_start + (i + 1) * (dist_end - dist_start) / dist_nbin
        for j in range(angle_nbin):
            th1 = angle_start + j * (angle_end - angle_start) / angle_nbin
            th2 = angle_start + (j + 1) * (angle_end - angle_start) / angle_nbin

            dist_angle_bin_vol[i, j] = 2.0 * np.pi * (0.5 * (r2 + r1))**2 * np.sin(0.5 * (th2 + th1) * np.pi / 180.0) * (r2 - r1) * (th2 - th1) * np.pi / 180.0

    for i, label, x in zip(range(len(data.T)), header, data.T):
        plt.figure()
        if np.mean(x) < 20:
            hist, bin_edges = np.histogram(x, bins=np.linspace(dist_start, dist_end, dist_nbin + 1))
            # print(hist)
            hist = np.divide(hist, dist_bin_vol) / len(x)
            # print(hist)
            bar_x = np.array([0.5 * (x1 + x2) for x1, x2 in zip(bin_edges[1:], bin_edges)])
            # plt.hist(x, bins=np.linspace(dist_start, dist_end, dist_nbin))
            plt.hist(bin_edges[:-1], bin_edges, weights=hist)
            if len(label.split("_")) == 2:
                plt.xlabel(r"dist %s, $\rm\AA$" % label)
            elif len(label.split("_")) == 3:
                plt.xlabel(r"angle %s, deg" % label)
            plt.xlim(dist_start, dist_end)
            plt.gca().set_xticks(range(dist_start, dist_end + 1))

        elif np.mean(x) > 20.0:
            hist, bin_edges = np.histogram(x, bins=np.linspace(angle_start, angle_end, angle_nbin + 1))
            hist = np.divide(hist, angle_bin_vol) / len(x)
            bar_x = np.array([0.5 * (x1 + x2) for x1, x2 in zip(bin_edges[1:], bin_edges)])
            # plt.hist(x, bins=np.linspace(angle_start, angle_end, angle_nbin))
            # plt.bar(bar_x, hist)
            plt.hist(bin_edges[:-1], bin_edges, weights=hist)
            if len(label.split("_")) == 2:
                plt.xlabel(r"dist %s, $\rm\AA$" % label)
            elif len(label.split("_")) == 3:
                plt.xlabel(r"angle %s, deg" % label)
            plt.xlim(angle_start, angle_end)
        plt.ylabel('Count')
        plt.title(label)
        plt.savefig("Figures/" + fn[:-4] + "_" + label + ".png", bbox_inches="tight")
        plt.close()

    # plot 2D
    fig, ax = plt.subplots()
    x_r = data.T[0]
    hist_r, bin_edges_r = np.histogram(x_r, bins=np.linspace(dist_start, dist_end, dist_nbin + 1))
    hist_r = np.divide(hist_r, dist_bin_vol) / len(x_r)
    x_a = data.T[1]
    hist_a, bin_edges_a = np.histogram(x_a, bins=np.linspace(angle_start, angle_end, angle_nbin + 1))
    hist_a = np.divide(hist_a, angle_bin_vol) / len(x_a)

    weights = []
    for (d, a) in zip(x_r, x_a):
        idx_r = 0
        idx_a = 0
        for i, edge_r in enumerate(bin_edges_r):
            if edge_r > d:
                idx_r = i - 1
                break
        for i, edge_a in enumerate(bin_edges_a):
            if edge_a > a:
                idx_a = i - 1
                break
        weights.append(1.0 / dist_angle_bin_vol[idx_r, idx_a] / len(data))

    if len(header[0].split("_")) == 2:
        x_max = 10
    elif len(header[0].split("_")) == 3:
        x_max = 180
    if len(header[1].split("_")) == 2:
        y_max = 10
    elif len(header[1].split("_")) == 3:
        y_max = 180

    counts, xedges, yedges, im = ax.hist2d(x_r, x_a, bins=[100, 180], range=[[0, x_max], [0, y_max]], weights=weights)

    if len(header[0].split("_")) == 2:
        plt.xlabel(r"dist %s, $\rm\AA$" % header[0])
        plt.xlim([0, 10])
    elif len(header[0].split("_")) == 3:
        plt.xlabel(r"angle %s, deg" % header[0])
        plt.xlim([0, 180])

    if len(header[1].split("_")) == 2:
        plt.ylabel(r"dist %s, $\rm\AA$" % header[1])
        plt.ylim([0, 10])
    elif len(header[1].split("_")) == 3:
        plt.ylabel(r"angle %s, deg" % header[1])
        plt.ylim([0, 180])

    plt.colorbar(im, ax=ax)
    plt.savefig("Figures/" + fn[:-4] + "_" + header[0] + "-" + header[1] + "_2D.png")
