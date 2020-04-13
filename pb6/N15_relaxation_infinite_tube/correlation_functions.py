import sys
from glob import glob
import numpy as np
from scipy.optimize import curve_fit


class Corfun:
    def __init__(self, path=None):
        self.data = None
        self.info = None
        self.resi = None
        self.pars = None
        self.s2s = None
        self.rmsd = None
        self.R1 = None
        self.R2 = None
        if path:
            self.read(path)

    def read(self, path):
        self.data = []
        self.info = []
        self.resi = []
        fns = [c.split('/')[-1] for c in sorted(glob(path + '/*.cor'))]
        fns.sort()
        for fn in fns:
            fds = fn.split('.')[0].split('_')
            if not fds[0].isdigit():
                fds = fds[1:]  # remove possible segid
            resi = int(fds[0])
            info = '-'.join(fds)
            data = np.loadtxt(path+'/'+fn)
            self.data.append(data)
            self.info.append(info)
            self.resi.append(resi)

    def write(self, path):
        for acf, res, inf in zip(self.data, self.resi, self.info):
            resname = inf.split("-")[1]
            fname = path + "/%04d_%3s.cor" % (res, resname)
            np.savetxt(fname, acf)

    @staticmethod
    def multiexponent(x, *args):
        N = int(len(args) / 2)
        amps, taus = args[:N], args[N:]
        assert len(amps) == len(taus), (amps, taus)
        res = np.zeros((N, len(x)))
        for i in range(N):
            res[i, :] = amps[i] * np.exp(-x / taus[i])
        return np.sum(res, axis=0)

    @staticmethod
    def multiexponent_s2(s2):
        def multiexponent_s2_(x, *args):
            N = int(len(args) / 2)
            amps, taus = args[:N], args[N:]
            assert len(amps) == len(taus), (amps, taus)
            res = np.zeros((N, len(x)))
            for i in range(N):
                res[i, :] = amps[i] * np.exp(-x / taus[i])
            return np.sum(res, axis=0) + s2

        return multiexponent_s2_

    @staticmethod
    def J(w, parms):
        N = int(len(parms)/2)
        amps, taus = parms[:N], parms[N:]
        assert len(amps) == len(taus)
        res = np.zeros((N,))
        for i in range(N):
            # 2.0 comes from "J(w)=2*FT(Corfun)"
            res[i] = 2.0/5.0*amps[i]*taus[i]/(1+(w*taus[i])**2)
        return np.sum(res, axis=0)

    def fit(self, n_exp=6, rep=5, s2=False, step=1, sp=None, lb=None, ub=None, fout=None):
        self.pars = []
        self.rmsd = []

        def randminmax(min, max):
            return np.random.rand(1) * (max - min) + min

        for M, info, res in zip(self.data, self.info, self.resi):
            xdata = M[:, 0]
            ydata = M[:, 1]
            normf = ydata[0]  # normalization factor
            ydata = ydata / normf

            rmsd_min = np.inf
            pars_min = [0] * n_exp * 2
            pars = pars_min

            for i in range(rep):
                sys.stdout.write('> fitting %s: %02d/%02d ' % (info, i+1, rep))
                guess = np.zeros(n_exp * 2)
                guess[:n_exp] = sp[:n_exp]
                for j in range(n_exp - 1):
                    guess[j + n_exp] = 10 ** (randminmax(np.log10(sp[j + n_exp]), np.log10(sp[j + n_exp + 1])))
                guess[-1] = 10 ** (randminmax(np.log10(sp[-1]), np.log10(sp[-1] ** 2 / sp[-2])))
                guess = guess.tolist()

                try:
                    if not s2:
                        popt, pcov = curve_fit(self.multiexponent,
                                               xdata, ydata,
                                               p0=guess,
                                               bounds=(lb, ub),
                                               max_nfev=10000,
                                               ftol=1e-6)
                        pars = [p for p in popt]
                        rmsd = np.sqrt(np.mean(np.power(np.subtract(ydata, self.multiexponent(xdata, *popt)), 2)))
                    else:
                        popt, pcov = curve_fit(self.multiexponent_s2(self.s2s[self.resi.index(res)]),
                                               xdata, ydata,
                                               p0=guess,
                                               bounds=(lb, ub),
                                               max_nfev=10000,
                                               ftol=1e-6)
                        pars = [p for p in popt]
                        rmsd = np.sqrt(np.mean(np.power(np.subtract(ydata, self.multiexponent(xdata, *popt)), 2)))
                except RuntimeError:
                    sys.stdout.write('[bad rep #%d] ' % i)


                if rmsd < rmsd_min:
                    rmsd_min = rmsd
                    pars_min = pars
                    idx_min = i

                sys.stdout.write('rmsd(%02d)=%.5e\n' % (idx_min + 1, rmsd_min))
                sys.stdout.write('\033[1A')  # cursor up by one line

            sys.stdout.write('\n')

            # recover normalization
            for i in range(n_exp):
                pars_min[i] = pars_min[i] * normf

            # convert tau units from points to seconds
            for i in range(n_exp):
                pars_min[i + n_exp] = pars_min[i + n_exp] * step

            self.pars.append(pars_min)
            self.rmsd.append(rmsd_min)

        if fout:
            f = open(fout, 'wt')
            for info, pars, rmsd in zip(self.info, self.pars, self.rmsd):
                n = int(len(pars) / 2)
                line = '%s ' % info + ('%.4e ' * n) % tuple(pars[:n]) + '  ' + \
                       ('%.4e ' * n) % tuple(pars[n:]) + ' %.3e\n' % rmsd
                f.write(line)
            f.close()

    def calc_R1(self, opfreq=600e6, X='N15',
                CSA=None, rHX=None, fnout=None):
        self.R1 = []
        for pars in self.pars:
            R1 = self.R1_fun(pars, opfreq, X, CSA, rHX)
            self.R1.append(R1)

        if fnout:
            f = open(fnout,'wt')
            for info,R1 in zip(self.info, self.R1):
                f.write('%s  %8.5e\n'%(info,R1))
            f.close()

    def calc_R2(self, opfreq=600e6, X='N15',
                CSA=None, rHX=None, fnout=None):
        self.R2 = []
        for pars in self.pars:
            R2 = self.R2_fun(pars, opfreq, X, CSA, rHX)
            self.R2.append(R2)

        if fnout:
            f = open(fnout,'wt')
            for info,R2 in zip(self.info, self.R2):
                f.write('%s  %8.5e\n'%(info,R2))
            f.close()

    def R1_fun(self, parms, opfreq=600e6, X='N15', CSA=None, rHX=None):
        # calculate R1 for fix-length inter-spin vector
        # INPUT: amps --- amplitute of the associated taus
        #        taus --- correlation time components
        #        opfreq in unit of Hz

        # Calculate dipolar coupling constant d2 -------------------------------
        u0 = 4 * np.pi * 1e-7
        h = 6.626069e-34
        gH = 267.522e6
        if X == 'N15':
            gX = -27.126e6
            if CSA == None: CSA = -172.0  # unit: ppm
            if rHX == None: rHX = 1.02e-10
        else:
            return None
        d2 = ((u0 / 4 / np.pi) * (h / 2 / np.pi) * (gH * gX) / (rHX ** 3)) ** 2
        # -----------------------------------------------------------------------
        gHX = gH / gX  # ratio gamma_H to gamma_X
        wH = 2 * np.pi * opfreq
        wX = wH / gHX
        c2 = (1. / 3.) * ((CSA * 1e-6 * wX) ** 2)

        R1 = 0.25 * d2 * (3 * self.J(wX, parms) + self.J(wH - wX, parms) + 6 * self.J(wH + wX, parms)) \
                  + c2 * self.J(wX, parms)

        return R1

    def R2_fun(self, parms, opfreq=600e6, X='N15', CSA=None, rHX=None):
        # calculate R1 for fix-length inter-spin vector
        # INPUT: amps --- amplitute of the associated taus
        #        taus --- correlation time components
        #        opfreq in unit of Hz

        # Calculate dipolar coupling constant d2 -------------------------------
        u0 = 4 * np.pi * 1e-7
        h = 6.626069e-34
        gH = 267.522e6
        if X == 'N15':
            gX = -27.126e6
            if CSA == None: CSA = -172.0  # unit: ppm
            if rHX == None: rHX = 1.02e-10
        else:
            return None
        d2 = ((u0 / 4 / np.pi) * (h / 2 / np.pi) * (gH * gX) / (rHX ** 3)) ** 2
        # -----------------------------------------------------------------------
        gHX = gH / gX  # ratio gamma_H to gamma_X
        wH = 2 * np.pi * opfreq
        wX = wH / gHX
        c2 = (1. / 3.) * ((CSA * 1e-6 * wX) ** 2)

        R2 = 0.125 * d2 * (4 * self.J(0, parms) + 3 * self.J(wX, parms) + self.J(wH - wX, parms)
                           + 6 * self.J(wH, parms) + 6 * self.J(wH + wX, parms)) + \
             (1. / 6.) * c2 * (4 * self.J(0, parms) + 3 * self.J(wX, parms))

        return R2


def cor(vec, grid, skip_frames):
    vectors = np.zeros((len(vec), 3))
    for i in range(len(vec)):
        vectors[i] = vec[i].to_np / np.linalg.norm(vec[i].to_np)
    acf = np.zeros(len(grid))
    for k, lag in enumerate(grid):
        res = []
        for i, j in zip(range(len(vectors)), range(lag, len(vectors), 1)):
            if i in skip_frames:
                continue
            if j in skip_frames:
                continue
            res.append((3.0 * np.dot(vectors[i], vectors[j])**2 - 1) / 2.0)
        acf[k] = np.mean(np.array(res))
    return acf


def gets2(vec):
    # vec - np.array of vectors
    n = len(vec)
    for i, ri in enumerate(vec):
        norm = np.linalg.norm(ri)
        if norm != 0:
            vec[i] = ri / norm
    xx = vec[:, 0].dot(vec[:, 0]) / n
    yy = vec[:, 1].dot(vec[:, 1]) / n
    zz = vec[:, 2].dot(vec[:, 2]) / n
    xy = vec[:, 0].dot(vec[:, 1]) / n
    xz = vec[:, 0].dot(vec[:, 2]) / n
    yz = vec[:, 1].dot(vec[:, 2]) / n

    return 0.25 * (3 * zz - 1) * (3 * zz - 1) + 3 * (xy * xy + xz * xz + yz * yz) + 0.75 * (xx - yy) * (xx - yy)
