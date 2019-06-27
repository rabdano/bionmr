import numpy as np

time, V = np.genfromtxt('summary.VOLUME', unpack=True)
time = time[999::1000]
V = V[999::1000]

inpcrd = '../1_build/box.inpcrd'

with open(inpcrd) as f:
    content = f.readlines()

content = content[-1].split()[:6]

a, b, c, alpha, beta, gamma = [float(x) for x in content]

# https://en.wikipedia.org/wiki/Parallelepiped
V0 = a*b*c*np.sqrt(1.0+2.0*np.cos(alpha/180*np.pi)*np.cos(beta/180*np.pi)*np.cos(gamma/180*np.pi)-np.cos(alpha/180*np.pi)**2-np.cos(beta/180*np.pi)**2-np.cos(gamma/180*np.pi)**2)

print(a, b, c, alpha, beta, gamma)
print(V0)

# XST file format
# http://www.ks.uiuc.edu/Research/namd/mailing_list/namd-l.2003-2004/0234.html

m_v1 = np.array([a, 0.0, 0.0])
m_v2 = np.array([b*np.cos(gamma/180*np.pi), b*np.sin(gamma/180*np.pi), 0])
m_v3 = np.zeros(3)

m_v3[0] = c*np.cos(beta/180*np.pi);
m_v3[1] = b/m_v2[1] * c * np.cos(alpha/180*np.pi) - m_v2[0]/m_v2[1]*m_v3[0];
m_v3[2] = np.sqrt(c*c - m_v3[0]*m_v3[0] - m_v3[1]*m_v3[1]);

# print(c*c)
# print(m_v3[0]*m_v3[0])
# print(m_v3[1]*m_v3[1])
# print(c*c - m_v3[0]*m_v3[0] - m_v3[1]*m_v3[1])
# print(m_v1, m_v2, m_v3)

xst = np.zeros((len(V), 13))
origin = np.zeros(3)

for i,f_V in enumerate(V):
    scaling_factor = np.cbrt(f_V / V0)
    # print(scaling_factor)
    xst[i, 0] = time[i]
    xst[i, 1:10] = np.array([m_v1[0], m_v1[1], m_v1[2], m_v2[0], m_v2[1], m_v2[2], m_v3[0], m_v3[1], m_v3[2]]) * scaling_factor
    xst[i, 10:13] = [origin[0], origin[1], origin[2]]

np.savetxt('traj.xst', xst, '%.5f')
