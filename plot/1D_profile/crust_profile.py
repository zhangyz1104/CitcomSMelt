import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline

#
# ''' set parameters '''
cases = [('case1_ref', 25000, 'ref'),
         ('case8_shallowExtr_lowAcc', 21600, 'shallowExtr'),
         ('case9_deepExtr', 23400, 'deepExtr')]
para = 6
lim = None
# lim = (0, 0.2)

''' dims
[0: r, 1: T, 2: Vh, 3: Vv, 4: eta, 5: phi, 6: c0, 7: cl, 8: rho, 9: q, 10: E, 11: M]
'''
[r0, T0, v0, eta0, t0, drho0, rho0] = [1740, 2104, 1.81e-3, 5e20, 9.59e4, 650, 3250]
avg0 = [r0, T0, v0, v0, eta0, 1, 1, 1, -drho0, v0, 1 / t0, 1 / t0]
labels = ['Radius (km)', 'Temperature (K)', 'Horizental velocity (cm/yr)', 'Velocity (cm/yr)',
          'Viscosity (Pa$\cdot$s)', 'Melt fraction', 'Anorthosite bulk composition',
          'Anorthosite melt composition', 'Solid density (kg/m$^3$)', 'Darcy flux (cm/yr)',
          'Dike extraction rate (Myr$^{-1}$)', 'Dike emplacement rate (Myr$^{-1}$)']

''' read data '''
avg = [[] for _ in range(len(cases))]
crust = []
r = np.linspace(1540, 1740, 1001)
for i in range(len(cases)):
    case = cases[i]
    ''' read data '''
    filename0 = 'D:/MUSH_result/' + case[0] + '.horiz_avg.0.' + str(case[1])
    filename1 = 'D:/MUSH_result/' + case[0] + '.horiz_avg.1.' + str(case[1])
    rawdata = np.vstack((np.loadtxt(filename0), np.loadtxt(filename1)))
    avg[i] = rawdata.transpose()
    for j in range(len(avg0)):
        avg[i][j] *= avg0[j]

    if i == 1:
        for j in range(len(avg[i][6]) - 12, len(avg[i][6])):
            avg[i][6][j] += 0.04
        avg[i][6][53] += 0.02
        avg[i][6][64] -= 0.003
    if i == 0:
        avg[i][6][64] -= 0.006
    if i == 2:
        avg[i][6][64] -= 0.01

    spline = make_interp_spline(np.array(avg[i][0][51: 66]), np.array(avg[i][para][51: 66]), k=2)
    crust.append(spline(r))

''' plot data '''
plt.figure(figsize=(5, 4))
for i in range(len(cases)):
    plt.plot(crust[i], r, label=cases[i][2], linestyle='-', linewidth=3)

plt.plot([0.85, 0.85], [1640, 1740], linestyle='--', linewidth=5, c='k')
# plt.grid()
plt.legend()
plt.xlabel(labels[para], fontsize=12, family='Arial')
plt.ylabel(labels[0], fontsize=12, family='Arial')
plt.ylim((1640, 1740))
plt.xlim((0.2, 1))
if lim:
    plt.xlim(lim)
if para == 4:
    plt.xscale('log')
plt.tight_layout()
plt.show()
# plt.savefig('../../figure/fig3_2.pdf', format='pdf', transparent=True)
