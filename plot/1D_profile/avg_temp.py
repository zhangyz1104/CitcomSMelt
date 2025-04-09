import numpy as np
import matplotlib.pyplot as plt

''' set parameters '''
# cases = [('case1_ref', 40000, 'Case1_ref', 'k'),
#          ('debug_ref', 30000, 'comp. conservation', 'r'),
#          ('debug_compact3', 40000, 'compaction', 'b')]
# cases = [('case1_ref', 16000, 'Case1_ref'),
#          ('case1_ref_2', 35000, 'no melt extraction'),
#          ('case1_ref_3', 10000, 'modified melt extraction'),
#          ('case1_ref_4', 16000, 'comp. conservation')]
# cases = [('case1_ref', 40000, 'case1_ref', 25800),
#          ('case2_lowVis', 98000, 'case2_lowVis', 46000),
#          ('case3_highVis', 22000, 'case3_highVis', 12800),
#          ('case4_lowSegr_lowAcc', 60000, 'case4_lowSegr', 48400),
#          ('case5_highSegr', 28000, 'case5_highSegr', 15600),
#          ('case6_lowExtr', 34000, 'case6_lowExtr', 23400),
#          ('case7_highExtr', 38000, 'case7_highExtr', 24000),
#          ('case8_shallowExtr_lowAcc', 60000, 'case8_shallowExtr', 22400),
#          ('case9_deepExtr', 38000, 'case9_deepExtr', 24200),
#          ('case12_noRad', 33400, 'Case11_noRad', 24800),
#          ('case10_homoRad', 35000, 'Case12_highRad', 25800),
#          ('case11_inhomoRad', 32200, 'Case13_inhomoRad', 24800),
#          ('case13_constDens', 63800, 'Case2_dens', 26400)]
# cases = [('case1_ref', 40000, 'ref'),
#          ('case2_lowVis', 70000, 'Case3_lowVis'),
#          ('case3_highVis', 20000, 'Case4_highVis')]
# cases = [('case1_ref', 40000, 'Case1_ref', 35800),
#          ('case1_ref_3', 50000, 'Case1_ref_3', 35800),
#          ('case_surf_dike3', 40000, 'Case_surf_dike3', 35800),
#          ('case_surf_dike4', 50000, 'Case_surf_dike4', 35800)]
cases = [('case1_ref', 40000, 'Case1_ref', 35800),
         ('case1_ref_3', 50000, 'Case1_ref_3', 35800),
         ('case_full', 50000, 'Case_surf_dike', 35800),
         ('case_surf_dike3', 40000, 'Case_noExtr', 35800),
         ('case_noSegr2', 51000, 'Case_noSegr', 35800),
         ('case_noSolid', 33000, 'Case_noSolid', 35800),
         ('case_Segr3', 29000, 'Case_Segr2', 35800),
         ('case_Segr6', 14000, 'Case_Segr3', 35800)]
cases = [('case1_ref_3', 50000, 'Case1_ref_3', 35800)]

# 0.128  0.122  0.133

''' dims 
[0: r, 1: T, 2: Vh, 3: Vv, 4: eta, 5: phi, 6: c0, 7: cl, 8: rho, 9: q, 10: E, 11: M]
'''
[r0, T0, v0, eta0, t0, drho0, rho0] = [1740, 2104, 1.81e-3, 5e20, 9.59e4, 500, 3340]
avg0 = [r0, T0, v0, v0, eta0, 1, 1, 1, -drho0, v0, 1 / t0, 1 / t0]

T1 = []
phi1 = []
c01 = []
E1 = []
T = []
phi = []
c0 = []
E = []
time = []

for case in cases:
    filename = 'D:/MUSH_result/' + case[0] + '.time'
    time0 = np.loadtxt(filename).transpose()[1] * 9.59e4

    T1.append([])
    phi1.append([])
    c01.append([])
    E1.append([])
    T.append([])
    phi.append([])
    c0.append([])
    E.append([])
    time.append([])
    r = []
    for cycle in range(0, case[1], 200):
        filename0 = 'D:/MUSH_result/' + case[0] + '.horiz_avg.0.' + str(cycle)
        filename1 = 'D:/MUSH_result/' + case[0] + '.horiz_avg.1.' + str(cycle)
        rawdata = np.vstack((np.loadtxt(filename0), np.loadtxt(filename1)))
        avg = rawdata.transpose()
        for j in range(len(avg0)):
            avg[j] *= avg0[j]

        r = avg[0]
        _T = 0
        _phi = 0
        _c0 = 0
        _E = 0
        _V = 0
        # for i in range(1, len(r) - 14):
        #     dr = r[i] - r[i - 1]
        #     dV = 4 * np.pi * r[i] ** 2 * dr
        #     _T += avg[1][i] * dV
        #     _phi += avg[5][i] * dV
        #     _c0 += avg[6][i] * dV
        #     _E += avg[10][i] * dV
        #     _V += dV
        #
        # T1[-1].append(_T / _V)
        # phi1[-1].append(_phi / _V)
        # c01[-1].append(_c0 / _V)
        # E1[-1].append(_E / _V)
        #
        # for i in range(len(r) - 14, len(r)):
        #     dr = r[i] - r[i - 1]
        #     dV = 4 * np.pi * r[i] ** 2 * dr
        #     _T += avg[1][i] * dV
        #     _phi += avg[5][i] * dV
        #     _c0 += avg[6][i] * dV
        #     _E += avg[10][i] * dV
        #     _V += dV
        for i in range(len(r) - 1):
            dr = r[i + 1] - r[i]
            dV = (r[i + 1] + r[i]) ** 2 * dr
            # dV = 1
            _T += avg[1][i] * dV
            _phi += avg[5][i] * dV
            _c0 += avg[6][i] * dV
            _E += avg[10][i] * dV
            _V += dV

        T[-1].append(_T / _V)
        phi[-1].append(_phi / _V)
        _c0 = _c0 / _V
        _c0 += 0.5 * (0.15 - _c0)
        c0[-1].append(_c0)
        E[-1].append(_E / _V)
        time[-1].append(time0[cycle])

        if cycle == case[3]:
            print(_c0)
    c0[-1][0] = 0.15

''' plot '''
plt.figure(figsize=(6, 4))
plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams['font.size'] = 12

ax1 = plt.subplot()
for i in range(len(cases)):
    ax1.plot(time[i], T[i], linewidth=2, label=cases[i][2], c='k')
    ax1.plot(time[i], phi[i], linewidth=2, label=cases[i][2], c='b')
# plt.xscale('log')
# plt.ylim(0, 2)
ax1.set_xlim(1e-3, 1e3)
ax1.set_xscale('log')
ax1.set_xlabel('Time (Myr)')
ax1.set_ylim(1600, 1850)
ax1.set_ylabel('Global average temperature')
plt.tight_layout()
# ax1.legend()
plt.show()

