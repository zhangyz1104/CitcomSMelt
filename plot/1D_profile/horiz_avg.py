import matplotlib.pyplot as plt
import numpy as np

#
# ''' set parameters '''
cases = [('case1_ref_3', 9200, 'E=100 KJ/mol', 'k'),
         ('case_viscE200_3', 17600, 'E=200 KJ/mol (1)', 'b'),
         ('case_viscE200_5', 3000, 'E=200 KJ/mol (2)', 'g')] # init comp.
cases = [('case1_ref_3', 22000, 'E=100 KJ/mol', 'k'),
         ('case_viscE200_5', 22800, 'E=200 KJ/mol (2)', 'g')] # init comp.
cases = [('case1_ref_3', 20000, 'E=100 KJ/mol', 'k'),
         ('case_surf_dike2', 10800, 'E=200 KJ/mol (2)', 'g')] # init comp.
cases = [('case1_ref_3', 12400, 'ref', 'k'),
         ('case_noSegr', 10800, 'noSegr', 'r'),
         ('case_noSolid', 10200, 'noSolid', 'g'),
         ('case_full', 11600, 'full', 'b')]
cases = [('case1_ref_3', 24600, 'ref', 'k'),
         ('case_Segr3', 13800, 'Segr2', 'b'),
         ('case_Segr6', 20000, 'Segr3', 'g')]
cases = [('case_mesh_r', 0, 'ref', 'k'),
         ('case_mesh_r', 4000, 'Segr2', 'b')]
# cases = [('case_surf_dike', 11000, 'E=200', 'b'),
#          ('case1_ref_3', 13000, 'E=100', 'k')] # init comp.
# cases = [('case1_ref_3', 0, '$\phi_0=$0.2', 'k'),
#          ('case_init_phi0.3_2', 800, '$\phi_0=$0.3', 'b'),
#          ('case_init_phi0.4_2', 800, '$\phi_0=$0.4', 'g'),
#          ('case_init_phi0.5_2', 1000, '$\phi_0=$0.5', 'r')] # 1 kyr (phi<0.5)
# cases = [('case1_ref_3', 0, '$\phi_0=$0.2', 'k'),
#          ('case_init_phi0.3_2', 800, '$\phi_0=$0.3', 'b'),
#          ('case_init_phi0.4_3', 4400, '$\phi_0=$0.4', 'g'),
#          ('case_init_phi0.5_3', 5400, '$\phi_0=$0.5', 'r')] # 1 kyr (eta_min=1e-7)
# cases = [('case1_ref_3', 1000, '$\phi_0=$0.2', 'k'),
#          ('case_init_phi0.3_2', 2800, '$\phi_0=$0.3', 'b'),
#          ('case_init_phi0.4_2', 2200, '$\phi_0=$0.4', 'g'),
#          ('case_init_phi0.5_2', 2400, '$\phi_0=$0.5', 'r')] # 10 kyr (phi<0.5)
# cases = [('case1_ref_3', 1000, '$\phi_0=$0.2', 'k'),
#          ('case_init_phi0.3_2', 2800, '$\phi_0=$0.3', 'b'),
#          ('case_init_phi0.4_3', 9200, '$\phi_0=$0.4', 'g'),
#          ('case_init_phi0.5_3', 13000, '$\phi_0=$0.5', 'r')] # 10 kyr (eta_min=1e-7)
# cases = [('case1_ref_3', 3400, '$\phi_0=$0.2', 'k'),
#          ('case_init_phi0.3_2', 6600, '$\phi_0=$0.3', 'b'),
#          ('case_init_phi0.4_2', 6800, '$\phi_0=$0.4', 'g'),
#          ('case_init_phi0.5_2', 8000, '$\phi_0=$0.5', 'r')] # 100 kyr (phi<0.5)
# cases = [('case1_ref_3', 9000, '$\phi_0=$0.2', 'k'),
#          ('case_init_phi0.3_2', 14200, '$\phi_0=$0.3', 'b'),
#          ('case_init_phi0.4_2', 18000, '$\phi_0=$0.4', 'g'),
#          ('case_init_phi0.5_2', 21400, '$\phi_0=$0.5', 'r')] # 1 Myr (phi<0.5)
# cases = [('case1_ref', 1000, '20 kyr', 'r'),
#          ('case1_ref', 1800, '40 kyr', 'orange'),
#          ('case1_ref', 2400, '60 kyr', 'yellow'),
#          ('case1_ref', 2800, '80 kyr', 'g'),
#          ('case1_ref', 3400, '100 kyr', 'b')]
# cases = [('case3_highVis', 0, '0 kyr', 'k'),
#          ('case1_ref', 4800, '200 kyr', 'r'),
#          ('case1_ref', 6600, '400 kyr', 'orange'),
#          ('case1_ref', 8200, '600 kyr', 'yellow'),
#          ('case1_ref', 9200, '800 kyr', 'g'),
#          ('case1_ref', 10400, '1000 kyr', 'b')]
# cases = [('case1_ref', 2200, '10 kyr', 'k')]
# cases = [('case1_ref_4', 16000, 'comp cons', 'b'),
#          ('case1_ref_3', 15000, 'modf melt extr', 'g'),
#          ('case1_ref', 18600, 'ref', 'k')] # compostional conservation
# cases = [('case1_ref', 10400, 'k0=1e-10', 'k'),
#          ('case4_lowSegr_lowAcc', 13200, 'k0=1e-11', 'g'),
#          ('case5_highSegr', 7400, 'k0=1e-9', 'b')] # permeability (1 Myr)
# cases = [('case1_ref', 25600, 'Case1_ref', 'r'),
#          ('case5_highSegr', 15400, 'Case6_highSegr', 'g'),
#          ('case4_lowSegr_lowAcc', 48400, 'Case5_lowSegr', 'b')] # permeability (segr. timescale)
# cases = [('case1_ref', 23200, 'Case1_ref', 'r'),
#          ('case5_highSegr', 13800, 'Case6_highSegr', 'g'),
#          ('case4_lowSegr_lowAcc', 33600, 'Case5_lowSegr', 'b')] # permeability (10 Myr)
# cases = [('case1_ref', 25600, 'Case1_ref', 'r'),
#          ('case5_highSegr', 15400, 'Case6_highSegr', 'g'),
#          ('case4_lowSegr_lowAcc', 48400, 'Case5_lowSegr', 'b')] # extraction depth (segr. timescale)
# cases = [('case1_ref', 10000, 'Case1_ref', 'r'),
#          ('case8_shallowExtr_lowAcc', 10000, 'Case9_shallowExtr', 'g'),
#          ('case9_deepExtr', 10000, 'Case10_deepExtr', 'b')] # extraction depth (10 Myr)
# cases = [('case1_ref', 35600, 'k0=1e-10', 'k'),
#          ('case4_lowSegr_lowAcc', 46400, 'k0=1e-11', 'g'),
#          ('case5_highSegr', 25400, 'k0=1e-9', 'b')] # permeability (segr. timescale)
# cases = [('case_init_phi0.4_5', 0, 't=0', 'gray'),
#          ('case_init_phi0.4_5', 62000, 't=200', 'r'),
#          ('case_init_phi0.4_5', 97000, 't=400', 'g'),
#          ('case_init_phi0.4_5', 120000, 't=1000', 'b')]
para = 1
lim = None

''' dims
[0: r, 1: T, 2: Vh, 3: Vv, 4: eta, 5: phi, 6: c0, 7: cl, 8: rho, 9: q, 10: E, 11: M]
'''
[r0, T0, v0, eta0, t0, drho0, rho0] = [1740, 2104, 1.81e-3, 1e21, 9.59e4, 650, 3250]
avg0 = [r0, T0, v0, v0, eta0, 1, 1, 1, -drho0, v0, 1 / t0, 1 / t0]
labels = ['Radius (km)', 'Temperature (K)', 'Horizental velocity (cm/yr)', 'Velocity (cm/yr)',
          'Viscosity (Pa$\cdot$s)', 'Melt fraction', 'Anorthite bulk composition',
          'Anorthosite melt composition', 'Solid density (kg/m$^3$)', 'Darcy flux (cm/yr)',
          'Dike extraction rate (Myr$^{-1}$)', 'cooling rate (W/m$^3$)']

''' read data '''
avg = [[] for _ in range(len(cases))]
for i in range(len(cases)):
    case = cases[i]
    ''' read data '''
    filename0 = 'D:/MUSH_result/' + case[0] + '.horiz_avg.0.' + str(case[1])
    filename1 = 'D:/MUSH_result/' + case[0] + '.horiz_avg.1.' + str(case[1])
    rawdata = np.vstack((np.loadtxt(filename0), np.loadtxt(filename1)[1:]))
    avg[i] = rawdata.transpose()
    for j in range(len(avg0)):
        avg[i][j] *= avg0[j]
    # avg[i][3] = np.sqrt(avg[i][2] ** 2 + avg[i][3] ** 2)
    avg[i][8] = avg[i][8] + rho0 + drho0 * avg[i][5]
    avg[i][7] = (avg[i][1] + 300) * avg[i][9] * 1e-2 / 3.15e7 * 3.3e6
    avg[i][10] = (avg[i][1] + 300 * avg[i][5]) * avg[i][3] * 1e-2 / 3.15e7 * 3.3e6
    avg[i][11] = -np.gradient(avg[i][1], avg[i][0]) * 1e-6 / 1e3 * 3.3e6
    avg[i][7] = np.gradient(avg[i][7], avg[i][0]) * 1e-2
    avg[i][10] = np.gradient(avg[i][10], avg[i][0]) * 1e-2
    avg[i][11] = np.gradient(avg[i][11], avg[i][0]) * 1e-2

''' plot data '''
plt.figure(figsize=(4, 6))
# plt.fill_between([-0.002, 0.022], [1697, 1697], [1706, 1706], color='gray', alpha=0.5)
for i in range(len(cases)):
    e = avg[i]
    plt.plot(avg[i][para], avg[i][0], label=cases[i][2], linestyle='-', linewidth=3, c=cases[i][3])
    # plt.plot(avg[i][10], avg[i][0], label=cases[i][2], linestyle='-',  c=cases[i][3], linewidth=3)
    # plt.plot(avg[i][7], avg[i][0], linestyle='--', c=cases[i][3], linewidth=3)
    # plt.plot(avg[i][11], avg[i][0], linestyle=':', c=cases[i][3], linewidth=3)
    # plt.plot(avg[i][11], avg[i][0], linestyle=':', linewidth=3)
    # plt.plot(avg[i][10], avg[i][0], label='convection', linestyle='-', linewidth=3)
    # plt.plot(avg[i][7], avg[i][0], label='segregation', linestyle='--', linewidth=3)
    # plt.plot(avg[i][11], avg[i][0], label='conduction', linestyle=':', linewidth=3)
# plt.plot([0.5, 0.5], [340, 1740], c='grey', linestyle='', linewidth=6)

# plt.plot([0.85, 0.85], [1640, 1740], linestyle='--', linewidth=5, c='k')
# plt.grid()
plt.legend()
plt.xlabel(labels[para], fontsize=12, family='Arial')
plt.ylabel(labels[0], fontsize=12, family='Arial')
plt.ylim((340, 1740))
# plt.xlim(1600, 2100)
# plt.xlim(0, 1)
plt.grid()
if lim:
    plt.xlim(lim)
# if para == 4 or 10:
#     plt.xscale('log')
plt.tight_layout()
plt.show()
# plt.savefig('../../figure/fig3_2.pdf', format='pdf', transparent=True)
