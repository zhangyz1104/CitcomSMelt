import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from copy import copy

''' set parameters '''
num = 8
cycles = [i for i in range(0, 60000, 200)]
timescale = 100

''' dims 
[0: r, 1: T, 2: Vh, 3: Vv, 4: eta, 5: phi, 6: c0, 7: cl, 8: rho, 9: q, 10: E, 11: M]
'''
cases = [0, 'case1_ref', 'case2_lowVis_lowAcc', 'case3_highVis', 'case4_lowSegr_lowAcc',
         'case5_highSegr', 'case6_lowExtr', 'case7_highExtr', 'case8_shallowExtr_lowAcc',
         'case9_deepExtr', 'case10_HomoRad_2', 'case11_inhomoRad', 'case12_noRad', 'case13_constDens']
[r0, T0, v0, eta0, t0, drho0, rho0] = [1740, 2104, 1.81e-3, 5e20, 9.59e4, 500, 3340]
avg0 = [r0, T0, v0, v0, eta0, 1, 1, 1, -drho0, v0, 1 / t0, 1 / t0]
ccr = 0.85
phicr = 1e-5

filename = 'D:/MUSH_result/' + cases[num] + '.time'
time0 = np.loadtxt(filename).transpose()[1] * 9.59e4

filename = 'D:/MUSH_result/' + cases[num] + '.qb.dat'
coreT = np.loadtxt(filename).transpose()[4] * 2104
time1 = np.loadtxt(filename).transpose()[0] * 9.59e4
t2 = 0
t3 = 0
for i in range(len(coreT)):
    if coreT[i] < 1700 and not t2 and time1[i] > 100:
        t2 = time1[i]
    if coreT[i] < 1600 and not t3:
        t3 = time1[i]

print(t2, t3)

T = []
phi = []
c0 = []
E = []
time = []
r = []

dike_time = 0
melt_time = 0
crust1 = 0
crust2 = 0
crust3 = 0
for cycle in cycles:
    filename0 = 'D:/MUSH_result/' + cases[num] + '.horiz_avg.0.' + str(cycle)
    filename1 = 'D:/MUSH_result/' + cases[num] + '.horiz_avg.1.' + str(cycle)
    rawdata = np.vstack((np.loadtxt(filename0), np.loadtxt(filename1)))
    avg = rawdata.transpose()
    for j in range(len(avg0)):
        avg[j] *= avg0[j]

    crust = 0
    for i in range(len(avg[6]) - 1, 0, -1):
        if len(avg[6]) - i > 1 and avg[6][i] < ccr and not crust:
            crust = (len(avg[6]) - i - 1 - (ccr - avg[6][i]) / (avg[6][i + 1] - avg[6][i])) * 10
    if np.max(avg[10]) < 1e-8 and not dike_time and time0[cycle] > 1:
        dike_time = (time[-1] + time0[cycle]) * 0.5
        # dike_time = time[-1]
    if np.max(avg[5]) < phicr and not melt_time:
        melt_time = (time[-1] + time0[cycle]) * 0.5
        # melt_time = time[-1]
        crust1 = crust

    r = avg[0]
    T.append(avg[1])
    phi.append(avg[5])
    c0.append(avg[6])
    E.append(avg[10])
    time.append(time0[cycle])

    if time[-1] > t2 and not crust2:
        crust2 = crust
    if time[-1] > t3 and not crust3:
        crust3 = crust

    if time[-1] > timescale:
        break

time0 = copy(time)
r, time = np.meshgrid(r, time)

''' plot '''
figure, (ax3) = plt.subplots(1, 1, figsize=(5, 3))

E = np.maximum(E, 1e-15)
subfig3 = ax3.pcolormesh(time, r, np.log10(E), vmin=0, vmax=-6, cmap='Reds', shading='gouraud')
ax3.contour(time, r, phi, [phicr], colors='green', linewidths=3, linestyles='dashed')
ax3.contour(time, r, c0, [ccr], colors=['blue'], linewidths=3)
# ax3.pcolormesh(time, r, c0, cmap='gray', vmin=0.85, vmax=1, shading='gouraud', alpha=0.5)
ax3.set_ylim(1600, 1740)
ax3.set_ylabel('Radius (km)', fontsize=12)
ax3.set_xlabel('Time (Myr)', fontsize=12)
# ax3.set_title('Melt fraction')
cbar3 = plt.colorbar(subfig3, ax=ax3, orientation='vertical',
                     ticks=[-6, -4, -2, 0],
                     label='Dike extraction rate (Myr$^{-1}$)')
cbar3.set_ticklabels(['$10^{-6}$', '$10^{-4}$', '$10^{-2}$', '$1$'])

plt.xlim(0, timescale)
plt.tight_layout()
plt.show()

print('dike timescale = %.2f' % dike_time)
print('melt timescale = %.2f' % melt_time)
print('max crust thickness = %.2f' % crust1)
print('min crust thickness = %.2f' % crust)
print('max2 crust thickness = %.2f' % crust2)
print('min2 crust thickness = %.2f' % crust3)

print('[%d, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]'
      % (num, dike_time, melt_time, crust, crust1, crust3, crust2))




