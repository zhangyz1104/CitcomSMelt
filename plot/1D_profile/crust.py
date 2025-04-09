import numpy as np
import matplotlib.pyplot as plt
#
# ''' set parameters '''
case = 'case1_ref_3'
cycles = [i for i in range(0, 45000, 200)]
crust = [i for i in range(28, 33)]
magma = [i for i in range(20, 28)]
colorbar = [0, 1]
form = (1, 33 * 33 * 12)
lim = None
# lim = (0, 0.2)

''' dims 
[0: r, 1: T, 2: Vh, 3: Vv, 4: eta, 5: phi, 6: c0, 7: cl, 8: rho, 9: q, 10: E, 11: M]
'''
[r0, T0, v0, eta0, t0, drho0, rho0] = [1740, 2104, 1.81e-3, 5e20, 9.59e4, 500, 3340]
avg0 = [r0, T0, v0, v0, eta0, 1, 1, 1, -drho0, v0, 1 / t0, 1 / t0]
Tcr = 1400
ccr = 0.9
phicr = 0.001

filename = 'D:/MUSH_result/' + case + '.time'
time0 = np.loadtxt(filename).transpose()[1] * 9.59e4

layer_T = []
layer_phi = []
layer_c0 = []
time = []
for cycle in cycles:
    ''' init '''
    findT = 0
    findc = 0
    findphi = 0

    ''' read data '''
    filename = 'D:/MUSH_result/' + case + '.horiz_avg.1.' + str(cycle)
    avg = np.loadtxt(filename).transpose()
    for j in range(len(avg0)):
        avg[j] *= avg0[j]

    ''' deal with data '''
    for i in range(32, -1, -1):
        r = avg[0][i]
        T = avg[1][i]
        phi = avg[5][i]
        c = avg[6][i]

        if T > Tcr and not findT:
            r1 = avg[0][i + 1]
            T1 = avg[1][i + 1]
            d = (r0 - r) - abs((T - Tcr) / (T - T1)) * (r1 - r)
            layer_T.append(d)
            findT = 1

        if c < ccr and not findc:
            if i == 32:
                layer_c0.append(0)
            else:
                r1 = avg[0][i + 1]
                c1 = avg[6][i + 1]
                d = (r0 - r) - abs((c - ccr) / (c - c1)) * (r1 - r)
                layer_c0.append(d)
            findc = 1

        if phi > phicr and not findphi:
            if i == 32:
                layer_phi.append(0)
            else:
                r1 = avg[0][i + 1]
                phi1 = avg[5][i + 1]
                d = (r0 - r) - abs((phi - phicr) / (phi - phi1)) * (r1 - r)
                layer_phi.append(d)
            findphi = 1

        if findc and findT and findphi:
            time.append(time0[cycle])
            break

''' plot '''
plt.figure(figsize=(5, 2))
plt.plot(time, layer_c0, label='composition layer', linestyle='-', linewidth=3, c='b')
# plt.plot(time, layer_T, label='temperature layer', linestyle='-', linewidth=3)
plt.plot(time, layer_phi, label='melt layer', linestyle='-', linewidth=3, c='g')
plt.xlabel('Time (Myr)', fontsize=10)
plt.ylabel('Thickness (km)', fontsize=10)
plt.xlim(0, 53)
plt.legend()
plt.grid()
plt.tight_layout()
plt.show()




