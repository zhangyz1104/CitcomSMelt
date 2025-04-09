import matplotlib.pyplot as plt
import numpy as np
import os
import imageio
import glob

''' dims
[0: r, 1: T, 2: Vh, 3: Vv, 4: eta, 5: phi, 6: c0, 7: cl, 8: rho, 9: q, 10: E, 11: M]
'''
[r0, T0, v0, eta0, t0, drho0, rho0] = [1740, 2104, 1.81e-3, 1e21, 9.59e4, 650, 3250]
avg0 = [r0, T0, v0, v0, eta0, 1, 1, 1, 3.6e-3, v0, 1 / t0, 3.6e-3]
labels = ['Radius (km)', 'Temperature (K)', 'Horizental velocity (cm/yr)', 'Velocity (cm/yr)',
          'Viscosity (Pa$\cdot$s)', 'Melt fraction', 'Anorthite bulk composition',
          'Anorthosite melt composition', 'Cooling rate (K/kyr)', 'Darcy flux (cm/yr)',
          'Dike extraction rate (Myr$^{-1}$)', 'Cooling rate (K/kyr)']

filename = 'D:/MUSH_result/debug_ref.time'
time0 = np.loadtxt(filename).transpose()[1] * 9.59e4

# for i in range(5000, 20000, 200):
#     ''' read data '''
#     filename0 = 'D:/MUSH_result/debug_ref.horiz_avg.0.' + str(i)
#     filename1 = 'D:/MUSH_result/debug_ref.horiz_avg.1.' + str(i)
#     rawdata = np.vstack((np.loadtxt(filename0), np.loadtxt(filename1)[1:]))
#     avg = rawdata.transpose()
#     for j in range(len(avg0)):
#         avg[j] *= avg0[j]
#
#     ''' plot data '''
#     plt.figure(figsize=(4, 6))
#     plt.plot(avg[1], avg[0], label='%.1eMyr' % time0[i], linestyle='-', linewidth=3, c='k')
#     plt.grid()
#     plt.legend()
#     plt.xlabel(labels[1], fontsize=12, family='Arial')
#     plt.ylabel(labels[0], fontsize=12, family='Arial')
#     plt.ylim((340, 1740))
#     plt.xlim(1600, 2100)
#     plt.tight_layout()
#     savename = '../../figure/debug_ref/T' + str(i) + '.png'
#     plt.savefig(savename)
#
#     print(i)


def createGif(duration=450):
    file = '../../figure/debug_ref/T'
    image_list = glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    gif_name = '../../figure/debug_ref/T.gif'

    dt = []
    t0 = 0
    totalt = 0
    time = [i for i in range(0, 20000, 200)]
    for i in time:
        t1 = time0[i] / duration
        _dt = t1 - t0
        dt.append(_dt)
        totalt += t1 - t0
        t0 = t1

    print(totalt)

    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))
    imageio.mimsave(gif_name, frames, 'GIF', duration=dt)

createGif()