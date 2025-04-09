import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from mpl_toolkits.basemap import Basemap

from mush_visual.global_data.extract import ReadCoordinates, ReadLayerData, FindCrust, FindLid, ncap

# ----------------------
case = 'case1_ref'
# cycles = [i for i in range(3200, 3300, 100)]  # 0.1 Myr
# cycles = [i for i in range(11400, 11500, 100)]  # 1 Myr
# cycles = [i for i in range(14600, 14700, 100)]  # 2 Myr
# cycles = [i for i in range(24200, 24300, 100)]  # 30 Myr
cycles = [i for i in range(24400, 24500, 100)]  # 100 Myr
# crust = [i for i in range(28, 33)]  # 50 km
magma = [i for i in range(20, 33)]  # 500 km
colorbar = [0, 1]


# ----------------------


def plot():
    files = []
    for c in range(ncap):
        _file = 'D:/MUSH_result/' + case + '.proc' + \
                str(2 * c + 1) + '.0.vts'
        files.append(_file)

    coord = ReadCoordinates(files)
    coord = mtri.Triangulation(coord[0], coord[1])
    lon = np.linspace(-180, 180, 360)
    lat = np.linspace(-90, 90, 180)
    lon, lat = np.meshgrid(lon, lat)

    for cycle in cycles:
        files = []
        for c in range(ncap):
            _file = 'D:/MUSH_result/' + case + '.proc' + \
                    str(2 * c + 1) + '.' + str(cycle) + '.vts'
            files.append(_file)

        lid = FindLid(files)
        # temp = ReadLayerData(files, magma, 'melt fraction')[0]
        temp = ReadLayerData(files, magma, 'dike extraction')[0]
        thickness = FindCrust(files, c0=0.85)
        # thickness = FindCrust2(files, dc0=0.1)
        tli = mtri.LinearTriInterpolator(coord, lid)
        lid = tli(lon, lat)
        tli = mtri.LinearTriInterpolator(coord, temp)
        temp = tli(lon, lat)
        tli = mtri.LinearTriInterpolator(coord, thickness)
        thickness = tli(lon, lat)

        m = Basemap(projection='moll', lon_0=0, resolution='c')
        m.drawparallels(np.arange(-90, 120, 30))
        m.drawmeridians(np.arange(0, 420, 45))
        # maps = m.pcolormesh(lon, lat, temp, cmap='magma', latlon=True, vmin=0, vmax=0.2)
        # maps = m.pcolormesh(lon, lat, lid, cmap='jet', latlon=True, vmin=30, vmax=100)
        # maps = m.pcolormesh(lon, lat, thickness, cmap='jet', latlon=True, vmin=0, vmax=10)
        # maps = m.pcolormesh(lon, lat, comp, cmap='jet', latlon=True, vmin=0.6, vmax=1)
        maps = m.pcolormesh(lon, lat, temp, cmap='jet', latlon=True, vmin=0, vmax=1e-5)
        # comp = m.contourf(lon, lat, comp, [0.93, 1], latlon=True, colors='gray', alpha=0.5)
        # temp = m.contour(lon, lat, temp, [0.1], latlon=True, colors='white')
        cbar = m.colorbar(maps, location='bottom', pad="5%")
        # cbar.set_label('Crust thickness (km)')
        # cbar.set_label('Melt fraction')
        # cbar.set_label('Anorthosite composition')

        plt.show()
        # filename = 'D:\\mush_image\\python\\' + case + '\\' + case + '.thickness.' + str(cycle) + '.png'
        # plt.savefig(filename, dpi=1000)
        # plt.clf()

        print(cycle)


plot()
