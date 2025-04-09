from math import sqrt, atan, pi, acos

# ----------------------
v_nnode_p = [33, 33, 33]  # number of nodes for (x, y, z) per process
ncap = 12  # number of caps
nlevel = 1  # number of levels
nnode_p = v_nnode_p[0] * v_nnode_p[1] * \
          v_nnode_p[2]  # number of nodes per process
nnode_l = v_nnode_p[0] * v_nnode_p[1] * ncap  # number of nodes per layer
ldataIndex = {
    'temperature': [6, 6 + nnode_p],
    'residual temperature': [8 + nnode_p, 8 + 2 * nnode_p],
    'melt fraction': [10 + 2 * nnode_p, 10 + 3 * nnode_p],
    'phase density constrast': [12 + 3 * nnode_p, 12 + 4 * nnode_p],
    'bulk density constrast': [14 + 4 * nnode_p, 14 + 5 * nnode_p],
    'latent heat': [16 + 5 * nnode_p, 16 + 6 * nnode_p],
    'dike extraction': [18 + 6 * nnode_p, 18 + 7 * nnode_p],
    'olivine composition': [20 + 7 * nnode_p, 20 + 8 * nnode_p],
    'anorthosite composition': [22 + 8 * nnode_p, 22 + 9 * nnode_p],
    'velocity': [24 + 9 * nnode_p, 24 + 10 * nnode_p],
    'viscosity': [26 + 10 * nnode_p, 26 + 11 * nnode_p],
    'coord': [34 + 11 * nnode_p, 34 + 12 * nnode_p]
}  # Layered dataset: [a + layerNum, b + layerNum, v_nnode_p(z)] (from bottom to top)
dataUnit = {
    'temperature': 2104,  # K
    'residual temperature': 2104,
    'velocity': 1.81e-3,  # cm/yr
    'viscosity': 5e20,
    'radius': 1740,  # km
    'time': 9.59e4,  # Myr
    'anorthosite composition': 1,
    'olivine composition': 1,
    'melt fraction': 1,
    'dike extraction': 1
}


# ----------------------


def ReadCoordinates(files):
    """
    Read coordinates from files
    Transform xyz to rt(-90,90)f(0,360)

    Parameters
    ----------
    files : [list(str)]
        [input files]

    Returns
    -------
    [coord_phi, coord_theta, coord_r]
     list(num)   list(num)  list(num)
        notes: the unit of coord_r is km
    """
    coord = [[], [], []]
    _id = ldataIndex['coord']
    step = v_nnode_p[2]

    for _file in files:
        fp = open(_file, 'r')
        rawdata = fp.readlines()

        for i in range(_id[0], _id[1], step):
            temp = rawdata[i].split()
            x = float(temp[0])
            y = float(temp[1])
            z = float(temp[2])

            r = sqrt(x ** 2 + y ** 2 + z ** 2)

            theta = 90 - acos(z / r) / pi * 180
            phi = atan(y / x) / pi * 180
            if x < 0:
                phi += 180
            elif y < 0:
                phi += 360

            coord[0].append(phi - 180)
            coord[1].append(theta)
            coord[2].append(r * dataUnit['radius'])

    return coord


def ReadTime(_file, cycles):
    """
    Read time(Myr) from files among specified cycles

    Parameters
    ----------
    _file : [str]
        [input files]
    cycles : [list(num)]
        [specified cycles]

    Returns
    -------
    time [list(num)]
    """
    time = []
    fp = open(_file, 'r')
    rawdata = fp.readlines()
    for i in cycles:
        temp = rawdata[i].split()[1]
        time.append(float(temp) * dataUnit['time'])

    return time


def ReadLayerData(files, layers, _type, dimension=1):
    """
    Read average layered data among specified layers from files
    note: can only deal with the files in the same level!

    Parameters
    ----------
    files : [list(str)]
        [input files]
    layers : [list(num)]
        [specified layers]
    _type : [str]
        [data type]
    dimension : [num]
        [default=1]

    Returns
    -------
    layered data [list(list(num))]
    """
    ldata = [[0 for i in range(nnode_l)] for j in range(dimension)]
    nlayer = len(layers)
    _id = ldataIndex[_type]
    unit = dataUnit[_type]
    step = v_nnode_p[2]
    node = 0

    for _file in files:
        fp = open(_file, 'r')
        rawdata = fp.readlines()

        node0 = node
        for layer in layers:
            node = node0
            for i in range(_id[0] + layer, _id[1] + layer, step):
                temp = rawdata[i].split()
                for d in range(dimension):
                    temp2 = float(temp[d]) * unit
                    ldata[d][node] += temp2 / nlayer
                node += 1

    return ldata


def FindCrust(files, c0=0.8):
    ldata = [0 for i in range(nnode_l)]
    _id = ldataIndex['anorthosite composition']
    step = v_nnode_p[2]
    node = 0

    for _file in files:
        fp = open(_file, 'r')
        rawdata = fp.readlines()

        node0 = node
        for i in range(_id[0], _id[1], step):
            for layer in range(32, 0, -1):
                comp = float(rawdata[i + layer])
                if comp <= c0:
                    if layer == 32:
                        ldata[node] = 0
                    else:
                        z1 = 32 - layer
                        comp0 = float(rawdata[i + layer + 1])
                        ldata[node] = (z1 - (c0 - comp) / (comp0 - comp)) * 10
                    break
            node += 1

    return ldata

def FindLid(files, phi0=0):
    ldata = [0 for i in range(nnode_l)]
    _id = ldataIndex['melt fraction']
    step = v_nnode_p[2]
    node = 0

    for _file in files:
        fp = open(_file, 'r')
        rawdata = fp.readlines()

        node0 = node
        for i in range(_id[0], _id[1], step):
            for layer in range(32, 0, -1):
                phi = float(rawdata[i + layer])
                if phi > phi0 or layer == 1:
                    if layer >= 22:
                        ldata[node] = (32 - layer) * 10
                    else:
                        ldata[node] = 100 + (22 - layer) * 24.0729
                    break
            node += 1

    return ldata

def FindCrust2(files, dc0=0.1):
    ldata = [0 for i in range(nnode_l)]
    _id = ldataIndex['anorthosite composition']
    step = v_nnode_p[2]
    node = 0

    for _file in files:
        fp = open(_file, 'r')
        rawdata = fp.readlines()

        node0 = node
        for i in range(_id[0], _id[1], step):
            comp0 = float(rawdata[i + 32])
            for layer in range(31, 0, -1):
                comp = float(rawdata[i + layer])
                dc = abs(comp - comp0)
                if dc > dc0:
                    z1 = 32 - layer
                    comp0 = float(rawdata[i + layer + 1])
                    ldata[node] = (z1 - dc0 / dc) * 10
                    break
            node += 1

    return ldata