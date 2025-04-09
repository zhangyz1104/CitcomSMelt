import os
import imageio
import glob
from PIL import Image, ImageDraw, ImageFont

case = "case1_ref"
output = "melt"


def createGif_300(duration=30):
    file = 'D:/mush_image/python/' + case + '(3)/*' + output + '.'
    timefile = 'D:/MUSH_result/' + case + '.time'
    # image_number = [0, 10400, 14600, 18600, 21000, 22000, 22600] + [i for i in range(23000, 25200, 200)]
    image_list = glob.glob(file + '?.png')
    image_list += glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    image_list = image_list[0:1] + image_list[int(22200 / 200): int(26000 / 200)]
    # print(image_list)
    gif_name = 'D:/mush_image/paraview/gif/' + case + '.' + output + '(300Myr).gif'

    dt = []
    t0 = 0
    totalt = 0
    fp = open(timefile)
    rawdata = fp.readlines()
    time = [0] + [i for i in range(22200, 26000, 200)]
    # time = [i for i in range(20400, 64000, 200)]
    for i in time:
        t1 = float(rawdata[i].split()[1]) * 9.59e4 / duration
        _dt = t1 - t0
        dt.append(_dt)
        totalt += t1 - t0
        t0 = t1

    print(totalt)

    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))
    imageio.mimsave(gif_name, frames, 'GIF', duration=dt)

def createGif_100(duration=10):
    file = 'D:/mush_image/python/' + case + '(3)/*' + output + '.'
    timefile = 'D:/MUSH_result/' + case + '.time'
    # image_number = [0, 10400, 14600, 18600, 21000, 22000, 22600] + [i for i in range(23000, 25200, 200)]
    image_list = glob.glob(file + '?.png')
    image_list += glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    image_list = image_list[0:1] + image_list[int(22200 / 200): int(25200 / 200)]
    # print(image_list)
    gif_name = 'D:/mush_image/paraview/gif/' + case + '.' + output + '(100Myr).gif'

    dt = []
    t0 = 0
    totalt = 0
    fp = open(timefile)
    rawdata = fp.readlines()
    time = [0] + [i for i in range(22200, 25200, 200)]
    # time = [i for i in range(20400, 64000, 200)]
    for i in time:
        t1 = float(rawdata[i].split()[1]) * 9.59e4 / duration
        _dt = t1 - t0
        dt.append(_dt)
        totalt += t1 - t0
        t0 = t1

    print(totalt)

    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))
    imageio.mimsave(gif_name, frames, 'GIF', duration=dt)

def createGif_5(duration=0.5):
    file = 'D:/mush_image/python/' + case + '(3)/*' + output + '.'
    timefile = 'D:/MUSH_result/' + case + '.time'
    # image_number = [0, 10400, 14600, 18600, 21000, 22000, 22600] + [i for i in range(23000, 25200, 200)]
    image_list = glob.glob(file + '?.png')
    image_list += glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    image_list = image_list[0:1] + image_list[int(10400 / 200): int(22200 / 200)]
    # print(image_list)
    gif_name = 'D:/mush_image/paraview/gif/' + case + '.' + output + '.gif'

    dt = []
    t0 = 0
    totalt = 0
    fp = open(timefile)
    rawdata = fp.readlines()
    time = [0] + [i for i in range(10400, 22200, 200)]
    # time = [i for i in range(20400, 64000, 200)]
    for i in time:
        t1 = float(rawdata[i].split()[1]) * 9.59e4 / duration
        _dt = t1 - t0
        dt.append(_dt)
        totalt += t1 - t0
        t0 = t1

    print(totalt)

    frames = []
    for image_name in image_list:
        frames.append(imageio.imread(image_name))
    imageio.mimsave(gif_name, frames, 'GIF', duration=dt)



def addTime():
    file = 'D:/mush_image/python/' + case + '(3)/*' + output + '.'
    timefile = 'D:/MUSH_result/' + case + '.time'
    image_list = glob.glob(file + '?.png')
    image_list += glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    selected_range = (0, 130)
    # print(image_list)

    fp = open(timefile)
    rawdata = fp.readlines()

    for i in range(selected_range[0], selected_range[1]):
        t1 = float(rawdata[i * 200].split()[1]) * 9.59e4
        img0 = Image.open(image_list[i])
        width = 50
        high = 70
        string = "%.2f Myr" % t1
        draw = ImageDraw.Draw(img0)
        fontStyle = ImageFont.truetype('C:/Windows/Fonts/Helvetica', size=25)
        draw.text((width, high), string, fill='black', font=fontStyle)
        img0.save(image_list[i], quality=100)


createGif_5()
# addTime()