import os
import imageio
import glob
from pv_control import cycles, gif, case, pvtype, output
from PIL import Image, ImageDraw, ImageFont


def createGif(duration=450):
    file = 'D:/mush_image/paraview/' + case + '/*' + pvtype + '*' + output + '.'
    timefile = 'D:/MUSH_result/' + case + '.time'
    image_list = glob.glob(file + '?.png')
    image_list += glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    # image_list = image_list[int(20400 / 200): int(64000 / 200)]
    print(image_list)
    gif_name = 'D:/mush_image/paraview/gif/' + case + \
               '.' + pvtype + '.' + output + '(2).gif'

    dt = []
    t0 = 0
    totalt = 0
    fp = open(timefile)
    rawdata = fp.readlines()
    # time = [0] + [i for i in range(1000, 22200, 200)]
    time = [i for i in range(20400, 64000, 200)]
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

def createMP4(duration=450):
    file = 'D:/mush_image/paraview/' + case + '/*' + pvtype + '*' + output + '.'
    timefile = 'D:/MUSH_result/' + case + '.time'
    image_list = glob.glob(file + '?.png')
    image_list += glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    print(image_list)
    gif_name = 'D:/mush_image/paraview/gif/' + case + \
               '.' + pvtype + '.' + output + '(2).mp4'

    dt = []
    t0 = 0
    totalt = 0
    fp = open(timefile)
    rawdata = fp.readlines()
    # time = [0] + [i for i in range(1000, 22200, 200)]
    time = [i for i in range(20400, 64000, 200)]
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
    imageio.mimsave(gif_name, frames, 'MP4', duration=dt)


def addTime():
    file = 'D:/mush_image/paraview/' + case + '/*' + pvtype + '*' + output + '.'
    timefile = 'D:/MUSH_result/' + case + '.time'
    image_list = glob.glob(file + '?.png')
    image_list += glob.glob(file + '???.png')
    image_list += glob.glob(file + '????.png')
    image_list += glob.glob(file + '?????.png')
    selected_range = (2, 53)
    print(image_list)

    fp = open(timefile)
    rawdata = fp.readlines()

    for i in range(selected_range[0], selected_range[1]):
        t1 = float(rawdata[i * 200].split()[1]) * 9.59e7
        img0 = Image.open(image_list[i])
        width = 700
        high = 90
        string = "%d kyr" % t1
        draw = ImageDraw.Draw(img0)
        fontStyle = ImageFont.truetype('C:/Windows/Fonts/Helvetica', size=50)
        draw.text((width, high), string, fill='black', font=fontStyle)
        img0.save(image_list[i], quality=100)


for i in range(cycles[0], cycles[1], 200):
    print(i)
    os.system("python pv_lib.py " + str(i))
# if gif:
#     createGif()

# createGif()
# addTime()