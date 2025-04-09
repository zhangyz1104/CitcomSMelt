import sys
from PIL import Image
import pv_control as pv
import matplotlib.pyplot as plt
import sys
sys.path.append('F:/app/paraview/ParaView 5.8.0-Windows-Python3.7-msvc2015-64bit/bin/Lib/site-packages')
from paraview.simple import *

cycle = sys.argv[1]
procNum = range(24)
procList = {}
blockdata = 0
renderView = 0


def init():
    global renderView
    LoadPalette(paletteName='WhiteBackground')
    renderView = GetActiveViewOrCreate('RenderView')
    renderView.ViewSize = [1000, 1100]

    return renderView


def inputVtsFile():
    global blockdata
    filename = 'D:\\MUSH_result\\' + pv.case + '.' + cycle + '.vtm'
    blockdata = XMLMultiBlockDataReader(FileName=[filename])


def setSlice():
    slice1 = Slice(Input=blockdata)
    slice1.SliceType.Origin = [0.0, 1e-5, 0.0]
    slice1.SliceType.Normal = pv.normal
    slice1Display = Show(slice1, renderView, 'GeometryRepresentation')
    ColorBy(slice1Display, ('POINTS', pv.output, 'X'))
    LUT = GetColorTransferFunction(pv.output)
    LUT.RescaleTransferFunction(pv.pvrange[0], pv.pvrange[1])
    if pv.output == 'melt_fraction':
        LUT.ApplyPreset('Magma (matplotlib)', True)
    elif pv.output == 'composition_2':
        LUT.ApplyPreset('Grayscale', True)
    else:
        LUT.ApplyPreset('Cool to Warm (Extended)', True)

    # colorbar setting
    colorbar = GetScalarBar(LUT, renderView)
    colorbar.AutoOrient = 0
    colorbar.Orientation = 'Horizontal'
    colorbar.WindowLocation = 'LowerCenter'
    # colorbar.Position = [0.3564940962761126, 0.1]
    if pv.output == 'melt_fraction':
        colorbar.Title = 'melt fraction'
    elif pv.output == 'composition_2':
        colorbar.Title = 'anorthosite composition'
    else:
        colorbar.Title = pv.output
    colorbar.ComponentTitle = ''
    colorbar.TitleFontFamily = 'Times'
    colorbar.TitleBold = 1
    colorbar.TitleFontSize = 21
    colorbar.LabelFontFamily = 'Times'
    colorbar.LabelBold = 1
    colorbar.LabelFontSize = 20
    colorbar.AutomaticLabelFormat = 1
    colorbar.UseCustomLabels = 1
    # colorbar.CustomLabels = [0, 0.05, 0.1, 0.15, 0.2]
    colorbar.CustomLabels = [0.7, 0.8, 0.9, 1]
    # colorbar.CustomLabels = [0, 0.1, 0.2, 0.3, 0.4]
    colorbar.AddRangeLabels = 0
    colorbar.DrawTickMarks = 1
    # colorbar.LabelFormat = '%-#6.1g'
    colorbar.TextPosition = 'Ticks left/bottom, annotations right/top'
    colorbar.ScalarBarThickness = 20
    colorbar.ScalarBarLength = 0.4

    # sphere
    sphere1 = Sphere()
    sphere1.Radius = 0.22
    sphere1.ThetaResolution = 20
    sphere1.PhiResolution = 20
    sphere1Display = Show(sphere1, renderView, 'GeometryRepresentation')

    return slice1Display


def setContour():
    for i in procNum:
        contour1 = Contour(Input=procList['proc' + str(i)])
        contour1.ContourBy = ['POINTS', pv.output]
        contour1.Isosurfaces = pv.pvrange
        contour1Display = Show(contour1, renderView, 'GeometryRepresentation')
        LUT = GetColorTransferFunction(pv.output)
        LUT.RescaleTransferFunction(pv.pvrange[0], pv.pvrange[1])


    sphere1 = Sphere()
    sphere1.Radius = 0.22
    sphere1.ThetaResolution = 20
    sphere1.PhiResolution = 20
    sphere1Display = Show(sphere1, renderView, 'GeometryRepresentation')
    sphere1 = Sphere()
    sphere1.Radius = 1
    sphere1.ThetaResolution = 100
    sphere1.PhiResolution = 100
    sphere1Display = Show(sphere1, renderView, 'GeometryRepresentation')
    sphere1Display.Opacity = 0.1
    # sphere1Display.DiffuseColor = [0.0, 0.0, 0.0]


def setGlyph():
    procNum1 = [i*2+1 for i in range(48)]
    for i in procNum1:
        r = Calculator(Input=procList['proc' + str(i)])
        r.ResultArrayName = 'r'
        r.Function = 'sqrt(coordsX^2+coordsY^2+coordsZ^2)'
        phi = Calculator(Input=r)
        phi.ResultArrayName = 'phi'
        phi.Function = 'atan(coordsY/coordsX) - 0.5*(sign(coordsX)-1)*sign(coordsY)*3.1415926'
        theta = Calculator(Input=phi)
        theta.ResultArrayName = 'theta'
        theta.Function = 'acos(coordsZ/r)'
        Vx = Calculator(Input=theta)
        Vx.ResultArrayName = 'Vx'
        Vx.Function = 'velocity_Z*sin(theta)*cos(phi)+velocity_X*cos(theta)*cos(phi)-velocity_Y*sin(phi)'
        Vy = Calculator(Input=Vx)
        Vy.ResultArrayName = 'Vy'
        Vy.Function = 'velocity_Z*sin(theta)*sin(phi)+velocity_X*cos(theta)*sin(phi)+velocity_Y*cos(phi)'
        Vz = Calculator(Input=Vy)
        Vz.ResultArrayName = 'Vz'
        Vz.Function = 'velocity_Z*sin(theta)-velocity_X*sin(phi)'
        Vxyz = Calculator(Input=Vz)
        Vxyz.ResultArrayName = 'Vxyz'
        Vxyz.Function = 'Vx*iHat+0*jHat+Vz*kHat'
        logV = Calculator(Input=Vz)
        logV.ResultArrayName = 'logV'
        logV.Function = 'log10(sqrt(Vx^2+Vy^2+Vz^2))'
        slice1 = Slice(Input=Vxyz)
        slice1.SliceType.Origin = [0.0, -0.02, 0.0]
        slice1.SliceType.Normal = pv.normal
        logV = Calculator(Input=slice1)
        logV.ResultArrayName = 'logV'
        logV.Function = 'log10(velocity)'
        glyph1 = Glyph(Input=logV)
        glyph1.OrientationArray = ['POINTS', 'Vxyz']
        glyph1.GlyphType = 'Arrow'
        glyph1.ScaleArray = ['POINTS', 'velocity']
        # glyph1.ScaleArray = ['POINTS', 'logV']
        glyph1.ScaleArray = ['POINTS', 'No scale array']
        glyph1.ScaleFactor = pv.arrowScale
        # glyph1.GlyphMode = 'Every Nth Point'
        glyph1.MaximumNumberOfSamplePoints = 5
        glyph1Display = Show(glyph1, renderView, 'GeometryRepresentation')
        ColorBy(glyph1Display, ('POINTS', 'temperature'))
        temperatureLUT = GetColorTransferFunction('temperature')
        temperatureLUT.RGBPoints = [0, 0, 0, 0.1, 0, 0]
        glyph1Display.AmbientColor = [0.0, 0.0, 0.0]
        glyph1Display.DiffuseColor = [0.0, 0.0, 0.0]


def saveImage():
    global renderView
    renderView.OrientationAxesVisibility = 0
    renderView.CameraViewUp = [0, 1, 0]
    renderView.CameraPosition =[0, 500, 0]
    renderView.ResetCamera()
    if pv.save:
        FileName = 'D:\\mush_image\\paraview\\' + pv.case + '\\c' + pv.case + \
            '.' + pv.pvtype + '.' + pv.output + '.' + cycle + '.png'
    else:
        FileName = 'D:\\temp\\temp.png'
    SaveScreenshot(FileName, renderView, ImageResolution=[1000, 1100])
    # SaveScreenshot(FileName, renderView, ImageResolution=[1000, 1100], TransparentBackground=1)
    if pv.show:
        plt.imshow(Image.open(FileName))
        plt.title(pv.case + '.' + cycle)
        plt.axis('off')
        plt.show()


def main():
    init()
    inputVtsFile()
    if pv.pvtype == 'slice':
        setSlice().SetScalarBarVisibility(renderView, True)
        # setContour()
        if pv.glyph:
            setGlyph()
    elif pv.pvtype == 'contour':
        setContour()
    saveImage()


main()
