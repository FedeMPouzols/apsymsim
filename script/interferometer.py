#
#
#    This file is part of APSYNSIM: A real-time Aperture Synthesis Simulator
#
#    Copyright (C) 2014  Ivan Marti-Vidal (Nordic ARC Node, OSO, Sweden)
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

from __future__ import (absolute_import, print_function)

import numpy as np
import pylab as pl
import os
import scipy.ndimage.interpolation as spndint
from ScrolledText import ScrolledText
import sys
import time

import matplotlib as mpl
from matplotlib.widgets import Slider, Button
import matplotlib.cm as cm
import matplotlib.image as plimg

try:
    import Tkinter as Tk
except ImportError:
    import tkinter as Tk

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
# although not used explicitly, required for projection='3d' in pl.axes()
from mpl_toolkits.mplot3d import Axes3D

import tkFileDialog
from tkMessageBox import showinfo

mpl.use('TkAgg')

import simple_clean_img
import cleaner

__help_text__ = """
     APSYNSIM, A REAL-TIME APERTURE SYNTHESIS SIMULATOR

                     IVAN MARTI-VIDAL
(ONSALA SPACE OBSERVATORY, NORDIC ALMA REGIONAL CENTER NODE)

You can click and drag the antennas in the plot called "ARRAY CONFIGURATION".
When you drag an antenna, all other plots (UV PLANE, DIRTY BEAM, and DIRTY
IMAGE) will be updated automatically (may need some time to refresh,
especially if working on Windows and/or with many antennas).

You can also click on any point of the DIRTY BEAM, MODEL IMAGE, or DIRTY
IMAGE plots, and the program will tell you the intensity value and the pixel
coordinates.

If you click on the UV PLANE image, the program will print the value of the
source Fourier transform at that point. If you click close to a point observed
with the interferometer, the program will tell you the baseline and hour
angle of observation.

You can also change the observing latitude, hour-angle coverage, source
declination, and observing wavelength by clicking on the blue sliders at
the bottom-right corner of the figure. The plots will be updated
automatically (may also need some time to refresh all plots).

The dirty beam is computed using Briggs weighting. The robustness parameter
can be changed by shifting the corresponding blue slider (robustness of -2
tends to uniform weighting, whereas +2 tends to natural weighting).

You can add and/or subtract antennas by pressing the "+ Antenna" and
"- Antenna" buttons. New antennas are inserted at the array origin (0,0).
If you add, drag, and subtract an antenna, the program will remember the
last antenna positions if you add them again.

You can save the current array, load a new array (for instance, from the
EXAMPLES folder), and/or load a new source model (for instance, from
the EXAMPLES folder) by pressing the corresponding buttons "Save array",
"Load array" and "Load model".

You can also zoom in/out by pressing "Z" or "z" (respectively). The program
will then zoom using the current cursor position as zooming center.

Pressing "c" will toggle the color code of the figures (from hue to grayscale).

Pressing "u" will pop-up a window with several plots in Fourier space.

Pressing the "Reduce data" button will open a new window, where you can
CLEAN your dirty image and apply corrupting gains to your antennas (see
help in that window for more details).

Enjoy!

"""


class Interferometer(object):

    def quit(self, event=None):

        self.tks.destroy()
        sys.exit()

    def __init__(self, antenna_file="", model_file="", tkroot=None):

    #  if tkroot is None:
    #    self.tks = Tk.Tk()
    #  else:
        self.tks = tkroot
        self.tks.protocol("WM_DELETE_WINDOW", self.quit)
        self.Hfac = np.pi / 180. * 15.
        self.deg2rad = np.pi / 180.
        self.curzoom = [0, 0, 0, 0]
        self.robust = 0.0
        self.deltaAng = 1. * self.deg2rad
        self.gamma = 0.5  # Gamma correction to plot model.
        self.lfac = 1.e6  # Lambda units (i.e., 1.e6 => Mlambda)
        self.ulab = r'U (M$\lambda$)'
        self.vlab = r'V (M$\lambda$)'
        self.W2W1 = 1.0   # Relative weighting for subarrays.
        self.currcmap = cm.jet

        self.GUIres = True  # Make some parts of the GUI respond to events
        self.antLock = False  # Lock antenna-update events

        self.my_cleaner = None   # Cleaner window instance (when initialized)

# Default of defaults!
        nH = 200
        Npix = 512   # Image pixel size. Must be a power of 2
        DefaultModel = 'Nebula.model'
        DefaultArray = 'Long_Golay_12.array'

# Overwrite defaults from config file:
        d1 = os.path.dirname(os.path.realpath(__file__))
        print(d1)

#   execfile(os.path.join(os.path.basename(d1),'apsynsim.config'))
        try:
            conf = open(os.path.join(d1, 'apsynsim.config'))
        except:
            d1 = os.getcwd()
            conf = open(os.path.join(d1, 'apsynsim.config'))

        for line in conf.readlines():
            temp = line.replace(' ', '')
            if len(temp) > 2:
                if temp[0:4] == 'Npix':
                    Npix = int(temp[5:temp.find('#')])
                if temp[0:2] == 'nH':
                    nH = int(temp[3:temp.find('#')])
                if temp[0:10] == 'DefaultMod':
                    DefaultModel = temp[12:temp.find('#')].replace(
                        '\'', '').replace('\"', '')
                if temp[0:12] == 'DefaultArray':
                    DefaultArray = temp[14:temp.find('#')].replace(
                        '\'', '').replace('\"', '')

        conf.close()

# Set instance configuration values:
        self.nH = nH
        self.Npix = Npix

        self.datadir = os.path.join(d1, '..', 'pictures')
        self.arraydir = os.path.join(d1, '..', 'arrays')
        self.modeldir = os.path.join(d1, '..', 'source_models')

     # Try to read a default initial array:
        if len(antenna_file) == 0:
            try:
                antenna_file = os.path.join(self.arraydir, DefaultArray)
            except:
                pass

     # Try to read a default initial model:
        if len(model_file) == 0:
            try:
                model_file = os.path.join(self.modeldir, DefaultModel)
            except:
                pass

        self.lock = False
        self._onSphere = False

        self.readModels(str(model_file))
        self.readAntennas(str(antenna_file))
        self.init_GUI()  # makefigs=makefigs)

    def showError(self, message):
        showinfo('ERROR!', message)
        raise Exception(message)

    def _getHelp(self):
        win = Tk.Toplevel(self.tks)
        win.title("Help")
        helptext = ScrolledText(win)
        helptext.config(state=Tk.NORMAL)
        helptext.insert('1.0', __help_text__)
        helptext.config(state=Tk.DISABLED)

        helptext.pack()
        Tk.Button(win, text='OK', command=win.destroy).pack()

    def init_GUI(self):  # ,makefigs=True):

        mpl.rcParams['toolbar'] = 'None'

        self.Nphf = self.Npix / 2
        self.robfac = 0.0
        self.figUV = pl.figure(figsize=(15, 8))

        if self.tks is None:
            self.canvas = self.figUV.canvas
        else:
            self.canvas = FigureCanvasTkAgg(self.figUV, master=self.tks)
            self.canvas.show()
            menubar = Tk.Menu(self.tks)
            menubar.add_command(label="Help", command=self._getHelp)
            menubar.add_command(label="Quit", command=self.quit)

            self.tks.config(menu=menubar)
            self.canvas.get_tk_widget().pack(
                side=Tk.TOP, fill=Tk.BOTH, expand=1)

        # UVPlot not visible. Note antPlot would overlaps UVPlot
        show_uvplot = False
        self.UVPlot = self.figUV.add_subplot(232, aspect='equal',
                                             axisbg=(0.4, 0.4, 0.4))
        sphere_dim = 0.12
        if not show_uvplot:
            sphere_dim = 0
        self.spherePlot = pl.axes([0.19, 0.82, sphere_dim, sphere_dim],
                                  projection='3d', aspect='equal')
        if not show_uvplot:
            self.UVPlot.set_visible(False)
            self.spherePlot.set_visible(False)

        self.antPlot = self.figUV.add_subplot(232, aspect='equal')
        # beamPlot (this is dirty-beam) not visible (bottom-right panel).
        # This space is currently used for the clean image
        self.show_beamPlot = False
        self.beamPlot = self.figUV.add_subplot(231, aspect='equal')
        self.beamText = self.beamPlot.text(0.05, 0.80, '',
                                           transform=self.beamPlot.transAxes,
                                           bbox=dict(facecolor='white', alpha=0.7))
        if not self.show_beamPlot:
            self.beamPlot.set_visible(False)

        self.modelPlot = self.figUV.add_subplot(234, aspect='equal')
        self.dirtyPlot = self.figUV.add_subplot(235, aspect='equal')

        self.cleanPlot = self.figUV.add_subplot(236, aspect='equal')
        self.cleanPlot.set_adjustable('box-forced')

        u = np.linspace(0, 2 * np.pi, 100)
        v = np.linspace(0, np.pi, 100)
        x = 10 * np.outer(np.cos(u), np.sin(v))
        y = 10 * np.outer(np.sin(u), np.sin(v))
        z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
        self.spherePlotPlot = self.spherePlot.plot_surface(x, y, z, rstride=4,
                                                           cstride=4, color=(0.8, 0.8, 1.0))
        self.spherePlot._axis3don = False
        self.spherePlot.patch.set_alpha(0.8)
        self.arrayPath = [
            np.zeros(self.nH), np.zeros(self.nH), np.zeros(self.nH)]
        self.sphereArray = self.spherePlot.plot([], [], [], 'y', linewidth=3)
        self.spherePlot.set_xlim3d((-6, 6))
        self.spherePlot.set_ylim3d((-6, 6))
        self.spherePlot.set_zlim3d((-6, 6))
        self.spherePlot.patch.set_alpha(0.8)
        self.spherePlot.elev = 45.

        self.figUV.subplots_adjust(
            left=0.05, right=0.99, top=0.95, bottom=0.07, hspace=0.25)
        self.canvas.mpl_connect('pick_event', self._onPick)
        self.canvas.mpl_connect('motion_notify_event', self._onAntennaDrag)
        self.canvas.mpl_connect('button_release_event', self._onRelease)
        self.canvas.mpl_connect('button_press_event', self._onPress)
        self.canvas.mpl_connect('key_press_event', self._onKeyPress)
        self.pickAnt = False

        self.fmtH = r'$\phi = $ %3.1f$^\circ$   $\delta = $ %3.1f$^\circ$' "\n" r'H = %3.1fh / %3.1fh'
        self.fmtBas = r'Bas %i $-$ %i  at  H = %4.2fh'
        self.fmtVis = r'Amp: %.1e Jy.   Phase: %5.1f deg.'
        self.fmtA = 'N = %i'
        self.fmtA2 = '  Picked Ant. #%i'
        self.fmtA3 = '\n%6.1fm | %6.1fm'
        fmtB1 = r'$\lambda = $ %4.1fmm  ' % (self.wavelength[2] * 1.e6)
        self.fmtB = fmtB1 + "\n" + r'% 4.2f Jy/beam' + "\n" + \
            r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f '
        self.fmtD = r'% .2e Jy/beam' "\n" r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f '
        self.fmtM = r'%.2e Jy/pixel' "\n"  r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f'

        self.wax = {}
        self.widget = {}

        # Place bars for latitude, dec, wave, robust, etc.
        # offsets wrt (original) bottom left panel
        bars_x = 0.63
        bars_y = 0.5 - 0.08  # -0.08 to re-center vertically for now
        self.wax['lat'] = pl.axes([bars_x + 0.07, bars_y + 0.45, 0.25, 0.04])
        self.wax['dec'] = pl.axes([bars_x + 0.07, bars_y + 0.40, 0.25, 0.04])
        self.wax['H0'] = pl.axes([bars_x + 0.07, bars_y + 0.35, 0.25, 0.04])
        self.wax['H1'] = pl.axes([bars_x + 0.07, bars_y + 0.30, 0.25, 0.04])
        self.wax['wave'] = pl.axes([bars_x + 0.07, bars_y + 0.25, 0.25, 0.04])
        self.wax['robust'] = pl.axes(
            [bars_x + 0.07, bars_y + 0.20, 0.25, 0.04])

        # Place buttons
        # offsets wrt (original) bottom left panel
        # 3 x rows: +0.07, 0.155, 0.24
        but_x = 0
        but_y = 0.5
        self.wax['loadarr'] = pl.axes([but_x + 0.07, but_y + 0.38, 0.10, 0.05])
        self.wax['save'] = pl.axes([but_x + 0.07, but_y + 0.32, 0.10, 0.05])
        self.wax['loadmod'] = pl.axes([but_x + 0.07, but_y + 0.26, 0.10, 0.05])
        self.wax['add'] = pl.axes([but_x + 0.28, but_y + 0.38, 0.08, 0.05])
        self.wax['rem'] = pl.axes([but_x + 0.28, but_y + 0.32, 0.08, 0.05])

        self.wax['reduce'] = pl.axes([but_x + 0.07, but_y + 0.20, 0.10, 0.05])
        self.wax['clean'] = pl.axes([but_x + 0.12, but_y + 0.08, 0.16, 0.05])
        have_quit = False
        if have_quit:
            self.wax['quit'] = pl.axes(
                [but_x + 0.155, but_y + 0.02, 0.08, 0.05])

        # Place plot output bars/labels
        # for the model image
        # offsets wrt (original) bottom center panel
        gamma_x = -0.33
        gamma_y = 0
        self.wax['gammacorr'] = pl.axes(
            [gamma_x + 0.46, gamma_y + 0.08, 0.13, 0.02], axisbg='white')
        # for the dirty image
        # offsets wrt (original) bottom center panel
        dirty_x = -0.33
        dirty_y = 0
        self.wax['diameter'] = pl.axes(
            [dirty_x + 0.825, dirty_y + 0.08, 0.10, 0.02], axisbg='white')

        # log(W1/W2) bar when visible for the array
        arr_x = 0.33
        arr_y = 0
        self.wax['subarrwgt'] = pl.axes([arr_x + 0.15, arr_y + 0.58,
                                         0.12, 0.02], axisbg='white')

        # create widgets for bars
        self.widget['lat'] = Slider(self.wax['lat'], r'Lat (deg)',
                                    -90., 90., valinit=self.lat / self.deg2rad)
        self.widget['dec'] = Slider(self.wax['dec'], r'Dec (deg)',
                                    -90., 90., valinit=self.dec / self.deg2rad)
        self.widget['H0'] = Slider(self.wax['H0'], r'H$_{0}$ (h)',
                                   -12., 12., valinit=self.Hcov[0] / self.Hfac)
        self.widget['H1'] = Slider(self.wax['H1'], r'H$_{1}$ (h)',
                                   -12., 12., valinit=self.Hcov[1] / self.Hfac)
        self.widget['wave'] = Slider(self.wax['wave'], r'$\lambda$ (mm)',
                                     self.wavelength[0] * 1.e6,
                                     self.wavelength[1] * 1.e6,
                                     valinit=self.wavelength[2] * 1.e6)
        self.widget['robust'] = Slider(self.wax['robust'], r'Robust',
                                       -2., 2., valinit=0.0)

        # create widgets for buttons
        self.widget['add'] = Button(self.wax['add'], r'+ Antenna')
        self.widget['rem'] = Button(self.wax['rem'], r'- Antenna')
        self.widget['reduce'] = Button(self.wax['reduce'], r'Adv. reduction')
        self.widget['clean'] = Button(self.wax['clean'], r'Clean image')
        clean_label = self.widget['clean'].label
        clean_label.set_fontsize(14)
        clean_label.set_weight('bold')
        self.widget['save'] = Button(self.wax['save'], 'Save array')
        self.widget['loadarr'] = Button(self.wax['loadarr'], 'Load array')
        self.widget['loadmod'] = Button(self.wax['loadmod'], 'Load model')
        if have_quit:
            self.widget['quit'] = Button(self.wax['quit'], 'Quit')

        # create widgets for output bars/labels
        self.widget['gammacorr'] = Slider(self.wax['gammacorr'],
                                          'gamma', 0.1, 1.0, valinit=self.gamma,
                                          color='red')
        self.widget['gammacorr'].label.set_color('white')
        self.widget['gammacorr'].valtext.set_color('white')

        self.widget['diameter'] = Slider(self.wax['diameter'], 'Dish size (m)',
                                         0, 100., valinit=0.0, color='red')
        self.widget['diameter'].label.set_color('white')
        self.widget['diameter'].valtext.set_color('white')

        self.widget['subarrwgt'] = Slider(self.wax['subarrwgt'], 'log(W1/W2)',
                                          -4, 4, valinit=0, color='red')
        self.widget['subarrwgt'].label.set_fontsize(9)
        self.widget['subarrwgt'].valtext.set_fontsize(9)

        # set on_ methods for bars
        self.widget['lat'].on_changed(self._onKeyLat)
        self.widget['dec'].on_changed(self._onKeyDec)
        self.widget['H0'].on_changed(self._onKeyH0)
        self.widget['H1'].on_changed(self._onKeyH1)
        self.widget['wave'].on_changed(self._changeWavelength)
        self.widget['robust'].on_changed(self._onRobust)

        # set on_ methods for buttons
        self.widget['add'].on_clicked(self._addAntenna)
        self.widget['rem'].on_clicked(self._removeAntenna)
        self.widget['save'].on_clicked(self.saveArray)
        self.widget['loadarr'].on_clicked(self.loadArray)
        self.widget['loadmod'].on_clicked(self.loadModel)
        self.widget['gammacorr'].on_changed(self._gammacorr)
        if have_quit:
            self.widget['quit'].on_clicked(self.quit)
        self.widget['reduce'].on_clicked(self._reduce)
        self.widget['clean'].on_clicked(self._clean_img)

        # set on_ methods for output bars/labels
        self.widget['subarrwgt'].on_changed(self._subarrwgt)
        self.widget['diameter'].on_changed(self._setDiameter)

        self._prepareBeam()
        self._prepareBaselines()
        self._setBaselines()
        self._setBeam()
        self._plotBeam()
        self._plotAntennas()
        self._prepareModel()
        self._plotModel()
        self._plotDirty()
        self._plotModelFFT()

        self._init_clean_img()

        self.canvas.draw()

    def _init_clean_img(self):
        # The cleaner with functionality separated from GUI
        self.my_clean_img = simple_clean_img.SimpleCleanImg(
            self, self.cleanPlot)
        # calculate and plot clean image - just 1 iteration
        self.my_clean_img.do_clean(pre_iter=0)

    def _setDiameter(self, diam):

        self.Diameters[0] = diam
        if self.GUIres:
            self._setPrimaryBeam(replotFFT=True)
            self._changeCoordinates(rescale=True)

    def _reduce(self, event):

        if self.tks is not None:
            self.my_cleaner = cleaner.Cleaner(self)

    def _clean_img(self, event):

        if self.my_clean_img is not None:
            self.cleanPlot.set_xlim((self.curzoom[1][0],
                                     self.curzoom[1][1]))
            self.cleanPlot.set_ylim((self.curzoom[1][2],
                                     self.curzoom[1][3]))
            self.my_clean_img.do_clean()

    def readAntennas(self, antenna_file):

        self.subarray = False
        self.Hcov = [-12.0 * self.Hfac, 12.0 * self.Hfac]
        self.Hmax = np.pi
        self.lat = 45. * self.deg2rad
        self.dec = 60. * self.deg2rad
        self.trlat = [np.sin(self.lat), np.cos(self.lat)]
        self.trdec = [np.sin(self.dec), np.cos(self.dec)]
        self.Xmax = 4.0
        self.Diameters = [0., 0.]
        self.wavelength = [3.e-6, 21.e-5, 6.e-5]  # in km.

        if len(antenna_file) == 0:
            self.Nant = 7
            self.antPos = [[0.0, 0.0], [0.0, 1.], [0.0, 2.0],
                           [1., -1.], [2.0, -2.0], [-1., -1.], [-2.0, -2.0]]
            self.antPos2 = []
            self.Nant2 = 0

        if len(antenna_file) > 0:
            if not os.path.exists(antenna_file):
                self.showError(
                    "\n\nAntenna file %s does not exist!\n\n" % antenna_file)
                return False

            else:
                antPos = []
                antPos2 = []
                Hcov = [0, 0]
                Nant = 0
                Nant2 = 0
                Xmax = 0.0
                fi = open(antenna_file)
                for li, l in enumerate(fi.readlines()):
                    comm = l.find('#')
                    if comm >= 0:
                        l = l[:comm]
                    it = l.split()
                    if len(it) > 0:

                        if it[0] == 'WAVELENGTH':
                            self.wavelength = [
                                float(it[1]) * 1.e-3, float(it[2]) * 1.e-3]
                            self.wavelength.append(
                                (self.wavelength[0] + self.wavelength[1]) / 2.)
                        elif it[0] == 'ANTENNA':
                            antPos.append(map(float, it[1:]))
                            Nant += 1
                            antPos[-1][0] *= 1.e-3
                            antPos[-1][1] *= 1.e-3
                            Xmax = np.max(np.abs(antPos[-1] + [Xmax]))
                        elif it[0] == 'ANTENNA2':
                            antPos2.append(map(float, it[1:]))
                            Nant2 += 1
                            antPos2[-1][0] *= 1.e-3
                            antPos2[-1][1] *= 1.e-3
                            Xmax = np.max(np.abs(antPos2[-1] + [Xmax]))
                        elif it[0] == 'DIAMETER':
                            Diams = map(float, it[1:])
                            self.Diameters[0] = Diams[0]
                            if len(Diams) > 1:
                                self.Diameters[1] = Diams[1]
                        elif it[0] == 'LATITUDE':
                            lat = float(it[1]) * self.deg2rad
                            trlat = [np.sin(lat), np.cos(lat)]
                        elif it[0] == 'DECLINATION':
                            dec = float(it[1]) * self.deg2rad
                            trdec = [np.sin(dec), np.cos(dec)]
                        elif it[0] == 'HOUR_ANGLE':
                            Hcov[0] = float(it[1]) * self.Hfac
                            Hcov[1] = float(it[2]) * self.Hfac
                        else:
                            self.showError("\n\nWRONG SYNTAX IN LINE %i:\n\n %s...\n\n" %
                                           (li + 1, l[:max(10, len(l))]))

                if Nant2 > 1:
                    self.subarray = True

                if np.abs(lat - dec >= np.pi / 2.):
                    self.showError(
                        "\n\nSource is either not observable or just at the horizon!\n\n")
                    return False
                if Nant < 2:
                    self.showError(
                        "\n\nThere should be at least 2 antennas!\n\n")
                    return False

                self.Nant = Nant
                self.antPos = antPos
                self.Nant2 = Nant2
                self.antPos2 = antPos2
                self.lat = lat
                self.dec = dec
                self.trlat = trlat
                self.trdec = trdec
                self.Hcov = Hcov
                self.Xmax = Xmax

                cosW = -np.tan(self.lat) * np.tan(self.dec)
                if np.abs(cosW) < 1.0:
                    Hhor = np.arccos(cosW)
                elif np.abs(self.lat - self.dec) > np.pi / 2.:
                    Hhor = 0
                else:
                    Hhor = np.pi

                if Hhor > 0.0:
                    if self.Hcov[0] < -Hhor:
                        self.Hcov[0] = -Hhor
                    if self.Hcov[1] > Hhor:
                        self.Hcov[1] = Hhor

                self.Hmax = Hhor
                self.Xmax = Xmax * 1.5
                fi.close()

        return True

    def readModels(self, model_file):

        self.imsize = 4.
        self.imfiles = []

        if len(model_file) == 0:
            self.models = [['G', 0., 0.4, 1.0, 0.1], ['D', 0., 0., 2., 0.5],
                           ['P', -0.4, -0.5, 0.1]]
            self.Xaxmax = self.imsize / 2.
            return True

        if len(model_file) > 0:
            if not os.path.exists(model_file):
                self.showError(
                    "\n\nModel file %s does not exist!\n\n" % model_file)
                return False

            else:
                fixsize = False
                models = []
                imfiles = []
                Xmax = 0.0
                fi = open(model_file)
                for li, l in enumerate(fi.readlines()):
                    comm = l.find('#')
                    if comm >= 0:
                        l = l[:comm]
                    it = l.split()
                    if len(it) > 0:
                        if it[0] == 'IMAGE':
                            imfiles.append([str(it[1]), float(it[2])])
                        elif it[0] in ['G', 'D', 'P']:
                            models.append([it[0]] + map(float, it[1:]))
                            if models[-1][0] != 'P':
                                models[-1][4] = np.abs(models[-1][4])
                                Xmax = np.max(
                                    [np.abs(models[-1][1]) + models[-1][4],
                                     np.abs(models[-1][2]) + models[-1][4], Xmax])
#          elif it[0] == 'WAVELENGTH':
#            wavelength = float(it[1])*1.e-3
                        elif it[0] == 'IMSIZE':
                            imsize = 2. * float(it[1])
                            fixsize = True
                        else:
                            self.showError("\n\nWRONG SYNTAX IN LINE %i:\n\n %s...\n\n" %
                                           (li + 1, l[:max(10, len(l))]))

                if len(models) + len(imfiles) == 0:
                    self.showError(
                        "\n\nThere should be at least 1 model component!\n\n")

                self.models = models
#      self.wavelength=wavelength
                self.imsize = imsize
                self.imfiles = imfiles

                if not fixsize:
                    self.imsize = Xmax * 1.1

                self.Xaxmax = self.imsize / 2.

                fi.close()

        return True

    def _changeWavelength(self, wave, redoUV=False):

        if not self.GUIres:
            return

        self.wavelength[2] = wave * 1.e-6
        fmtB1 = r'$\lambda = $ %4.1fmm  ' % (self.wavelength[2] * 1.e6)
        self.fmtB = fmtB1 + \
            "\n" r'% 4.2f Jy/beam' "\n" r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f '

      #  self._plotAntennas(redo=False)
        self._setPrimaryBeam(replotFFT=True)
        self._changeCoordinates(rescale=True, redoUV=redoUV)
      #  self._plotModelFFT(redo=False)

    def _changeCoordinates(self, rescale=False, redoUV=False):

        if self.lat > np.pi / 2.:
            self.lat = np.pi / 2.
            return
        elif self.lat < -np.pi / 2.:
            self.lat = -np.pi / 2.
            return

        if self.dec > np.pi / 2.:
            self.dec = np.pi / 2.
            return
        elif self.dec < -np.pi / 2.:
            self.dec = -np.pi / 2.
            return

        self.trlat = [np.sin(self.lat), np.cos(self.lat)]
        self.trdec = [np.sin(self.dec), np.cos(self.dec)]

        cosW = -np.tan(self.lat) * np.tan(self.dec)
        if np.abs(cosW) < 1.0:
            Hhor = np.arccos(cosW)
        elif np.abs(self.lat - self.dec) > np.pi / 2.:
            Hhor = 0
        else:
            Hhor = np.pi

        self.Hmax = Hhor

        if self.Hmax > 0.0:
            if self.Hcov[0] < -self.Hmax:
                self.Hcov[0] = -self.Hmax
                self.lock = True
                self.widget['H0'].set_val(self.Hcov[0] / self.Hfac)
                self.lock = False
            if self.Hcov[1] > self.Hmax:
                self.Hcov[1] = self.Hmax
                self.lock = True
                self.widget['H1'].set_val(self.Hcov[1] / self.Hfac)
                self.lock = False

        if redoUV:
            self.UVPlot.cla()
            self._plotModelFFT(redo=True)
            self._plotAntennas(redo=True, rescale=True)

        newtext = self.fmtH % (
            self.lat / self.deg2rad, self.dec / self.deg2rad,
                               self.Hcov[0] / self.Hfac, self.Hcov[1] / self.Hfac)
        self.latText.set_text(newtext)
        self.Horig = np.linspace(self.Hcov[0], self.Hcov[1], self.nH)
        H = self.Horig[np.newaxis, :]
        self.H = [np.sin(H), np.cos(H)]
        self._setBaselines()
        self._setBeam()
        self._plotBeam(redo=False)
        self._plotAntennas(redo=False, rescale=True)
        self._plotDirty(redo=False)

        self.beamText.set_text(self.fmtB % (1.0, 0.0, 0.0))
        newtext = self.fmtVis % (self.totflux, 0.0)
        self.visText.set_text(newtext)
        dirflux = self.dirtymap[self.Nphf, self.Nphf]
        modflux = self.modelimTrue[self.Nphf, self.Nphf]
        self.dirtyText.set_text(self.fmtD % (dirflux, 0.0, 0.0))
        self.modelText.set_text(self.fmtM % (modflux, 0.0, 0.0))
        self.basText.set_text(self.fmtBas % (0, 0, 0.0))
        self.antPlotBas.set_data([[0], [0]])

        pl.draw()
        self.canvas.draw()

    def _setNoise(self, noise):
        if noise == 0.0:
            self.Noise[:] = 0.0
        else:
            self.Noise[:] = (np.random.normal(loc=0.0, scale=noise,
                                              size=np.shape(self.Noise)) +
                             1.j * np.random.normal(loc=0.0, scale=noise,
                                                    size=np.shape(self.Noise)))
        self._setBaselines()
        self._setBeam()
        self._plotBeam(redo=False)
        self._plotDirty(redo=False)
        self.canvas.draw()

    def _setGains(self, An1, An2, H0, H1, G):

        self.Gains[:] = 1.0

        for nb in range(self.Nbas):
            if An1 == self.antnum[nb][0]:
                if An2 == -1 or An2 == self.antnum[nb][1]:
                    self.Gains[nb, H0:H1] *= G
            if An1 == self.antnum[nb][1]:
                if An2 == -1 or An2 == self.antnum[nb][0]:
                    self.Gains[nb, H0:H1] *= np.conjugate(G)

        self._setBaselines()
        self._setBeam()
        self._plotBeam(redo=False)
        self._plotDirty(redo=False)
        self.canvas.draw()

    def _prepareBeam(self):

        self.beam = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.totsampling = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.dirtymap = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.noisemap = np.zeros((self.Npix, self.Npix), dtype=np.complex64)
        self.robustsamp = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.Gsampling = np.zeros((self.Npix, self.Npix), dtype=np.complex64)
        self.Grobustsamp = np.zeros((self.Npix, self.Npix), dtype=np.complex64)
        self.GrobustNoise = np.zeros(
            (self.Npix, self.Npix), dtype=np.complex64)

        self.beam2 = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.totsampling2 = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.dirtymap2 = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.robustsamp2 = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        # self.Gsampling2 = np.zeros((self.Npix,self.Npix),dtype=np.complex64)
        # self.Grobustsamp2 =
        # np.zeros((self.Npix,self.Npix),dtype=np.complex64)

    def _prepareBaselines(self):

        self.Nbas = int(self.Nant * (self.Nant - 1) / 2)
        NBmax = self.Nbas
        self.B = np.zeros((NBmax, self.nH), dtype=np.float32)
        self.basnum = np.zeros((self.Nant, self.Nant - 1), dtype=np.int8)
        self.basidx = np.zeros((self.Nant, self.Nant), dtype=np.int8)
        self.antnum = np.zeros((NBmax, 2), dtype=np.int8)
        self.Gains = np.ones((self.Nbas, self.nH), dtype=np.complex64)
        self.Noise = np.zeros((self.Nbas, self.nH), dtype=np.complex64)
        self.Horig = np.linspace(self.Hcov[0], self.Hcov[1], self.nH)
        H = self.Horig[np.newaxis, :]
        self.H = [np.sin(H), np.cos(H)]

        bi = 0
        nii = [0 for n in range(self.Nant)]
        for n1 in range(self.Nant - 1):
            for n2 in range(n1 + 1, self.Nant):
                self.basnum[n1, nii[n1]] = bi
                self.basnum[n2, nii[n2]] = bi
                self.basidx[n1, n2] = bi
                self.antnum[bi] = [n1, n2]
                nii[n1] += 1
                nii[n2] += 1
                bi += 1

        self.u = np.zeros((NBmax, self.nH))
        self.v = np.zeros((NBmax, self.nH))
        self.ravelDims = (NBmax, self.nH)

        if self.Nant2 > 1:
            self.Nbas2 = int(self.Nant2 * (self.Nant2 - 1) / 2)
            NBmax2 = self.Nbas2
            self.B2 = np.zeros((NBmax2, self.nH), dtype=np.float32)
            self.basnum2 = np.zeros(
                (self.Nant2, self.Nant2 - 1), dtype=np.int8)
            self.basidx2 = np.zeros((self.Nant2, self.Nant2), dtype=np.int8)
            self.antnum2 = np.zeros((NBmax2, 2), dtype=np.int8)
            self.Gains2 = np.ones((self.Nbas2, self.nH), dtype=np.complex64)
            self.H = [np.sin(H), np.cos(H)]

            bi = 0
            nii = [0 for n in range(self.Nant2)]
            for n1 in range(self.Nant2 - 1):
                for n2 in range(n1 + 1, self.Nant2):
                    self.basnum2[n1, nii[n1]] = bi
                    self.basnum2[n2, nii[n2]] = bi
                    self.basidx2[n1, n2] = bi
                    self.antnum2[bi] = [n1, n2]
                    nii[n1] += 1
                    nii[n2] += 1
                    bi += 1

            self.u2 = np.zeros((NBmax2, self.nH))
            self.v2 = np.zeros((NBmax2, self.nH))
            self.ravelDims2 = (NBmax2, self.nH)

    def _setBaselines(self, antidx=-1):

        if antidx == -1:
            bas2change = range(self.Nbas)
        elif antidx < self.Nant:
            bas2change = self.basnum[antidx].flatten()
        else:
            bas2change = []

        for currBas in bas2change:
            n1, n2 = self.antnum[currBas]
            self.B[currBas, 0] = (-(self.antPos[n2][1] - self.antPos[n1][1]) *
                                  self.trlat[0] / self.wavelength[2])
            self.B[currBas, 1] = (
                self.antPos[n2][0] - self.antPos[n1][0]) / self.wavelength[2]
            self.B[currBas, 2] = ((self.antPos[n2][1] - self.antPos[n1][1]) *
                                  self.trlat[1] / self.wavelength[2])
            self.u[currBas, :] = - \
                (self.B[currBas, 0] * self.H[
                 0] + self.B[currBas, 1] * self.H[1])
            self.v[currBas, :] = (-self.B[currBas, 0] * self.trdec[0] * self.H[1] + self.B[currBas, 1] *
                                  self.trdec[0] * self.H[0] + self.trdec[1] * self.B[currBas, 2])

        if self.Nant2 > 1:

            if antidx == -1:
                bas2change = range(self.Nbas2)
            elif antidx >= self.Nant:
                bas2change = self.basnum2[antidx - self.Nant].flatten()
            else:
                bas2change = []

            for currBas in bas2change:
                n1, n2 = self.antnum2[currBas]
                self.B2[currBas, 0] = (-(self.antPos2[n2][1] - self.antPos2[n1][1]) *
                                       self.trlat[0] / self.wavelength[2])
                self.B2[currBas, 1] = ((self.antPos2[n2][0] - self.antPos2[n1][0]) /
                                       self.wavelength[2])
                self.B2[currBas, 2] = ((self.antPos2[n2][1] - self.antPos2[n1][1]) *
                                       self.trlat[1] / self.wavelength[2])
                self.u2[currBas, :] = (-(self.B2[currBas, 0] * self.H[0] +
                                         self.B2[currBas, 1] * self.H[1]))
                self.v2[currBas, :] = ((-self.B2[currBas, 0] * self.trdec[0] * self.H[1] +
                                        self.B2[currBas, 1] * self.trdec[0] * self.H[0] +
                                        self.trdec[1] * self.B2[currBas, 2]))

    def _gridUV(self, antidx=-1):

        if antidx == -1:
            bas2change = range(self.Nbas)
            self.pixpos = [[] for nb in bas2change]
            self.totsampling[:] = 0.0
            self.Gsampling[:] = 0.0
            self.noisemap[:] = 0.0
        elif antidx < self.Nant:
            bas2change = map(int, list(self.basnum[antidx].flatten()))
        else:
            bas2change = []

        self.UVpixsize = 2. / (self.imsize * np.pi / 180. / 3600.)

        for nb in bas2change:
            pixU = np.rint(
                self.u[nb] / self.UVpixsize).flatten().astype(np.int32)
            pixV = np.rint(
                self.v[nb] / self.UVpixsize).flatten().astype(np.int32)
            goodpix = np.where(
                np.logical_and(np.abs(pixU) < self.Nphf, np.abs(pixV) < self.Nphf))[0]
            pU = pixU[goodpix] + self.Nphf
            pV = pixV[goodpix] + self.Nphf
            mU = -pixU[goodpix] + self.Nphf
            mV = -pixV[goodpix] + self.Nphf

            if not antidx == -1:
            # print bas2change
            # print np.shape(goodpix), np.shape(self.Gains),
            # np.shape(self.pixpos[nb][0]), nb
                self.totsampling[self.pixpos[nb][1], self.pixpos[nb][2]] -= 1.0
                self.totsampling[self.pixpos[nb][3], self.pixpos[nb][0]] -= 1.0
                self.Gsampling[
                    self.pixpos[nb][1], self.pixpos[nb][2]] -= self.Gains[nb, goodpix]
                self.Gsampling[self.pixpos[nb][3], self.pixpos[
                    nb][0]] -= np.conjugate(self.Gains[nb, goodpix])
                self.noisemap[self.pixpos[nb][1], self.pixpos[nb][2]] -= (self.Noise[nb, goodpix] *
                                                                          np.abs(self.Gains[nb, goodpix]))
                self.noisemap[self.pixpos[nb][3], self.pixpos[nb][0]] -= (np.conjugate(self.Noise[nb, goodpix]) *
                                                                          np.abs(self.Gains[nb, goodpix]))

            self.pixpos[nb] = [
                np.copy(pU), np.copy(pV), np.copy(mU), np.copy(mV)]
            for pi, gp in enumerate(goodpix):
                gabs = np.abs(self.Gains[nb, gp])
                pVi = pV[pi]
                mUi = mU[pi]
                mVi = mV[pi]
                pUi = pU[pi]
                self.totsampling[pVi, mUi] += 1.0
                self.totsampling[mVi, pUi] += 1.0
                self.Gsampling[pVi, mUi] += self.Gains[nb, gp]
                self.Gsampling[mVi, pUi] += np.conjugate(self.Gains[nb, gp])
                self.noisemap[pVi, mUi] += self.Noise[nb, gp] * gabs
                self.noisemap[mVi, pUi] += np.conjugate(
                    self.Noise[nb, gp]) * gabs

        self.robfac = (5. * 10.**(-self.robust))**2. * (
            2. * self.Nbas * self.nH) / np.sum(self.totsampling**2.)

        if self.Nant2 > 1:

            if antidx == -1:
                bas2change = range(self.Nbas2)
                self.pixpos2 = [[] for nb in bas2change]
                self.totsampling2[:] = 0.0
            #   self.Gsampling2[:] = 0.0
            elif antidx >= self.Nant:
                bas2change = map(
                    int, list(self.basnum2[antidx - self.Nant].flatten()))
            else:
                bas2change = []

            for nb in bas2change:
                pixU = np.rint(
                    self.u2[nb] / self.UVpixsize).flatten().astype(np.int32)
                pixV = np.rint(
                    self.v2[nb] / self.UVpixsize).flatten().astype(np.int32)
                goodpix = np.logical_and(
                    np.abs(pixU) < self.Nphf, np.abs(pixV) < self.Nphf)
                pU = pixU[goodpix] + self.Nphf
                pV = pixV[goodpix] + self.Nphf
                mU = -pixU[goodpix] + self.Nphf
                mV = -pixV[goodpix] + self.Nphf
                if not antidx == -1:
                    self.totsampling2[
                        self.pixpos2[nb][1], self.pixpos2[nb][2]] -= 1.0
                    self.totsampling2[
                        self.pixpos2[nb][3], self.pixpos2[nb][0]] -= 1.0
                    # self.Gsampling2[self.pixpos2[nb][1],self.pixpos2[nb][2]] -= self.Gains[nb,goodpix]
                    # self.Gsampling2[self.pixpos2[nb][3],self.pixpos2[nb][0]]
                    # -= np.conjugate(self.Gains[nb,goodpix])
                self.pixpos2[nb] = [
                    np.copy(pU), np.copy(pV), np.copy(mU), np.copy(mV)]

                self.totsampling2[pV, mU] += 1.0
                self.totsampling2[mV, pU] += 1.0
        #    self.Gsampling2[pV,mU] += self.Gains[nb,goodpix]
        #    self.Gsampling2[mV,pU] += np.conjugate(self.Gains[nb,goodpix])

            self.robfac2 = (5. * 10.**(-self.robust))**2. * (
                2. * self.Nbas2 * self.nH) / np.sum(self.totsampling2**2.)

    def _setBeam(self, antidx=-1):

        self._gridUV(antidx=antidx)

        denom = 1. + self.robfac * self.totsampling
        self.robustsamp[:] = self.totsampling / denom
        self.Grobustsamp[:] = self.Gsampling / denom
        self.GrobustNoise[:] = self.noisemap / denom

        self.beam[:] = np.fft.ifftshift(
            np.fft.ifft2(np.fft.fftshift(self.robustsamp))).real / (1. + self.W2W1)
        # self.beamScale =
        # np.max(self.beam[self.Nphf:self.Nphf+1,self.Nphf:self.Nphf+1])

        if self.Nant2 > 1:
            self.robustsamp2[:] = self.totsampling2 / (
                1. + self.robfac2 * self.totsampling2)
            self.beam[:] += np.fft.ifftshift(
                np.fft.ifft2(np.fft.fftshift(self.robustsamp2))).real * (self.W2W1 / (1. + self.W2W1))
            self.beamScale2 = np.max(
                self.beam[self.Nphf:self.Nphf + 1, self.Nphf:self.Nphf + 1])
            self.beam[:] /= self.beamScale2
        else:
            self.beamScale = np.max(
                self.beam[self.Nphf:self.Nphf + 1, self.Nphf:self.Nphf + 1])
            self.beam[:] /= self.beamScale

    def _prepareModel(self):

        pixsize = float(self.imsize) / self.Npix
        xx = np.linspace(-self.imsize / 2., self.imsize / 2., self.Npix)
        yy = np.ones(self.Npix, dtype=np.float32)
        distmat = np.zeros((self.Npix, self.Npix), dtype=np.float32)
        self.modelim = [np.zeros((self.Npix, self.Npix), dtype=np.float32)
                        for i in [0, 1]]
        self.modelimTrue = np.zeros((self.Npix, self.Npix), dtype=np.float32)

        for model in self.models:
            xsh = -model[1]
            ysh = -model[2]
            xpix = np.rint(xsh / pixsize).astype(np.int32)
            ypix = np.rint(ysh / pixsize).astype(np.int32)
            centy = np.roll(xx, ypix)
            centx = np.roll(xx, xpix)
            distmat[:] = np.outer(centy**2., yy) + np.outer(yy, centx**2.)
            if model[0] == 'D':
                mask = np.logical_or(
                    distmat <= model[4]**2., distmat == np.min(distmat))
                self.modelimTrue[mask] += float(model[3]) / np.sum(mask)
            elif model[0] == 'G':
                gauss = np.exp(-distmat / (2. * model[4]**2.))
                self.modelimTrue[:] += float(model[3]) * gauss / np.sum(gauss)
            elif model[0] == 'P':
                if np.abs(xpix + self.Nphf) < self.Npix and np.abs(ypix + self.Nphf) < self.Npix:
                    yint = ypix + self.Nphf
                    xint = xpix + self.Nphf
                    self.modelimTrue[yint, xint] += float(model[3])

        for imfile in self.imfiles:
            if not os.path.exists(imfile[0]):
                imfile[0] = os.path.join(self.datadir, imfile[0])
                if not os.path.exists(imfile[0]):
                    self.showError(
                        'File %s does NOT exist. Cannot read the model!' % imfile[0])
                    return

            Np4 = self.Npix / 4
            img = plimg.imread(imfile[0]).astype(np.float32)
            dims = np.shape(img)
            d3 = min(2, dims[2])
            d1 = float(max(dims))
            avimg = np.average(img[:, :, :d3], axis=2)
            avimg -= np.min(avimg)
            avimg *= imfile[1] / np.max(avimg)
            if d1 == self.Nphf:
                sh0 = (self.Nphf - dims[0]) / 2
                sh1 = (self.Nphf - dims[1]) / 2
                self.modelimTrue[sh0 + Np4:sh0 + Np4 + dims[
                    0], sh1 + Np4:sh1 + Np4 + dims[1]] += zoomimg
            else:
                zoomimg = spndint.zoom(avimg, float(self.Nphf) / d1)
                zdims = np.shape(zoomimg)
                zd0 = min(zdims[0], self.Nphf)
                zd1 = min(zdims[1], self.Nphf)
                sh0 = (self.Nphf - zdims[0]) / 2
                sh1 = (self.Nphf - zdims[1]) / 2
                self.modelimTrue[
                    sh0 + Np4:sh0 + Np4 + zd0, sh1 + Np4:sh1 + Np4 + zd1] += zoomimg[:zd0, :zd1]

        self.modelimTrue[self.modelimTrue < 0.0] = 0.0
        xx = np.linspace(-self.imsize / 2., self.imsize / 2., self.Npix)
        yy = np.ones(self.Npix, dtype=np.float32)
        self.distmat = (
            -np.outer(xx**2., yy) - np.outer(yy, xx**2.)) * pixsize**2.
        self._setPrimaryBeam(replotFFT=True)

    def _setPrimaryBeam(self, replotFFT=False):

        if self.Diameters[0] > 0.0:
            PB = 2. * \
                (1220. * 180. / np.pi * 3600. * self.wavelength[
                 2] / self.Diameters[0] / 2.3548)**2.  # 2*sigma^2
            #  print PB, np.max(self.distmat),self.wavelength
            beamImg = np.exp(self.distmat / PB)
            self.modelim[0][:] = self.modelimTrue * beamImg
        else:
            self.modelim[0][:] = self.modelimTrue

        if self.Nant2 > 1:
            if self.Diameters[1] > 0.0:
                PB = 2. * \
                    (1220. * 180. / np.pi * 3600. * self.wavelength[
                     2] / self.Diameters[1] / 2.3548)**2.  # 2*sigma^2
                beamImg = np.exp(self.distmat / PB)
                self.modelim[1][:] = self.modelimTrue * beamImg
            else:
                self.modelim[1][:] = self.modelimTrue

        self.modelfft = np.fft.fft2(np.fft.fftshift(self.modelim[0]))
        self.modelfft2 = np.fft.fft2(np.fft.fftshift(self.modelim[1]))
        if replotFFT:
            self._plotModelFFT(redo=True)

    def _plotModel(self, redo=True):

        Np4 = self.Npix / 4

        if redo:
            self.modelPlot.cla()
            self.modelPlotPlot = self.modelPlot.imshow(
                np.power(
                    self.modelimTrue[
                        Np4:self.Npix - Np4, Np4:self.Npix - Np4], self.gamma),
              picker=True, interpolation='nearest', vmin=0.0,
              vmax=np.max(self.modelimTrue)**self.gamma, cmap=self.currcmap)

            modflux = self.modelimTrue[self.Nphf, self.Nphf]
            self.modelText = self.modelPlot.text(0.05, 0.87,
                                                 self.fmtM % (
                                                     modflux, 0.0, 0.0),
                                                 transform=self.modelPlot.transAxes,
                                                 bbox=dict(facecolor='white', alpha=0.7))
            pl.setp(self.modelPlotPlot,
                    extent=(self.Xaxmax / 2., -self.Xaxmax / 2., -self.Xaxmax / 2., self.Xaxmax / 2.))
            self.modelPlot.set_ylabel('Dec offset (as)')
            self.modelPlot.set_xlabel('RA offset (as)')
            self._plotAntennas(redo=False)
        else:
            self.modelPlotPlot.set_data(
                np.power(
                    self.modelimTrue[Np4:self.Npix - Np4, Np4:self.Npix - Np4],
                       self.gamma))
            extr = [0.0, np.max(self.modelimTrue)**self.gamma]
            self.modelPlotPlot.norm.vmin = extr[0]
            self.modelPlotPlot.norm.vmax = extr[1]
            pl.setp(self.modelPlotPlot,
                    extent=(self.Xaxmax / 2., -self.Xaxmax / 2., -self.Xaxmax / 2., self.Xaxmax / 2.))

        self.totflux = np.sum(
            self.modelimTrue[Np4:self.Npix - Np4, Np4:self.Npix - Np4])
        self.modelPlot.set_title('MODEL IMAGE: %.2e Jy' % self.totflux)

    def _plotModelFFT(self, redo=True):

        self.UVmax = self.Npix / 2. / self.lfac * self.UVpixsize
        self.UVSh = -self.UVmax / self.Npix
        self.FFTtoplot = np.fft.fftshift(self.modelfft)
        toplot = np.abs(self.FFTtoplot)
        mval = np.min(toplot)
        Mval = np.max(toplot)
        dval = (Mval - mval) / 2.

        if redo:
            mymap = pl.gray()
            self.UVPlotFFTPlot = self.UVPlot.imshow(
                toplot, cmap=mymap, vmin=0.0,
                                                    vmax=Mval + dval, picker=5)
            pl.setp(self.UVPlotFFTPlot,
                    extent=(-self.UVmax + self.UVSh, self.UVmax + self.UVSh,
                            -self.UVmax - self.UVSh, self.UVmax - self.UVSh))
        else:
            self.UVPlotFFTPlot.set_data(toplot)
            self.UVPlotFFTPlot.norm.vmin = mval - dval
            self.UVPlotFFTPlot.norm.vmax = Mval + dval

    def _plotDirty(self, redo=True):
        Np4 = self.Npix / 4

        self.dirtymap[:] = (np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(self.GrobustNoise) +
                                                         self.modelfft * np.fft.ifftshift(self.Grobustsamp)))).real / (1. + self.W2W1)

    # print 'RMS:
    # ',np.std(np.abs(self.dirtymap[:])),np.max(np.abs(self.GrobustNoise)),np.max(np.abs(self.totsampling))

        if self.Nant2 > 1:
            self.dirtymap[:] += (np.fft.fftshift(np.fft.ifft2(self.modelfft2 * np.fft.ifftshift(self.robustsamp2)))).real * (
                self.W2W1 / (1. + self.W2W1))
            self.dirtymap /= self.beamScale2
        else:
            self.dirtymap /= self.beamScale

    # print 'RMS2: ',np.std(np.abs(self.dirtymap[:]))  #, self.beamScale

        extr = [np.min(self.dirtymap), np.max(self.dirtymap)]
        if redo:
            self.dirtyPlot.cla()
            self.dirtyPlotPlot = self.dirtyPlot.imshow(
                self.dirtymap[Np4:self.Npix - Np4,
                              Np4:self.Npix - Np4],
                                                       interpolation='nearest', picker=True,
                                                       cmap=self.currcmap)
            modflux = self.dirtymap[self.Nphf, self.Nphf]
            self.dirtyText = self.dirtyPlot.text(
                0.05, 0.87, self.fmtD % (modflux, 0.0, 0.0),
                                                 transform=self.dirtyPlot.transAxes,
                                                 bbox=dict(facecolor='white', alpha=0.7))
            pl.setp(self.dirtyPlotPlot,
                    extent=(self.Xaxmax / 2., -self.Xaxmax / 2., -self.Xaxmax / 2., self.Xaxmax / 2.))
            self.curzoom[1] = (
                self.Xaxmax / 2., -self.Xaxmax / 2., -self.Xaxmax / 2., self.Xaxmax / 2.)
            self.dirtyPlot.set_ylabel('Dec offset (as)')
            self.dirtyPlot.set_xlabel('RA offset (as)')
            self.dirtyPlot.set_title('DIRTY IMAGE')
        else:
            self.dirtyPlotPlot.set_data(
                self.dirtymap[Np4:self.Npix - Np4, Np4:self.Npix - Np4])
            self.dirtyPlotPlot.norm.vmin = extr[0]
            self.dirtyPlotPlot.norm.vmax = extr[1]

    def _plotAntennas(self, redo=True, rescale=False):

        mw = 2. * self.Xmax / self.wavelength[2] / self.lfac
        if mw < 0.1 and self.lfac == 1.e6:
            self.lfac = 1.e3
            self.ulab = r'U (k$\lambda$)'
            self.vlab = r'V (k$\lambda$)'
        elif mw >= 100. and self.lfac == 1.e3:
            self.lfac = 1.e6
            self.ulab = r'U (M$\lambda$)'
            self.vlab = r'V (M$\lambda$)'

        if redo:

            toplot = np.array(self.antPos[:self.Nant])
            self.antPlot.cla()
            if self.Nant2 > 1:
                pl.setp(self.wax['subarrwgt'], visible=True)
            else:
                pl.setp(self.wax['subarrwgt'], visible=False)
            self.antPlotBas = self.antPlot.plot([0], [0], '-b')[0]
            self.antPlotPlot = self.antPlot.plot(toplot[:, 0], toplot[:, 1],
                                                 'o', color='lime', picker=5)[0]
            if self.Nant2 > 1:
                toplot2 = np.array(self.antPos2[:self.Nant2])
                self.antPlotPlot2 = self.antPlot.plot(
                    toplot2[:, 0], toplot2[:, 1], 'or', picker=5)[0]

            self.antPlot.set_xlim((-self.Xmax, self.Xmax))
            self.antPlot.set_ylim((-self.Xmax, self.Xmax))
            self.curzoom[3] = (-self.Xmax, self.Xmax, -self.Xmax, self.Xmax)
            self.antPlot.set_xlabel('E-W offset (km)')
            self.antPlot.set_ylabel('N-S offset (km)')
            self.antPlot.set_title('ARRAY CONFIGURATION')
            self.antText = self.antPlot.text(0.05, 0.88,
                                             self.fmtA % (
                                                 self.Nant + self.Nant2),
                                             transform=self.antPlot.transAxes)
            self.UVPlotPlot = []
            toplotu = self.u.flatten() / self.lfac
            toplotv = self.v.flatten() / self.lfac
            self.UVPlotPlot.append(self.UVPlot.plot(toplotu, toplotv, '.',
                                                    color='lime', markersize=1, picker=2)[0])
            self.UVPlotPlot.append(self.UVPlot.plot(-toplotu, -toplotv, '.',
                                                    color='lime', markersize=1, picker=2)[0])
            if self.Nant2 > 1:
                self.UVPlotPlot2 = []
                toplotu = self.u2.flatten() / self.lfac
                toplotv = self.v2.flatten() / self.lfac
                self.UVPlotPlot2.append(self.UVPlot.plot(toplotu, toplotv,
                                                         '.r', markersize=1, picker=2)[0])
                self.UVPlotPlot2.append(self.UVPlot.plot(-toplotu, -toplotv,
                                                         '.r', markersize=1, picker=2)[0])
            self.UVPlot.set_xlim(
                (2. * self.Xmax / self.wavelength[2] / self.lfac,
                 -2. * self.Xmax / self.wavelength[2] / self.lfac))
            self.UVPlot.set_ylim(
                (2. * self.Xmax / self.wavelength[2] / self.lfac,
                 -2. * self.Xmax / self.wavelength[2] / self.lfac))
            self.curzoom[2] = (
                2. * self.Xmax / self.lfac, -2. * self.Xmax / self.lfac,
                               2. * self.Xmax / self.lfac, -2. * self.Xmax / self.lfac)
            self.latText = self.UVPlot.text(0.05, 0.87,
                                            self.fmtH % (
                                                self.lat / self.deg2rad,
                                                         self.dec /
                                                             self.deg2rad,
                                                         self.Hcov[
                                                             0] / self.Hfac,
                                                         self.Hcov[
                                                             1] / self.Hfac),
                                            transform=self.UVPlot.transAxes)
            self.latText.set_color('orange')
            self.basText = self.UVPlot.text(0.05, 0.02,
                                            self.fmtBas % (0, 0, 0.0),
                                            transform=self.UVPlot.transAxes)
            self.antPlotBas.set_data([[0], [0]])

            self.visText = self.UVPlot.text(0.05, 0.08,
                                            self.fmtVis % (0.0, 0.0),
                                            transform=self.UVPlot.transAxes)
            self.visText.set_color('orange')

            self.basText.set_color('orange')
            self.UVPlot.set_xlabel(self.ulab)
            self.UVPlot.set_ylabel(self.vlab)
            self.UVPlot.set_title('UV PLANE')

            self.antLabelPlot = []
            self.antLabelPlot2 = []

            for i in range(self.Nant):
                self.antLabelPlot.append(
                    self.antPlot.annotate(
                        str(i + 1), textcoords='offset points',
                                                               xy=(toplot[
                                                                   i, 0], toplot[
                                                                       i, 1]),
                                                               xytext=(-7, 4)))

            if self.Nant2 > 1:
                for i in range(self.Nant2):
                    self.antLabelPlot2.append(
                        self.antPlot.annotate(str(i + 1 + self.Nant),
                                              textcoords='offset points',
                                              xy=(toplot[
                                                  i, 0], toplot[
                                                  i, 1]),
                                              xytext=(-7, 4)))

        else:

            if rescale:
                self.antPlot.set_xlim((-self.Xmax, self.Xmax))
                self.antPlot.set_ylim((-self.Xmax, self.Xmax))
                self.UVPlot.set_xlim(
                    (2. * self.Xmax / self.wavelength[2] / self.lfac,
                     -2. * self.Xmax / self.wavelength[2] / self.lfac))
                self.UVPlot.set_ylim(
                    (2. * self.Xmax / self.wavelength[2] / self.lfac,
                     -2. * self.Xmax / self.wavelength[2] / self.lfac))
                self.curzoom[2] = (
                    2. * self.Xmax / self.lfac, -2. * self.Xmax / self.lfac,
                                   2. * self.Xmax / self.lfac, -2. * self.Xmax / self.lfac)
                self.curzoom[3] = (
                    -self.Xmax, self.Xmax, -self.Xmax, self.Xmax)

            if len(self.antLabelPlot) > self.Nant:
                for i in range(self.Nant, len(self.antLabelPlot)):
                    self.antLabelPlot[i].set_visible(False)

            toplot = np.array(self.antPos[:self.Nant])
            self.antPlotPlot.set_data(toplot[:, 0], toplot[:, 1])
            toplotu = self.u.flatten() / self.lfac
            toplotv = self.v.flatten() / self.lfac
            for i in range(self.Nant):
                if i > len(self.antLabelPlot) - 1:
                    self.antLabelPlot.append(self.antPlot.annotate(str(i + 1),
                                                                   textcoords='offset points',
                                                                   xy=(toplot[
                                                                       i, 0], toplot[
                                                                           i, 1]),
                                                                   xytext=(-7, 4)))
                else:
                    self.antLabelPlot[i].set_visible(True)
                    self.antLabelPlot[i].xy = (toplot[i, 0], toplot[i, 1])
            self.UVPlotPlot[0].set_data(toplotu, toplotv)
            self.UVPlotPlot[1].set_data(-toplotu, -toplotv)

            if self.Nant2 > 1:
                toplot = np.array(self.antPos2[:self.Nant2])
                self.antPlotPlot2.set_data(toplot[:, 0], toplot[:, 1])
                toplotu = self.u2.flatten() / self.lfac
                toplotv = self.v2.flatten() / self.lfac
                for i in range(self.Nant2):
                    self.antLabelPlot2[i].xy = (toplot[i, 0], toplot[i, 1])
                self.UVPlotPlot2[0].set_data(toplotu, toplotv)
                self.UVPlotPlot2[1].set_data(-toplotu, -toplotv)

            self.UVPlot.set_xlabel(self.ulab)
            self.UVPlot.set_ylabel(self.vlab)

    def _plotBeam(self, redo=True):

        if not self.show_beamPlot:
            return

        Np4 = self.Npix / 4
        if redo:
            self.beamPlot.cla()
            self.beamPlotPlot = self.beamPlot.imshow(
                self.beam[Np4:self.Npix - Np4,
                          Np4:self.Npix - Np4],
                                                     picker=True, interpolation='nearest',
                                                     cmap=self.currcmap)
            self.beamText = self.beamPlot.text(
                0.05, 0.80, self.fmtB % (1.0, 0.0, 0.0),
                                               transform=self.beamPlot.transAxes,
                                               bbox=dict(facecolor='white', alpha=0.7))
            self.beamPlot.set_ylabel('Dec offset (as)')
            self.beamPlot.set_xlabel('RA offset (as)')
            pl.setp(self.beamPlotPlot, extent=(self.Xaxmax / 2.,
                                               -self.Xaxmax / 2.,
                                               -self.Xaxmax / 2.,
                                               self.Xaxmax / 2.))
            self.curzoom[0] = (self.Xaxmax / 2., -self.Xaxmax / 2.,
                               -self.Xaxmax / 2., self.Xaxmax / 2.)
            self.beamPlot.set_title('DIRTY BEAM')
            pl.draw()
            self.canvas.draw()
        else:
            self.beamPlotPlot.set_data(self.beam[Np4:self.Npix - Np4,
                                                 Np4:self.Npix - Np4])
            self.beamText.set_text(self.fmtB % (1.0, 0.0, 0.0))

        self.nptot = np.sum(self.totsampling[:])
        self.beamPlotPlot.norm.vmin = np.min(self.beam)
        self.beamPlotPlot.norm.vmax = 1.0

        if np.sum(self.totsampling[self.Nphf - 4:self.Nphf + 4,
                                   self.Nphf - 4:self.Nphf + 4]) == self.nptot:
            warn = ('WARNING!\nToo short baselines for such a small image\n'
                    'PLEASE, INCREASE THE IMAGE SIZE!\nAND/OR DECREASE THE WAVELENGTH')
            self.beamText.set_text(warn)

        self.spherePlot.view_init(elev=self.dec / self.deg2rad, azim=0)
        self.arrayPath[0][:] = 10. * self.H[1] * np.cos(self.lat)
        self.arrayPath[1][:] = 10. * self.H[0] * np.cos(self.lat)
        self.arrayPath[2][:] = 10. * np.sin(self.lat)
        self.sphereArray[0].set_data(self.arrayPath[0], self.arrayPath[1])
        self.sphereArray[0].set_3d_properties(self.arrayPath[2])

    def _onPick(self, event):

        if event.mouseevent.inaxes == self.UVPlot:

            Up = event.mouseevent.xdata - self.UVSh
            Vp = event.mouseevent.ydata + self.UVSh
            yi = np.floor((self.UVmax + Up) / (self.UVmax) * self.Npix / 2.)
            xi = np.floor((self.UVmax - Vp) / (self.UVmax) * self.Npix / 2.)
            Flux = self.FFTtoplot[xi, yi]
            Phas, Amp = np.angle(Flux, deg=True), np.abs(Flux)
            newtext = self.fmtVis % (Amp, Phas)
            self.visText.set_text(newtext)

            if event.artist in self.UVPlotPlot:
                idata = np.unravel_index(event.ind, self.ravelDims)
                if event.artist == self.UVPlotPlot[0]:
                    n1, n2 = self.antnum[idata[0][0]]
                else:
                    n2, n1 = self.antnum[idata[0][0]]

                H = self.Horig[idata[1][0]] / self.Hfac
                newtext = self.fmtBas % (n1 + 1, n2 + 1, H)
                self.basText.set_text(newtext)
                self.antPlotBas.set_data([[self.antPos[n1][0],
                                           self.antPos[n2][0]],
                                          [self.antPos[n1][1],
                                           self.antPos[n2][1]]])

            elif self.Nant2 > 1 and event.artist in self.UVPlotPlot2:
                idata = np.unravel_index(event.ind, self.ravelDims2)
                if event.artist == self.UVPlotPlot2[0]:
                    n1, n2 = self.antnum2[idata[0][0]]
                else:
                    n2, n1 = self.antnum2[idata[0][0]]

                H = self.Horig[idata[1][0]] / self.Hfac
                newtext = self.fmtBas % (
                    n1 + 1 + self.Nant, n2 + 1 + self.Nant, H)
                self.basText.set_text(newtext)
                self.antPlotBas.set_data([[self.antPos2[n1][0],
                                           self.antPos2[n2][0]],
                                          [self.antPos2[n1][1],
                                           self.antPos2[n2][1]]])

            pl.draw()
            self.canvas.draw()
            return

        #  self.antPlotBas.set_data([[0],[0]])
        #  pl.draw()
        #  self.canvas.draw()

        elif event.mouseevent.inaxes == self.beamPlot:

            RA = event.mouseevent.xdata
            Dec = event.mouseevent.ydata
            yi = np.floor((self.Xaxmax - RA) / (2. * self.Xaxmax) * self.Npix)
            xi = np.floor((self.Xaxmax - Dec) / (2. * self.Xaxmax) * self.Npix)
            Flux = self.beam[xi, yi]
            self.beamText.set_text(self.fmtB % (Flux, RA, Dec))
            pl.draw()
            self.canvas.draw()

        elif event.mouseevent.inaxes == self.dirtyPlot:

            RA = event.mouseevent.xdata
            Dec = event.mouseevent.ydata
            yi = np.floor((self.Xaxmax - RA) / (2. * self.Xaxmax) * self.Npix)
            xi = np.floor((self.Xaxmax - Dec) / (2. * self.Xaxmax) * self.Npix)
            Flux = self.dirtymap[xi, yi]
            self.dirtyText.set_text(self.fmtD % (Flux, RA, Dec))
            pl.draw()
            self.canvas.draw()

        elif event.mouseevent.inaxes == self.modelPlot:

            RA = event.mouseevent.xdata
            Dec = event.mouseevent.ydata
            yi = np.floor((self.Xaxmax - RA) / (2. * self.Xaxmax) * self.Npix)
            xi = np.floor((self.Xaxmax - Dec) / (2. * self.Xaxmax) * self.Npix)
            Flux = self.modelimTrue[xi, yi]
            self.modelText.set_text(self.fmtM % (Flux, RA, Dec))
            pl.draw()
            self.canvas.draw()

        elif event.mouseevent.inaxes == self.antPlot:

            if event.artist is self.antPlotPlot:
                self.pickSub = 0
            elif self.Nant2 > 1 and event.artist is self.antPlotPlot2:
                self.pickSub = 1

#   else:
            if event.mouseevent.button == 1 and not self.pickAnt:
                self.antidx = event.ind
                if len(self.antidx) > 1:
                    self.antidx = self.antidx[-1]
                self.pickAnt = True
                if self.pickSub == 0:
                    self.antText.set_text(self.fmtA % (self.Nant + self.Nant2) +
                                          self.fmtA2 % (self.antidx + 1) +
                                          self.fmtA3 % tuple([1000 * a for a in self.antPos[self.antidx]]))
                else:
                    self.antText.set_text(self.fmtA % (self.Nant + self.Nant2) +
                                          self.fmtA2 % (self.antidx + self.Nant + 1) +
                                          self.fmtA3 %
                                          tuple([1000 * a for a in self.antPos2[self.antidx]]))

                pl.draw()
                self.canvas.draw()

    def _onAntennaDrag(self, event):
        if self.pickAnt:
            if self.pickSub == 0:
                self.antPos[self.antidx] = [event.xdata, event.ydata]
                self.antText.set_text(self.fmtA % (self.Nant + self.Nant2) +
                                      self.fmtA2 % (self.antidx + 1) +
                                      self.fmtA3 %
                                      tuple([1000 * a for a in self.antPos[self.antidx]]))
                self._setBaselines(-1)  # antidx=self.antidx)
                self._plotAntennas(redo=False)
                self._setBeam(-1)  # antidx=self.antidx)
            else:
                self.antPos2[self.antidx] = [event.xdata, event.ydata]
                self.antText.set_text(self.fmtA % (self.Nant + self.Nant2) +
                                      self.fmtA2 % (self.antidx + self.Nant + 1) +
                                      self.fmtA3 %
                                      tuple([1000 * a for a in self.antPos2[self.antidx]]))
                self._setBaselines(-1)  # antidx=self.antidx+self.Nant)
                self._plotAntennas(redo=False)
                self._setBeam(-1)  # antidx=self.antidx+self.Nant)

            self._plotBeam(redo=False)
            self._plotDirty(redo=False)

        #    pl.draw()
            self.canvas.draw()

# Drag the sphere plot (to change source position)
        if self._onSphere:
            oldDec = self.dec / self.deg2rad
            newDec, _newH0 = self.spherePlot.elev, self.spherePlot.azim

        # Limits on declination:
            if np.abs(newDec) > 90.:
                self.spherePlot.view_init(elev=oldDec, azim=0.0)
                return

            newDec *= self.deg2rad
            if not self.lock:
                self.lock = True
                if newDec != self.dec and np.abs(newDec - self.lat) < np.pi / 2.:
                    self.widget['dec'].set_val(newDec / self.deg2rad)
                    self.spherePlot.view_init(
                        elev=newDec / self.deg2rad, azim=0.0)
                    self.dec = newDec
                    self._changeCoordinates()
                    self.lock = False
                else:
                    self.spherePlot.view_init(elev=oldDec, azim=0.0)
                    self.lock = False
            else:
                self.spherePlot.view_init(elev=oldDec, azim=0.0)

    def _onRelease(self, event):
        self._onSphere = False
        if self.pickAnt:
            self.pickAnt = False
            self.antText.set_text(self.fmtA % self.Nant)
            pl.draw()
            self.canvas.draw()

    def _onRobust(self, newrob):
        self.robust = newrob
        self._changeCoordinates()

    def _onKeyLat(self, newlat):
        if not self.GUIres:
            return

        newlat *= self.deg2rad
        if not self.lock:
            self.lock = True
            if newlat != self.lat:
                if np.abs(newlat - self.dec) < np.pi / 2.:
                    self.lat = newlat
                    self._changeCoordinates()
                else:
                    self.widget['lat'].set_val(self.lat / self.deg2rad)
            self.lock = False

    def _subarrwgt(self, w1w2):
        self.W2W1 = 10.**(-w1w2)
        self._changeCoordinates()

    def _gammacorr(self, gamma):
        self.gamma = gamma
        self._plotModel(redo=False)
        pl.draw()
        self.canvas.draw()

    def _onKeyDec(self, newdec):
        if not self.GUIres:
            return

        newdec *= self.deg2rad
        if not self.lock:
            self.lock = True
            if newdec != self.dec:
                if np.abs(newdec - self.lat) < np.pi / 2.:
                    self.dec = newdec
                    self._changeCoordinates()
                else:
                    self.widget['dec'].set_val(self.dec / self.deg2rad)
            self.lock = False

    def _onKeyH0(self, newH0):
        if not self.GUIres:
            return

        newH0 *= self.Hfac
        if not self.lock:
            self.lock = True
            if np.abs(newH0) < self.Hmax:
                self.Hcov[0] = newH0
            else:
                self.Hcov[0] = -self.Hmax
                self.widget['H0'].set_val(self.Hcov[0] / self.Hfac)
            if self.Hcov[1] < self.Hcov[0]:
                self.Hcov[1] = self.Hcov[0]
                self.widget['H1'].set_val(self.Hcov[1] / self.Hfac)
            self._changeCoordinates()
            self.lock = False

    def _onKeyH1(self, newH1):
        if not self.GUIres:
            return

        newH1 *= self.Hfac
        if not self.lock:
            self.lock = True
            if np.abs(newH1) < self.Hmax:
                self.Hcov[1] = newH1
            else:
                self.Hcov[1] = self.Hmax
                self.widget['H1'].set_val(self.Hcov[1] / self.Hfac)
            if self.Hcov[0] > self.Hcov[1]:
                self.Hcov[0] = self.Hcov[1]
                self.widget['H0'].set_val(self.Hcov[0] / self.Hfac)
            self._changeCoordinates()
            self.lock = False

    def _addAntenna(self, antenna):
        if not self.antLock:

            self.antLock = True

            if self.Nant >= len(self.antPos):
                self.antPos.append([0., 0.])
                self.Nant += 1
                self.antLabelPlot.append(
                    self.antPlot.annotate(str(self.Nant),
                                          textcoords='offset points',
                                          xy=(0, 0), xytext=(-7, 4)))
            else:
                self.antLabelPlot[self.Nant].xy = (self.antPos[self.Nant][0],
                                                   self.antPos[self.Nant][1])
                self.Nant += 1

            self.antLabelPlot[self.Nant - 1].set_visible(True)

            newtext = self.fmtA % self.Nant
            self.antText.set_text(newtext)
            self._prepareBaselines()
            self._changeCoordinates()

            self.antLock = False

    def _removeAntenna(self, antenna):

        if not self.antLock:

            self.antLock = True

            if self.Nant > 2:
                self.Nant -= 1
                for i in range(self.Nant, len(self.antLabelPlot)):
                    self.antLabelPlot[i].set_visible(False)
            newtext = self.fmtA % self.Nant
            self.antText.set_text(newtext)
            self._prepareBaselines()
            self._changeCoordinates()

            self.antLock = False

    def _onKeyPress(self, event):
        from uvplotter2 import ApSynSim_UV_Plotter2

        if event.key == 'u' or event.key == 'U':
            if self.tks is not None:
                self.myUVPLOT2 = ApSynSim_UV_Plotter2(self)

        if event.key == 'c' or event.key == 'C':
            if self.currcmap == cm.jet:
                self.currcmap = cm.Greys_r
            else:
                self.currcmap = cm.jet

            self._plotBeam(redo=True)
            self._plotModel(redo=True)
            self._plotDirty(redo=True)
            self._plotModelFFT(redo=True)

            pl.draw()
            self.canvas.draw()

            if self.my_cleaner:
                self.my_cleaner.ResidPlotPlot.set_cmap(self.currcmap)
                self.my_cleaner.CLEANPlotPlot.set_cmap(self.currcmap)
                self.my_cleaner.canvas1.draw()

        if event.key == 'Z':
            event.button = 1
            event.dblclick = True
            self._onPress(event)
        if event.key == 'z':
            event.button = 3
            event.dblclick = True
            self._onPress(event)

    def _onPress(self, event):
        if event.inaxes == self.spherePlot:
            self._onSphere = True

        if not hasattr(event, 'dblclick'):
            event.dblclick = False

        if event.dblclick:
            mpl = [self.modelPlot, self.dirtyPlot]
            if self.my_cleaner:
                mpl += [self.my_cleaner.ResidPlot, self.my_cleaner.CLEANPlot]

                # ZOOM IN:
            if event.inaxes == self.beamPlot:
                toZoom = [self.beamPlot]
                cz = 0
                inv = True
                inv2 = False
                scal = 1.0
            elif event.inaxes in mpl:  # [self.modelPlot, self.dirtyPlot]:
                toZoom = mpl  # [self.modelPlot, self.dirtyPlot]
                cz = 1
                inv = True
                inv2 = False
                scal = 1.0
            elif event.inaxes == self.UVPlot:
                toZoom = [self.UVPlot]
                cz = 2
                inv = True
                inv2 = True
                scal = self.wavelength[2]
            elif event.inaxes == self.antPlot:
                toZoom = [self.antPlot]
                cz = 3
                inv = False
                inv2 = False
                scal = 1.0
            else:
                cz = -1
                inv = False
                inv2 = False
                scal = 1.0

            if cz >= 0:
                if event.button == 1 and cz >= 0:
                    RA = event.xdata
                    Dec = event.ydata
                    xL = np.abs(
                        self.curzoom[cz][1] - self.curzoom[cz][0]) / 4. / scal
                    yL = np.abs(
                        self.curzoom[cz][3] - self.curzoom[cz][2]) / 4. / scal
                    x0 = RA - xL
                    x1 = RA + xL
                    y0 = Dec - xL
                    y1 = Dec + xL
                    if cz in [0, 1]:
                        if x0 < -self.Xaxmax / 2.:
                            x0 = -self.Xaxmax / 2.
                            x1 = x0 + 2. * xL
                        if x1 > self.Xaxmax / 2.:
                            x1 = self.Xaxmax / 2.
                            x0 = x1 - 2. * xL
                        if y0 < -self.Xaxmax / 2.:
                            y0 = -self.Xaxmax / 2.
                            y1 = y0 + 2. * xL
                        if y1 > self.Xaxmax / 2.:
                            y1 = self.Xaxmax / 2.
                            y0 = y1 - 2. * xL
                            # ZOOM OUT:
                if event.button == 3:
                    RA = event.xdata
                    Dec = event.ydata
                    xL = np.abs(
                        self.curzoom[cz][1] - self.curzoom[cz][0]) / scal
                    yL = np.abs(
                        self.curzoom[cz][3] - self.curzoom[cz][2]) / scal
                    if cz in [0, 1]:
                        if xL > self.Xaxmax / 2.:
                            xL = self.Xaxmax / 2.
                        if yL > self.Xaxmax / 2.:
                            yL = self.Xaxmax / 2.

                    x0 = RA - xL
                    x1 = RA + xL
                    y0 = Dec - xL
                    y1 = Dec + xL
                    if cz in [0, 1]:
                        if x0 < -self.Xaxmax / 2.:
                            x0 = -self.Xaxmax / 2.
                            x1 = x0 + 2. * xL
                        if x1 > self.Xaxmax / 2.:
                            x1 = self.Xaxmax / 2.
                            x0 = x1 - 2. * xL
                        if y0 < -self.Xaxmax / 2.:
                            y0 = -self.Xaxmax / 2.
                            y1 = y0 + 2. * xL
                        if y1 > self.Xaxmax / 2.:
                            y1 = self.Xaxmax / 2.
                            y0 = y1 - 2. * xL

                if inv:
                    xx0 = x1
                    xx1 = x0
                else:
                    xx0 = x0
                    xx1 = x1
                if inv2:
                    yy0 = y1
                    yy1 = y0
                else:
                    yy0 = y0
                    yy1 = y1

                for ax in toZoom:
                    ax.set_xlim((xx0, xx1))
                    ax.set_ylim((yy0, yy1))

                self.curzoom[cz] = (
                    xx0 * scal, xx1 * scal, yy0 * scal, yy1 * scal)

                pl.draw()
                self.canvas.draw()

                if self.my_cleaner:
                    self.my_cleaner.canvas1.draw()

    def saveArray(self, array):

        fname = tkFileDialog.asksaveasfilename(defaultextension='.array',
                                               title='Save current array...')
        iff = open(fname, 'w')

        print >> iff, 'LATITUDE % 3.1f' % (self.lat / self.deg2rad)
        print >> iff, 'DECLINATION % 3.1f' % (self.dec / self.deg2rad)
        toprint = tuple([l / self.Hfac for l in self.Hcov])
        print >> iff, 'HOUR_ANGLE % 3.1f  % 3.1f' % toprint

        if self.Diameters[0] != 0.0 or self.Diameters[1] != 0.0:
            print >> iff, 'DIAMETER % 3.1f  % 3.1f' % tuple(self.Diameters)

        for ant in self.antPos:
            toprint = tuple([p * 1.e3 for p in ant])
            print >> iff, 'ANTENNA  % .3e   % .3e' % toprint

        self.antText.set_text('SAVED: %s' % os.path.basename(fname))
        pl.draw()
        self.canvas.draw()

        time.sleep(3)
        self.antText.set_text(self.fmtA % self.Nant)
        pl.draw()
        self.canvas.draw()

        iff.close()

    def loadArray(self, array):

        antenna_file = tkFileDialog.askopenfilename(title='Load array...',
                                                    initialdir=self.arraydir)
        self.lock = False

        if len(antenna_file) > 0:
            goodread = self.readAntennas(str(antenna_file))

            if goodread:
                self.GUIres = False
                newtext = self.fmtA % self.Nant
                self.antText.set_text(newtext)
                self.widget['diameter'].set_val(self.Diameters[0])
                self.widget['lat'].set_val(self.lat / self.deg2rad)
                self.widget['dec'].set_val(self.dec / self.deg2rad)
                self.widget['H0'].set_val(self.Hcov[0] / self.Hfac)
                self.widget['H1'].set_val(self.Hcov[1] / self.Hfac)
                self.wax['wave'].cla()
                self.widget['wave'] = Slider(self.wax['wave'],
                                             r'$\lambda$ (mm)',
                                             self.wavelength[0] * 1.e6,
                                             self.wavelength[1] * 1.e6,
                                             valinit=self.wavelength[2] * 1.e6)
                self.widget['wave'].on_changed(self._changeWavelength)
                self.widget['wave'].set_val(self.wavelength[2] * 1.e6)
                self.GUIres = True

                self._prepareBaselines()
                self._setBaselines()
                self._plotModelFFT(redo=True)
                self._plotAntennas(redo=True, rescale=True)
                self._plotModel(redo=True)
                self._changeCoordinates(redoUV=True)
                self.widget['wave'].set_val(self.wavelength[2] * 1.e6)

                self._init_clean_img()

                pl.draw()
                self.canvas.draw()

    def loadModel(self, model):

        model_file = tkFileDialog.askopenfilename(title='Load model...',
                                                  initialdir=self.modeldir)
        self.lock = False

        if len(model_file) > 0:
            goodread = self.readModels(str(model_file))
            if goodread:
                self._prepareModel()
                self._plotModel(redo=True)
                self._setBaselines()
                self._setBeam()
                self._changeCoordinates()
                self._plotModelFFT(redo=True)
                self._plotBeam(redo=True)
                self._plotDirty(redo=True)

                self._init_clean_img()

                pl.draw()
                self.canvas.draw()
