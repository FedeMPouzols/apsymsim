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

import scipy.optimize as spfit

from ScrolledText import ScrolledText

try:
    import Tkinter as Tk
except:
    import tkinter as Tk

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)

from tkMessageBox import showinfo


__CLEAN_help_text__ = """
APSYNSIM - CLEAN GUI

Here you can experiment with CLEAN deconvolution on (noise-free)
visibilities. You can also corrupt the visibilities by adding a complex
gain to one of your antennas (or baselines).

Clicking and dragging, with the LEFT mouse button, on the RESIDUALS image
creates new CLEAN mask regions. Clicking and dragging with the RIGHT mouse
button removes CLEAN mask regions. You can add as many CLEAN mask regions
as you want.

The CLEAN gain and number of iterations can be changed in the text boxes.
Pressing CLEAN executes the iterations, refreshing all images in real time.
You can further click on CLEAN, to continue deconvolving. The box "Thres"
is the CLEAN threshold (in Jy per beam). Setting it to negative values will
allow CLEANing negative components.


Pressing RELOAD will undo all CLEANING and update the images from the
main window. That is, if you change anything in the main program window
(e.g., observing wavelength, antenna positions, etc.), pressing RELOAD
will apply such changes to the images in the CLEAN interface.

TIP: You can load more than one CLEAN GUI, change anything in the main
window and press "RELOAD" just in one of the GUIs. This way, you can
compare directly how the changes you made in the main window affect the
CLEANing!

Pressing "+/- Resid" will add (or remove) the residuals from the CLEANed
image.  By default, the residuals are NOT added (i.e., only the restored
CLEAN components are shown in the CLEAN image).

Pressing "(Un)restore" will restore (or unrestore) the CLEAN model with
the CLEAN beam when plotting. Default status is to apply the restore.

Pressing "Rescale" will rescale the color palette (e.g., to see better
the structure of the residuals).

Pressing "True source (conv.)" will show the true source structure
convolved with the CLEAN beam. This is to compare the fidelity of the
CLEAN deconvolution algorithm, by comparing the CLEAN image to the
true source brightness distribution (downgraded to the CLEAN resolution).

You can add random noise to your visibilities by setting a sensitivity
(in the "Sensit." text) and pressing "Redo Noise". Any time that you
press this button, a new realisation of the random noise will be
computed. "Sensit." is the expected rms that you would get from of a
source-free observation, using natural weighting. Basically, the noise
added to each visibility is proportional to Sensit.*sqrt(Nbas*Nt), where
Nbas is the number of baselines and Nt is the number of integration
times per baseline.

-----------------------------
HOW TO ADD A CORRUPTING GAIN
-----------------------------

Just select an antenna from the "Ant. 1" list to corrupt it. If you select
a different antenna from the "Ant. 2" list, only the baseline between
the two antennas will be corrupted. But if the two antennas are the same,
then ALL the baselines to that antenna will be corrupted.

The two first sliders ("From integration" and "to integration") mark the
first and last observing scans where the corruption term will be applied.
By default, the whole duration of the experiment is selected.

The last two sliders ("Amplitude gain" and "phase gain") define the gain
that will be applied to the corrupted antenna.

The button "APPLY GAIN" actually applies the gain and reloads the new
images.

The button "RESET GAIN", undoes the gain correction (so the data become
perfectly calibrated again).

NOTICE THAT if a new antenna is added, or subtracted, the gains are
reset automatically (but you will need to refresh the images in this
window, by pressing the "RESET" button, just below the "CLEAN"
button, to load the correct images).

"""


class Cleaner(object):

    def quit(self):

        self.parent.my_cleaner = None
        self.parent._setNoise(0.0)
        self.parent._setGains(-1, -1, 0, 0, 1.0)
        self.me.destroy()
        self.residuals[:] = 0.0
        self.cleanmod[:] = 0.0
        self.cleanmodd[:] = 0.0
        self.cleanBeam[:] = 0.0

    def __init__(self, parent):

        self.parent = parent
        self.me = Tk.Toplevel(parent.tks)

        menubar = Tk.Menu(self.me)
        menubar.add_command(label="Help", command=self._getHelp)
        menubar.add_command(label="Quit", command=self.quit)

        self.me.config(menu=menubar)
        self.me.protocol("WM_DELETE_WINDOW", self.quit)
        self.Np4 = self.parent.Npix / 4

        self.figCL1 = pl.figure(figsize=(12, 6))
        # self.figCL2 = pl.figure(figsize=(6,6))

        self.ResidPlot = self.figCL1.add_subplot(121, aspect='equal')
        # pl.axes([0.01,0.43,0.5,0.5],aspect='equal')
        self.CLEANPlot = self.figCL1.add_subplot(122, aspect='equal',
                                                 sharex=self.ResidPlot,
                                                 sharey=self.ResidPlot)
        # pl.axes([0.55,0.43,0.5,0.5],aspect='equal')
        self.ResidPlot.set_adjustable('box-forced')
        self.CLEANPlot.set_adjustable('box-forced')

        self.frames = {}
        self.frames['FigFr'] = Tk.Frame(self.me)
        self.frames['GFr'] = Tk.Frame(self.me)

        self.canvas1 = FigureCanvasTkAgg(
            self.figCL1, master=self.frames['FigFr'])
        # self.canvas2 = FigureCanvasTkAgg(self.figCL2,
        # master=self.frames['FigFr'])

        self.canvas1.show()
        # self.canvas2.show()

        self.frames['FigFr'].pack(side=Tk.TOP)
        self.frames['GFr'].pack(side=Tk.TOP)

        self.frames['CLOpt'] = Tk.Frame(self.frames['FigFr'])

        self.frames['Gain'] = Tk.Frame(self.frames['CLOpt'])
        self.frames['Niter'] = Tk.Frame(self.frames['CLOpt'])
        self.frames['Thres'] = Tk.Frame(self.frames['CLOpt'])
        self.frames['Sensit'] = Tk.Frame(self.frames['CLOpt'])

        Gtext = Tk.Label(self.frames['Gain'], text="Gain:  ")
        Ntext = Tk.Label(self.frames['Niter'], text="# iter:")
        Ttext = Tk.Label(self.frames['Thres'], text="Thres (Jy/b):")
        Stext = Tk.Label(self.frames['Sensit'], text="Sensit. (Jy/b):")

        self.entries = {}
        self.entries['Gain'] = Tk.Entry(self.frames['Gain'])
        self.entries['Gain'].insert(0, "0.1")
        self.entries['Gain'].config(width=5)

        self.entries['Niter'] = Tk.Entry(self.frames['Niter'])
        self.entries['Niter'].insert(0, "100")
        self.entries['Niter'].config(width=5)

        self.entries['Thres'] = Tk.Entry(self.frames['Thres'])
        self.entries['Thres'].insert(0, "0.0")
        self.entries['Thres'].config(width=5)

        self.entries['Sensit'] = Tk.Entry(self.frames['Sensit'])
        self.entries['Sensit'].insert(0, "0.0")
        self.entries['Sensit'].config(width=5)

        GTitle = Tk.Label(self.frames['GFr'], text="CALIBRATION ERROR:")
        GTitle.pack(side=Tk.TOP)
        self.frames['Ant1L'] = Tk.Frame(self.frames['GFr'])
        Ant1T = Tk.Label(self.frames['Ant1L'], text="Ant. 1:")
        Ant1T.pack(side=Tk.TOP)
        self.entries['Ant1'] = Tk.Listbox(self.frames['Ant1L'],
                                          exportselection=False,
                                          width=5)
        self.frames['Ant2L'] = Tk.Frame(self.frames['GFr'])
        Ant2T = Tk.Label(self.frames['Ant2L'], text="Ant. 2:")
        Ant2T.pack(side=Tk.TOP)
        self.entries['Ant2'] = Tk.Listbox(self.frames['Ant2L'],
                                          exportselection=False,
                                          width=5)
        self.entries['Ant1'].pack(side=Tk.TOP)
        self.entries['Ant2'].pack(side=Tk.TOP)

        self.frames['Ant1L'].pack(side=Tk.LEFT)
        self.frames['Ant2L'].pack(side=Tk.LEFT)

        self.frames['space1'] = Tk.Frame(self.frames['GFr'], width=40)
        self.frames['space1'].pack(side=Tk.LEFT)

        self.frames['H0Fr'] = Tk.Frame(self.frames['GFr'])
        self.frames['H1Fr'] = Tk.Frame(self.frames['GFr'])
        self.frames['AmpFr'] = Tk.Frame(self.frames['GFr'])
        self.frames['PhasFr'] = Tk.Frame(self.frames['GFr'])

        self.frames['H0Fr'].pack(side=Tk.TOP)
        self.frames['H1Fr'].pack(side=Tk.TOP)
        self.frames['AmpFr'].pack(side=Tk.TOP)
        self.frames['PhasFr'].pack(side=Tk.TOP)

        self.entries['H0'] = Tk.Scale(self.frames['H0Fr'], from_=0,
                                      to=self.parent.nH, orient=Tk.HORIZONTAL,
                                      length=200)
        self.entries['H1'] = Tk.Scale(self.frames['H1Fr'], from_=0,
                                      to=self.parent.nH, orient=Tk.HORIZONTAL,
                                      length=200)
        H0Text = Tk.Label(self.frames['H0Fr'], text="From integration #: ",
                          width=15)
        H1Text = Tk.Label(self.frames['H1Fr'], text="To integration #: ",
                          width=15)
        H0Text.pack(side=Tk.LEFT)
        self.entries['H0'].pack(side=Tk.RIGHT)
        H1Text.pack(side=Tk.LEFT)
        self.entries['H1'].pack(side=Tk.RIGHT)
        self.entries['H1'].set(self.parent.nH)

        self.entries['Amp'] = Tk.Scale(self.frames['AmpFr'],
                                       from_=0.1, to=1000.,
                                       orient=Tk.HORIZONTAL, length=200)
        self.entries['Phas'] = Tk.Scale(self.frames['PhasFr'],
                                        from_=-180., to=180.,
                                        orient=Tk.HORIZONTAL, length=200)
        AmpText = Tk.Label(self.frames['AmpFr'], text="Amplitude gain (%): ",
                           width=15)
        PhasText = Tk.Label(self.frames['PhasFr'], text="Phase Gain: ",
                            width=15)
        AmpText.pack(side=Tk.LEFT)
        self.entries['Amp'].pack(side=Tk.RIGHT)
        PhasText.pack(side=Tk.LEFT)
        self.entries['Phas'].pack(side=Tk.RIGHT)

        Gtext.pack(side=Tk.LEFT)
        self.entries['Gain'].pack(side=Tk.RIGHT)

        Ntext.pack(side=Tk.LEFT)
        self.entries['Niter'].pack(side=Tk.RIGHT)

        Ttext.pack(side=Tk.LEFT)
        self.entries['Thres'].pack(side=Tk.RIGHT)

        Stext.pack(side=Tk.LEFT)
        self.entries['Sensit'].pack(side=Tk.RIGHT)

        self.frames['CLOpt'].pack(side=Tk.LEFT)
        # self.canvas2.get_tk_widget().pack(side=Tk.LEFT) # , fill=Tk.BOTH,
        # expand=1)
        self.canvas1.get_tk_widget().pack(
            side=Tk.LEFT)  # , fill=Tk.BOTH, expand=1)

        self.buttons = {}
        self.buttons['Noise'] = Tk.Button(self.frames['CLOpt'],
                                          text="Redo Noise",
                                          command=self._ReNoise)
        self.buttons['clean'] = Tk.Button(self.frames['CLOpt'],
                                          text="CLEAN",
                                          command=self._CLEAN)
        self.buttons['reset'] = Tk.Button(self.frames['CLOpt'],
                                          text="RELOAD",
                                          command=self._reset)
        self.buttons['addres'] = Tk.Button(self.frames['CLOpt'],
                                           text="+/- Resid",
                                           command=self._AddRes)
        self.buttons['dorestore'] = Tk.Button(self.frames['CLOpt'],
                                              text="(Un)restore",
                                              command=self._doRestore)
        self.buttons['dorescale'] = Tk.Button(self.frames['CLOpt'],
                                              text="Rescale",
                                              command=self._doRescale)

        self.buttons['showfft'] = Tk.Button(self.frames['CLOpt'],
                                            text="Show FFT",
                                            command=self._showFFT)
        self.buttons['convsource'] = Tk.Button(self.frames['CLOpt'],
                                               text="True source (conv.)",
                                               command=self._convSource)

        self.buttons['apply'] = Tk.Button(self.frames['GFr'],
                                          text="APPLY GAIN",
                                          command=self._ApplyGain)
        self.buttons['apply'].pack(side=Tk.RIGHT)
        self.buttons['recal'] = Tk.Button(self.frames['GFr'],
                                          text="RESET GAIN",
                                          command=self._reCalib)
        self.buttons['recal'].pack(side=Tk.RIGHT)

        self.frames['Gain'].pack(side=Tk.TOP)
        self.frames['Niter'].pack(side=Tk.TOP)
        self.frames['Thres'].pack(side=Tk.TOP)

        self.buttons['clean'].pack(side=Tk.TOP)
        self.buttons['reset'].pack(side=Tk.TOP)
        self.buttons['addres'].pack(side=Tk.TOP)
        self.buttons['dorestore'].pack(side=Tk.TOP)
        self.buttons['dorescale'].pack(side=Tk.TOP)
        self.buttons['showfft'].pack(side=Tk.TOP)
        self.buttons['convsource'].pack(side=Tk.TOP)

        separator = Tk.Frame(self.frames['CLOpt'], height=4,
                             bd=5, relief=Tk.SUNKEN)
        separator.pack(fill=Tk.X, padx=10, pady=20, side=Tk.TOP)

        self.frames['Sensit'].pack(side=Tk.TOP)
        self.buttons['Noise'].pack(side=Tk.TOP)

        self.canvas1.mpl_connect('pick_event', self._onPick)
        self.canvas1.mpl_connect('motion_notify_event', self._doMask)
        self.canvas1.mpl_connect('button_release_event', self._onRelease)
        self.canvas1.mpl_connect('button_press_event', self._onPress)
        self.canvas1.mpl_connect('key_press_event', self.parent._onKeyPress)

        # toolbar_frame = Tk.Frame(self.me)
        # toolbar = NavigationToolbar2TkAgg(self.canvas1, toolbar_frame)
        # toolbar_frame.pack(side=Tk.LEFT)

        self.pressed = -1
        self.xy0 = [0, 0]
        self.moved = False
        self.resadd = False
        self.dorestore = True
        self._makeMask()
        self._reCalib()

    def _ReNoise(self):
        try:
            sensit = float(self.entries['Sensit'].get())
        except:
            showinfo(
                'ERROR!', 'Please, check the content of Sensit!\nIt should be a number!')
            return

        if sensit < 0.0:
            showinfo('ERROR!', 'The sensitivity should be >= 0!')
            return

         # Get the number of baselines and the number of integration times:

        Nsamples = float(self.parent.Nbas * self.parent.nH)
        sensPerSamp = sensit * np.sqrt(Nsamples) / np.sqrt(2.)
        self.parent._setNoise(sensPerSamp)
        self._reset(donoise=False)

    def _doRestore(self):

        if self.dorestore:
            self.dorestore = False
            toadd = self.cleanmodd[self.Np4:self.parent.Npix - self.Np4,
                                   self.Np4:self.parent.Npix - self.Np4]

        else:
            self.dorestore = True
            if self.resadd:
                toadd = (
                    self.cleanmod + self.residuals)[self.Np4:self.parent.Npix - self.Np4,
                                                    self.Np4:self.parent.Npix - self.Np4]
            else:
                toadd = self.cleanmod[self.Np4:self.parent.Npix - self.Np4,
                                      self.Np4:self.parent.Npix - self.Np4]

        self.CLEANPlotPlot.set_array(toadd)
        self.CLEANPlotPlot.norm.vmin = np.min(toadd)
        self.CLEANPlotPlot.norm.vmax = np.max(toadd)
        self.canvas1.draw()
        del toadd

    def _doRescale(self):

            # if True:
        clarr = self.CLEANPlotPlot.get_array()
        self.CLEANPlotPlot.norm.vmin = np.min(clarr)
        self.CLEANPlotPlot.norm.vmax = np.max(clarr)
        self.CLEANPlotPlot.set_array(clarr)
        rsarr = self.ResidPlotPlot.get_array()
        self.ResidPlotPlot.norm.vmin = np.min(rsarr)
        self.ResidPlotPlot.norm.vmax = np.max(rsarr)
        self.ResidPlotPlot.set_array(rsarr)

        del clarr, rsarr

        # self.CLEANPlot.autoscale() #.norm.vmin = np.min(clarr)
        self.canvas1.draw()

    def _ApplyGain(self):

        try:
            an1 = int(self.entries['Ant1'].curselection()[0])
        except:
            showinfo('WARNING!', 'No antenna selected!')
            return

        try:
            an2 = int(self.entries['Ant2'].curselection()[0])
        except:
            an2 = an1

        if an2 == an1:
            an2 = -1

        G = float((self.entries['Amp'].get()) / 100. *
                  np.exp(1.j * float(self.entries['Phas'].get()) * np.pi / 180.))
        H0 = int(self.entries['H0'].get())
        H1 = int(self.entries['H1'].get())

        self.parent._setGains(an1, an2, H0, H1, G)
        self._reset()

    def _makeMask(self):
        self.mask = np.zeros(np.shape(self.parent.beam))
        self.bmask = np.zeros(np.shape(self.parent.beam)).astype(np.bool)

    def _onPick(self, event):
        RA = event.mouseevent.xdata
        Dec = event.mouseevent.ydata
        yi = np.floor((self.Xaxmax - RA) / (
            2. * self.Xaxmax) * self.parent.Npix)
        xi = np.floor((self.Xaxmax - Dec) / (
            2. * self.Xaxmax) * self.parent.Npix)
        Flux = self.residuals[xi, yi]
        self.pickcoords = [xi, yi, RA, Dec]
        self.ResidText.set_text(self.fmtD2 % (Flux, RA, Dec,
                                              self.PEAK, self.RMS))
        if self.dorestore:
            if self.resadd:
                Flux = self.cleanmod[xi, yi] + self.residuals[xi, yi]
            else:
                Flux = self.cleanmod[xi, yi]
        else:
            Flux = self.cleanmodd[xi, yi]

        self.CLEANText.set_text(self.fmtDC % (Flux, RA, Dec,
                                              self.CLEANPEAK,
                                              self.CLEANPEAK / self.RMS) +
                                '\n' + self.Beamtxt)

        self.canvas1.draw()
        # self.canvas2.draw()

    def _onPress(self, event):
        self.canvas1._tkcanvas.focus_set()
        if event.inaxes == self.ResidPlot:
            self.pressed = int(event.button)
            RA = event.xdata
            Dec = event.ydata
            self.xydata = [RA, Dec]
            self.xy0[1] = np.floor(
                (self.Xaxmax - RA) / (2. * self.Xaxmax) * self.parent.Npix)
            self.xy0[0] = np.floor(
                (self.Xaxmax - Dec) / (2. * self.Xaxmax) * self.parent.Npix)
            self.moved = False

    def _onRelease(self, event):
        if event.inaxes != self.ResidPlot:
            self.moved = False

        if self.moved:
            RA = event.xdata
            Dec = event.ydata
            y1 = np.floor(
                (self.Xaxmax - RA) / (2. * self.Xaxmax) * self.parent.Npix)
            x1 = np.floor(
                (self.Xaxmax - Dec) / (2. * self.Xaxmax) * self.parent.Npix)
            xi, xf = [min(self.xy0[0], x1), max(self.xy0[0], x1)]
            yi, yf = [min(self.xy0[1], y1), max(self.xy0[1], y1)]
            if self.pressed == 1:
                self.mask[xi:xf, yi:yf] = 1.0
                self.bmask[xi:xf, yi:yf] = True
            else:
                self.mask[xi:xf, yi:yf] = 0.0
                self.bmask[xi:xf, yi:yf] = False

            for coll in self.MaskPlot.collections:
                self.ResidPlot.collections.remove(coll)

            self.MaskPlot = self.ResidPlot.contour(
                np.linspace(self.parent.Xaxmax / 2.,
                            -self.parent.Xaxmax /
                            2.,
                            self.parent.Npix / 2),
                                                   np.linspace(
                                                       self.parent.Xaxmax / 2.,
                                                               -self.parent.Xaxmax /
                                                                   2.,
                                                               self.parent.Npix / 2),
                                                   self.mask[
                                                       self.Np4:self.parent.Npix - self.Np4,
                                                             self.Np4:self.parent.Npix - self.Np4],
                                                   levels=[0.5])

            # self.ResidPlot.set_xlim((self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.))
            # self.ResidPlot.set_ylim((-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))
            self.CLEANPlot.set_xlim((self.parent.curzoom[1][0],
                                     self.parent.curzoom[1][1]))
            self.CLEANPlot.set_ylim((self.parent.curzoom[1][2],
                                     self.parent.curzoom[1][3]))
            self.canvas1.draw()

            self.Box.set_data([0., 0., 0., 0., 0.],
                              [0., 0., 0., 0., 0.])

        self.moved = False
        self.pressed = -1
        self.canvas1.draw()

    def _doMask(self, event):
        if self.pressed >= 0 and event.inaxes == self.ResidPlot:
            self.moved = True
            RA = event.xdata
            Dec = event.ydata
            self.Box.set_data([self.xydata[0], self.xydata[0], RA, RA,
                               self.xydata[0]], [self.xydata[1], Dec, Dec,
                                                 self.xydata[1],
                                                 self.xydata[1]])
            self.canvas1.draw()

    def _AddRes(self):
        if not self.dorestore:
            showinfo('ERROR',
                     'Cannot add residual to the (unrestored) CLEAN model!\nRestore first!')

        if self.resadd:
            self.resadd = False
            toadd = self.cleanmod[self.Np4:self.parent.Npix - self.Np4,
                                  self.Np4:self.parent.Npix - self.Np4]
        else:
            self.resadd = True
            toadd = (
                self.cleanmod + self.residuals)[self.Np4:self.parent.Npix - self.Np4,
                                                self.Np4:self.parent.Npix - self.Np4]

        self.CLEANPlotPlot.set_array(toadd)
        self.CLEANPlotPlot.norm.vmin = np.min(toadd)
        self.CLEANPlotPlot.norm.vmax = np.max(toadd)

        self.canvas1.draw()
        # self.canvas2.draw()
        del toadd

    def _reCalib(self):
        self.entries['Ant1'].delete(0, Tk.END)
        self.entries['Ant2'].delete(0, Tk.END)
        self.entries['H0'].set(0)
        self.entries['H1'].set(self.parent.nH)
        self.entries['Amp'].set(100)
        self.entries['Phas'].set(0)

        for i in range(self.parent.Nant):
            self.entries['Ant1'].insert(Tk.END, str(i + 1))
            self.entries['Ant2'].insert(Tk.END, str(i + 1))

        self.parent._setGains(-1, -1, 0, 0, 1.0)
        self._reset()

    def _reset(self, donoise=False):

        self.ResidPlot.cla()
        self.dorestore = True

        self.fmtD2 = r'% .2e Jy/beam at point' "\n" r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f ' "\n" r'Peak: % 4.2f Jy/beam ; rms: % 4.2f Jy/beam'
        self.fmtDC = r'Model: % .2e Jy/beam at point' "\n" r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f ' "\n" r'Peak: % 4.2f Jy/beam ; Dyn. Range: % 4.2f'

        dslice = self.parent.dirtymap[self.Np4:self.parent.Npix - self.Np4,
                                      self.Np4:self.parent.Npix - self.Np4]
        self.ResidPlotPlot = self.ResidPlot.imshow(dslice,
                                                   interpolation='nearest',
                                                   picker=True,
                                                   cmap=self.parent.currcmap)
        modflux = self.parent.dirtymap[self.parent.Nphf, self.parent.Nphf]
        self.RMS = np.sqrt(np.var(dslice) + np.average(dslice)**2.)
        self.PEAK = np.max(dslice)
        self.CLEANPEAK = 0.0
        self.pickcoords = [self.parent.Nphf, self.parent.Nphf, 0., 0.]
        self.ResidText = self.ResidPlot.text(0.05, 0.87,
                                             self.fmtD2 % (modflux, 0.0, 0.0,
                                                           self.PEAK, self.RMS),
                                             transform=self.ResidPlot.transAxes,
                                             bbox=dict(facecolor='white',
                                                       alpha=0.7))
        pl.setp(self.ResidPlotPlot, extent=(self.parent.Xaxmax / 2.,
                                            -self.parent.Xaxmax / 2.,
                                            -self.parent.Xaxmax / 2.,
                                            self.parent.Xaxmax / 2.))

        self.Xaxmax = float(self.parent.Xaxmax)

        self.Box = self.ResidPlot.plot([0., 0., 0., 0., 0.],
                                       [0., 0., 0., 0., 0.],
                                       lw=2, color='w')[0]

        self.ResidPlot.set_ylabel('Dec offset (as)')
        self.ResidPlot.set_xlabel('RA offset (as)')
        self.ResidPlot.set_title('RESIDUALS')

        self.MaskPlot = self.ResidPlot.contour(
            np.linspace(self.parent.Xaxmax / 2.,
                        -self.parent.Xaxmax /
                        2.,
                        self.parent.Npix / 2),
                                               np.linspace(
                                                   self.parent.Xaxmax / 2.,
                                                           -self.parent.Xaxmax /
                                                               2.,
                                                           self.parent.Npix / 2),
                                               self.mask[
                                                   self.Np4:self.parent.Npix - self.Np4,
                                                         self.Np4:self.parent.Npix - self.Np4],
                                               levels=[0.5])
        # pl.setp(self.MaskPlot,
        # extent=(self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))

        # self.ResidPlot.set_xlim((self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.))
        # self.ResidPlot.set_ylim((-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))

        self.residuals = np.copy(self.parent.dirtymap)
        self.cleanmod = np.zeros(np.shape(self.parent.dirtymap))
        self.cleanmodd = np.zeros(np.shape(self.parent.dirtymap))

        self.CLEANPlot.cla()
        self.CLEANPlotPlot = self.CLEANPlot.imshow(
            self.parent.dirtymap[self.Np4:self.parent.Npix - self.Np4,
                                 self.Np4:self.parent.Npix - self.Np4],
                                                   interpolation='nearest',
                                                   picker=True,
                                                   cmap=self.parent.currcmap)
        modflux = self.parent.dirtymap[self.parent.Nphf,
                                       self.parent.Nphf]
        self.CLEANText = self.CLEANPlot.text(0.05, 0.83,
                                             self.fmtDC % (
                                                 0.0, 0.0, 0.0, 0., 0.),
                                             transform=self.CLEANPlot.transAxes,
                                             bbox=dict(facecolor='white',
                                                       alpha=0.7))
        pl.setp(self.CLEANPlotPlot, extent=(self.parent.Xaxmax / 2.,
                                            -self.parent.Xaxmax / 2.,
                                            -self.parent.Xaxmax / 2.,
                                            self.parent.Xaxmax / 2.))
        self.CLEANPlot.set_ylabel('Dec offset (as)')
        self.CLEANPlot.set_xlabel('RA offset (as)')
        self.CLEANPlot.set_title('CLEAN (0 ITER)')
        self.CLEANPlotPlot.set_array(
            self.cleanmod[self.Np4:self.parent.Npix - self.Np4,
                          self.Np4:self.parent.Npix - self.Np4])

        # self.CLEANPlot.set_xlim((self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.))
        # self.CLEANPlot.set_ylim((-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))
        self.CLEANPlot.set_xlim((self.parent.curzoom[1][0],
                                 self.parent.curzoom[1][1]))
        self.CLEANPlot.set_ylim((self.parent.curzoom[1][2],
                                 self.parent.curzoom[1][3]))

        self.totiter = 0

        # DERIVE THE CLEAN BEAM
        MainLobe = np.where(self.parent.beam > 0.6)
        self.cleanBeam = np.zeros(np.shape(self.residuals))

        if len(MainLobe[0]) < 5:
            showinfo('ERROR!',
                     'The main lobe of the PSF is too narrow!\n CLEAN model will not be restored')
            self.cleanBeam[:] = 0.0
            self.cleanBeam[self.parent.Npix / 2, self.parent.Npix / 2] = 1.0
        else:
            dX = MainLobe[0] - self.parent.Npix / 2
            dY = MainLobe[1] - self.parent.Npix / 2

            try:
                fit = spfit.leastsq(lambda x:
                                    np.exp(-(dX * dX * x[0] + dY * dY * x[1] + dX * dY * x[2])) - self.parent.beam[
                                    MainLobe],
                                    [1., 1., 0.])
                Pang = 180. / np.pi * (np.arctan2(fit[0][2],
                                                  (fit[0][0] - fit[0][1])) / 2.)
                AmB = fit[0][2] / np.sin(2. * np.pi / 180. * Pang)
                ApB = fit[0][0] + fit[0][1]
                A = 2.355 * \
                    (2. / (ApB + AmB))**0.5 * \
                    self.parent.imsize / self.parent.Npix
                B = 2.355 * \
                    (2. / (ApB - AmB))**0.5 * \
                    self.parent.imsize / self.parent.Npix
                if A < B:
                    A, B = B, A
                    Pang = Pang - 90.
                if Pang < -90.:
                    Pang += 180.
                if Pang > 90.:
                    Pang -= 180.

                if B > 0.1:
                    self.Beamtxt = '%.1f x %.1f as (PA = %.1f deg.)' % (
                        A, B, Pang)
                else:
                    self.Beamtxt = '%.1f x %.1f mas (PA = %.1f deg.)' % (1000. * A,
                                                                         1000. *
                                                                         B,
                                                                         Pang)

                self.CLEANText.set_text(self.fmtDC % (0., 0., 0., 0., 0.) +
                                        '\n' + self.Beamtxt)
                # print 'BEAM FIT: ',fit[0], A, B, Pang
                ddX = np.outer(np.ones(self.parent.Npix),
                               np.arange(-self.parent.Npix / 2,
                                         self.parent.Npix / 2).astype(np.float64))
                ddY = np.outer(np.arange(-self.parent.Npix / 2,
                                         self.parent.Npix / 2).astype(
                                             np.float64),
                               np.ones(self.parent.Npix))

                self.cleanBeam[:] = np.exp(-(ddY * ddY * fit[0][0] +
                                             ddX * ddX * fit[0][1] +
                                             ddY * ddX * fit[0][2]))

                del ddX, ddY

            except:
                showinfo('ERROR!',
                         'Problems fitting the PSF main lobe!\n CLEAN model will not be restored')
                self.cleanBeam[:] = 0.0
                self.cleanBeam[
                    self.parent.Npix / 2, self.parent.Npix / 2] = 1.0

        self.resadd = False
        self.dorestore = True
        self.ffti = False

        self.totalClean = 0.0

        if donoise:
            self._ReNoise()
            # self.canvas1.mpl_connect('key_press_event',
            # self.parent._onKeyPress)
        self.canvas1.draw()
        # self.canvas2.draw()

        del modflux

    def _CLEAN(self):
        if np.sum(self.bmask) == 0:
            goods = np.ones(np.shape(self.bmask)).astype(np.bool)
            tempres = self.residuals
        else:
            goods = self.bmask
            tempres = self.residuals * self.mask

        psf = self.parent.beam

        try:
            gain = float(self.entries['Gain'].get())
            niter = int(self.entries['Niter'].get())
            thrs = float(self.entries['Thres'].get())
        except:
            showinfo('ERROR!',
                     'Please, check the content of Gain, # Iter, and Thres!\nShould be numbers!')
            return

        for i in range(niter):
            self.totiter += 1

            if thrs != 0.0:
                tempres[tempres < thrs] = 0.0
                if thrs < 0.0:
                    tempres = np.abs(tempres)

                if np.sum(tempres) == 0.0:
                    showinfo('INFO', 'Threshold reached in CLEAN masks!')
                    break

            rslice = self.residuals[self.Np4:self.parent.Npix - self.Np4,
                                    self.Np4:self.parent.Npix - self.Np4]
            peakpos = np.unravel_index(np.argmax(tempres),
                                       np.shape(self.residuals))
            peakval = self.residuals[peakpos[0], peakpos[1]]
            self.residuals -= gain * peakval * np.roll(np.roll(psf,
                                                               peakpos[
                                                               0] - self.parent.Npix / 2,
                                                               axis=0),
                                                       peakpos[
                                                       1] - self.parent.Npix / 2,
                                                       axis=1)
            tempres[goods] = self.residuals[goods]
            # MODIFY CLEAN MODEL!!
            self.cleanmodd[peakpos[0], peakpos[1]] += gain * peakval
            self.cleanmod += gain * peakval * np.roll(np.roll(self.cleanBeam,
                                                              peakpos[
                                                              0] - self.parent.Npix / 2,
                                                              axis=0),
                                                      peakpos[
                                                      1] - self.parent.Npix / 2,
                                                      axis=1)
            self.ResidPlotPlot.set_array(rslice)

            self.CLEANPEAK = np.max(self.cleanmod)
            self.totalClean += gain * peakval
            self.CLEANPlot.set_title(
                'CLEAN (%i ITER): %.2e Jy' % (self.totiter,
                                              self.totalClean))

            xi, yi, RA, Dec = self.pickcoords

            if self.dorestore:
                if self.resadd:
                    toadd = (self.cleanmod + self.residuals)
                else:
                    toadd = self.cleanmod
            else:
                toadd = self.cleanmodd

            clFlux = toadd[xi, yi]

            self.CLEANPlotPlot.set_array(
                toadd[self.Np4:self.parent.Npix - self.Np4,
                      self.Np4:self.parent.Npix - self.Np4])
            self.CLEANPlotPlot.norm.vmin = np.min(
                toadd[self.Np4:self.parent.Npix - self.Np4,
                      self.Np4:self.parent.Npix - self.Np4])
            self.CLEANPlotPlot.norm.vmax = np.max(
                toadd[self.Np4:self.parent.Npix - self.Np4,
                      self.Np4:self.parent.Npix - self.Np4])

            self.RMS = np.sqrt(np.var(rslice) + np.average(rslice)**2.)
            self.PEAK = np.max(rslice)
            # self.RMS = np.std(self.residuals)
            self.ResidText.set_text(
                self.fmtD2 % (self.residuals[xi, yi], RA, Dec,
                              self.PEAK, self.RMS))
            self.CLEANText.set_text(
                self.fmtDC % (clFlux, RA, Dec, self.CLEANPEAK,
                              self.CLEANPEAK / self.RMS) + '\n' + self.Beamtxt)

            self.canvas1.draw()

        # Re-draw if threshold reached:
        self.canvas1.draw()
        del tempres, psf, goods
        try:
            del toadd
        except:
            pass

    def _getHelp(self):
        win = Tk.Toplevel(self.me)
        win.title("Help")
        helptext = ScrolledText(win)
        helptext.config(state=Tk.NORMAL)
        helptext.insert('1.0', __CLEAN_help_text__)
        helptext.config(state=Tk.DISABLED)

        helptext.pack()
        Tk.Button(win, text='OK', command=win.destroy).pack()

    def _showFFT(self):
        from uvplotter1 import ApSynSim_UV_Plotter1

        try:
            self.parent.myUVPLOT.destroy()
        except:
            self.parent.myUVPLOT = ApSynSim_UV_Plotter1(self.parent)

#    try:
#      self.FFTwin.destroy()
#    except:
#      pass

    def _convSource(self):

        try:
            self.convSource.destroy()
        except:
            pass

        self.convSource = Tk.Toplevel(self.me)
        self.convSource.title("True source image")

        self.figCS1 = pl.figure(figsize=(6, 6))

        self.CS1 = self.figCS1.add_subplot(111, aspect='equal')
        # pl.axes([0.55,0.43,0.5,0.5],aspect='equal')

        self.figCS1.subplots_adjust(left=0.05, right=0.98)

        self.CSText = self.CS1.text(0.05, 0.87,
                                    self.parent.fmtD % (0.0, 0., 0.),
                                    transform=self.CS1.transAxes,
                                    bbox=dict(facecolor='white',
                                              alpha=0.7))

        self.frames = {}
        self.frames['FigCS'] = Tk.Frame(self.convSource)
        self.frames['CSFr'] = Tk.Frame(self.convSource)

        self.canvasCS1 = FigureCanvasTkAgg(self.figCS1,
                                           master=self.frames['FigCS'])

        self.canvasCS1.show()
        self.frames['FigCS'].pack(side=Tk.TOP)
        self.frames['CSFr'].pack(side=Tk.TOP)

        self.buttons['reloadCS'] = Tk.Button(self.frames['CSFr'],
                                             text="Reload",
                                             command=self._CSRead)
        self.buttons['reloadCS'].pack()

        self.canvasCS1.mpl_connect('pick_event',
                                   self._onCSPick)
        self.canvasCS1.get_tk_widget().pack(side=Tk.LEFT)
        # , fill=Tk.BOTH, expand=1)
        toolbar_frame = Tk.Frame(self.convSource)
        _toolbar = NavigationToolbar2TkAgg(self.canvasCS1, toolbar_frame)
        toolbar_frame.pack(side=Tk.LEFT)

        self._CSRead()

    def _CSRead(self):
        self.CSImage = (np.fft.fftshift(np.fft.ifft2(
                                        self.parent.modelfft * np.fft.fft2(
                                        np.fft.fftshift(self.cleanBeam))))).real

        self.CS1Plot = self.CS1.imshow(
            self.CSImage[self.Np4:self.parent.Npix - self.Np4,
                         self.Np4:self.parent.Npix - self.Np4],
                                       interpolation='nearest', picker=True,
                                       cmap=self.parent.currcmap)
        modflux = self.parent.dirtymap[self.parent.Nphf, self.parent.Nphf]
        self.CSText.set_text(self.parent.fmtD % (modflux, 0.0, 0.0))
     # transform=self.CS1.transAxes,bbox=dict(facecolor='white', alpha=0.7))
        pl.setp(self.CS1Plot, extent=(self.Xaxmax / 2.,
                                      -self.Xaxmax / 2.,
                                      -self.Xaxmax / 2.,
                                      self.Xaxmax / 2.))

        self.CS1.set_xlabel('RA offset (as)')
        self.CS1.set_ylabel('Dec offset (as)')

        self.CS1.set_title('TRUE SOURCE - CONVOLVED')

        _to_plot = np.abs(np.fft.fftshift(np.fft.fft2(
                                          np.fft.fftshift(self.residuals))))

        self.canvasCS1.draw()

    def _onCSPick(self, event):
        RA = event.mouseevent.xdata
        Dec = event.mouseevent.ydata
        yi = np.floor((self.Xaxmax - RA) / (
            2. * self.Xaxmax) * self.parent.Npix)
        xi = np.floor((self.Xaxmax - Dec) / (
            2. * self.Xaxmax) * self.parent.Npix)
        Flux = self.CSImage[xi, yi]
        self.CSText.set_text(self.parent.fmtD % (Flux, RA, Dec))

        self.canvasCS1.draw()
