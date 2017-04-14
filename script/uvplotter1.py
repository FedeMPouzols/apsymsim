#############################################################################
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
#############################################################################

import numpy as np
import pylab as pl

try:
    import Tkinter as Tk
except:
    import tkinter as Tk

from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)


class ApSynSim_UV_Plotter1(object):

  def quit(self):
    self.FFTwin.destroy()

  def __init__(self, parent):

    self.parent = parent

    self.FFTwin = Tk.Toplevel(self.parent.tks)
    self.FFTwin.title("UV space")

    menubar = Tk.Menu(self.FFTwin)
    menubar.add_command(label="Quit", command=self.quit)
    self.FFTwin.config(menu=menubar)

    self.figUV1 = pl.figure(figsize=(8.5, 7))

    self.UVPSF = self.figUV1.add_subplot(231, aspect='equal')
    pl.setp(self.UVPSF.get_xticklabels(), visible=False)
    self.UVCLEANMOD = self.figUV1.add_subplot(
      232, sharex=self.UVPSF, sharey=self.UVPSF, aspect='equal')
    pl.setp(self.UVCLEANMOD.get_xticklabels(), visible=False)
    pl.setp(self.UVCLEANMOD.get_yticklabels(), visible=False)

    self.UVResid = self.figUV1.add_subplot(
      234, sharex=self.UVPSF, sharey=self.UVPSF, aspect='equal')
    self.UVCLEAN = self.figUV1.add_subplot(
      235, sharex=self.UVPSF, sharey=self.UVPSF, aspect='equal')
    pl.setp(self.UVCLEAN.get_yticklabels(), visible=False)

    self.UVSOURCE = self.figUV1.add_subplot(
      233, sharex=self.UVPSF, sharey=self.UVPSF, aspect='equal')
    pl.setp(self.UVSOURCE.get_xticklabels(), visible=False)
    pl.setp(self.UVSOURCE.get_yticklabels(), visible=False)

    self.UVSOURCECONV = self.figUV1.add_subplot(
      236, sharex=self.UVPSF, sharey=self.UVPSF, aspect='equal')
    pl.setp(self.UVSOURCECONV.get_yticklabels(), visible=False)

    self.figUV1.subplots_adjust(left=0.1, right=0.98,
                                top=0.90, bottom=0.15,
                                wspace=0.02, hspace=0.15)

    self.UVfmt = '%.2e Jy'
    self.PSFfmt = '%.2e'

    self.PSFText = self.UVPSF.text(0.05, 0.87, self.PSFfmt % (0.0),
                                   transform=self.UVPSF.transAxes,
                                   bbox=dict(facecolor='white',
                                             alpha=0.7))

    self.ResidText = self.UVResid.text(0.05, 0.87, self.UVfmt % (0.0),
                                       transform=self.UVResid.transAxes,
                                       bbox=dict(facecolor='white',
                                                 alpha=0.7))

    self.CLEANText = self.UVCLEAN.text(0.05, 0.87, self.UVfmt % (0.0),
                                       transform=self.UVCLEAN.transAxes,
                                       bbox=dict(facecolor='white',
                                                 alpha=0.7))

    self.CLEANMODText = self.UVCLEANMOD.text(0.05, 0.87, self.UVfmt % (0.0),
                                             transform=self.UVCLEANMOD.transAxes,
                                             bbox=dict(facecolor='white',
                                                       alpha=0.7))

    self.UVSOURCEText = self.UVSOURCE.text(0.05, 0.87, self.UVfmt % (0.0),
                                           transform=self.UVSOURCE.transAxes,
                                           bbox=dict(facecolor='white',
                                                     alpha=0.7))

    self.UVSOURCECONVText = self.UVSOURCECONV.text(0.05, 0.87, self.UVfmt % (0.0),
                                                   transform=self.UVSOURCECONV.transAxes,
                                                   bbox=dict(facecolor='white',
                                                             alpha=0.7))

    self.frames = {}
    self.frames['FigUV'] = Tk.Frame(self.FFTwin)
    self.frames['BFr'] = Tk.Frame(self.FFTwin)

    self.canvasUV1 = FigureCanvasTkAgg(self.figUV1, master=self.frames['FigUV'])

    self.canvasUV1.show()
    self.frames['FigUV'].pack(side=Tk.TOP)
    self.frames['BFr'].pack(side=Tk.TOP)

    self.buttons = {}
    self.buttons['reload'] = Tk.Button(self.frames['BFr'], text="Reload",
                                       command=self._FFTRead)
    self.buttons['reload'].pack()

    self.canvasUV1.mpl_connect('pick_event', self._onUVPick)
    self.canvasUV1.get_tk_widget().pack(side=Tk.LEFT)  # , fill=Tk.BOTH, expand=1)
    toolbar_frame = Tk.Frame(self.FFTwin)
    __toolbar = NavigationToolbar2TkAgg(self.canvasUV1, toolbar_frame)
    toolbar_frame.pack(side=Tk.LEFT)

    self._FFTRead()

  def _FFTRead(self):

    Toplot = np.abs(np.fft.fftshift(
      np.fft.fft2(np.fft.fftshift(self.parent.beam))))

    vmax = np.max(Toplot)
    __vmin = np.min(Toplot)
    self.UVPSFPlot = self.UVPSF.imshow(Toplot, vmin=0.0, vmax=vmax,
                                       cmap=self.parent.currcmap,
                                       picker=True, interpolation='nearest')
    pl.setp(self.UVPSFPlot,
            extent=(-self.parent.UVmax + self.parent.UVSh,
                    self.parent.UVmax + self.parent.UVSh,
                    -self.parent.UVmax - self.parent.UVSh,
                    self.parent.UVmax - self.parent.UVSh))
    self.UVPSF.set_ylabel(self.parent.vlab)

    self.UVPSF.set_title('UV - PSF')

    try:
      Toplot = np.abs(np.fft.fftshift(
        np.fft.fft2(np.fft.fftshift(self.parent.myCLEAN.residuals))))
    except:
      Toplot = np.zeros(np.shape(self.parent.beam))

    vmax = np.max(np.abs(np.fft.fft2(self.parent.dirtymap)))
    self.UVResidPlot = self.UVResid.imshow(Toplot, vmin=0.0, vmax=vmax,
                                           cmap=self.parent.
                                           currcmap, picker=True,
                                           interpolation='nearest')
    pl.setp(self.UVResidPlot,
            extent=(-self.parent.UVmax + self.parent.UVSh,
                    self.parent.UVmax + self.parent.UVSh,
                    -self.parent.UVmax - self.parent.UVSh,
                    self.parent.UVmax - self.parent.UVSh))
    self.UVResid.set_xlabel(self.parent.ulab)
    self.UVResid.set_ylabel(self.parent.vlab)

    self.UVResid.set_title('UV - RESIDUALS (REST.)')

    try:
      Toplot = np.abs(np.fft.fftshift(
        np.fft.fft2(np.fft.fftshift(self.parent.myCLEAN.cleanmod))))
    except:
      Toplot = np.zeros(np.shape(self.parent.beam))

    vmax = np.max(Toplot)
    self.UVCLEANPlot = self.UVCLEAN.imshow(Toplot, vmin=0.0, vmax=vmax,
                                           cmap=self.parent.currcmap,
                                           picker=True,
                                           interpolation='nearest')
    pl.setp(self.UVCLEANPlot,
            extent=(-self.parent.UVmax + self.parent.UVSh,
                    self.parent.UVmax + self.parent.UVSh,
                    -self.parent.UVmax - self.parent.UVSh,
                    self.parent.UVmax - self.parent.UVSh))
    self.UVCLEAN.set_xlabel(self.parent.ulab)

    self.UVCLEAN.set_title('UV - CLEAN (REST.)')

    try:
      Toplot = np.abs(np.fft.fftshift(
        np.fft.fft2(np.fft.fftshift(self.parent.myCLEAN.cleanmodd))))
    except:
      Toplot = np.zeros(np.shape(self.parent.beam))

    vmax = np.max(Toplot)
    self.UVCLEANMODPlot = self.UVCLEANMOD.imshow(Toplot, vmin=0.0, vmax=vmax,
                                                 cmap=self.parent.currcmap,
                                                 picker=True,
                                                 interpolation='nearest')
    pl.setp(self.UVCLEANMODPlot,
            extent=(-self.parent.UVmax + self.parent.UVSh,
                    self.parent.UVmax + self.parent.UVSh,
                    -self.parent.UVmax - self.parent.UVSh,
                    self.parent.UVmax - self.parent.UVSh))

    self.UVCLEANMOD.set_title('UV - CLEAN (MODEL)')

    try:
      Toplot = np.fft.fftshift(np.abs(
        self.parent.modelfft*np.fft.fft2(self.parent.myCLEAN.cleanBeam)))
    except:
      Toplot = np.zeros(np.shape(self.parent.beam))

    vmax = np.max(Toplot)
    self.UVSOURCECONVPlot = self.UVSOURCECONV.imshow(Toplot, vmin=0.0, vmax=vmax,
                                                     cmap=self.parent.currcmap,
                                                     picker=True,
                                                     interpolation='nearest')
    pl.setp(self.UVSOURCECONVPlot,
            extent=(-self.parent.UVmax+self.parent.UVSh,
                    self.parent.UVmax+self.parent.UVSh,
                    -self.parent.UVmax-self.parent.UVSh,
                    self.parent.UVmax-self.parent.UVSh))

    self.UVSOURCECONV.set_title('UV - SOURCE (REST.)')

    Toplot = np.fft.fftshift(np.abs(self.parent.modelfft))

    vmax = np.max(Toplot)
    self.UVSOURCEPlot = self.UVSOURCE.imshow(Toplot, vmin=0.0, vmax=vmax,
                                             cmap=self.parent.currcmap,
                                             picker=True,
                                             interpolation='nearest')
    pl.setp(self.UVSOURCEPlot,
            extent=(-self.parent.UVmax + self.parent.UVSh,
                    self.parent.UVmax + self.parent.UVSh,
                    -self.parent.UVmax - self.parent.UVSh,
                    self.parent.UVmax - self.parent.UVSh))

    self.UVSOURCE.set_title('UV - SOURCE')

    self.canvasUV1.draw()

  def _onUVPick(self, event):

    Up = event.mouseevent.xdata - self.parent.UVSh
    Vp = event.mouseevent.ydata + self.parent.UVSh

    yi = np.floor((self.parent.UVmax+Up)/(self.parent.UVmax)*self.parent.Npix/2.)
    xi = np.floor((self.parent.UVmax-Vp)/(self.parent.UVmax)*self.parent.Npix/2.)

    if xi > 0 and yi > 0 and xi < self.parent.Npix and yi < self.parent.Npix:
      self.PSFText.set_text(self.PSFfmt % self.UVPSFPlot.get_array()[xi, yi])
      self.ResidText.set_text(self.UVfmt % self.UVResidPlot.get_array()[xi, yi])
      self.CLEANText.set_text(self.UVfmt % self.UVCLEANPlot.get_array()[xi, yi])
      self.CLEANMODText.set_text(self.UVfmt % self.UVCLEANMODPlot.get_array()[xi, yi])
      self.UVSOURCEText.set_text(self.UVfmt % self.UVSOURCEPlot.get_array()[xi, yi])
      self.UVSOURCECONVText.set_text(self.UVfmt % self.UVSOURCECONVPlot.get_array()[xi, yi])
    else:
      self.PSFText.set_text(self.PSFfmt % 0.0)
      self.ResidText.set_text(self.UVfmt % 0.0)
      self.CLEANText.set_text(self.UVfmt % 0.0)
      self.CLEANMODText.set_text(self.UVfmt % 0.0)
      self.UVSOURCEText.set_text(self.UVfmt % 0.0)
      self.UVSOURCECONVText.set_text(self.UVfmt % 0.0)

    self.canvasUV1.draw()
