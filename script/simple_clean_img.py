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
import matplotlib as mpl
import scipy.optimize as spfit

from tkMessageBox import showinfo


class SimpleCleanImg(object):

    """ Simplified clean image. Does the clean-image calculations and shows
    results onto one of the subaxes of the 'parent'. Assumes too much about
    the 'parent'. Requires among other things parent.cleanPlot."""

    DEFAULT_NITER = 400

    def __init__(self, parent, plot_axis):
        self.parent = parent
        self.plot_axis = plot_axis
        self.frames = {}
        self.cleanmod = None
        self.cleanmodd = None
        self.cleanBeam = None
        self.Beamtxt = ''
        self._recalib()

    def __del__(self):
        self.residuals[:] = 0.0
        self.cleanmod[:] = 0.0
        self.cleanmodd[:] = 0.0
        self.cleanBeam[:] = 0.0

    def _make_mask(self):
        self.mask = np.zeros(np.shape(self.parent.beam))
        self.bmask = np.zeros(np.shape(self.parent.beam)).astype(np.bool)

    def _recalib(self):
        self.entries = {}
        self.entries['Gain'] = 0.1
        self.entries['Niter'] = self.DEFAULT_NITER
        self.entries['Thres'] = 0.0

        self.entries['Ant1'] = 0
        self.entries['Ant2'] = 0
        self.entries['H0'] = 0
        self.entries['H1'] = self.parent.nH
        self.entries['Amp'] = 100
        self.entries['Phas'] = 0

        self.parent._setGains(-1, -1, 0, 0, 1.0)
        self._reset()

    def _init_values(self):
        self.totiter = 0

        dslice = self.parent.dirtymap[self.Np4:self.parent.Npix - self.Np4,
                                      self.Np4:self.parent.Npix - self.Np4]
        self.RMS = np.sqrt(np.var(dslice) + np.average(dslice)**2.)
        self.PEAK = np.max(dslice)
        self.CLEANPEAK = 0.0
        self.pickcoords = [self.parent.Nphf, self.parent.Nphf, 0., 0.]

        self.Xaxmax = float(self.parent.Xaxmax)
        self.residuals = np.copy(self.parent.dirtymap)
        self.cleanmod = np.zeros(np.shape(self.parent.dirtymap))
        self.cleanmodd = np.zeros(np.shape(self.parent.dirtymap))

    def _reset(self, donoise=False):

        self.Np4 = self.parent.Npix / 4

        self.dorestore = True

        self.fmtD2 = (r'% .2e Jy/beam at point' "\n"
                      r'$\Delta\alpha = $ % 4.2f / '
                      '$\Delta\delta = $ % 4.2f ' "\n"
                      r'Peak: % 4.2f Jy/beam ; rms: % 4.2f Jy/beam')
        self.fmtDC = (r'Model: % .2e Jy/beam at point' "\n"
                      r'$\Delta\alpha = $ % 4.2f / '
                      '$\Delta\delta = $ % 4.2f ' "\n"
                      r'Peak: % 4.2f Jy/beam ; Dyn. Range: % 4.2f')

        self._init_values()

        self.parent.cleanPlot.cla()
        self.CLEANPlotPlot = self.parent.cleanPlot.imshow(
            self.parent.dirtymap[self.Np4:self.parent.Npix - self.Np4,
                                 self.Np4:self.parent.Npix - self.Np4],
            interpolation='nearest',
            picker=True,
            cmap=self.parent.currcmap)

        self.CLEANText = self.parent.cleanPlot.text(
            0.05, 0.78, self.fmtDC % (0.0, 0.0, 0.0, 0., 0.),
            transform=self.parent.cleanPlot.transAxes,
            bbox=dict(facecolor='white', alpha=0.7),
            fontsize=10, verticalalignment='bottom')

        pl.setp(self.CLEANPlotPlot, extent=(self.parent.Xaxmax / 2.,
                                            -self.parent.Xaxmax / 2.,
                                            -self.parent.Xaxmax / 2.,
                                            self.parent.Xaxmax / 2.))
        self.parent.cleanPlot.set_ylabel('Dec offset (as)')
        self.parent.cleanPlot.set_xlabel('RA offset (as)')
        # more dynamic than self.parent.cleanPlot.set_title('CLEAN (0 ITER)')
        self.CLEANTitle = self.parent.cleanPlot.text(
            0.5, 1.02, 'CLEAN ( 0 ITER)',
            transform=self.parent.cleanPlot.transAxes,
            backgroundcolor=self.parent.figUV.get_facecolor(),
            verticalalignment='bottom',
            horizontalalignment='center',
            fontsize=14)

        self.CLEANPlotPlot.set_array(
            self.cleanmod[self.Np4:self.parent.Npix - self.Np4,
                          self.Np4:self.parent.Npix - self.Np4])

        self.parent.cleanPlot.set_xlim((self.parent.curzoom[1][0],
                                        self.parent.curzoom[1][1]))
        self.parent.cleanPlot.set_ylim((self.parent.curzoom[1][2],
                                        self.parent.curzoom[1][3]))

        # DERIVE THE CLEAN BEAM
        MainLobe = np.where(self.parent.beam > 0.6)
        self.cleanBeam = np.zeros(np.shape(self.residuals))

        if len(MainLobe[0]) < 5:
            showinfo('ERROR!',
                     'The main lobe of the PSF is too narrow!\n'
                     'CLEAN model will not be restored')
            self.cleanBeam[:] = 0.0
            self.cleanBeam[self.parent.Npix / 2, self.parent.Npix / 2] = 1.0
        else:
            dX = MainLobe[0] - self.parent.Npix / 2
            dY = MainLobe[1] - self.parent.Npix / 2

            try:
                fit = spfit.leastsq(lambda x:
                                    np.exp(-(dX * dX * x[0] +
                                             dY * dY * x[1] + dX * dY * x[2])) -
                                    self.parent.beam[MainLobe],
                                    [1., 1., 0.])
                Pang = 180. / np.pi * (np.arctan2(fit[0][2],
                                                  (fit[0][0] - fit[0][1])) / 2.)
                AmB = fit[0][2] / np.sin(2. * np.pi / 180. * Pang)
                ApB = fit[0][0] + fit[0][1]
                A = (2.355 * (2. / (ApB + AmB))**0.5 *
                     self.parent.imsize / self.parent.Npix)
                B = (2.355 * (2. / (ApB - AmB))**0.5 *
                     self.parent.imsize / self.parent.Npix)
                if A < B:
                    A, B = B, A
                    Pang = Pang - 90.
                if Pang < -90.:
                    Pang += 180.
                if Pang > 90.:
                    Pang -= 180.

                if B > 0.1:
                    self.Beamtxt = ('%.1f x %.1f as (PA = %.1f deg.)'
                                    % (A, B, Pang))
                else:
                    self.Beamtxt = ('%.1f x %.1f mas (PA = %.1f deg.)' %
                                    (1000. * A, 1000. * B, Pang))

                self.CLEANText.set_text(self.fmtDC % (0., 0., 0., 0., 0.) +
                                        '\n' + self.Beamtxt)
                # print('BEAM FIT: ',fit[0], A, B, Pang)
                ddX = np.outer(np.ones(self.parent.Npix),
                               np.arange(-self.parent.Npix / 2,
                                         self.parent.Npix / 2).
                               astype(np.float64))
                ddY = np.outer(np.arange(-self.parent.Npix / 2,
                                         self.parent.Npix / 2).
                               astype(np.float64),
                               np.ones(self.parent.Npix))

                self.cleanBeam[:] = np.exp(-(ddY * ddY * fit[0][0] +
                                             ddX * ddX * fit[0][1] +
                                             ddY * ddX * fit[0][2]))

                del ddX, ddY

            except Exception as exc:
                showinfo('ERROR!',
                         'Problems fitting the PSF main lobe!\n'
                         'CLEAN model will not be restored. Details: {}'.
                         format(exc))
                self.cleanBeam[:] = 0.0
                self.cleanBeam[
                    self.parent.Npix / 2, self.parent.Npix / 2] = 1.0

        self.resadd = False
        self.dorestore = True
        self.ffti = False

        self.totalClean = 0.0

        if donoise:
            self._ReNoise()

        pl.draw()

    def do_clean(self, pre_iter=None):
        self._make_mask()

        if np.sum(self.bmask) == 0:
            goods = np.ones(np.shape(self.bmask)).astype(np.bool)
            tempres = self.residuals
        else:
            goods = self.bmask
            tempres = self.residuals * self.mask

        psf = self.parent.beam

        try:
            gain = float(self.entries['Gain'])
            if pre_iter is not None:
                niter = pre_iter
            else:
                niter = int(self.entries['Niter'])
            thrs = float(self.entries['Thres'])
        except:
            showinfo('ERROR!',
                     'Please, check the content of Gain, '
                     '# Iter, and Thres!\nShould be numbers!')
            return

        self.parent.figUV.canvas.draw()

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
            self.residuals -= gain * peakval * np.roll(
                np.roll(psf, peakpos[0] - self.parent.Npix / 2, axis=0),
                peakpos[1] - self.parent.Npix / 2, axis=1)
            tempres[goods] = self.residuals[goods]
            # MODIFY CLEAN MODEL!!
            self.cleanmodd[peakpos[0], peakpos[1]] += gain * peakval
            self.cleanmod += gain * peakval * np.roll(
                np.roll(self.cleanBeam, peakpos[0] - self.parent.Npix / 2,
                        axis=0),
                peakpos[1] - self.parent.Npix / 2, axis=1)
            self.CLEANPEAK = np.max(self.cleanmod)
            self.totalClean += gain * peakval
            # self.parent.cleanPlot.set_title('CLEAN (%i ITER): %.2e Jy' %
            #                        (self.totiter, self.totalClean))
            self.CLEANTitle.set_text('CLEAN (%i ITER): %.2e Jy' %
                                     (self.totiter, self.totalClean))

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

            self.CLEANText.set_text(self.fmtDC %
                                    (clFlux, RA, Dec, self.CLEANPEAK,
                                     self.CLEANPEAK / self.RMS) +
                                    '\n' + self.Beamtxt)

            # This is to avoid a full pl.draw() which is much slower as
            # it refreshes the whole pictures, all panels
            # Consider using matplotlib.animations in the future
            self.CLEANTitle.set_text('CLEAN (%i ITER): %.2e Jy' %
                                     (self.totiter, self.totalClean))

            self.parent.cleanPlot.draw_artist(self.CLEANPlotPlot)
            self.parent.cleanPlot.draw_artist(self.CLEANText)
            self.parent.cleanPlot.draw_artist(self.CLEANTitle)
            # blit does not update axes, titles, legends, etc.
            # Also, blit with cleanPlot.bbox would leave out the title
            # self.parent.figUV.canvas.blit(self.parent.cleanPlot.bbox)
            bbox_title = mpl.transforms.TransformedBbox(
                mpl.transforms.Bbox.from_extents([0, 0, 1, 1.1]),
                self.parent.cleanPlot.transAxes)
            self.parent.figUV.canvas.blit(bbox_title)

        # Re-draw all if threshold reached:
        pl.draw()

        del tempres, psf, goods
        try:
            del toadd
        except:
            pass
