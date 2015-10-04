#############################################################################
#
#    APSYNSIM: A real-time Aperture Synthesis Simulator
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

#import Tkinter
import FileDialog
import matplotlib as mpl
mpl.use('TkAgg')
import numpy as np
import pylab as pl
import scipy.ndimage.interpolation as spndint
import scipy.optimize as spfit
from matplotlib.widgets import Slider, Button
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
import matplotlib.image as plimg
from ScrolledText import ScrolledText


from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
try:
  import Tkinter as Tk
except:
  import tkinter as Tk


from matplotlib.backend_bases import NavigationToolbar2
import tkFileDialog
from tkMessageBox import showinfo
import os
import time
import sys

__version__ = '0.4-b1'



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
 
Pressing the "Reduce data" button will open a new window, where you can
CLEAN your dirty image and apply corrupting gains to your antennas (see 
help in that window for more details).

Enjoy!


"""


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
You can further click on CLEAN, to continue deconvolving.

Pressing RESET will undo all CLEANING. 

You also need to press RESET to refresh the images for CLEANing, if you 
change anything in the main program window (e.g., observing wavelength, 
antenna positions, etc.).

Pressing "+/- Resid" will add (or remove) the residuals from the CLEANed 
image.  By default, the residuals are NOT added (i.e., only the restored
CLEAN components are shown in the CLEAN image).


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


class Interferometer(object):

  def quit(self,event=None):
  
    self.tks.destroy()
    sys.exit()

  def __init__(self,antenna_file="",model_file="",tkroot=None):

    self.__version__ = __version__

  #  if tkroot is None:
   #   self.tks = None  #Tkinter.Tk()
  #  else:
    self.tks = tkroot
   # self.tks.withdraw()

    self.Hfac = np.pi/180.*15.
    self.deg2rad = np.pi/180.
    self.curzoom = [0,0,0,0]
    self.robust = 0.0
    self.deltaAng = 1.*self.deg2rad
    self.gamma = 0.5  # Gamma correction to plot model.
    self.lfac = 1.e6   # Lambda units (i.e., 1.e6 => Mlambda)
    self.ulab = r'U (M$\lambda$)'
    self.vlab = r'V (M$\lambda$)'

    self.currcmap = cm.jet

# Default of defaults!
    nH = 200
    Npix = 512   # Image pixel size. Must be a power of 2
    DefaultMod = 'Nebula.model'
    DefaultArray = 'Golay_12.array'

# Overwrite defaults from config file:
    d1 = os.path.dirname(os.path.realpath(__file__))
    print d1

#   execfile(os.path.join(os.path.basename(d1),'APSYNSIM.config'))
    try:
      conf = open(os.path.join(d1,'APSYNSIM.config'))
    except:
      d1 = os.getcwd()
      conf = open(os.path.join(d1,'APSYNSIM.config'))
      
    for line in conf.readlines():
      temp=line.replace(' ','')
      if len(temp)>2:
         if temp[0:4] == 'Npix':
           Npix = int(temp[5:temp.find('#')])
         if temp[0:2] == 'nH':
           nH = int(temp[3:temp.find('#')])
         if temp[0:10] == 'DefaultMod':
           DefaultModel = temp[12:temp.find('#')].replace('\'','').replace('\"','')
         if temp[0:12] == 'DefaultArray':
           DefaultArray = temp[14:temp.find('#')].replace('\'','').replace('\"','')

    conf.close()

# Set instance configuration values:
    self.nH = nH
    self.Npix = Npix

    self.datadir = os.path.join(d1,'..','PICTURES')
    self.arraydir  = os.path.join(d1,'..','ARRAYS')
    self.modeldir  = os.path.join(d1,'..','SOURCE_MODELS')

 # Try to read a default initial array:
    if len(antenna_file)==0:
      try:
        antenna_file = os.path.join(self.arraydir,DefaultArray)
      except:
        pass

 # Try to read a default initial model:
    if len(model_file)==0:
      try:
        model_file = os.path.join(self.modeldir,DefaultModel)
      except:
        pass


    self.lock=False
    self._onSphere = False

    self.readAntennas(str(antenna_file))
    self.readModels(str(model_file))
    self.GUI() #makefigs=makefigs)


  def showError(self,message):
    showinfo('ERROR!', message)
    raise Exception(message)


  def _getHelp(self):
    win = Tk.Toplevel(self.tks)
    win.title("Help")
    helptext = ScrolledText(win)
    helptext.config(state=Tk.NORMAL)
   # about = __help_text__
    helptext.insert('1.0',__help_text__)
    helptext.config(state=Tk.DISABLED)

   # Tk.Label(win, text=about, width=100, height=10).pack()
    helptext.pack()
    Tk.Button(win, text='OK', command=win.destroy).pack()


  def GUI(self): # ,makefigs=True):

    mpl.rcParams['toolbar'] = 'None'

    self.Nphf = self.Npix/2
    self.robfac = 0.0
    self.figUV = pl.figure(figsize=(15,8))

    if self.tks is None:
       self.canvas = self.figUV.canvas
    else:
      self.canvas = FigureCanvasTkAgg(self.figUV, master=self.tks)
      self.canvas.show()
      menubar = Tk.Menu(self.tks)
      menubar.add_command(label="Help", command=self._getHelp)
      menubar.add_command(label="Quit", command=self.quit)

      self.tks.config(menu=menubar)
    #  menubar.pack(side=Tk.TOP)
      self.canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)




#    self.figUV.canvas.toolbar = CustomToolbar(self.figUV.canvas, self.figUV)
#    self.figUV.canvas.toolbar.update()

  #  self.figUV.canvas.set_window_title('Aperture Synthesis Simulator (I. Marti-Vidal, Onsala Space Observatory) - version  %s'%__version__)
    self.antPlot = self.figUV.add_subplot(231,aspect='equal')
    self.UVPlot = self.figUV.add_subplot(232,aspect='equal',axisbg=(0.4,0.4,0.4))
    self.beamPlot = self.figUV.add_subplot(233,aspect='equal')
    self.modelPlot = self.figUV.add_subplot(235,aspect='equal')
    self.dirtyPlot = self.figUV.add_subplot(236,aspect='equal')

    self.spherePlot = pl.axes([0.53,0.82,0.12,0.12],projection='3d',aspect='equal')

    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = 10 * np.outer(np.cos(u), np.sin(v))
    y = 10 * np.outer(np.sin(u), np.sin(v))
    z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))
    self.spherePlotPlot = self.spherePlot.plot_surface(x, y, z,  rstride=4, cstride=4, color=(0.8,0.8,1.0))
    self.spherePlot._axis3don = False
    self.spherePlot.patch.set_alpha(0.8)
    beta = np.zeros(100)
    self.arrayPath = [np.zeros(self.nH), np.zeros(self.nH), np.zeros(self.nH)]
    self.sphereArray = self.spherePlot.plot([],[],[],'y',linewidth=3)
    self.spherePlot.set_xlim3d((-6,6))
    self.spherePlot.set_ylim3d((-6,6))
    self.spherePlot.set_zlim3d((-6,6))
    self.spherePlot.patch.set_alpha(0.8)
    self.spherePlot.elev = 45.


    self.figUV.subplots_adjust(left=0.05,right=0.99,top=0.95,bottom=0.07,hspace=0.25)
    self.canvas.mpl_connect('pick_event', self._onPick)
    self.canvas.mpl_connect('motion_notify_event', self._onAntennaDrag)
    self.canvas.mpl_connect('button_release_event',self._onRelease)
    self.canvas.mpl_connect('button_press_event',self._onPress)
    self.canvas.mpl_connect('key_press_event', self._onKeyPress)
    self.pickAnt = False

    self.fmtH = r'$\phi = $ %3.1f$^\circ$   $\delta = $ %3.1f$^\circ$' "\n" r'H = %3.1fh / %3.1fh'
    self.fmtBas = r'Bas %i $-$ %i  at  H = %4.2fh'
    self.fmtVis = r'Amp: %.1e Jy.   Phase: %5.1f deg.' 
    self.fmtA = 'N = %i'
    self.fmtA2 = '  Picked Ant. #%i' 
    self.fmtA3 = '\n%6.1fm | %6.1fm'
    fmtB1 = r'$\lambda = $ %4.1fmm  '%(self.wavelength*1.e6)
    self.fmtB = fmtB1 + "\n" + r'% 4.2f Jy/beam' + "\n" + r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f '
    self.fmtD = r'% .2e Jy/beam' "\n" r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f '
    self.fmtM = r'%.2e Jy/pixel' "\n"  r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f'

    self.wax = {}
    self.widget = {}
    self.wax['lat'] = pl.axes([0.07,0.45,0.25,0.04])
    self.wax['dec'] = pl.axes([0.07,0.40,0.25,0.04])
    self.wax['H0'] = pl.axes([0.07,0.35,0.25,0.04])
    self.wax['H1'] = pl.axes([0.07,0.30,0.25,0.04])
    self.wax['wave'] = pl.axes([0.07,0.25,0.25,0.04])
    self.wax['robust'] = pl.axes([0.07,0.20,0.25,0.04])
    self.wax['add'] = pl.axes([0.07,0.14,0.08,0.05])
    self.wax['rem'] = pl.axes([0.155,0.14,0.08,0.05])
    self.wax['reduce'] = pl.axes([0.24,0.14,0.08,0.05])
    self.wax['save'] =  pl.axes([0.07,0.08,0.08,0.05]) 
    self.wax['loadarr']=pl.axes([0.155,0.08,0.08,0.05])
    self.wax['quit']=pl.axes([0.155,0.02,0.08,0.05])
    self.wax['loadmod']=pl.axes([0.24,0.08,0.08,0.05])
    self.wax['gammacorr']=pl.axes([0.46,0.08,0.13,0.02],axisbg='white')
   # self.wax['gammacorr'].yaxis.label.set_color('white')
   # self.wax['gammacorr'].label.set_color('white')

    self.widget['robust'] = Slider(self.wax['robust'],r'Robust',-2.,2.,valinit=0.0)
    self.widget['lat'] = Slider(self.wax['lat'],r'Lat (deg)',-90.,90.,valinit=self.lat/self.deg2rad)
    self.widget['dec'] = Slider(self.wax['dec'],r'Dec (deg)',-90.,90.,valinit=self.dec/self.deg2rad)
    self.widget['H0'] = Slider(self.wax['H0'],r'H$_{0}$ (h)',-12.,12.,valinit=self.Hcov[0]/self.Hfac)
    self.widget['H1'] = Slider(self.wax['H1'],r'H$_{1}$ (h)',-12.,12.,valinit=self.Hcov[1]/self.Hfac)
    self.widget['wave'] = Slider(self.wax['wave'],r'$\lambda$ (mm)',self.wavelength*1.e6/10.,10.*self.wavelength*1.e6,valinit=self.wavelength*1.e6)
    self.widget['add'] = Button(self.wax['add'],r'+ Antenna')
    self.widget['rem'] = Button(self.wax['rem'],r'$-$ Antenna')
    self.widget['reduce'] = Button(self.wax['reduce'],r'Reduce data')
    self.widget['save'] = Button(self.wax['save'],'Save array')
    self.widget['loadarr'] = Button(self.wax['loadarr'],'Load array')
    self.widget['loadmod'] = Button(self.wax['loadmod'],'Load model')
    self.widget['quit'] = Button(self.wax['quit'],'Quit')
    self.widget['gammacorr'] = Slider(self.wax['gammacorr'],'gamma',0.1,1.0,valinit=self.gamma,color='red')
    self.widget['gammacorr'].label.set_color('white')
    self.widget['gammacorr'].valtext.set_color('white')

    self.widget['robust'].on_changed(self._onRobust)
    self.widget['lat'].on_changed(self._onKeyLat)
    self.widget['dec'].on_changed(self._onKeyDec)
    self.widget['H0'].on_changed(self._onKeyH0)
    self.widget['H1'].on_changed(self._onKeyH1)
    self.widget['wave'].on_changed(self._changeWavelength)
    self.widget['add'].on_clicked(self._addAntenna)
    self.widget['rem'].on_clicked(self._removeAntenna)
    self.widget['save'].on_clicked(self.saveArray)
    self.widget['loadarr'].on_clicked(self.loadArray)
    self.widget['loadmod'].on_clicked(self.loadModel)
    self.widget['gammacorr'].on_changed(self._gammacorr)
    self.widget['quit'].on_clicked(self.quit)
    self.widget['reduce'].on_clicked(self.reduce)


   # if makefigs:
   #   self._makefigs()
 # def _makefigs(self):

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


#    def new_home(self, *args, **kwargs):
#       print 'home pressed'
#    pl.get_current_fig_manager().toolbar.home = new_home #.bind('<Button-1>', new_home)


 #   pl.show()
    self.canvas.draw()
 #   Tk.mainloop()



  def reduce(self,event):

    if self.tks is not None:
      self.myCLEAN = CLEANer(self)
   

  def readAntennas(self,antenna_file):
    
    self.Hcov = [-12.0*self.Hfac,12.0*self.Hfac]
    self.Hmax = np.pi
    self.lat = 45.*self.deg2rad
    self.dec = 60.*self.deg2rad
    self.trlat = [np.sin(self.lat),np.cos(self.lat)]
    self.trdec = [np.sin(self.dec), np.cos(self.dec)]
    self.Xmax = 4.0

    if len(antenna_file)==0:
      self.Nant = 7
      self.antPos=[[0.0,0.0],[0.0,1.],[0.0,2.0],[1.,-1.],[2.0,-2.0],[-1.,-1.],[-2.0,-2.0]]
    

    if len(antenna_file)>0: 
     if not os.path.exists(antenna_file):
      self.showError("\n\nAntenna file %s does not exist!\n\n"%antenna_file)
      return False
    
     else:
      antPos=[]
      Hcov = [0,0]
      Nant = 0
      Xmax = 0.0
      fi = open(antenna_file)
      for li,l in enumerate(fi.readlines()):
        comm = l.find('#')
        if comm>=0:
          l = l[:comm]     
        it = l.split()
        if len(it)>0:
          if it[0]=='ANTENNA':   
            antPos.append(map(float,it[1:]))
            Nant += 1
            antPos[-1][0] *= 1.e-3 ; antPos[-1][1] *= 1.e-3
            Xmax = np.max(np.abs(antPos[-1]+[Xmax]))
          elif it[0]=='LATITUDE':
            lat = float(it[1])*self.deg2rad
            trlat = [np.sin(lat),np.cos(lat)]
          elif it[0]=='DECLINATION':
            dec = float(it[1])*self.deg2rad
            trdec = [np.sin(dec),np.cos(dec)]
          elif it[0]=='HOUR_ANGLE':
            Hcov[0] = float(it[1])*self.Hfac
            Hcov[1] = float(it[2])*self.Hfac
          else:
            self.showError("\n\nWRONG SYNTAX IN LINE %i:\n\n %s...\n\n"%(li+1,l[:max(10,len(l))]))


      if np.abs(lat-dec>=np.pi/2.):
         self.showError("\n\nSource is either not observable or just at the horizon!\n\n")
         return False
      if Nant<2:
         self.showError("\n\nThere should be at least 2 antennas!\n\n")
         return False



      self.Nant = Nant
      self.antPos = antPos
      self.lat = lat
      self.dec = dec
      self.trlat = trlat
      self.trdec = trdec
      self.Hcov = Hcov
      self.Xmax = Xmax


      cosW = -np.tan(self.lat)*np.tan(self.dec)
      if np.abs(cosW) < 1.0:
        Hhor = np.arccos(cosW)
      elif np.abs(self.lat-self.dec)>np.pi/2.:
        Hhor = 0
      else:
        Hhor = np.pi

      if Hhor>0.0:
        if self.Hcov[0]< -Hhor:
          self.Hcov[0] = -Hhor
        if self.Hcov[1]>  Hhor:
          self.Hcov[1] =  Hhor

      self.Hmax = Hhor
      H = np.linspace(self.Hcov[0],self.Hcov[1],self.nH)[np.newaxis,:]
      self.Xmax = Xmax*1.5
      fi.close()
   
    return True




  def readModels(self,model_file):

    self.wavelength = 1.e-6  # in km.
    self.imsize = 4. 
    self.imfiles = []

    if len(model_file)==0:
      self.models = [['G',0.,0.4,1.0,0.1],['D',0.,0.,2.,0.5],['P',-0.4,-0.5,0.1]]
      self.Xaxmax = self.imsize/2.
      return True

    if len(model_file)>0: 
     if not os.path.exists(model_file):
      self.showError("\n\nModel file %s does not exist!\n\n"%model_file)
      return False

     else:
      fixsize = False
      models = []
      imfiles = []
      Xmax = 0.0
      fi = open(model_file)
      for li,l in enumerate(fi.readlines()):
        comm = l.find('#')
        if comm>=0:
          l = l[:comm]     
        it = l.split()
        if len(it)>0:
          if it[0]=='IMAGE':
            imfiles.append([str(it[1]),float(it[2])])
          elif it[0] in ['G','D','P']:   
            models.append([it[0]]+map(float,it[1:]))
            if models[-1][0] != 'P':
              models[-1][4] = np.abs(models[-1][4])
              Xmax = np.max([np.abs(models[-1][1])+models[-1][4],
                 np.abs(models[-1][2])+models[-1][4],Xmax])
          elif it[0] == 'WAVELENGTH':
            wavelength = float(it[1])*1.e-3
          elif it[0] == 'IMSIZE':
            imsize = 2.*float(it[1])
            fixsize = True
          else:
            self.showError("\n\nWRONG SYNTAX IN LINE %i:\n\n %s...\n\n"%(li+1,l[:max(10,len(l))]))

      if len(models)+len(imfiles)==0:
         self.showError("\n\nThere should be at least 1 model component!\n\n")

      self.models=models
      self.wavelength=wavelength
      self.imsize=imsize
      self.imfiles = imfiles

      if not fixsize:
        self.imsize = Xmax*1.1


      self.Xaxmax = self.imsize/2.

      fi.close()



    return True





  def _changeWavelength(self,wave):

     self.wavelength = wave*1.e-6
     fmtB1 = r'$\lambda = $ %4.1fmm  '%(self.wavelength*1.e6)
     self.fmtB = fmtB1 + "\n" r'% 4.2f Jy/beam' "\n" r'$\Delta\alpha = $ % 4.2f / $\Delta\delta = $ % 4.2f '

     self._plotAntennas(redo=False)
     self._changeCoordinates(rescale=True)
     self._plotModelFFT(redo=False) 



  def _changeCoordinates(self,rescale=False):

    if self.lat>np.pi/2.:
      self.lat = np.pi/2.
      return
    elif self.lat<-np.pi/2.:
      self.lat = -np.pi/2.
      return

    if self.dec>np.pi/2.:
      self.dec = np.pi/2.
      return
    elif self.dec<-np.pi/2.:
      self.dec = -np.pi/2.
      return

    self.trlat = [np.sin(self.lat),np.cos(self.lat)]
    self.trdec = [np.sin(self.dec), np.cos(self.dec)]


    cosW = -np.tan(self.lat)*np.tan(self.dec)
    if np.abs(cosW) < 1.0:
      Hhor = np.arccos(cosW)
    elif np.abs(self.lat-self.dec)>np.pi/2.:
      Hhor = 0
    else:
      Hhor = np.pi

    self.Hmax = Hhor

    if self.Hmax>0.0:
      if self.Hcov[0]< -self.Hmax:
        self.Hcov[0] = -self.Hmax
        self.lock=True
        self.widget['H0'].set_val(self.Hcov[0]/self.Hfac)
        self.lock=False
      if self.Hcov[1]>  self.Hmax:
        self.Hcov[1] =  self.Hmax
        self.lock=True
        self.widget['H1'].set_val(self.Hcov[1]/self.Hfac)
        self.lock=False


    newtext = self.fmtH%(self.lat/self.deg2rad,self.dec/self.deg2rad,self.Hcov[0]/self.Hfac,self.Hcov[1]/self.Hfac)
    self.latText.set_text(newtext)
    self.Horig =  np.linspace(self.Hcov[0],self.Hcov[1],self.nH)
    H = self.Horig[np.newaxis,:]
    self.H = [np.sin(H),np.cos(H)] 
    self._setBaselines()
    self._setBeam()
    self._plotBeam(redo=False)
    self._plotAntennas(redo=False,rescale=True)
    self._plotDirty(redo=False)

    pl.draw()
    self.canvas.draw()



  def _setGains(self,An1,An2,H0,H1,G):

    self.Gains[:] = 1.0
 #   print 'FROM: ',An1,An2,H0,H1,G,np.sum(self.Gains)

    for nb in range(self.Nbas):
  #    print self.antnum[nb],An1,An2
      if An1 == self.antnum[nb][0]:
        if An2 == -1 or An2 == self.antnum[nb][1]:
          self.Gains[nb,H0:H1] *= G
      if An1 == self.antnum[nb][1]:
        if An2 == -1 or An2 == self.antnum[nb][0]:
         self.Gains[nb,H0:H1] *= np.conjugate(G)


    self._setBaselines()
  #  self._plotAntennas(redo=False)
    self._setBeam()
    self._plotBeam(redo=False)
    self._plotDirty(redo=False)
   #    pl.draw()
    self.canvas.draw()

   # self.myCLEAN._reset()

 #   print 'TO: ',np.sum(self.Gains)



  def _prepareBeam(self):
    
    self.beam = np.zeros((self.Npix,self.Npix),dtype=np.float32)
    self.totsampling = np.zeros((self.Npix,self.Npix),dtype=np.float32)
    self.dirtymap = np.zeros((self.Npix,self.Npix),dtype=np.float32)
#    self.totsamplingU = np.zeros((self.Npix,self.Npix),dtype=np.float32)
    self.robustsamp = np.zeros((self.Npix,self.Npix),dtype=np.float32)
    self.Gsampling = np.zeros((self.Npix,self.Npix),dtype=np.complex64)
    self.Grobustsamp = np.zeros((self.Npix,self.Npix),dtype=np.complex64)


  def _prepareBaselines(self):

    self.Nbas = self.Nant*(self.Nant-1)/2
    NBmax = self.Nbas
    self.B = np.zeros((NBmax,self.nH),dtype=np.float32)
    self.basnum = np.zeros((self.Nant,self.Nant-1),dtype=np.int8)
  #  self.refant = np.zeros((self.Nant,self.Nant-1),dtype=np.int8)
    self.basidx = np.zeros((self.Nant,self.Nant),dtype=np.int8)
    self.antnum = np.zeros((NBmax,2),dtype=np.int8)

    self.Gains = np.ones((self.Nbas,self.nH),dtype=np.complex64)

    self.Horig = np.linspace(self.Hcov[0],self.Hcov[1],self.nH)
    H = self.Horig[np.newaxis,:]
    self.H = [np.sin(H),np.cos(H)]

    bi = 0
    nii = [0 for n in range(self.Nant)]
    for n1 in range(self.Nant-1):
      for n2 in range(n1+1,self.Nant):
        self.basnum[n1,nii[n1]] = bi
        self.basnum[n2,nii[n2]] = bi
        self.basidx[n1,n2] = bi
        self.antnum[bi] = [n1,n2]
   #     self.refant[n1,nii[n1]] = 0
   #     self.refant[n2,nii[n2]] = 1
        nii[n1] += 1; nii[n2] += 1
        bi += 1

    self.u = np.zeros((NBmax,self.nH))
    self.v = np.zeros((NBmax,self.nH))
    self.ravelDims = (NBmax,self.nH)

#    self.sampling = np.zeros((self.Nant,self.Npix,self.Npix),dtype=np.float32)



  def _setBaselines(self,antidx=-1):

   if antidx==-1:
     bas2change = range(self.Nbas)
   else:
     bas2change = self.basnum[antidx].flatten()

   for currBas in bas2change:
     n1,n2 = self.antnum[currBas]
     self.B[currBas,0] = -(self.antPos[n2][1]-self.antPos[n1][1])*self.trlat[0]/self.wavelength
     self.B[currBas,1] = (self.antPos[n2][0]-self.antPos[n1][0])/self.wavelength
     self.B[currBas,2] = (self.antPos[n2][1]-self.antPos[n1][1])*self.trlat[1]/self.wavelength
     self.u[currBas,:] = -(self.B[currBas,0]*self.H[0] + self.B[currBas,1]*self.H[1])
     self.v[currBas,:] = -self.B[currBas,0]*self.trdec[0]*self.H[1]+self.B[currBas,1]*self.trdec[0]*self.H[0]+self.trdec[1]*self.B[currBas,2]




  def _gridUV(self,antidx=-1):

 #  antidx=-1

   if antidx==-1:
     bas2change = range(self.Nbas)
 #    dummy = np.zeros(self.Npix*self.Npix,dtype=np.int8)
     self.pixpos = [[] for nb in bas2change]
     self.totsampling[:] = 0.0
     self.Gsampling[:] = 0.0
   else:
     bas2change = map(int,list(self.basnum[antidx].flatten()))

  # Nphf = int(self.Npix/2)
   self.UVpixsize = 2./(self.imsize*np.pi/180./3600.)

  # print bas2change

   for nb in bas2change:
     pixU = np.rint(self.u[nb]/self.UVpixsize).flatten().astype(np.int32)
     pixV = np.rint(self.v[nb]/self.UVpixsize).flatten().astype(np.int32)
     goodpix = np.logical_and(np.abs(pixU)<self.Nphf,np.abs(pixV)<self.Nphf)
     pU = pixU[goodpix] + self.Nphf
     pV = pixV[goodpix] + self.Nphf
     mU = -pixU[goodpix] + self.Nphf
     mV = -pixV[goodpix] + self.Nphf
     if not antidx==-1:
     # for pi in range(len(self.pixU[nb])):
     #  self.totsampling[self.pixV[nb][pi]+self.Nphf,-self.pixU[nb][pi]+self.Nphf] -= 1.0
     #  self.totsampling[-self.pixV[nb][pi]+self.Nphf,self.pixU[nb][pi]+self.Nphf] -= 1.0
     #  self.Gsampling[self.pixV[nb][pi]+self.Nphf,-self.pixU[nb][pi]+self.Nphf] -= self.Gains[nb,pi]
     #  self.Gsampling[-self.pixV[nb][pi]+self.Nphf,self.pixU[nb][pi]+self.Nphf] -= np.conjugate(self.Gains[nb,pi])
   #    tempG = float(self.Gsampling[pV[2],mU[2]])
       self.totsampling[self.pixpos[nb][1],self.pixpos[nb][2]] -= 1.0
       self.totsampling[self.pixpos[nb][3],self.pixpos[nb][0]] -= 1.0
       self.Gsampling[self.pixpos[nb][1],self.pixpos[nb][2]] -= self.Gains[nb,goodpix]
       self.Gsampling[self.pixpos[nb][3],self.pixpos[nb][0]] -= np.conjugate(self.Gains[nb,goodpix])
   #    print tempG, self.totsampling[pV[2],mU[2]]

    # self.pixU[nb] = np.copy(pixU[goodpix])
    # self.pixV[nb] = np.copy(pixV[goodpix])
 #    for pi in range(len(self.pixU[nb])):
 #      self.totsampling[self.pixV[nb][pi]+self.Nphf,-self.pixU[nb][pi]+self.Nphf] += 1.0
 #      self.totsampling[-self.pixV[nb][pi]+self.Nphf,self.pixU[nb][pi]+self.Nphf] += 1.0
 #      self.Gsampling[self.pixV[nb][pi]+self.Nphf,-self.pixU[nb][pi]+self.Nphf] += self.Gains[nb,pi]
 #      self.Gsampling[-self.pixV[nb][pi]+self.Nphf,self.pixU[nb][pi]+self.Nphf] += np.conjugate(self.Gains[nb,pi])
     self.pixpos[nb] = [np.copy(pU),np.copy(pV),np.copy(mU),np.copy(mV)]

     self.totsampling[pV,mU] += 1.0
     self.totsampling[mV,pU] += 1.0
     self.Gsampling[pV,mU] += self.Gains[nb,goodpix]
     self.Gsampling[mV,pU] += np.conjugate(self.Gains[nb,goodpix])



  # self.totsamplingU[:] = 0.0
  # self.GsamplingU[:] = 0.0
  # isdata = self.totsampling>0.0
  # self.totsamplingU[isdata] = 1.0
  # self.GsamplingU[isdata] /= self.totsampling[isdata]

   self.robfac = (5.*10.**(-self.robust))**2.*(2.*self.Nbas*self.nH)/np.sum(self.totsampling**2.)


  def _setBeam(self,antidx=-1):

  # print 'A:',self.robfac
   self._gridUV(antidx=antidx) 
  # print 'B:',self.robfac

   self.robustsamp[:] = self.totsampling/(1.+self.robfac*self.totsampling)
   self.Grobustsamp[:] = self.Gsampling/(1.+self.robfac*self.totsampling)


   self.beam[:] = np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(self.robustsamp))).real 
  # Nphf = int(self.Npix/2)
   self.beamScale = np.max(self.beam[self.Nphf:self.Nphf+1,self.Nphf:self.Nphf+1])
   self.beam[:] /= self.beamScale

  # print 'C:',np.std(self.beam)


  def _prepareModel(self):

    pixsize = float(self.imsize)/self.Npix
    xx = np.linspace(-self.imsize/2.,self.imsize/2.,self.Npix)
    yy = np.ones(self.Npix,dtype=np.float32)
    distmat = np.zeros((self.Npix,self.Npix),dtype=np.float32)
    self.modelim = np.zeros((self.Npix,self.Npix),dtype=np.float32)
   # Nphf = int(self.Npix/2)
    for model in self.models:
      xsh = -model[1]
      ysh = -model[2]
      xpix = np.rint(xsh/pixsize).astype(np.int32)
      ypix = np.rint(ysh/pixsize).astype(np.int32)
      centy = np.roll(xx,ypix)
      centx = np.roll(xx,xpix)
      distmat[:] = np.outer(centy**2.,yy) + np.outer(yy,centx**2.)
      if model[0]=='D':
        mask = np.logical_or(distmat<=model[4]**2.,distmat==np.min(distmat))
        self.modelim[mask] += float(model[3])/np.sum(mask)
      elif model[0]=='G':
        gauss = np.exp(-distmat/(2.*model[4]**2.))
        self.modelim[:] += float(model[3])*gauss/np.sum(gauss)
      elif model[0]=='P':
        if np.abs(xpix+self.Nphf)<self.Npix and np.abs(ypix+self.Nphf)<self.Npix:
          yint = ypix+self.Nphf
          xint = xpix+self.Nphf
       #   print ypix+Nphf, xpix+Nphf
          self.modelim[yint,xint] += float(model[3])

    for imfile in self.imfiles:
      if not os.path.exists(imfile[0]):
       # showinfo('File missing!','File %s does NOT exist. \nWill try to look in the %s directory'%(imfile[0],self.datadir)) 
        imfile[0] = os.path.join(self.datadir,imfile[0])
        if not os.path.exists(imfile[0]):
          self.showError('File %s does NOT exist. Cannot read the model!'%imfile[0]) 
          return

      Np4 = self.Npix/4
      img = plimg.imread(imfile[0]).astype(np.float32)
      dims = np.shape(img)
      d3 = min(2,dims[2])
      d1 = float(max(dims))
      avimg = np.average(img[:,:,:d3],axis=2)
      avimg -= np.min(avimg)
      avimg *= imfile[1]/np.max(avimg)
      if d1 == self.Nphf:
        sh0 = (self.Nphf-dims[0])/2
        sh1 = (self.Nphf-dims[1])/2
        self.modelim[sh0+Np4:sh0+Np4+dims[0], sh1+Np4:sh1+Np4+dims[1]] += zoomimg
      else:
        zoomimg = spndint.zoom(avimg,float(self.Nphf)/d1)
        zdims = np.shape(zoomimg)
        zd0 = min(zdims[0],self.Nphf)
        zd1 = min(zdims[1],self.Nphf)
    #    print zd0, zd1, zdims, dims
        sh0 = (self.Nphf-zdims[0])/2
        sh1 = (self.Nphf-zdims[1])/2
        self.modelim[sh0+Np4:sh0+Np4+zd0, sh1+Np4:sh1+Np4+zd1] += zoomimg[:zd0,:zd1]

    #  figg = pl.figure()
    #  pl.imshow(self.modelim[Np4:self.Npix-Np4,Np4:self.Npix-Np4])
    #  pl.show()

    self.modelim[self.modelim<0.0] = 0.0
    self.modelfft = np.fft.fft2(np.fft.fftshift(self.modelim))


  def _plotModel(self,redo=True):
    
    Np4 = self.Npix/4

    if redo:
      self.modelPlot.cla()
      self.modelPlotPlot = self.modelPlot.imshow(np.power(self.modelim[Np4:self.Npix-Np4,Np4:self.Npix-Np4],self.gamma),picker=True,interpolation='nearest',vmin=0.0,vmax=np.max(self.modelim)**self.gamma,cmap=self.currcmap)
 #     self.modelPlotPlot = self.modelPlot.imshow(self.modelim[Np4:self.Npix-Np4,Np4:self.Npix-Np4],picker=True,interpolation='nearest') #,vmin=0.0,vmax=np.max(self.modelim))

      modflux = self.modelim[self.Nphf,self.Nphf]
      self.modelText = self.modelPlot.text(0.05,0.87,self.fmtM%(modflux,0.0,0.0),
         transform=self.modelPlot.transAxes,bbox=dict(facecolor='white', 
         alpha=0.7))
      pl.setp(self.modelPlotPlot, extent=(self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.))
      self.modelPlot.set_ylabel('Dec offset (as)')
      self.modelPlot.set_xlabel('RA offset (as)')
      self.modelPlot.set_title('MODEL IMAGE')
      self._plotAntennas(redo=False)
    else:
      self.modelPlotPlot.set_data(np.power(self.modelim[Np4:self.Npix-Np4,Np4:self.Npix-Np4],self.gamma))
  #    self.modelPlotPlot.set_data(self.modelim[Np4:self.Npix-Np4,Np4:self.Npix-Np4])
      extr = [0.0,np.max(self.modelim)**self.gamma]
      self.modelPlotPlot.norm.vmin = extr[0]
      self.modelPlotPlot.norm.vmax = extr[1]
      pl.setp(self.modelPlotPlot, extent=(self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.))



  def _plotModelFFT(self,redo=True):

    self.UVmax = self.Npix/2./self.lfac*self.UVpixsize
    self.UVSh = -self.UVmax/self.Npix
    self.FFTtoplot = np.fft.fftshift(self.modelfft)
    toplot = np.abs(self.FFTtoplot)
    mval = np.min(toplot)
    Mval = np.max(toplot)
    dval = (Mval-mval)/2.
 #   print 'MVALS: ',mval,Mval

    if redo:
       mymap = pl.gray()
       self.UVPlotFFTPlot = self.UVPlot.imshow(toplot,cmap=mymap,vmin=mval-dval,vmax=Mval+dval,picker=5)
  #     pl.setp(self.UVPlotFFTPlot, extent=(-self.UVmax/2.,self.UVmax/2.,-self.UVmax/2.,self.UVmax/2.))
       pl.setp(self.UVPlotFFTPlot, extent=(-self.UVmax+self.UVSh,self.UVmax+self.UVSh,-self.UVmax-self.UVSh,self.UVmax-self.UVSh))
    else:
       self.UVPlotFFTPlot.set_data(toplot)
    #   pl.setp(self.UVPlotFFTPlot, extent=(-self.UVmax/2.,self.UVmax/2.,-self.UVmax/2.,self.UVmax/2.))
       self.UVPlotFFTPlot.norm.vmin = mval-dval
       self.UVPlotFFTPlot.norm.vmax = Mval+dval

    


  def _plotDirty(self,redo=True):
    Np4 = self.Npix/4

    self.dirtymap[:] = (np.fft.fftshift(np.fft.ifft2(self.modelfft*np.fft.ifftshift(self.Grobustsamp)))/self.beamScale).real
    extr = [np.min(self.dirtymap),np.max(self.dirtymap)]
   # Nphf = int(self.Npix/2)
    if redo:
      self.dirtyPlot.cla()
      self.dirtyPlotPlot = self.dirtyPlot.imshow(self.dirtymap[Np4:self.Npix-Np4,Np4:self.Npix-Np4],interpolation='nearest',picker=True, cmap=self.currcmap)
      modflux = self.dirtymap[self.Nphf,self.Nphf]
      self.dirtyText = self.dirtyPlot.text(0.05,0.87,self.fmtD%(modflux,0.0,0.0),
         transform=self.dirtyPlot.transAxes,bbox=dict(facecolor='white', 
         alpha=0.7))
      pl.setp(self.dirtyPlotPlot, extent=(self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.))
      self.curzoom[1] = (self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.)
      self.dirtyPlot.set_ylabel('Dec offset (as)')
      self.dirtyPlot.set_xlabel('RA offset (as)')
      self.dirtyPlot.set_title('DIRTY IMAGE')
    else:
      self.dirtyPlotPlot.set_data(self.dirtymap[Np4:self.Npix-Np4,Np4:self.Npix-Np4])
      self.dirtyPlotPlot.norm.vmin = extr[0]
      self.dirtyPlotPlot.norm.vmax = extr[1]
  #    pl.setp(self.dirtyPlotPlot, extent=(self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.))






  def _plotAntennas(self,redo=True,rescale=False):

   mw = 2.*self.Xmax/self.wavelength/self.lfac
   if mw < 0.1 and self.lfac == 1.e6:
      self.lfac = 1.e3
    #  self.Xmax *= 1.e-3
      self.ulab = r'U (k$\lambda$)'
      self.vlab = r'V (k$\lambda$)'
   elif mw >=100. and self.lfac == 1.e3:
      self.lfac = 1.e6
    #  self.Xmax *= 1.e3
      self.ulab = r'U (M$\lambda$)'
      self.vlab = r'V (M$\lambda$)'
   

   if redo:
    toplot = np.array(self.antPos[:self.Nant])
    self.antPlot.cla()
    self.antPlotPlot = self.antPlot.plot(toplot[:,0],toplot[:,1],'og', picker=5)[0]
    self.antPlot.set_xlim((-self.Xmax,self.Xmax))
    self.antPlot.set_ylim((-self.Xmax,self.Xmax))
    self.curzoom[3] = (-self.Xmax,self.Xmax,-self.Xmax,self.Xmax)
    self.antPlot.set_xlabel('E-W offset (km)')
    self.antPlot.set_ylabel('N-S offset (km)')
    self.antPlot.set_title('ARRAY CONFIGURATION')
    self.antText = self.antPlot.text(0.05,0.88,self.fmtA%self.Nant,transform=self.antPlot.transAxes)
    self.UVPlotPlot = []
    toplotu = self.u.flatten()/self.lfac ;  toplotv = self.v.flatten()/self.lfac ; 
    self.UVPlotPlot.append(self.UVPlot.plot(toplotu, toplotv,'.b',markersize=1,picker=2)[0])
    self.UVPlotPlot.append(self.UVPlot.plot(-toplotu,-toplotv,'.r',markersize=1,picker=2)[0])
    self.UVPlot.set_xlim((2.*self.Xmax/self.wavelength/self.lfac,-2.*self.Xmax/self.wavelength/self.lfac))
    self.UVPlot.set_ylim((2.*self.Xmax/self.wavelength/self.lfac,-2.*self.Xmax/self.wavelength/self.lfac))
    self.curzoom[2] = (2.*self.Xmax/self.lfac,-2.*self.Xmax/self.lfac,2.*self.Xmax/self.lfac,-2.*self.Xmax/self.lfac)
    self.latText = self.UVPlot.text(0.05,0.87,self.fmtH%(self.lat/self.deg2rad,self.dec/self.deg2rad,self.Hcov[0]/self.Hfac,self.Hcov[1]/self.Hfac),transform=self.UVPlot.transAxes)
    self.latText.set_color('orange')
    self.basText = self.UVPlot.text(0.05,0.02,self.fmtBas%(0,0,0.0),transform=self.UVPlot.transAxes)
    self.visText = self.UVPlot.text(0.05,0.08,self.fmtVis%(0.0,0.0),transform=self.UVPlot.transAxes)
    self.visText.set_color('orange')

    self.basText.set_color('orange')
    self.UVPlot.set_xlabel(self.ulab)
    self.UVPlot.set_ylabel(self.vlab) 
    self.UVPlot.set_title('UV PLANE')

    self.antLabelPlot = []

    for i in range(self.Nant):
      self.antLabelPlot.append(self.antPlot.annotate(str(i+1),textcoords = 'offset points',xy=(toplot[i,0],toplot[i,1]),xytext=(-7,4)))
  

   else:

    toplot = np.array(self.antPos[:self.Nant])
    self.antPlotPlot.set_data(toplot[:,0],toplot[:,1])
    toplotu = self.u.flatten()/self.lfac ;  toplotv = self.v.flatten()/self.lfac ; 
    if rescale:
      self.antPlot.set_xlim((-self.Xmax,self.Xmax))
      self.antPlot.set_ylim((-self.Xmax,self.Xmax))
      self.UVPlot.set_xlim((2.*self.Xmax/self.wavelength/self.lfac,-2.*self.Xmax/self.wavelength/self.lfac))
      self.UVPlot.set_ylim((2.*self.Xmax/self.wavelength/self.lfac,-2.*self.Xmax/self.wavelength/self.lfac))
      self.curzoom[2] = (2.*self.Xmax/self.lfac,-2.*self.Xmax/self.lfac,2.*self.Xmax/self.lfac,-2.*self.Xmax/self.lfac)
      self.curzoom[3] = (-self.Xmax,self.Xmax,-self.Xmax,self.Xmax)


    if len(self.antLabelPlot)>self.Nant:
      for i in range(self.Nant,len(self.antLabelPlot)):
         self.antLabelPlot[i].set_visible(False)


    for i in range(self.Nant):
      if i>len(self.antLabelPlot)-1:
        self.antLabelPlot.append(self.antPlot.annotate(str(i+1),textcoords = 'offset points',xy=(toplot[i,0],toplot[i,1]),xytext=(-7,4)))
      else:
        self.antLabelPlot[i].set_visible(True)
        self.antLabelPlot[i].xy = (toplot[i,0],toplot[i,1])
    

    self.UVPlotPlot[0].set_data(toplotu,toplotv)
    self.UVPlotPlot[1].set_data(-toplotu,-toplotv)
    self.UVPlot.set_xlabel(self.ulab)
    self.UVPlot.set_ylabel(self.vlab) 

    

  def _plotBeam(self,redo=True):

    Np4 = self.Npix/4
    if redo:
      self.beamPlot.cla()
      self.beamPlotPlot = self.beamPlot.imshow(self.beam[Np4:self.Npix-Np4,Np4:self.Npix-Np4],picker=True,interpolation='nearest', cmap=self.currcmap)
      self.beamText = self.beamPlot.text(0.05,0.80,self.fmtB%(1.0,0.0,0.0),
           transform=self.beamPlot.transAxes,bbox=dict(facecolor='white', alpha=0.7))
      self.beamPlot.set_ylabel('Dec offset (as)')
      self.beamPlot.set_xlabel('RA offset (as)')
      pl.setp(self.beamPlotPlot, extent=(self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.))
      self.curzoom[0] = (self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.)
      self.beamPlot.set_title('DIRTY BEAM')
      pl.draw()
      self.canvas.draw()
    else:
      self.beamPlotPlot.set_data(self.beam[Np4:self.Npix-Np4,Np4:self.Npix-Np4])
      self.beamText.set_text(self.fmtB%(1.0,0.0,0.0))
   #   pl.setp(self.beamPlotPlot, extent=(self.Xaxmax/2.,-self.Xaxmax/2.,-self.Xaxmax/2.,self.Xaxmax/2.))


    self.nptot = np.sum(self.totsampling[:])
    self.beamPlotPlot.norm.vmin = np.min(self.beam)
    self.beamPlotPlot.norm.vmax = 1.0

    if np.sum(self.totsampling[self.Nphf-4:self.Nphf+4,self.Nphf-4:self.Nphf+4])==self.nptot:
      warn = 'WARNING!\nToo short baselines for such a small image\nPLEASE, INCREASE THE IMAGE SIZE!\nAND/OR DECREASE THE WAVELENGTH'
      self.beamText.set_text(warn)

    self.spherePlot.view_init(elev=self.dec/self.deg2rad,azim=0)
    self.arrayPath[0][:] = 10.*self.H[1]*np.cos(self.lat)
    self.arrayPath[1][:] = 10.*self.H[0]*np.cos(self.lat)
    self.arrayPath[2][:] = 10.*np.sin(self.lat)
    self.sphereArray[0].set_data(self.arrayPath[0], self.arrayPath[1])
    self.sphereArray[0].set_3d_properties(self.arrayPath[2])



  def _onPick(self,event):

   if event.mouseevent.inaxes == self.UVPlot:

     if event.artist in self.UVPlotPlot:
       idata = np.unravel_index(event.ind,self.ravelDims)
       if event.artist == self.UVPlotPlot[0]:
          n1,n2 = self.antnum[idata[0][0]]
       else:
          n2,n1 = self.antnum[idata[0][0]]

       H = self.Horig[idata[1][0]]/self.Hfac

       newtext = self.fmtBas%(n1+1,n2+1,H)
       self.basText.set_text(newtext)
       pl.draw()
       self.canvas.draw()
     else:
       Up = event.mouseevent.xdata-self.UVSh
       Vp = event.mouseevent.ydata+self.UVSh
#       yi = np.floor((self.UVmax/2.+Up)/(self.UVmax)*self.Npix)
#       xi = np.floor((self.UVmax/2.-Vp)/(self.UVmax)*self.Npix)
       yi = np.floor((self.UVmax+Up)/(self.UVmax)*self.Npix/2.)
       xi = np.floor((self.UVmax-Vp)/(self.UVmax)*self.Npix/2.)
     #  print Up, Vp, yi, xi
       Flux = self.FFTtoplot[xi,yi]
       Phas, Amp = np.angle(Flux,deg=True), np.abs(Flux)
       newtext = self.fmtVis%(Amp,Phas)
       self.visText.set_text(newtext)
       pl.draw()
       self.canvas.draw()
   #    print Up,Vp
#       self.auxVar = event

   elif event.mouseevent.inaxes == self.beamPlot:

     RA = event.mouseevent.xdata
     Dec = event.mouseevent.ydata
   #  print RA, Dec
     yi = np.floor((self.Xaxmax-RA)/(2.*self.Xaxmax)*self.Npix)
     xi = np.floor((self.Xaxmax-Dec)/(2.*self.Xaxmax)*self.Npix)
   #  print xi, yi
     Flux = self.beam[xi,yi]
     self.beamText.set_text(self.fmtB%(Flux,RA,Dec))
     pl.draw()
     self.canvas.draw()

   elif event.mouseevent.inaxes == self.dirtyPlot:

     RA = event.mouseevent.xdata
     Dec = event.mouseevent.ydata
     yi = np.floor((self.Xaxmax-RA)/(2.*self.Xaxmax)*self.Npix)
     xi = np.floor((self.Xaxmax-Dec)/(2.*self.Xaxmax)*self.Npix)
     Flux = self.dirtymap[xi,yi]
     self.dirtyText.set_text(self.fmtD%(Flux,RA,Dec))
     pl.draw()
     self.canvas.draw()

   elif event.mouseevent.inaxes == self.modelPlot:

     RA = event.mouseevent.xdata
     Dec = event.mouseevent.ydata
     yi = np.floor((self.Xaxmax-RA)/(2.*self.Xaxmax)*self.Npix)
     xi = np.floor((self.Xaxmax-Dec)/(2.*self.Xaxmax)*self.Npix)
     Flux = self.modelim[xi,yi]
     self.modelText.set_text(self.fmtM%(Flux,RA,Dec))
     pl.draw()
     self.canvas.draw()


   elif event.mouseevent.inaxes == self.antPlot:

    if event.mouseevent.button==1 and not self.pickAnt:
     self.antidx = event.ind
     if len(self.antidx) > 1:
       self.antidx = self.antidx[-1]
     self.pickAnt = True
     self.antText.set_text(self.fmtA%self.Nant + self.fmtA2%self.antidx + self.fmtA3%tuple([1000*a for a in self.antPos[self.antidx]]))
     pl.draw()
     self.canvas.draw()




  def _onAntennaDrag(self,event):   
     if self.pickAnt:
       self.antPos[self.antidx] = [event.xdata,event.ydata]
       self.antText.set_text(self.fmtA%self.Nant + self.fmtA2%self.antidx + self.fmtA3%tuple([1000*a for a in self.antPos[self.antidx]]))
       self._setBaselines(antidx=self.antidx)
       self._plotAntennas(redo=False)
       self._setBeam(antidx=self.antidx)
       self._plotBeam(redo=False)
       self._plotDirty(redo=False)
   #    pl.draw()
       self.canvas.draw()


# Drag the sphere plot (to change source position)
     if self._onSphere:
       oldDec = self.dec/self.deg2rad
       newDec, newH0 = self.spherePlot.elev, self.spherePlot.azim

   # Limits on declination:
       if np.abs(newDec)>90.:
         self.spherePlot.view_init(elev = oldDec, azim=0.0)
         return

       newDec *= self.deg2rad
       if not self.lock:
         self.lock=True
         if newDec != self.dec and np.abs(newDec-self.lat)<np.pi/2.:
           self.widget['dec'].set_val(newDec/self.deg2rad)
           self.spherePlot.view_init(elev = newDec/self.deg2rad, azim=0.0)
           self.dec = newDec
           self._changeCoordinates()
           self.lock = False
         else:
           self.spherePlot.view_init(elev = oldDec, azim=0.0)
           self.lock = False
       else:
           self.spherePlot.view_init(elev = oldDec, azim=0.0)





  def _onRelease(self,event):
     self._onSphere = False
     if self.pickAnt:
       self.pickAnt = False
       self.antText.set_text(self.fmtA%self.Nant)
       pl.draw()
       self.canvas.draw()


  def _onRobust(self,newrob):
    self.robust = newrob
    self._changeCoordinates()


  def _onKeyLat(self,newlat):
   newlat *= self.deg2rad
   if not self.lock:
    self.lock = True
    if newlat != self.lat: 
     if np.abs(newlat-self.dec)<np.pi/2.:
      self.lat = newlat
      self._changeCoordinates()
     else:
      self.widget['lat'].set_val(self.lat/self.deg2rad)
    self.lock = False


  def _gammacorr(self,gamma):
    self.gamma = gamma
    self._plotModel(redo=False)
    pl.draw()
    self.canvas.draw()

  def _onKeyDec(self,newdec):
   newdec *= self.deg2rad
   if not self.lock:
    self.lock=True
    if newdec != self.dec: 
     if np.abs(newdec-self.lat)<np.pi/2.:
      self.dec = newdec
      self._changeCoordinates()
     else:
      self.widget['dec'].set_val(self.dec/self.deg2rad)
    self.lock = False



  def _onKeyH0(self,newH0):
   newH0 *= self.Hfac
   if not self.lock:
    self.lock = True
    if np.abs(newH0) < self.Hmax:
      self.Hcov[0] = newH0
    else:
      self.Hcov[0] = -self.Hmax
      self.widget['H0'].set_val(self.Hcov[0]/self.Hfac) 
    if self.Hcov[1]<self.Hcov[0]:
      self.Hcov[1] = self.Hcov[0]
      self.widget['H1'].set_val(self.Hcov[1]/self.Hfac) 
    self._changeCoordinates()
    self.lock=False



  def _onKeyH1(self,newH1):
   newH1 *= self.Hfac
   if not self.lock:
    self.lock = True
    if np.abs(newH1) < self.Hmax:
      self.Hcov[1] = newH1
    else:
      self.Hcov[1] = self.Hmax
      self.widget['H1'].set_val(self.Hcov[1]/self.Hfac) 
    if self.Hcov[0]>self.Hcov[1]:
      self.Hcov[0] = self.Hcov[1]
      self.widget['H0'].set_val(self.Hcov[0]/self.Hfac) 
    self._changeCoordinates()
    self.lock=False


  def _addAntenna(self,antenna):

    if self.Nant >= len(self.antPos):
      self.antPos.append([0.,0.])
      self.Nant += 1
      self.antLabelPlot.append(self.antPlot.annotate(str(self.Nant),textcoords = 'offset points',xy=(0,0),xytext=(-7,4)))

    else:
      self.antLabelPlot[self.Nant].xy = (self.antPos[self.Nant][0],self.antPos[self.Nant][1])
      self.Nant += 1



    self.antLabelPlot[self.Nant-1].set_visible(True)

    newtext = self.fmtA%self.Nant 
    self.antText.set_text(newtext)
    self._prepareBaselines()
    self._changeCoordinates()


  def _removeAntenna(self,antenna):
    if self.Nant > 2:
        self.Nant -= 1
        for i in range(self.Nant,len(self.antLabelPlot)):
           self.antLabelPlot[i].set_visible(False)
    newtext = self.fmtA%self.Nant 
    self.antText.set_text(newtext)
    self._prepareBaselines()
    self._changeCoordinates()


  def _onKeyPress(self,event):

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


    if event.key == 'Z':
      event.button = 1
      event.dblclick = True
      self._onPress(event)
    if event.key == 'z':
      event.button = 3
      event.dblclick = True
      self._onPress(event)



  def _onPress(self,event):
    if event.inaxes == self.spherePlot:
      self._onSphere = True


    if not hasattr(event,'dblclick'):
       event.dblclick = False

    if event.dblclick:
# ZOOM IN:
      if event.inaxes == self.beamPlot:
         toZoom = [self.beamPlot]
         cz = 0; inv = True; inv2 = False; scal = 1.0
      elif event.inaxes in [self.modelPlot, self.dirtyPlot]:
         toZoom = [self.modelPlot, self.dirtyPlot]
         cz = 1; inv = True; inv2 = False; scal = 1.0
      elif event.inaxes == self.UVPlot:
         toZoom = [self.UVPlot]
         cz = 2; inv = True; inv2 = True; scal = self.wavelength
      elif event.inaxes == self.antPlot:
         toZoom = [self.antPlot]
         cz = 3; inv = False; inv2 = False; scal = 1.0
      else:
         cz = -1; inv = False;inv2 = False; scal = 1.0

      if cz >= 0:
       if event.button == 1 and cz >=0:
         RA = event.xdata
         Dec = event.ydata
         xL = np.abs(self.curzoom[cz][1]-self.curzoom[cz][0])/4./scal
         yL = np.abs(self.curzoom[cz][3]-self.curzoom[cz][2])/4./scal
         x0 = RA-xL
         x1 = RA+xL
         y0 = Dec-xL
         y1 = Dec+xL
         if cz in [0,1]:
          if x0 < -self.Xaxmax/2.:
           x0 = -self.Xaxmax/2.
           x1 = x0 + 2.*xL
          if x1 > self.Xaxmax/2.:
           x1 = self.Xaxmax/2.
           x0 = x1 - 2.*xL
          if y0 < -self.Xaxmax/2.:
           y0 = -self.Xaxmax/2.
           y1 = y0 + 2.*xL
          if y1 > self.Xaxmax/2.:
           y1 = self.Xaxmax/2.
           y0 = y1 - 2.*xL
# ZOOM OUT:
       if event.button == 3:
         RA = event.xdata
         Dec = event.ydata
         xL = np.abs(self.curzoom[cz][1]-self.curzoom[cz][0])/scal
         yL = np.abs(self.curzoom[cz][3]-self.curzoom[cz][2])/scal
         if cz in [0,1]:
          if xL > self.Xaxmax/2.:
           xL = self.Xaxmax/2.
          if yL > self.Xaxmax/2.:
           yL = self.Xaxmax/2.
     #    if cz == 2:
     #     if xL > 4.*self.Xmax/self.lfac/scal:
     #      xL = 4.*self.Xmax/self.lfac/scal
     #     if yL > 4.*self.Xmax/self.lfac/scal:
     #      yL = 4.*self.Xmax/self.lfac/scal

         x0 = RA-xL
         x1 = RA+xL
         y0 = Dec-xL
         y1 = Dec+xL
         if cz in [0,1]:
          if x0 < -self.Xaxmax/2.:
           x0 = -self.Xaxmax/2.
           x1 = x0 + 2.*xL
          if x1 > self.Xaxmax/2.:
           x1 = self.Xaxmax/2.
           x0 = x1 - 2.*xL
          if y0 < -self.Xaxmax/2.:
           y0 = -self.Xaxmax/2.
           y1 = y0 + 2.*xL
          if y1 > self.Xaxmax/2.:
           y1 = self.Xaxmax/2.
           y0 = y1 - 2.*xL

       if inv:
         xx0 = x1 ; xx1 = x0
       else:
         xx0 = x0 ; xx1 = x1
       if inv2:
         yy0 = y1 ; yy1 = y0
       else:
         yy0 = y0 ; yy1 = y1

       for ax in toZoom:
           ax.set_xlim((xx0,xx1))
           ax.set_ylim((yy0,yy1))

       self.curzoom[cz] = (xx0*scal,xx1*scal,yy0*scal,yy1*scal)

       pl.draw()
       self.canvas.draw()

  def saveArray(self,array):

    fname = tkFileDialog.asksaveasfilename(defaultextension='.array',title='Save current array...')
    iff = open(fname,'w')

    print >> iff,'LATITUDE % 3.1f'%(self.lat/self.deg2rad)
    print >> iff,'DECLINATION % 3.1f'%(self.dec/self.deg2rad)
    toprint = tuple([l/self.Hfac for l in self.Hcov])
    print >> iff,'HOUR_ANGLE % 3.1f  % 3.1f'%toprint
    for ant in self.antPos:
      toprint = tuple([p*1.e3 for p in ant])
      print >> iff,'ANTENNA  % .3e   % .3e'%toprint

    self.antText.set_text('SAVED: %s'%os.path.basename(fname))
    pl.draw()
    self.canvas.draw()

    time.sleep(3)
    self.antText.set_text(self.fmtA%self.Nant) 
    pl.draw()
    self.canvas.draw()

    iff.close()

          
    return


  def loadArray(self,array):

    antenna_file = tkFileDialog.askopenfilename(title='Load array...',initialdir=self.arraydir)
    self.lock=False

    if len(antenna_file)>0:
      goodread = self.readAntennas(str(antenna_file))

      if goodread:
        self._prepareBaselines()
        self._changeCoordinates()
        self._plotModel(redo=True)
        self._plotModelFFT(redo=True) 
        newtext = self.fmtA%self.Nant 
        self.antText.set_text(newtext)
        self.widget['lat'].set_val(self.lat/self.deg2rad)
        self.widget['dec'].set_val(self.dec/self.deg2rad)
        self.widget['H0'].set_val(self.Hcov[0]/self.Hfac)
        self.widget['H1'].set_val(self.Hcov[1]/self.Hfac)
        self.widget['wave'].set_val(self.wavelength*1.e6)
        self._plotAntennas(redo=False,rescale=True)
      #  self._setBaselines()
        self._plotBeam(redo=True)
        self._plotDirty(redo=True)
        pl.draw()
        self.canvas.draw()

    return


  def loadModel(self,model):

    model_file = tkFileDialog.askopenfilename(title='Load model...',initialdir=self.modeldir)
    self.lock=False

    if len(model_file)>0:
      goodread = self.readModels(str(model_file))
    #  print goodread
      if goodread:
        self._prepareModel()
        self._plotModel(redo=True)
        self.widget['wave'].set_val(self.wavelength*1.e6)
        self._setBaselines()
        self._setBeam()
        self._changeCoordinates()
        self._plotModelFFT(redo=True) 
        self.wax['wave'].cla()
        self.widget['wave'] = Slider(self.wax['wave'],r'$\lambda$ (mm)',self.wavelength*1.e6/10.,10.*self.wavelength*1.e6,valinit=self.wavelength*1.e6)
        self.widget['wave'].on_changed(self._changeWavelength)
        self._plotBeam(redo=True)
        self._plotDirty(redo=True)
        pl.draw()
        self.canvas.draw()


    return














class CLEANer(object):

  def quit(self):

    self.me.destroy()
      

  def __init__(self,parent):

    self.parent = parent
    self.me = Tk.Toplevel(parent.tks)

    menubar = Tk.Menu(self.me)
    menubar.add_command(label="Help", command=self._getHelp)
    menubar.add_command(label="Quit", command=self.quit)

    self.me.config(menu=menubar)

    self.Np4 = self.parent.Npix/4

    self.figCL1 = pl.figure(figsize=(5,5))    
    self.figCL2 = pl.figure(figsize=(5,5))    

  #  self.ResidPlot = self.figCL.add_subplot(231,aspect='equal')
  #  self.CLEANPlot = self.figCL.add_subplot(233,aspect='equal')

    self.ResidPlot = self.figCL1.add_subplot(111) #pl.axes([0.01,0.43,0.5,0.5],aspect='equal')
    self.CLEANPlot = self.figCL2.add_subplot(111) #pl.axes([0.55,0.43,0.5,0.5],aspect='equal')

    self.frames = {}
    self.frames['FigFr'] = Tk.Frame(self.me)
    self.frames['GFr'] = Tk.Frame(self.me)

    self.canvas1 = FigureCanvasTkAgg(self.figCL1, master=self.frames['FigFr'])
    self.canvas2 = FigureCanvasTkAgg(self.figCL2, master=self.frames['FigFr'])

    self.canvas1.show()
    self.canvas2.show()

    self.frames['FigFr'].pack(side=Tk.TOP)
    self.frames['GFr'].pack(side=Tk.TOP)


    self.frames['CLOpt'] = Tk.Frame(self.frames['FigFr'])
    self.frames['Gain'] = Tk.Frame(self.frames['CLOpt'])
    self.frames['Niter'] = Tk.Frame(self.frames['CLOpt'])
    Gtext = Tk.Label(self.frames['Gain'],text="Gain:  ")
    Ntext = Tk.Label(self.frames['Niter'],text="# iter:")

    self.entries = {}
    self.entries['Gain'] = Tk.Entry(self.frames['Gain'])
    self.entries['Gain'].insert(0,"0.1")
    self.entries['Gain'].config(width=5)
    self.entries['Niter'] = Tk.Entry(self.frames['Niter'])
    self.entries['Niter'].insert(0,"100")
    self.entries['Niter'].config(width=5)


    GTitle = Tk.Label(self.frames['GFr'],text="CALIBRATION ERROR:")
    GTitle.pack(side=Tk.TOP)
    self.frames['Ant1L'] = Tk.Frame(self.frames['GFr'])
    Ant1T = Tk.Label(self.frames['Ant1L'],text="Ant. 1:")
    Ant1T.pack(side=Tk.TOP)
    self.entries['Ant1'] = Tk.Listbox(self.frames['Ant1L'],exportselection=False,width=5)
    self.frames['Ant2L'] = Tk.Frame(self.frames['GFr'])
    Ant2T = Tk.Label(self.frames['Ant2L'],text="Ant. 2:")
    Ant2T.pack(side=Tk.TOP)
    self.entries['Ant2'] = Tk.Listbox(self.frames['Ant2L'],exportselection=False,width=5)
    self.entries['Ant1'].pack(side=Tk.TOP)
    self.entries['Ant2'].pack(side=Tk.TOP)

    self.frames['Ant1L'].pack(side=Tk.LEFT)
    self.frames['Ant2L'].pack(side=Tk.LEFT)

    self.frames['space1'] = Tk.Frame(self.frames['GFr'],width=40)
    self.frames['space1'].pack(side=Tk.LEFT)

    self.frames['H0Fr'] = Tk.Frame(self.frames['GFr'])
    self.frames['H1Fr'] = Tk.Frame(self.frames['GFr'])
    self.frames['AmpFr'] = Tk.Frame(self.frames['GFr'])
    self.frames['PhasFr'] = Tk.Frame(self.frames['GFr'])

    self.frames['H0Fr'].pack(side=Tk.TOP)
    self.frames['H1Fr'].pack(side=Tk.TOP)
    self.frames['AmpFr'].pack(side=Tk.TOP)
    self.frames['PhasFr'].pack(side=Tk.TOP)


    self.entries['H0'] = Tk.Scale(self.frames['H0Fr'],from_=0,to=self.parent.nH,orient=Tk.HORIZONTAL,length=200)
    self.entries['H1'] = Tk.Scale(self.frames['H1Fr'],from_=0,to=self.parent.nH,orient=Tk.HORIZONTAL,length=200)
    H0Text = Tk.Label(self.frames['H0Fr'],text="From integration #: ",width=15)
    H1Text = Tk.Label(self.frames['H1Fr'],text="To integration #: ",width=15)
    H0Text.pack(side=Tk.LEFT)
    self.entries['H0'].pack(side=Tk.RIGHT)
    H1Text.pack(side=Tk.LEFT)
    self.entries['H1'].pack(side=Tk.RIGHT)
    self.entries['H1'].set(self.parent.nH)

    self.entries['Amp'] = Tk.Scale(self.frames['AmpFr'],from_=0.1,to=10.,orient=Tk.HORIZONTAL,length=200)
    self.entries['Phas'] = Tk.Scale(self.frames['PhasFr'],from_=-180.,to=180.,orient=Tk.HORIZONTAL,length=200)
    AmpText = Tk.Label(self.frames['AmpFr'],text="Amplitude gain: ",width=15)
    PhasText = Tk.Label(self.frames['PhasFr'],text="Phase Gain: ",width=15)
    AmpText.pack(side=Tk.LEFT)
    self.entries['Amp'].pack(side=Tk.RIGHT)
    PhasText.pack(side=Tk.LEFT)
    self.entries['Phas'].pack(side=Tk.RIGHT)


    Gtext.pack(side=Tk.LEFT)
    self.entries['Gain'].pack(side=Tk.RIGHT)

    Ntext.pack(side=Tk.LEFT)
    self.entries['Niter'].pack(side=Tk.RIGHT)

    self.canvas1.get_tk_widget().pack(side=Tk.LEFT) #, fill=Tk.BOTH, expand=1)
    self.frames['CLOpt'].pack(side=Tk.LEFT)
    self.canvas2.get_tk_widget().pack(side=Tk.LEFT) #, fill=Tk.BOTH, expand=1)

    self.buttons = {}
    self.buttons['clean'] = Tk.Button(self.frames['CLOpt'],text="CLEAN",command=self._CLEAN)
    self.buttons['reset'] = Tk.Button(self.frames['CLOpt'],text="RESET",command=self._reset)
    self.buttons['addres'] = Tk.Button(self.frames['CLOpt'],text="+/- Resid",command=self._AddRes)

    self.buttons['apply'] = Tk.Button(self.frames['GFr'],text="APPLY GAIN",command=self._ApplyGain)
    self.buttons['apply'].pack(side=Tk.RIGHT)
    self.buttons['recal'] = Tk.Button(self.frames['GFr'],text="RESET GAIN",command=self._reCalib)
    self.buttons['recal'].pack(side=Tk.RIGHT)

    self.frames['Gain'].pack(side=Tk.TOP)
    self.frames['Niter'].pack(side=Tk.TOP)
    self.buttons['clean'].pack(side=Tk.TOP)
    self.buttons['reset'].pack(side=Tk.TOP)
    self.buttons['addres'].pack(side=Tk.TOP)

    self.canvas1.mpl_connect('pick_event', self._onPick)
    self.canvas2.mpl_connect('pick_event', self._onPick)
    self.canvas1.mpl_connect('motion_notify_event', self._doMask)
    self.canvas1.mpl_connect('button_release_event',self._onRelease)
    self.canvas1.mpl_connect('button_press_event',self._onPress)

    self.pressed = -1
    self.xy0 = [0,0]
    self.moved = False
    self.resadd = False
    self._makeMask()
    self._reCalib()


  def _ApplyGain(self):

    try:
      an1 = int(self.entries['Ant1'].curselection()[0])
    except:
      return

    try:
      an2 = int(self.entries['Ant2'].curselection()[0])
    except:
      an2=an1


    if an2==an1:
      an2=-1

    G = float(self.entries['Amp'].get())*np.exp(1.j*float(self.entries['Phas'].get())*np.pi/180.)
    H0 = int(self.entries['H0'].get())
    H1 = int(self.entries['H1'].get())

    self.parent._setGains(an1,an2,H0,H1,G)
    self._reset()

  def _makeMask(self):

    self.mask = np.zeros(np.shape(self.parent.beam))
    self.bmask = np.zeros(np.shape(self.parent.beam)).astype(np.bool)


  def _onPick(self,event):

     RA = event.mouseevent.xdata
     Dec = event.mouseevent.ydata
  #   print RA, Dec
     yi = np.floor((self.Xaxmax-RA)/(2.*self.Xaxmax)*self.parent.Npix)
     xi = np.floor((self.Xaxmax-Dec)/(2.*self.Xaxmax)*self.parent.Npix)
  #   print xi, yi
     Flux = self.residuals[xi,yi]
     self.ResidText.set_text(self.parent.fmtM%(Flux,RA,Dec))
     if self.resadd:
       Flux = self.cleanmod[xi,yi] + self.residuals[xi,yi]
     else:
       Flux = self.cleanmod[xi,yi]

     self.CLEANText.set_text(self.parent.fmtM%(Flux,RA,Dec))


     self.canvas1.draw()
     self.canvas2.draw()


  def _onPress(self,event):

    if event.inaxes == self.ResidPlot:
      self.pressed = int(event.button)
 #   print 'WAS PRESSED:',self.pressed, event.xdata, event.ydata
 #   self.xy0 = [int(event.xdata),int(event.ydata)]
      RA = event.xdata
      Dec = event.ydata
      self.xydata = [RA,Dec]
      self.xy0[1] = np.floor((self.Xaxmax-RA)/(2.*self.Xaxmax)*self.parent.Npix)
      self.xy0[0] = np.floor((self.Xaxmax-Dec)/(2.*self.Xaxmax)*self.parent.Npix)
    #  self.Box.set_size=[0.,0.]
    #  self.Box.set_xy(self.xy0)
      self.moved = False 

  def _onRelease(self,event):

    if event.inaxes != self.ResidPlot:
      self.moved=False

    if self.moved:
      RA = event.xdata
      Dec = event.ydata
      y1 = np.floor((self.Xaxmax-RA)/(2.*self.Xaxmax)*self.parent.Npix)
      x1 = np.floor((self.Xaxmax-Dec)/(2.*self.Xaxmax)*self.parent.Npix)
      xi,xf = [min(self.xy0[0],x1),max(self.xy0[0],x1)]
      yi,yf = [min(self.xy0[1],y1),max(self.xy0[1],y1)]
      if self.pressed==1:
        self.mask[xi:xf,yi:yf] = 1.0
        self.bmask[xi:xf,yi:yf] = True

   #     print 'DONE'
      else:
        self.mask[xi:xf,yi:yf] = 0.0
        self.bmask[xi:xf,yi:yf] = False

   #     print 'DONE2'

      for coll in self.MaskPlot.collections:
         self.ResidPlot.collections.remove(coll)

      self.MaskPlot = self.ResidPlot.contour(np.linspace(self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Npix/2),np.linspace(self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Npix/2),self.mask[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4],levels=[0.5])

#      self.MaskPlot.set_array(self.mask[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4])
      self.ResidPlot.set_xlim((self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.))
      self.ResidPlot.set_ylim((-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))

      self.canvas1.draw()
      self.canvas2.draw()

 #     print 'WAS RELEASED ',self.moved,np.sum(self.mask),self.xy0,x1,y1
      self.Box.set_data([0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.])


    self.moved = False
    self.pressed = -1
    self.canvas1.draw()

  def _doMask(self,event):
    if self.pressed>=0 and event.inaxes==self.ResidPlot:
      self.moved = True
      RA = event.xdata
      Dec = event.ydata
      y1 = np.floor((self.Xaxmax-RA)/(2.*self.Xaxmax)*self.parent.Npix)
      x1 = np.floor((self.Xaxmax-Dec)/(2.*self.Xaxmax)*self.parent.Npix)
      self.Box.set_data([self.xydata[0],self.xydata[0],RA,RA,self.xydata[0]],[self.xydata[1],Dec,Dec,self.xydata[1],self.xydata[1]])
      self.canvas1.draw()

  def _AddRes(self):

    if self.resadd:
      self.resadd = False
      toadd = self.cleanmod[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4]
    else:
      self.resadd = True
      toadd = (self.cleanmod + self.residuals)[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4]

    self.CLEANPlotPlot.set_array(toadd)
    self.CLEANPlotPlot.norm.vmin = np.min(toadd)
    self.CLEANPlotPlot.norm.vmax = np.max(toadd)

    self.canvas1.draw()
    self.canvas2.draw()
    del toadd


  def _reCalib(self):

    self.entries['Ant1'].delete(0,Tk.END)
    self.entries['Ant2'].delete(0,Tk.END)
    self.entries['H0'].set(0)
    self.entries['H1'].set(self.parent.nH)
    self.entries['Amp'].set(1)
    self.entries['Phas'].set(0)

    for i in range(self.parent.Nant):
      self.entries['Ant1'].insert(Tk.END,str(i+1))
      self.entries['Ant2'].insert(Tk.END,str(i+1))

    self.parent._setGains(-1,-1,0,0,1.0)
    self._reset()


  def _reset(self):


   # self.parent.Gains[:] = 1.0

  #  self.parent._prepareBaselines()
  #  self.parent._changeCoordinates()


   # self.dirtymap[:] = np.fft.fftshift(np.fft.ifft2(self.modelfft*np.fft.ifftshift(self.robustsamp))).real/self.beamScale
    extr = [np.min(self.parent.dirtymap),np.max(self.parent.dirtymap)]
   # Nphf = int(self.Npix/2)

    self.ResidPlot.cla()

    self.ResidPlotPlot = self.ResidPlot.imshow(self.parent.dirtymap[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4],interpolation='nearest',picker=True, cmap=self.parent.currcmap)
    modflux = self.parent.dirtymap[self.parent.Nphf,self.parent.Nphf]
    self.ResidText = self.ResidPlot.text(0.05,0.87,self.parent.fmtD%(modflux,0.0,0.0),
         transform=self.ResidPlot.transAxes,bbox=dict(facecolor='white', alpha=0.7))
    pl.setp(self.ResidPlotPlot, extent=(self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))

 #   self.curzoom[1] = (self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.)
    self.Xaxmax = float(self.parent.Xaxmax)

    self.Box = self.ResidPlot.plot([0.,0.,0.,0.,0.],[0.,0.,0.,0.,0.],lw=2,color='w')[0]

    self.ResidPlot.set_ylabel('Dec offset (as)')
    self.ResidPlot.set_xlabel('RA offset (as)')
    self.ResidPlot.set_title('RESIDUALS')

    self.MaskPlot = self.ResidPlot.contour(np.linspace(self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Npix/2),np.linspace(self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Npix/2),self.mask[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4],levels=[0.5])

    self.ResidPlot.set_xlim((self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.))
    self.ResidPlot.set_ylim((-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))


    self.residuals = np.copy(self.parent.dirtymap)
    self.cleanmod = np.zeros(np.shape(self.parent.dirtymap))

    self.CLEANPlot.cla()
    self.CLEANPlotPlot = self.CLEANPlot.imshow(self.parent.dirtymap[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4],interpolation='nearest',picker=True, cmap=self.parent.currcmap)
    modflux = self.parent.dirtymap[self.parent.Nphf,self.parent.Nphf]
    self.CLEANText = self.CLEANPlot.text(0.05,0.87,self.parent.fmtD%(0.0,0.0,0.0),
         transform=self.CLEANPlot.transAxes,bbox=dict(facecolor='white', alpha=0.7))
    pl.setp(self.CLEANPlotPlot, extent=(self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,-self.parent.Xaxmax/2.,self.parent.Xaxmax/2.))
    self.CLEANPlot.set_ylabel('Dec offset (as)')
    self.CLEANPlot.set_xlabel('RA offset (as)')
    self.CLEANPlot.set_title('CLEAN (0 ITER)')
    self.CLEANPlotPlot.set_array(self.cleanmod[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4])

    self.totiter = 0

    # DERIVE THE CLEAN BEAM
    MainLobe = np.where(self.parent.beam>0.5)
    dX = MainLobe[0]-self.parent.Npix/2 ; dY = MainLobe[1]-self.parent.Npix/2
    fit = spfit.leastsq(lambda x: np.exp(-(dX*dX*x[0]+dY*dY*x[1]+dX*dY*x[2]))-self.parent.beam[MainLobe],[0.,0.,0.])
    print 'BEAM FIT: ',fit[0]
    self.cleanBeam = np.zeros(np.shape(self.residuals))
    ddX = np.outer(np.ones(self.parent.Npix),np.arange(-self.parent.Npix/2,self.parent.Npix/2).astype(np.float64))
    ddY = np.outer(np.arange(-self.parent.Npix/2,self.parent.Npix/2).astype(np.float64),np.ones(self.parent.Npix))

   # self.cleanBeam[:] = np.exp(-(ddX*ddX*fit[0][0]+ddY*ddY*fit[0][1]+ddY*ddX*fit[0][2]))
    self.cleanBeam[:] = np.exp(-(ddY*ddY*fit[0][0]+ddX*ddX*fit[0][1]+ddY*ddX*fit[0][2]))

    self.resadd = False
    self.ffti = False

    self.canvas1.draw()
    self.canvas2.draw()

    del ddX, ddY, modflux

  def _CLEAN(self):

     if np.sum(self.bmask)==0:
       goods = np.ones(np.shape(self.bmask)).astype(np.bool)
       tempres = self.residuals
     else:
       goods = self.bmask
       tempres = self.residuals*self.mask

     psf = self.parent.beam
     gain = float(self.entries['Gain'].get())
     niter = int(self.entries['Niter'].get())
     for i in range(niter):
       self.totiter += 1
       peakpos = np.unravel_index(np.argmax(tempres),np.shape(self.residuals))
       peakval = self.residuals[peakpos[0],peakpos[1]]
       self.residuals -= gain*peakval*np.roll(np.roll(psf,peakpos[0]-self.parent.Npix/2,axis=0), peakpos[1]-self.parent.Npix/2,axis=1)
       tempres[goods] = self.residuals[goods]
       # MODIFY CLEAN MODEL!!
       self.cleanmod += gain*peakval*np.roll(np.roll(self.cleanBeam,peakpos[0]-self.parent.Npix/2,axis=0), peakpos[1]-self.parent.Npix/2,axis=1)
       self.ResidPlotPlot.set_array(self.residuals[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4])
       self.CLEANPlot.set_title('CLEAN (%i ITER)'%self.totiter)

       if self.resadd:
         toadd = (self.cleanmod + self.residuals)[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4]
       else:
         toadd = self.cleanmod[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4]

       self.CLEANPlotPlot.set_array(toadd)
       self.CLEANPlotPlot.norm.vmin = np.min(toadd)
       self.CLEANPlotPlot.norm.vmax = np.max(toadd)

       self.canvas1.draw()
       self.canvas2.draw()

     del psf, toadd,tempres, goods


  def _getHelp(self):
    win = Tk.Toplevel(self.me)
    win.title("Help")
    helptext = ScrolledText(win)
    helptext.config(state=Tk.NORMAL)
   # about = __help_text__
    helptext.insert('1.0',__CLEAN_help_text__)
    helptext.config(state=Tk.DISABLED)

   # Tk.Label(win, text=about, width=100, height=10).pack()
    helptext.pack()
    Tk.Button(win, text='OK', command=win.destroy).pack()



#  def FFT(self):
#     Np4 = self.parent.Npix/4
#     if self.ffti:
#        self.ffti=False
#        self.CLEANPlot.set_ylabel('Dec offset (as)')
#        self.CLEANPlot.set_xlabel('RA offset (as)')
#        self.ResidPlot.set_ylabel('Dec offset (as)')
#        self.ResidPlot.set_xlabel('RA offset (as)')
#        self.ResidPlotPlot.set_array(self.residuals[self.Np4:self.parent.Npix-self.Np4,self.Np4:self.parent.Npix-self.Np4])


  #   print 'HILIDILI'











#oldtb = mpl.backends.backend_gtkagg.NavigationToolbar2GTKAgg

#class CustomToolbar(oldtb):
#   def __init__(self,canvas_,parent_):
#     self.toolitems = (
#         ('Home', 'Lorem ipsum dolor sit amet', 'home', 'home'),
#         ('Back', 'consectetuer adipiscing elit', 'back', 'back'),
#         ('Forward', 'sed diam nonummy nibh euismod', 'forward', 'forward'),
#         (None, None, None, None),
#         ('Pan', 'tincidunt ut laoreet', 'move', 'pan'),
#         ('Zoom', 'dolore magna aliquam', 'zoom_to_rect', 'zoom'),
#         (None, None, None, None),
#         ('Subplots', 'putamus parum claram', 'subplots', 'configure_subplots'),
#         ('Save', 'sollemnes in futurum', 'filesave', 'save_figure'),
#         )
#     oldtb.__init__(self,canvas_,parent_)


#class MyApp(object):
#    def __init__(self,root):
#        self.root = root
#        self._init_app()
#
#    # here we embed the figure in the Tk GUI
#    def _init_app(self):
#        self.interferometer = Interferometer()
#      #  self.figure = mpl.figure.Figure()
#      #  self.ax = self.figure.add_subplot(111)
#        self.canvas = FigureCanvasTkAgg(self.interferometer.figUV,self.root)
#        self.toolbar = CustomToolbar(self.canvas,self.root)
#        self.toolbar.update()
#        self.plot_widget = self.canvas.get_tk_widget()
#        self.plot_widget.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
#        self.toolbar.pack(side=tk.TOP, fill=tk.BOTH, expand=1)
#        self.canvas.show()
##
##
##
##
#def main():
#    root = Tkinter.Tk()
#    app = MyApp(root)
#    root.mainloop()
#
#if __name__ == "__main__":
#    main()

if __name__ == "__main__":

  root = Tk.Tk()
  TITLE = 'Aperture Synthesis Simulator (I. Marti-Vidal, Onsala Space Observatory) - version  %s'%__version__
  root.wm_title(TITLE)


  myint = Interferometer(tkroot=root)
  Tk.mainloop()

