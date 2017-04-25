# -*- mode: python -*-

# SPEC FILE FOR PYINSTALLER:

# RUN PYINSTALLER WITH THE -w FLAG!!

a = Analysis(['./script/apsynsim.py', './script/cleaner.py',
              './script/interferometer.py', './script/simple_clean_img.py',
              './script/uvplotter1.py', './script/uvplotter2.py'],
             pathex=['./apsinsym_mac'],
             hiddenimports=['scipy.special._ufuncs_cxx','mpl_toolkits','mpl_toolkits.mplot3d'],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='apsynsim.app',
          debug=False,
          strip=None,
          upx=True,
          console=True , icon='./compile/apsynsim_icon_small.ico')

import mpl_toolkits.mplot3d
import os


coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas, 
               strip=None,
               upx=True,
               name='apsynsim')

