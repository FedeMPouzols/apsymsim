# -*- mode: python -*-
a = Analysis(['z:\\APSYNSIM-v0.2\\SOURCE\\APSYNSIM-v0.2.py'],
             pathex=['Z:\\APSYNSIM-v0.2\\COMPILE'],
             hiddenimports=[],
             hookspath=None,
             runtime_hooks=None)
pyz = PYZ(a.pure)
exe = EXE(pyz,
          a.scripts,
          exclude_binaries=True,
          name='APSYNSIM-v0.2.exe',
          debug=False,
          strip=None,
          upx=True,
          console=True , icon='z:\\APSYNSIM-v0.2\\COMPILE\\APSYNSIM_icon_small.ico')
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=None,
               upx=True,
               name='APSYNSIM-v0.2')
