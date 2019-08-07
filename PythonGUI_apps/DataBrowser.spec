# -*- mode: python -*-

block_cipher = None


a = Analysis(['C:\\Users\\lindat18\\Dropbox\\Ginger_Lab\\Data_Analysis\\PythonGUI_apps\\DataBrowser.py'],
             pathex=['C:\\Users\\lindat18\\Dropbox\\Ginger_Lab\\Data_Analysis\\PythonGUI_apps'],
             binaries=[],
             datas=[],
             hiddenimports=['ipykernel.datapub'],
             hookspath=[],
             runtime_hooks=[],
             excludes=[],
             win_no_prefer_redirects=False,
             win_private_assemblies=False,
             cipher=block_cipher,
             noarchive=False)
pyz = PYZ(a.pure, a.zipped_data,
             cipher=block_cipher)
exe = EXE(pyz,
          a.scripts,
          [],
          exclude_binaries=True,
          name='DataBrowser',
          debug=False,
          bootloader_ignore_signals=False,
          strip=False,
          upx=True,
          console=True )
coll = COLLECT(exe,
               a.binaries,
               a.zipfiles,
               a.datas,
               strip=False,
               upx=True,
               name='DataBrowser')
