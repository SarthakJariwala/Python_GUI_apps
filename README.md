# GLabViz
GUI Python apps written in python and qt for quick analysis of custom data. It also includes the ability to convert data to *H5* if need to be analysed in other software such as MATLAB.

_**Python is not required to use GLabViz**_ (see **How to use?**)

The primary users for this Python package are Ginger Lab members at the University of Washington, Seattle but is licensed under MIT License and open for everyone to use.

## Includes
* Fluorescence Lifetime Analysis
* Spectra Analysis
* Fluorescence Lifetime Imaging Microscopy (FLIM) Data Analysis
* Photluminescence Quantum Efficiency (PLQE) Analysis
* UV-Vis Data Analysis
* General *H5* View and Plot
* *H5* and *PKL* File Viewer

## How to use?
### Standalone App - without Python or any dependencies (_only for Windows users_)
* Under the releases page, download the latest release of the _**DataBrowser**_ zip file
* Extract the zip file and run _**DataBrowser.exe**_
### With Python and its dependencies
```
git clone https://github.com/SarthakJariwala/Python_GUI_apps.git
```
* Install all dependencies
* Run the application by double-clicking DataBrowser.py.
* OR Run it from command-line while in the PythonGUI_apps folder:
```
python DataBrowser.py
```

#### Dependencies
* ScopeFoundry
* pyqtgraph 
* numpy
* pyqt
* qtpy
* h5py
* matplotlib
* scipy
* lmfit
* customplotting

#### Installing dependencies from command-line
```
conda install numpy pyqt qtpy h5py pyqtgraph
pip install ScopeFoundry
pip install matplotlib scipy lmfit customplotting==0.1.4.dev0
```

## Screenshots
### Welcome Screen
![Welcome Screen](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_interface_1.png)
