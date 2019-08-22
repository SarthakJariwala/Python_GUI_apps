# GLabViz
GUI Python apps written in python and qt for quick analysis of custom data. It also includes the ability to convert data to *H5* if need to be analysed in other software such as MATLAB.

_**Python is not required to use GLabViz**_ (see **How to use?**)

[_**DOWNLOAD HERE**_](https://github.com/SarthakJariwala/Python_GUI_apps/releases)

The primary users for this Python package application are Ginger Lab members at the University of Washington, Seattle but is licensed under MIT License and open for everyone to use.

## Includes
* Fluorescence Lifetime Analysis
    * analyze lifetime
    * fit data with or without IRF
    * fit with stretched, single, or double exponential functions by diff_ev or fmin_tnc
    * calculate surface recombination velocity
    * export graph and fit results
* Spectra Analysis
    * analyze single spectrum
        * fit with or without background and white light
        * fit with single Lorentzian, single Gaussian, double Gaussian, triple Gaussian models
        * export graph and fit results
    * analyze spectra scan
        * load spectra scan data in .h5 or .pkl files
        * plot raw scan data
        * plot scan intensity sums
        * plot fitted scan by pk_pos, fwhm, sigma, or height
        * export fitted scan
    * .pkl to .txt, .pkl to .h5 converters
* Fluorescence Lifetime Imaging Microscopy (FLIM) Data Analysis
    * load lifetime scans in .h5 or .pkl files
    * plot histogram intensity sums and analyze PSF
    * export intensities array and intensities image
    * plot raw histogram data and analyze lifetime
    * compare lifetime in two different regions
* Photluminescence Quantum Efficiency (PLQE) Analysis
    * plot PLQE data
    * calculate PLQE
* UV-Vis Data Analysis
    * plot UV-Vis data
    * correct UV-Vis data for scattering
    * plot Tauc data
    * calculate bandgap
    * export UV-Vis and Tauc plots
* General *H5* View and Plot
    * load .h5 file to view file structure
    * plot datasets as a graph or an image
* *H5* and *PKL* File Viewer
    * load .h5 or .pkl file to view file structure
* Image Analysis
    * load image on SPOT or Pixera settings, or specify pixel size
    * handle RGB and greyscale images 
    * select magnification
    * color profile horizontally or vertically

## Screenshots
### Welcome Screen
![Welcome Screen](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_interface_1.png)
### Lifetime Analysis
![Lifetime Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Lifetime_analysis_2.PNG)
### Spectra Analysis
![Single Spectrum](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Spectrum_analysis_1.PNG)
![Scan Data](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Spectrum_analysis_2.PNG)
### FLIM Analysis
![FLIM Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_FLIM_analysis_2.png)
### UV-Vis Analysis
![UV-Vis Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_UVvis_analysis_1.PNG)
### H5 & Pkl View
![H5-pkl-viewer](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_h5_ViewPlot_analysis_1.PNG)
### H5 Quick Plot
![h5- quick plot](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_h5_ViewPlot_analysis_2.PNG)
### Image Analysis
![Image Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Image_analysis_1.PNG.png)

## How to use?
### Standalone App - without Python or any dependencies (_only for Windows users_)
* Under the [releases](https://github.com/SarthakJariwala/Python_GUI_apps/releases) page, download the latest release of the _**DataBrowser**_ zip file
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
* [ScopeFoundry](https://github.com/ScopeFoundry/ScopeFoundry)
* [pyqtgraph](http://www.pyqtgraph.org/) 
* numpy
* pyqt
* qtpy
* h5py
* matplotlib
* scipy
* [lmfit](https://lmfit.github.io/lmfit-py/)
* [customplotting](https://github.com/SarthakJariwala/Custom-Plotting)

#### Installing dependencies from command-line
```
conda install numpy pyqt qtpy h5py pyqtgraph
pip install ScopeFoundry
pip install matplotlib scipy lmfit customplotting==0.1.4.dev0
```
