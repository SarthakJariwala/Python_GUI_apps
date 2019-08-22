# GLabViz
GUI Python apps written in python and qt for quick analysis of custom data. It also includes the ability to convert data to *H5* if need to be analysed in other software such as MATLAB.

_**Python is not required to use GLabViz**_ (see **How to use?**)

[_**DOWNLOAD HERE**_](https://github.com/SarthakJariwala/Python_GUI_apps/releases)

The primary users for this Python package application are Ginger Lab members at the University of Washington, Seattle but is licensed under MIT License and open for everyone to use.

## Includes
* Fluorescence Lifetime Analysis
    * Analyze lifetime
    * Fit data with or without IRF
    * Fit with stretched, single, or double exponential functions by diff_ev or fmin_tnc
    * Calculate surface recombination velocity
    * Export graph and fit results
* Spectra Analysis
    * Analyze single spectrum
        * Fit with or without background and white light
        * Fit with single Lorentzian, single Gaussian, double Gaussian, triple Gaussian models
        * Export graph and fit results
    * Analyze spectra scan
        * Load spectra scan data in .h5 or .pkl files
        * Plot raw scan data
        * Plot scan intensity sums
        * Plot fitted scan by pk_pos, fwhm, sigma, or height
        * Export fitted scan
    * .pkl to .txt, .pkl to .h5 converters
* Fluorescence Lifetime Imaging Microscopy (FLIM) Data Analysis
    * Load lifetime scans in .h5 or .pkl files
    * Plot histogram intensity sums and analyze PSF
    * Export intensities array and intensities image
    * Plot raw histogram data and analyze lifetime
    * Compare lifetime in two different regions
* Photluminescence Quantum Efficiency (PLQE) Analysis
    * Plot PLQE data
    * Calculate PLQE
* UV-Vis Data Analysis
    * Plot UV-Vis data
    * Correct UV-Vis data for scattering
    * Plot Tauc data
    * Calculate bandgap
    * Export UV-Vis and Tauc plots
* General *H5* View and Plot
    * Load .h5 file to view file structure
    * Plot datasets as a graph or an image
* *H5* and *PKL* File Viewer
    * Load .h5 or .pkl file to view file structure
* Image Analysis
    * Load image on SPOT or Pixera settings, or specify pixel size
    * Handle RGB and greyscale images 
    * Select magnification
    * Color profile horizontally or vertically

## Screenshots
### Welcome Screen
![Welcome Screen](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_interface_1.png)
### Lifetime Analysis
![Lifetime Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Lifetime_analysis_2.png)
### Spectra Analysis
![Single Spectrum](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Spectrum_analysis_1.png)
![Scan Data](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Spectrum_analysis_2.png)
### FLIM Analysis
![FLIM Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_FLIM_analysis_2.png)
### UV-Vis Analysis
![UV-Vis Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_UVvis_analysis_1.PNG)
### H5 & Pkl View
![H5-pkl-viewer](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_h5_ViewPlot_analysis_1.PNG)
### H5 Quick Plot
![h5- quick plot](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_h5_ViewPlot_analysis_2.PNG)
### Image Analysis
![Image Analysis](https://github.com/SarthakJariwala/Python_GUI_apps/blob/master/Screenshots/GLabViz_Image_analysis_1.png)

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
