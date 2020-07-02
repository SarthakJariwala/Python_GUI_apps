# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

# system imports
import sys
import h5py
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import pickle
import lmfit
from lmfit.models import GaussianModel, LinearModel
from scipy import interpolate
import customplotting.mscope as cpm

sys.path.append(os.path.abspath('../H5_Pkl'))
from H5_Pkl import h5_pkl_view
sys.path.append(os.path.abspath('../Export_Windows'))
try:
    from Export_window import ExportFigureWindow, ExportPlotWindow
except:
    from Export_Windows.Export_window import ExportFigureWindow, ExportPlotWindow
# local modules
try:
    from Spectra_fit_funcs import Spectra_Fit, Single_Gaussian, Single_Lorentzian, Double_Gaussian, Multi_Gaussian
    from pyqtgraph_MATPLOTLIBWIDGET import MatplotlibWidget
except:
    from Spectrum_analysis.Spectra_fit_funcs import Spectra_Fit, Single_Gaussian, Single_Lorentzian, Double_Gaussian, Multi_Gaussian
    from Spectrum_analysis.pyqtgraph_MATPLOTLIBWIDGET import MatplotlibWidget

matplotlib.use('Qt5Agg')

"""Recylce params for plotting"""
plt.rc('xtick', labelsize = 10)
plt.rc('xtick.major', pad = 1)
plt.rc('ytick', labelsize = 10)
plt.rc('lines', lw = 1.5, markersize = 3.5)
plt.rc('legend', fontsize = 10)
plt.rc('axes', linewidth=1.5)

pg.mkQApp()
pg.setConfigOption('background', 'w')


base_path = Path(__file__).parent
file_path = (base_path / "Spectra_plot_fit_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

def updateDelay(scale, time):
    """ Hack fix for scalebar inaccuracy"""
    QtCore.QTimer.singleShot(time, scale.updateBar)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        pg.setConfigOption('imageAxisOrder', 'row-major')
        super(TemplateBaseClass, self).__init__()
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        # self.ui.fitFunc_comboBox.addItems(["Single Gaussian","Single Lorentzian", "Double Gaussian", "Multiple Gaussians"])
#        self.ui.actionExit.triggered.connect(self.close_application)
        
        ##setup ui signals
        self.ui.importSpec_pushButton.clicked.connect(self.open_file)
        self.ui.importBck_pushButton.clicked.connect(self.open_bck_file)
        self.ui.importWLRef_pushButton.clicked.connect(self.open_wlref_file)

        self.ui.load_spectra_scan_pushButton.clicked.connect(self.open_spectra_scan_file)
        self.ui.load_bck_file_pushButton.clicked.connect(self.open_spectra_bck_file)
        self.ui.load_fitted_scan_pushButton.clicked.connect(self.open_fit_scan_file)
        
        self.ui.plot_pushButton.clicked.connect(self.plot)
        self.ui.plot_fit_scan_pushButton.clicked.connect(self.plot_fit_scan)
        self.ui.plot_raw_scan_pushButton.clicked.connect(self.plot_raw_scan)
        self.ui.plot_intensity_sums_pushButton.clicked.connect(self.plot_intensity_sums)
        
        self.ui.fit_pushButton.clicked.connect(self.fit_and_plot)
        self.ui.fit_scan_pushButton.clicked.connect(self.fit_and_plot_scan)
        # self.ui.config_fit_params_pushButton.clicked.connect(self.configure_fit_params)
        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        self.ui.add_to_mem_pushButton.clicked.connect(self.add_trace_to_mem)
        self.ui.export_single_figure_pushButton.clicked.connect(self.export_plot_window)
        self.ui.export_scan_figure_pushButton.clicked.connect(self.export_window)
        self.ui.analyze_spectra_fits_pushButton.clicked.connect(self.analyze_spectra_fits)

        self.ui.import_pkl_pushButton.clicked.connect(self.open_pkl_file)
        self.ui.data_txt_pushButton.clicked.connect(self.pkl_data_to_txt)
        self.ui.scan_params_txt_pushButton.clicked.connect(self.pkl_params_to_txt)

        self.ui.pkl_to_h5_pushButton.clicked.connect(self.pkl_to_h5)

        self.ui.tabWidget.currentChanged.connect(self.switch_overall_tab)
        self.ui.fitFunc_comboBox.currentTextChanged.connect(self.switch_bounds_and_guess_tab)
        self.ui.adjust_param_checkBox.stateChanged.connect(self.switch_adjust_param)

        self.ui.export_data_pushButton.clicked.connect(self.export_data)
        self.ui.clear_export_data_pushButton.clicked.connect(self.clear_export_data)

        # for i in reversed(range(self.ui.bounds_groupBox.layout().count())):
        #     self.ui.bounds_groupBox.layout().itemAt(i).widget().deleteLater()
        #self.ui.single_bounds_page.layout().addWidget(QtWidgets.QPushButton("test"))

        self.ui.plot = pg.PlotWidget()
        self.ui.plot.setTitle(title="Spectrum Plot")
        #self.ui.plot.show()
        
        self.file = None
        self.bck_file = None
        self.wlref_file = None
        self.x = None
        self.y = None
        self.out = None # output file after fitting
        
        # Peak parameters if adjust params is selected
        self.center_min = None
        self.center_max = None

        #variables accounting for data received from FLIM analysis
        self.opened_from_flim = False #switched to True in FLIM_plot when "analyze lifetime" clicked
        self.sum_data_from_flim = []
        
        # fit scan file variable set to None for analyze_spectra_fits
        self.fit_scan_file = None

        #container for data to append to txt file
        self.data_list = []

        #for adding traces to memory for plotting/exporting all at once
        self.x_mem = []
        self.y_mem = []
        self.best_fit_mem = []
        self.legend = []
        self.single_spec_fit_called = False
        
        self.show()
    
    """ Open Single Spectrum files """
    def open_file(self):
        try:
            self.single_spec_filename = QtWidgets.QFileDialog.getOpenFileName(self)
            try:
                try:
                    self.file = np.loadtxt(self.single_spec_filename[0], skiprows = 1, delimiter=" ")
                except:        
                    self.file = np.loadtxt(self.single_spec_filename[0], skiprows = 16, delimiter='\t')
            except:
                self.file = np.genfromtxt(self.single_spec_filename[0], skip_header=1, skip_footer=3, delimiter='\t')
            self.opened_from_flim = False
        except:
            pass
    
    def open_bck_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            try:
                try:
                    self.bck_file = np.loadtxt(filename[0], skiprows=1, delimiter=" ")
                except:
                     self.bck_file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
            except:
                self.bck_file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')    
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass
    
    def open_wlref_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            try:
                try:
                    self.wlref_file = np.loadtxt(filename[0], skiprows=1, delimiter= " ")
                except:
                    self.wlref_file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
            except:
                self.wlref_file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')
        except:
            pass
        
    """Open Scan Files"""
    def open_spectra_scan_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self, filter="Scan files (*.pkl *.h5)")
            self.filename_for_viewer_launch = filename[0]
            if ".pkl" in filename[0]:
                self.spec_scan_file = pickle.load(open(filename[0], 'rb'))
                if self.ui.launch_data_viewer_checkBox.isChecked():
                    self.launch_h5_pkl_viewer()
                self.scan_file_type = "pkl"
            elif ".h5" in filename[0]:
                self.spec_scan_file = h5py.File(filename[0], 'r')
                if self.ui.launch_data_viewer_checkBox.isChecked():
                    self.launch_h5_pkl_viewer()
                self.scan_file_type = "h5"
            self.get_data_params()
            self.ui.result_textBrowser2.append("Done Loading - Spectra Scan File")
        except Exception as e:
            self.ui.result_textBrowser2.append(str(e))
            pass
    
    def open_spectra_bck_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            self.bck_file = np.loadtxt(filename[0])#, skiprows=1, delimiter=None)
            self.ui.result_textBrowser2.append("Done Loading - Background File")
        except Exception as e:
            self.ui.result_textBrowser2.append(str(e))
            pass
       
    def open_fit_scan_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            with pg.BusyCursor():
                self.filename_for_viewer_launch = filename[0]
                self.fit_scan_file = pickle.load(open(filename[0], 'rb'))
                if self.ui.launch_data_viewer_checkBox.isChecked():
                    self.launch_h5_pkl_viewer() # TODO Needs to implement reading the fit result datatype in PKL Viewer
                self.ui.result_textBrowser2.append("Done Loading - Scan Fit File")
        except Exception as e:
            self.ui.result_textBrowser2.append(str(e))
            pass

    def open_pkl_file(self):
        """ Open pkl file to convert to txt """
        try:
            self.pkl_to_convert = QtWidgets.QFileDialog.getOpenFileNames(self)
            files = self.pkl_to_convert[0]
            for i in range(len(files)):
                self.filename_for_viewer_launch = files[i]
                if self.ui.launch_data_viewer_checkBox_2.isChecked():
                    self.launch_h5_pkl_viewer()
        except:
            pass
    
    def launch_h5_pkl_viewer(self):
        """ Launches H5/PKL viewer to give an insight into the data and its structure"""
        viewer_window = h5_pkl_view.H5PklView(sys.argv)
        viewer_window.settings['data_filename'] = self.filename_for_viewer_launch
    
    def analyze_spectra_fits(self):
        """Analyze fits to the individual spectrum within a spectra scan fit file"""
        if self.fit_scan_file is None:
            self.open_fit_scan_file()
        
        #result_no = int(self.ui.result_spinBox.value())
        #self.matplotlibwidget = MatplotlibWidget(size=(12,8), dpi=300)
        #self.fit_scan_file['result_'+str(0)].plot(fig=self.matplotlibwidget.getFigure().add_subplot(111))
        #self.matplotlibwidget.draw()
        #self.matplotlibwidget.show()
        analyze_window = Analyze(scan_fit_file=self.fit_scan_file)
        analyze_window.run()
        
    def switch_overall_tab(self):
        """ Enable/disable fit settings on right depending on current tab """
        if self.ui.tabWidget.currentIndex() == 0:
            self.ui.fitting_settings_groupBox.setEnabled(True)
            self.ui.fit_pushButton.setEnabled(True)
            self.ui.fit_scan_pushButton.setEnabled(False)
            self.ui.save_all_checkBox.setEnabled(False)
            self.ui.scan_fit_settings_groupBox.setEnabled(False)
        elif self.ui.tabWidget.currentIndex() == 1:
            self.ui.fitting_settings_groupBox.setEnabled(False)
            self.ui.fit_pushButton.setEnabled(False)
            self.ui.fit_scan_pushButton.setEnabled(True)
            self.ui.save_all_checkBox.setEnabled(True)
            self.ui.scan_fit_settings_groupBox.setEnabled(True)
        elif self.ui.tabWidget.currentIndex() == 2:
            self.ui.fitting_settings_groupBox.setEnabled(False)
            self.ui.fit_pushButton.setEnabled(False)
            self.ui.fit_scan_pushButton.setEnabled(False)
            self.ui.save_all_checkBox.setEnabled(False)
            self.ui.scan_fit_settings_groupBox.setEnabled(False)

    """ Single spectrum functions """
    def switch_bounds_and_guess_tab(self):
        """ Show the appropriate bounds and initial guess params based on fit function """
        fit_func = self.ui.fitFunc_comboBox.currentText()
        if fit_func == "Single Gaussian" or fit_func == "Single Lorentzian":
            self.ui.n_label.setEnabled(False)
            self.ui.n_spinBox.setEnabled(False)
            self.ui.bounds_stackedWidget.setCurrentIndex(0)
            self.ui.guess_stackedWidget.setCurrentIndex(0)
            self.ui.plot_components_checkBox.setEnabled(False)
            self.ui.n_spinBox.setValue(1)
        elif fit_func == "Double Gaussian":
            self.ui.n_label.setEnabled(False)
            self.ui.n_spinBox.setEnabled(False)
            self.ui.bounds_stackedWidget.setCurrentIndex(1)
            self.ui.guess_stackedWidget.setCurrentIndex(1)
            self.ui.plot_components_checkBox.setEnabled(True)
            self.ui.n_spinBox.setValue(2)
        elif fit_func == "Triple Gaussian":
            self.ui.n_label.setEnabled(False)
            self.ui.n_spinBox.setEnabled(False)
            self.ui.bounds_stackedWidget.setCurrentIndex(2)
            self.ui.guess_stackedWidget.setCurrentIndex(2)
            self.ui.plot_components_checkBox.setEnabled(True)
            self.ui.n_spinBox.setValue(3)

    def switch_adjust_param(self):
        """ Enable bounds and initial guess only when adjust parameters is checked """
        checked = self.ui.adjust_param_checkBox.isChecked()
        self.ui.bounds_groupBox.setEnabled(checked)
        self.ui.guess_groupBox.setEnabled(checked)

    def check_loaded_files(self):
        """ 
        Check if 'subtract background' or 'white light correction' is checked 
        and if required files have been loaded. 
        """
        if self.ui.subtract_bck_radioButton.isChecked() and self.bck_file is None:
            self.ui.result_textBrowser.setText("You need to load a background file.")
        elif self.wlref_file is not None and self.ui.WLRef_checkBox.isChecked() == False:
            self.ui.result_textBrowser.setText("You need to check the White Light Correction option!")
        elif self.wlref_file is None and self.ui.WLRef_checkBox.isChecked():
            self.ui.result_textBrowser.setText("You need to load a White Light Ref file.")
        else:
            return True

    def plot(self):
        try:
            if self.opened_from_flim:
                flim_data = self.sum_data_from_flim.T
                interp = interpolate.interp1d(flim_data[:,0], flim_data[:,1])
                x_range = [flim_data[:,0][0], flim_data[:,0][-1]]
                xnew = np.linspace(x_range[0], x_range[1], 100 )
                ynew = interp(xnew)
                self.file = np.zeros((xnew.shape[0], 2))
                self.file[:,0] = xnew
                self.file[:,1] = ynew
                self.x = xnew
                self.y = ynew

            elif self.file is None: #elif
                self.ui.result_textBrowser.setText("You need to load a data file.")    
            else:
                self.x = self.file[:,0]
                self.y = self.file[:,1]
                
            self.check_loaded_files()

            if self.check_loaded_files() == True: #check the following conditions if all required files have been provided
                if self.ui.subtract_bck_radioButton.isChecked() == True and self.ui.WLRef_checkBox.isChecked() == False:
                    bck_y = self.bck_file[:,1]
                    self.y = self.y - bck_y
                elif self.ui.subtract_bck_radioButton.isChecked() == False and self.ui.WLRef_checkBox.isChecked() == True:
                    wlref_y = self.wlref_file[:,1]
                    self.y = (self.y)/wlref_y
                
                elif self.ui.subtract_bck_radioButton.isChecked() == True and self.ui.WLRef_checkBox.isChecked() == True:
                    bck_y = self.bck_file[:,1]
                    wlref_y = self.wlref_file[:,1]
                    self.y = (self.y-bck_y)/wlref_y
            
            
            if self.ui.norm_checkBox.isChecked():
                self.normalize()
            
            self.check_eV_state()
            self.ui.plot.show()
            self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
            
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        if self.ui.fit_in_eV.isChecked():
            self.ui.plot.setLabel('bottom', 'Energy (eV)')
        else:    
            self.ui.plot.setLabel('bottom', 'Wavelength (nm)')

        self.single_spec_fit_called = False
    
    def normalize(self):
        self.y = (self.y) / np.amax(self.y)
    
    def check_eV_state(self):
        if self.ui.fit_in_eV.isChecked():
            self.x = np.sort(1240/self.file[:,0])
            self.y = [self.y[i] for i in np.argsort(1240/self.file[:,0])]
        else:
            self.x = self.file[:,0]
            self.y = self.file[:,1]
    
    def clear_plot(self):
        self.ui.plot.clear()
        self.ui.result_textBrowser.clear()
        
    def clear_check(self):
        if self.ui.clear_checkBox.isChecked() == True:
            return True
        elif self.ui.clear_checkBox.isChecked() == False:
            return False

    def fit_and_plot(self):
        fit_func = self.ui.fitFunc_comboBox.currentText()
        
        try:
            self.plot()
            if self.opened_from_flim:
                self.file = np.zeros((self.x.shape[0], 2))
                self.file[:,0] = self.x
                self.file[:,1] = self.y
                bck = lmfit.models.LinearModel(prefix='line_')
                gmodel = GaussianModel(prefix='g1_')
                pars = bck.make_params(intercept=self.y.min(), slope=0)
                pars += gmodel.guess(self.y, x=self.x)
                comp_model = gmodel + bck
                self.result = comp_model.fit(self.y, pars, x=self.x, nan_policy='propagate')
                self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
                self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                self.ui.result_textBrowser.setText(self.result.fit_report())
                

            if self.ui.plot_without_bck_radioButton.isChecked(): #if plot w/o bck, create dummy bck_file
                self.bck_file = np.zeros(shape=(self.file.shape[0], 2))
                self.bck_file[:,0] = self.file[:,0]
 
            # if self.ui.subtract_bck_radioButton.isChecked() == False:
            #     self.ui.result_textBrowser.setText("You need to check the subtract background option!")
            if self.check_loaded_files is None:
                pass
            elif self.opened_from_flim:
                pass
            else:
                self.check_eV_state()
                if fit_func == "Single Gaussian":
                    single_gauss = Single_Gaussian(self.file, self.bck_file, wlref=self.wlref_file, fit_in_eV=self.ui.fit_in_eV.isChecked())
                    if self.ui.adjust_param_checkBox.isChecked():
                        center1_min = self.ui.single_peakcenter1_min_spinBox.value()
                        center1_max = self.ui.single_peakcenter1_max_spinBox.value()
                        center1_guess = self.ui.single_peakcenter1_guess_spinBox.value()
                        sigma1_guess = self.ui.single_sigma1_guess_spinBox.value()
                        self.result = single_gauss.gaussian_model_w_lims(center1_guess, sigma1_guess,
                            [center1_min, center1_max])
                    else:
                        self.result = single_gauss.gaussian_model()
                    #self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
                    self.plot()
                    self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                    self.ui.result_textBrowser.setText(self.result.fit_report())
                
                elif fit_func == "Single Lorentzian": #and self.ui.subtract_bck_radioButton.isChecked() == True:
                    single_lorentzian = Single_Lorentzian(self.file, self.bck_file, wlref=self.wlref_file, fit_in_eV=self.ui.fit_in_eV.isChecked())
                    
                    if self.ui.adjust_param_checkBox.isChecked():
                        center1_min = self.ui.single_peakcenter1_min_spinBox.value()
                        center1_max = self.ui.single_peakcenter1_max_spinBox.value()
                        center1_guess = self.ui.single_peakcenter1_guess_spinBox.value()
                        sigma1_guess = self.ui.single_sigma1_guess_spinBox.value()
                        self.result = single_lorentzian.lorentzian_model_w_lims(center1_guess, sigma1_guess,
                                [center1_min, center1_max])
                    else:
                        self.result = single_lorentzian.lorentzian_model()
                    #self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
                    self.plot()
                    self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                    self.ui.result_textBrowser.setText(self.result.fit_report())
                
                elif fit_func == "Double Gaussian": #and self.ui.subtract_bck_radioButton.isChecked() == True:
                    double_gauss = Double_Gaussian(self.file, self.bck_file, wlref=self.wlref_file, fit_in_eV=self.ui.fit_in_eV.isChecked())
                    if self.ui.adjust_param_checkBox.isChecked():
                        center1_min = self.ui.double_peakcenter1_min_spinBox.value()
                        center1_max = self.ui.double_peakcenter1_max_spinBox.value()
                        center2_min = self.ui.double_peakcenter2_min_spinBox.value()
                        center2_max = self.ui.double_peakcenter2_max_spinBox.value()
                        center1_guess = self.ui.double_peakcenter1_guess_spinBox.value()
                        sigma1_guess = self.ui.double_sigma1_guess_spinBox.value()
                        center2_guess = self.ui.double_peakcenter2_guess_spinBox.value()
                        sigma2_guess = self.ui.double_sigma2_guess_spinBox.value()

                        peak_pos = [center1_guess, center2_guess]
                        sigma = [sigma1_guess, sigma2_guess]
                        min_max_range = [ [center1_min, center1_max], [center2_min, center2_max] ]
                        self.result = double_gauss.gaussian_model_w_lims(peak_pos, sigma, min_max_range)

                    else:
                        self.result = double_gauss.gaussian_model()

                    #self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
                    self.plot()
                    self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                    if self.ui.plot_components_checkBox.isChecked():
                        comps = self.result.eval_components(x=self.x)
                        self.ui.plot.plot(self.x, comps['g1_'], pen='b', clear=False)
                        self.ui.plot.plot(self.x, comps['g2_'], pen='g', clear=False)

                    self.ui.result_textBrowser.setText(self.result.fit_report())
                
                elif fit_func == "Triple Gaussian": #and self.ui.subtract_bck_radioButton.isChecked() == True:
                    #currently only works for triple gaussian (n=3)
                    multiple_gauss = Multi_Gaussian(self.file, self.bck_file, 3, wlref=self.wlref_file, fit_in_eV=self.ui.fit_in_eV.isChecked())
                    if self.ui.adjust_param_checkBox.isChecked():
                        center1_min = self.ui.multi_peakcenter1_min_spinBox.value()
                        center1_max = self.ui.multi_peakcenter1_max_spinBox.value()
                        center2_min = self.ui.multi_peakcenter2_min_spinBox.value()
                        center2_max = self.ui.multi_peakcenter2_max_spinBox.value()
                        center3_min = self.ui.multi_peakcenter3_min_spinBox.value()
                        center3_max = self.ui.multi_peakcenter3_max_spinBox.value()
                        center1_guess = self.ui.multi_peakcenter1_guess_spinBox.value()
                        sigma1_guess = self.ui.multi_sigma1_guess_spinBox.value()
                        center2_guess = self.ui.multi_peakcenter2_guess_spinBox.value()
                        sigma2_guess = self.ui.multi_sigma2_guess_spinBox.value()
                        center3_guess = self.ui.multi_peakcenter3_guess_spinBox.value()
                        sigma3_guess = self.ui.multi_sigma3_guess_spinBox.value()
#                        num_gaussians = 3
                        peak_pos = [center1_guess, center2_guess, center3_guess]
                        sigma = [sigma1_guess, sigma2_guess, sigma3_guess]
                        min_max_range = [ [center1_min, center1_max], [center2_min, center2_max], [center3_min, center3_max] ]
                        
                        self.result = multiple_gauss.gaussian_model_w_lims(peak_pos, sigma, min_max_range)
                    else:
                        self.result = multiple_gauss.gaussian_model()

                    #self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
                    self.plot()
                    self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                    if self.ui.plot_components_checkBox.isChecked():
                        comps = self.result.eval_components(x=self.x)
                        self.ui.plot.plot(self.x, comps['g1_'], pen='b', clear=False)
                        self.ui.plot.plot(self.x, comps['g2_'], pen='g', clear=False)
                        self.ui.plot.plot(self.x, comps['g3_'], pen='c', clear=False)
                    self.ui.result_textBrowser.setText(self.result.fit_report())

                self.data_list.append(self.ui.result_textBrowser.toPlainText())
                self.single_spec_fit_called = True
        
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
    
    def export_window(self):
        self.export_window = ExportImages()
        self.export_window.export_fig_signal.connect(self.pub_ready_plot_export)
    
    def export_plot_window(self):
        self.exportplotwindow = ExportPlotWindow()
        self.exportplotwindow.export_fig_signal.connect(self.pub_ready_plot_export)

    def pub_ready_plot_export(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
        """Recylce params for plotting"""
        plt.rc('xtick', labelsize = 20)
        plt.rc('xtick.major', pad = 3)
        plt.rc('ytick', labelsize = 20)
        plt.rc('lines', lw = 1.5, markersize = 7.5)
        plt.rc('legend', fontsize = 20)
        plt.rc('axes', linewidth=3.5)
        try:
            try:
                try:
                    data = self.spec_scan_file
                except:
                    data = self.fit_scan_file
                if self.export_window.ui.reverse_checkBox.isChecked():
                    colormap = str(self.export_window.ui.cmap_comboBox.currentText())+"_r"
                else:
                    colormap = str(self.export_window.ui.cmap_comboBox.currentText())
                if str(self.export_window.ui.dataChannel_comboBox.currentText()) == "Fitted":
                    param_selection = str(self.ui.comboBox.currentText())
                    if param_selection == 'pk_pos': label = 'PL Peak Position (n.m.)'
                    elif param_selection == 'fwhm': label = 'PL FWHM (n.m.)'
                    cpm.plot_confocal(self.img, FLIM_adjust = False, stepsize = data['Scan Parameters']['X step size (um)'], cmap=colormap, cbar_label=label,
                                      vmin=self.export_window.ui.vmin_spinBox.value(), vmax=self.export_window.ui.vmax_spinBox.value())
                elif str(self.export_window.ui.dataChannel_comboBox.currentText()) == "Raw":
                    cpm.plot_confocal(self.sums, FLIM_adjust = False, figsize=(10,10), stepsize = data['Scan Parameters']['X step size (um)'], cmap=colormap, 
                                      vmin=self.export_window.ui.vmin_spinBox.value(), vmax=self.export_window.ui.vmax_spinBox.value())
                plt.tick_params(direction='out', length=8, width=3.5)
                plt.tight_layout()
                plt.savefig(filename[0],bbox_inches='tight', dpi=300)
                plt.close()
            except:
                if self.x_mem == []:
                    self.ui.result_textBrowser.setText("Add traces to memory first!")
                else:
                    plt.figure(figsize=(8,6))
                    plt.tick_params(direction='out', length=8, width=3.5)
                    for i in range(len(self.x_mem)):
                        plt.plot(self.x_mem[i], self.y_mem[i], label=str(self.legend[i]))
                        if self.single_spec_fit_called == True:
                            plt.plot(self.x_mem[i], self.best_fit_mem[i],'k')
                    
                    if self.ui.fit_in_eV.isChecked():
                        plt.xlabel("Energy (eV)", fontsize=20, fontweight='bold')
                    else:
                        plt.xlabel("Wavelength (nm)", fontsize=20, fontweight='bold')
                    plt.ylabel("Intensity (a.u.)", fontsize=20, fontweight='bold')
                    plt.xlim([self.exportplotwindow.ui.lowerX_spinBox.value(),self.exportplotwindow.ui.upperX_spinBox.value()])
                    plt.ylim([self.exportplotwindow.ui.lowerY_spinBox.value(),self.exportplotwindow.ui.upperY_doubleSpinBox.value()])
                    plt.legend()
                    plt.tight_layout()
                    
                    plt.savefig(filename[0],bbox_inches='tight', dpi=300)
                    plt.close()
            
        except AttributeError:
            self.ui.result_textBrowser.setText("Need to fit the data first!")

    def export_data(self):
        """ Save fit params and srv calculations stored in data_list as .txt """
        folder = os.path.dirname(self.single_spec_filename[0])
        filename_ext = os.path.basename(self.single_spec_filename[0])
        filename = os.path.splitext(filename_ext)[0] #get filename without extension

        path = folder + "/" + filename + "_fit_results.txt"
        if not os.path.exists(path):
            file = open(path, "w+")
        else:
            file = open(path, "a+")

        for i in range(len(self.data_list)):
            file.write(self.data_list[i] + "\n\n")

        self.data_list = []
        file.close()

    def clear_export_data(self):
        self.data_list = []
        self.x_mem = []
        self.y_mem = []
        self.legend = []
        self.best_fit_mem = []

    def add_trace_to_mem(self):
        try:
            self.x_mem.append(self.x)
            self.y_mem.append(self.y)
            if self.single_spec_fit_called == True:
                self.best_fit_mem.append(self.result.best_fit)
            self.legend.append(self.ui.lineEdit.text())
        except Exception as e:
            print(e)


    """ Scan spectra functions """
    def get_data_params(self):
        data = self.spec_scan_file
        if self.scan_file_type == "pkl":
            self.intensities = data['Intensities']
            self.wavelengths = data['Wavelengths']
            # try:
            self.x_scan_size = data['Scan Parameters']['X scan size (um)']
            self.y_scan_size = data['Scan Parameters']['Y scan size (um)']
            self.x_step_size = data['Scan Parameters']['X step size (um)']
            self.y_step_size = data['Scan Parameters']['Y step size (um)']
            # except: # TODO test and debug loading pkl file w/o scan parameters
            #     self.configure_scan_params()
            #     while not hasattr(self, "scan_params_entered"):
            #         pass
            #     self.x_scan_size = self.param_window.ui.x_scan_size_spinBox.value()
            #     self.y_scan_size = self.param_window.ui.y_scan_size_spinBox.value()
            #     self.x_step_size = self.param_window.ui.x_step_size_spinBox.value()
            #     self.y_step_size = self.param_window.ui.y_step_size_spinBox.value()

        else: #run this if scan file is h5
            self.x_scan_size = data['Scan Parameters'].attrs['X scan size (um)']
            self.y_scan_size = data['Scan Parameters'].attrs['Y scan size (um)']
            self.x_step_size = data['Scan Parameters'].attrs['X step size (um)']
            self.y_step_size = data['Scan Parameters'].attrs['Y step size (um)']
            self.intensities = data['Intensities'][()] #get dataset values
            self.wavelengths = data['Wavelengths'][()]

        self.numb_x_pixels = int(np.ceil(self.x_scan_size/self.x_step_size))
        self.numb_y_pixels = int(np.ceil(self.y_scan_size/self.y_step_size))

        """Open param window and get peak center range values and assign it to variables to use later"""
    # def configure_scan_params(self):
    #      self.param_window = ParamWindow()
    #     self.param_window.peak_range.connect(self.peak_range)
    
    # def peak_range(self, peaks):
    #     self.center_min = peaks[0]
    #     self.center_max = peaks[1]

    def plot_fit_scan(self):
        try:
            if self.ui.use_raw_scan_settings.isChecked():
                num_x = self.numb_x_pixels
                num_y  =self.numb_y_pixels
            else:
                num_x = self.ui.num_x_spinBox.value()
                num_y = self.ui.num_y_spinBox.value()
            
            numb_of_points = num_x * num_y #75*75
            
            fwhm = np.zeros(shape=(numb_of_points,1))
            pk_pos = np.zeros(shape=(numb_of_points,1))
#            sigma = np.zeros(shape=(numb_of_points,1))
#            height = np.zeros(shape=(numb_of_points,1))
            
            if type(self.fit_scan_file['result_0']) == dict:
                for i in range(numb_of_points):
                    fwhm[i, 0] = 2.3548200*self.fit_scan_file['result_'+str(i)]['g1_sigma']
                    pk_pos[i, 0] = self.fit_scan_file['result_'+str(i)]['g1_center']
#                    sigma[i, 0] = self.fit_scan_file['result_'+str(i)]['g1_sigma']
#                    height[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_height']
            
            elif type(self.fit_scan_file['result_0']) == lmfit.model.ModelResult:
                for i in range(numb_of_points):
                    fwhm[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_fwhm']
                    pk_pos[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_center']
#                    sigma[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_sigma']
#                    height[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_height']
            
            newshape = (num_x, num_y)
            
            param_selection = str(self.ui.comboBox.currentText())
            self.img = np.reshape(eval(param_selection), newshape)

            if num_y == 1:
                x = np.linspace(0, self.x_scan_size, num_x)
                self.graph_layout=pg.GraphicsLayoutWidget()
                self.plot = self.graph_layout.addPlot(title="Line Scan")
                self.plot.plot(x, self.img[:,0], pen="r")
                self.graph_layout.show()
            elif num_x == 1:
                y = np.linspace(0, self.y_scan_size, num_y)
                self.graph_layout=pg.GraphicsLayoutWidget()
                self.plot = self.graph_layout.addPlot(title="Line Scan")
                self.plot.plot(y, self.img[0,:], pen="r")
                self.graph_layout.show()
            
            else:
                self.fit_scan_viewbox = pg.ImageView()
                if self.ui.use_raw_scan_settings.isChecked():
                    self.fit_scan_viewbox.setImage(self.img, scale=
                                                    (self.x_step_size,
                                                    self.y_step_size))
                    scale = pg.ScaleBar(size=2,suffix='um')
                    scale.setParentItem(self.fit_scan_viewbox.view)
                    scale.anchor((1, 1), (1, 1), offset=(-30, -30))
                    self.fit_scan_viewbox.view.sigRangeChanged.connect(lambda: updateDelay(scale, 10))
                else:
                    self.fit_scan_viewbox.setImage(self.img)
                
                self.fit_scan_viewbox.view.invertY(False)
                self.fit_scan_viewbox.show()

        except Exception as e:
            self.ui.result_textBrowser2.append(str(e))
            pass
            
    def plot_raw_scan(self):
        try:
            # TODO test line scan plots

            intensities = self.intensities.T #this is only there because of how we are saving the data in the app
            intensities = np.reshape(intensities, newshape=(2048,self.numb_x_pixels, self.numb_y_pixels))
            self.raw_scan_viewbox = pg.ImageView() 
            self.raw_scan_viewbox.setImage(intensities, scale=
                                                  (self.x_step_size,
                                                   self.y_step_size), xvals=self.wavelengths)
            
            #roi_plot = self.ui.raw_scan_viewBox.getRoiPlot()
            #roi_plot.plot(data['Wavelengths'], intensities)
            self.raw_scan_viewbox.view.invertY(False)
            scale = pg.ScaleBar(size=2,suffix='um')
            scale.setParentItem(self.raw_scan_viewbox.view)
            scale.anchor((1, 1), (1, 1), offset=(-30, -30))
            self.raw_scan_viewbox.view.sigRangeChanged.connect(lambda: updateDelay(scale, 10))
            self.raw_scan_viewbox.show()
            
        except Exception as e:
            self.ui.result_textBrowser2.append(str(e))

    def plot_intensity_sums(self):
        try:
            # TODO test line scan plots

            #intensities = np.reshape(intensities, newshape=(2048, numb_pixels_X*numb_pixels_Y))
            
            sums = np.sum(self.intensities, axis=-1)
            self.sums = np.reshape(sums, newshape=(self.numb_x_pixels, self.numb_y_pixels))
            self.intensity_sums_viewBox = pg.ImageView()
            
            self.intensity_sums_viewBox.setImage(self.sums, scale=
                                                  (self.x_step_size,
                                                   self.y_step_size))
            self.intensity_sums_viewBox.view.invertY(False)
            
            scale = pg.ScaleBar(size=2,suffix='um')
            scale.setParentItem(self.intensity_sums_viewBox.view)
            scale.anchor((1, 1), (1, 1), offset=(-30, -30))
            self.intensity_sums_viewBox.view.sigRangeChanged.connect(lambda: updateDelay(scale, 10))
            self.intensity_sums_viewBox.show()

        except Exception as e:
            self.ui.result_textBrowser2.append(str(e))

    
    def fit_and_plot_scan(self):
#        self.ui.result_textBrowser.append("Starting Scan Fitting")
        print("Starting Scan Fitting")
        print("Using Single Gaussian to Fit\nThis is the only fitting functions implemented")
        
        with pg.BusyCursor():
            try:
                """Define starting and stopping wavelength values here"""
                start_nm = int(self.ui.start_nm_spinBox.value())
                stop_nm = int(self.ui.stop_nm_spinBox.value())
                
                if self.bck_file is None:
                    print("Load Background file!")
                ref = self.bck_file
                index = (ref[:,0]>start_nm) & (ref[:,0]<stop_nm)
                
                x = self.wavelengths
                x = x[index]
                
                data_array = self.intensities
                
                result_dict = {}
                result_dict["Scan Parameters"] = self.spec_scan_file['Scan Parameters']
                result_dict["OceanOptics Parameters"] = self.spec_scan_file["OceanOptics Parameters"]
                
                for i in range(data_array.shape[0]):
                    
                    y = data_array[i, index] # intensity
                    yref = ref[index, 1]
                    
                    y = y - yref # background correction
                    y = y - np.mean(y[(x>start_nm) & (x<start_nm + 25)]) # removing any remaining bckgrnd
                    
                    gmodel = GaussianModel(prefix = 'g1_') # calling gaussian model
                    pars = gmodel.guess(y, x=x) # parameters - center, width, height
                    result = gmodel.fit(y, pars, x=x, nan_policy='propagate')
                    if self.ui.save_all_checkBox.isChecked():
                        result_dict["result_"+str(i)] = result
                    else:
                        result_dict["result_"+str(i)] = result.best_values
                                
    #            self.ui.result_textBrowser.append("Scan Fitting Complete!")
                print("Scan Fitting Complete!")
    
                filename = QtWidgets.QFileDialog.getSaveFileName(self)
                pickle.dump(result_dict, open(filename[0]+"_fit_result_dict.pkl", "wb"))
                
    #            self.ui.result_textBrowser.append("Data Saved!")
                print("Data Saved!")
            
            except Exception as e:
                self.ui.result_textBrowser2.append(str(e))
                pass
            
    #        self.ui.result_textBrowser.append("Loading Fit Data and Plotting")
            print("Loading Fit Data and Plotting")
            try:
                self.fit_scan_file = pickle.load(open(filename[0]+"_fit_result_dict.pkl", 'rb'))
                self.plot_fit_scan()
                
            except Exception as e:
                self.ui.result_textBrowser2.append(str(e))
                pass


    """ Pkl conversion functions """
    def pkl_data_to_txt(self):
        """ Get data from ocean optics scan pkl file, convert to txt"""
        for i in range(len(self.pkl_to_convert[0])):
            folder = os.path.dirname(self.pkl_to_convert[0][i])
            filename_ext = os.path.basename(self.pkl_to_convert[0][i])
            filename = os.path.splitext(filename_ext)[0] #get filename without extension
            pkl_file = pickle.load(open(self.pkl_to_convert[0][i], 'rb'))
    
            txt_file = np.zeros(shape=(2048,pkl_file['Intensities'].shape[0] + 1))
    
            data_array = pkl_file['Intensities']
            data_array = np.transpose(data_array)
            wavelength = pkl_file['Wavelengths']
    
            txt_file[:,0] = wavelength
    
            for i in range(pkl_file['Intensities'].shape[0]):
                txt_file[:,i+1] = data_array[:,i]
    
            np.savetxt(folder +"/"+ filename +"_data.txt", txt_file, fmt = '%.2f', delimiter= "\t", header="wavelength(nm), Intensities at different points")

    def pkl_params_to_txt(self):
        """ Get scan parameters from ocean optics scan pkl file, convert to txt """
        for i in range(len(self.pkl_to_convert[0])):
            folder = os.path.dirname(self.pkl_to_convert[0][i])
            filename_ext = os.path.basename(self.pkl_to_convert[0][i])
            filename = os.path.splitext(filename_ext)[0] #get filename without extension
            pkl_file = pickle.load(open(self.pkl_to_convert[0][i], 'rb'))
    
            pkl_scan = pkl_file['Scan Parameters']
            pkl_oo = pkl_file['OceanOptics Parameters']
            
            param_list = []
            param_list.append(['X scan start (um)', 'Y scan start (um)', 'X scan size (um)', 'Y scan size (um)',
                'X step size (um)', 'Y step size (um)', 'Integration Time (ms)', 'Scans to Average', 'Correct Dark Counts']) #list of param names
            param_list.append([ pkl_scan['X scan start (um)'], pkl_scan['Y scan start (um)'], pkl_scan['X scan size (um)'],
                pkl_scan['Y scan size (um)'], pkl_scan['X step size (um)'], pkl_scan['Y step size (um)'],
                pkl_oo['Integration Time (ms)'], pkl_oo['Scans Averages'], pkl_oo['Correct Dark Counts'] ]) #list of param values
    
            param_list = list(zip(*param_list)) #transpose so names and values are side-by-side
            save_to = folder +"/"+ filename +"_scan_parameters.txt"
            
            with open(save_to, 'w') as f:
                for item in param_list:
                    f.write("%s\t" % str(item[0])) #write name
                    f.write("%s\n" % str(item[1])) #write value

    def pkl_to_h5(self):
        """ Convert raw scan .pkl file to h5 """
        for i in range(len(self.pkl_to_convert[0])):                
            folder = os.path.dirname(self.pkl_to_convert[0][i])
            filename_ext = os.path.basename(self.pkl_to_convert[0][i])
            filename = os.path.splitext(filename_ext)[0] #get filename without extension
            pkl_file = pickle.load(open(self.pkl_to_convert[0][i], 'rb'))
    
            h5_filename = folder + "/" + filename + ".h5"
            h5_file = h5py.File(h5_filename, "w")
            self.traverse_dict_into_h5(pkl_file, h5_file)

    def traverse_dict_into_h5(self, dictionary, h5_output):
        """ 
        Create an h5 file using .pkl with scan data and params

        dictionary -- dictionary to convert
        h5_output -- h5 file or group to work in
        """
        for key in dictionary:
            if type(dictionary[key]) == dict:
                group = h5_output.create_group(key)
                previous_dict = dictionary[key]
                self.traverse_dict_into_h5(dictionary[key], group)
            else:
                if key == "Wavelengths" or key == "Intensities":
                    h5_output.create_dataset(key, data=dictionary[key])
                else:
                    h5_output.attrs[key] = dictionary[key]


    def close_application(self):
        choice = QtGui.QMessageBox.question(self, 'EXIT!',
                                            "Do you want to exit the app?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass
        
        
"""Analyze Window GUI and Functions"""
base_path = Path(__file__).parent
file_path = (base_path / "analyze_fit_results.ui").resolve()

uiFile = file_path

Analyze_WindowTemplate, Analyze_TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class Analyze(Analyze_TemplateBaseClass):  
    
    def __init__(self, scan_fit_file):
        pg.setConfigOption('imageAxisOrder', 'row-major')
        super(Analyze_TemplateBaseClass, self).__init__()
        
        # Create the main window
        self.ui = Analyze_WindowTemplate()
        self.ui.setupUi(self)
        
        self.ui.plot_pushButton.clicked.connect(self.plot)
        
        self.fit_scan_file = scan_fit_file
        self.show()
        
    def plot(self):
        matplotlib.use('Qt5Agg')
        try:
            result_no = int(self.ui.result_spinBox.value())
            if type(self.fit_scan_file['result_0']) == lmfit.model.ModelResult:
                self.matplotlibwidget = MatplotlibWidget(size=(12,8), dpi=300)
                self.fit_scan_file['result_'+str(result_no)].plot(fig=self.matplotlibwidget.getFigure().add_subplot(111))
                plt.tick_params(length=8, width=3)
                plt.xlabel("Wavelength (nm)", fontsize=25, fontweight='bold')
                plt.ylabel("Intensity (a.u.)", fontsize=25, fontweight='bold')
                plt.show()
                plt.tight_layout()
        except Exception as e:
            print(str(e))
            pass
    
    def run(self):
        win = Analyze()
        QtGui.QApplication.instance().exec_()
        return win

"""Export Images GUI"""
#ui_file_path = (base_path / "export_fig_gui.ui").resolve()
#export_WindowTemplate, export_TemplateBaseClass = pg.Qt.loadUiType(ui_file_path)

class ExportImages(ExportFigureWindow):
    
    def __init__(self):
        ExportFigureWindow.__init__(self)
        
        self.ui.dataChannel_comboBox.addItems(['Raw', 'Fitted'])

    
"""Run the Main Window"""    
def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win

#Uncomment below if you want to run this as standalone
#run()