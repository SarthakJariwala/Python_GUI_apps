# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

# system imports
import sys
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import matplotlib.pyplot as plt
import pickle
import time
from lmfit.models import GaussianModel
import customplotting.mscope as cpm
# local modules
try:
    from Spectra_fit_funcs import Spectra_Fit, Single_Gaussian, Single_Lorentzian
except:
    from Spectrum_analysis.Spectra_fit_funcs import Spectra_Fit, Single_Gaussian, Single_Lorentzian    


"""Recylce params for plotting"""
plt.rc('xtick', labelsize = 20)
plt.rc('xtick.major', pad = 3)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', lw = 1.5, markersize = 7.5)
plt.rc('legend', fontsize = 20)
plt.rc('axes', linewidth=3.5)

pg.mkQApp()
pg.setConfigOption('background', 'w')
pg.setConfigOption('imageAxisOrder', 'row-major')

base_path = Path(__file__).parent
file_path = (base_path / "Spectra_plot_fit_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        super(TemplateBaseClass, self).__init__()
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        self.ui.fitFunc_comboBox.addItems(["Single Gaussian","Single Lorentzian", "Double Gaussian", "Multiple Gaussians"])
        
#        self.ui.actionSave.triggered.connect(self.save_file)
#        self.ui.actionExit.triggered.connect(self.close_application)
        
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
        self.ui.config_fit_params_pushButton.clicked.connect(self.configure_fit_params)
        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        self.ui.export_fig_pushButton.clicked.connect(self.pub_ready_plot_export)

        self.ui.import_pkl_pushButton.clicked.connect(self.open_pkl_file)
        self.ui.data_txt_pushButton.clicked.connect(self.pkl_data_to_txt)
        self.ui.scan_params_txt_pushButton.clicked.connect(self.pkl_params_to_txt)

        
        self.file = None
        self.bck_file = None
        self.wlref_file = None
        self.x = None
        self.y = None
        self.out = None # output file after fitting
        
        # Peak parameters if adjust params is selected
        self.center_min = None
        self.center_max = None
        
        self.show()
    
    """Open Single Spectrum Files"""    
    def open_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            try:
                self.file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
            except:
                self.file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')
        except:
            pass
    
    def open_bck_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            try:
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
                self.wlref_file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
            except:
                self.wlref_file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')
        except:
            pass
        
    """Open Scan Files"""
    def open_spectra_scan_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            self.spec_scan_file = pickle.load(open(filename[0], 'rb'))
            self.ui.result_textBrowser.append("Done Loading - Spectra Scan File")
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass
    
    def open_spectra_bck_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            self.bck_file = np.loadtxt(filename[0])#, skiprows=1, delimiter=None)
            self.ui.result_textBrowser.append("Done Loading - Background File")
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass
       
    def open_fit_scan_file(self):
        try:
            filename = QtWidgets.QFileDialog.getOpenFileName(self)
            self.fit_scan_file = pickle.load(open(filename[0], 'rb'))
            self.ui.result_textBrowser.append("Done Loading - Scan Fit File")
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass

    def open_pkl_file(self):
        try:
            self.pkl_to_convert = QtWidgets.QFileDialog.getOpenFileName(self)
            self.ui.result_textBrowser.append("Done Loading - .pkl to convert")
        except:
            pass
    
    def save_file(self):# not used yet!
        try:
            filename = QtWidgets.QFileDialog.getSaveFileName(self)
            np.savetxt(filename[0], self.out, fmt = '%.5f', header = 'Time, Raw_PL, Sim_PL', delimiter = ' ')
        except:
            pass

    def plot(self):
        try:
            self.x = self.file[:,0]
            self.y = self.file[:,1]
            
            if self.ui.subtract_bck_checkBox.isChecked() == True and self.ui.WLRef_checkBox.isChecked() == False:
                bck_y = self.bck_file[:,1]
                self.y = self.y - bck_y
            
            elif self.ui.subtract_bck_checkBox.isChecked() == False and self.ui.WLRef_checkBox.isChecked() == True:
                wlref_y = self.wlref_file[:,1]
                self.y = (self.y)/wlref_y
            
            elif self.ui.subtract_bck_checkBox.isChecked() == True and self.ui.WLRef_checkBox.isChecked() == True:
                bck_y = self.bck_file[:,1]
                wlref_y = self.wlref_file[:,1]
                self.y = (self.y-bck_y)/wlref_y
            
            
            if self.ui.norm_checkBox.isChecked():
                self.normalize()
                
            self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
            
        except:
            pass
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Wavelength (nm)')
    
    def plot_fit_scan(self):
        try:
            if self.ui.use_raw_scan_settings.isChecked():
                data = self.spec_scan_file
                num_x = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
                num_y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
            else:
                num_x = self.ui.num_x_spinBox.value()
                num_y = self.ui.num_y_spinBox.value()
            
            numb_of_points = num_x * num_y #75*75
            
            fwhm = np.zeros(shape=(numb_of_points,1))
            pk_pos = np.zeros(shape=(numb_of_points,1))
#            pk_pos_plus = np.zeros(shape=(numb_of_points,1))
#            pk_pos_minus = np.zeros(shape=(numb_of_points,1))
            sigma = np.zeros(shape=(numb_of_points,1))
            height = np.zeros(shape=(numb_of_points,1))
            
            for i in range(numb_of_points):
                fwhm[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_fwhm']
                pk_pos[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_center']
                sigma[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_sigma']
                height[i, 0] = self.fit_scan_file['result_'+str(i)].values['g1_height']
            
            newshape = (num_x, num_y)
            
            param_selection = str(self.ui.comboBox.currentText())
            self.img = np.reshape(eval(param_selection), newshape)

            if self.ui.use_raw_scan_settings.isChecked():
                self.ui.fit_scan_viewbox.setImage(self.img, scale=
                                                  (data['Scan Parameters']['X step size (um)'],
                                                   data['Scan Parameters']['Y step size (um)']))
                scale = pg.ScaleBar(size=2,suffix='um')
                scale.setParentItem(self.ui.fit_scan_viewbox.view)
                scale.anchor((1, 1), (1, 1), offset=(-30, -30))
            else:
                self.ui.fit_scan_viewbox.setImage(self.img)
            
            self.ui.fit_scan_viewbox.view.invertY(False)
                
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass
            
    def plot_raw_scan(self):
        try:
            data = self.spec_scan_file
            numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
            numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
            # TODO test line scan plots

            intensities = data['Intensities'].T #this is only there because of how we are saving the data in the app
            
            intensities = np.reshape(intensities, newshape=(2048,numb_pixels_X,numb_pixels_Y))
            
            wavelengths = data['Wavelengths']
            
            self.ui.raw_scan_viewbox.view.invertY(False)
            self.ui.raw_scan_viewbox.setImage(intensities, scale=
                                                  (data['Scan Parameters']['X step size (um)'],
                                                   data['Scan Parameters']['Y step size (um)']), xvals=wavelengths)
            

            #roi_plot = self.ui.raw_scan_viewBox.getRoiPlot()
            #roi_plot.plot(data['Wavelengths'], intensities)
            scale = pg.ScaleBar(size=2,suffix='um')
            scale.setParentItem(self.ui.raw_scan_viewbox.view)
            scale.anchor((1, 1), (1, 1), offset=(-30, -30))
            
        except:
            pass

    def plot_intensity_sums(self):
        try:
            data = self.spec_scan_file
            numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
            numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
            # TODO test line scan plots

            intensities = data['Intensities']

            #intensities = np.reshape(intensities, newshape=(2048, numb_pixels_X*numb_pixels_Y))
            
            sums = np.sum(intensities, axis=-1)
            sums = np.reshape(sums, newshape=(numb_pixels_X, numb_pixels_Y))
            
            self.ui.intensity_sums_viewBox.setImage(sums, scale=
                                                  (data['Scan Parameters']['X step size (um)'],
                                                   data['Scan Parameters']['Y step size (um)']))
            self.ui.intensity_sums_viewBox.view.invertY(False)
            
            scale = pg.ScaleBar(size=2,suffix='um')
            scale.setParentItem(self.ui.intensity_sums_viewBox.view)
            scale.anchor((1, 1), (1, 1), offset=(-30, -30))

        except:
            pass

    def normalize(self):
        self.y = (self.y) / np.amax(self.y)
    
    def clear_plot(self):
        self.ui.plot.clear()
        self.ui.result_textBrowser.clear()
        
    def clear_check(self):
        if self.ui.clear_checkBox.isChecked() == True:
            return True
        elif self.ui.clear_checkBox.isChecked() == False:
            return False
    
    """Open param window and get peak center range values and assign it to variables to use later"""
    def configure_fit_params(self):
        self.param_window = ParamWindow()
        self.param_window.peak_range.connect(self.peak_range)
    
    def peak_range(self, peaks):
        self.center_min = peaks[0]
        self.center_max = peaks[1]
        
    
    def fit_and_plot(self):
        fit_func = self.ui.fitFunc_comboBox.currentText()
        
        try:
            
            if self.ui.subtract_bck_checkBox.isChecked() == False:
                self.ui.result_textBrowser.setText("You need to check the subtract background option!")
            
            elif self.wlref_file is not None and self.ui.WLRef_checkBox.isChecked() == False:
                self.ui.result_textBrowser.setText("You need to check the White Light Correction option!")
                
            else:
                if fit_func == "Single Gaussian" and self.ui.subtract_bck_checkBox.isChecked() == True:
                    
                    single_gauss = Single_Gaussian(self.file, self.bck_file, wlref=self.wlref_file)
                    
                    if self.ui.adjust_param_checkBox.isChecked():
                        self.result = single_gauss.gaussian_model_w_lims(
                                center_min=self.center_min, center_max=self.center_max)
                    else:
                        self.result = single_gauss.gaussian_model()
                    self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
                    self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                    self.ui.result_textBrowser.setText(self.result.fit_report())
                
                elif fit_func == "Single Lorentzian" and self.ui.subtract_bck_checkBox.isChecked() == True:
                    
                    single_lorentzian = Single_Lorentzian(self.file, self.bck_file, wlref=self.wlref_file)
                    
                    if self.ui.adjust_param_checkBox.isChecked():
                        self.result = single_lorentzian.lorentzian_model_w_lims(
                                center_min = self.center_min, center_max = self.center_max)
                    else:
                        self.result = single_lorentzian.lorentzian_model()
                    self.ui.plot.plot(self.x, self.y, clear=self.clear_check(), pen='r')
                    self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                    self.ui.result_textBrowser.setText(self.result.fit_report())
                
                elif fit_func == "Double Gaussian" and self.ui.subtract_bck_checkBox.isChecked() == True:
                    self.ui.result_textBrowser.setText("Not Implemented Yet!")
                
                elif fit_func == "Multiple Gaussians" and self.ui.subtract_bck_checkBox.isChecked() == True:
                    self.ui.result_textBrowser.setText("Not Implemented Yet!")
        
        except Exception as e:
            self.ui.result_textBrowser.setText(str(e))
    
    def fit_and_plot_scan(self):
#        self.ui.result_textBrowser.append("Starting Scan Fitting")
        print("Starting Scan Fitting")
        
        try:
            """Define starting and stopping wavelength values here"""
            start_nm = int(self.ui.start_nm_spinBox.value())
            stop_nm = int(self.ui.stop_nm_spinBox.value())
            
            ref = self.bck_file
            index = (ref[:,0]>start_nm) & (ref[:,0]<stop_nm)
            
            x = self.spec_scan_file['Wavelengths']
            x = x[index]
            
            data_array = self.spec_scan_file['Intensities']
            
            result_dict = {}
            
            for i in range(data_array.shape[0]):
                
                y = data_array[i, index] # intensity
                yref = ref[index, 1]
                
                y = y - yref # background correction
                y = y - np.mean(y[(x>start_nm) & (x<start_nm + 25)]) # removing any remaining bckgrnd
                
                gmodel = GaussianModel(prefix = 'g1_') # calling gaussian model
                pars = gmodel.guess(y, x=x) # parameters - center, width, height
                result = gmodel.fit(y, pars, x=x, nan_policy='propagate')
                result_dict["result_"+str(i)] = result
            
#            self.ui.result_textBrowser.append("Scan Fitting Complete!")
            print("Scan Fitting Complete!")

            filename = QtWidgets.QFileDialog.getSaveFileName(self)
            pickle.dump(result_dict, open(filename[0]+"_fit_result_dict.pkl", "wb"))
            
#            self.ui.result_textBrowser.append("Data Saved!")
            print("Data Saved!")
        
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass
        
#        self.ui.result_textBrowser.append("Loading Fit Data and Plotting")
        print("Loading Fit Data and Plotting")
        try:
            self.fit_scan_file = pickle.load(open(filename[0]+"_fit_result_dict.pkl", 'rb'))
            self.plot_fit_scan()
            
        except Exception as e:
            self.ui.result_textBrowser.append(str(e))
            pass

    def pub_ready_plot_export(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
        try:
            try:
                data = self.spec_scan_file
                param_selection = str(self.ui.comboBox.currentText())
                if param_selection == 'pk_pos': label = 'PL Peak Position (n.m.)'
                elif param_selection == 'fwhm': label = 'PL FWHM (n.m.)'
                cpm.plot_confocal(self.img, figsize=(10,10), stepsize = data['Scan Parameters']['X step size (um)'], cmap="seismic", cbar_label=label)
                plt.savefig(filename[0],bbox_inches='tight', dpi=300)
                plt.close()
            except:
                plt.figure(figsize=(8,6))
                plt.tick_params(direction='out', length=8, width=3.5)
                plt.plot(self.x, self.y)
                plt.plot(self.x, self.result.best_fit,'k')
                plt.xlabel("Wavelength (nm)", fontsize=20, fontweight='bold')
                plt.ylabel("Intensity (a.u.)", fontsize=20, fontweight='bold')
                plt.tight_layout()
                
                plt.savefig(filename[0],bbox_inches='tight', dpi=300)
                plt.close()
            
        except AttributeError:
            self.ui.result_textBrowser.setText("Need to fit the data first!")

    #Get data from ocean optics scan pkl file, and convert to txt
    def pkl_data_to_txt(self):
        folder = os.path.dirname(self.pkl_to_convert[0])
        filename_ext = os.path.basename(self.pkl_to_convert[0])
        filename = os.path.splitext(filename_ext)[0] #get filename without extension
        pkl_file = pickle.load(open(self.pkl_to_convert[0], 'rb'))

        txt_file = np.zeros(shape=(2048,pkl_file['Intensities'].shape[0] + 1))

        data_array = pkl_file['Intensities']
        data_array = np.transpose(data_array)
        wavelength = pkl_file['Wavelengths']

        txt_file[:,0] = wavelength

        for i in range(pkl_file['Intensities'].shape[0]):
            txt_file[:,i+1] = data_array[:,i]

        np.savetxt(folder +"/"+ filename +"_data.txt", txt_file, fmt = '%.2f', delimiter= "\t", header="wavelength(nm), Intensities at different points")
        self.ui.result_textBrowser.append("Data from .pkl saved as .txt")

    #Get scan parameters from ocean optics scan pkl file, and convert to txt
    def pkl_params_to_txt(self):
        folder = os.path.dirname(self.pkl_to_convert[0])
        filename_ext = os.path.basename(self.pkl_to_convert[0])
        filename = os.path.splitext(filename_ext)[0] #get filename without extension
        pkl_file = pickle.load(open(self.pkl_to_convert[0], 'rb'))

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

        self.ui.result_textBrowser.append("Scan parameters from .pkl saved as .txt")

    def close_application(self):
        choice = QtGui.QMessageBox.question(self, 'EXIT!',
                                            "Do you want to exit the app?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass
        
        
"""Parameter Window GUI and Functions"""
param_file_path = (base_path / "peak_bounds_input.ui").resolve()

param_uiFile = param_file_path

param_WindowTemplate, param_TemplateBaseClass = pg.Qt.loadUiType(param_uiFile)

class ParamWindow(param_TemplateBaseClass):
    
    peak_range = QtCore.pyqtSignal(list)
    
    def __init__(self):
#        super(param_TemplateBaseClass, self).__init__()
        param_TemplateBaseClass.__init__(self)
        
        # Create the param window
        self.pui = param_WindowTemplate()
        self.pui.setupUi(self)
        
        self.pui.pushButton.clicked.connect(self.done)
        
        self.show()
    
    def current_peak_range(self):
        center_min = self.pui.cent_min_doubleSpinBox.value()
        center_max = self.pui.cent_max_doubleSpinBox.value()
        return center_min, center_max
    
    def done(self):
        center_min, center_max = self.current_peak_range()
        self.peak_range.emit([center_min, center_max])
        self.close()
    
"""Run the Main Window"""    
def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win

#Uncomment below if you want to run this as standalone
#run()