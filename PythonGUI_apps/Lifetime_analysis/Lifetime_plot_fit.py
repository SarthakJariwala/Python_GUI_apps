# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

# system imports
import sys
from pathlib import Path

# module imports
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import numpy as np
import matplotlib.pyplot as plt

# local module imports
try:
    from Lifetime_analysis.Fit_functions import stretch_exp_fit, double_exp_fit, single_exp_fit
    from Lifetime_analysis.picoharp_phd import read_picoharp_phd
    from Lifetime_analysis.Fit_functions_with_irf import fit_exp_stretch_diffev, fit_exp_stretch_fmin_tnc, fit_multi_exp_diffev, fit_multi_exp_fmin_tnc
except:
    from Fit_functions import stretch_exp_fit, double_exp_fit, single_exp_fit
    from Fit_functions_with_irf import fit_exp_stretch_diffev, fit_exp_stretch_fmin_tnc, fit_multi_exp_diffev, fit_multi_exp_fmin_tnc
    from picoharp_phd import read_picoharp_phd

"""Recylce params for plotting"""
plt.rc('xtick', labelsize = 20)
plt.rc('xtick.major', pad = 3)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', lw = 2.5, markersize = 7.5)
plt.rc('legend', fontsize = 20)
plt.rc('axes', linewidth=3.5)

pg.mkQApp()
pg.setConfigOption('background', 'w')
#pg.setConfigOption('crashWarning', True)

base_path = Path(__file__).parent
file_path = (base_path / "Lifetime_analysis_gui_layout.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        TemplateBaseClass.__init__(self)
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        self.ui.Res_comboBox.addItems(["0.004","0.008","0.016","0.032","0.064","0.128","0.256","0.512"])
        #self.ui.FittingFunc_comboBox.addItems(["Stretched Exponential","Double Exponential", "Single Exponential"])
        
        self.ui.actionOpen.triggered.connect(self.open_file)
        self.ui.actionOpen_IRF_File.triggered.connect(self.open_irf_file)
        self.ui.actionSave.triggered.connect(self.save_file)
        self.ui.actionExit.triggered.connect(self.close_application)
        
        self.ui.plot_pushButton.clicked.connect(self.plot)
        self.ui.fit_pushButton.clicked.connect(self.call_fit_and_plot)
        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        self.ui.export_plot_pushButton.clicked.connect(self.pub_ready_plot_export)
        self.ui.calculate_srv_pushButton.clicked.connect(self.calculate_srv)

        self.ui.log_checkBox.stateChanged.connect(self.make_semilog)
        self.ui.fit_with_irf_checkBox.stateChanged.connect(self.switch_fit_settings)
        self.ui.FittingFunc_comboBox.currentTextChanged.connect(self.switch_function_tab)
        self.ui.FittingMethod_comboBox.currentTextChanged.connect(self.switch_init_params_groupBox)
        self.ui.separate_irf_checkBox.stateChanged.connect(self.switch_open_irf)

        self.plot_color_button = pg.ColorButton(color=(255,0,0))
        self.ui.plot_color_button_container.layout().addWidget(self.plot_color_button)
        self.plot_color = self.plot_color_button.color()
        self.plot_color_button.sigColorChanged.connect(self.plot_color_changed)
        self.file = None
        self.out = None # output file after fitting
        
        self.show()
        
    def open_file(self):
#        try:
        filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            self.file = np.loadtxt(filename[0], skiprows=10)
#        except ValueError:
#            self.file = np.loadtxt(filename[0], skiprows=10)
        except UnicodeDecodeError:
            self.file = read_picoharp_phd(filename[0])
        except:
            pass

    def open_irf_file(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            self.irf_file = np.loadtxt(filename[0], skiprows=10)
        except UnicodeDecodeError:
            self.file = read_picoharp_phd(filename[0])
        except:
            pass
    
    def save_file(self):
        try:
            filename = QtWidgets.QFileDialog.getSaveFileName(self)
            np.savetxt(filename[0], self.out, fmt = '%.5f', header = 'Time, Raw_PL, Sim_PL', delimiter = ' ')
        except:
            pass

    def switch_open_irf(self):
        self.ui.actionOpen_IRF_File.setEnabled(self.ui.separate_irf_checkBox.isChecked())


    def switch_fit_settings(self):
        checked = self.ui.fit_with_irf_checkBox.isChecked()
        for func in "str de se".split(" "):
            boundsGb = eval("self.ui."+func+"_bounds_groupBox")
            #initGb = eval("self.ui."+func+"_init_groupBox")
            boundsGb.setEnabled(checked)
            #initGb.setEnabled(checked)
            if checked == True:
                self.switch_init_params_groupBox()
            else:
                initGb = eval("self.ui."+func+"_init_groupBox")
                initGb.setEnabled(checked)
        self.ui.FittingMethod_comboBox.setEnabled(checked)

    def switch_function_tab(self):
        fitting_func = self.ui.FittingFunc_comboBox.currentText()
        if fitting_func == "Stretched Exponential":
            self.ui.fitting_params_stackedWidget.setCurrentIndex(0)
        elif fitting_func == "Double Exponential":
            self.ui.fitting_params_stackedWidget.setCurrentIndex(1)
        elif fitting_func == "Single Exponential":
            self.ui.fitting_params_stackedWidget.setCurrentIndex(2)
        

    def switch_init_params_groupBox(self):
        if self.ui.FittingMethod_comboBox.currentText() == "diff_ev":
            for func in "str de se".split(" "):
                initGb = eval("self.ui."+func+"_init_groupBox")
                initGb.setEnabled(False)
            #initGb.setEnabled(checked)
        elif self.ui.FittingMethod_comboBox.currentText() == "fmin_tnc":
            for func in "str de se".split(" "):
                initGb = eval("self.ui."+func+"_init_groupBox")
                initGb.setEnabled(True)

    def plot_color_changed(self):
        self.plot_color = self.plot_color_button.color()
    
    def acquire_settings(self, mode="data"):# mode --looks whether argument is data or irf
        self.resolution = float(self.ui.Res_comboBox.currentText())
        if mode == "data":
            channel = int(self.ui.Data_channel_spinBox.value())
        elif mode == "irf":
            channel = int(self.ui.irf_channel_spinBox.value())
        try:
            try:
                if self.ui.separate_irf_checkBox.isChecked() and mode=="irf":
                    y = self.irf_file[:,channel]
                else:
                    y = self.file[:,channel]
            except:
                res, y = self.file.get_curve(channel)
                # TO DO - check if res read in is the same as selected
                time_window = int(np.floor(self.file.get_time_window_in_ns(channel)))
                y = y[0:time_window]
            
            length = np.shape(y)[0]
            x = np.arange(0, length, 1) * self.resolution
            return x,y
        
        except Exception as e:
            self.ui.Result_textBrowser.setText(str(e))

    def plot(self):
        try:
            x,y = self.acquire_settings()
            self.ui.plot.plot(x, y, clear=False, pen=pg.mkPen(self.plot_color))

            try:
                self.ui.Result_textBrowser.setText("Integral Counts :\n" "{:.2E}".format(
                        self.file.get_integral_counts(int(self.ui.Channel_spinBox.value()))))
            except:
                self.ui.Result_textBrowser.setText("Integral Counts :\n" "{:.2E}".format(np.sum(y)))
        except:
            pass
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Time (ns)')
        
    
    def make_semilog(self):
        self.ui.plot.setLogMode(False,self.ui.log_checkBox.isChecked())
    
    def clear_plot(self):
        self.ui.plot.clear()
        self.ui.Result_textBrowser.clear()
    
    def fit_and_plot(self):
        try:
            x,y = self.acquire_settings()
            
            y_norm = y/np.max(y)
            # find the max intensity in the array and start things from there
            find_max_int = np.nonzero(y_norm == 1)
            y = y[np.asscalar(find_max_int[0]):]
            x = x[np.asscalar(find_max_int[0]):]
            
            t = x
            
            time_fit = t 
            TRPL_interp = np.interp(time_fit, t, y)
            
            fit_func = self.ui.FittingFunc_comboBox.currentText()
            self.ui.plot.plot(t, y, clear=True, pen=pg.mkPen(self.plot_color))
            
            if fit_func == "Stretched Exponential": #stretch exponential tab
                tc, beta, a, avg_tau, PL_fit = stretch_exp_fit(TRPL_interp, t)
                self.out = np.empty((len(t), 3))
                self.out[:,0] = t #time
                self.out[:,1] = TRPL_interp #Raw PL 
                self.out[:,2] = PL_fit # PL fit
                self.ui.plot.plot(t, PL_fit, clear=False, pen='k')
                self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Stretched Exponential"
                                                   "\nFit Method: " + "diff_ev" + #TODO : change when diff_ev and fmin_tnc implemented for non-irf
                                                   "\nAverage Lifetime = " + str(avg_tau)+ " ns"
                                                   "\nCharacteristic Tau = " + str(tc)+" ns"
                                                   "\nBeta = "+str(beta))
            
            elif fit_func == "Double Exponential": #double exponential tab
                tau1, a1, tau2, a2, avg_tau, PL_fit = double_exp_fit(TRPL_interp, t)
                self.out = np.empty((len(t), 3))
                self.out[:,0] = t #time
                self.out[:,1] = TRPL_interp #Raw PL 
                self.out[:,2] = PL_fit # PL fit
                self.ui.plot.plot(t, PL_fit, clear=False, pen='k')
                self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Double Exponential"
                                                   "\nFit Method: " + "diff_ev" +
                                                   "\nAverage Lifetime = " + str(avg_tau)+ " ns"
                                                   "\nTau 1 = " + str(tau1)+" ns"
                                                   "\nA 1 = " + str(a1)+
                                                   "\nTau 2 = " + str(tau2)+" ns"
                                                   "\nA 2 = " + str(a2))
            
            elif fit_func == "Single Exponential": #single exponential tab
                tau, a, PL_fit = single_exp_fit(TRPL_interp, t)
                self.out = np.empty((len(t), 3))
                self.out[:,0] = t #time
                self.out[:,1] = TRPL_interp #Raw PL 
                self.out[:,2] = PL_fit # PL fit
                self.ui.plot.plot(t, PL_fit, clear=False, pen='k')
                self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Single Exponential"
                                                   "\nFit Method: " + "diff_ev" +
                                                   "\nLifetime = " + str(tau)+ " ns"
                                                   "\nA = " + str(a))
                
            self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
            self.ui.plot.setLabel('bottom', 'Time (ns)')
            return self.out
        
        except Exception as err:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print(exc_type, exc_tb.tb_lineno)

    def fit_and_plot_with_irf(self):
        try:
            x,y = self.acquire_settings()
            _, irf_counts = self.acquire_settings(mode="irf")

            #make sure Irf and data have the same length
            if len(y) != len(irf_counts):
                y = y[0:min(len(y), len(irf_counts))]
                irf_counts = irf_counts[0:min(len(y), len(irf_counts))]
                x = x[0:min(len(y), len(irf_counts))]

            y_norm = y/np.max(y)
            irf_norm = irf_counts/np.max(irf_counts)
            
            t = x
            time_fit = t 
            y = y_norm
            irf_counts = irf_norm
            
            TRPL_interp = np.interp(time_fit, t, y)

            fit_func = self.ui.FittingFunc_comboBox.currentText()
            self.ui.plot.plot(t, y, clear=True, pen=pg.mkPen(self.plot_color))
            if fit_func == "Stretched Exponential": #stretched exponential tab
                tc_bounds = (self.ui.str_tc_min_spinBox.value(), self.ui.str_tc_max_spinBox.value()) #(0, 10000)
                a_bounds = (self.ui.str_a_min_spinBox.value(), self.ui.str_a_max_spinBox.value())#(0.9, 1.1)
                beta_bounds = (self.ui.str_beta_min_spinBox.value(), self.ui.str_beta_max_spinBox.value())#(0,1)
                noise_bounds = (self.ui.str_noise_min_spinBox.value(), self.ui.str_noise_max_spinBox.value())#(0, 1e4)
                stretch_exp_bounds = [tc_bounds, beta_bounds, a_bounds, noise_bounds]
                stretch_exp_init_params = [self.ui.str_tc_init_spinBox.value(), self.ui.str_a_init_spinBox.value(), self.ui.str_beta_init_spinBox.value(), self.ui.str_noise_init_spinBox.value()]

                #tc, beta, a, avg_tau, PL_fit = stretch_exp_fit(TRPL_interp, t)
#                resolution = float(self.ui.Res_comboBox.currentText())
                if self.ui.FittingMethod_comboBox.currentText() == "diff_ev":
                    bestfit_params, t_avg, bestfit_model, data_array, time_array, irf = fit_exp_stretch_diffev(t, self.resolution, TRPL_interp, irf_counts, stretch_exp_bounds)
                else: #if fmin_tnc fitting method selected
                    bestfit_params, t_avg, bestfit_model, data_array, time_array, irf  = fit_exp_stretch_fmin_tnc(t, self.resolution, TRPL_interp, irf_counts, stretch_exp_init_params, stretch_exp_bounds)
                self.out = np.empty((len(t), 3))
                self.out[:,0] = t #time
                self.out[:,1] = TRPL_interp #Raw PL 
                self.out[:,2] = bestfit_model # PL fit
                self.ui.plot.plot(t, bestfit_model, clear=False, pen='k')
                self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Stretched Exponential"
                    "\nFit Method: "+ self.ui.FittingMethod_comboBox.currentText() +
                    "\ntau_avg = %.5f ns"
                    "\nbeta = %.5f"
                    "\ntau_c = %.5f ns"
                    "\na = %.5f \nnoise = %.5f counts" %(t_avg, bestfit_params[1], bestfit_params[0], bestfit_params[2], bestfit_params[3]))
                #self.effective_lifetime = t_avg
                self.ui.average_lifetime_spinBox.setValue(t_avg)
            
            elif fit_func == "Double Exponential": #double exponential tab
                a1_bounds = (self.ui.de_a1_min_spinBox.value(), self.ui.de_a1_max_spinBox.value())
                tau1_bounds = (self.ui.de_tau1_min_spinBox.value(), self.ui.de_tau1_max_spinBox.value())
                a2_bounds = (self.ui.de_a2_min_spinBox.value(), self.ui.de_a2_max_spinBox.value())
                tau2_bounds = (self.ui.de_tau2_min_spinBox.value(), self.ui.de_tau2_max_spinBox.value())
                noise_bounds = (self.ui.de_noise_min_spinBox.value(), self.ui.de_noise_max_spinBox.value())
                double_exp_bounds = [a1_bounds, tau1_bounds, a2_bounds, tau2_bounds, noise_bounds]
                double_exp_init_params = [self.ui.de_a1_init_spinBox.value(), self.ui.de_tau1_init_spinBox.value(), self.ui.de_a2_init_spinBox.value(), 
                    self.ui.de_tau2_init_spinBox.value(), self.ui.de_noise_init_spinBox.value()]

                if self.ui.FittingMethod_comboBox.currentText() == "diff_ev":
                    bestfit_params, bestfit_model, data_array, time_array, irf = fit_multi_exp_diffev(t, self.resolution, TRPL_interp, irf_counts,  double_exp_bounds, 2)
                    #bestfit_params, bestfit_model, data_array, time_array, irf = fit_multi_exp_diffev(t, resolution, TRPL_interp, irf_counts, double_exp_init_bounds, 2)
                else:
                    bestfit_params, bestfit_model, data_array, time_array, irf = fit_multi_exp_fmin_tnc(t, self.resolution, TRPL_interp, irf_counts, double_exp_init_params, double_exp_bounds, 2)
                self.out = np.empty((len(t), 3))
                self.out[:,0] = t #time
                self.out[:,1] = TRPL_interp #Raw PL 
                self.out[:,2] = bestfit_model # PL fit
                self.ui.plot.plot(t, bestfit_model, clear=False, pen='k')
                self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Double Exponential"
                    "\nFit Method: "+ self.ui.FittingMethod_comboBox.currentText() +
                    "\na1 = %.5f"
                    "\ntau1 = %.5f ns"
                    "\na2 = %.5f"
                    "\ntau2 = %.5f ns"
                    "\nnoise = %.5f counts" %(bestfit_params[0], bestfit_params[1], bestfit_params[2], bestfit_params[3], bestfit_params[4]))

            elif fit_func == "Single Exponential": #single exponential tab
                a_bounds = (self.ui.se_a_min_spinBox.value(), self.ui.se_a_max_spinBox.value())
                tau_bounds = (self.ui.se_tau_min_spinBox.value(), self.ui.se_tau_max_spinBox.value())
                noise_bounds = (self.ui.se_noise_min_spinBox.value(), self.ui.se_noise_max_spinBox.value())
                single_exp_bounds = [a_bounds, tau_bounds, noise_bounds]
                single_exp_init_params = [self.ui.se_a_init_spinBox.value(), self.ui.se_tau_init_spinBox.value(), self.ui.se_noise_init_spinBox.value()]

                if self.ui.FittingMethod_comboBox.currentText() == "diff_ev":
                    bestfit_params, bestfit_model, data_array, time_array, irf = fit_multi_exp_diffev(t, self.resolution, TRPL_interp, irf_counts, single_exp_bounds, 1)
                else:
                    bestfit_params, bestfit_model, data_array, time_array, irf = fit_multi_exp_fmin_tnc(t, self.resolution, TRPL_interp, irf_counts, single_exp_init_params, single_exp_bounds, 1)
                self.out = np.empty((len(t), 3))
                self.out[:,0] = t #time
                self.out[:,1] = TRPL_interp #Raw PL 
                self.out[:,2] = bestfit_model # PL fit
                self.ui.plot.plot(t, bestfit_model, clear=False, pen='k')
                self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Single Exponential"
                    "\nFit Method: "+ self.ui.FittingMethod_comboBox.currentText() +
                    "\na = %.5f"
                    "\ntau = %.5f ns"
                    "\nnoise = %.5f counts" %(bestfit_params[0], bestfit_params[1], bestfit_params[2]))

        
        except Exception as err:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            print(exc_type, exc_tb.tb_lineno)

    def call_fit_and_plot(self):
        if self.ui.fit_with_irf_checkBox.isChecked():
            self.fit_and_plot_with_irf()
            if self.ui.calculate_srv_groupBox.isChecked() and self.ui.FittingFunc_comboBox.currentText() == "Stretched Exponential":
                #self.calculate_surface_lifetime()
                self.calculate_srv()
        else:
            self.fit_and_plot()
    
    def calculate_surface_lifetime(self):
        effective_lifetime = self.ui.average_lifetime_spinBox.value()
        self.bulk_lifetime = self.ui.bulk_lifetime_spinBox.value() # in ns
        self.surface_lifetime = (effective_lifetime * self.bulk_lifetime)/(self.bulk_lifetime - effective_lifetime)
        self.ui.surface_lifetime_label.setText(str(self.surface_lifetime))
        
    def calculate_srv (self):
        self.calculate_surface_lifetime()
        self.thickness = self.ui.thickness_spinBox.value()*1e-7 # convert to cm
        self.diffusion_coeffecient = self.ui.diffusion_coefficient_spinBox.value() # in cm2/s
        
        if self.ui.srv1_srv2_checkBox.isChecked():
            self.srv = self.thickness / (2*((1e-9*self.surface_lifetime) - ((1/self.diffusion_coeffecient)*((self.thickness/np.pi)**2)) ))
        else:
            self.srv = self.thickness / ((1e-9*self.surface_lifetime) - ((4/self.diffusion_coeffecient)*((self.thickness/np.pi)**2)) )
        
        self.ui.srv_label.setText(str(self.srv))
    
    def pub_ready_plot_export(self):
        try:
            filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
            
            plt.figure(figsize=(8,6))
            plt.tick_params(direction='out', length=8, width=3.5)
            if self.ui.save_w_fit_checkBox.isChecked():
                plt.plot(self.out[:,0],self.out[:,1]/np.max(self.out[:,1]))
                plt.plot(self.out[:,0],self.out[:,2]/np.max(self.out[:,1]),'k')
            else:
                plt.plot(self.acquire_settings()[0],self.acquire_settings()[1]/np.max(self.acquire_settings()[1]))
            plt.yscale('log')
            plt.xlabel("Time (ns)", fontsize=20, fontweight='bold')
            plt.ylabel("Intensity (norm.)", fontsize=20, fontweight='bold')
            plt.tight_layout()
            
            plt.savefig(filename[0],bbox_inches='tight', dpi=300)
            plt.close()
        
        except:
            pass
            
    def close_application(self):
        choice = QtGui.QMessageBox.question(self, 'EXIT!',
                                            "Do you want to exit the app?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass

def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win

#Uncomment below if you want to run this as standalone
#run()