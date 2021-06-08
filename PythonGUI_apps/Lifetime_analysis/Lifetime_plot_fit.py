# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

# system imports
import sys
import os
from pathlib import Path

# module imports
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.abspath('../Export_Windows'))
try:
    from Export_window import ExportPlotWindow
except:
    from Export_Windows.Export_window import ExportPlotWindow
    
# local module imports
try:
    from Lifetime_analysis.Fit_functions import stretch_exp_fit, double_exp_fit, single_exp_fit
    from Lifetime_analysis.read_ph_phd import read_picoharp_phd
    from Lifetime_analysis.Fit_functions_with_irf import fit_exp_stretch_diffev, fit_exp_stretch_fmin_tnc, fit_multi_exp_diffev, fit_multi_exp_fmin_tnc
except:
    from Fit_functions import stretch_exp_fit, double_exp_fit, single_exp_fit
    from Fit_functions_with_irf import fit_exp_stretch_diffev, fit_exp_stretch_fmin_tnc, fit_multi_exp_diffev, fit_multi_exp_fmin_tnc
    from read_ph_phd import read_picoharp_phd

"""Recylce params for plotting"""
plt.rc('xtick', labelsize = 20)
plt.rc('xtick.major', pad = 3)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', lw = 2.5, markersize = 7.5)
plt.rc('legend', fontsize = 20)
plt.rc('axes', linewidth=3.5)

pg.mkQApp()
pg.setConfigOption('background', 'w')
##pg.setConfigOption('crashWarning', True)

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
        self.ui.FittingFunc_comboBox.addItems(["Stretched Exponential","Double Exponential", "Single Exponential"])
        self.ui.FittingMethod_comboBox.addItems(["diff_ev", "fmin_tnc"])
        
        #set up file menu
        self.ui.actionOpen.triggered.connect(self.open_file)
        self.ui.actionOpen_IRF_File.triggered.connect(self.open_irf_file)
        self.ui.actionSave.triggered.connect(self.save_file)
        self.ui.actionExit.triggered.connect(self.close_application)
        
        #set up ui signals
        self.ui.plot_pushButton.clicked.connect(self.plot)
        self.ui.fit_pushButton.clicked.connect(self.call_fit_and_plot)
        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        self.ui.export_plot_pushButton.clicked.connect(self.export_window)#pub_ready_plot_export
        self.ui.calculate_srv_pushButton.clicked.connect(self.calculate_srv)

        self.ui.log_checkBox.stateChanged.connect(self.make_semilog)
        self.ui.fit_with_irf_checkBox.stateChanged.connect(self.switch_fit_settings)
        self.ui.FittingFunc_comboBox.currentTextChanged.connect(self.switch_function_tab)
        self.ui.FittingMethod_comboBox.currentTextChanged.connect(self.switch_init_params_groupBox)
        self.ui.separate_irf_checkBox.stateChanged.connect(self.switch_open_irf)
        self.ui.add_to_mem_pushButton.clicked.connect(self.add_trace_to_mem)
        self.ui.export_data_pushButton.clicked.connect(self.export_data)
        self.ui.clear_export_data_pushButton.clicked.connect(self.clear_export_data)
        self.ui.smoothData_checkBox.stateChanged.connect(self.smooth_trace_enabled)

        #set up plot color button
        self.plot_color_button = pg.ColorButton(color=(255,0,0))
        self.ui.plot_color_button_container.layout().addWidget(self.plot_color_button)
        self.plot_color = self.plot_color_button.color()
        self.plot_color_button.sigColorChanged.connect(self.plot_color_changed)

        self.file = None
        self.out = None # output file after fitting
        self.data_list = []
        self.fit_lifetime_called_w_irf = False
        self.fit_lifetime_called_wo_irf = False
        self.x_mem = [] # containers for adding x data to memory
        self.y_mem = [] # containers for adding y data to memory
        self.best_fit_mem = [] # containers for adding best fit data to memory
        self.best_fit_mem_x = [] # containers for adding best fit data to memory
        self.legend = [] # containers for adding legend to memory

        #variables accounting for data received from FLIM analysis
        self.opened_from_flim = False #switched to True in FLIM_plot when "analyze lifetime" clicked
        self.hist_data_from_flim = [] #container for flim roi data

        self.show()
        
    def open_file(self):
        """ Open data file """
#        try:
        self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            if ".csv" in self.filename[0] or ".txt" in self.filename[0]: #if txt or csv, prompt user to enter # of rows to skip
                self.skip_rows_window = SkipRowsWindow()
                self.skip_rows_window.skip_rows_signal.connect(self.open_with_skip_rows_window)
                self.ui.Res_comboBox.setEnabled(True)
            else:
                self.file = read_picoharp_phd(self.filename[0])
            self.opened_from_flim = False
        except:
            pass

    def open_with_skip_rows_window(self):
        """ Prompts user to enter how many rows to skip """
        skip_rows = self.skip_rows_window.ui.skip_rows_spinBox.value()
        if ".txt" in self.filename[0]:
            self.file = np.loadtxt(self.filename[0], skiprows=skip_rows)
            
            if self.file.ndim == 1: # if there is only one trace, reshape to 2D
                self.file = self.file.reshape(self.file.shape[0], 1)
    
        elif ".csv" in self.filename[0]:
            self.file = np.genfromtxt(self.filename[0], skip_header=skip_rows, delimiter=",")

    def open_irf_file(self):
        """ Open file with irf - enabled if 'load separate irf' is checled """
        self.irf_filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            if ".txt" in self.irf_filename[0] or ".csv" in self.irf_filename[0]:
                self.irf_skip_rows_window = SkipRowsWindow()
                self.irf_skip_rows_window.skip_rows_signal.connect(self.open_irf_with_skip_rows_window)
                self.ui.Res_comboBox.setEnabled(True)
            else:
                self.irf_file = read_picoharp_phd(self.irf_filename[0])
        except:
            pass
    
    def open_irf_with_skip_rows_window(self):
        irf_skip_rows = self.irf_skip_rows_window.ui.skip_rows_spinBox.value()
        if ".txt" in self.irf_filename[0]:
            self.irf_file = np.loadtxt(self.irf_filename[0], skiprows=irf_skip_rows)
            
            if self.irf_file.ndim == 1: # if there is only one trace, reshape to 2d array
                self.irf_file = self.irf_file.reshape(self.irf_file.shape[0], 1)
                        
        elif ".csv" in self.irf_filename[0]:
            self.irf_file = np.genfrontxt(self.irf_filename[0], skip_header=irf_skip_rows, delimiter=",")

    def save_file(self):
        try:
            filename = QtWidgets.QFileDialog.getSaveFileName(self)
            np.savetxt(filename[0], self.out, fmt = '%.5f', header = 'Time, Raw_PL, Sim_PL', delimiter = ' ')
        except:
            pass

    def switch_open_irf(self):
        """ Handle 'load separate irf' checkbox """
        self.ui.actionOpen_IRF_File.setEnabled(self.ui.separate_irf_checkBox.isChecked())

    def switch_fit_settings(self):
        """ Enable bounds/initial guess groupboxes only when 'Fit with IRF' is checked """
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
        """ Switch bounds groupbox contents depending on selected fit function """
        fitting_func = self.ui.FittingFunc_comboBox.currentText()
        if fitting_func == "Stretched Exponential":
            self.ui.fitting_params_stackedWidget.setCurrentIndex(0)
        elif fitting_func == "Double Exponential":
            self.ui.fitting_params_stackedWidget.setCurrentIndex(1)
        elif fitting_func == "Single Exponential":
            self.ui.fitting_params_stackedWidget.setCurrentIndex(2)
        
    def switch_init_params_groupBox(self):
        """ Enable initial guess groupbox only when fmin_tnc fit method selected """
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
        """ Grab new plot_color when color button value is changed """
        self.plot_color = self.plot_color_button.color()
    
    def smooth_trace_enabled(self):
        """Enable smooth spin box when smooth data is checked"""
        if self.ui.smoothData_checkBox.isChecked():
            self.ui.smoothData_spinBox.setEnabled(True)
        else:
            self.ui.smoothData_spinBox.setEnabled(False)
    
    def acquire_settings(self, mode="data"):
        """
        Acquire data or irf from channel specified in spinbox.

        mode -- string specifying whether to use data or irf channel (default "data")
        """
        if mode == "data":
            channel = int(self.ui.Data_channel_spinBox.value())
        elif mode == "irf":
            channel = int(self.ui.irf_channel_spinBox.value())
        try:
            try: #
                if self.ui.separate_irf_checkBox.isChecked() and mode=="irf": #if separate irf, get from irf file
                    try:
                        y = self.irf_file[:,channel]
                    except:
                        y = self.irf_file.get_curve(channel)[1]
                else: #otherwise, get data/irf from data file
                    y = self.file[:,channel]

                self.resolution = float(self.ui.Res_comboBox.currentText())
            except:
                res, y = self.file.get_curve(channel)
                time_window = int(np.floor(self.file.get_time_window_in_ns(channel)))
                y = y[0:time_window]
                self.resolution = res
            
            length = np.shape(y)[0]
            x = np.arange(0, length*self.resolution, self.resolution, np.float)
            
            if self.ui.smoothData_checkBox.isChecked() and mode=="data":
                y = np.convolve(y, np.ones(self.ui.smoothData_spinBox.value())/self.ui.smoothData_spinBox.value(), mode="same")
            
            if self.ui.normalize_checkBox.isChecked():
                y = y / np.amax(y)

            return x,y
        
        except Exception as e:
            self.ui.Result_textBrowser.setText(str(e))

    def plot(self):
        try:
            if self.opened_from_flim:
                x, y = self.hist_data_from_flim
            else:
                x,y = self.acquire_settings() #get data
                
            self.ui.plot.plot(x, y, clear=self.ui.clear_plot_checkBox.isChecked(), pen=pg.mkPen(self.plot_color))
            self.fit_lifetime_called_w_irf = False
            self.fit_lifetime_called_wo_irf = False

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
        """ Switch y-log on/off """
        self.ui.plot.setLogMode(False,self.ui.log_checkBox.isChecked())
    
    def clear_plot(self):
        self.ui.plot.clear()
        self.ui.Result_textBrowser.clear()
    
    def fit_and_plot(self):
        """ Fit and plot without IRF """
        try:
            if not hasattr(self, "file"):
                self.ui.Result_textBrowser.setText("You need to load a data file.")
            else:
                if self.opened_from_flim:
                    x, y = self.hist_data_from_flim
                else:
                    x,y = self.acquire_settings() #get data
                y_norm = y/np.max(y) #normalized y

                # find the max intensity in the array and start things from there
                find_max_int = np.nonzero(y_norm == 1)
                y = y[np.asscalar(find_max_int[0]):]
                x = x[np.asscalar(find_max_int[0]):]

                t = x
                time_fit = t
                TRPL_interp = np.interp(time_fit, t, y)
                
                fit_func = self.ui.FittingFunc_comboBox.currentText()
                self.ui.plot.plot(t, y, clear=self.ui.clear_plot_checkBox.isChecked(), pen=pg.mkPen(self.plot_color))
                
                if fit_func == "Stretched Exponential": #stretch exponential tab
                    tc, beta, a, avg_tau, PL_fit, noise = stretch_exp_fit(TRPL_interp, t)
                    self.out = np.empty((len(t), 3))
                    self.out[:,0] = t #time
                    self.out[:,1] = TRPL_interp #Raw PL 
                    self.out[:,2] = PL_fit # PL fit
                    self.ui.plot.plot(t, PL_fit, clear=self.ui.clear_plot_checkBox.isChecked(), pen='k')
                    self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Stretched Exponential"
                                                       "\nFit Method: " + "diff_ev" + #TODO : change when diff_ev and fmin_tnc implemented for non-irf
                                                       "\nAverage Lifetime = " + str(avg_tau)+ " ns"
                                                       "\nCharacteristic Tau = " + str(tc)+" ns"
                                                       "\nBeta = "+str(beta)+
                                                       "\nNoise = "+ str(noise))
                    self.ui.average_lifetime_spinBox.setValue(avg_tau)
                
                elif fit_func == "Double Exponential": #double exponential tab
                    tau1, a1, tau2, a2, avg_tau, PL_fit, noise = double_exp_fit(TRPL_interp, t)
                    self.out = np.empty((len(t), 3))
                    self.out[:,0] = t #time
                    self.out[:,1] = TRPL_interp #Raw PL 
                    self.out[:,2] = PL_fit # PL fit
                    self.ui.plot.plot(t, PL_fit, clear=self.ui.clear_plot_checkBox.isChecked(), pen='k')
                    self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Double Exponential"
                                                       "\nFit Method: " + "diff_ev" +
                                                       "\nAverage Lifetime = " + str(avg_tau)+ " ns"
                                                       "\nTau 1 = " + str(tau1)+" ns"
                                                       "\nA 1 = " + str(a1)+
                                                       "\nTau 2 = " + str(tau2)+" ns"
                                                       "\nA 2 = " + str(a2)+
                                                       "\nNoise = "+ str(noise))
                    #TODO - once tau_avg implemented, set average lifetime spinbox to tau_avg value
                
                elif fit_func == "Single Exponential": #single exponential tab
                    tau, a, PL_fit, noise = single_exp_fit(TRPL_interp, t)
                    self.out = np.empty((len(t), 3))
                    self.out[:,0] = t #time
                    self.out[:,1] = TRPL_interp #Raw PL 
                    self.out[:,2] = PL_fit # PL fit
                    self.ui.plot.plot(t, PL_fit, clear=self.ui.clear_plot_checkBox.isChecked(), pen='k')
                    self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Single Exponential"
                                                       "\nFit Method: " + "diff_ev" +
                                                       "\nLifetime = " + str(tau)+ " ns"
                                                       "\nA = " + str(a)+
                                                       "\nNoise = "+ str(noise))
                    self.ui.average_lifetime_spinBox.setValue(tau)
                
                #add fit params to data_list
                self.data_list.append("Data Channel: " + str(self.ui.Data_channel_spinBox.value()) + "\n" + self.ui.Result_textBrowser.toPlainText())
                self.fit_lifetime_called_wo_irf = True
                self.fit_lifetime_called_w_irf = False

                self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
                self.ui.plot.setLabel('bottom', 'Time (ns)')
                return self.out
        
        except Exception as e:
            self.ui.Result_textBrowser.append(format(e))

    def fit_and_plot_with_irf(self):
        """ Fit and plot with IRF """
        try:
            self.ui.Result_textBrowser.clear()
            if not hasattr(self, "file"):
                self.ui.Result_textBrowser.append("You need to load a data file.")
            if not hasattr(self, "irf_file") and self.ui.separate_irf_checkBox.isChecked():
                self.ui.Result_textBrowser.append("You need to load an IRF file.")
            else:
                if self.opened_from_flim:
                    x,y = self.hist_data_from_flim
                else:
                    x,y = self.acquire_settings() #get data
                _, irf_counts = self.acquire_settings(mode="irf") #get irf counts

                #make sure Irf and data have the same length
                if len(y) != len(irf_counts):
                    y = y[0:min(len(y), len(irf_counts))]
                    irf_counts = irf_counts[0:min(len(y), len(irf_counts))]
                    x = x[0:min(len(y), len(irf_counts))]

                y_norm = y/np.max(y) #normalized y
                irf_norm = irf_counts/np.amax(irf_counts) #normalized irf
                
                t = x
                time_fit = t 
                y = y_norm
                irf_counts = irf_norm
                
                TRPL_interp = np.interp(time_fit, t, y)

                fit_func = self.ui.FittingFunc_comboBox.currentText()
                self.ui.plot.plot(t, y, clear=self.ui.clear_plot_checkBox.isChecked(), pen=pg.mkPen(self.plot_color))
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
                    self.ui.plot.plot(t, bestfit_model, clear=self.ui.clear_plot_checkBox.isChecked(), pen='k')
                    self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Stretched Exponential with IRF"
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
                    self.ui.plot.plot(t, bestfit_model, clear=self.ui.clear_plot_checkBox.isChecked(), pen='k')
                    self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Double Exponential with IRF"
                        "\nFit Method: "+ self.ui.FittingMethod_comboBox.currentText() +
                        "\na1 = %.5f"
                        "\ntau1 = %.5f ns"
                        "\na2 = %.5f"
                        "\ntau2 = %.5f ns"
                        "\nnoise = %.5f counts" %(bestfit_params[0], bestfit_params[1], bestfit_params[2], bestfit_params[3], bestfit_params[4]))
                    #TODO - once tau_avg implemented, set average lifetime spinbox to tau_avg value
                    if bestfit_params[3] > bestfit_params[1]:
                        self.ui.average_lifetime_spinBox.setValue(bestfit_params[3])
                    elif bestfit_params[1] > bestfit_params[3]:
                        self.ui.average_lifetime_spinBox.setValue(bestfit_params[1])

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
                    self.ui.plot.plot(t, bestfit_model, clear=self.ui.clear_plot_checkBox.isChecked(), pen='k')
                    self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Single Exponential with IRF"
                        "\nFit Method: "+ self.ui.FittingMethod_comboBox.currentText() +
                        "\na = %.5f"
                        "\ntau = %.5f ns"
                        "\nnoise = %.5f counts" %(bestfit_params[0], bestfit_params[1], bestfit_params[2]))
                    self.ui.average_lifetime_spinBox.setValue(bestfit_params[1]) #set spinbox to tau value

                #add fit params to data_list
                self.data_list.append("Data Channel: " + str(self.ui.Data_channel_spinBox.value()) + "\n" + self.ui.Result_textBrowser.toPlainText())
                self.fit_lifetime_called_w_irf = True
                self.fit_lifetime_called_wo_irf = False
        except Exception as e:
            self.ui.Result_textBrowser.append(format(e))

    def call_fit_and_plot(self):
        if self.ui.fit_with_irf_checkBox.isChecked():
            self.fit_and_plot_with_irf()
        else:
            self.fit_and_plot()
        if self.ui.calculate_srv_groupBox.isChecked():
            self.calculate_srv() #calculate srv on plot
            self.data_list.append(self.get_srv_string()) #add srv params to data_list
    
    def calculate_surface_lifetime(self):
        effective_lifetime = self.ui.average_lifetime_spinBox.value()
        self.bulk_lifetime = self.ui.bulk_lifetime_spinBox.value() # in ns
        self.surface_lifetime = (effective_lifetime * self.bulk_lifetime)/(self.bulk_lifetime - effective_lifetime)
        self.ui.surface_lifetime_label.setText(str(self.surface_lifetime))
        
    def calculate_srv (self):
        self.calculate_surface_lifetime()
        self.thickness = self.ui.thickness_spinBox.value()*1e-7 # convert to cm
        self.diffusion_coeffecient = self.ui.diffusion_coefficient_spinBox.value() # in cm2/s
        
        self.srv1_srv2_equal = self.thickness / (2*((1e-9*self.surface_lifetime) - ((1/self.diffusion_coeffecient)*((self.thickness/np.pi)**2)) ))
        self.srv1_zero = self.thickness / ((1e-9*self.surface_lifetime) - ((4/self.diffusion_coeffecient)*((self.thickness/np.pi)**2)) )
        
        self.ui.srv1_srv2_equal_label.setText(str(self.srv1_srv2_equal))
        self.ui.srv1_zero_label.setText(str(self.srv1_zero))

    def get_srv_string(self):
        """ Get info from SRV Calculation groupbox as string """
        srv_string = "SRV Calculation:"\
                + "\nAverage Lifetime (ns): " + str(self.ui.average_lifetime_spinBox.value()) \
                + "\nBulk Lifetime (ns): " + str(self.ui.bulk_lifetime_spinBox.value()) \
                + "\nThickness (nm): " + str(self.ui.thickness_spinBox.value()) \
                + "\nDiffusion Coefficient (cm2/s): " + str(self.ui.diffusion_coefficient_spinBox.value())
        srv_string += "\nSurface Lifetime (ns): " + self.ui.surface_lifetime_label.text()
        
        srv_string += "\nSRV1 = SRV2"\
                + "\nSRV (cm/s): " + self.ui.srv1_srv2_equal_label.text()
        
        srv_string += "\nSRV1 = 0"\
                + "\nSRV (cm/s): " + self.ui.srv1_zero_label.text()
        return srv_string
    
    def export_data(self):
        """ Save fit params and srv calculations stored in data_list as .txt """
        folder = os.path.dirname(self.filename[0])
        filename_ext = os.path.basename(self.filename[0])
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
        self.clean_up_after_fig_export()
    
    def clean_up_after_fig_export(self):
        self.x_mem = []
        self.y_mem = []
        self.legend = []
        self.best_fit_mem = []
        self.best_fit_mem_x = []
    
    def add_trace_to_mem(self):
        try:
            if self.fit_lifetime_called_w_irf == True:
                self.x_mem.append(self.out[:,0])
                self.y_mem.append(self.out[:,1])
                self.best_fit_mem_x.append(self.out[:,0])
                self.best_fit_mem.append(self.out[:,2])
            elif self.fit_lifetime_called_wo_irf == True:
                self.x_mem.append(self.acquire_settings()[0])
                self.y_mem.append(self.acquire_settings()[1])
                self.best_fit_mem_x.append(self.out[:,0])
                self.best_fit_mem.append(self.out[:,2])
            else:
                self.x_mem.append(self.acquire_settings()[0])
                self.y_mem.append(self.acquire_settings()[1])
            self.legend.append(self.ui.lineEdit.text())
        except Exception as e:
            print(e)
    
    def export_window(self):
        self.exportplotwindow = ExportPlotWindow()
        self.exportplotwindow.export_fig_signal.connect(self.pub_ready_plot_export)

    def pub_ready_plot_export(self):
        try:
            if self.x_mem == []:
                self.ui.result_textBrowser.setText("Add traces to memory first!")
            
            else:
                filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
                
                plt.figure(figsize=(8,6))
                plt.tick_params(direction='out', length=8, width=3.5)
                for i in range(len(self.x_mem)):
                    plt.plot(self.x_mem[i], self.y_mem[i], label=str(self.legend[i]))
                    if self.fit_lifetime_called_w_irf == True or self.fit_lifetime_called_wo_irf == True:
                        plt.plot(self.best_fit_mem_x[i], self.best_fit_mem[i],'k--')

                plt.yscale('log')
                plt.xlabel("Time (ns)", fontsize=20, fontweight='bold')
                plt.ylabel("Intensity (norm.)", fontsize=20, fontweight='bold')
                plt.legend()
                plt.tight_layout()
                plt.xlim([self.exportplotwindow.ui.lowerX_spinBox.value(),self.exportplotwindow.ui.upperX_spinBox.value()])
                plt.ylim([self.exportplotwindow.ui.lowerY_spinBox.value(),self.exportplotwindow.ui.upperY_doubleSpinBox.value()])
                
                plt.savefig(filename[0],bbox_inches='tight', dpi=300)
                plt.close()
                self.clean_up_after_fig_export()
        
        except Exception as e:
            self.ui.Result_textBrowser.append(format(e))
            pass
            
    def close_application(self):
        choice = QtGui.QMessageBox.question(self, 'EXIT!',
                                            "Do you want to exit the app?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass

"""Skip rows GUI"""
ui_file_path = (base_path / "skip_rows.ui").resolve()
skiprows_WindowTemplate, skiprows_TemplateBaseClass = pg.Qt.loadUiType(ui_file_path)

class SkipRowsWindow(skiprows_TemplateBaseClass):
    
    skip_rows_signal = QtCore.pyqtSignal() #signal to help with pass info back to MainWindow
    
    def __init__(self):
        skiprows_TemplateBaseClass.__init__(self)

        # Create the param window
        self.ui = skiprows_WindowTemplate()
        self.ui.setupUi(self)
        self.ui.done_pushButton.clicked.connect(self.done)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, False)
        self.show()
    
    def done(self):
        self.skip_rows_signal.emit()
        self.close()


def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win

#Uncomment below if you want to run this as standalone
#run()