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
except:
    from Fit_functions import stretch_exp_fit, double_exp_fit, single_exp_fit
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
        self.ui.Channel_comboBox.addItems(["0","1","2","3","4","5","6","7"])
        self.ui.FittingFunc_comboBox.addItems(["Stretched Exponential","Double Exponential", "Single Exponential"])
        
        self.ui.actionOpen.triggered.connect(self.open_file)
        self.ui.actionSave.triggered.connect(self.save_file)
        self.ui.actionExit.triggered.connect(self.close_application)
        
        self.ui.plot_pushButton.clicked.connect(self.plot)
        self.ui.log_pushButton.clicked.connect(self.make_semilog)
        self.ui.fit_pushButton.clicked.connect(self.fit_and_plot)
        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        self.ui.export_plot_pushButton.clicked.connect(self.pub_ready_plot_export)
        
        self.file = None
        self.out = None # output file after fitting
        
        self.show()
        
    def open_file(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            self.file = np.loadtxt(filename[0], skiprows=0)
        except ValueError:
            self.file = np.loadtxt(filename[0], skiprows=10)
        except:
            self.file = read_picoharp_phd(filename[0])
    
    def save_file(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self)
        np.savetxt(filename[0], self.out, fmt = '%.5f', header = 'Time, Raw_PL, Sim_PL', delimiter = ' ')
    
    def acquire_settings(self):
        resolution = float(self.ui.Res_comboBox.currentText())
        channel = int(self.ui.Channel_comboBox.currentText())
        try:
            y = self.file[:,channel]
        except:
            res, y = self.file.get_curve(channel)
            # TO DO - check if res read in is the same as selected
            time_window = int(np.floor(self.file.get_time_window_in_ns(channel)))
            y = y[0:time_window]
        
        length = np.shape(y)[0]
        x = np.arange(0, length, 1) * resolution
        return x,y

    def plot(self):
        x,y = self.acquire_settings()
        self.ui.plot.plot(x, y, clear=False, pen='r')
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Time', units='ns')
        try:
            self.ui.Result_textBrowser.setText("Integral Counts :\n" "{:.2E}".format(
                    self.file.get_integral_counts(int(self.ui.Channel_comboBox.currentText()))))
        except:
            self.ui.Result_textBrowser.setText("Integral Counts :\n" "{:.2E}".format(np.sum(y)))
    
    def make_semilog(self):
        self.ui.plot.setLogMode(False,True)
    
    def clear_plot(self):
        self.ui.plot.clear()
        self.ui.Result_textBrowser.clear()
    
    def fit_and_plot(self):
        
        x,y = self.acquire_settings()
        
        y_norm = y/np.max(y)
        # find the max intensity in the array and start things from there
        find_max_int = np.nonzero(y_norm == 1)
        y = y[np.asscalar(find_max_int[0]):]
        
        resolution = float(self.ui.Res_comboBox.currentText())
#        x = x[np.asscalar(find_max_int[0]):]
        x = np.arange(0, len(y), 1) * resolution
        
        t = x
        
        time_fit = t 
        TRPL_interp = np.interp(time_fit, t, y)
        
        fit_func = self.ui.FittingFunc_comboBox.currentText()
        self.ui.plot.plot(t, y, clear=True, pen='r')
        
        if fit_func == "Stretched Exponential":
            tc, beta, a, avg_tau, PL_fit = stretch_exp_fit(TRPL_interp, t)
            self.out = np.empty((len(t), 3))
            self.out[:,0] = t #time
            self.out[:,1] = TRPL_interp #Raw PL 
            self.out[:,2] = PL_fit # PL fit
            self.ui.plot.plot(t, PL_fit, clear=False, pen='k')
            self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Stretched Exponential"
                                               "\nAverage Lifetime = " + str(avg_tau)+ " ns"
                                               "\nCharacteristic Tau = " + str(tc)+" ns"
                                               "\nBeta = "+str(beta))
        
        elif fit_func == "Double Exponential":
            tau1, a1, tau2, a2, avg_tau, PL_fit = double_exp_fit(TRPL_interp, t)
            self.out = np.empty((len(t), 3))
            self.out[:,0] = t #time
            self.out[:,1] = TRPL_interp #Raw PL 
            self.out[:,2] = PL_fit # PL fit
            self.ui.plot.plot(t, PL_fit, clear=False, pen='k')
            self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Double Exponential"
                                               "\nAverage Lifetime = " + str(avg_tau)+ " ns"
                                               "\nTau 1 = " + str(tau1)+" ns"
                                               "\nA 1 = " + str(a1)+
                                               "\nTau 2 = " + str(tau2)+" ns"
                                               "\nA 2 = " + str(a2))
        
        elif fit_func == "Single Exponential":
            tau, a, PL_fit = single_exp_fit(TRPL_interp, t)
            self.out = np.empty((len(t), 3))
            self.out[:,0] = t #time
            self.out[:,1] = TRPL_interp #Raw PL 
            self.out[:,2] = PL_fit # PL fit
            self.ui.plot.plot(t, PL_fit, clear=False, pen='k')
            self.ui.Result_textBrowser.setText("Fit Results:\n\nFit Function: Single Exponential"
                                               "\nLifetime = " + str(tau)+ " ns"
                                               "\nA = " + str(a))
            
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Time', units='ns')
        return self.out
    
    def pub_ready_plot_export(self):
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