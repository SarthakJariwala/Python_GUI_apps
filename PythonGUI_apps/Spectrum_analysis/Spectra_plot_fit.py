# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

# system imports
import sys
from pathlib import Path

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import matplotlib.pyplot as plt

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

base_path = Path(__file__).parent
file_path = (base_path / "Spectra_plot_fit_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        TemplateBaseClass.__init__(self)
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        self.ui.fitFunc_comboBox.addItems(["Single Gaussian","Single Lorentzian", "Double Gaussian", "Multiple Gaussians"])
        
#        self.ui.actionSave.triggered.connect(self.save_file)
#        self.ui.actionExit.triggered.connect(self.close_application)
        
        self.ui.importSpec_pushButton.clicked.connect(self.open_file)
        self.ui.importBck_pushButton.clicked.connect(self.open_bck_file)
        self.ui.importWLRef_pushButton.clicked.connect(self.open_wlref_file)
        self.ui.plot_pushButton.clicked.connect(self.plot)
        self.ui.fit_pushButton.clicked.connect(self.fit_and_plot)
        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        self.ui.export_fig_pushButton.clicked.connect(self.pub_ready_plot_export)
        
        self.file = None
        self.bck_file = None
        self.wlref_file = None
        self.x = None
        self.y = None
        self.out = None # output file after fitting
        
        self.show()
        
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
        except:
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
    
    def save_file(self):
        try:
            filename = QtWidgets.QFileDialog.getSaveFileName(self)
            np.savetxt(filename[0], self.out, fmt = '%.5f', header = 'Time, Raw_PL, Sim_PL', delimiter = ' ')
        except:
            pass

    def plot(self):
        
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
            
        self.ui.plot.plot(self.x, self.y, clear=False, pen='r')
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Wavelength (nm)')
        
    
    def normalize(self):
        self.y = (self.y) / np.amax(self.y)
    
    def clear_plot(self):
        self.ui.plot.clear()
#        self.ui.Result_textBrowser.clear()
    
    def fit_and_plot(self):
        fit_func = self.ui.fitFunc_comboBox.currentText()
        
        if self.ui.subtract_bck_checkBox.isChecked() == False:
            self.ui.result_textBrowser.setText("You need to check the subtract background option!")
        
        elif self.wlref_file is not None and self.ui.WLRef_checkBox.isChecked() == False:
            self.ui.result_textBrowser.setText("You need to check the White Light Correction option!")
            
        else:
            if fit_func == "Single Gaussian" and self.ui.subtract_bck_checkBox.isChecked() == True:
                
                single_gauss = Single_Gaussian(self.file, self.bck_file, wlref=self.wlref_file)
                self.result = single_gauss.gaussian_model()
                self.ui.plot.plot(self.x, self.y, clear=True, pen='r')
                self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                self.ui.result_textBrowser.setText(self.result.fit_report())
            
            elif fit_func == "Single Lorentzian" and self.ui.subtract_bck_checkBox.isChecked() == True:
                
                single_lorentzian = Single_Lorentzian(self.file, self.bck_file, wlref=self.wlref_file)
                self.result = single_lorentzian.lorentzian_model()
                self.ui.plot.plot(self.x, self.y, clear=True, pen='r')
                self.ui.plot.plot(self.x, self.result.best_fit, clear=False, pen='k')
                self.ui.result_textBrowser.setText(self.result.fit_report())
            
            elif fit_func == "Double Gaussian" and self.ui.subtract_bck_checkBox.isChecked() == True:
                self.ui.result_textBrowser.setText("Not Implemented Yet!")
            
            elif fit_func == "Multiple Gaussians" and self.ui.subtract_bck_checkBox.isChecked() == True:
                self.ui.result_textBrowser.setText("Not Implemented Yet!")
            
    
    def pub_ready_plot_export(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
        try:
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