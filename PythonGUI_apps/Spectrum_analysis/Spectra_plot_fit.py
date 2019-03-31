# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import sys
#from Fit_functions import stretch_exp_fit, double_exp_fit, single_exp_fit

pg.mkQApp()
pg.setConfigOption('background', 'w')

uiFile = "Spectra_plot_fit_gui.ui"

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        TemplateBaseClass.__init__(self)
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
#        self.ui.FittingFunc_comboBox.addItems(["Strected Exponential","Double Exponential", "Single Exponential"])
        
#        self.ui.actionSave.triggered.connect(self.save_file)
#        self.ui.actionExit.triggered.connect(self.close_application)
        
        self.ui.importSpec_pushButton.clicked.connect(self.open_file)
        self.ui.importBck_pushButton.clicked.connect(self.open_bck_file)
        self.ui.importWLRef_pushButton.clicked.connect(self.open_wlref_file)
        self.ui.plot_pushButton.clicked.connect(self.plot)
        
#        self.ui.fit_pushButton.clicked.connect(self.fit_and_plot)
        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        
        self.file = None
        self.bck_file = None
        self.wlref_file = None
        self.x = None
        self.y = None
        self.out = None # output file after fitting
        
        self.show()
        
    def open_file(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            self.file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
        except:
            self.file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')
    
    def open_bck_file(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            self.bck_file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
        except:
            self.bck_file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')
    
    def open_wlref_file(self):
        filename = QtWidgets.QFileDialog.getOpenFileName(self)
        try:
            self.wlref_file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
        except:
            self.wlref_file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')
    
    def save_file(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self)
        np.savetxt(filename[0], self.out, fmt = '%.5f', header = 'Time, Raw_PL, Sim_PL', delimiter = ' ')

    def plot(self):
#        x,y = self.acquire_settings()
        self.x = self.file[:,0]
        self.y = self.file[:,1]
        
        if self.ui.subtract_bck_checkBox.isChecked() == True:
            bck_y = self.bck_file[:,1]
            self.y = self.y - bck_y
        
        elif self.ui.subtract_bck_checkBox.isChecked() == True and self.ui.WLRef_checkBox.isChecked() ==True:
            bck_y = self.bck_file[:,1]
            wlref_y = self.wlref_file[:,1]
            self.y = (self.y-bck_y)/wlref_y
        
        elif self.ui.WLRef_checkBox.isChecked() == True:
            wlref_y = self.wlref_file[:,1]
            self.y = (self.y)/wlref_y
        
        if self.ui.norm_checkBox.isChecked():
            self.normalize()
            
        self.ui.plot.plot(self.x, self.y, clear=False, pen='r')
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Wavelength', units='nm')
        
    
    def normalize(self):
        self.y = (self.y) / np.amax(self.y)
        
    
    def clear_plot(self):
        self.ui.plot.clear()
#        self.ui.Result_textBrowser.clear()
    
    def close_application(self):
        choice = QtGui.QMessageBox.question(self, 'EXIT!',
                                            "Do you want to exit the app?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass
        
        
win = MainWindow()


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
