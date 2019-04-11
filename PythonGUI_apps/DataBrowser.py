# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:41:08 2019

@author: Sarthak
"""

# system imports
from pathlib import Path

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

from Lifetime_analysis import Lifetime_plot_fit
from Spectrum_analysis import Spectra_plot_fit

pg.mkQApp()
pg.setConfigOption('background', 'w')

base_path = Path(__file__).parent
file_path = (base_path / "DataBrowser_GUI.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        TemplateBaseClass.__init__(self)
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        self.ui.select_comboBox.addItems(["Lifetime Analysis", "Spectrum Analysis", "UV-Vis Analysis"])
        self.ui.load_pushButton.clicked.connect(self.load_app)
        
        self.show()

    
    def load_app(self):
        
        analysis_software = self.ui.select_comboBox.currentText()
        
        if analysis_software == "Lifetime Analysis":
            self.lifetime_window = Lifetime_plot_fit.MainWindow()
            self.lifetime_window.show()
        
        elif analysis_software == "Spectrum Analysis":
            self.spectrum_window = Spectra_plot_fit.MainWindow()
            self.spectrum_window.show()
            
        else:
            print("not yet linked -- coming soon")

def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return

run()