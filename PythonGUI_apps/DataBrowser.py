# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 14:41:08 2019

@author: Sarthak
"""

# system imports
import sys
from pathlib import Path

import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore

from Lifetime_analysis import Lifetime_plot_fit
from Spectrum_analysis import Spectra_plot_fit
from FLIM_analysis import FLIM_plot
from UV_Vis_analysis import uv_vis_analysis
from PLQE_analysis import plqe_analysis
from H5_Pkl import h5_pkl_view, h5_view_and_plot

pg.mkQApp()
#pg.setConfigOption('background', 'w')

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
        self.ui.select_comboBox.addItems(["Lifetime Analysis", "Spectrum Analysis", "FLIM Analysis", 
            "UV-Vis Analysis", "PLQE Analysis", "H5 View/Plot", "H5/PKL Viewer"])
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
        elif analysis_software == "FLIM Analysis":
            self.flim_window = FLIM_plot.MainWindow()
            self.flim_window.show()
        elif analysis_software == "UV-Vis Analysis":
            self.uv_vis_window = uv_vis_analysis.MainWindow()
            self.uv_vis_window.show()
        elif analysis_software == "PLQE Analysis":
            self.plqe_window = plqe_analysis.MainWindow()
            self.plqe_window.show()
        elif analysis_software == "H5 View/Plot":
            app  = h5_view_and_plot.H5ViewPlot(sys.argv)
            #sys.exit(app.exec_())
        elif analysis_software == "H5/PKL Viewer":
            app = h5_pkl_view.H5PklView(sys.argv)
            #sys.exit(app.exec_())
        
        

def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return

run()