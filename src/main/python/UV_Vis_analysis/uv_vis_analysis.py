# system imports
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph import exporters
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import matplotlib.pyplot as plt

import numpy as np
import time

# local modules

pg.mkQApp()
pg.setConfigOption('background', 'w')

base_path = Path(__file__).parent
file_path = (base_path / "uv_vis_analysis_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

"""params for plotting"""
plt.rc('xtick', labelsize = 20)
plt.rc('xtick.major', pad = 3)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', lw = 1.5, markersize = 7.5)
plt.rc('legend', fontsize = 20)
plt.rc('axes', linewidth = 3.5)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        super(TemplateBaseClass, self).__init__()
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        #setup uv vis plot
        self.absorbance_plot_layout = pg.GraphicsLayoutWidget()
        self.ui.absorbance_plot_container.layout().addWidget(self.absorbance_plot_layout)
        self.absorbance_plot = self.absorbance_plot_layout.addPlot(title="Wavelengths vs. Absorbance")
        self.absorbance_plot.setLabel('bottom', 'Wavelength', unit='nm')
        self.absorbance_plot.setLabel('left', 'Absorbance', unit='a.u.')

        #setup correction region for uv vis
        self.correction_region = pg.LinearRegionItem()
        self.correction_region_min = 600
        self.correction_region_max = 900
        self.correction_region.setRegion((self.correction_region_min, self.correction_region_max))
        
        #setup uv vis ui signals
        self.ui.correct_for_scattering_checkBox.stateChanged.connect(self.scattering_checkBox_state)
        self.scatter_corrected = False
        self.ui.actionLoad_data.triggered.connect(self.open_data_file)
        self.ui.plot_absorbance_pushButton.clicked.connect(self.plot_absorbance)
        self.ui.clear_uvvis_pushButton.clicked.connect(self.clear_uvvis)
        self.ui.export_uv_vis_pushButton.clicked.connect(self.export_uv_vis)
        self.correction_region.sigRegionChanged.connect(self.update_correction_region)

        #setup tauc plot
        self.tauc_plot_layout = pg.GraphicsLayoutWidget()
        self.ui.tauc_plot_container.layout().addWidget(self.tauc_plot_layout)
        self.tauc_plot = self.tauc_plot_layout.addPlot(title="Tauc plot fit")
        self.tauc_plot.setLabel('bottom', 'hv', unit='ev')
        y_label = '(ahv)' + chr(0x00B2) #char is superscripted 2
        self.tauc_plot.setLabel('left', y_label)

        #setup tauc ui signals
        self.ui.plot_tauc_pushButton.clicked.connect(self.plot_tauc)
        self.ui.clear_tauc_pushButton.clicked.connect(self.clear_tauc)
        self.ui.export_tauc_pushButton.clicked.connect(self.export_tauc)

        self.show()

    def open_data_file(self):
        try:
            self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
            self.data = np.loadtxt(self.filename[0], delimiter = ',', skiprows = 2)
            self.Wavelength = self.data[:,0] # in nm
            self.Absorbance = self.data[:,1]
        except Exception as err:
            print(format(err))

    def update_correction_region(self):
        """ Update correction region variables from region """
        self.correction_region_min, self.correction_region_max = self.correction_region.getRegion()
    
    def scattering_checkBox_state(self):
        if self.ui.correct_for_scattering_checkBox.isChecked():
            self.scatter_corrected = True
            self.ui.mean_radioButton.setEnabled(True)
            self.ui.fourth_orderpoly_radioButton.setEnabled(True)
        else:
            self.scatter_corrected = False
            self.ui.mean_radioButton.setEnabled(False)
            self.ui.fourth_orderpoly_radioButton.setEnabled(False)

    def plot_absorbance(self):
        try:
            self.plotted_absorbance = self.Absorbance #by default set to original absorbance data
            if self.scatter_corrected == True:
                if self.ui.fourth_orderpoly_radioButton.isChecked():
                    p = (np.polyfit(self.Wavelength[(self.Wavelength>self.correction_region_min) & (self.Wavelength<self.correction_region_max)],
                                                    self.Absorbance[(self.Wavelength>self.correction_region_min) & (self.Wavelength<self.correction_region_max)],4))
                    p_val = np.polyval(p,self.Wavelength[(self.Wavelength>self.correction_region_min) & (self.Wavelength<self.correction_region_max)])
                    self.plotted_absorbance[(self.Wavelength>self.correction_region_min) & (self.Wavelength<self.correction_region_max)] = self.Absorbance[(self.Wavelength>self.correction_region_min) & (self.Wavelength<self.correction_region_max)] - p_val
                elif self.ui.mean_radioButton.isChecked():
                    self.plotted_absorbance = self.Absorbance - np.mean(self.Absorbance[(self.Wavelength>self.correction_region_min) & (self.Wavelength<self.correction_region_max)])
            self.absorbance_plot.plot(self.Wavelength, self.plotted_absorbance, pen='r', clear=True)
            self.absorbance_plot.setYRange(0, 4)
            self.absorbance_plot.addItem(self.correction_region, ignoreBounds=True)
        except Exception as err:
            print(format(err))

    def clear_uvvis(self):
        self.absorbance_plot.clear()

    def plot_tauc(self):
        try:
            self.hv_min = self.ui.hv_min_spinBox.value()
            self.hv_max = self.ui.hv_max_spinBox.value()
            self.hv = 1240/self.Wavelength
            self.Alpha_hv = (self.plotted_absorbance * self.hv)**2.0
            self.index = (self.hv > self.hv_min) & (self.hv < self.hv_max)
            model = np.polyfit(self.hv[self.index], self.Alpha_hv[self.index], 1) 
            self.Alpha_hv_fit = self.hv * model[0] + model[1] #This is the linear fit
            self.tauc_plot.plot(self.hv, self.Alpha_hv, pen='r')
            self.tauc_plot.plot(self.hv, self.Alpha_hv_fit, pen='k')
            self.tauc_plot.setXRange(1,2)
            self.tauc_plot.setYRange(0, np.max(self.Alpha_hv[self.index]) + 1)

            self.Eg = - model[1]/model[0]
            self.ui.bandgap_label.setText(str(self.Eg))
        except:
            pass

    def clear_tauc(self):
        self.tauc_plot.clear()

    def export_uv_vis(self):
        """ Export publication ready uv vis figure """
        try:
            filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
            plt.figure(figsize=(8,6))
            plt.tick_params(direction='out', length=8, width=3.5)
            plt.plot(self.Wavelength, self.plotted_absorbance, linewidth = 3, color = 'r')
            if self.scatter_corrected:
                plt.xlim(self.correction_region_min, self.correction_region_max)
                plt.ylim(0, np.max(self.plotted_absorbance[(self.Wavelength>self.correction_region_min)]) +0.5)
            else:
                plt.xlim(self.correction_region_min, self.correction_region_max)
                plt.ylim(0, np.max(self.plotted_absorbance[(self.Wavelength>self.correction_region_min)] +0.5))
            plt.xlabel('Wavelength (nm)', fontsize = 20)
            plt.ylabel('Absorbance (a.u.)', fontsize = 20)
            plt.savefig(filename[0],bbox_inches='tight', dpi=300)
            plt.close()
        except:
            pass

    def export_tauc(self):
        """ Export publication ready tauc figure"""
        try:
            filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
            plt.figure(figsize=(8,6))
            plt.tick_params(direction='out', length=8, width=3.5)
            plt.plot(self.hv, self.Alpha_hv, linewidth = 3, color = 'r')
            plt.plot(self.hv, self.Alpha_hv_fit, linewidth = 2, color = 'black')
            plt.xlim(1,2)
            plt.ylim(0, np.max(self.Alpha_hv[self.index]) + 1)
            plt.xlabel('h$\\nu$ (eV)', fontsize = 20)
            plt.ylabel('($\\alpha$h$\\nu$)$^2$', fontsize = 20)
            #plt.title(Plot_title, fontsize = 20)

            plt.text(1.2, 1.2, r'E$_{g}$ = %.2f eV'%self.Eg, fontsize = 15)
            plt.tight_layout()
            plt.savefig(filename[0],bbox_inches='tight', dpi=300)
            plt.close()
        except:
            pass

"""Run the Main Window"""
def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win