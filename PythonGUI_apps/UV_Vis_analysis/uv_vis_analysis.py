# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

# system imports
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

pg.mkQApp()
pg.setConfigOption('background', 'w')

base_path = Path(__file__).parent
file_path = (base_path / "uv_vis_analysis_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
	
	def __init__(self):
		super(TemplateBaseClass, self).__init__()
		
		# Create the main window
		self.ui = WindowTemplate()
		self.ui.setupUi(self)
		
		self.ui.actionLoad_data.triggered.connect(self.open_data_file)
		self.ui.plot_absorbance_pushButton.clicked.connect(self.plot_absorbance)
		self.ui.plot_tauc_pushButton.clicked.connect(self.plot_tauc)
		self.ui.export_tauc_pushButton.clicked.connect(self.export_tauc)

		self.ui.correct_for_scattering_checkBox.stateChanged.connect(self.switch_correct_for_scattering)


		self.absorbance_plot_layout = pg.GraphicsLayoutWidget()
		self.ui.absorbance_plot_container.layout().addWidget(self.absorbance_plot_layout)
		self.absorbance_plot = self.absorbance_plot_layout.addPlot(title="Wavelengths vs. Absorbance")
		self.absorbance_plot.setLabel('left', 'Wavelength', unit='nm')
		self.absorbance_plot.setLabel('bottom', 'Absorbance', unit='a.u.')

		self.tauc_plot_layout = pg.GraphicsLayoutWidget()
		self.ui.tauc_plot_container.layout().addWidget(self.tauc_plot_layout)
		self.tauc_plot = self.tauc_plot_layout.addPlot(title="Tauc plot fit")
		self.tauc_plot.setLabel('left', '($\\alpha$h$\\nu$)$^2$') #fix formatting
		self.tauc_plot.setLabel('bottom', 'h$\\nu', unit='ev')

		self.show()

	"""Open Scan Files"""
	def open_data_file(self):
		try:
			self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
			self.data = np.loadtxt(self.filename[0], delimiter = ',', skiprows = 1)
			self.Wavelength = self.data[:,0] # in nm
			self.Absorbance = self.data[:,1] 
		except Exception as err:
			print(format(err))

	def switch_correct_for_scattering(self):
		uncorrected = self.data[:,1]
		if self.ui.correct_for_scattering_checkBox.isChecked():
			self.Absorbance = uncorrected - np.mean(uncorrected[(self.Wavelength>800) & (self.Wavelength<1000)])
		else:
			self.Absorbance = uncorrected

	def plot_absorbance(self):
		self.absorbance_plot.plot(self.Wavelength, self.Absorbance, pen='r', clear=True)

	def plot_tauc(self):
		hv = 1240/Wavelength
		hv_min = self.ui.tauc_start_spinBox.value()
		hv_max = self.ui.tauc_end_spinBox.value()
		index = (hv > hv_min) & (hv < hv_max)
		Alpha_hv = (self.Absorbance * hv)**2.0
		model = np.polyfit(hv[index], Alpha_hv[index], 1) 
		
		self.tauc_plot.plot(hv, Alpha_hv, color = 'r')
		self.tauc_plot.plot(hv, Alpha_hv_fit, color = 'k')

	def export_tauc(self):
		filename_ext = os.path.basename(self.filename[0])
		filename = os.path.splitext(filename_ext)[0] #get filename without extension
		save_to = os.getcwd() + "\\" + filename + "_TaucPlot.tiff"
		cpm.plt.savefig(save_to, bbox_inches='tight', dpi = 300)

"""Run the Main Window"""    
def run():
	win = MainWindow()
	QtGui.QApplication.instance().exec_()
	return win