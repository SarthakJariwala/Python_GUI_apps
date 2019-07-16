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

		self.show()
	
	"""Open Scan Files"""
	def open_data_file(self):
		try:
			self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
			self.data_file = np.loadtxt(self.filename[0])
		except Exception as err:
			print(format(err))

	def plot_intensity_sums(self):
		try:
			data = self.pkl_file
			numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			# TODO test line scan plots

			hist_data = data["Histogram data"]
			hist_data = np.reshape(hist_data, newshape=(hist_data.shape[0], numb_pixels_X*numb_pixels_Y))
			self.intensity_sums = np.sum(hist_data, axis=0)
			self.intensity_sums = np.reshape(self.intensity_sums, newshape=(numb_pixels_X, numb_pixels_Y))
			self.ui.intensity_sums_viewBox.setImage(self.intensity_sums, scale=
												  (data['Scan Parameters']['X step size (um)'],
												   data['Scan Parameters']['Y step size (um)']))
			scale = pg.ScaleBar(size=2,suffix='um')
			scale.setParentItem(self.ui.intensity_sums_viewBox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
		except Exception as err:
			print(format(err))

	def plot_raw_scan(self):
		try:
			data = self.pkl_file
			numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			# TODO test line scan plots
			hist_data = data['Histogram data']
			
			hist_data = np.reshape(hist_data, newshape=(hist_data.shape[0],numb_pixels_X,numb_pixels_Y))
			
			time_data = data['Time data']
			times = time_data[:, 0, 0]*1e-3
			self.ui.raw_hist_data_viewBox.setImage(hist_data, scale=
												  (data['Scan Parameters']['X step size (um)'],
												   data['Scan Parameters']['Y step size (um)']), xvals=times)
			# time_data = data['Time data']
			# print(time_data.shape)
			# print(time_data)
			# roiPlot = self.ui.raw_hist_data_viewBox.getRoiPlot()
			# roiPlot.tVals = time_data#(np.min(time_data), np.max(time_data))
			
			scale = pg.ScaleBar(size=2,suffix='um')
			scale.setParentItem(self.ui.raw_hist_data_viewBox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
			
		except Exception as err:
			print(format(err))

	def save_intensities_image(self):
		try:
			filename_ext = os.path.basename(self.filename[0])
			filename = os.path.splitext(filename_ext)[0] #get filename without extension
			save_to = os.getcwd() + "\\" + filename + "_intensity_sums.png"
			cpm.plot_confocal(self.intensity_sums, stepsize=np.abs(self.pkl_file['Scan Parameters']['X step size (um)']))
			cpm.plt.savefig(save_to, bbox_inches='tight', dpi=300)
		except:
			pass

	def save_intensities_array(self):
		try:
			filename_ext = os.path.basename(self.filename[0])
			filename = os.path.splitext(filename_ext)[0] #get filename without extension
			save_to = os.getcwd() + "\\" + filename + "_intensity_sums.txt"
			np.savetxt(save_to, self.intensity_sums, fmt='%f')
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

"""Run the Main Window"""    
def run():
	win = MainWindow()
	QtGui.QApplication.instance().exec_()
	return win

#Uncomment below if you want to run this as standalone
#run()