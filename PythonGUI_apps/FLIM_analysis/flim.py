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
file_path = (base_path / "flim_analysis.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
	
	def __init__(self):
		super(TemplateBaseClass, self).__init__()
		
		# Create the main window
		self.ui = WindowTemplate()
		self.ui.setupUi(self)
		
		# self.file = None
		# self.bck_file = None
		# self.wlref_file = None
		# self.x = None
		# self.y = None
		# self.out = None # output file after fitting
		
		# Peak parameters if adjust params is selected		
		self.show()
	
	"""Open Scan Files"""
	def open_pkl_file(self):
		try:
			self.pkl_to_convert = QtWidgets.QFileDialog.getOpenFileName(self, "*.pkl")
			self.pkl_file = pickle.load(open(self.pkl_to_convert[0], 'rb'))
			self.ui.result_textBrowser.append("Done Loading - .pkl to convert")
		except:
			pass
	
	def plot_raw_scan(self):
		try:
			data = self.pkl_file
			numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			# TODO test line scan plots

			hist_data = data['Histogram Data'].T
			
			hist_data = np.reshape(hist_data, newshape=(2048,numb_pixels_X,numb_pixels_Y))
			
			self.ui.raw_hist_data_viewbox.setImage(hist_data, scale=
												  (data['Scan Parameters']['X step size (um)'],
												   data['Scan Parameters']['Y step size (um)']))
			
			scale = pg.ScaleBar(size=2,suffix='um')
			scale.setParentItem(self.ui.raw_hist_data_viewbox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
			
		except:
			pass

	def plot_intensity_sums(self):
		try:
			data = self.pkl_file
			numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			# TODO test line scan plots

			hist_data = data['Histogram Data']

			#intensities = np.reshape(intensities, newshape=(2048, numb_pixels_X*numb_pixels_Y))
			
			sums = np.sum(hist_data, axis=-1)
			sums = np.reshape(sums, newshape=(numb_pixels_X, numb_pixels_Y))
			self.ui.intensity_sums_viewbox.setImage(sums, scale=
												  (data['Scan Parameters']['X step size (um)'],
												   data['Scan Parameters']['Y step size (um)']))
			
			scale = pg.ScaleBar(size=2,suffix='um')
			scale.setParentItem(self.ui.intensity_sums_viewbox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
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
run()