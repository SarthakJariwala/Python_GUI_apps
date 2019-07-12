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
	
	"""Open Single Spectrum Files"""    
	# def open_file(self):
	# 	try:
	# 		filename = QtWidgets.QFileDialog.getOpenFileName(self)
	# 		try:
	# 			self.file = np.loadtxt(filename[0], skiprows = 16, delimiter='\t')
	# 		except:
	# 			self.file = np.genfromtxt(filename[0], skip_header=1, skip_footer=3, delimiter='\t')
	# 	except:
	# 		pass
		
	"""Open Scan Files"""
	def open_pkl_file(self):
		try:
			self.pkl_file = QtWidgets.QFileDialog.getOpenFileName(self, "*.pkl")
			#self.ui.result_textBrowser.append("Done Loading - .pkl")
		except:
			pass
	
	# def save_file(self):# not used yet!
	# 	try:
	# 		filename = QtWidgets.QFileDialog.getSaveFileName(self)
	# 		np.savetxt(filename[0], self.out, fmt = '%.5f', header = 'Time, Raw_PL, Sim_PL', delimiter = ' ')
	# 	except:
	# 		pass

	def plot_raw_scan(self):
		try:
			data = self.spec_scan_file
			numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			# TODO test line scan plots

			intensities = data['Intensities'].T
			
			intensities = np.reshape(intensities, newshape=(2048,numb_pixels_X,numb_pixels_Y))
			
			self.ui.raw_scan_viewbox.setImage(intensities, scale=
												  (data['Scan Parameters']['X step size (um)'],
												   data['Scan Parameters']['Y step size (um)']))
			
			scale = pg.ScaleBar(size=2,suffix='um')
			scale.setParentItem(self.ui.raw_scan_viewbox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
			
		except:
			pass

	def plot_intensity_sums(self):
		try:
			data = self.spec_scan_file
			numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			# TODO test line scan plots

			intensities = data['Intensities']

			#intensities = np.reshape(intensities, newshape=(2048, numb_pixels_X*numb_pixels_Y))
			
			sums = np.sum(intensities, axis=-1)
			sums = np.reshape(sums, newshape=(numb_pixels_X, numb_pixels_Y))
			print(sums)
			print(sums.shape)
			self.ui.intensity_sums_viewbox.setImage(sums, scale=
												  (data['Scan Parameters']['X step size (um)'],
												   data['Scan Parameters']['Y step size (um)']))
			
			scale = pg.ScaleBar(size=2,suffix='um')
			scale.setParentItem(self.ui.intensity_sums_viewbox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
		except:
			pass

	# def pub_ready_plot_export(self):
	# 	filename = QtWidgets.QFileDialog.getSaveFileName(self,caption="Filename with EXTENSION")
	# 	try:
	# 		try:
	# 			data = self.spec_scan_file
	# 			param_selection = str(self.ui.comboBox.currentText())
	# 			if param_selection == 'pk_pos': label = 'PL Peak Position (n.m.)'
	# 			elif param_selection == 'fwhm': label = 'PL FWHM (n.m.)'
	# 			cpm.plot_confocal(self.img, figsize=(10,10), stepsize = data['Scan Parameters']['X step size (um)'], cmap="seismic", cbar_label=label)
	# 			plt.savefig(filename[0],bbox_inches='tight', dpi=300)
	# 			plt.close()
	# 		except:
	# 			plt.figure(figsize=(8,6))
	# 			plt.tick_params(direction='out', length=8, width=3.5)
	# 			plt.plot(self.x, self.y)
	# 			plt.plot(self.x, self.result.best_fit,'k')
	# 			plt.xlabel("Wavelength (nm)", fontsize=20, fontweight='bold')
	# 			plt.ylabel("Intensity (a.u.)", fontsize=20, fontweight='bold')
	# 			plt.tight_layout()
				
	# 			plt.savefig(filename[0],bbox_inches='tight', dpi=300)
	# 			plt.close()
			
	# 	except AttributeError:
	# 		self.ui.result_textBrowser.setText("Need to fit the data first!")

	# #Get data from ocean optics scan pkl file, and convert to txt
	# def pkl_data_to_txt(self):
	# 	folder = os.path.dirname(self.pkl_to_convert[0])
	# 	filename_ext = os.path.basename(self.pkl_to_convert[0])
	# 	filename = os.path.splitext(filename_ext)[0] #get filename without extension
	# 	pkl_file = pickle.load(open(self.pkl_to_convert[0], 'rb'))

	# 	txt_file = np.zeros(shape=(2048,pkl_file['Intensities'].shape[0] + 1))

	# 	data_array = pkl_file['Intensities']
	# 	data_array = np.transpose(data_array)
	# 	wavelength = pkl_file['Wavelengths']

	# 	txt_file[:,0] = wavelength

	# 	for i in range(pkl_file['Intensities'].shape[0]):
	# 		txt_file[:,i+1] = data_array[:,i]

	# 	np.savetxt(folder +"/"+ filename +"_data.txt", txt_file, fmt = '%.2f', delimiter= "\t", header="wavelength(nm), Intensities at different points")
	# 	self.ui.result_textBrowser.append("Data from .pkl saved as .txt")

	# #Get scan parameters from ocean optics scan pkl file, and convert to txt
	# def pkl_params_to_txt(self):
	# 	folder = os.path.dirname(self.pkl_to_convert[0])
	# 	filename_ext = os.path.basename(self.pkl_to_convert[0])
	# 	filename = os.path.splitext(filename_ext)[0] #get filename without extension
	# 	pkl_file = pickle.load(open(self.pkl_to_convert[0], 'rb'))

	# 	pkl_scan = pkl_file['Scan Parameters']
	# 	pkl_oo = pkl_file['OceanOptics Parameters']
		
	# 	param_list = []
	# 	param_list.append(['X scan start (um)', 'Y scan start (um)', 'X scan size (um)', 'Y scan size (um)',
	# 		'X step size (um)', 'Y step size (um)', 'Integration Time (ms)', 'Scans to Average', 'Correct Dark Counts']) #list of param names
	# 	param_list.append([ pkl_scan['X scan start (um)'], pkl_scan['Y scan start (um)'], pkl_scan['X scan size (um)'],
	# 		pkl_scan['Y scan size (um)'], pkl_scan['X step size (um)'], pkl_scan['Y step size (um)'],
	# 		pkl_oo['Integration Time (ms)'], pkl_oo['Scans Averages'], pkl_oo['Correct Dark Counts'] ]) #list of param values

	# 	param_list = list(zip(*param_list)) #transpose so names and values are side-by-side
	# 	save_to = folder +"/"+ filename +"_scan_parameters.txt"
		
	# 	with open(save_to, 'w') as f:
	# 		for item in param_list:
	# 			f.write("%s\t" % str(item[0])) #write name
	# 			f.write("%s\n" % str(item[1])) #write value

	# 	self.ui.result_textBrowser.append("Scan parameters from .pkl saved as .txt")

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