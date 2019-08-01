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
pg.setConfigOption('imageAxisOrder', 'row-major')

base_path = Path(__file__).parent
file_path = (base_path / "flim_plot_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
	
	def __init__(self):
		super(TemplateBaseClass, self).__init__()
		
		# Create the main window
		self.ui = WindowTemplate()
		self.ui.setupUi(self)
		
		self.ui.load_scan_pushButton.clicked.connect(self.open_pkl_file)
		self.ui.plot_intensity_sums_pushButton.clicked.connect(self.plot_intensity_sums)
		self.ui.plot_raw_hist_data_pushButton.clicked.connect(self.plot_raw_scan)
		self.ui.save_intensities_image_pushButton.clicked.connect(self.save_intensities_image)
		self.ui.save_intensities_array_pushButton.clicked.connect(self.save_intensities_array)
		self.ui.compare_checkBox.stateChanged.connect(self.switch_compare)
		self.ui.intensity_sums_viewBox.roi.sigRegionChanged.connect(self.line_profile_update_plot)
		self.show()

	"""Open Scan Files"""
	def open_pkl_file(self):
		try:
			self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
			self.pkl_file = pickle.load(open(self.filename[0], 'rb'))
		except Exception as err:
			print(format(err))

	def plot_intensity_sums(self):
		try:
			data = self.pkl_file
			self.numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			self.numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			self.x_step_size = float(data['Scan Parameters']['X step size (um)'])
			self.x_scan_size = float(data['Scan Parameters']['X scan size (um)'])
			self.y_step_size = float(data['Scan Parameters']['Y step size (um)'])
			self.y_scan_size = float(data['Scan Parameters']['Y scan size (um)'])
			# TODO test line scan plots

			hist_data = data["Histogram data"]
			hist_data = np.reshape(hist_data, newshape=(hist_data.shape[0], self.numb_pixels_X*self.numb_pixels_Y))
			self.intensity_sums = np.sum(hist_data, axis=0)
			self.intensity_sums = np.reshape(self.intensity_sums, newshape=(self.numb_pixels_X, self.numb_pixels_Y))
			self.ui.intensity_sums_viewBox.view.invertY(False) # stop y axis invert
			self.ui.intensity_sums_viewBox.setImage(self.intensity_sums, scale=
												  (data['Scan Parameters']['X step size (um)'],
												   data['Scan Parameters']['Y step size (um)']))
			self.ui.intensity_sums_viewBox.roi.setSize([self.x_scan_size, self.y_step_size]) #line roi
			scale = pg.ScaleBar(size=1,suffix='um')
			scale.setParentItem(self.ui.intensity_sums_viewBox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
		except Exception as err:
			print(format(err))

	def line_profile_update_plot(self):
		if hasattr(self, "intensity_sums"):
			roiPlot = self.ui.intensity_sums_viewBox.getRoiPlot()
			roiPlot.clear()
			roi = self.ui.intensity_sums_viewBox.roi

			image = self.ui.intensity_sums_viewBox.getProcessedImage()

			# Extract image data from ROI
			axes = (self.ui.intensity_sums_viewBox.axes['x'], self.ui.intensity_sums_viewBox.axes['y'])
			data, coords = roi.getArrayRegion(image.view(np.ndarray), self.ui.intensity_sums_viewBox.imageItem, axes, returnMappedCoords=True)

			#calculate sums along columns in region
			sums_to_plot = np.sum(data, axis=0)

			#get scan x-coordinates in region
			x_values = coords[1][0]

			try:
				roiPlot.plot(x_values, sums_to_plot)
			except:
				pass

	def plot_raw_scan(self):
		try:
			data = self.pkl_file
			self.numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			self.numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			self.x_step_size = float(data['Scan Parameters']['X step size (um)'])
			self.x_scan_size = float(data['Scan Parameters']['X scan size (um)'])
			self.y_step_size = float(data['Scan Parameters']['Y step size (um)'])
			self.y_scan_size = float(data['Scan Parameters']['Y scan size (um)'])
			# TODO test line scan plots
			hist_data = data['Histogram data']
			
			self.hist_image = np.reshape(hist_data, newshape=(hist_data.shape[0],self.numb_pixels_X,self.numb_pixels_Y))
			time_data = data['Time data']
			self.times = time_data[:, 0, 0]*1e-3
			self.ui.raw_hist_data_viewBox.view.invertY(False) # stops y-axis invert
			self.ui.raw_hist_data_viewBox.setImage(self.hist_image, scale=
												(data['Scan Parameters']['X step size (um)'],
												data['Scan Parameters']['Y step size (um)']), xvals=self.times)
			self.ui.raw_hist_data_viewBox.roi.setSize([self.x_scan_size, self.y_scan_size])
			# if self.ui.compare_checkBox.isChecked():
			# 	self.ui.imv2.setImage(self.hist_image, scale= (data['Scan Parameters']['X step size (um)'],
			# 									data['Scan Parameters']['Y step size (um)']), xvals=self.times)
			self.switch_compare()
			self.ui.raw_hist_data_viewBox.ui.roiBtn.clicked.connect(self.switch_compare)
			scale = pg.ScaleBar(size=1,suffix='um')
			scale.setParentItem(self.ui.raw_hist_data_viewBox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
			
		except Exception as err:
			print(format(err))

	def switch_compare(self):
		"""
		Handles compare checkbox. If checked, show second ROI that user can use for comparison to first ROI.
		"""
		if self.ui.compare_checkBox.isChecked() and hasattr(self, "hist_image"):
			if not hasattr(self, "roi2"):
				self.roi2 = pg.ROI(pos=[0,0], size=[int(self.x_scan_size/2), int(self.y_scan_size/2)], movable=True, pen='r')
				self.roi2.addScaleHandle([1, 1], [0, 0])
				self.roi2.addRotateHandle([0, 0], [1, 1])
				self.roi2.sigRegionChanged.connect(self.update_roi2_plot)
				self.ui.raw_hist_data_viewBox.addItem(self.roi2)
				self.update_roi2_plot()
				self.roi2.hide()
				self.roi2_plot.hide()
			if self.ui.raw_hist_data_viewBox.ui.roiBtn.isChecked():
				self.roi2.show()
				self.roi2_plot.show()
			else:
				self.roi2.hide()
				self.roi2_plot.hide()
		else:
			if hasattr(self, "roi2"):
				self.roi2.hide()
				self.roi2_plot.hide()

	def update_roi2_plot(self):
		#Adapted from pyqtgraph imageview sourcecode
		
		image = self.ui.raw_hist_data_viewBox.getProcessedImage()

		# Extract image data from ROI
		axes = (self.ui.raw_hist_data_viewBox.axes['x'], self.ui.raw_hist_data_viewBox.axes['y'])
		data, coords = self.roi2.getArrayRegion(image.view(np.ndarray), self.ui.raw_hist_data_viewBox.imageItem, axes, returnMappedCoords=True)
		if data is None:
			return
		# Average data within entire ROI for each frame
		data = data.mean(axis=max(axes)).mean(axis=min(axes))
		xvals = self.ui.raw_hist_data_viewBox.tVals
		if hasattr(self, "roi2_plot"):
			self.roi2_plot.clear()
		self.roi2_plot = self.ui.raw_hist_data_viewBox.getRoiPlot().plot(xvals, data, pen='r')

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
			np.savetxt(save_to, self.intensity_sums.T, fmt='%f') #save transpoed intensity sums, as original array handles x in cols and y in rows
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