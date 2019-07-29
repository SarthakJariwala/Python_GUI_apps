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
			self.ui.intensity_sums_viewBox.roi.setSize([self.x_scan_size, self.y_scan_size])
			scale = pg.ScaleBar(size=1,suffix='um')
			scale.setParentItem(self.ui.intensity_sums_viewBox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
		except Exception as err:
			print(format(err))

	def custom_round(self, a, round_step):
		#Round a to the nearest interval of round_step
		return round(float(a)/round_step) * round_step

	def line_profile_update_plot(self):
		if self.ui.line_profile_checkBox.isChecked() and hasattr(self, "intensity_sums"):
			roiPlot = self.ui.intensity_sums_viewBox.getRoiPlot()
			roiPlot.clear()
			intensity_sums = self.intensity_sums
			roi = self.ui.intensity_sums_viewBox.roi

			##select relevant data in x direction by looking at roi
			roi_width = roi.size()[0]
			x_roi_min = self.custom_round(roi.pos()[0], self.x_step_size)
			x_roi_max = self.custom_round(roi.pos()[0] + roi_width, self.x_step_size)
			x_values = np.arange(start=x_roi_min, stop=x_roi_max+self.x_step_size, step=self.x_step_size) #get x values within x range, in x_step_size increments
			x_index_min = 0 
			if x_roi_min != 0: #to avoid divide by 0 error
				x_index_min = round(x_roi_min/self.x_step_size) -1
			x_index_max = round(x_roi_max/self.x_step_size)+1
			column_slice = intensity_sums[:, x_index_min:x_index_max] #get array of data in x range
		
			##select relevant data in y direction by looking at roi
			roi_height = roi.size()[1]
			y_roi_min = self.custom_round(roi.pos()[1], self.y_step_size)
			y_roi_max = self.custom_round(roi.pos()[1] + roi_height, self.y_step_size)
			y_range = [y_roi_min, y_roi_max] #get y range of roi
			y_index_min = 0
			if y_range[0] != 0: #to avoid divide by 0 error
				y_index_min = round(y_roi_min/self.y_step_size) -1 
			y_index_max = round(y_range[1]/self.y_step_size)+1
			row_slice = column_slice[y_index_min:y_index_max] #get array of data in y range

			#sum intensities along columns
			sums_to_plot = np.sum(row_slice, axis=0)

 			#even with rounding, x_values and sums_to_plot sometimes have slightly uneven lengths
 			#handle this issue here
			if x_values.shape[0] != sums_to_plot.shape[0]:
				newshape = min(x_values.shape[0], sums_to_plot.shape[0])
				x_values = x_values[0:newshape]
				sums_to_plot = sums_to_plot[0:newshape]
			
			try:
				roiPlot.plot(x_values, sums_to_plot)
			except:
				pass

	def plot_raw_scan(self):
		try:
			data = self.pkl_file
			self.numb_pixels_X = int((data['Scan Parameters']['X scan size (um)'])/(data['Scan Parameters']['X step size (um)']))
			self.numb_pixels_Y = int((data['Scan Parameters']['Y scan size (um)'])/(data['Scan Parameters']['Y step size (um)']))
			# TODO test line scan plots
			hist_data = data['Histogram data']
			
			self.hist_image = np.reshape(hist_data, newshape=(hist_data.shape[0],self.numb_pixels_X,self.numb_pixels_Y))
			time_data = data['Time data']
			self.times = time_data[:, 0, 0]*1e-3
			self.ui.raw_hist_data_viewBox.view.invertY(False) # stops y-axis invert
			self.ui.raw_hist_data_viewBox.setImage(self.hist_image, scale=
												(data['Scan Parameters']['X step size (um)'],
												data['Scan Parameters']['Y step size (um)']), xvals=self.times)
			if self.ui.compare_checkBox.isChecked():
				self.ui.imv2.setImage(self.hist_image, scale= (data['Scan Parameters']['X step size (um)'],
												data['Scan Parameters']['Y step size (um)']), xvals=self.times)
			scale = pg.ScaleBar(size=1,suffix='um')
			scale.setParentItem(self.ui.raw_hist_data_viewBox.view)
			scale.anchor((1, 1), (1, 1), offset=(-30, -30))
			
		except Exception as err:
			print(format(err))

	def switch_compare(self):
		"""
		Handles compare checkbox. If checked, show second ROI that user can use for comparison to first ROI.
		"""
		if self.ui.compare_checkBox.isChecked():
			self.imv2 = pg.ImageView()
			if hasattr(self, "hist_image"):
				self.imv2.setImage(self.hist_image, scale= (self.pkl_file['Scan Parameters']['X step size (um)'],
					self.pkl_file['Scan Parameters']['Y step size (um)']), xvals=self.times)
				self.imv2.view.invertY(False) # stop y-axis invert
			self.ui.gridLayout.addWidget(self.imv2, 10, 4)
		else:
			self.ui.gridLayout.removeWidget(self.imv2)
			self.imv2.hide()

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