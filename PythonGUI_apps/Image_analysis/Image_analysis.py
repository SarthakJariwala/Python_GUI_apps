import sys
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

# local modules
 
pg.mkQApp()

base_path = Path(__file__).parent
file_path = (base_path / "image_analysis_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  

	def __init__(self):
		super(TemplateBaseClass, self).__init__()
		
		# Create the main window
		self.ui = WindowTemplate()
		self.ui.setupUi(self)

		#setup imageview
		self.imv = pg.ImageView()
		self.imv.getView().setAspectLocked(lock=False, ratio=1)
		self.imv.getView().setMouseEnabled(x=True, y=True)
		self.imv.getView().invertY(False)
		self.roi = self.imv.roi
		self.roi.translateSnap = True
		self.roi.scaleSnap = True
		self.update_camera() #initialize camera pixel size
		self.update_scaling_factor() #initialize scaling_factor

		self.roi_plot = self.imv.getRoiPlot().getPlotItem() #get roi plot
		self.ui.image_groupBox.layout().addWidget(self.imv)

		#set up ui signals
		self.roi.sigRegionChanged.connect(self.line_profile_update_plot)
		self.ui.load_image_pushButton.clicked.connect(self.load_image)
		self.ui.custom_pixel_size_checkBox.stateChanged.connect(self.switch_custom_pixel_size)
		self.ui.update_scaling_factor_pushButton.clicked.connect(self.update_scaling_factor)
		self.ui.spot_radioButton.toggled.connect(self.update_camera)

		self.num_plots = 0
		self.show()

	def load_image(self):
		"""
		Prompts the user to select a text file containing image data.
		"""
		try:
			file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', os.getcwd())
			self.original_image = Image.open(file[0])
			#self.original_image = self.original_image.rotate(-90, expand=True)
			self.resize_to_scaling_factor(self.original_image)
		except Exception as err:
			print(format(err))

	def resize_to_scaling_factor(self, image):
		"""
		Handles loading of image according to scaling_factor
		"""
		image = image.resize((round(image.size[0]*self.scaling_factor), round(image.size[1]*self.scaling_factor)))
		if self.ui.greyscale_checkBox.isChecked():
			image = image.convert("L") #convert to greyscale
		image_array = np.asarray(image).T #correct numpy array auto-flip
		width = image_array.shape[0]
		height = image_array.shape[1]
		try:
			x_vals = np.arange(width) #imv x-axis
			self.imv.setImage(img=image_array, xvals= x_vals)
			self.roi.setPos((0,0))
			self.roi.setSize([width, height * self.scaling_factor]) #set line roi
			self.line_profile_update_plot()
		except:
			pass

	def line_profile_update_plot(self):
		""" Handle line profile for intensity sum viewbox """
		self.roi_plot.clear()
		image = self.imv.getProcessedImage()

		# Extract image data from ROI
		axes = (self.imv.axes['x'], self.imv.axes['y'])
		data, coords = self.roi.getArrayRegion(image.view(np.ndarray), self.imv.imageItem, axes, returnMappedCoords=True)
		if data is None:
			return

		x_values = coords[1][0]

		#calculate sums along columns in region
		if len(data.shape) == 2: #if grayscale, average intensities 
			sums_to_plot = np.mean(data, axis=0)
			try:
				self.roi_plot.plot(x_values, sums_to_plot)
			except:
				pass
		elif len(data.shape) > 2: #if rgb arrays, plot individual components
			r_values = data[:,:,0]
			g_values = data[:,:,1]
			b_values = data[:,:,2]
			r_avg = np.mean(r_values, axis=0) #average red values across columns
			g_avg = np.mean(g_values, axis=0) #average green values
			b_avg = np.mean(b_values, axis=0) #average blue values
			try:
				self.roi_plot.plot(x_values, r_avg, pen='r')
				self.roi_plot.plot(x_values, g_avg, pen='g')
				self.roi_plot.plot(x_values, b_avg, pen='b')
			except:
				pass

	def update_scaling_factor(self):
		"""
		Calculate scaling factor
		"""
		if self.ui.custom_pixel_size_checkBox.isChecked():
			self.scaling_factor = self.ui.custom_pixel_size_spinBox.value()
		else:
			self.scaling_factor = self.camera_pixel_size/int(self.ui.magnification_comboBox.currentText())
		self.roi.snapSize = self.scaling_factor #roi snaps to multiples of scaling_factor
		if hasattr(self, "original_image"):
			self.resize_to_scaling_factor(self.original_image) #resize image, sets up roi

	def switch_custom_pixel_size(self):
		checked = self.ui.custom_pixel_size_checkBox.isChecked()
		self.ui.custom_pixel_size_spinBox.setEnabled(checked)
		self.ui.magnification_comboBox.setEnabled(not checked)

	def update_camera(self):
		if self.ui.spot_radioButton.isChecked():
			self.camera_pixel_size = 7.4
		elif self.ui.pixera_radioButton.isChecked():
			self.camera_pixel_size = 3

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