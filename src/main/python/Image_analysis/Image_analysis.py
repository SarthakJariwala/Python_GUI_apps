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

def updateDelay(scale, time):
	""" Hack fix for scalebar inaccuracy"""
	QtCore.QTimer.singleShot(time, scale.updateBar)

class MainWindow(TemplateBaseClass):  

	def __init__(self):
		pg.setConfigOption('imageAxisOrder', 'col-major') 
		super(TemplateBaseClass, self).__init__()
		
		# Create the main window
		self.ui = WindowTemplate()
		self.ui.setupUi(self)

		#setup image plot
		self.image_plot_layout=pg.GraphicsLayoutWidget()
		self.ui.image_groupBox.layout().addWidget(self.image_plot_layout)
		self.image_plot = self.image_plot_layout.addPlot()
		self.img_item = pg.ImageItem()
		self.image_plot.addItem(self.img_item)
		self.image_plot_view = self.image_plot.getViewBox()

		#setup lookup table
		self.hist_lut = pg.HistogramLUTItem()
		self.image_plot_layout.addItem(self.hist_lut)

		#region of interest - allows user to select scan area
		self.roi = pg.ROI([0,0],[10, 10], movable=True)
		self.roi.addScaleHandle([1, 1], [0, 0])
		self.roi.addRotateHandle([0, 0], [1, 1])
		self.roi.translateSnap = True
		self.roi.scaleSnap = True
		self.roi.sigRegionChanged.connect(self.line_profile_update_plot)
		self.image_plot.addItem(self.roi)

		#setup rgb plot
		self.rgb_plot_layout=pg.GraphicsLayoutWidget()
		self.ui.rgb_plot_groupBox.layout().addWidget(self.rgb_plot_layout)
		self.rgb_plot = self.rgb_plot_layout.addPlot()

		#set up ui signals
		self.ui.load_image_pushButton.clicked.connect(self.load_image)
		self.ui.custom_pixel_size_checkBox.stateChanged.connect(self.switch_custom_pixel_size)
		self.ui.update_settings_pushButton.clicked.connect(self.reload_image)
		self.ui.spot_radioButton.toggled.connect(self.update_camera)

		self.update_camera() #initialize camera pixel size
		self.update_scaling_factor() #initialize scaling_factor
		self.show()

		#row major. invert y false, rotate false
	def load_image(self):
		"""
		Prompts the user to select a text file containing image data.
		"""
		try:
			file = QtWidgets.QFileDialog.getOpenFileName(self, 'Open file', os.getcwd())
			self.original_image = Image.open(file[0])
			self.original_image = self.original_image.rotate(-90, expand=True) #correct image orientation
			self.resize_to_scaling_factor(self.original_image)
		except Exception as err:
			print(format(err))

	def resize_to_scaling_factor(self, image):
		"""
		Handles loading of image according to scaling_factor
		"""
		self.update_scaling_factor()

		if self.ui.spot_radioButton.isChecked() and self.ui.resize_image_checkBox.isChecked():
			image = self.original_image.resize((round(image.size[0]*self.scaling_factor), round(image.size[1]*self.scaling_factor)))
			self.image_plot.getAxis("bottom").setScale(scale = 1)
			self.image_plot.getAxis("left").setScale(scale = 1)
		else:
			image = self.original_image
			self.image_plot.getAxis("bottom").setScale(scale = self.scaling_factor)
			self.image_plot.getAxis("left").setScale(scale = self.scaling_factor)
			
		if self.ui.greyscale_checkBox.isChecked():
			image = image.convert("L") #convert to greyscale

		self.image_array = np.array(image)
		width = self.image_array.shape[0]
		height = self.image_array.shape[1]
		
		try:
			#set image bounds with qrect
			self.img_item_rect = QtCore.QRectF(0, 0, width, height)
			self.img_item.setImage(image=self.image_array)
			self.img_item.setRect(self.img_item_rect)

			# if self.ui.greyscale_checkBox.isChecked():
			# 	self.hist_lut.setImageItem(self.img_item)
			
			if self.ui.vertical_radioButton.isChecked():
				roi_height = self.scaling_factor * height
				self.roi.setSize([width, roi_height])
			elif self.ui.horizontal_radioButton.isChecked():
				roi_height = self.scaling_factor * width
				self.roi.setSize([roi_height, height])
			self.roi.setAngle(0)
			self.roi.setPos(0, 0)
			self.line_profile_update_plot()
		except:
			pass

	def line_profile_update_plot(self):
		""" Handle line profile for intensity sum viewbox """
		self.rgb_plot.clear()

		# Extract image data from ROI
		data, coords = self.roi.getArrayRegion(self.image_array, self.img_item, returnMappedCoords=True)
		if data is None:
			return

		if self.ui.vertical_radioButton.isChecked():
			x_values = coords[0,:,0]
		elif self.ui.horizontal_radioButton.isChecked():
			x_values = coords[1,0,:]

		if self.ui.pixera_radioButton.isChecked() or (self.ui.spot_radioButton.isChecked() and not self.ui.resize_image_checkBox.isChecked()):
			x_values = x_values * self.scaling_factor
		
		#calculate average along columns in region
		if len(data.shape) <= 2: #if grayscale, average intensities 
			if self.ui.vertical_radioButton.isChecked():
				avg_to_plot = np.mean(data, axis=-1)
			elif self.ui.horizontal_radioButton.isChecked():
				avg_to_plot = np.mean(data, axis=0)
			try:
				self.rgb_plot.plot(x_values, avg_to_plot)
			except:
				pass
		elif len(data.shape) > 2: #if rgb arrays, plot individual components
			r_values = data[:,:,0]
			g_values = data[:,:,1]
			b_values = data[:,:,2]
			if self.ui.vertical_radioButton.isChecked():
				r_avg = np.mean(r_values, axis=-1) #average red values across columns
				g_avg = np.mean(g_values, axis=-1) #average green values
				b_avg = np.mean(b_values, axis=-1) #average blue values
			elif self.ui.horizontal_radioButton.isChecked():
				r_avg = np.mean(r_values, axis=0)
				g_avg = np.mean(g_values, axis=0)
				b_avg = np.mean(b_values, axis=0)
			try:
				self.rgb_plot.plot(x_values, r_avg, pen='r')
				self.rgb_plot.plot(x_values, g_avg, pen='g')
				self.rgb_plot.plot(x_values, b_avg, pen='b')
			except Exception as e:
				pass

	def update_scaling_factor(self):
		"""
		Calculate scaling factor
		"""
		if self.ui.custom_pixel_size_checkBox.isChecked():
			self.camera_pixel_size = self.ui.custom_pixel_size_spinBox.value()
			self.scaling_factor = self.camera_pixel_size
		else:
			self.scaling_factor = self.camera_pixel_size/int(self.ui.magnification_comboBox.currentText())
		self.roi.snapSize = self.scaling_factor #roi snaps to multiples of scaling_factor

	def reload_image(self):
		if hasattr(self, "original_image"):
			self.resize_to_scaling_factor(self.original_image) #resize image, sets up roi

	def switch_custom_pixel_size(self):
		checked = self.ui.custom_pixel_size_checkBox.isChecked()
		self.ui.custom_pixel_size_spinBox.setEnabled(checked)
		self.ui.magnification_comboBox.setEnabled(not checked)

	def update_camera(self):
		if self.ui.spot_radioButton.isChecked():
			self.camera_pixel_size = 7.4
			self.ui.greyscale_checkBox.setChecked(False)
			self.ui.resize_image_checkBox.setEnabled(True)
			self.update_scaling_factor()
		elif self.ui.pixera_radioButton.isChecked():
			self.camera_pixel_size = 3
			self.ui.greyscale_checkBox.setChecked(True)
			self.ui.resize_image_checkBox.setEnabled(False)
			self.update_scaling_factor()

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