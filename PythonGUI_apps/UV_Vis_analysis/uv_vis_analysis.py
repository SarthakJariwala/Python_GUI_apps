# system imports
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph import exporters
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
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
		
		#setup uv vis plot
		self.absorbance_plot_layout = pg.GraphicsLayoutWidget()
		self.ui.absorbance_plot_container.layout().addWidget(self.absorbance_plot_layout)
		self.absorbance_plot = self.absorbance_plot_layout.addPlot(title="Wavelengths vs. Absorbance")
		self.absorbance_plot.setLabel('bottom', 'Wavelength', unit='nm')
		self.absorbance_plot.setLabel('left', 'Absorbance', unit='a.u.')

		self.correction_region = pg.LinearRegionItem()
		self.correction_region.setZValue(10)
		
		self.ui.actionLoad_data.triggered.connect(self.open_data_file)
		self.ui.plot_absorbance_pushButton.clicked.connect(self.plot_absorbance)
		self.ui.correct_for_scattering_checkBox.stateChanged.connect(self.switch_correction_region)
		self.correction_region.sigRegionChanged.connect(self.update_correction_region)

		#setup tauc plot
		self.tauc_plot_layout = pg.GraphicsLayoutWidget()
		self.ui.tauc_plot_container.layout().addWidget(self.tauc_plot_layout)
		self.tauc_plot = self.tauc_plot_layout.addPlot(title="Tauc plot fit")
		self.tauc_plot.setLabel('bottom', 'hv', unit='ev')
		y_label = '(ahv)' + chr(0x00B2)
		self.tauc_plot.setLabel('left', y_label)

		self.ui.plot_tauc_pushButton.clicked.connect(self.plot_tauc)
		self.ui.export_tauc_pushButton.clicked.connect(self.export_tauc)

		self.show()

	"""Open Scan Files"""
	def open_data_file(self):
		try:
			self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
			self.data = np.loadtxt(self.filename[0], delimiter = ',', skiprows = 1)
			self.Wavelength = self.data[:,0] # in nm
			self.Absorbance = self.data[:,1] 
			self.correction_region.setRegion((np.min(self.Wavelength), np.max(self.Wavelength)))
		except Exception as err:
			print(format(err))

	def update_correction_region(self):
		self.correction_region_min, self.correction_region_max = self.correction_region.getRegion()

	def switch_correction_region(self):
		if self.ui.correct_for_scattering_checkBox.isChecked():
			self.absorbance_plot.addItem(self.correction_region, ignoreBounds=True)
			self.correction_region.show()
		else:
			self.correction_region.hide()

	def plot_absorbance(self):
		try:
			Absorbance = self.Absorbance
			if self.ui.correct_for_scattering_checkBox.isChecked() and hasattr(self, "correction_region_min"):
				Absorbance = Absorbance - np.mean(Absorbance[(self.Wavelength>self.correction_region_min) & (self.Wavelength<self.correction_region_max)])
			self.absorbance_plot.plot(self.Wavelength, Absorbance, pen='r', clear=True)
			if self.ui.correct_for_scattering_checkBox.isChecked():
				self.absorbance_plot.addItem(self.correction_region, ignoreBounds=True)
				self.correction_region.show()
		except Exception as err:
			print(format(err))

	def plot_tauc(self):
		try:
			hv_min = self.ui.hv_min_spinBox.value()
			hv_max = self.ui.hv_max_spinBox.value()
			hv = 1240/self.Wavelength
			Alpha_hv = (self.Absorbance * hv)**2.0
			index = (hv > hv_min) & (hv < hv_max)
			model = np.polyfit(hv[index], Alpha_hv[index], 1) 
			Alpha_hv_fit = hv * model[0] + model[1] #This is the linear fit
			self.tauc_plot.plot(hv, Alpha_hv, pen='r')
			self.tauc_plot.plot(hv, Alpha_hv_fit, pen='k')
		except:
			pass

	def export_tauc(self):
		try:
			filename_ext = os.path.basename(self.filename[0])
			filename = os.path.splitext(filename_ext)[0] #get filename without extension
			save_to = os.getcwd() + "\\" + filename + "_TaucPlot.tiff"

			exporter = pg.exporters.ImageExporter(self.tauc_plot)
			exporter.params.param('width').setValue(800, blockSignal=exporter.widthChanged)
			exporter.params.param('height').setValue(600, blockSignal=exporter.heightChanged)
			# save to file
			exporter.export(save_to)
		except Exception as err:
			print(format(err))

"""Run the Main Window"""
def run():
	win = MainWindow()
	QtGui.QApplication.instance().exec_()
	return win