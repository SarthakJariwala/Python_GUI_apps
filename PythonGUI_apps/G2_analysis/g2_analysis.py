"""
this app has not been added to the main databrowser as
some features do not work as they should

see blocks of code under ###comments
"""


# system imports
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph import exporters
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import matplotlib.pyplot as plt
import numpy as np
import read_ptu
import pandas as pd
import sys
from helper_funcs import refresh_text_box, find_nearest


pg.mkQApp()
pg.setConfigOption("background", "w")

base_path = Path(__file__).parent
file_path = (base_path / "g2_analysis_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)



class MainWindow(TemplateBaseClass):  
	
	def __init__(self):
		"""
		Sets up window ui 
		"""
		super(TemplateBaseClass, self).__init__()
		
		# Create the main window
		self.ui = WindowTemplate()
		self.ui.setupUi(self)

		#setup time trace plot
		self.time_counts_plot = self.ui.counts_plotWidget.getPlotItem()
		self.time_counts_plot.enableAutoRange(enable=True)
		self.time_counts_plot.setTitle(title="Time vs. Occurrences")
		self.time_counts_plot.setLabel("bottom", "Time", unit="s")
		self.time_counts_plot.setLabel("left", "Occurrences")
		
		#setup g2 plot
		self.tau_g2_plot = self.ui.g2_plotWidget.getPlotItem()
		self.tau_g2_plot.enableAutoRange(enable=True)
		self.tau_g2_plot.setTitle(title="Time vs. g2(t)")
		self.tau_g2_plot.setLabel("bottom", "Time", unit="ps")
		self.tau_g2_plot.setLabel("left", "g2(t)")
	
		#setup g2 plot integration region selection
		self.integration_region = pg.LinearRegionItem(brush=QtGui.QBrush \
			(QtGui.QColor(255, 0, 0, 50)))
		self.integration_region.sigRegionChanged.connect \
			(self.update_integration_spinBoxes)
		
		#setup ui signals
		self.ui.load_data_pushButton.clicked.connect(self.open_data_file)
		self.ui.plot_pushButton.clicked.connect(self.on_plot)
		self.ui.calculate_integration_pushButton.clicked.connect \
			(self.calculate_integration)
		
		self.show()

	def open_data_file(self):
		"""
		Load ptu file and extract timetag data.
		"""
		try:
			self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
			self.data = None #this will hold data from ptu file
			ptu_txt = "" #this will hold data from 
			if ".ptu" in self.filename[0][-4:]:
				read_ptu.ptufile = self.filename[0]
				read_ptu.setup_conversion() #open file and read settings
				self.glob_res = read_ptu.globRes
				self.num_records = read_ptu.numRecords
				
				#get readable file - ie. convert ptu to txt 
				ptu_txt = read_ptu.generate_txt(self.ui.textBrowser) 
			else:
				refresh_text_box(self.ui.textBrowser, "Invalid file.\n")
			
			refresh_text_box(self.ui.textBrowser, "Retrieving dataframe...")
			self.data = pd.read_csv(ptu_txt, delim_whitespace=True, \
				encoding="utf-16le", skiprows=77) #read txt file into array
			refresh_text_box(self.ui.textBrowser, "Dataframe acquired.\n")
			self.chan0_data = self.data.loc[self.data['chan'] == 0] \
				['truetime/ps'].values
			self.chan1_data = self.data.loc[self.data['chan'] == 1] \
				['truetime/ps'].values
		except Exception as err:
			refresh_text_box(self.ui.textBrowser, format(err) + "\n")
			
	def update_integration_spinBoxes(self):
		"""
		Update laser spinboxes based on line rois
		"""
		self.integration_start, self.integration_stop = \
			self.integration_region.getRegion()
		self.ui.start_integration_doubleSpinBox.setValue(self.integration_start)
		self.ui.stop_integration_doubleSpinBox.setValue(self.integration_stop)

	def update_integration_region(self):
		"""
		Update laser line rois based on spinboxes
		"""
		self.integration_start = self.ui.start_integration_doubleSpinBox.value()
		self.integration_stop = self.ui.stop_integration_doubleSpinBox.value()
		self.integration_region.setRegion((integration_start, integration_stop))

					
	def calculate_timecorr(self):
		#user input bin size: ms to s
		self.counts_bin_size = self.ui.counts_bin_size_doubleSpinBox.value() * 1e-3
		
		#store time_tags
		both_time_tags = np.array(self.data['truetime/ps'])
		chan0_time_tags = np.array(self.chan0_data)
		chan1_time_tags = np.array(self.chan1_data)

		#generate even-width bins
		num_bins = np.abs(int((np.max(both_time_tags * 1e-12) - np.min(both_time_tags * 1e-12)) / self.counts_bin_size))
		
		#generate histograms
		self.counts_hist = np.histogram(both_time_tags, num_bins)
		self.chan0_counts_hist = np.histogram(chan0_time_tags, num_bins)
		self.chan1_counts_hist = np.histogram(chan1_time_tags, num_bins)
		
		#get min and max counts from ui
		self.min_counts = self.ui.min_counts_spinBox.value()
		self.max_counts = self.ui.max_counts_spinBox.value()
		
	def plot_timecorr(self):
		"""
		Plot time correlation event of detected photons
		"""
		self.calculate_timecorr()
		
		#adjust y view range based on user input min and max
		self.time_counts_plot.setYRange(self.min_counts, self.max_counts)
		
		#plot based on selected radio button
		try:
			if self.ui.both_channels_radioButton.isChecked():
				self.time_counts_plot.plot(self.counts_hist[1] * 1e-12, self.counts_hist[0], stepMode=True, fillLevel=0, \
					fillOutline=True, brush=(0,0,255,150))
			elif self.ui.channel0_radioButton.isChecked():
				self.time_counts_plot.plot(self.chan0_counts_hist[1] * 1e-12, self.chan0_counts_hist[0], stepMode=True, fillLevel=0, \
					fillOutline=True, brush=(0,0,255,150))
			elif self.ui.channel1_radioButton.isChecked():
				self.time_counts_plot.plot(self.chan1_counts_hist[1] * 1e-12, self.chan1_counts_hist[0], stepMode=True, fillLevel=0, \
					fillOutline=True, brush=(0,0,255,150))
		except Exception as err:
			refresh_text_box(self.ui.textBrowser, format(err) + "\n")
			
	def calc_corr(self, first_channel_data, second_channel_data):
		"""
		Calculate correlation values in relation to lag times.
		Algorithm based on Wahl, Gregor, et al. Fast calculation 
		of fluorescence correlation data with async TCSPC. Section 2.4.
		
		Procedure:
		1. Traverse first channel until we find greater value than 
			first value in second channel. Add 1 to correlation value.
		2. Traverse second channel until we find greater value than 
			what we found in step 1. Add 1 to correlation value.
		3. Starting with value in step 1, traverse first channel until 
			we find greater value than what we found in step 2. 
			Add 1 to correlation value.
		4. Starting with value in step 2, raverse second channel until 
			we find greater value than what we found in step 3.
			Add 1 to correlation value/
		5. And so on... continue until we reach the end of either 
			channel array.
			
		:param first_channel_data: time tags of first channel to start traversal with
		:type first_channel_data: array-like
		
		:param second_channel_data: time tags of second channel in correlation calculation
		:type second_channel_data: array-like
		
		:returns: correlation array in format (n, 2)
			first column is lag time calculated by subtracting current first channel 
			point from corresponding second channel point
			second column is the correlation value
		:rtype: numpy.ndarray
		"""
		
		first_channel_pos = 0
		second_channel_pos = 0
		corr = 0
		no_val_found = 0
		corr_array = []
		while first_channel_pos < len(first_channel_data) and second_channel_pos < len(second_channel_data):

			#run this loop until we find the next value in channel 0 
			#greater than second_channel_pos
			for i in range(first_channel_pos, len(first_channel_data)):
				if first_channel_data[i] >= second_channel_data[second_channel_pos]:
					corr += 1
					
					#save index that holds greater value for first channel traversal
					first_channel_pos = i 
					
					#record corresponding lag time 
					#along with current correlation value
					try:
						corr_array.append([second_channel_data[i] - first_channel_data[i] , corr])
					except:
						break
					
					break 
				
				else:
					no_val_found = i
					
			#if last index where no greater value was found is end of array,
			#leave for loop, as we have satisfied the while loop condition
			if no_val_found > len(first_channel_data) - 1 or no_val_found > len(second_channel_data) - 1: 
				break

			#same as above, but traverse through channel 1
			for j in range(second_channel_pos, len(second_channel_data)):
				if second_channel_data[j] >= first_channel_data[first_channel_pos]:
					corr += 1

					second_channel_pos = j
					
					try:
						corr_array.append([second_channel_data[j] - first_channel_data[j] , corr])
					except:
						break
						
					break
					
				else:
					no_val_found = j
					
			#if last index where no greater value was found is end of array,
			#leave for loop, as we have satisfied the while loop condition
			if no_val_found > len(second_channel_data) - 1 or no_val_found > len(first_channel_data) - 1:
				break
				
		return np.array(corr_array)
	
	def plot_g2(self):
		"""
		Plot lag times vs. correlations.
		"""
		#convert ps to microseconds
		self.chan0_data = self.chan0_data * 1e-6
		self.chan1_data = self.chan1_data * 1e-6
		
		###only consider time tags at which the number of photons is within user defined count range 
		self.counts_hist = np.array([np.append(self.counts_hist[0], 0), self.counts_hist[1]]).T
		not_counts_range = self.counts_hist[(self.counts_hist[:, 0] < self.min_counts) | (self.counts_hist[:, 0] > self.max_counts)]
		times_to_exclude = not_counts_range[:, 1] * 1e-6 #convert ps to microseconds
		#for loop checks chan0 and chan1 arrays and deletes data that does not fit user range
		for i in range(len(times_to_exclude)):
			current_value = times_to_exclude[i]
			idx0 = np.where(self.chan0_data == current_value)
			for j in range(len(idx0)):
				self.chan0_data = np.delete(self.chan0_data, j, axis = 0)
			idx1 = np.where(self.chan1_data == current_value)
			for j in range (len(idx1)):
				self.chan1_data = np.delete(self.chan1_data, j, axis = 0)

		###only get data within user defined time range
		max_x = self.ui.max_time_doubleSpinBox.value() * 1e-3 #nanoseconds to microseconds
		self.chan0_data = self.chan0_data[(self.chan0_data > max_x * -1) | (self.chan0_data < max_x)]
		self.chan1_data = self.chan1_data[(self.chan1_data > max_x * -1) | (self.chan1_data < max_x)]
		
		#user input num samples
		num_bins = self.ui.sampling_points_spinBox.value()
		
		#calculate correlation values!
		corr_array0 = self.calc_corr(self.chan0_data, self.chan1_data) #a-b
		corr_array1 = self.calc_corr(self.chan1_data, self.chan0_data) #b-a

		###average correlation data of both channels
		#find shorter array
		if len(corr_array1) > len(corr_array0):
			longer = corr_array1
			shorter = corr_array0
		else:
			longer = corr_array0
			shorter = corr_array1	
		diff = len(longer) - len(shorter)
		#add placeholder 0s to match shorter array to longer array
		for i in range(diff):
			shorter = np.append(shorter, np.array([[0, 0]]), axis=0)
		#get averages
		self.time_avg = (shorter[:,0] + longer[:,0])/2
		self.corr_avg = (shorter[:,1] + longer[:,1])/2
		
		#generate correlation histogram
		self.g2_y, self.g2_x = np.histogram([self.time_avg, self.corr_avg], num_bins)
			
		try:
			#setup integration region selection
			quarter = (np.max(self.time_avg) - np.min(self.time_avg))/4
			self.integration_region.setRegion((quarter, quarter * 3))
			self.tau_g2_plot.addItem(self.integration_region)
			
			#plot g2 data
			self.tau_g2_plot.plot(self.g2_x, self.g2_y, stepMode=True, \
				brush=(0,0,255,150), fillLevel=0, fillOutline=True)
		except Exception as err:
			refresh_text_box(self.ui.textBrowser, format(err) + "\n")
				
	def on_plot(self):
		"""
		On plot button press, clear and re-plot data.
		"""
		self.ui.counts_plotWidget.clear()
		self.ui.g2_plotWidget.clear()
		self.plot_timecorr()	
		self.plot_g2()
		
	def calculate_integration(self):
		"""
		Integrate over region and display on ui.
		"""
		start_index = find_nearest(self.corr_avg, self.integration_start)
		stop_index = find_nearest(self.corr_avg, self.integration_stop)
		integration_region_x = self.corr_avg[start_index : stop_index]
		integration_region_y = self.corr_avg[start_index : stop_index]
		area = np.trapz(x = integration_region_x, y = integration_region_y)
		self.ui.integration_result_label.setText("Result: %d" % area)

"""Run the Main Window"""
def run():
	win = MainWindow()
	QtGui.QApplication.instance().exec_()
	return win
	
#uncomment if you want to run directly from this py file
run()