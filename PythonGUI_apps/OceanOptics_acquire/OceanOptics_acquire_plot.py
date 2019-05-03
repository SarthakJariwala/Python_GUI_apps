# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import pickle
import sys
import seabreeze.spectrometers as sb
from pipython import GCSDevice

pg.mkQApp()
pg.setConfigOption('background', 'w')

#uiFile = "OO_acquire_gui.ui" "OO_PZstageScan_acquire_gui"
uiFile = "OO_PZstageScan_acquire_gui.ui"

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        TemplateBaseClass.__init__(self)
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        self.ui.live_pushButton.clicked.connect(self.live)
        self.ui.close_connection_pushButton.clicked.connect(self.close_connection)
        self.ui.config_save_pushButton.clicked.connect(self.save_file_location)
        self.ui.save_single_spec_pushButton.clicked.connect(self.save_single_spec)
        
        self.ui.init_stage_pushButton.clicked.connect(self.init_piezo_stage)
        self.ui.show_curr_pos_pushButton.clicked.connect(self.current_piezo_stage_pos)
        self.ui.center_stage_pushButton.clicked.connect(self.center_piezo)
        self.ui.abs_mov_pushButton.clicked.connect(self.abs_mov)
        self.ui.rel_mov_pushButton.clicked.connect(self.rel_mov)
        self.ui.estimate_scan_time_pushButton.clicked.connect(self.estimate_scan_time)
        self.ui.path_to_folder_pushButton.clicked.connect(self.save_file_location)
        self.ui.start_x_y_sacan_pushButton.clicked.connect(self.x_y_scan)
        
#        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        
        self.spec = sb.Spectrometer(sb.list_devices()[0])
        self.ui.status_textBrowser.setText("Ocean Optics Device Connection Initiated\n\n Device:\n\n"+str(sb.list_devices()[0]))
        
        self.save_filename = None
        self.y = None
        
        self.show()
    
    def connect_spectrometer(self):
        if self.spec is not None:
            devices = sb.list_devices()
            self.spec = sb.Spectrometer(devices[0])
        else:
            self.ui.status_textBrowser.setText("Already Connected")
        
    def init_piezo_stage(self):
        self.pi_device = GCSDevice("E-710")	# Creates a Controller instant
        self.pi_device.ConnectNIgpib(board=0,device=4) # Connect to GPIB board
        self.ui.status_textBrowser.setText('Connected: {}'.format(self.pi_device.qIDN().strip()))
        print('connected: {}'.format(self.pi_device.qIDN().strip()))
        
        self.axes = self.pi_device.axes[0:2] # selecting x and y axis of the stage
        
        self.pi_device.INI()
        self.pi_device.REF(axes=self.axes)
        
        self.pi_device.SVO(axes=self.axes, values=[True,True])	# Turn on servo control for both axes
        self.ui.status_textBrowser.setText("Current Stage Position:\n"+self.pi_device.qPOS(axes=self.axes))
        print(self.pi_device.qPOS(axes=self.axes))
        
    def center_piezo(self):
        self.pi_device.MOV(axes=self.axes, values=[50,50])
        self.ui.status_textBrowser.setText("Piezo Stage Centered: [50x,50y]")
    
    def current_piezo_stage_pos(self):
        self.ui.status_textBrowser.setText("Current Stage Position:\n"+self.pi_device.qPOS(axes=self.axes))
    
    def abs_mov(self):
        x_abs_pos = self.ui.x_abs_doubleSpinBox.value()
        y_abs_pos = self.ui.y_abs_doubleSpinBox.value()
        self.pi_device.MOV(axes=self.axes, values=[x_abs_pos,y_abs_pos])
        
    def rel_mov(self):
        x_rel_pos = self.ui.x_rel_doubleSpinBox.value()
        y_rel_pos = self.ui.y_rel_doubleSpinBox.value()
        self.pi_device.MVR(axes=self.axes, values=[x_rel_pos,y_rel_pos])
        
    def estimate_scan_time(self):
        x_scan_size = self.ui.x_size_doubleSpinBox.value()
        y_scan_size = self.ui.y_size_doubleSpinBox.value()
        
        x_step = self.ui.x_step_doubleSpinBox.value()
        y_step = self.ui.y_step_doubleSpinBox.value()
        
        y_range = int(np.ceil(y_scan_size/y_step))
        x_range = int(np.ceil(x_scan_size/x_step))
        
        total_points = x_range*y_range
        
        intg_time_ms = self.ui.intg_time_spinBox.value()
        scans_to_avg = self.ui.scan_to_avg_spinBox.value()
        
        total_time = total_points*(intg_time_ms*1e-3)*(scans_to_avg) # in seconds
        
        self.ui.status_textBrowser.setText("Estimated scan time: "+str(total_time/60)+" mins")
    
    def x_y_scan(self):
        x_start = self.ui.x_start_doubleSpinBox.value()
        y_start = self.ui.y_start_doubleSpinBox.value()
        
        x_scan_size = self.ui.x_size_doubleSpinBox.value()
        y_scan_size = self.ui.y_size_doubleSpinBox.value()
        
        x_step = self.ui.x_step_doubleSpinBox.value()
        y_step = self.ui.y_step_doubleSpinBox.value()
        
        if y_scan_size == 0:
            y_scan_size = 1
        
        if x_scan_size == 0:
            x_scan_size = 1
        
        if y_step == 0:
            y_step = 1
            
        if x_step == 0:
            x_step = 1
            
        y_range = int(np.ceil(y_scan_size/y_step))
        x_range = int(np.ceil(x_scan_size/x_step))
        
        # Define empty array for saving intensities
        data_array = np.zeros(shape=(x_range*y_range,2048))
        
        self.ui.status_textBrowser.setText("Starting Scan...")
        # Move to the starting position
        self.pi_device.MOV(axes=self.axes, values=[x_start,y_start])
        
        self.ui.status_textBrowser.setText("Scan in Progress...")
        
        # This for loop works perfectly for raster x-y scan
        # Needs a little modification for line scans! --- np.ceil fixed the error!
        k = 0 
        for i in range(y_range):
            for j in range(x_range):
                print(self.pi_device.qPOS(axes=self.axes))
                
                self._read_spectrometer()
                data_array[k,:] = self.y
                self.ui.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
                
                self.pi_device.MVR(axes=self.axes[0], values=[x_step])
                
                self.ui.progressBar.setValue(((k+1)/(x_range*y_range)))
                k+=1
                
            if i == y_range-1: # this if statement is there to keep the stage at the finish position (in x) and not bring it back like we were doing during the scan 
                self.pi_device.MVR(axes=self.axes[1], values=[y_step])
            else:
                self.pi_device.MVR(axes=self.axes, values=[-x_scan_size, y_step])
        
        self.ui.status_textBrowser.setText("Scan Complete!\nSaving Data...")
        
        save_dict = {"Wavelengths": self.spec.wavelengths(), "Intensities": data_array,
                     "Scan Parameters":{"X scan start (um)": x_start, "Y scan start (um)": y_start,
                                        "X scan size (um)": x_scan_size, "Y scan size (um)": y_scan_size,
                                        "X step size (um)": x_step, "Y step size (um)": y_step}
                     }
        
        pickle.dump(save_dict, open(self.save_folder+"/"+self.ui.lineEdit.text()+"_raw_PL_spectra_data.pkl", "wb"))
        
        self.ui.status_textBrowser.setText("Data saved!")
        
        
    
    def save_file_location(self):
        self.save_folder = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Select Folder")
#        filename = QtWidgets.QFileDialog.getSaveFileName(self,
#                                                         caption="NO EXTENSION IN FILENAME")
#        self.save_filename = filename[0]
    
    def save_single_spec(self):
        save_array = np.zeros(shape=(2048,2))
        save_array[:,0] = self.spec.wavelengths()
        save_array[:,1] = self.y
        np.savetxt(self.save_folder+"/"+self.ui.lineEdit.text()+".txt", save_array, fmt = '%.5f', 
                   header = 'Wavelength (nm), Intensity (counts)', delimiter = ' ')

    def live(self):
            
#        intg_time_ms = self.ui.intg_time_spinBox.value()
#        self.spec.integration_time_micros(intg_time_ms*1e3)
#        
#        scans_to_avg = self.ui.scan_to_avg_spinBox.value()
#        Int_array = np.zeros(shape=(2048,scans_to_avg))
        save_array = np.zeros(shape=(2048,2))
        save_array[:,0] = self.spec.wavelengths()
        
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Wavelength', units='nm')
        j = 0
        while self.ui.connect_checkBox.isChecked(): # this while loop works!
#            for i in range(scans_to_avg):
#        
#                data = self.spec.spectrum(correct_dark_counts=self.ui.correct_dark_counts_checkBox.isChecked())
#                Int_array[:,i] = data[1]
#                self.y = np.mean(Int_array, axis=-1)
#                save_array[:,0] = self.spec.wavelengths()
#                save_array[:,1] = self.y
            
            self._read_spectrometer() #substitute for above 'for' loop
            save_array[:,1] = self.y
            
            self.ui.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
            
            if self.ui.save_every_spec_checkBox.isChecked():
                np.savetxt(self.save_folder+"/"+self.ui.lineEdit.text()+str(j)+".txt", save_array, fmt = '%.5f', 
                           header = 'Wavelength (nm), Intensity (counts)', delimiter = ' ')
            
            pg.QtGui.QApplication.processEvents()
            j += 1
    
    def _read_spectrometer(self):
        if self.spec is not None:
            
            intg_time_ms = self.ui.intg_time_spinBox.value()
            self.spec.integration_time_micros(intg_time_ms*1e3)
            
            scans_to_avg = self.ui.scan_to_avg_spinBox.value()
            Int_array = np.zeros(shape=(2048,scans_to_avg))
            
            for i in range(scans_to_avg): #software average
                data = self.spec.spectrum(correct_dark_counts=self.ui.correct_dark_counts_checkBox.isChecked())
                Int_array[:,i] = data[1]
                self.y = np.mean(Int_array, axis=-1)
        
        else:
            self.ui.status_textBrowser.setText("Connect to Spectrometer!")
            
            
    def close_connection(self):
         self.spec.close()
         self.ui.status_textBrowser.setText("Ocean Optics Device Disconnected")
         del self.spec
         self.spec = None
    
    def close_application(self):
        choice = QtGui.QMessageBox.question(self, 'EXIT!',
                                            "Do you want to exit the app?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass
        
        
win = MainWindow()


## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
