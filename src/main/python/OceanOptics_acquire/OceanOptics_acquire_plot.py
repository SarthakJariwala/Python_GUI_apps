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
import time

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
        
        self.ui.connect_checkBox.stateChanged.connect(self._handle_spectrometer_connection)
        self.ui.live_pushButton.clicked.connect(self.live)
        
        self.ui.path_to_folder_pushButton.clicked.connect(self.save_file_location)
        self.ui.save_single_spec_pushButton.clicked.connect(self.save_single_spec)
        
        self.ui.init_stage_pushButton.clicked.connect(self.init_piezo_stage)
        self.ui.show_curr_pos_pushButton.clicked.connect(self.current_piezo_stage_pos)
        self.ui.center_stage_pushButton.clicked.connect(self.center_piezo)
        self.ui.abs_mov_pushButton.clicked.connect(self.abs_mov)
        self.ui.rel_mov_pushButton.clicked.connect(self.rel_mov)
        self.ui.estimate_scan_time_pushButton.clicked.connect(self.estimate_scan_time)
        self.ui.start_x_y_sacan_pushButton.clicked.connect(self.x_y_scan)
        
#        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        self.pi_device = None
        self.spec = None
        
        self.ui.status_textBrowser.setText("Welcome!\nGUI Initiated!")
        
        self.show()
    
    def _handle_spectrometer_connection(self):
        if self.ui.connect_checkBox.isChecked():
            self.connect_spectrometer()
        else:
            self.close_connection()
        
    def connect_spectrometer(self):
        if self.spec is None:
            devices = sb.list_devices()
            self.spec = sb.Spectrometer(devices[0])
            self.ui.status_textBrowser.append("Ocean Optics Device Connected!\n\n Device:\n\n"+str(sb.list_devices()[0]))
        else:
            self.ui.status_textBrowser.append("Already Connected")
        
    def init_piezo_stage(self):
        if self.pi_device is None:
            self.pi_device = GCSDevice("E-710")	# Creates a Controller instant
            self.pi_device.ConnectNIgpib(board=0,device=4) # Connect to GPIB board
            self.ui.status_textBrowser.append('Connected: {}'.format(self.pi_device.qIDN().strip()))
    #        print('connected: {}'.format(self.pi_device.qIDN().strip()))
            
            self.axes = self.pi_device.axes[0:2] # selecting x and y axis of the stage
            
            self.pi_device.INI()
            self.pi_device.REF(axes=self.axes)
            
            self.pi_device.SVO(axes=self.axes, values=[True,True])	# Turn on servo control for both axes
            self.ui.status_textBrowser.append("Current Stage Position:\n{}".format(self.pi_device.qPOS(axes=self.axes)))
    #        print(self.pi_device.qPOS(axes=self.axes))
        else:
            self.ui.status_textBrowser.append("Piezo Stage Is Already Initialized!")
        
    def center_piezo(self):
        if self.pi_device is None:
            self.init_piezo_stage()
        self.pi_device.MOV(axes=self.axes, values=[50,50])
        self.ui.status_textBrowser.append("Piezo Stage Centered: [50x,50y]")
    
    def current_piezo_stage_pos(self):
        if self.pi_device is None:
            self.init_piezo_stage()
        self.ui.status_textBrowser.append("Current Stage Position:\n{}".format(self.pi_device.qPOS(axes=self.axes)))
    
    def abs_mov(self):
        if self.pi_device is None:
            self.init_piezo_stage()
        x_abs_pos = self.ui.x_abs_doubleSpinBox.value()
        y_abs_pos = self.ui.y_abs_doubleSpinBox.value()
        self.pi_device.MOV(axes=self.axes, values=[x_abs_pos,y_abs_pos])
        
    def rel_mov(self):
        if self.pi_device is None:
            self.init_piezo_stage()
        x_rel_pos = self.ui.x_rel_doubleSpinBox.value()
        y_rel_pos = self.ui.y_rel_doubleSpinBox.value()
        self.pi_device.MVR(axes=self.axes, values=[x_rel_pos,y_rel_pos])
        
    def estimate_scan_time(self):
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
        
        total_points = x_range*y_range
        
        intg_time_ms = self.ui.intg_time_spinBox.value()
        scans_to_avg = self.ui.scan_to_avg_spinBox.value()
        
        total_time = total_points*(intg_time_ms*1e-3)*(scans_to_avg) # in seconds
        
        self.ui.status_textBrowser.append("Estimated scan time: "+str(np.float16(total_time/60))+" mins")
    
    def x_y_scan(self):
        if self.pi_device is None:
            self.init_piezo_stage()
            
        if self.spec is None:
            self.ui.status_textBrowser.append("Spectrometer not connected!\nForce connecting the spectrometer...")
            self.connect_spectrometer()
            
        start_time = time.time()
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
        
        self.ui.status_textBrowser.append("Starting Scan...")
        # Move to the starting position
        self.pi_device.MOV(axes=self.axes, values=[x_start,y_start])
        
        self.ui.status_textBrowser.append("Scan in Progress...")
        
        k = 0 
        for i in range(y_range):
            for j in range(x_range):
#                print(self.pi_device.qPOS(axes=self.axes))
                
                self._read_spectrometer()
                data_array[k,:] = self.y
                self.ui.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
                pg.QtGui.QApplication.processEvents()
                
                self.pi_device.MVR(axes=self.axes[0], values=[x_step])
                
                self.ui.progressBar.setValue(100*((k+1)/(x_range*y_range)))
                k+=1
            # TODO
            # if statement needs to be modified to keep the stage at the finish y-pos for line scans in x, and same for y
            if i == y_range-1: # this if statement is there to keep the stage at the finish position (in x) and not bring it back like we were doing during the scan 
                self.pi_device.MVR(axes=self.axes[1], values=[y_step])
            else:
                self.pi_device.MVR(axes=self.axes, values=[-x_scan_size, y_step])
        
        self.ui.status_textBrowser.append("Scan Complete!\nSaving Data...")

        save_dict = {"Wavelengths": self.spec.wavelengths(), "Intensities": data_array,
                     "Scan Parameters":{"X scan start (um)": x_start, "Y scan start (um)": y_start,
                                        "X scan size (um)": x_scan_size, "Y scan size (um)": y_scan_size,
                                        "X step size (um)": x_step, "Y step size (um)": y_step},
                                        "OceanOptics Parameters":{"Integration Time (ms)": self.ui.intg_time_spinBox.value(),
                                                                  "Scans Averages": self.ui.scan_to_avg_spinBox.value(),
                                                                  "Correct Dark Counts": self.ui.correct_dark_counts_checkBox.isChecked()}
                     }
        
        pickle.dump(save_dict, open(self.save_folder+"/"+self.ui.lineEdit.text()+"_raw_PL_spectra_data.pkl", "wb"))
        
        self.ui.status_textBrowser.append("Data saved!\nTotal time taken:"+str(np.float16((time.time()-start_time)/60))+" mins")
          
    def save_file_location(self):
        self.save_folder = QtWidgets.QFileDialog.getExistingDirectory(self, caption="Select Folder")
    
    def save_single_spec(self):
        save_array = np.zeros(shape=(2048,2))
        self._read_spectrometer()
        save_array[:,1] = self.y
        save_array[:,0] = self.spec.wavelengths()

        np.savetxt(self.save_folder+"/"+self.ui.lineEdit.text()+".txt", save_array, fmt = '%.5f', 
                   header = 'Wavelength (nm), Intensity (counts)', delimiter = ' ')

    def live(self):
        save_array = np.zeros(shape=(2048,2))
        
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Wavelength', units='nm')
        j = 0
        while self.spec is not None:#self.ui.connect_checkBox.isChecked(): # this while loop works!
            self._read_spectrometer()
            save_array[:,1] = self.y
            
            self.ui.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
            
            if self.ui.save_every_spec_checkBox.isChecked():
                save_array[:,0] = self.spec.wavelengths()
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
            self.ui.status_textBrowser.append("Connect to Spectrometer!")
            raise Exception("Must connect to spectrometer first!")
            
            
    def close_connection(self):
        if self.spec is not None:
             self.spec.close()
             self.ui.status_textBrowser.append("Ocean Optics Device Disconnected")
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
