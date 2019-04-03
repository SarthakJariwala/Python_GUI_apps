# -*- coding: utf-8 -*-
"""
Created on Wed Mar 27 16:50:26 2019

@author: Sarthak
"""

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import sys
import seabreeze.spectrometers as sb

pg.mkQApp()
pg.setConfigOption('background', 'w')

uiFile = "OO_acquire_gui.ui"

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
#        self.ui.clear_pushButton.clicked.connect(self.clear_plot)
        
        self.spec = sb.Spectrometer(sb.list_devices()[0])
        self.ui.status_textBrowser.setText("Ocean Optics Device Connection Initiated\n\n Device:\n\n"+str(sb.list_devices()[0]))
        
        self.save_filename = None
        self.y = None
        
        self.show()
    
    def connect_spectrometer(self):
        devices = sb.list_devices()
        self.spec = sb.Spectrometer(devices[0])
        
    
    def save_file_location(self):
        filename = QtWidgets.QFileDialog.getSaveFileName(self,
                                                         caption="NO EXTENSION IN FILENAME")
        self.save_filename = filename[0]
    
    def save_single_spec(self):
        save_array = np.zeros(shape=(2048,2))
        save_array[:,0] = self.spec.wavelengths()
        save_array[:,1] = self.y
        np.savetxt(self.save_filename+".txt", save_array, fmt = '%.5f', 
                   header = 'Wavelength (nm), Intensity (counts)', delimiter = ' ')

    def live(self):
            
        intg_time_ms = self.ui.intg_time_spinBox.value()
        self.spec.integration_time_micros(intg_time_ms*1e3)
        
        scans_to_avg = self.ui.scan_to_avg_spinBox.value()
        Int_array = np.zeros(shape=(2048,scans_to_avg))
        save_array = np.zeros(shape=(2048,2))
        
        self.ui.plot.setLabel('left', 'Intensity', units='a.u.')
        self.ui.plot.setLabel('bottom', 'Wavelength', units='nm')
        j = 0
        while self.ui.connect_checkBox.isChecked(): # this while loop works!
            for i in range(scans_to_avg):
        
                data = self.spec.spectrum(correct_dark_counts=self.ui.correct_dark_counts_checkBox.isChecked())
                Int_array[:,i] = data[1]
                self.y = np.mean(Int_array, axis=-1)
                save_array[:,0] = self.spec.wavelengths()
                save_array[:,1] = self.y
            
            self.ui.plot.plot(self.spec.wavelengths(), self.y, pen='r', clear=True)
            
            if self.ui.save_every_spec_checkBox.isChecked():
                np.savetxt(self.save_filename+str(j)+".txt", save_array, fmt = '%.5f', 
                           header = 'Wavelength (nm), Intensity (counts)', delimiter = ' ')
            
            pg.QtGui.QApplication.processEvents()
            j += 1
            
            
    def close_connection(self):
         self.spec.close()
         self.ui.status_textBrowser.setText("Ocean Optics Device Disconnected")
    
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
