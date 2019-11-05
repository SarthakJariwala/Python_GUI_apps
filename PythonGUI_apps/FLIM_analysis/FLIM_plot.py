import sys
import h5py
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import matplotlib.pyplot as plt
import pickle
#import time
from lmfit.models import GaussianModel
import customplotting.mscope as cpm

sys.path.append(os.path.abspath('../Lifetime_analysis'))
sys.path.append(os.path.abspath('../Spectrum_analysis'))
sys.path.append(os.path.abspath('../H5_Pkl'))
sys.path.append(os.path.abspath('../Export_Windows'))
from Lifetime_analysis import Lifetime_plot_fit
from Spectrum_analysis import Spectra_plot_fit
from H5_Pkl import h5_pkl_view
try:
    from Export_window import ExportFigureWindow
except:
    from Export_Windows.Export_window import ExportFigureWindow
# local modules
 
pg.mkQApp()
pg.setConfigOption('background', 'w')


base_path = Path(__file__).parent
file_path = (base_path / "flim_plot_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

def updateDelay(scale, time):
    """ Hack fix for scalebar inaccuracy """
    QtCore.QTimer.singleShot(time, scale.updateBar)

class MainWindow(TemplateBaseClass):  

    hist_data_signal = QtCore.pyqtSignal()

    def __init__(self):
        pg.setConfigOption('imageAxisOrder', 'row-major')
        super(TemplateBaseClass, self).__init__()
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        #set up ui signals
        self.ui.load_scan_pushButton.clicked.connect(self.open_file)
        self.ui.plot_intensity_sums_pushButton.clicked.connect(self.plot_intensity_sums)
        self.ui.plot_raw_hist_data_pushButton.clicked.connect(self.plot_raw_scan)
        self.ui.save_intensities_image_pushButton.clicked.connect(self.export_window)
        self.ui.save_intensities_array_pushButton.clicked.connect(self.save_intensities_array)
        self.ui.compare_checkBox.stateChanged.connect(self.switch_compare)
        self.ui.intensity_sums_viewBox.roi.sigRegionChanged.connect(self.line_profile_update_plot)
        self.ui.import_pkl_pushButton.clicked.connect(self.import_pkl_to_convert)
        self.ui.pkl_to_h5_pushButton.clicked.connect(self.pkl_to_h5)
        self.ui.analyze_lifetime_pushButton.clicked.connect(self.on_analyze_lifetime)
        self.ui.analyze_psf_pushButton.clicked.connect(self.on_analyze_psf)

        self.show()

    def open_file(self):
        """ Open FLIM scan file """
        try:
            self.filename = QtWidgets.QFileDialog.getOpenFileName(self, filter="Scan files (*.pkl *.h5 *.txt)")
            if ".pkl" in self.filename[0]:
                self.flim_scan_file = pickle.load(open(self.filename[0], 'rb'))
                self.scan_file_type = "pkl"
                self.launch_h5_pkl_viewer()
                self.get_data_params()
            elif ".h5" in self.filename[0]:
                self.flim_scan_file = h5py.File(self.filename[0], 'r')
                self.scan_file_type = "h5"
                self.launch_h5_pkl_viewer()
                self.get_data_params()
            elif ".txt" in self.filename[0]:
                self.intensity_sums = np.loadtxt(self.filename[0]).T
                self.stepsize_window = StepSizeWindow()
                self.stepsize_window.stepsize_signal.connect(self.get_stepsize)
                self.scan_file_type = "txt"
            # self.pkl_file = pickle.load(open(self.filename[0], 'rb'))
        except Exception as err:
            print(format(err))
    
    def launch_h5_pkl_viewer(self):
        """ Launches H5/PKL viewer to give an insight into the data and its structure"""
        viewer_window = h5_pkl_view.H5PklView(sys.argv)
        viewer_window.settings['data_filename'] = self.filename[0]

    def import_pkl_to_convert(self):
        """ Open pkl file to convert to h5 """
        try:
            self.pkl_to_convert = QtWidgets.QFileDialog.getOpenFileName(self)
            self.ui.result_textBrowser.append("Done Loading - .pkl to convert")
        except:
            pass
    
    def get_stepsize(self):
        """ Get step size from user input -- specfically written for loading 
        txt files from legacy labview code, but can also be run on txt file 
        saved using the new FLIM acquistion code """
        self.stepsize = self.stepsize_window.ui.stepsize_doubleSpinBox.value()
        self.x_step_size = self.stepsize
        self.y_step_size = self.stepsize

    def get_data_params(self):

        data = self.flim_scan_file
        if self.scan_file_type == "pkl":
            self.x_scan_size = data['Scan Parameters']['X scan size (um)']
            self.y_scan_size = data['Scan Parameters']['Y scan size (um)']
            self.x_step_size = data['Scan Parameters']['X step size (um)']
            self.y_step_size = data['Scan Parameters']['Y step size (um)']
            self.hist_data = data['Histogram data']
            self.time_data = data['Time data']
        else: #run this if scan file is h5
            self.x_scan_size = data['Scan Parameters'].attrs['X scan size (um)']
            self.y_scan_size = data['Scan Parameters'].attrs['Y scan size (um)']
            self.x_step_size = data['Scan Parameters'].attrs['X step size (um)']
            self.y_step_size = data['Scan Parameters'].attrs['Y step size (um)']
            self.hist_data = data['Histogram data'][()] #get dataset values
            self.time_data = data['Time data'][()]

        self.numb_x_pixels = int(self.x_scan_size/self.x_step_size)
        self.numb_y_pixels = int(self.y_scan_size/self.y_step_size)


    def plot_intensity_sums(self):
        try:
            if self.scan_file_type is "pkl" or self.scan_file_type is "h5":
                pg.setConfigOption('imageAxisOrder', 'row-major')
                self.hist_data = np.reshape(self.hist_data, newshape=(self.hist_data.shape[0], self.numb_x_pixels*self.numb_y_pixels))
                self.intensity_sums = np.sum(self.hist_data, axis=0) #sum intensities for each pixel
                self.intensity_sums = np.reshape(self.intensity_sums, newshape=(self.numb_x_pixels, self.numb_y_pixels))
            else:
                pg.setConfigOption('imageAxisOrder', 'col-major')
            self.ui.intensity_sums_viewBox.view.invertY(False) # stop y axis invert
            self.ui.intensity_sums_viewBox.setImage(self.intensity_sums, scale=
                                                  (self.x_step_size,
                                                   self.y_step_size))
            if self.scan_file_type is "pkl" or self.scan_file_type is "h5":
                self.ui.intensity_sums_viewBox.roi.setSize([self.x_scan_size, self.y_step_size]) #line roi
            scale = pg.ScaleBar(size=1,suffix='um')
            scale.setParentItem(self.ui.intensity_sums_viewBox.view)
            scale.anchor((1, 1), (1, 1), offset=(-30, -30))
            self.ui.intensity_sums_viewBox.view.sigRangeChanged.connect(lambda: updateDelay(scale, 10))
        except Exception as err:
            print(format(err))

    def line_profile_update_plot(self):
        """ Handle line profile for intensity sum viewbox """
        if hasattr(self, "intensity_sums"):
            roiPlot = self.ui.intensity_sums_viewBox.getRoiPlot()
            roiPlot.clear()
            roi = self.ui.intensity_sums_viewBox.roi

            image = self.ui.intensity_sums_viewBox.getProcessedImage()

            # Extract image data from ROI
            axes = (self.ui.intensity_sums_viewBox.axes['x'], self.ui.intensity_sums_viewBox.axes['y'])
            data, coords = roi.getArrayRegion(image.view(np.ndarray), self.ui.intensity_sums_viewBox.imageItem, axes, returnMappedCoords=True)

            #calculate sums along columns in region
            sums_to_plot = np.mean(data, axis=0)

            #get scan x-coordinates in region
            x_values = coords[1][0]

            try:
                roiPlot.plot(x_values, sums_to_plot)
            except:
                pass

    def on_analyze_psf(self):
        self.spectrum_window = Spectra_plot_fit.MainWindow()
        self.spectrum_window.show()
        self.spectrum_window.opened_from_flim = True
        sum_data = self.ui.intensity_sums_viewBox.getRoiPlot().getPlotItem().curves[0].getData()
        self.spectrum_window.sum_data_from_flim = np.asarray(sum_data)
        self.spectrum_window.ui.plot_without_bck_radioButton.setChecked(True)
        self.spectrum_window.ui.result_textBrowser.setText("Data successfully loaded from FLIM analysis.")

    def plot_raw_scan(self):
        try:
            self.hist_image = np.reshape(self.hist_data, newshape=(self.hist_data.shape[0],self.numb_x_pixels,self.numb_y_pixels))
            self.times = self.time_data[:, 0, 0]*1e-3
            self.ui.raw_hist_data_viewBox.view.invertY(False) # stops y-axis invert
            self.ui.raw_hist_data_viewBox.setImage(self.hist_image, scale=
                                                (self.x_step_size,
                                                self.y_step_size), xvals=self.times)
            self.ui.raw_hist_data_viewBox.roi.setSize([self.x_scan_size, self.y_scan_size])
            # if self.ui.compare_checkBox.isChecked():
            #     self.ui.imv2.setImage(self.hist_image, scale= (data['Scan Parameters']['X step size (um)'],
            #                                     data['Scan Parameters']['Y step size (um)']), xvals=self.times)
            self.switch_compare()
            self.ui.raw_hist_data_viewBox.ui.roiBtn.clicked.connect(self.switch_compare)
            scale = pg.ScaleBar(size=1,suffix='um')
            scale.setParentItem(self.ui.raw_hist_data_viewBox.view)
            scale.anchor((1, 1), (1, 1), offset=(-30, -30))
            self.ui.raw_hist_data_viewBox.view.sigRangeChanged.connect(lambda: updateDelay(scale, 10))
            
        except Exception as err:
            print(format(err))

    def switch_compare(self):
        """
        Handles compare checkbox. If checked, show second ROI on raw histogram data that user can use for comparison to first ROI.
        """
        if self.ui.compare_checkBox.isChecked() and hasattr(self, "hist_image"):
            if not hasattr(self, "roi2"): #create roi if doesn't exist yet
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
        else: #if not checked, hide roi
            if hasattr(self, "roi2"):
                self.roi2.hide()
                self.roi2_plot.hide()

    def update_roi2_plot(self):
        """ Update plot corresponding to second roi """
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
        if hasattr(self, "roi2_plot"): #make sure second plot is properly cleared everytime
            self.roi2_plot.clear()
            c = self.ui.raw_hist_data_viewBox.getRoiPlot().getPlotItem().curves.pop()
            c.scene().removeItem(c)
        self.roi2_plot = self.ui.raw_hist_data_viewBox.getRoiPlot().plot(xvals, data, pen='r')

    def get_raw_hist_curve(self, curve_index):
        #curve_index = 0 for original roi
        #curve_index = 1 for second comparison roi
        curves = self.ui.raw_hist_data_viewBox.getRoiPlot().getPlotItem().curves
        return curves[curve_index].getData()

    def on_analyze_lifetime(self):
        self.lifetime_window = Lifetime_plot_fit.MainWindow()
        self.lifetime_window.show()
        self.lifetime_window.opened_from_flim = True
        self.lifetime_window.hist_data_from_flim = np.asarray(self.get_raw_hist_curve(0))
        self.lifetime_window.ui.Result_textBrowser.setText("Data successfully loaded from FLIM analysis.")
    
    def export_window(self):
        self.export_window = ExportFigureWindow()
        self.export_window.ui.vmin_spinBox.setValue(np.min(self.intensity_sums))
        self.export_window.ui.vmax_spinBox.setValue(np.max(self.intensity_sums))
        self.export_window.export_fig_signal.connect(self.save_intensities_image)

    def save_intensities_image(self):
        try:
            folder = os.path.dirname(self.filename[0])
            filename_ext = os.path.basename(self.filename[0])
            filename = os.path.splitext(filename_ext)[0] #get filename without extension
            save_to = folder + "\\" + filename + "_intensity_sums.png"
            if self.export_window.ui.reverse_checkBox.isChecked():
                colormap = str(self.export_window.ui.cmap_comboBox.currentText())+"_r"
            else:
                colormap = str(self.export_window.ui.cmap_comboBox.currentText())
            if self.export_window.ui.cbar_checkBox.isChecked():
                label = str(self.export_window.ui.cbar_label.text())
            else:
                label = "PL Intensity (a.u.)"
            cpm.plot_confocal(self.intensity_sums, FLIM_adjust=False, 
                              stepsize=np.abs(self.x_step_size),cmap=colormap, 
                              cbar_label=label, vmin=self.export_window.ui.vmin_spinBox.value(), 
                              vmax=self.export_window.ui.vmax_spinBox.value())
            plt.savefig(save_to, bbox_inches='tight', dpi=300)
        except Exception as e:
            print(format(e))

    def save_intensities_array(self):
        try:
            folder = os.path.dirname(self.filename[0])
            filename_ext = os.path.basename(self.filename[0])
            filename = os.path.splitext(filename_ext)[0] #get filename without extension
            save_to = folder + "\\" + filename + "_intensity_sums.txt"
            np.savetxt(save_to, self.intensity_sums.T, fmt='%f') #save transposed intensity sums, as original array handles x in cols and y in rows
        except:
            pass

    def pkl_to_h5(self):
        #Convert scan .pkl file into h5
        try:
            folder = os.path.dirname(self.pkl_to_convert[0])
            filename_ext = os.path.basename(self.pkl_to_convert[0])
            filename = os.path.splitext(filename_ext)[0] #get filename without extension
            pkl_file = pickle.load(open(self.pkl_to_convert[0], 'rb'))

            h5_filename = folder + "/" + filename + ".h5"
            h5_file = h5py.File(h5_filename, "w")
            self.traverse_dict_into_h5(pkl_file, h5_file)
        except Exception as err:
            print(format(err))

    def traverse_dict_into_h5(self, dictionary, h5_output):
        #Create an h5 file using .pkl with scan data and params
        for key in dictionary:
            if type(dictionary[key]) == dict: #if subdictionary, create a group
                group = h5_output.create_group(key)
                previous_dict = dictionary[key]
                self.traverse_dict_into_h5(dictionary[key], group) #traverse subdictionary
            else:
                if key == "Histogram data" or key == "Time data":
                    h5_output.create_dataset(key, data=dictionary[key])
                else:
                    h5_output.attrs[key] = dictionary[key] #if not dataset, create attribute

    def close_application(self):
        choice = QtGui.QMessageBox.question(self, 'EXIT!',
                                            "Do you want to exit the app?",
                                            QtGui.QMessageBox.Yes | QtGui.QMessageBox.No)
        if choice == QtGui.QMessageBox.Yes:
            sys.exit()
        else:
            pass

"""Skip rows GUI"""
ui_file_path = (base_path / "step_size_labview_files.ui").resolve()
stepsize_WindowTemplate, stepsize_TemplateBaseClass = pg.Qt.loadUiType(ui_file_path)

class StepSizeWindow(stepsize_TemplateBaseClass):
    
    stepsize_signal = QtCore.pyqtSignal() #signal to help with pass info back to MainWindow
    
    def __init__(self):
        stepsize_TemplateBaseClass.__init__(self)

        # Create the param window
        self.ui = stepsize_WindowTemplate()
        self.ui.setupUi(self)
        self.ui.done_pushButton.clicked.connect(self.done)
        self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, False)
        self.show()
    
    def done(self):
        self.stepsize_signal.emit()
        self.close()

"""Run the Main Window"""    
def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win

#Uncomment below if you want to run this as standalone
#run()