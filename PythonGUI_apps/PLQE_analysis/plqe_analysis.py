# system imports
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph import exporters
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import matplotlib.pyplot as plt
import numpy as np

# local modules

pg.mkQApp()
pg.setConfigOption('background', 'w')

base_path = Path(__file__).parent
file_path = (base_path / "plqe_analysis_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

"""params for plotting"""
plt.rc('xtick', labelsize = 20)
plt.rc('xtick.major', pad = 3)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', lw = 1.5, markersize = 7.5)
plt.rc('legend', fontsize = 20)
plt.rc('axes', linewidth = 3.5)

class MainWindow(TemplateBaseClass):  
    
    def __init__(self):
        super(TemplateBaseClass, self).__init__()
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        #setup uv vis plot
        self.plot = self.ui.plotWidget.getPlotItem()
        self.plot.setTitle(title="Wavelength vs. Intensity")
        self.plot.setLabel('bottom', 'Wavelength', unit='nm')
        self.plot.setLabel('left', 'Intensity', unit='a.u.')
        self.plot.setLogMode(x=None, y=1)

        #setup line rois for laser and emission
        self.laser_region = pg.LinearRegionItem(brush=QtGui.QBrush(QtGui.QColor(255, 0, 0, 50)))
        self.laser_region.sigRegionChanged.connect(self.update_laser_spinBoxes)
        self.emission_region = pg.LinearRegionItem()
        self.emission_region.sigRegionChanged.connect(self.update_emission_spinBoxes)
        self.laser_region.setRegion((200, 400))
        self.emission_region.setRegion((700, 800))

        #setup ui signals
        self.ui.load_data_pushButton.clicked.connect(self.open_data_file)
        self.ui.plot_pushButton.clicked.connect(self.plot_intensity)
        self.ui.clear_pushButton.clicked.connect(self.clear)
        self.ui.calculate_plqe_pushButton.clicked.connect(self.calculate_plqe)
        self.ui.laser_start_spinBox.valueChanged.connect(self.update_laser_region)
        self.ui.laser_stop_spinBox.valueChanged.connect(self.update_laser_region)
        self.ui.emission_start_spinBox.valueChanged.connect(self.update_emission_region)
        self.ui.emission_stop_spinBox.valueChanged.connect(self.update_emission_region)

        self.show()

    def open_data_file(self):
        """ Open data file """
        try:
            self.filename = QtWidgets.QFileDialog.getOpenFileName(self)
            #self.data = np.loadtxt(self.filename[0], delimiter = '\t', skiprows = 1)
            if ".txt" in self.filename[0]:
                self.data = np.loadtxt(self.filename[0], delimiter = '\t', skiprows = 1)
            elif ".csv" in self.filename[0]:
                self.data = np.loadtxt(self.filename[0], delimiter = ',', skiprows = 1)
            elif ".qua" in self.filename[0]:#TODO: Include a Pop-up window for input for skipping header
                self.data = np.genfromtxt(self.filename[0], delimiter = '\t', skip_header = 28)
            self.cs_window = ColSelectionWindow(self.data)
            self.cs_window.col_selection_signal.connect(self.open_with_col_selection)
            self.nm = np.copy(self.data[:,0])
            #self.ref_data = np.copy(self.data[:,1])
            #self.inpath_data = np.copy(self.data[:,2])
            #self.outpath_data = np.copy(self.data[:,3])
        except Exception as err:
            print(format(err))

    def open_with_col_selection(self):
        ref_data_col = self.cs_window.ui.ref_spinBox.value() - 1 #subtract since spinboxes refer to column num and not index
        inpath_data_col = self.cs_window.ui.inpath_spinBox.value() - 1
        outpath_data_col = self.cs_window.ui.outpath_spinBox.value() - 1
        self.ref_data = np.copy(self.data[:,ref_data_col])
        self.inpath_data = np.copy(self.data[:,inpath_data_col])
        self.outpath_data = np.copy(self.data[:,outpath_data_col])

    def update_laser_spinBoxes(self):
        """ Update laser spinboxes based on line rois """
        self.laser_start, self.laser_stop = self.laser_region.getRegion()
        self.ui.laser_start_spinBox.setValue(self.laser_start)
        self.ui.laser_stop_spinBox.setValue(self.laser_stop)


    def update_emission_spinBoxes(self):
        """ Update emission spinboxes based on line rois """
        self.emission_start, self.emission_stop = self.emission_region.getRegion()
        self.ui.emission_start_spinBox.setValue(self.emission_start)
        self.ui.emission_stop_spinBox.setValue(self.emission_stop)

    def update_laser_region(self):
        """ Update laser line rois based on spinboxes """
        laser_start = self.ui.laser_start_spinBox.value()
        laser_stop = self.ui.laser_stop_spinBox.value()
        self.laser_region.setRegion((laser_start, laser_stop))

    def update_emission_region(self):
        """ Update emission line rois based on spinboxes """
        emission_start = self.ui.emission_start_spinBox.value()
        emission_stop = self.ui.emission_stop_spinBox.value()
        self.emission_region.setRegion((emission_start, emission_stop))

    def plot_intensity(self):
        try:
            self.plot.plot(self.nm, self.inpath_data, pen='r')
            self.plot.addItem(self.laser_region, ignoreBounds=True)
            self.plot.addItem(self.emission_region, ignoreBounds=True)
        except Exception as err:
            print(format(err))

    def find_nearest(self,array,value):
        idx = (np.abs(array-value)).argmin()
        return idx

    def calculate_plqe(self):
        
        nm_interp_step = 1
        nm_interp_start = np.ceil(self.nm[0] / nm_interp_step) * nm_interp_step
        nm_interp_stop = np.floor(self.nm[len(self.nm) - 1] / nm_interp_step) * nm_interp_step
        nm_interp = np.arange(nm_interp_start, nm_interp_stop + nm_interp_step, nm_interp_step)    
        
        ref_interp = np.interp(nm_interp, self.nm, self.ref_data)
        
        
        inpath_interp = np.interp(nm_interp, self.nm, self.inpath_data)
        outpath_interp = np.interp(nm_interp,  self.nm, self.outpath_data)
        
        
        """L_x is area under laser profile for experiment x"""
        """P_x_ is area under emission profile for experiment x"""
        
        
        #plt.semilogy(nm, a1_outpath_data[:,1])
        
        emission_start_idx = self.find_nearest(nm_interp, self.emission_start)
        emission_stop_idx = self.find_nearest(nm_interp, self.emission_stop)
        
        laser_start_idx = self.find_nearest(nm_interp, self.laser_start)
        laser_stop_idx = self.find_nearest(nm_interp, self.laser_stop)
        
        la = np.trapz(ref_interp[laser_start_idx: laser_stop_idx], x = nm_interp[laser_start_idx:laser_stop_idx])
        lb = np.trapz(outpath_interp[laser_start_idx: laser_stop_idx], x = nm_interp[laser_start_idx:laser_stop_idx])
        lc = np.trapz(inpath_interp[laser_start_idx: laser_stop_idx], x = nm_interp[laser_start_idx:laser_stop_idx])
        
        pa = np.trapz(ref_interp[emission_start_idx:emission_stop_idx], x = nm_interp[emission_start_idx:emission_stop_idx])
        pb = np.trapz(outpath_interp[emission_start_idx:emission_stop_idx], x = nm_interp[emission_start_idx:emission_stop_idx])
        pc = np.trapz(inpath_interp[emission_start_idx:emission_stop_idx], x = nm_interp[emission_start_idx:emission_stop_idx])
        
        absorb = 1.0 - (lc / lb)
        
        plqe = 100 * (pc - ((1.0 - absorb) * pb)) / (la * absorb)
        #print('PLQE Percent = %.3f' %(plqe))
        #return plqe
        self.ui.plqe_label.setText("%.3f" %(plqe))

    def clear(self):
        self.plot.clear()

"""Table view GUI"""
ui_file_path = (base_path / "column_selection_gui.ui").resolve()
col_selection_WindowTemplate, col_selection_TemplateBaseClass = pg.Qt.loadUiType(ui_file_path)

class ColSelectionWindow(col_selection_TemplateBaseClass):
    
    col_selection_signal = QtCore.pyqtSignal() #signal to help with pass info back to MainWindow
    
    def __init__(self, data):
        col_selection_TemplateBaseClass.__init__(self)
        # Create the param window
        self.ui = col_selection_WindowTemplate()
        self.ui.setupUi(self)
        self.ui.done_pushButton.clicked.connect(self.done)
        
        self.table_widget = pg.TableWidget()
        self.ui.data_preview_groupBox.layout().addWidget(self.table_widget)

        self.table_widget.setData(data)

        #self.setWindowFlag(QtCore.Qt.WindowCloseButtonHint, False)
        self.show()
    
    def done(self):
        self.col_selection_signal.emit()
        self.ui.textBrowser.setText("Data successfully loaded.")
        #self.close()

    def closeEvent(self, event):
        self.col_selection_signal.emit()

"""Run the Main Window"""
def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win

#run()
