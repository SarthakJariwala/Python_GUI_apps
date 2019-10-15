from __future__ import division, print_function, absolute_import
from ScopeFoundry import BaseApp
from ScopeFoundry.helper_funcs import load_qt_ui_file, sibling_path
from collections import OrderedDict
import os
from qtpy import QtCore, QtWidgets, QtGui
import pyqtgraph as pg
import pyqtgraph.dockarea as dockarea
import numpy as np
from ScopeFoundry.logged_quantity import LQCollection
from scipy.stats import spearmanr
import argparse
from .h5_tree import H5TreeSearchView
from .pkl_tree import PklTreeSearchView



class H5ViewPlot(BaseApp):
    
    name = "h5_view_plot"
    
    def __init__(self, argv):
        pg.setConfigOption('imageAxisOrder', 'row-major')
        BaseApp.__init__(self, argv)
        self.setup()
        parser = argparse.ArgumentParser()
        for lq in self.settings.as_list():
            parser.add_argument("--" + lq.name)
        args = parser.parse_args()
        for lq in self.settings.as_list():
            if lq.name in args:
                val = getattr(args,lq.name)
                if val is not None:
                    lq.update_value(val)

    def setup(self):
        self.ui_filename = sibling_path(__file__, "h5_view_and_plot_gui.ui")
        self.ui = load_qt_ui_file(self.ui_filename)
        self.ui.show()
        self.ui.raise_()

        self.settings.New('data_filename', dtype='file')
        
        self.settings.data_filename.add_listener(self.on_change_data_filename)

        self.settings.New('view_name', dtype=str, initial='0', choices=('0',))
        
        # UI Connections
        self.settings.data_filename.connect_to_browse_widgets(self.ui.data_filename_lineEdit, 
                                                              self.ui.data_filename_browse_pushButton)
        self.ui.plot_pushButton.clicked.connect(self.plot_dataset)
        self.ui.dataset_listWidget.currentItemChanged.connect(self.on_data_selection)
        self.ui.plot_radioButton.toggled.connect(self.update_data_widget)
        self.ui.image_radioButton.toggled.connect(self.update_data_widget)

        #set up image item for 2d array
        self.data_img_layout = pg.GraphicsLayoutWidget()
        self.ui.imageItem_page.layout().addWidget(self.data_img_layout)
        self.data_img_layout = self.data_img_layout.addViewBox()
        self.data_img  = pg.ImageItem()
        self.data_img_layout.addItem(self.data_img)

        #set up image view for 3d array
        self.ui.data_imageView.getView().invertY(False)

        self.h5treeview = H5TreeSearchView(self)
        self.ui.dataview_page.layout().addWidget(self.h5treeview.ui)
        self.h5treeview.ui.hide()
        self.ui.show()

    def on_change_data_filename(self):
        """ Handle file change """
        try:
            fname = self.settings.data_filename.val 
            if os.path.isfile(fname):
                self.f = self.h5treeview.on_change_data_filename(fname)
                self.ui.dataview_placeholder.hide()
                self.h5treeview.ui.show()
        except:
            pass

    def plot_dataset(self):
        """ Plot data set depending on dataset shape and plot type option. """
        self.plot = self.ui.data_plotWidget.getPlotItem()
        self.plot.clear()

        data = self.dataset[()]
        if self.dataset_shape == 1:
            x_start = self.ui.plotWidget_x_start_spinBox.value()
            x_end = self.ui.plotWidget_x_end_spinBox.value()
            num_points = self.dataset.shape[0]
            x_values = np.linspace(x_start, x_end, num_points) 
            self.plot.plot(x_values, data)
        elif self.dataset_shape == 2 and self.ui.plot_radioButton.isChecked():
            self.plot.plot(data[0], data[1]) # TODO check and test this
        elif self.dataset_shape == 2 and self.ui.image_radioButton.isChecked():
            self.data_img.setImage(data)
        elif self.dataset_shape == 3:
            if self.f['Cube/Info/Cube'].attrs['AcqMode'] == b'Hyperspectral Acquisition': # This works for our PhotonEtc. Hyperspectral Camera output
                x_start = int(self.f['Cube/Info/Cube'].attrs['LowerWavelength'])
                x_end = int(self.f['Cube/Info/Cube'].attrs['UpperWavelength'])
            else:
                x_start = self.ui.imageView_x_start_spinBox.value()
                x_end = self.ui.imageView_x_end_spinBox.value()
            num_points = self.dataset.shape[0]
            x_values = np.linspace(x_start, x_end, num_points) #scale x axis
            self.ui.data_imageView.setImage(data, xvals=x_values)

    def on_data_selection(self):
        """ Handle dataset selection """
        try:
            dataset_name = self.ui.dataset_listWidget.currentItem().text()
            self.dataset = self.h5treeview.dataset_dict[dataset_name]
            self.dataset_shape = len(self.dataset[()].shape)
            self.update_data_widget()
            if self.dataset_shape == 1:
                self.ui.plot_type_groupBox.setEnabled(False)
                self.ui.plot_radioButton.setChecked(True)
            elif self.dataset_shape == 2:
                self.ui.plot_type_groupBox.setEnabled(True)
            elif self.dataset_shape == 3:
                self.ui.plot_type_groupBox.setEnabled(False)
                self.ui.image_radioButton.setChecked(True)
        except:
            pass

    def update_data_widget(self):
        """ Decide which widget to display based on dataset shape and plot type option. """
        if self.dataset_shape == 1:
            self.ui.data_stackedWidget.setCurrentIndex(0)
        elif self.dataset_shape == 2 and self.ui.plot_radioButton.isChecked():
            self.ui.data_stackedWidget.setCurrentIndex(0)
        elif self.dataset_shape == 2 and self.ui.image_radioButton.isChecked():
            self.ui.data_stackedWidget.setCurrentIndex(1)
        elif self.dataset_shape == 3:
            self.ui.data_stackedWidget.setCurrentIndex(2)

# if __name__ == '__main__':
#     import sys
#     app = H5ViewPlot(sys.argv)
#     sys.exit(app.exec_())    