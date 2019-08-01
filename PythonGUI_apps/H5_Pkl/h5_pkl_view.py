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

pg.setConfigOption('imageAxisOrder', 'row-major')

class H5PklView(BaseApp):
    
    name = "h5_pkl_view"
    
    def __init__(self, argv):
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
        self.ui_filename = sibling_path(__file__, "h5_pkl_view_gui.ui")
        self.ui = load_qt_ui_file(self.ui_filename)
        self.ui.show()
        self.ui.raise_()
        
        self.views = OrderedDict()        

        self.settings.New('data_filename', dtype='file')
        self.settings.New('auto_select_view',dtype=bool, initial=True)
        self.settings.New('view_name', dtype=str, initial='0', choices=('0',))

        self.settings.data_filename.add_listener(self.on_change_data_filename)

        # UI Connections/
        self.settings.data_filename.connect_to_browse_widgets(self.ui.data_filename_lineEdit, 
                                                              self.ui.data_filename_browse_pushButton)

        # set views
        self.h5treeview = H5TreeSearchView(self)
        self.load_view(self.h5treeview)
        self.pkltreeview = PklTreeSearchView(self)
        self.load_view(self.pkltreeview)

        self.settings.view_name.add_listener(self.on_change_view_name)

        self.current_view = None

        self.ui.show()

    def load_view(self, new_view): 
        # add to views dict
        self.views[new_view.name] = new_view
        
        self.ui.dataview_page.layout().addWidget(new_view.ui)
        new_view.ui.hide()
        
        # update choices for view_name
        self.settings.view_name.change_choice_list(list(self.views.keys()))
        return new_view

    def on_change_data_filename(self):
        #Handle file change
        try:
            fname = self.settings.data_filename.val 
            if not self.settings['auto_select_view']:
                self.current_view.on_change_data_filename(fname)
            else:
                view_name = self.auto_select_view(fname)
                if self.current_view is None or view_name != self.current_view.name:
                    # update view (automatically calls on_change_data_filename)
                    self.settings['view_name'] = view_name
                else:
                    # force update
                    if os.path.isfile(fname): 
                        self.current_view.on_change_data_filename(fname)
        except:
            pass
             
    def on_change_view_name(self):
        #Handle view change - happens when filetype changes
        self.ui.dataview_placeholder.hide()
        previous_view = self.current_view
        
        self.current_view = self.views[self.settings['view_name']]
        # hide current view 
        # (handle the initial case where previous_view is None )
        if previous_view:
            previous_view.ui.hide() 

        # show new view
        self.current_view.ui.show()
        
        # set datafile for new (current) view
        fname = self.settings['data_filename']
        if  os.path.isfile(fname):
            self.current_view.on_change_data_filename(self.settings['data_filename'])

    def auto_select_view(self, fname):
        #return the name of the last supported view for the given fname
        for view_name, view in list(self.views.items())[::-1]:
            if view.is_file_supported(fname):
                return view_name
        # return default file_info view if no others work
        return "h5_tree_search"

# if __name__ == '__main__':
#     import sys
#     app = H5PklView(sys.argv)
#     sys.exit(app.exec_())
