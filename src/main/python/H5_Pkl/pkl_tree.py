from ScopeFoundry.data_browser import DataBrowserView
from qtpy import QtWidgets
import h5py
import pickle
import numpy as np
import lmfit
            
class PklTreeSearchView(DataBrowserView):
    
    name = 'pkl_tree_search'
    
    def is_file_supported(self, fname):
        return ('.pkl' in fname)
    
    def setup(self):
        #self.settings.New('search_text', dtype=str, initial="")
        
        self.ui = QtWidgets.QWidget()
        self.ui.setLayout(QtWidgets.QVBoxLayout())
        self.search_lineEdit = QtWidgets.QLineEdit()
        self.search_lineEdit.setPlaceholderText("Search")
        self.tree_textEdit = QtWidgets.QTextEdit("")
        
        self.ui.layout().addWidget(self.search_lineEdit)
        self.ui.layout().addWidget(self.tree_textEdit)
         
        #self.settings.search_text.connect_to_widget(self.search_lineEdit)
        #self.settings.search_text.add_listener(self.on_new_search_text)
        self.search_text = ""
        
        self.search_lineEdit.textChanged.connect(self.on_new_search_text)
        
    def on_change_data_filename(self, fname=None):
        """ Handle file change """
        self.tree_textEdit.setText("loading {}".format(fname))
        try:
            self.fname = fname        
            #self.f = h5py.File(fname, 'r')
            self.dictionary = pickle.load(open(self.fname, 'rb'))
            self.on_new_search_text()
            self.databrowser.ui.statusbar.showMessage("")
            
        except Exception as err:
            msg = "Failed to load %s:\n%s" %(fname, err)
            self.databrowser.ui.statusbar.showMessage(msg)
            self.tree_textEdit.setText(msg)
            raise(err)

    def on_new_search_text(self, x=None):
        if x is not None:
            self.search_text = x.lower()
        old_scroll_pos = self.tree_textEdit.verticalScrollBar().value()
        self.tree_str = ""  
        #self.f.visititems(self._visitfunc)
        self.traverse_dict(self.dictionary, self.dictionary, 0)

        
        self.tree_text_html = \
        """<html><b>{}</b><hr/>
        <div style="font-family: Courier;">
         {} 
         </div>
         </html>""".format(self.fname, self.tree_str)
        
        self.tree_textEdit.setText(self.tree_text_html)
        self.tree_textEdit.verticalScrollBar().setValue(old_scroll_pos)

    def traverse_dict(self, dictionary, previous_dict, level):
        """
        Visit all values in the dictionary and its subdictionaries.

        dictionary -- dictionary to traverse
        previous_dict -- dictionary one level up
        level -- track how far to indent 
        """
        for key in dictionary:
            if key not in previous_dict:
                level -=1
            indent = "&nbsp;"*4*(level)

            if type(dictionary[key]) == dict:
                print_string = key
                if self.search_text and self.search_text in print_string:
                    self.tree_str += indent + """<span style="color: red;">{}</span>""".format(print_string)
                else:
                    self.tree_str += indent + "|> <b>{}/</b><br/>".format(print_string)
                level += 1
                previous_dict = dictionary[key]
                self.traverse_dict(dictionary[key], previous_dict, level)                
            else:
                value = dictionary[key]
                if type(value) == np.ndarray or type(value)==np.memmap:
                    value = str(value.shape) + " " + str(value.dtype)
                elif type(value) == lmfit.model.ModelResult:
                    value = "lmfit.model.ModelResult"
                # if type(value) == list and len(value) > 5: ##account for data stored in lists
                #     value = str(np.asarray(value).shape) + " " + str(type(value[0]))

                print_string = key + " = " + str(value)
                if self.search_text and self.search_text in print_string:
                    self.tree_str += indent + """<span style="color: red;">{}</span>""".format(print_string)
                else:
                    self.tree_str += indent + "|- {}<br/>".format(print_string)