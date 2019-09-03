# -*- coding: utf-8 -*-
"""
Created on Mon Sep  2 17:04:49 2019

@author: Sarthak
"""

from pathlib import Path
import pyqtgraph as pg
from pyqtgraph.python2_3 import asUnicode
from pyqtgraph.Qt import QtCore, QtGui


pg.mkQApp()
pg.setConfigOption('background', 'w')


base_path = Path(__file__).parent
file_path = (base_path / "Table_widget_gui.ui").resolve()

uiFile = file_path

WindowTemplate, TemplateBaseClass = pg.Qt.loadUiType(uiFile)

class MainWindow(TemplateBaseClass):  

    def __init__(self):
        super(TemplateBaseClass, self).__init__()
        
        # Create the main window
        self.ui = WindowTemplate()
        self.ui.setupUi(self)
        
        self.clear()
        
        self.ui.clear_pushButton.clicked.connect(self.clear)
        self.ui.add_row_pushButton.clicked.connect(self.add_row)
        self.ui.add_column_pushButton.clicked.connect(self.add_column)
        
        """Saving and Copying --- implemented from pyqtgraph TableWidget"""
        self.contextMenu = QtGui.QMenu()
        self.contextMenu.addAction('Copy Selection').triggered.connect(self.copySel)
        self.contextMenu.addAction('Copy All').triggered.connect(self.copyAll)
        self.contextMenu.addAction('Save Selection').triggered.connect(self.saveSel)
        self.contextMenu.addAction('Save All').triggered.connect(self.saveAll)
        
        self.show()
        
    def clear(self):
        self.ui.tableWidget.clear()
        self.verticalHeadersSet = False
        self.horizontalHeadersSet = False
        
    def add_row(self):
        row_position = self.ui.tableWidget.rowCount()
        self.ui.tableWidget.insertRow(row_position)
        
    def add_column(self):
        column_position = self.ui.tableWidget.columnCount()
        self.ui.tableWidget.insertColumn(column_position)
        
    def save_table(self):# Needs to be implemented
        print(self.ui.tableWidget.currentItem().text())
    
    def serialize(self, useSelection=False):
        """Convert entire table (or just selected area) into tab-separated text values"""
        if useSelection:
            selection = self.ui.tableWidget.selectedRanges()[0]
            rows = list(range(selection.topRow(),
                              selection.bottomRow() + 1))
            columns = list(range(selection.leftColumn(),
                                 selection.rightColumn() + 1))        
        else:
            rows = list(range(self.ui.tableWidget.rowCount()))
            columns = list(range(self.ui.tableWidget.columnCount()))
    
        data = []
        if self.horizontalHeadersSet:
            row = []
            if self.verticalHeadersSet:
                row.append(asUnicode(''))
            
            for c in columns:
                row.append(asUnicode(self.ui.tableWidget.horizontalHeaderItem(c).text()))
            data.append(row)
        
        for r in rows:
            row = []
            if self.verticalHeadersSet:
                row.append(asUnicode(self.ui.tableWidget.verticalHeaderItem(r).text()))
            for c in columns:
                item = self.ui.tableWidget.item(r, c)
                if item is not None:
                    row.append(asUnicode(item.text()))
                else:
                    row.append(asUnicode(''))
            data.append(row)
            
        s = ''
        for row in data:
            s += ('\t'.join(row) + '\n')
        return s


    def copySel(self):
        """Copy selected data to clipboard."""
        QtGui.QApplication.clipboard().setText(self.serialize(useSelection=True))

    def copyAll(self):
        """Copy all data to clipboard."""
        QtGui.QApplication.clipboard().setText(self.serialize(useSelection=False))


    def saveSel(self):
        """Save selected data to file."""
        self.save(self.serialize(useSelection=True))


    def saveAll(self):
        """Save all data to file."""
        self.save(self.serialize(useSelection=False))


    def save(self, data):
        fileName = QtGui.QFileDialog.getSaveFileName(self, "Save As..", "", "Tab-separated values (*.tsv)")
        if fileName == '':
            return
        open(fileName[0], 'w').write(data)

    def contextMenuEvent(self, ev):
        self.contextMenu.popup(ev.globalPos())
        
    def keyPressEvent(self, ev):
        if ev.key() == QtCore.Qt.Key_C and ev.modifiers() == QtCore.Qt.ControlModifier:
            ev.accept()
            self.copySel()
#        else:
#            QtGui.QTableWidget.keyPressEvent(self, ev)
"""Run the Main Window"""    
def run():
    win = MainWindow()
    QtGui.QApplication.instance().exec_()
    return win

#run()