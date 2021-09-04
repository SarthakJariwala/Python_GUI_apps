from pyqtgraph.Qt import QtCore, QtGui, QtWidgets
import numpy as np

def refresh_text_box(textBrowser, string_update): 
	textBrowser.append(string_update) #append string
	QtGui.QApplication.processEvents() #update gui for pyqt

def find_nearest(array, value):
	"""
	Find index with nearest value

	:param array: array-like in which to look for value
	:type array: array-like
	
	:param value: value to look for in array
	:type value: float
	
	:returns: index of nearest value
	:rtype: int
	"""
	idx = (np.abs(array - value)).argmin()
	return idx