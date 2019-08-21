import sys
from pathlib import Path
import os.path
import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui, QtWidgets#, QColorDialog
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image

scale = 7.4/50
image = Image.open(r"C:\Users\lindat18\Desktop\BREAD.png")
image = image.rotate(-90, expand=True)
image = image.resize((round(image.size[0]*scale), round(image.size[1]*scale)))
image_array = np.array(image)
width = image_array.shape[0]
height = image_array.shape[1]
print(image_array.shape)
arr = image_array[:,0,0]
print(arr.shape)
