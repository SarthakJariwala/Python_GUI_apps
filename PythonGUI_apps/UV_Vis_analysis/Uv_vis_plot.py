# -*- coding: utf-8 -*-
"""
Created on Thu May 10 13:45:12 2018

@author: Sarthak
"""

import numpy as np
import matplotlib.pyplot as plt

"""Load data files here - """
data = np.loadtxt('E:/IonExchange/Irikas_data/6_14_2018/MAI.csv', delimiter = ',', skiprows = 1)
data2 = np.loadtxt('E:/IonExchange/Irikas_data/6_14_2018/MAIFAI.csv', delimiter = ',', skiprows = 1)
#data3 = np.loadtxt('E:/PassivationBeyondMAPI/UV-vis/2_20_2018_abs_bef_aftertreatment/Cs17-Br15-PbSCN-CON.csv', delimiter = ',', skiprows = 1)


image_name = 'Comibed' + '.tiff'

"""Enter plot legends here ---"""
legend1 = 'MAPbI$_3$ - Covered Side'
legend2 = 'FA - Exchanged Side'
#legend3 = 'Cs17/Br15'

Wavelength = data[:,0] # in nm
Absorbance = data[:,1]
Absorbance = Absorbance - np.mean(Absorbance[(Wavelength>800) & (Wavelength<1000)])

Wavelength2 = data2[:,0] # in nm
Absorbance2 = data2[:,1]
Absorbance2 = Absorbance2 - np.mean(Absorbance2[(Wavelength2>900) & (Wavelength2<1000)])

#Wavelength3 = data3[:,0] # in nm
#Absorbance3 = data3[:,1]

"""Recylce params for plotting"""
plt.rc('xtick', labelsize = 20)
plt.rc('xtick.major', pad = 3)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', lw = 1.5, markersize = 7.5)
plt.rc('legend', fontsize = 20)
plt.rc('axes', linewidth = 3.5)
""" Plotting UV-Vis Data """

plt.figure(figsize=(8,6))
plt.tick_params(direction='out', length=8, width=3.5)
#plt.plot(Wavelength3, Absorbance3, linewidth = 3)
plt.plot(Wavelength, Absorbance, linewidth = 3, color = 'b')
plt.plot(Wavelength2, Absorbance2, linewidth = 3, color = 'r')
plt.xlim(600,900)
plt.ylim(-0.1,2)
plt.xlabel('Wavelength (nm)', fontsize = 20)
plt.ylabel('Absorbance (a.u.)', fontsize = 20)
#plt.legend([legend1, legend2])

#plt.savefig('E:/IonExchange/Irikas_data/6_14_2018/Scattering_corrected_Covered_Exchanged_side_no_legend.tiff', bbox_inches='tight', dpi = 300)
