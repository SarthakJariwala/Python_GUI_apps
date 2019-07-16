# -*- coding: utf-8 -*-
"""
Created on Thu May 10 15:47:59 2018

@author: Sarthak
"""

import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('E:/IonExchange/Irikas_data/6_14_2018/MAIFAI.csv', delimiter = ',', skiprows = 1)

Plot_title = 'FA - Exchanged Side'

#image_name = 'Cs25Br20'

hv_min = 1.55 # in eV
hv_max = 1.62 # in eV

"""No user input needed"""

Wavelength = data[:,0] # in nm
Absorbance = data[:,1] 
Absorbance = Absorbance - np.mean(Absorbance[(Wavelength>900) & (Wavelength<1000)])# scattering correction

hv = 1240/Wavelength # in eV

Alpha_hv = (Absorbance * hv)**2.0

"""Fitting here---"""

index = (hv > hv_min) & (hv < hv_max)

model = np.polyfit(hv[index], Alpha_hv[index], 1) 

Alpha_hv_fit = hv * model[0] + model[1] #This is the linear fit

"""Calculating Bandgap (in eV)"""

Eg = - model[1]/model[0]
print('Bandgap (eV) - ')
print(Eg)

"""Recylce params for plotting"""
plt.rc('xtick', labelsize = 20)
plt.rc('xtick.major', pad = 3)
plt.rc('ytick', labelsize = 20)
plt.rc('lines', lw = 1.5, markersize = 7.5)
plt.rc('legend', fontsize = 20)
plt.rc('axes', linewidth = 3.5)

"""Plotting Tauc plot with the Linear fit"""

plt.figure(figsize = (8,6))
plt.tick_params(direction='out', length=8, width=3.5)
plt.plot(hv, Alpha_hv, linewidth = 3, color = 'r')
plt.plot(hv, Alpha_hv_fit, linewidth = 2, color = 'black')
plt.xlim(1,2)
plt.ylim(0, np.max(Alpha_hv[index]) + 1)
plt.xlabel('h$\\nu$ (eV)', fontsize = 20)
plt.ylabel('($\\alpha$h$\\nu$)$^2$', fontsize = 20)
plt.title(Plot_title, fontsize = 20)

plt.text(1.2, 1.2, r'E$_{g}$ = %.2f eV'%Eg, fontsize = 15)
plt.tight_layout()

#plt.savefig('E:/IonExchange/Irikas_data/6_14_2018/Scattering_corrected_FA_exch_side_TaucPlot.tiff', bbox_inches='tight', dpi = 300)
#plt.savefig('E:/PassivationBeyondMAPI/UV-vis/2_28_2018_CsBr_mixed/' + image_name, bbox_inches='tight', dpi = 300)
