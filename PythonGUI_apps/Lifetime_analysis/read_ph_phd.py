# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:52:05 2019

@author: Sarthak
"""

import picoharp_phd
import numpy
import matplotlib.pyplot as plt
#import sys

#datafile = "E:/Python code/APTES_APTMS.phd"


def read_picoharp_phd(datafile):
    parser = picoharp_phd.PicoharpParser(datafile)
    return parser

def phd_to_csv(datafile):
    parser = read_picoharp_phd(datafile)
    name, ext = datafile.rsplit('.', 1)
    
    total_curves = parser.no_of_curves()
    y = []
    for i in range(total_curves):
        res, curve = parser.get_curve(i)
        time_window = int(numpy.floor(parser.get_time_window_in_ns(curve_no)))
        curve = curve[0:time_window]
        y.append(curve)

parser = picoharp_phd.PicoharpParser(datafile)
name, ext = datafile.rsplit('.', 1)

#for i in range(parser.no_of_curves()):
#    
#    res, curve1 = parser.get_curve(0)
#    res, curve2 = parser.get_curve(1)
#    res, curve3 = parser.get_curve(2)
#    res, curve4 = parser.get_curve(3)
#    res, curve5 = parser.get_curve(4)
##    res, curve2 = parser.get_curve(1)
#    size = len(curve1)

curve_no = 3

res, curve1 = parser.get_curve(curve_no)
time_window = int(numpy.floor(parser.get_time_window_in_ns(curve_no)))
curve1 = curve1[0:time_window]
size = len(curve1)
X = numpy.arange(0, size*res, res, numpy.float)

plt.figure()
plt.plot(X,curve1)
plt.yscale('log')
plt.ylim([1,1e4])
#csvname = '%s.csv' % name
#csv = open(csvname, 'w')
#
#for x, y1, y2, y3, y4, y5 in zip(X, curve1, curve2, curve3, curve4, curve5):
#    csv.write('%f,%d,%d,%d,%d,%d\n' % (x, y1, y2, y3,y4,y5))
#
#csv.close()
#
#print('Saved %s.' % csvname)


#if __name__ == '__main__':
#    main()
