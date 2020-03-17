# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:52:05 2019

@author: Sarthak
"""

try:
    from Lifetime_analysis import picoharp_phd
except:
    import picoharp_phd
import numpy as np
import matplotlib.pyplot as plt
#import sys

#datafile = "E:/Python code/APTES_APTMS.phd"


def read_picoharp_phd(datafile):
    parser = picoharp_phd.PicoharpParser(datafile)
    return parser

#def phd_to_csv(datafile, return_df = False):
#    parser = read_picoharp_phd(datafile)
#    name, ext = datafile.rsplit('.', 1)
#    
#    total_curves = parser.no_of_curves()
#    y = []
#    for i in range(total_curves):
#        res, curve = parser.get_curve(i)
#        time_window = int(np.floor(parser.get_time_window_in_ns(curve_no)))
#        curve = curve[0:time_window]
#        y.append(curve)
#    
#    df = pd.DataFrame(y)
#    df.T.to_csv(name+".csv", index=False, header=False)
#    if return_df == True:
#        return df.T

def smooth(curve, boxwidth):
    sm_curve = np.convolve(curve, np.ones(boxwidth)/boxwidth, mode="same")
    return sm_curve

def get_x_y(curve_no, parser, smooth_trace = False, boxwidth = 3):
    
    assert type(parser) == picoharp_phd.PicoharpParser, 'must be picoharp parser'
    res, curve = parser.get_curve(curve_no)
    time_window = int(np.floor(parser.get_time_window_in_ns(curve_no)))
    curve = curve[0:time_window]
    size = len(curve)
    x = np.arange(0, size*res, res, np.float)
    if smooth_trace == True:
        curve = smooth(curve, boxwidth=boxwidth)
    return x,curve
