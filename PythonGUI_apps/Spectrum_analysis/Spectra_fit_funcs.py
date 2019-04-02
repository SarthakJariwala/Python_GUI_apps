# -*- coding: utf-8 -*-
"""
Created on Sun Mar 31 15:46:54 2019

@author: Sarthak
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Jan 31 20:48:03 2019

@author: Sarthak
"""
import numpy as np
from lmfit.models import GaussianModel

class Spectra_Fit(object):
    """
    Fit a spectrum.
    
    Attributes:
        data: spectrum data (x-axis and y-axis)
        ref: reference spectrum (both x and y-axis) for background correction
    """
    
    def __init__(self, data, ref):
        self.data = data
        self.ref = ref
        
    def background_correction(self):
        """Return the background corrected spectrum"""
        x = self.data[:, 0] # wavelengths
        y = self.data[:, 1] # intensity
        yref = self.ref[:, 1]
        y = y - yref # background correction
        # y = y - np.mean(y[(x>600) & (x<700)]) # removing any remaining bckgrnd
        return [x,y]

class Single_Gaussian(Spectra_Fit):
    """Fit a single gaussian to the spectrum, plot it and save it
    
    Attributes:
        data: spectrum data (x-axis and y-axis)
        ref: reference spectrum (both x and y-axis) for background correction
    """
    
    def gaussian_model(self):
        x,y = self.background_correction()
        gmodel = GaussianModel(prefix = 'g1_') # calling gaussian model
        pars = gmodel.guess(y, x=x) # parameters - center, width, height
        result = gmodel.fit(y, pars, x=x, nan_policy='propagate')
        return result

class Double_Gaussian(Spectra_Fit):
    """Fit two gaussians to the spectrum
    
    Attributes:
        data: spectrum data (x-axis and y-axis)
        ref: reference spectrum (both x and y-axis) for background correction
    """
    
    def gaussian_model(self):
        x,y = self.background_correction()
        gmodel_1 = GaussianModel(prefix='g1_') # calling gaussian model
        pars = gmodel_1.guess(y, x=x) # parameters - center, width, height
        pars['g1_center'].set(800, min = 795, max = 820)
        pars['g1_sigma'].set(15)
        pars['g1_amplitude'].set(min=0)

        gmodel_2 = GaussianModel(prefix='g2_')
        pars.update(gmodel_2.make_params()) # update parameters - center, width, height
        pars['g2_center'].set(767, min = 760, max = 775)
        pars['g2_amplitude'].set(min = 0)
        pars['g2_sigma'].set(min = pars['g1_sigma'].value)

        gmodel = gmodel_1 + gmodel_2
        result = gmodel.fit(y, pars, x=x, nan_policy='propagate')
        return result

class Multi_Gaussian(Spectra_Fit):
    
    def __init__(self, data, ref, num_of_gaussians, peak_pos, min_max_range):
        Spectra_Fit.__init__(self, data, ref)
        self.num_of_gaussians = num_of_gaussians
        self.peak_pos = peak_pos
        self.min_max_range = min_max_range
        
    def multi_gaussian(self):
        composite_model = None
        composite_pars = None
        
        x,y = self.background_correction()
        
        assert self.num_of_gaussians == len(self.peak_pos), ("Number of gaussians must be equal to the number of peak positions")
        assert len(self.min_max_range) == len(self.peak_pos), ("Number of bounds on the range must be equal to the number of peak positions")

        
        for i in range(self.num_of_gaussians):

            model = GaussianModel(prefix='g'+str(i+1)+'_')

            if composite_pars is None:
                composite_pars = model.guess(y, x=x)
#                 composite_pars = model.make_params()
                composite_pars['g'+str(i+1)+'_center'].set(self.peak_pos[i], 
                                                           min = self.min_max_range[0][0], max = self.min_max_range[0][1])
                composite_pars['g'+str(i+1)+'_sigma'].set(15)
                composite_pars['g'+str(i+1)+'_amplitude'].set(min = 0)
                                                
            else:
                composite_pars.update(model.make_params())
                composite_pars['g'+str(i+1)+'_center'].set(self.peak_pos[i],
                                                          min = self.min_max_range[1][0], max = self.min_max_range[1][1])
                composite_pars['g'+str(i+1)+'_sigma'].set(min = composite_pars['g1_sigma'].value)
                composite_pars['g'+str(i+1)+'_amplitude'].set(min = 0)


            if composite_model is None:
                composite_model = model
            else:
                composite_model += model
        
        result = composite_model.fit(y, composite_pars, x=x, nan_policy='propagate')
        return result