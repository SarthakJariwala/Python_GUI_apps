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
from lmfit.models import GaussianModel, LorentzianModel

class Spectra_Fit(object):
    """
    Fit a spectrum.
    
    Attributes:
        data: spectrum data (x-axis and y-axis)
        ref: reference spectrum (both x and y-axis) for background correction
    """
    
    def __init__(self, data, ref, wlref = None):
        self.data = data
        self.ref = ref
        self.wlref = wlref
        
    def background_correction(self):
        """Return the background corrected spectrum"""
        x = self.data[:, 0] # wavelengths
        y = self.data[:, 1] # intensity
        yref = self.ref[:, 1]
        y = y - yref # background correction
        if self.wlref is not None:
            wlref = self.wlref[:,1]
            y = y/wlref
            
        # y = y - np.mean(y[(x>600) & (x<700)]) # removing any remaining bckgrnd
        return [x,y]

class Single_Gaussian(Spectra_Fit):
    """Fit a single gaussian to the spectrum
    
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
    
    def gaussian_model_w_lims(self, center_initial_guess=None, sigma_initial_guess=None, center_min=None, center_max=None):
        x,y = self.background_correction()
        gmodel = GaussianModel(prefix = 'g1_') # calling gaussian model
        pars = gmodel.guess(y, x=x) # parameters - center, width, height
        pars['g1_center'].set(center_initial_guess, min=center_min, max=center_max)
        pars['g1_sigma'].set(sigma_initial_guess)
        result = gmodel.fit(y, pars, x=x, nan_policy='propagate')
        return result #770 760 780   sigma 15 

class Single_Lorentzian(Spectra_Fit):
    """Fit a single Lorentzian to the spectrum
    
    Attributes:
        data: spectrum data (x-axis and y-axis)
        ref: reference spectrum (both x and y-axis) for background correction
    """
    
    def lorentzian_model(self):
        x,y = self.background_correction()
        lmodel = LorentzianModel(prefix = 'l1_') # calling lorentzian model
        pars = lmodel.guess(y, x=x) # parameters - center, width, height
        result = lmodel.fit(y, pars, x=x, nan_policy='propagate')
        return result
    
    # def lorentzian_model_w_lims(self, center_min = None, center_max = None):
    #     x,y = self.background_correction()
    #     lmodel = LorentzianModel(prefix = 'l1_') # calling lorentzian model
    #     pars = lmodel.guess(y, x=x) # parameters - center, width, height
    #     pars['l1_center'].set(min = center_min, max = center_max)
    #     result = lmodel.fit(y, pars, x=x, nan_policy='propagate')
    #     return result

    def lorentzian_model_w_lims(self, center_initial_guess=None, sigma_initial_guess=None, center_min = None, center_max = None):
        x,y = self.background_correction()
        lmodel = LorentzianModel(prefix = 'l1_') # calling lorentzian model
        pars = lmodel.guess(y, x=x) # parameters - center, width, height
        pars['l1_center'].set(center_initial_guess, min = center_min, max = center_max)
        pars['l1_sigma'].set(sigma_initial_guess)
        result = lmodel.fit(y, pars, x=x, nan_policy='propagate')
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

        gmodel_2 = GaussianModel(prefix='g2_')
        pars.update(gmodel_2.make_params()) # update parameters - center, width, height

        gmodel = gmodel_1 + gmodel_2
        result = gmodel.fit(y, pars, x=x, nan_policy='propagate')
        return result

    def gaussian_model_w_lims(self, center_initial_guesses=None, sigma_initial_guesses=None, min_max_range=None):
        #center_initial_guesses - list containing initial guesses for peak centers. [center_guess1, center_guess2]
        #sigma_initial_guesses - list containing initial guesses for sigma. [sigma1, sigma2]
        #min_max_range - list containing lists of min and max for peak center. [ [min1, max1], [min2, max2] ] 
        
        x,y = self.background_correction()
        gmodel_1 = GaussianModel(prefix='g1_') # calling gaussian model
        pars = gmodel_1.guess(y, x=x) # parameters - center, width, height
        pars['g1_center'].set(center_initial_guesses[0], min = min_max_range[0][0], max = min_max_range[0][1])
        pars['g1_sigma'].set(sigma_initial_guesses[0])
        pars['g1_amplitude'].set(min=0)

        gmodel_2 = GaussianModel(prefix='g2_')
        pars.update(gmodel_2.make_params()) # update parameters - center, width, height
        pars['g2_center'].set(center_initial_guesses[1], min = min_max_range[1][0], max = min_max_range[1][1])
        pars['g2_sigma'].set(sigma_initial_guesses[1])
        pars['g2_amplitude'].set(min = 0)

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
                                                          min = self.min_max_range[i][0], max = self.min_max_range[i][1])
                composite_pars['g'+str(i+1)+'_sigma'].set(min = composite_pars['g1_sigma'].value)
                composite_pars['g'+str(i+1)+'_amplitude'].set(min = 0)


            if composite_model is None:
                composite_model = model
            else:
                composite_model += model
        
        result = composite_model.fit(y, composite_pars, x=x, nan_policy='propagate')
        return result