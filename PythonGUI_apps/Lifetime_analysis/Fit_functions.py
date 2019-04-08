# -*- coding: utf-8 -*-
"""
Created on Fri Mar 29 11:21:59 2019

@author: Sarthak
"""

import numpy as np
from scipy.optimize import differential_evolution
from scipy.special import gamma

def stretch_exp_fit(TRPL, t, Tc = (0,1e5), Beta = (0,1), A = (0,1e6)):

    def exp_stretch(t, tc, beta, a):
        return (a * np.exp(-((1.0 / tc) * t) ** beta)) 
    
    def avg_tau_from_exp_stretch(tc, beta):
        return (tc / beta) * gamma(1.0 / beta)
    
    def Diff_Ev_Fit_SE(TRPL):
        TRPL = TRPL
        
        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tc = params[0]
            beta = params[1]
            a = params[2]
            
            PL_sim = exp_stretch(t,tc,beta,a)
    
            Resid= np.sqrt(np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data
        
        bounds = [Tc, Beta, A] 
    
        result = differential_evolution(residuals, bounds)
        return result.x
    
    p = Diff_Ev_Fit_SE(TRPL)

    tc = p[0]
    beta = p[1]
    a = p[2]
    
    PL_fit = exp_stretch(t,tc,beta,a)
    
    avg_tau = avg_tau_from_exp_stretch(tc,beta)
    
    return tc, beta, a, avg_tau, PL_fit

def double_exp_fit(TRPL, t, tau1_bounds=(0,100), a1_bounds=(0,1e5), tau2_bounds=(0,100), a2_bounds=(0,1e5)):

    def single_exp(t, tau, a):
        return (a * np.exp(-((1.0 / tau)*t) ))
    
    def double_exp(t, tau1, a1, tau2, a2):
        return ((single_exp(t, tau1, a1)) + (single_exp(t, tau2, a2)))
    
    def avg_tau_from_double_exp(tau1, a1, tau2, a2):
        return (((tau1*a1) + (tau2*a2))/(a1+a2))
    
    def Diff_Ev_Fit_DE(TRPL):
        TRPL = TRPL
        
        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tau1 = params[0]
            a1 = params[1]
            tau2 = params[2]
            a2 = params[3]
            
            PL_sim = double_exp(t,tau1, a1, tau2, a2)
    
            Resid= np.sqrt(np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data
        
        bounds = [tau1_bounds, a1_bounds, tau2_bounds, a2_bounds] 
    
        result = differential_evolution(residuals, bounds)
        return result.x
    
    p = Diff_Ev_Fit_DE(TRPL)

    tau1 = p[0]
    a1 = p[1]
    tau2 = p[2]
    a2 = p[3]
    
    PL_fit = double_exp(t, tau1, a1, tau2, a2)
    
    avg_tau = avg_tau_from_double_exp(tau1, a1, tau2, a2)
    
    return tau1, a1, tau2, a2, avg_tau, PL_fit

def single_exp_fit(TRPL, t, tau_bounds=(0,1000), a_bounds=(0,1e5)):
    
    def single_exp(t, tau, a):
        return (a * np.exp(-((1.0 / tau)*t) ))
    
    def Diff_Ev_Fit_singleExp(TRPL):
        TRPL = TRPL
        
        def residuals(params):#params are the parameters to be adjusted by differential evolution or leastsq, interp is the data to compare to the model.
            #Variable Rates
            tau = params[0]
            a = params[1]
            
            PL_sim = single_exp(t, tau, a)
    
            Resid= np.sqrt(np.sum(((PL_sim-TRPL)**2)/(np.sqrt(PL_sim)**2)))
            return Resid #returns the difference between the PL data and simulated data
        
        bounds = [tau_bounds, a_bounds] 
    
        result = differential_evolution(residuals, bounds)
        return result.x
    
    p = Diff_Ev_Fit_singleExp(TRPL)

    tau = p[0]
    a = p[1]
    
    PL_fit = single_exp(t, tau, a)
    
    return tau, a, PL_fit
    