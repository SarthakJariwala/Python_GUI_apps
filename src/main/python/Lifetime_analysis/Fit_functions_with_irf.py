import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fmin_tnc, differential_evolution
from scipy.special import gamma
from scipy.signal import fftconvolve
from scipy.integrate import odeint


"""Fit TCSPC data to a model by reconvolution with the IRF
Convolution is done in time domain with np.convolve()
Convolution can be done in the frequency domain of np.convolve() is replaced by scipy.signal.fftconvolve()

For a good tutorial on numerical convolution, see section 13.1 of Numerical Recipes:

Press, W. H.; Teukolsky, S. A.; Vetterling, W. T.; Flannery, B. P., 
Numerical Recipes 3rd Edition: The Art of Scientific Computing. 3 ed.; 
Cambridge University Press: New York, 2007

***Note that algorithm given in Numerical Recipes does convolution in the frequency
domain using the FFT. However the discussion of convolution in 13.1 applies to the time
domain and can be used to understand this Python code.***

-MZ, 2/2017
""" 

def convolve_sig_resp(signal_array, response_array, t_array, tstep):

    def normalize_response(response_array, t_array):
        area = np.trapz(response_array, x = t_array)
        return response_array / area        
    
    def array_zeropad_neg(array, pad_length):
    
        return np.pad(array, (pad_length, 0), 'constant', constant_values = (0,0))

    def array_zeropad_pos(array, pad_length):
        return np.pad(array, (0, pad_length), 'constant', constant_values = (0,0))

#    def array_symmetricpad_neg(array, pad_length):
#    
#        return np.pad(array, (pad_length, 0), 'symmetric')

    def signal_and_resp_forconv(signal_array, response_array):
        resp_pad_negtime = array_zeropad_neg(normalize_response(response_array, t_array), len(response_array) - 1)
        sig_pad_negtime = array_zeropad_neg(signal_array, len(signal_array) - 1)
        sig_pad_postime = array_zeropad_pos(sig_pad_negtime, len(response_array))
        return [resp_pad_negtime, sig_pad_postime]    
    
    resp, sig = signal_and_resp_forconv(signal_array, response_array)    
    convolution = tstep * fftconvolve(sig, resp, mode = 'same')#np.convolve(resp, sig, mode = 'same') 

    return convolution[len(signal_array) - 1 : (2*len(signal_array)) - 1]

def convolution_plusnoise(signal_array, response_array, t_array, tstep, noiselevel):
    return convolve_sig_resp(signal_array, response_array, t_array, tstep) + noiselevel

def herz_ode(t, n0, params):
    a = params[0]
    k1 = params[1]
    k2 = params[2]
    k3 = params[3]
    def odefun(n, t, k1, k2, k3):    
        dndt = -(k1 * n) - (k2 * (n ** 2.0)) - (k3 * (n ** 3.0)) 
        return dndt
    ode_sol = odeint(odefun, n0, t, args = (k1, k2, k3))[:,0]
    pl = k2 * (ode_sol ** 2.0)
    return a*pl
       
def fit_herz_ode_global_3traces_fmin_tnc(t1, t2, t3, tstep, d1, d2, d3, irf, init_params, bounds, n0array):
    time_array1 = t1
    time_array2 = t2
    time_array3 = t3
    data_array1 = d1
    data_array2 = d2
    data_array3 = d3
    n0 = n0array[0]
    n1 = n0array[1]
    n2 = n0array[2]    
    def min_fit_decay(params):
        #Minimize chi-squre for data set with Poisson distribution ()        

        a0 = params[0]
        a1 = params[1]
        a2 = params[2]
        k1 = params[3]
        k2 = params[4]
        k3 = params[5]
        noise1 = params[6]
        noise2 = params[7]
        noise3 = params[8]            
        decaymodel1 = herz_ode(time_array1, n0, np.array([a0,k1,k2,k3]))
        decaymodel2 = herz_ode(time_array2, n1, np.array([a1,k1,k2,k3]))
        decaymodel3 = herz_ode(time_array3, n2, np.array([a2,k1,k2,k3]))
        model1 = convolution_plusnoise(decaymodel1, irf, time_array1, tstep, noise1)
        model2 = convolution_plusnoise(decaymodel2, irf, time_array2, tstep, noise2)

        model3 = convolution_plusnoise(decaymodel3, irf, time_array3, tstep, noise3)

        data_fit_idx1 = np.nonzero(data_array1)
        data_fit_idx2 = np.nonzero(data_array2)
        data_fit_idx3 = np.nonzero(data_array3)
        data_array_fit1 = data_array1[data_fit_idx1]
        data_array_fit2 = data_array2[data_fit_idx2]
        data_array_fit3 = data_array3[data_fit_idx3]
        
        model_fit1 = model1[data_fit_idx1]
        model_fit2 = model2[data_fit_idx2]
        model_fit3 = model3[data_fit_idx3]
        
        min1 = np.sum(((data_array_fit1 - model_fit1)** 2.0) / (np.sqrt(data_array_fit1) ** 2.0))
        min2 = np.sum(((data_array_fit2 - model_fit2)** 2.0) / (np.sqrt(data_array_fit2) ** 2.0))
        min3 = np.sum(((data_array_fit3 - model_fit3)** 2.0) / (np.sqrt(data_array_fit3) ** 2.0))
        return np.sqrt((min1 ** 2.0) + (min2 ** 2.0) + (min3 ** 2.0)) 
    bestfit_params = fmin_tnc(min_fit_decay, init_params, approx_grad = True, bounds = bounds)[0]
    def bestfit_decay(params):
        a0 = params[0]
        a1 = params[1]
        a2 = params[2]
        k1 = params[3]
        k2 = params[4]
        k3 = params[5]
        noise1 = params[6]
        # noise2 = params[7]
        # noise3 = params[8]            
        decaymodel1 = herz_ode(time_array, n0, np.array([a0,k1,k2,k3]))
        decaymodel2 = herz_ode(time_array, n1, np.array([a1,k1,k2,k3]))
        decaymodel3 = herz_ode(time_array, n2, np.array([a2,k1,k2,k3]))
        
        model1 = convolution_plusnoise(decaymodel1, irf, time_array, tstep, noise1)
        model2 = convolution_plusnoise(decaymodel2, irf, time_array, tstep, noise2)
        model3 = convolution_plusnoise(decaymodel3, irf, time_array, tstep, noise3)
        return [model1, model2, model3]
            
    bestfit_model = bestfit_decay(bestfit_params)
    # plt.figure()
    # plt.ylabel('PL Counts')
    # plt.xlabel('Time (ns)')
    # plt.semilogy(time_array1, data_array1,'b', label = 'Data')
    # plt.semilogy(time_array1, bestfit_model[0], 'r', label = 'Fit')
    # plt.semilogy(time_array2, data_array2,'b', label = 'Data')
    # plt.semilogy(time_array2, bestfit_model[1], 'r', label = 'Fit')
    # plt.semilogy(time_array3, data_array3,'b', label = 'Data')
    # plt.semilogy(time_array3, bestfit_model[2], 'r', label = 'Fit')
    # plt.legend(loc = 'best')
    return bestfit_params, bestfit_model, data_array, time_array, irf 

def multi_exp(t, params, num_exp):
    exp_array = np.empty((len(t), num_exp))

    i = 0
    while (i<len(params)):
        remain = i % 2
        if remain == 0:
            exp = params[i] * np.exp(-(1.0 / params[i + 1]) * t)
            exp_array[:, int(i / 2)] = exp
        else:
            pass
        i = i + 1
    
    return np.sum(exp_array, axis = 1)
    
def exp_stretch(t, tc, beta, a):
    return (a * np.exp(-((1.0 / tc) * t) ** beta)) 

def avg_tau_from_exp_stretch(tc, beta):
    return (tc / beta) * gamma(1.0 / beta)
#modtest = np.remainder(2,2)        
#print np.remainder(5,2), np.divide(5,2)#2 % 2        
    
def fit_exp_stretch_fmin_tnc(t, tstep, data, irf, init_params, bounds):
    time_array = t
    data_array = data
    def min_fit_decay(params):
        #Minimize chi-squre for data set with Poisson distribution ()        
        tau_c = params[0]
        beta = params[1]
        a = params[2]
        noise = params[3]            
        decaymodel = exp_stretch(time_array, tau_c, beta, a)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        data_fit_idx = np.nonzero(data_array)
        data_array_fit = data_array[data_fit_idx]
        model_fit = model[data_fit_idx]
        return np.sum(((data_array_fit - model_fit)** 2.0) / (np.sqrt(data_array_fit) ** 2.0))
    bestfit_params = fmin_tnc(min_fit_decay, init_params, approx_grad = True, bounds = bounds)[0]
    def bestfit_decay(params):
        tau_c = params[0]
        beta = params[1]
        a = params[2]
        noise = params[3]            
        decaymodel = exp_stretch(time_array, tau_c, beta, a)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        return model
            
    bestfit_model = bestfit_decay(bestfit_params)
    t_avg = avg_tau_from_exp_stretch(bestfit_params[0], bestfit_params[1])
    # print ('--Stretched Exponential Best Fit Parameters--\ntau_avg = %.5f ns\nbeta = %.5f \ntau_c = %.5f ns \na = %.5f \nnoise = %.5f counts' %(t_avg, bestfit_params[1], bestfit_params[0], bestfit_params[2], bestfit_params[3]))
    # plt.figure()
    # plt.ylabel('PL (au)', fontsize = 25)
    # plt.xlabel('Time (ns)', fontsize = 25)
    # plt.semilogy(time_array, data_array,'b', label = 'Data')
    # plt.semilogy(time_array, bestfit_model, 'r', label = 'Fit')
    # plt.legend(loc = 'best')
    return bestfit_params, t_avg, bestfit_model, data_array, time_array, irf  

def fit_exp_stretch_diffev(t, tstep, data, irf,  bounds):
    time_array = t
    data_array = data
    def min_fit_decay(params):
        #Minimize chi-squre for data set with Poisson distribution ()        
        tau_c = params[0]
        beta = params[1]
        a = params[2]
        noise = params[3]            
        decaymodel = exp_stretch(time_array, tau_c, beta, a)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        data_fit_idx = np.nonzero(data_array)
        data_array_fit = data_array[data_fit_idx]
        model_fit = model[data_fit_idx]
        return np.sum(((data_array_fit - model_fit)** 2.0) / (np.sqrt(data_array_fit) ** 2.0))
    bestfit = differential_evolution(min_fit_decay, bounds = bounds, polish = True)
    bestfit_params = bestfit.x    
    def bestfit_decay(params):
        tau_c = params[0]
        beta = params[1]
        a = params[2]
        noise = params[3]            
        decaymodel = exp_stretch(time_array, tau_c, beta, a)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        return model
    bestfit_model = bestfit_decay(bestfit_params)
    t_avg = avg_tau_from_exp_stretch(bestfit_params[0], bestfit_params[1])
    # print ('--Stretched Exponential Best Fit Parameters--\ntau_avg = %.5f ns\nbeta = %.5f \ntau_c = %.5f ns \na = %.5f \nnoise = %.5f counts' %(t_avg, bestfit_params[1], bestfit_params[0], bestfit_params[2], bestfit_params[3]))
    # plt.figure()
    # plt.ylabel('PL (au)', fontsize = 25)
    # plt.xlabel('Time (ns)', fontsize = 25)
    # plt.semilogy(time_array, data_array,'b', label = 'Data')
    # plt.semilogy(time_array, bestfit_model, 'r', label = 'Fit')
    # plt.legend(loc = 'best')
    return bestfit_params, t_avg, bestfit_model, data_array, time_array, irf

def fit_multi_exp_fmin_tnc(t, tstep, data, irf, init_params, bounds, n):
    time_array = t
    data_array = data
    def min_fit_decay(params):
        #Minimize chi-squre for data set with Poisson distribution ()        

        noise = params[len(params) - 1]            
        decaymodel = multi_exp(time_array, params[:-1], n)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        data_fit_idx = np.nonzero(data_array)
        data_array_fit = data_array[data_fit_idx]
        model_fit = model[data_fit_idx]
        return np.sum(((data_array_fit - model_fit)** 2.0) / (np.sqrt(data_array_fit) ** 2.0))
    bestfit_params = fmin_tnc(min_fit_decay, init_params, approx_grad = True, bounds = bounds)[0]    
    def bestfit_decay(params):

        noise = params[len(params) - 1]                
        decaymodel = decaymodel = multi_exp(time_array, params[:-1], n)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        return model
    bestfit_model = bestfit_decay(bestfit_params)
    
    # print ('--Multi Exponential Best Fit Parameters--')
    # print ('# Exponentials = %.0f' %(n))
    
    # for i in range(len(bestfit_params) - 1):
    #     remain = i % 2
          
    #     if remain == 0:
    #         print ('a%.0f = %.5f \ntau%.0f = %.5f ns' %((i / 2) + 1, bestfit_params[i], (i / 2) + 1, bestfit_params[i + 1]))
            
    #     else:
    #         pass
        
    # print ('Noise = %.5f Counts'%(bestfit_model[len(bestfit_model) - 1]))
   
    # plt.figure()
    # plt.ylabel('PL Counts')
    # plt.xlabel('Time (ns)')
    # plt.semilogy(time_array, data_array,'b', label = 'Data')
    # plt.semilogy(time_array, bestfit_model, 'r', label = 'Fit')
    # plt.legend(loc = 'best')
    return bestfit_params, bestfit_model, data_array, time_array, irf  

def fit_multi_exp_diffev(t, tstep, data, irf,  bounds, n):
    time_array = t
    data_array = data
    def min_fit_decay(params):
        #Minimize chi-squre for data set with Poisson distribution ()        

        noise = params[len(params) - 1]            
        decaymodel = multi_exp(time_array, params[:-1], n)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        data_fit_idx = np.nonzero(data_array)
        data_array_fit = data_array[data_fit_idx]
        model_fit = model[data_fit_idx]
        return np.sum(((data_array_fit - model_fit)** 2.0) / (np.sqrt(data_array_fit) ** 2.0))
    bestfit = differential_evolution(min_fit_decay, bounds = bounds, polish = True)
    bestfit_params = bestfit.x    
    def bestfit_decay(params):

        noise = params[len(params) - 1]                
        decaymodel = decaymodel = multi_exp(time_array, params[:-1], n)
        model = convolution_plusnoise(decaymodel, irf, time_array, tstep, noise)
        return model
    bestfit_model = bestfit_decay(bestfit_params)
    
    # print ('--Multi Exponential Best Fit Parameters--')
    # print ('# Exponentials = %.0f' %(n))
    
    # for i in range(len(bestfit_params) - 1):
    #     remain = i % 2
          
    #     if remain == 0:
    #         print ('a%.0f = %.5f \ntau%.0f = %.5f ns' %((i / 2) + 1, bestfit_params[i], (i / 2) + 1, bestfit_params[i + 1]))
            
    #     else:
    #         pass
        
    # print ('Noise = %.5f Counts'%(bestfit_model[len(bestfit_model) - 1]))
   
    # plt.figure()
    # plt.ylabel('PL Counts')
    # plt.xlabel('Time (ns)')
    # plt.semilogy(time_array, data_array,'b', label = 'Data')
    # plt.semilogy(time_array, bestfit_model, 'r', label = 'Fit')
    # plt.legend(loc = 'best')
    return bestfit_params, bestfit_model, data_array, time_array, irf  