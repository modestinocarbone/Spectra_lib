
"""
Created on Fri March 21 2025

@author: Modestino Carbone


"""

import numpy as np
import scipy.signal as signal
import scipy.fft as fourier


def fft(x, fs, window, N):
    
    #window genaration
    win = signal.windows.get_window(window,N)
    
    #scipy fft algorithm 
    sp = fourier.fft(x[0:N]*win)
    freq = np.linspace(0,fs,N)
    
    #normalization of the scale 2/N factor
    fft = 2*np.abs(sp)/(win.sum())
    
    #normalizing DC
    fft[0]=fft[0]/2

    return freq[0:int(N/2)], fft[0:int(N/2)]
    
 