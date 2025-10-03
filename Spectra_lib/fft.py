"""
 *  File:        fft.py
 *  Author:      Modestino Carbone
 *  Created on:  21/03/2025
 *  Last update: 01/10/2025
 *
 *  Description:
 *      A simple fast Fourier transform (FFT) implementation, just for kids.
 *
"""

import numpy as np
import scipy.signal as signal
import scipy.fft as fourier

def fft(x, fs, window, N):
    
    # Window genaration
    win = signal.windows.get_window(window, N)
    
    # Check length
    if(N > len(x)):
        raise Exception(f"Window length N: {N} is greater than actual input dimension: {len(x)}")
      
    # 'scipy' FFT algorithm
    # 'scipy' is faster than 'numpy'
    sp = fourier.fft(x[0:N]*win)
    
    # Frequency axis
    freq = np.linspace(0, fs/2, N//2)
    
    # Normalization of the scale with the 2/N factor
    # This function retrieves only the magnitude of FFT
    # Change this part to convert the FFT output to an array of complex values
    fft = 2*sp/(win.sum())
    fft = fft[0:N//2]
    
    # Normalizing DC
    fft[0]=fft[0]/2

    return freq, np.abs(fft)