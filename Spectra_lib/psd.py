"""
 *  File:        psd.py
 *  Author:      Modestino Carbone
 *  Created on:  21/03/2025
 *  Last update: 01/10/2025
 *
 *  Description:
 *      - 'psd_bartlett' and 'psd_welch': The main basics of Power spectral density (PSD).
 *      - 'psd_Min_Max': Power spectral density (CSD) with
 *      advanced averaging options.
 *
"""

import numpy as np
import scipy.signal as signal
import scipy.fft as fourier

def psd_bartlett(x, fs, avg, window, N):

    # Check length
    if(N > len(x)):
        raise Exception(f"Window length N: {N} is greater than actual input dimension: {len(x)}")
      
    # Verify that N is a multiple of len(x)
    if(len(x) % N != 0):
        raise Exception(f"input length: {len(x)} is not a multiple of N: {N}")
        
    # verify number of averages
    if (len(x) // N < avg):
       raise Exception(f"Too many averages is: {avg}; actual maximum is: {len(x)/N}")

    # Window genaration
    win = signal.windows.get_window(window,N)

    # Buffer of PSD that later will be averaged
    psd_buff = []
    
    # Averaging loop
    for i in range(avg):
        
        X = fourier.fft(x[i*N:i*N+N]*win)
        app = 2*(1/fs)*X*np.conj(X)/((win*win).sum())
        psd_buff.append(app)
    
    # Converts out buffer
    buff = np.array(psd_buff)

    # Final average
    PSD = (np.sum(buff , axis=0)/avg)
    
    # Frequency axis
    freq = np.linspace(0,fs,N)
    
    return freq[1:N//2], PSD[1:N//2]

def psd_welch(x, fs, avg, window, N):
    
    # Check length
    if(N > len(x)):
        raise Exception(f"Window length N: {N} is greater than actual input dimension: {len(x)}")
      
    # Verify that len(x) is a multiple of N
    if(len(x) % N != 0):
        raise Exception(f"input length: {len(x)} is not a multiple of N: {N}")
        
    # verify number of averages
    if (len(x) // N < avg):
       raise Exception(f"Too many averages is: {avg}; actual maximum is: {len(x)/N}")

    # Window genaration
    win = signal.windows.get_window(window,N)

    # Buffer of PSD that later will be averaged
    psd_buff = []
    
    # Averaging loop
    for i in range(avg*2 -1):
        
        X = fourier.fft(x[i*int(N/2):i*int(N/2)+N]*win) 
        app = 2*(1/fs)*X*np.conj(X)/((win*win).sum())
        psd_buff.append(app)
        
    # Converts out buffer
    buff = np.array(psd_buff)

    # Final average
    PSD = (np.sum(buff , axis=0)/(avg*2 -1))
    
    # Frequency axis
    freq = np.linspace(0,fs,N)
    
    return freq[1:N//2], PSD[1:N//2]

def psd_Min_Max( x, fs, Min_avg, Max_avg, win_type): 

    # Verify that len(x) is a multiple of Min_avg
    if(len(x) % Min_avg != 0):
        raise Exception(f"input length: {len(x)} is not a multiple of the minimum averages requested: {Min_avg}")
    
    # Verify that len(x) is a multiple of Max_avg
    if(len(x) % Max_avg != 0):
        raise Exception(f"input length: {len(x)} is not a multiple of the maximum averages requested: {Max_avg}")
    
    # check average coherence
    if(Max_avg < Min_avg):
        raise Exception(f"Minimum averages parameter: {Min_avg} is greater than maximum limit: {Max_avg}")
    
    # final vector which contains computed spectrum
    PSD = []
    # final vector of frequency
    fr = []
    # loop index
    j=0
    # Vectors of averages and decades boundaries
    averages=[]
    decades=[]
    
    while(True):
        
        # cration of window
        Window = int(len(x)/(Max_avg/(10**j)))
        win    = signal.windows.get_window(win_type,Window)
        
        buff      = []
        buff2     = []
        csd_cross = []
        
        averages.append((Max_avg/(10**j)))
        
        for i in range(int(Max_avg/(10**j))):

            X = fourier.fft(x[i*Window:i*Window+Window]*win)
            app= 2*(1/fs)*X*np.conj(X)/((win*win).sum())
            csd_cross.append(app)

        buff = np.array(csd_cross)
        buff2 = (np.sum(buff , axis=0)/(Max_avg/(10**j)))

        freq =[]
        freq = np.linspace(0,fs,int(len(x)/(Max_avg/(10**j))))
        
        # compute max/min frequency
        f_min = (fs/2)/(10**(j+1))
        f_max = (fs/2)/(10**(j))
        
        # append the max frequency value in decades vector
        decades.append(f_max)

        # increment j       
        j=j+1
        
        # check if window dimension is equal to the smallest window possible
        if  Window == len(x)//(Min_avg):
            
            f_min = 0
            mask = (freq >= f_min) & (freq <= f_max)
            
            PSD.extend(buff2[mask])
            fr.extend(freq[mask])
            
            # append the last frequncy value when the maximum is reached
            decades.append(0)
            
            # final index sorting
            sorted_indices = np.argsort(fr)
            fr = np.array(fr)[sorted_indices]
            PSD = np.array(PSD)[sorted_indices]
            
            return fr[1:], PSD[1:], np.array(averages), np.array(decades)
       
        mask = (freq >= f_min) & (freq <= f_max)
        
        PSD.extend(buff2[mask])
        fr.extend(freq[mask])

def psd_composite(x, fs, avg, win_type, type):

    # buffer output array
    fr = []
    PSD = []

    for i in range(len(x)):

        # check if average is possible    
        if(len(x[i]) % avg != 0):
            raise Exception(f"input length: {len(x[i])} is not a multiple of avg: {avg}")
        
        # just remenber that resolution is fs/W 
        # where W is the length of each subwindow
        W = len(x[i])//avg

        # just for initialization
        f   = [] 
        csd = [] 

        # selected operation depends on the string type
        if type == "bartlett":
            f, csd  = psd_bartlett(x[i], fs[i], avg, win_type, W)  
        if type == "welch":
            f, csd  = psd_welch(x[i], fs[i], avg, win_type, W)  
        
        # mask is from 10**(int(np.log10(fs[i]/2))-1) 10**int(np.log10(fs[i]/2))
        mask = (f >= 10**(int(np.log10(fs[i]/2))-1)) & (f <= 10**int(np.log10(fs[i]/2)))
        fr.extend(f[mask])
        PSD.extend(csd[mask])
    
    return fr, PSD