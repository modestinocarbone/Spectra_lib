
"""
Created on Fri March 21 2025

@author: Modestino Carbone


"""

import numpy as np
import scipy.signal as signal
import scipy.fft as fourier

def psd_bartlett(x, fs, avg, window, N):
    
    #window genaration
    win = signal.windows.get_window(window,N)
    
    if (len(x)/N >= avg):
        
        psd_buff = []
        
        for i in range(avg):
            
            X = fourier.fft(x[i*N:i*N+N]*win) 
            app = 2*(1/fs)*X*np.conj(X)/((win*win).sum())
            psd_buff.append(app)
            
        buff = np.array(psd_buff)

        PSD = (np.sum(buff , axis=0)/avg)
        
    else:
        
        X = fourier.fft(x)
        PSD = 2*(1/fs)*X*np.conj(X)/((win*win).sum())
        
    freq = np.linspace(0,fs,N)
    
    return freq[0:int(N/2)], PSD[0:int(N/2)]  



def psd_welch(x, fs, avg, window, N):
    
    #window genaration
    win = signal.windows.get_window(window,N)
    
    if (len(x)/N >= avg):
        
        psd_buff = []
        
        for i in range(avg*2 -1):
            
            X = fourier.fft(x[i*int(N/2):i*int(N/2)+N]*win) 
            app = 2*(1/fs)*X*np.conj(X)/((win*win).sum())
            psd_buff.append(app)
            
            
        buff = np.array(psd_buff)
        
        PSD = (np.sum(buff , axis=0)/(avg*2 -1))
        
    else:
        
        X = fourier.fft(x)
        PSD = 2*(1/fs)*X*np.conj(X)/((win*win).sum())
        
    freq = np.linspace(0,fs,N)
    
    return freq[0:int(N/2)], PSD[0:int(N/2)]  





def psd_MA(x, fs, Max_avg, N, win_type): 
       
    #Power Specral Density Bertlett with Maximum Average
    PSD = []
    fr = []
    j=0
    averages=[]
    decades=[]
    
    while(True):
        
        Window=int(N/(Max_avg/(10**j)))
        win = signal.windows.get_window(win_type,Window)
        
        buff = []
        buff2 = []
        csd_cross= []
        
        averages.append((Max_avg/(10**j)))
        
        for i in range(int(Max_avg/(10**j))):

            X = fourier.fft(x[i*Window:i*Window+Window]*win)
            app= 2*(1/fs)*X*np.conj(X)/((win*win).sum())
            csd_cross.append(app)

        buff = np.array(csd_cross)
        buff2 = (np.sum(buff , axis=0)/(Max_avg/(10**j)))
        freq =[]

        freq = np.linspace(0,fs,int(len(x)/(Max_avg/(10**j))))

        f_min = (fs/2)/(10**(j+1))
        f_max = (fs/2)/(10**(j))
        
        decades.append(f_max)
        
        j=j+1
        
        if  Window == N:
            
            f_min = 0
            mask = (freq >= f_min) & (freq <= f_max)
            
            PSD.extend(buff2[mask])
            fr.extend(freq[mask])
            
            #append the last frequncy value when the maximum is reached
            decades.append(0)
            
            #final index sorting
            sorted_indices = np.argsort(fr)
            fr = np.array(fr)[sorted_indices]
            PSD = np.array(PSD)[sorted_indices]
            
            return fr[0:int(N/2)], np.array(PSD[0:int(N/2)]), np.array(averages), np.array(decades)
        
        mask = (freq >= f_min) & (freq <= f_max)
              
        PSD.extend(buff2[mask])
        fr.extend(freq[mask])
       
    #final index sorting
    sorted_indices = np.argsort(fr)
    fr = np.array(fr)[sorted_indices]
    PSD = np.array(PSD)[sorted_indices]
    
    return fr[0:int(N/2)], np.array(PSD[0:int(N/2)]), np.array(averages), np.array(decades)




def psd_Min_Max( x, fs, Min_avg, Max_avg, N, win_type): 
       
    PSD = []
    fr = []
    j=0
    averages=[]
    decades=[]
    
    while(True):
        
        Window=int(len(x)/(Max_avg/(10**j)))
        win = signal.windows.get_window(win_type,Window)
        
        buff = []
        buff2 = []
        psd_cross= []
        
        averages.append((Max_avg/(10**j)))
        
        for i in range(int(Max_avg/(10**j))):

            X = fourier.fft(x[i*Window:i*Window+Window]*win)
            app= 2*(1/fs)*X*np.conj(X)/((win*win).sum())
            psd_cross.append(app)

        buff = np.array(psd_cross)
        buff2 = (np.sum(buff , axis=0)/(Max_avg/(10**j)))
        freq =[]

        freq = np.linspace(0,fs,int(len(x)/(Max_avg/(10**j))))
        
        f_min = (fs/2)/(10**(j+1))
        f_max = (fs/2)/(10**(j))
        
        decades.append(f_max)
                
        j=j+1
        
        if  Window == int(len(x)/(Min_avg)):
            
            f_min = 0
            mask = (freq >= f_min) & (freq <= f_max)
            
            PSD.extend(buff2[mask])
            fr.extend(freq[mask])
            
            #append the last frequncy value when the maximum is reached
            decades.append(0)
            
            #final index sorting
            sorted_indices = np.argsort(fr)
            fr = np.array(fr)[sorted_indices]
            PSD = np.array(PSD)[sorted_indices]
            
            return fr[0:int(N/2)], np.array(PSD[0:int(N/2)]), np.array(averages), np.array(decades)
       
        mask = (freq >= f_min) & (freq <= f_max)
        
        PSD.extend(buff2[mask])
        fr.extend(freq[mask])
       
    #final index sorting
    sorted_indices = np.argsort(fr)
    fr = np.array(fr)[sorted_indices]
    PSD = np.array(PSD)[sorted_indices]
    
    return fr[0:int(N/2)], np.array(PSD[0:int(N/2)]), np.array(averages), np.array(decades)


