
"""
Created on Fri Jun 6 2024

@author: Modestino Carbone

"""

import numpy as np
from scipy import signal
import matplotlib.pyplot as plt


class Spectra(object):

    def __init__(self):
        super(Spectra,self).__init__()

#general purpose spectrum functions################################################################################################

    def Fft_custom(self,x,fs,N = None,window = 'hann'):
        
        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided
        T=N*1/fs
        win = signal.windows.get_window(window,N)
        
        t = np.arange(0,T,1/fs)
        sp = np.fft.fft(x)
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)
        fft = 2*np.abs(sp)/((win*win).sum())

        return freq[0:int(N/2)], fft[0:int(N/2)]


###### PSD (POWER SPECTRAL DENSISITY) ################# V/rad(Hz)

    def Psd_basic(self,x,fs,N = None):

        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        X = np.fft.fft(x)
        PSD= 2*(1/fs)*X*np.conj(X)/(N)
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)

        return freq[0:int(N/2)], PSD[0:int(N/2)]


    def Psd_windowed(self,x,fs,N = None,window = 'hann'):
        
        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        win = signal.windows.get_window(window,N)
        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        X = np.fft.fft(x*win)
        PSD= 2*(1/fs)*X*np.conj(X)/( (win*win).sum())
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)

        return freq[0:int(N/2)], PSD[0:int(N/2)]

    def Psd_Bertlett(self,x,fs,avg,N = None,window = 'hann'):
    
        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        win = signal.windows.get_window(window,N)
        
        if (len(x)/N >= avg):
            
            psd_buff = []
            
            for i in range(avg):
                
                X = np.fft.fft(x[i*N:i*N+N]*win)
                app = 2*(1/fs)*X*np.conj(X)/((win*win).sum())
                psd_buff.append(app)
                
            buff = np.array(psd_buff)
            PSD = (np.sum(buff , axis=0)/avg)
            
        else:
            
            X = np.fft.fft(x)
            PSD = 2*(1/fs)*X*np.conj(X)/(N)
            
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)
        
        return freq[0:int(N/2)], PSD[0:int(N/2)]

    def Asd_Bertlett(self,x,fs,avg,N = None,window = 'hann'): #amplitude spectral density

        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        freq, PSD = self.Psd_Bertlett(x, fs, N, avg, window)
        ASD = np.sqrt(PSD)
        
        return freq, ASD



###### CSD (POWER SPECTRAL DENSISITY) ################# V**2/Hz

    def Csd_basic(self,x,y,fs,avg,N = None,window = 'hann'):

        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        win = signal.windows.get_window(window,N)
        
        if (len(x)/N >= avg):
        
            csd_cross=[]
            
            for i in range(avg):
            
                X = np.fft.fft(x[i*N:i*N+N]*win)
                Y = np.fft.fft(y[i*N:i*N+N]*win)
                app= 2*(1/fs)*Y*np.conj(X)/((win*win).sum())
                csd_cross.append(app)
                
            buff = np.array(csd_cross)
            CSD = (np.sum(buff , axis=0)/avg)
            
        else:
        
            X = np.fft.fft(x)
            Y = np.fft.fft(y)
            CSD= 2*(1/fs)*Y*np.conj(X)/(N)
            
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)

        return freq[0:int(N/2)], CSD[0:int(N/2)] 
 
        

    def Coherence(self, x, y, fs, N=None, avg=1, window='hann'):
       
        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        freq, CSD = self.Csd_basic(x, y, fs, avg,  N,window)
        
        _, PSD_x = self.Psd_Bertlett(x, fs, avg, N, window)

        _, PSD_y = self.Psd_Bertlett(y, fs, avg, N, window)

        coherence = (np.abs(CSD) ** 2) / (PSD_x * PSD_y)
        
        return freq, coherence 
        
        
        
#advanced functions##############################################################################

    def PSD_avg_decades(self,x,fs,avg,N = None):
    
        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        dec = 5
        PSD = []
        fr = []
        
        if (len(x)/N >= avg ):
    
            print("performing cross spectrum with "+str(avg)+" means")
            
            for j in range(dec):
                
                Window=int(N/(2**j))
                
                buff = []
                buff2 = []
                psd_cross= []
                
                for i in range(avg):

                    X = np.fft.fft(x[i*Window:i*Window+Window])
                    Y = X
                    app= 2*(1/fs)*Y*np.conj(X)/(Window)
                    psd_cross.append(app)
        
                buff = np.array(psd_cross)
                buff2 = (np.sum(buff , axis=0)/avg)
                freq =[]
                T=Window*(1/fs)
                t = np.arange(0,T,1/fs)
                freq = np.fft.fftfreq(t.shape[-1],d=1/fs)
                
                f_min = (10**j)-1
                f_max = (10**(j+1))-10**j
                
                mask = (freq >= f_min) & (freq <= f_max)
                
           
                PSD.extend(buff2[mask])
                fr.extend(freq[mask])
                print('between'+str(10**(j+1))+'and'+str(10**j)+'there are point of freq n:'+str(len(freq[mask])))
                print('resolution is'+str( (10**(j+1))-10**j/len(freq[mask])))
            
        else:
            X = np.fft.fft(x)
            Y = X
            PSD= 2*(1/fs)*Y*np.conj(X)/(N)
        
        print(fr)
        return fr[0:int(N/2)], np.array(PSD[0:int(N/2)])




    def cross_spectrum_avg_decades(self,x,y,fs,avg,N = None):
    
        if N is None:
            N = len(x)  # Use the length of the signal if N is not provided

        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        dec = 5
        CSD = []
        fr = []
        
        if (len(x)/N >= avg ):
    
            print("performing cross spectrum with "+str(avg)+" means")
            
            for j in range(dec):
                
                Window=int(N/(2**j))
                
                buff = []
                buff2 = []
                csd_cross= []
                
                for i in range(avg):

                    X = np.fft.fft(x[i*Window:i*Window+Window])
                    Y = np.fft.fft(y[i*Window:i*Window+Window])
                    app= 2*(1/fs)*Y*np.conj(X)/(Window)
                    csd_cross.append(app)
        
                buff = np.array(csd_cross)
                buff2 = (np.sum(buff , axis=0)/avg)
                freq =[]
                T=Window*(1/fs)
                t = np.arange(0,T,1/fs)
                freq = np.fft.fftfreq(t.shape[-1],d=1/fs)
                
                f_min = (10**j)-1
                f_max = (10**(j+1))-10**j
                
                mask = (freq >= f_min) & (freq <= f_max)
                
           
                CSD.extend(buff2[mask])
                fr.extend(freq[mask])
                    
            
        else:
            print("averaging not possible inpute len "+ str(len(x))+" divided by the number of window"+str(N)+" excedes avg "+ str(avg))
            X = np.fft.fft(x)
            Y = np.fft.fft(y)
            CSD= 2*(1/fs)*Y*np.conj(X)/(N)
        
        print(fr)
        return fr[0:int(N/2)], np.array(CSD[0:int(N/2)])


# maximum mean functions ##############################################################################


    def Psd_Bertlett_MA(self,x,fs,Max_avg,N = None): 
           
        #Power Specral Density Bertlett with Maximum Average
    
        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        PSD = []
        fr = []
        j=0
        averages=[]
        decades=[]
        
    
        while(True):
            
            Window=int(N/(Max_avg/(10**j)))
            
            buff = []
            buff2 = []
            psd_cross= []
            
            averages.append((Max_avg/(10**j)))
            
            for i in range(int(Max_avg/(10**j))):
    
                X = np.fft.fft(x[i*Window:i*Window+Window])
                Y = X
                app= 2*(1/fs)*Y*np.conj(X)/(Window)
                psd_cross.append(app)
    
            buff = np.array(psd_cross)
            buff2 = (np.sum(buff , axis=0)/(Max_avg/(10**j)))
            freq =[]
            T=Window*(1/fs)
            t = np.arange(0,T,1/fs)
            freq = np.fft.fftfreq(t.shape[-1],d=1/fs)
            
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
    
    
    
    def Csd_MA(self,x,y,fs,Max_avg,N = None): 
           
        #Power Specral Density Bertlett with Maximum Average
    
        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        PSD = []
        fr = []
        j=0
        averages=[]
        decades=[]
        
    
        while(True):
            
            Window=int(N/(Max_avg/(10**j)))
            
            buff = []
            buff2 = []
            psd_cross= []
            
            averages.append((Max_avg/(10**j)))
            
            for i in range(int(Max_avg/(10**j))):
    
                X = np.fft.fft(x[i*Window:i*Window+Window])
                Y = np.fft.fft(y[i*Window:i*Window+Window])
                app= 2*(1/fs)*Y*np.conj(X)/(Window)
                psd_cross.append(app)
    
            buff = np.array(psd_cross)
            buff2 = (np.sum(buff , axis=0)/(Max_avg/(10**j)))
            freq =[]
            T=Window*(1/fs)
            t = np.arange(0,T,1/fs)
            freq = np.fft.fftfreq(t.shape[-1],d=1/fs)
            
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




    def marker_plotter_MA(self, avgerages, decades):
    
        colors = [
        'lightblue', 'lightgreen', 'lightcoral',  'lightpink', 'lavender', 'lightcyan', 'peachpuff', 'mintcream', 'honeydew', 'thistle',
        'powderblue', 'moccasin', 'palegoldenrod', 'lemonchiffon', 'lightgoldenrodyellow', 'mistyrose'
        ]
        
        
        # Loop per tracciare le linee tratteggiate e colorare le porzioni
        for i  in range(len(decades)-1):
            
            # Aggiungere linee tratteggiate
            plt.axvline(x=decades[i-1], color='gray', linestyle='--', lw=1)
            plt.axvline(x=decades[i], color='gray', linestyle='--', lw=1)
            # Colorare la porzione
            plt.axvspan(decades[i],decades[i+1], color=colors[i], alpha=0.5, label=str(int(avgerages[i]))+' average')
 

