
"""
Created on Fri Jun 6 2024

@author: Modestino Carbone

"""


import matplotlib.pyplot as plt
import numpy as np
from scipy import signal
from scipy.interpolate import interp1d

class Spectra(object):

    def __init__(self):
        super(Spectra,self).__init__()

#general purpose fft

    def fft_custom(self,x,fs,N):

        T=N*1/fs
        t = np.arange(0,T,1/fs)
        sp = np.fft.fft(x)
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)
        fft = 2*np.abs(sp)/N

        return freq[0:int(N/2)], fft[0:int(N/2)]


#unoptimized spectrum functions

    def cross_spectrum_avg(self,x,y,fs,N,avg):

        T=N*(1/fs)
        t = np.arange(0,T,1/fs)

        if (len(x)/N >= avg):
    
            csd_cross=[]
            print("performing cross spectrum with "+str(avg)+" means")

            for i in range(avg):

                X = np.fft.fft(x[i*N:i*N+N-1])
                Y = np.fft.fft(y[i*N:i*N+N-1])
                app= 2*(1/fs)*Y*np.conj(X)/(N)
                csd_cross.append(app)

            buff = np.array(csd_cross)
            CSD = (np.sum(buff , axis=0)/avg)
            
        else:
            print("averaging not possible inpute len "+ str(len(x))+" divided by the number of window"+str(N)+" excedes avg "+ str(avg))
            X = np.fft.fft(x)
            Y = np.fft.fft(y)
            CSD= 2*(1/fs)*Y*np.conj(X)/(N)
        
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)

        return freq[0:int(N/2)], CSD[0:int(N/2)]

    def Psd_custom(self,x,fs,N):

        T=N*(1/fs)
        t = np.arange(0,T,1/fs)
        X = np.fft.fft(x)
        CSD= 2*(1/fs)*X*np.conj(X)/(N)
        freq = np.fft.fftfreq(t.shape[-1],d=1/fs)

        return freq[0:int(N/2)], CSD[0:int(N/2)]


#optimized PSD functin
#4 decades fast algorithm

    def PSD_avg_decades(self,x,fs,N,avg):
    
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
                    cpsd_cross.append(app)
        
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
             print("averaging not possible inpute len "+ str(len(x))+" divided by the number of window"+str(N)+" excedes avg "+ str(avg))
            X = np.fft.fft(x)
            Y = X
            PSD= 2*(1/fs)*Y*np.conj(X)/(N)
        
        print(fr)
        return fr[0:int(N/2)], np.array(PSD[0:int(N/2)])


#optimized CSD functin
#4 decades fast algorithm

    def cross_spectrum_avg_decades(self,x,y,fs,N,avg):
    
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

