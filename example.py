import numpy as np
import matplotlib.pyplot as plt
from Spectra_lib import Spectra
spec = Spectra()

fs= 100000
N = 100000
T = N/fs

# nyquist freq = fs/2 = 5000 Hz
# freq of the signal at least nyquist_freq/2.5 => freq_sig <= 2500 Hz
# at least 10 sample per period are needed to sample a real signal to avoid distortion

freq = 2000

t = np.linspace(0,T,N)
noise1 = np.random.uniform(-0.01, 0.01, N)
noise2 = np.random.uniform(-0.01, 0.01, N)
x = np.sin(2*np.pi*freq*t)+noise1
y = np.sin(2*np.pi*freq*t)+noise2

#%%

plt.figure(1)

plt.grid()
plt.plot(t,x)
plt.xlabel("time /s")
plt.show()

#%%

plt.figure(2)
plt.title("FFT")
f, ff = spec.Fft(x, fs, "hann", N)
#logarithmic grid
plt.grid(True, which="both", ls="--", color='gray') 
#fft plot
plt.semilogx(f,ff)
plt.xlabel("Frequency /Hz")
plt.legend()
plt.show()

#%%
  
plt.figure(3)
plt.title("PSD")

f, ps = spec.Psd_Bartlett(x, fs, avg = 1, window = "boxcar", N = int(N/1))
plt.semilogx(f,10*np.log10(abs(ps)),label='Bertlett Psd')

f, ps = spec.Psd_Welch(x, fs, avg = 1, window = "boxcar", N = int(N/1))
plt.semilogx(f,10*np.log10(abs(ps)),label='Welch Psd')

#peak encountered at -> Vrms^2/bin; bin = 1Hz 

plt.grid(True, which="both", ls="--", color='gray')
plt.ylabel(r"dB , $V^2$ 1/Hz")
plt.xlabel("Frequency /Hz")
plt.xlim(1,max(f))
plt.legend()
plt.show()

#%%
  
plt.figure(4)
plt.title("PSD Maximum average")

f, ps, avg, dec = spec.Psd_MA(x, fs, 1000, window = "boxcar", N = int(N/1))
plt.semilogx(f,10*np.log10(abs(ps)),label='PSD')

spec.Avg_dec_col(avg, dec, label = True)

plt.grid(True, which="both", ls="--", color='gray')
plt.ylabel("dB , $V^2$ 1/Hz")
plt.xlabel("Frequency /Hz")
plt.xlim(1,max(f))
plt.legend()
plt.show()

#%%
  
plt.figure(5)
plt.title("PSD Min-Max average")


f, ps, avg, dec = spec.Psd_Min_Max(x,fs, 10, 1000, window = "boxcar", N = int(N/1))
plt.semilogx(f,10*np.log10(abs(ps)),label='PSD')

spec.Avg_dec_col(avg, dec,  label = True)

plt.grid(True, which="both", ls="--", color='gray')
plt.ylabel("dB , $V^2$ 1/Hz")
plt.xlabel("Frequency /Hz")
plt.xlim(1,max(f))
plt.legend()
plt.show()

#%%

plt.figure(6)
plt.title("CSD")

f, ps = spec.Csd(x, y, fs, avg = 1, window = "boxcar", N = int(N/1))
plt.semilogx(f,10*np.log10(abs(ps)),label='Csd')

plt.grid(True, which="both", ls="--", color='gray')
plt.ylabel("dB , $V^2$ 1/Hz")
plt.xlabel("Frequency /Hz")
plt.xlim(1,max(f))
plt.legend()
plt.show()

#%%

plt.figure(7)
plt.title("CSD Maximum average")

f, ps, avg, dec = spec.Csd_MA(x, y, fs, 1000, window = "boxcar", N = int(N/1))
plt.semilogx(f,10*np.log10(abs(ps)),label='Csd')

spec.Avg_dec_col(avg, dec,  label = True)

plt.grid(True, which="both", ls="--", color='gray')
plt.ylabel("dB , $V^2$ 1/Hz")
plt.xlabel("Frequency /Hz")
plt.xlim(1,max(f))
plt.legend()
plt.show()

#%%

plt.figure(8)
plt.title("CSD Min-Max average")

f, ps, avg, dec = spec.Csd_Min_Max(x, y, fs, 10, 1000, window = "boxcar", N = int(N/1))
plt.semilogx(f,10*np.log10(abs(ps)),label='Csd')

spec.Avg_dec_col(avg, dec, label = True)

plt.grid(True, which="both", ls="--", color='gray')
plt.ylabel("dB , $V^2$ 1/Hz")
plt.xlabel("Frequency /Hz")
plt.xlim(1,max(f))
plt.legend()
plt.show()