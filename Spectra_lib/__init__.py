

"""
Created on Fri March 21 2025

@author: Modestino Carbone


"""


from .fft import fft
from .psd import psd_bartlett, psd_welch, psd_MA, psd_Min_Max
from .csd import csd_basic, csd_MA, csd_Min_Max
from .graph import marker_plotter_MA


class Spectra:

    def __init__(self):
        print("Spectra functions imported")
        
    def Fft(self,x, fs, window , N):
        """
        Description:\n
        This function computes the Fast Fourier Transform (FFT) of a given signal `x`, applying a specified window function (`window`) and normalizing the result.\n
        
        Parameters:\n
        - x (numpy array): The input signal to be transformed. It can be a 1D array or list.\n
        - fs (float): The sampling frequency of the signal `x` (in Hz).\n
        - window (str): The type of window to apply to the signal before computing the FFT (e.g., `"hann"`, `"hamming"`, etc.).\n
        - N (int): The number of points to use for the FFT. If not provided, the function will use the length of `x` as `N`.\n
        
        Returns:\n
        - freq (numpy array): A 1D array containing the frequency values corresponding to the FFT bins, ranging from 0 to `fs/2`.\n
        - FFT (numpy array): A 1D array containing the magnitude of the FFT for each frequency, normalized according to the window function used.\n
        
        """
        return fft(x, fs, window,N )
    
    def Psd_Bartlett(self,x, fs, avg, window, N):
        
        """
        Description:\n
        This function computes the **Power Spectral Density (PSD)** of a signal `x` using the **Bartlett's method** (also known as the **modified periodogram**). It applies a specified window function (`window`) and averages the results over multiple segments. This method is commonly used for estimating the PSD of a signal, especially when the signal is not stationary over time.
        The Bartlett method splits the signal into overlapping segments, applies the FFT on each segment, and averages the individual periodograms to reduce variance. If the length of `x` is larger than the segment length multiplied by `avg`, it will split the signal into segments, otherwise, it will compute the PSD on the entire signal.\n

        Parameters:\n
        -  x (numpy array): The input signal to be analyzed. It can be a 1D array or list containing the signal data.\n
        -  fs (float): The sampling frequency of the signal `x` (in Hz).\n
        -  avg (int): The number of segments over which to average the PSD. This helps in reducing variance in the estimation.\n
        -  window (str): The type of window to apply to the signal before computing the FFT (e.g., `"hann"`, `"hamming"`, etc.).\n
        -  N (int): The number of points to use for the FFT in each segment. If not provided, the function will use the length of `x` as `N`.\n

        Returns:\n
        - freq (numpy array): A 1D array containing the frequency values corresponding to the PSD bins, ranging from 0 to `fs/2`.\n
        - PSD (numpy array): A 1D array containing the estimated Power Spectral Density of the signal for each frequency.\n
        """
        
        return psd_bartlett(x, fs, avg, window, N )
    
    def Psd_Welch(self, x, fs, avg, window, N):
        
        """
        Description:\n
        This function computes the **Power Spectral Density (PSD)** of a signal `x` using the **Welch's method**. It applies a specified window function (`window`) and averages the results over multiple and overlapped segments. This method is commonly used for estimating the PSD of a signal, especially when the signal is stationary over time.
        If the length of `x` is larger than the segment length multiplied by `avg`, it will split the signal into segments, otherwise, it will compute the PSD on the entire signal.\n

        Parameters:\n
        - x (numpy array): The input signal to be analyzed. It can be a 1D array or list containing the signal data.\n
        - fs (float): The sampling frequency of the signal `x` (in Hz).\n
        - avg (int): The number of segments over which to average the PSD. This helps in reducing variance in the estimation.\n
        - window (str): The type of window to apply to the signal before computing the FFT (e.g., `"hann"`, `"hamming"`, etc.).\n
        - N (int): The number of points to use for the FFT in each segment. If not provided, the function will use the length of `x` as `N`.\n

        Returns:\n
        - freq (numpy array): A 1D array containing the frequency values corresponding to the PSD bins, ranging from 0 to `fs/2`.\n
        - PSD (numpy array): A 1D array containing the estimated Power Spectral Density of the signal for each frequency.\n
        """
        
        return psd_welch(x, fs, avg, window, N)
    
    def Psd_MA(self, x, fs, Max_avg, window, N):
        
        """
        Description:\n
        This function computes the Power Spectral Density (PSD) using Bartlett's method fixing the Maximum Averaging. \n

        Parameters:\n
        - x (numpy array): Input signal .\n
        - fs (float): Sampling frequency of the signal.\n
        - Max_avg (int): Maximum number of averages for PSD estimation.\n
        - window (str): Window type for FFT processing.\n
        - N (int) : Number of points for the Fourier Transform.\n

        Returns:\n
        - freq (numpy array): Frequency vector corresponding to the PSD.\n
        - PSD (numpy array): Estimated Power Spectral Density.\n
        - averages (list of int): Number of segments used for averaging in each frequency range.\n
        - decades (list of int): Frequency range divisions in logarithmic scale.\n
        """
        
        return psd_MA(x, fs, Max_avg, N, window)
    
    def Psd_Min_Max(self, x, fs, Min_avg, Max_avg,  window, N): 
        
        """
        Description:\n
        This function computes the Power Spectral Density (PSD) using Bartlett's method fixing the Maximum Averaging and the Minimum one. \n

        Parameters:\n
        - x (numpy array): Input signal .\n
        - fs (float): Sampling frequency of the signal.\n
        - Max_avg (int): Maximum number of averages for PSD estimation.\n
        - Min_avg (int): Maximum number of averages for PSD estimation.\n
        - window (str): Window type for FFT processing.\n
        - N (int) : Number of points for the Fourier Transform.\n

        Returns:\n
        - freq (numpy array): Frequency vector corresponding to the PSD.\n
        - PSD (numpy array): Estimated Power Spectral Density.\n
        - averages (list of int): Number of segments used for averaging in each frequency range.\n
        - decades (list of int): Frequency range divisions in logarithmic scale.\n
        """
        
        return psd_Min_Max(x, fs, Min_avg, Max_avg, N, window)
    
    
    def Csd(self, x, y, fs, avg, window, N):
        
        """
        Description:\n
        This function computes the Cross-Spectral Density (CSD) between two signals.\n
    
        Parameters:\n
        - x (numpy array):  The first input signal to be analyzed. It can be a 1D array or list containing the signal data.\n
        - y (numpy array): The second input signal to be analyzed. It can be a 1D array or list containing the signal data.\n
        - fs (float): Sampling frequency of the signals.\n
        - avg (int): Number of averages to improve spectral estimation.\n
        - window (str): Window type to be applied before FFT.\n
        - N (int): Number of points for the Fourier Transform.\n
    
        Returns:\n
        - freq (numpy array): Frequency vector corresponding to the CSD.\n
        - CSD (numpy array): Estimated Cross-Spectral Density.\n
        """
        
        return csd_basic(x, y, fs, avg, window, N)
    
    def Csd_MA(self, x, y, fs, Max_avg, window, N): 
        
        """
        Description:\n
        This function computes the Cross-Spectral Density (CSD) fixing the Maximum Average. \n
    
        Parameters:\n
        - x (numpy array): First input signal.\n
        - y (numpy array): Second input signal.\n
        - fs (float): Sampling frequency.\n
        - Max_avg (int): Maximum number of averages.\n
        - window (str): Window type for FFT processing.\n
        - N (int): Number of points for the Fourier Transform.\n
    
        Returns:\n
        - freq (numpy array): Frequency vector.\n
        - CSD (numpy array): Computed Cross-Spectral Density.\n
        - averages (list of int): List of averages used per frequency range.\n
        - decades (list of int): Frequency range divisions.\n
        """
        
        return csd_MA(x, y, fs, Max_avg, N, window)
    
    def Csd_Min_Max(self, x, y, fs, Min_avg, Max_avg, window, N):
        
        """
        Description:\n
        This function computes the Cross-Spectral Density (CSD) fixing the Maximum Average and the Minimum one.\n
    
        Parameters:\n
        - x (numpy array): First input signal.\n
        - y (numpy array): Second input signal.\n
        - fs (float): Sampling frequency.\n
        - Max_avg (int): Maximum number of averages.\n
        - Min_avg (int): Minimum number of averages.\n
        - window (str): Window type for FFT processing.\n
        - N (int): Number of points for the Fourier Transform.\n
    
        Returns:\n
        - freq (numpy array): Frequency vector.\n
        - CSD (numpy array): Computed Cross-Spectral Density.\n
        - averages (list of int): List of averages used per frequency range.\n
        - decades (list of int): Frequency range divisions.\n
        """
        
        return  csd_Min_Max(x, y, fs, Min_avg, Max_avg, N, window)
    
    def Avg_dec_col(self, averages, decades, label):
        
        """
        Description:\n
        This function computes plots frequency ranges with shaded regions and dashed vertical lines based on averaging parameters.\n

        Parameters:\n
        - averages (list of int): List of averaging values used for each frequency range.\n
        - decades (list of int): List of frequency range divisions in logarithmic scale.\n
        - label (boolean): Boolean indicating whether to add labels to the shaded regions.\n

        The function highlights different frequency ranges using different colors and marks decade boundaries with dashed lines.\n
        """

        return marker_plotter_MA(averages, decades, label);