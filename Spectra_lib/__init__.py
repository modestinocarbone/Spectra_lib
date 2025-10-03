"""
 *  File:        _init_.py
 *  Author:      Modestino Carbone
 *  Created on:  21/03/2025
 *  Last update: 01/10/2025
 *
 *  Description:
 *      Top hierarchy; this file contains defintions and guide to the use of library 
 *      All function are described in the docstring format
 *
"""

from .fft import fft
from .psd import psd_bartlett, psd_welch, psd_Min_Max, psd_composite
from .csd import csd_bartlett, csd_welch, csd_Min_Max, csd_composite
from .graph import marker_plotter_MA
from .save_load import save, load, savewstyle, loadwstyle

class Spectra:

    def __init__(self):
        print("Spectra functions imported")

# FFT ------------------------------------------      
    def Fft(self, x, fs, window , N):
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

# PSD ------------------------------------------ 
    def Psd(self, x, fs, avg, window, N):
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

    def Psd_Min_Max(self, x, fs, Min_avg, Max_avg,  window):
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
        return psd_Min_Max(x, fs, Min_avg, Max_avg, window)
    
    def psd_composite(self, x, fs, avg, win_type, type):
        """
        Computes the composite Power Spectral Density (CSD) of the 'x' input array.\n

        It computes the PSD using the specified method ("bartlett" or "welch")\n
        and aggregates only the frequency components within a specific frequency range\n
        based on the Nyquist frequency.\n

        Parameters:\n
            x (list of array-like): List of input signals.\n
            fs (list of float): Sampling frequencies for each pair of signals.\n
            avg (int): Number of segments to divide each signal into (for averaging).\n
            win_type (str): Type of window to apply to each segment (e.g., 'hann', 'hamming').\n
            type (str): Spectral estimation method; must be either "bartlett" or "welch".\n

        Returns:\n
            tuple:\n
                fr (list): Concatenated frequency values after masking.\n
                PSD (list): Corresponding composite PSD values.\n
        Notes:\n
            - The resolution of the spectral estimate is fs / W, where W = len(signal) / avg.\n
            - Frequencies are masked to retain only those within a logarithmic decade below and up to the Nyquist frequency.\n
            - Make sure that `csd_bartlett` and `csd_welch` functions are implemented and return correct outputs.\n
        """
        return psd_composite(x, fs, avg, win_type, type)
    
# CSD ------------------------------------------    
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
        return csd_bartlett(x, y, fs, avg, window, N)
    
    def Csd_Welch(self, x, y, fs, avg, window, N):
        """
        Description:\n
        This function computes the Cross-Spectral Density (CSD) between two signals using the **Welch's method**.\n
    
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
        
        return csd_welch(x, y, fs, avg, window, N)
 
    def Csd_Min_Max(self, x, y, fs, Min_avg, Max_avg, window):
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
        return  csd_Min_Max(x, y, fs, Min_avg, Max_avg, window)
    
    def csd_composite(self, x, y, fs, avg, win_type, type):
        """
        Computes the composite Cross Spectral Density (CSD) between corresponding elements of two input arrays.\n

        For each pair of signals in the input lists, it computes the CSD using the specified method\n
        ("bartlett" or "welch") and aggregates only the frequency components within a specific frequency range\n
        based on the Nyquist frequency.\n

        Parameters:\n
            x (list of array-like): List of input signals (first set).\n
            y (list of array-like): List of input signals (second set), must match `x` in length and shape.\n
            fs (list of float): Sampling frequencies for each pair of signals.\n
            avg (int): Number of segments to divide each signal into (for averaging).\n
            win_type (str): Type of window to apply to each segment (e.g., 'hann', 'hamming').\n
            type (str): Spectral estimation method; must be either "bartlett" or "welch".\n

        Returns:\n
            tuple:\n
                fr (list): Concatenated frequency values after masking.\n
                CSD (list): Corresponding composite CSD values.\n
        Notes:\n
            - The resolution of the spectral estimate is fs / W, where W = len(signal) / avg.\n
            - Frequencies are masked to retain only those within a logarithmic decade below and up to the Nyquist frequency.\n
            - Make sure that `csd_bartlett` and `csd_welch` functions are implemented and return correct outputs.\n
        """
        return csd_composite(x, y, fs, avg, win_type, type)
# LOAD & STORE OPERATION   

    def Save(self, f, spec, file_name):
        """
        Saves four arrays (frequency, spectrum, average, and decimation) into a tab-delimited text file.\n

        Parameters:\n
            f (array-like): Frequency data.\n
            spec (array-like): Spectrum data.\n
            avg (array-like): Averaged data.\n
            dec (array-like): Decimated data.\n
            file_name (str): Path to the output file.\n

        Notes:\n
            If the file already exists, it will be overwritten.\n
        """
        save(f, spec, file_name)

    def Savewstyle(self, f, spec, avg, dec, file_name):
        """
        Saves four arrays (frequency, spectrum, averages, and decades) into a tab-delimited text file.\n

        Parameters:\n
            f (array-like): Frequency data.\n
            spec (array-like): Spectrum data.\n
            avg (array-like): Average array.\n
            dec (array-like): decades array.\n
            file_name (str): Path to the output file.\n

        Notes:\n
            If the file already exists, it will be overwritten.\n
        """
        savewstyle(f, spec, avg, dec, file_name)

    def Load(self, file_name):
        """
        Loads frequency and spectrum data from a tab-delimited text file.\n

        Parameters:\n
            file_name (str): Path to the file to load.\n

        Returns:\n
            tuple: Two arrays (f, spec) representing frequency and spectrum.\n

        Raises:\n
            Exception: If the file does not exist.\n
        """
        return load(file_name)

    def Loadwstyle(self, file_name):
        """
        Loads frequency, spectrum, average, and decimation data from a tab-delimited text file.\n

        Parameters:\n
            file_name (str): Path to the file to load.\n

        Returns:\n
            tuple: Four arrays (f, spec, avg, dec) representing frequency, spectrum,\n
                average and decades arrays.\n

        Raises:\n
            Exception: If the file does not exist.\n
        
        Note:\n
            This function currently returns incorrect data for avg and dec.\n
            It assigns them to the same columns as f and spec.\n
            Ensure the file contains four columns and correct indexing is used.\n
        """
        return loadwstyle(file_name)


# GRAPHICS  ------------------------------------ 
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
