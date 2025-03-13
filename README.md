# Spectra Library

A Python library for custom FFT, cross-spectrum, and PSD calculations, with support for optimized algorithms across multiple decades.
In some signal processing computation specific spectrum tools can be really useful like a multipurpose Swiss Army knife.

## Features
- **FFT**: Compute the FFT Fast Fourier transform of a signal.
- **PSD**: Compute the PSD Power spectral density of a signal.
- **ASD**:  Compute the ASD Amplitude Spectral density  of a signal.
- **CSD**:  Compute the CSD Cross-Spectral density  between two signals.
- **Coherence**:  Compute the Coherence between two signals.
- **Fixed average**: PSD and CSD  is computed with a smaller frequency resolution for high frequency.
- **Maximum average**: PSD and CSD  is computed with a smaller frequency resolution for higher decades and an discend average for every lower slot. 

## Installation

Clone the repository:
```bash
git clone https://github.com/lunaticmercury/Spectra_lib.git
```
Install dependencies using pip:
```bash
pip install -r requirements.txt
```

## Usage
clone the python file Spectra_lib.py in the desired workspace and put the following lines at the head of the main.py file

```python
from Spectra_lib import Spectra 
spec = Spectra()
```
