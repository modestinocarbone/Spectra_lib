# Spectra Library

A Python library for custom FFT, cross-spectrum, and PSD calculations, with support for optimized algorithms across multiple decades.
In some signal processing computation specific spectrum tools can be really useful like a multipurpose Swiss Army knife.

## Features
- **FFT**: Compute the Fast Fourier Transform (FFT) of a signal.
- **CSD**: Cross-spectral density estimation  between two signals with support for averaging.
- **PSD**: Power Spectral Density estimation for a signal.

## Installation

Clone the repository:
```bash
git clone https://github.com/modestinocarbone/Spectra_lib
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

The file example.py is a perfect proof of concept of various kind of function implementations offered by the library
