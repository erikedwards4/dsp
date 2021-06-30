# dsp

dsp: functions for DSP (digital signal processing)

================================================

This is a set of C functions, and associated command-line tools in C++,  
that implement many of the usual DSP functions found in Octave, etc.  
The purpose is to support speech and neural work.  
This does not, therefore, have every possible DSP function.  
In particular, only 1D signal processing is supported,  
including matrices and tensors of multivariate signals,  
but not 2D (image) processing.  

The command-line programs are written in C++ with a consistent style and interface.  
The low-level functions themselves are written in C for fastest performance (e.g., openBLAS).  

The C functions are meant for the developer; the C++ command-line tools are meant for the end-user.  
The interface to each C function is BLAS-like, meaning that one specifies the input and/or output dimensions,  
the matrix order as row-major or column-major, and so on.  

The C++ command-line programs are written in a consistent style that was developed for command-line tools in general.  
All of these command-line tools use argtable2 (http://argtable.sourceforge.net/) for parsing inputs and option flags.  
All of them allow -h (--help) as a flag to give description and usage info.  

Input/output is supported for NumPy tensors (https://numpy.org/)  
and several C++ tensor formats: Armadillo (http://arma.sourceforge.net/),  
ArrayFire (https://arrayfire.com/), and a minimal format for Eigen (http://eigen.tuxfamily.org/).  


## Dependencies
Requires argtable2, openBLAS, LAPACKE, FFTW.  
For Ubuntu, these are available by apt-get:  
```
sudo apt-get install libargtable2-0 libblas3 libopenblas-base liblapack3 liblapacke fftw3  
```

You must first install the util library:  
https://github.com/erikedwards4/util  
And install dsp into the same parent directory as util.  
Preferably: /opt/codee/util and /opt/codee/dsp  
For full examples and support functions, also install math:  
https://github.com/erikedwards4/math  


## Installation
```
cd /opt/codee  
git clone https://github.com/erikedwards4/dsp  
cd /opt/codee/dsp  
make  
```

Each C function can also be compiled and used separately; see c subdirectory Makefile for details.  


## Usage
See each resulting command-line tool for help (use -h or --help option).  
For example:  
```
/opt/codee/dsp/bin/fir --help
```


## List of functions
All: Generate Interp Filter Conv Transform AR_Poly AC_LP  
Generate: Noise Waves Wins  
Noise: white pink  
Waves: sinewave squarewave triwave sawwave pulsewave  
Wins: rectangular triangular bartlett hann hamming blackman blackmanharris flattop povey gauss tukey planck  
Interp: interp1q  
Transform: fft ifft dct idct dst idst hilbert  
Filter: fir iir filter filtfilt fftfilt medfilt integrate spencer  
AR_Poly: poly2roots roots2poly poly2ar ar2poly ar2rc rc2ar poly2rc rc2poly ar2psd poly2psd  
AC_LP: sig2ac sig2ac_fft ac2ar_levdurb ac2poly_levdurb sig2poly_levdurb sig2ar_levdurb sig2ar_burg sig2poly_burg ac2rc ac2cc ac2mvdr  
ZCs_LCs: zcs lcs mcs zcr lcr mcr zcr_windowed lcr_windowed mcr_windowed  
Freqs: convert_freqs get_cfs get_stft_freqs get_cfs_T  
Frame: frame_univar apply_win window_univar  
STFT: fft_hc hc_square fft_squared stft  
Wavelets: gabor analytic  


## Contributing
This is currently only to view the project in progress.


## License
[BSD 3-Clause](https://choosealicense.com/licenses/bsd-3-clause/)
