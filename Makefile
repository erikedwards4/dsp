#@author Erik Edwards
#@date 2018-present
#@license BSD 3-clause

#dsp is my own library of functions for digital signal processing (DSP) in C and C++.
#This is the makefile for the C++ command-line tools.

#This is meant primarily to support work in audio, speech, NNs, and some general time-series analyses.
#Thus, not every possible DSP function is included.
#Also, 2D DSP (which is essentially image processing) is not included.
#This does support matrices and tensors up to 4D, i.e. for multivariate 1D signals.

SHELL=/bin/bash
ss=../util/bin/srci2src
CC=clang++

ifeq ($(CC),clang++)
	STD=-std=c++11
	WFLAG=-Weverything -Wno-c++98-compat -Wno-old-style-cast -Wno-gnu-imaginary-constant
else
	STD=-std=gnu++14
	WFLAG=-Wall -Wextra
endif

INCLS=-Ic -I../util
CFLAGS=$(WFLAG) $(STD) -O2 -ffast-math -march=native $(INCLS)


All: all
all: Dirs Generate Wins Transform Conv Filter Interp ZCs_LCs AR_Poly AC_LP Freqs Frame STFT Clean
	rm -f 7 obj/*.o

Dirs:
	mkdir -pm 777 bin obj


#Generate: generate noise, simple waveforms, windows
Generate: Noise Waves Wins

#Noise: generate colored noise
Noise: #white pink brown
white: srci/white.cpp c/white.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-c++98-compat-pedantic; $(CC) obj/$@.o -obin/$@ -largtable2
pink: srci/pink.cpp c/pink.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-c++98-compat-pedantic; $(CC) obj/$@.o -obin/$@ -largtable2
brown: srci/brown.cpp c/brown.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-c++98-compat-pedantic; $(CC) obj/$@.o -obin/$@ -largtable2

#Waves: generate periodic waveforms
Waves: #sinewave squarewave triwave sawwave pulsewave
sinewave: srci/sinewave.cpp c/sinewave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
squarewave: srci/squarewave.cpp c/squarewave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
triwave: srci/triwave.cpp c/triwave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
sawwave: srci/sawwave.cpp c/sawwave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
pulsewave: srci/pulsewave.cpp c/pulsewave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Wins: generate common 1-D windows
Wins: rectangular triangular bartlett hann hamming blackman blackmanharris flattop povey gauss tukey planck
rectangular: srci/rectangular.cpp c/rectangular.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
triangular: srci/triangular.cpp c/triangular.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
bartlett: srci/bartlett.cpp c/bartlett.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
hann: srci/hann.cpp c/hann.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
hamming: srci/hamming.cpp c/hamming.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
blackman: srci/blackman.cpp c/blackman.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
blackmanharris: srci/blackmanharris.cpp c/blackmanharris.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
flattop: srci/flattop.cpp c/flattop.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
povey: srci/povey.cpp c/povey.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
gauss: srci/gauss.cpp c/gauss.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
tukey: srci/tukey.cpp c/tukey.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
planck: srci/planck.cpp c/planck.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm


#Transform: common 1-D signal transforms
Transform: fft ifft dct idct dst idst hilbert
fft: srci/fft.cpp c/fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3
ifft: srci/ifft.cpp c/ifft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3
dct: srci/dct.cpp c/dct.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
idct: srci/idct.cpp c/idct.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
dst: srci/dst.cpp c/dst.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3
idst: srci/idst.cpp c/idst.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3
hilbert: srci/hilbert.cpp c/hilbert.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3


#Filter: FIR and IIR filters
Filter: fir iir #filter filtfilt fftfilt medfilt integrate spencer
integrate: srci/integrate.cpp c/integrate.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
fir: srci/fir.cpp c/fir.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
iir: srci/iir.cpp c/iir.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm


#Conv: 1-D convolution
Conv: #conv1 conv1_fft
conv1: srci/conv1.cpp c/conv1.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
conv1_fft: srci/conv1_fft.cpp c/conv1_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm



#Interp: 1-D interpolation
Interp: #interp1q interp1t
interp1q: srci/interp1q.cpp c/interp1q.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
interp1t: srci/interp1t.cpp c/interp1t.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2


#ZC_LCs: zero-crossings, level-crossings and mean-crossings
ZCs_LCs: zcs lcs mcs zcr lcr mcr zcr_windowed lcr_windowed mcr_windowed
zcs: srci/zcs.cpp c/zcs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
lcs: srci/lcs.cpp c/lcs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
mcs: srci/mcs.cpp c/mcs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
zcr: srci/zcr.cpp c/zcr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
lcr: srci/lcr.cpp c/lcr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
mcr: srci/mcr.cpp c/mcr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
zcr_windowed: srci/zcr_windowed.cpp c/zcr_windowed.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
lcr_windowed: srci/lcr_windowed.cpp c/lcr_windowed.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
mcr_windowed: srci/mcr_windowed.cpp c/mcr_windowed.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm


#AR_Poly: conversion between AR (autoregressive), poly (polynomial),
#RC (reflection coeff), and PSD (power spectral density) representations
AR_Poly: poly2roots roots2poly poly2ar ar2poly ar2rc rc2ar poly2rc rc2poly ar2psd poly2psd
poly2roots: srci/poly2roots.cpp c/poly2roots.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -llapacke -lm
roots2poly: srci/roots2poly.cpp c/roots2poly.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
poly2ar: srci/poly2ar.cpp c/poly2ar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
ar2poly: srci/ar2poly.cpp c/ar2poly.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
ar2rc: srci/ar2rc.cpp c/ar2rc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
rc2ar: srci/rc2ar.cpp c/rc2ar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
poly2rc: srci/poly2rc.cpp c/poly2rc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
rc2poly: srci/rc2poly.cpp c/rc2poly.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
ar2psd: srci/ar2psd.cpp c/ar2psd.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
poly2psd: srci/poly2psd.cpp c/poly2psd.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm


#AC_LP: conversions between sig (signal), AC (autocorrelation), LP (linear prediction), and related.
#For example, ac2rc converts from AC to RCs (reflection coeffs).
AC_LP: autocorr autocorr_fft sig2ac sig2ac_fft ac2ar_levdurb ac2poly_levdurb sig2poly_levdurb sig2ar_levdurb sig2ar_burg sig2poly_burg ac2rc ac2cc ac2mvdr
sig2ac: srci/sig2ac.cpp c/sig2ac.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
sig2ac_fft: srci/sig2ac_fft.cpp c/sig2ac_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lfftw3f -lfftw3 -lm
ac2ar_levdurb: srci/ac2ar_levdurb.cpp c/ac2ar_levdurb.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
ac2poly_levdurb: srci/ac2poly_levdurb.cpp c/ac2poly_levdurb.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
sig2ar_levdurb: srci/sig2ar_levdurb.cpp c/sig2ar_levdurb.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lfftw3f -lfftw3 -lm
sig2poly_levdurb: srci/sig2poly_levdurb.cpp c/sig2poly_levdurb.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lfftw3f -lfftw3 -lm
sig2ar_burg: srci/sig2ar_burg.cpp c/sig2ar_burg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
sig2poly_burg: srci/sig2poly_burg.cpp c/sig2poly_burg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
ac2rc: srci/ac2rc.cpp c/ac2rc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
ac2cc: srci/ac2cc.cpp c/ac2cc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
ac2mvdr: srci/ac2mvdr.cpp c/ac2mvdr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lfftw3f -lfftw3 -lm


#Frame: get frames and apply windows for univariate signal to put into matrix
Frame: frame_univar apply_win window_univar
frame_univar: srci/frame_univar.cpp c/frame_univar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
apply_win: srci/apply_win.cpp c/apply_win.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas
window_univar: srci/window_univar.cpp c/window_univar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm


#Freqs: gets center frequencies or STFT frequencies, and converts btwn frequency scales
Freqs: convert_freqs get_cfs get_stft_freqs get_cfs_T #get_cns
convert_freqs: srci/convert_freqs.cpp c/convert_freqs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
get_cfs: srci/get_cfs.cpp c/get_cfs.c c/convert_freqs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
get_stft_freqs: srci/get_stft_freqs.cpp c/get_stft_freqs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
get_cfs_T: srci/get_cfs_T.cpp c/get_cfs_T.c c/convert_freqs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm


#STFT: steps to do the STFT (short-term Fourier transform)
STFT: fft_hc hc_square fft_squared stft
fft_hc: srci/fft_hc.cpp c/fft_hc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lfftw3f -lfftw3 -lm
hc_square: srci/hc_square.cpp c/hc_square.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
fft_squared: srci/fft_squared.cpp c/fft_squared.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lfftw3f -lfftw3 -lm
stft: srci/stft.cpp c/stft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lfftw3f -lfftw3 -lm


#Wavelets: a couple of my most often-used wavelets
Wavelets: #gabor analytic
gabor: srci/gabor.cpp c/gabor.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
analytic: srci/analytic.cpp c/analytic.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm



#make clean
Clean: clean
clean:
	find ./obj -type f -name *.o | xargs rm -f
	rm -f 7 X* Y* x* y*
