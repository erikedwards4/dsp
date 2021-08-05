#@author Erik Edwards
#@date 2017-present
#@license BSD 3-clause

#DSP is my own library of C functions for digital signal processing.

#This is the Makefile used during development of the C functions.
#Since these are not stand-alone programs (no main functions), the only outputs are object files to obj dir.
#It is not necessary to run this when making the C++ command-line tools.

SHELL=/bin/bash
CC=clang++

ifeq ($(CC),$(filter $(CC),clang gcc))
	STD=-std=c11
else
	STD=-std=gnu++14
endif

ifeq ($(CC),clang++)
	STD=-std=c++11
	WFLAG=-Weverything -Wno-c++98-compat -Wno-padded -Wno-gnu-imaginary-constant -Wno-deprecated -Wno-old-style-cast
else
	WFLAG=-Wall -Wextra
endif

CFLAGS=$(WFLAG) $(STD) -O3 -march=native -fPIC


All: all
all: Dirs Generate Transform Filter Conv Interp ZCs_LCs AR_Poly AC_LP Frame STFT Wavelets Clean

Dirs:
	mkdir -pm 777 obj


#Generate: generate noise, pulses, waveforms, windows
Generate: Noise Pulses Waves Wins

#Noise: generate colored noise
Noise: white pink red brown blue violet
white: white.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pink: pink.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
red: red.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
brown: brown.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
blue: blue.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
violet: violet.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Pulses: generate pulse functions (unit impulse is a.k.a. dirac delta)
Pulses: unit_impulse delta_impulse rect_pulse
unit_impulse: unit_impulse.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
delta_impulse: delta_impulse.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rect_pulse: rect_pulse.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Waves: generate periodic waveforms
Waves: sinewave cosinewave squarewave triwave sawwave pulsewave
sinewave: sinewave.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
cosinewave: cosinewave.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
squarewave: squarewave.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
triwave: triwave.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sawwave: sawwave.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pulsewave: pulsewave.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Wins: generate common 1-D windows
Wins: rectangular triangular bartlett hann hamming blackman blackmanharris flattop povey gauss tukey planck
rectangular: rectangular.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
triangular: triangular.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
bartlett: bartlett.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hann: hann.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hamming: hamming.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
blackman: blackman.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
blackmanharris: blackmanharris.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
flattop: flattop.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
povey: povey.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
gauss: gauss.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
tukey: tukey.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
planck: planck.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Transform: common 1-D signal transforms
Transform: fft.fftw fft.fftw.r2hc fft.rad2 ifft.fftw ifft.rad2 dct idct dst idst hilbert
fft.fftw: fft.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fft.fftw.r2hc: fft.fftw.r2hc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fft.rad2: fft.rad2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ifft.fftw: ifft.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ifft.rad2: ifft.rad2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dct: dct.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idct: idct.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dst: dst.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idst: idst.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hilbert: hilbert.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Filter: FIR and IIR filters
Filter: #fir iir #fir_tau iir_tau filter filter_tau filtfilt spencer
fir: fir.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iir: iir.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Conv: convolution
Conv: #conv1d conv2d conv1d_fft conv2d_fft
conv1d: conv1d.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv2d: conv2d.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv1d_fft: conv1d_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv2d_fft: conv2d_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Interp: interpolation
Interp: interp1q #interp1t
interp1q: interp1q.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
interp1t: interp1t.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#ZC_LCs: zero-crossings, level-crossings and mean-crossings
ZCs_LCs: #zcs lcs mcs zcr lcr mcr zcr_windowed lcr_windowed mcr_windowed
zcs: zcs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lcs: lcs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mcs: mcs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
zcr: zcr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lcr: lcr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mcr: mcr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
zcr_windowed: zcr_windowed.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
lcr_windowed: lcr_windowed.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
mcr_windowed: mcr_windowed.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#AR_Poly: conversion between AR (autoregressive), poly (polynomial),
#RC (reflection coeff), and PSD (power spectral density) representations.
#Currently, the RC-related functions are not working, which are less important,
#and the Octave tsa functions for RCs are also not working.
AR_Poly: poly2roots roots2poly poly2ar ar2poly poly2psd ar2psd #ar2rc rc2ar poly2rc rc2poly
poly2roots: poly2roots.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
roots2poly: roots2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2ar: poly2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2poly: ar2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2psd: poly2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2psd: ar2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
#ar2rc: ar2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
#rc2ar: rc2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
#poly2rc: poly2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
#rc2poly: rc2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#AC_LP: conversions related to AC (autocovariance) and LP (linear prediction).
#These work from the signal (sig) or autocovariance (ac) to get
#ac, rc (reflection coeff), poly (polynomial params), ar (autoregressive),
#cc (cepstral coeff), or mvdr (minimum-variance distortionless-response) representations.
#Here sig is signal, ac is autocovariance, rc is reflection coefficients,
#poly is polynomial params, and ar is autoregressive params.
AC_LP: sig2ac sig2ac_fft ac2rc_levdurb ac2ar_levdurb ac2poly_levdurb sig2rc_levdurb sig2ar_levdurb sig2poly_levdurb sig2ar_burg sig2poly_burg #ac2rc ac2cc ac2mvdr
sig2ac: sig2ac.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ac_fft: sig2ac_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2rc_levdurb: ac2rc_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2ar_levdurb: ac2ar_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2poly_levdurb: ac2poly_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2rc_levdurb: sig2rc_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ar_levdurb: sig2ar_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2poly_levdurb: sig2poly_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ar_burg: sig2ar_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2poly_burg: sig2poly_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

ac2cc_levdurb: ac2cc_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2cc_levdurb: sig2cc_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2mvdr: ac2mvdr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2mvdr: sig2mvdr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Frame: get frames and apply windows for univariate signal to put into matrix
Frame: frame_univar frame_univar_flt apply_win window_univar window_univar_flt
frame_univar: frame_univar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
frame_univar_flt: frame_univar_flt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
apply_win: apply_win.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
window_univar: window_univar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
window_univar_flt: window_univar_flt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#STFT: steps to do the STFT (short-term Fourier transform)
STFT: fft_power get_stft_freqs stft stft_flt
fft_power: fft_power.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
get_stft_freqs: get_stft_freqs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
stft: stft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
stft_flt: stft_flt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Spectrogram: steps to make a spectrogram (frequency-reweighted STFT)
Spectrogram: convert_freqs
convert_freqs: convert_freqs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Wavelets: a couple of my most often-used wavelets
Wavelets: #gabor analytic
gabor: gabor.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
analytic: analytic.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#If want to make a library archive for C++ use, this works with clang++
#ar crs ../lib/libdsp.a obj/*.o


#make clean
Clean: clean
clean:
	find ./obj -type f -name *.o | xargs rm -f
	rm -f *.s tmp.*
