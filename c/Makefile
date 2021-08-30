#@author Erik Edwards
#@date 2017-present
#@license BSD 3-clause

#DSP is my own library of C functions for digital signal processing.

#This is the Makefile used during development of the C functions.
#Since these are not stand-alone programs (no main functions), the only outputs are object files to obj dir.
#It is not necessary to run this when making the C++ command-line tools.

SHELL=/bin/bash
CC=clang

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

CFLAGS=$(WFLAG) $(STD) -O3 -ffast-math -march=native -mfpmath=sse -fPIC


All: all
all: Dirs Generate Transform Filter Conv Interp ZCs_LCs AR_Poly Linear_Pred Frame STFT Wavelets Nonlinear Clean

Dirs:
	mkdir -pm 777 obj


#Generate: generate noise, pulses, waveforms, windows
Generate: Noise Pulses Waves Wins Filts

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

#Filts: generate some common/useful filter coefficients (vectors A and/or B to use with iir, fir, filter, etc.)
Filts: spencer #rc butter smooth_diff smooth_diffdiff
spencer: spencer.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Transform: common 1-D signal transforms
Transform: FFT DFT DCT DST Hilbert

#FFT: fast Fourier transforms
FFT: fft.rad2 ifft.rad2 fft.fftw ifft.fftw fft.fftw.r2hc fft.ffts ifft.ffts fft.kiss ifft.kiss
fft.rad2: fft.rad2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ifft.rad2: ifft.rad2.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fft.fftw: fft.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ifft.fftw: ifft.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fft.fftw.r2hc: fft.fftw.r2hc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fft.ffts: fft.ffts.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ifft.ffts: ifft.ffts.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fft.kiss: fft.kiss.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ifft.kiss: ifft.kiss.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#DFT: discrete Fourier transforms (same output as FFT, but uses matrix multiplies)
DFT: dft idft dft.cblas idft.cblas
dft: dft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idft: idft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dft.cblas: dft.cblas.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idft.cblas: idft.cblas.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#DCT: discrete cosine transforms
DCT: dct idct dct.cblas idct.cblas dct.fftw idct.fftw dct.ffts #idct.ffts
dct: dct.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idct: idct.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dct.cblas: dct.cblas.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idct.cblas: idct.cblas.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dct.fftw: dct.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idct.fftw: idct.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dct.ffts: dct.ffts.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
#idct.ffts: idct.ffts.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#DST: discrete sine transforms
DST: dst idst dst.cblas idst.cblas dst.fftw idst.fftw
dst: dst.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idst: idst.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dst.cblas: dst.cblas.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idst.cblas: idst.cblas.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dst.fftw: dst.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idst.fftw: idst.fftw.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Hilbert: Hilbert transform and related
Hilbert: hilbert analytic_sig analytic_amp analytic_pow inst_phase inst_freq
hilbert: hilbert.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
analytic_sig: analytic_sig.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
analytic_amp: analytic_amp.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
analytic_pow: analytic_pow.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
inst_phase: inst_phase.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
inst_freq: inst_freq.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Filter: FIR and IIR filters
Filter: fir iir filter fir_fft #fir_tau iir_tau filter filter_tau filtfilt spencer
fir: fir.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iir: iir.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
filter: filter.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
fir_fft: fir_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
filtfilt: filtfilt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Conv: convolution
Conv: conv xcorr conv1d xcorr1d conv_fft xcorr_fft conv1d_fft xcorr1d_fft
conv: conv.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
xcorr: xcorr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv1d: conv1d.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
xcorr1d: xcorr1d.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv_fft: conv_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
xcorr_fft: xcorr_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv1d_fft: conv1d_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
xcorr1d_fft: xcorr1d_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Interp: interpolation
Interp: interp1q #interp1ft
interp1q: interp1q.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
interp1ft: interp1ft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#ZC_LCs: zero-crossings, level-crossings and mean-crossings
ZCs_LCs: zcs lcs mcs zcr lcr mcr zcr_windowed lcr_windowed mcr_windowed
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
AR_Poly: poly2roots roots2poly poly2ar ar2poly rc2ar poly2ar ar2psd poly2psd ar2rc rc2ar poly2rc rc2poly #ar2cc poly2cc
poly2roots: poly2roots.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS) -Wno-c99-extensions
roots2poly: roots2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2ar: poly2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2poly: ar2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2rc: poly2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rc2poly: rc2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2rc: ar2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rc2ar: rc2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2psd: ar2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2psd: poly2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Linear_Pred: Linear prediction and related (autocorrelation, Levinson-Durbin, Burg, MVDR)
Linear_Pred: AC_LP Burg MVDR

#AC_LP: conversions related to AC (autocovariance) and LP (linear prediction).
#These work from the signal (sig) or autocovariance (AC) to get
#AC, RC (reflection coeff), poly (polynomial), or AR (autoregressive) representations.
#These use time-domain computations, except the _fft versions.
#The LP is done by Levinson-Durbin recursion from the AC.
AC_LP: sig2ac sig2ac_fft ac2rc ac2ar ac2poly ac2psd sig2rc sig2ar sig2poly sig2psd
sig2ac: sig2ac.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ac_fft: sig2ac_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2rc: ac2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2ar: ac2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2poly: ac2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2psd: ac2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2rc: sig2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ar: sig2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2poly: sig2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2psd: sig2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#Burg: Linear prediction by time-domain maximum-entropy method of Burg
Burg: #sig2ar_burg sig2poly_burg sig2psd_burg sig2cc_burg
sig2ar_burg: sig2ar_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2poly_burg: sig2poly_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2psd_burg: sig2psd_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
# sig2ar_fft_burg: sig2ar_fft_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
# sig2poly_fft_burg: sig2poly_fft_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
# sig2psd_fft_burg: sig2psd_fft_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)

#MVDR: minimum variance distortionless response
MVDR: #ac2mvdr sig2mvdr
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
Spectrogram: convert_freqs pow_compress
convert_freqs: convert_freqs.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
pow_compress: pow_compress.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Wavelets: a couple of my most often-used wavelets
Wavelets: #gabor analytic
gabor: gabor.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
analytic: analytic.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Nonlinear: various nonlinear DSP methods
Nonlinear: tkeo #medfilt
tkeo: tkeo.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
medfilt: medfilt.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#If want to make a library archive for C++ use, this works with clang++
#ar crs ../lib/libdsp.a obj/*.o


#make clean
Clean: clean
clean:
	find ./obj -type f -name *.o | xargs rm -f
	rm -f *.s tmp.*
