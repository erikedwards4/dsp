#@author Erik Edwards
#@date 2017-2020
#@license BSD 3-clause

#DSP is my own library of C functions for digital signal processing in C.

#This is the Makefile used during development of the C functions.
#Since these are not stand-alone programs (no main functions), the only outputs are object files to obj.


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

CFLAGS=$(WFLAG) $(STD) -O3 -march=native -fPIC


All: all
all: Wins Interp Filter Conv Transform AR_Poly AC_LP Clean


#Wins: common 1-D windows
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


Interp: #interp1q
interp1q: interp1q.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Filter: FIR and IIR filters
Filter: fir iir #fir_tau iir_tau filter filter_tau filtfilt spencer
fir: fir.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
iir: iir.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Conv: convolution
Conv: #conv1d conv2d conv1d_fft conv2d_fft
conv1d: conv1d.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv2d: conv2d.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv1d_fft: conv1d_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
conv2d_fft: conv2d_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#Transform: common 1-D signal transforms
Transforms: fft ifft dct idct dst idst hilbert
fft: fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ifft: ifft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dct: dct.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idct: idct.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
dst: dst.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
idst: idst.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
hilbert: hilbert.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#AR_Poly: conversion between autoregressive (AR), polynomial, reflection coeff (RC)
#and power spectral density (PSD) representations
AR_Poly: poly2roots roots2poly ar2poly poly2ar ar2rc rc2ar poly2rc rc2poly ar2psd poly2psd
poly2roots: poly2roots.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
roots2poly: roots2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2poly: ar2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2ar: poly2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2rc: ar2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rc2ar: rc2ar.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2rc: poly2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
rc2poly: rc2poly.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ar2psd: ar2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
poly2psd: poly2psd.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)


#AC_LP: conversions involving autocorrelation (AC), linear prediction (LP), and related,
#for example, ac2cc converts from AC to cepstral coeffs (CCs).
AC_LP: autocorr autocorr_fft sig2ac sig2ac_fft ac2ar_levdurb ac2poly_levdurb sig2ar_levdurb sig2poly_levdurb sig2ar_burg sig2poly_burg ac2rc ac2cc ac2mvdr
sig2ac: sig2ac.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ac_fft: sig2ac_fft.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2ar_levdurb: ac2ar_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2poly_levdurb: ac2poly_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ar_levdurb: sig2ar_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2poly_levdurb: sig2poly_levdurb.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2ar_burg: sig2ar_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
sig2poly_burg: sig2poly_burg.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2rc: ac2rc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2cc: ac2cc.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)
ac2mvdr: ac2mvdr.c; $(CC) -c $@.c -o obj/$@.o $(CFLAGS)



#If want to make a library archive for C++ use, this works with clang++
#ar crs ../lib/libdsp.a obj/*.o


#make clean
Clean: clean
clean:
	find ./obj -type f -name *.o | xargs rm -f
	rm -f *.s tmp.*
