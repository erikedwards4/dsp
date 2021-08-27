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
CC=g++

ifeq ($(CC),clang++)
	STD=-std=c++11
	WFLAG=-Weverything -Wno-c++98-compat -Wno-old-style-cast -Wno-gnu-imaginary-constant
else
	STD=-std=gnu++14
	WFLAG=-Wall -Wextra
endif

INCLS=-Ic -I../util
CFLAGS=$(WFLAG) $(STD) -O3 -ffast-math -march=native -mfpmath=sse $(INCLS)


All: all
all: Dirs Generate Transform Filter Conv Interp ZCs_LCs AR_Poly AC_LP Frame STFT Clean
	rm -f 7 obj/*.o

Dirs:
	mkdir -pm 777 bin obj


#Generate: generate noise, simple waveforms, windows
Generate: Noise Pulses Waves Wins

#Noise: generate colored noise
Noise: white pink red brown blue violet
white: srci/white.cpp c/white.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
pink: srci/pink.cpp c/pink.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
red: srci/red.cpp c/red.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
brown: srci/brown.cpp c/brown.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
blue: srci/blue.cpp c/blue.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
violet: srci/violet.cpp c/violet.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

#Pulses: generate pulse functions (unit impulse is a.k.a. dirac delta)
Pulses: unit_impulse delta_impulse rect_pulse
unit_impulse: srci/unit_impulse.cpp c/unit_impulse.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
delta_impulse: srci/delta_impulse.cpp c/delta_impulse.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
rect_pulse: srci/rect_pulse.cpp c/rect_pulse.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2

#Waves: generate periodic waveforms
Waves: sinewave cosinewave squarewave triwave sawwave pulsewave
sinewave: srci/sinewave.cpp c/sinewave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
cosinewave: srci/cosinewave.cpp c/cosinewave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
squarewave: srci/squarewave.cpp c/squarewave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
triwave: srci/triwave.cpp c/triwave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
sawwave: srci/sawwave.cpp c/sawwave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
pulsewave: srci/pulsewave.cpp c/pulsewave.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

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
Transform: FFT DCT DST Hilbert

#FFT: fast Fourier transforms
FFT: fft ifft fft.rad2 ifft.rad2 fft.fftw ifft.fftw fft.fftw.r2hc fft.ffts ifft.ffts fft.kiss ifft.kiss
fft: srci/fft.cpp c/fft.fftw.c c/fft.rad2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
ifft: srci/ifft.cpp c/ifft.fftw.c c/ifft.rad2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
fft.rad2: srci/fft.rad2.cpp c/fft.rad2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ifft.rad2: srci/ifft.rad2.cpp c/ifft.rad2.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
fft.fftw: srci/fft.fftw.cpp c/fft.fftw.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
ifft.fftw: srci/ifft.fftw.cpp c/ifft.fftw.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
fft.fftw.r2hc: srci/fft.fftw.r2hc.cpp c/fft.fftw.r2hc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
fft.ffts: srci/fft.ffts.cpp c/fft.ffts.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lffts -lm
ifft.ffts: srci/ifft.ffts.cpp c/ifft.ffts.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lffts -lm
fft.kiss: srci/fft.kiss.cpp c/fft.kiss.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ifft.kiss: srci/ifft.kiss.cpp c/ifft.kiss.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm

#DCT: discrete cosine transforms
DCT: dct idct dct.cblas idct.cblas dct.fftw idct.fftw dct.ffts
dct: srci/dct.cpp c/dct.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
idct: srci/idct.cpp c/idct.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
dct.cblas: srci/dct.cblas.cpp c/dct.cblas.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
idct.cblas: srci/idct.cblas.cpp c/idct.cblas.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
dct.fftw: srci/dct.fftw.cpp c/dct.fftw.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
idct.fftw: srci/idct.fftw.cpp c/idct.fftw.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
dct.ffts: srci/dct.ffts.cpp c/dct.ffts.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lffts -lm

#DST: discrete sine transforms
DST: dst idst dst.cblas idst.cblas dst.fftw idst.fftw dst.ffts
dst: srci/dst.cpp c/dst.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
idst: srci/idst.cpp c/idst.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
dst.cblas: srci/dst.cblas.cpp c/dst.cblas.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
idst.cblas: srci/idst.cblas.cpp c/idst.cblas.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
dst.fftw: srci/dst.fftw.cpp c/dst.fftw.c
	$(ss) -tvd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
idst.fftw: srci/idst.fftw.cpp c/idst.fftw.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm


#Hilbert: Hilbert transform and related
Hilbert: hilbert analytic_sig analytic_amp analytic_pow inst_phase inst_freq
hilbert: srci/hilbert.cpp c/hilbert.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
analytic_sig: srci/analytic_sig.cpp c/analytic_sig.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
analytic_amp: srci/analytic_amp.cpp c/analytic_amp.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
analytic_pow: srci/analytic_pow.cpp c/analytic_pow.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
inst_phase: srci/inst_phase.cpp c/inst_phase.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
inst_freq: srci/inst_freq.cpp c/inst_freq.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm


#Filter: FIR and IIR filters
Filter: fir fir_fft iir filter filtfilt #fir_ola medfilt spencer
fir: srci/fir.cpp c/fir.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
fir_fft: srci/fir_fft.cpp c/fir_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
iir: srci/iir.cpp c/iir.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
filter: srci/filter.cpp c/filter.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm
filtfilt: srci/filtfilt.cpp c/filtfilt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lopenblas -lm


#Conv: 1-D convolution and cross-correlation
Conv: conv xcorr conv1d xcorr1d conv_fft xcorr_fft conv1d_fft xcorr1d_fft
conv: srci/conv.cpp c/conv.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
xcorr: srci/xcorr.cpp c/xcorr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
conv1d: srci/conv1d.cpp c/conv1d.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
xcorr1d: srci/xcorr1d.cpp c/xcorr1d.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
conv_fft: srci/conv_fft.cpp c/conv_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
xcorr_fft: srci/xcorr_fft.cpp c/xcorr_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
conv1d_fft: srci/conv1d_fft.cpp c/conv1d_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
xcorr1d_fft: srci/xcorr1d_fft.cpp c/xcorr1d_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm


#Interp: 1-D interpolation
Interp: interp1q #interp1ft
interp1q: srci/interp1q.cpp c/interp1q.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
interp1ft: srci/interp1ft.cpp c/interp1ft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2


#ZC_LCs: zero-crossings, level-crossings and mean-crossings
ZCs_LCs: zcs lcs mcs zcr lcr mcr zcr_windowed lcr_windowed mcr_windowed
zcs: srci/zcs.cpp c/zcs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
lcs: srci/lcs.cpp c/lcs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
mcs: srci/mcs.cpp c/mcs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
zcr: srci/zcr.cpp c/zcr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
lcr: srci/lcr.cpp c/lcr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
mcr: srci/mcr.cpp c/mcr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
zcr_windowed: srci/zcr_windowed.cpp c/zcr_windowed.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
lcr_windowed: srci/lcr_windowed.cpp c/lcr_windowed.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
mcr_windowed: srci/mcr_windowed.cpp c/mcr_windowed.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm


#AR_Poly: conversion between AR (autoregressive), poly (polynomial),
#RC (reflection coeff), and PSD (power spectral density) representations
#poly2roots uses LAPACKE, so for complex case requires -Wno-c99-extensions.
AR_Poly: #poly2roots roots2poly poly2ar ar2poly ar2psd poly2psd #ar2rc rc2ar poly2rc rc2poly
poly2roots: srci/poly2roots.cpp c/poly2roots.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS) -Wno-c99-extensions; $(CC) obj/$@.o -obin/$@ -largtable2 -llapacke -lm
roots2poly: srci/roots2poly.cpp c/roots2poly.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
poly2ar: srci/poly2ar.cpp c/poly2ar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
ar2poly: srci/ar2poly.cpp c/ar2poly.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
ar2psd: srci/ar2psd.cpp c/ar2psd.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
poly2psd: srci/poly2psd.cpp c/poly2psd.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
#ar2rc: srci/ar2rc.cpp c/ar2rc.c
#	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
#rc2ar: srci/rc2ar.cpp c/rc2ar.c
#	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
#poly2rc: srci/poly2rc.cpp c/poly2rc.c
#	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
#rc2poly: srci/rc2poly.cpp c/rc2poly.c
#	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2


#AC_LP: conversions between sig (signal), AC (autocorrelation), LP (linear prediction), and related.
#For example, ac2rc converts from AC to RCs (reflection coeffs).
AC_LP: #sig2ac sig2ac_fft ac2rc ac2ar ac2poly sig2rc sig2ar sig2poly #sig2ar_burg sig2poly_burg ac2cc ac2mvdr
ac2rc: srci/ac2rc.cpp c/ac2rc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ac2ar: srci/ac2ar.cpp c/ac2ar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ac2poly: srci/ac2poly.cpp c/ac2poly.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
ac2mvdr: srci/ac2mvdr.cpp c/ac2mvdr.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
sig2ac: srci/sig2ac.cpp c/sig2ac.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
sig2ac_fft: srci/sig2ac_fft.cpp c/sig2ac_fft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
sig2rc: srci/sig2rc.cpp c/sig2rc.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
sig2ar: srci/sig2ar.cpp c/sig2ar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
sig2poly: srci/sig2poly.cpp c/sig2poly.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
sig2ar_burg: srci/sig2ar_burg.cpp c/sig2ar_burg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
sig2poly_burg: srci/sig2poly_burg.cpp c/sig2poly_burg.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm


#Frame: get frames and apply windows for univariate signal to put into matrix
Frame: frame_univar frame_univar_flt apply_win window_univar window_univar_flt
frame_univar: srci/frame_univar.cpp c/frame_univar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
frame_univar_flt: srci/frame_univar_flt.cpp c/frame_univar_flt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
apply_win: srci/apply_win.cpp c/apply_win.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
window_univar: srci/window_univar.cpp c/window_univar.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
window_univar_flt: srci/window_univar_flt.cpp c/window_univar_flt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2


#STFT: steps to do the STFT (short-term Fourier transform)
STFT: fft_power get_stft_freqs stft stft_flt
fft_power: srci/fft_power.cpp c/fft_power.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
get_stft_freqs: srci/get_stft_freqs.cpp c/get_stft_freqs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2
stft: srci/stft.cpp c/stft.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm
stft_flt: srci/stft_flt.cpp c/stft_flt.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lfftw3f -lfftw3 -lm


#Spectrogram: steps to make a spectrogram (frequency-reweighted STFT)
Spectrogram: convert_freqs pow_compress
convert_freqs: srci/convert_freqs.cpp c/convert_freqs.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm
pow_compress: srci/pow_compress.cpp c/pow_compress.c
	$(ss) -vd srci/$@.cpp > src/$@.cpp; $(CC) -c src/$@.cpp -oobj/$@.o $(CFLAGS); $(CC) obj/$@.o -obin/$@ -largtable2 -lm


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
