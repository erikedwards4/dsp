#pragma once

//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <complex.h>
//#include <lapacke.h>
//#include <cblas.h>
#include <time.h>


#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ac2ar_s (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2ar_d (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2ar_c (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2ar_z (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int ac2mvdr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const float preg);
int ac2mvdr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const double preg);
int ac2mvdr_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const float preg);
int ac2mvdr_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const double preg);

int ac2poly_s (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2poly_d (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2poly_c (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2poly_z (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int ac2psd_s (float *Y, const float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2psd_d (double *Y, const double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2psd_c (float *Y, const float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2psd_z (double *Y, const double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int ac2rc_s (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_d (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_c (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_z (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int analytic_amp_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);
int analytic_amp_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);

int analytic_pow_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);
int analytic_pow_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);

int analytic_sig_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);
int analytic_sig_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);

int ar2psd_s (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2psd_d (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2psd_c (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2psd_z (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int ar2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int bartlett_s (float *Y, const size_t L, const size_t norm);
int bartlett_d (double *Y, const size_t L, const size_t norm);
int bartlett_c (float *Y, const size_t L, const size_t norm);
int bartlett_z (double *Y, const size_t L, const size_t norm);

int blackman_s (float *Y, const size_t L, const int exact, const size_t norm);
int blackman_d (double *Y, const size_t L, const int exact, const size_t norm);
int blackman_c (float *Y, const size_t L, const int exact, const size_t norm);
int blackman_z (double *Y, const size_t L, const int exact, const size_t norm);

int blackmanharris_s (float *Y, const size_t L, const size_t norm);
int blackmanharris_d (double *Y, const size_t L, const size_t norm);
int blackmanharris_c (float *Y, const size_t L, const size_t norm);
int blackmanharris_z (double *Y, const size_t L, const size_t norm);

int blue_s (float *Y, const size_t N, const float std, const int zmn);
int blue_d (double *Y, const size_t N, const double std, const int zmn);
int blue_c (float *Y, const size_t N, const float std, const int zmn);
int blue_z (double *Y, const size_t N, const double std, const int zmn);

int brown_s (float *Y, const size_t N, const float std, const int zmn);
int brown_d (double *Y, const size_t N, const double std, const int zmn);
int brown_c (float *Y, const size_t N, const float std, const int zmn);
int brown_z (double *Y, const size_t N, const double std, const int zmn);

int conv1_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);
int conv1_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);
int conv1_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);
int conv1_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);

int conv1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);

int conv1d_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_fft_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_fft_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int conv1d_fft_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);

int conv_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);

int convert_freqs_s (float *frqs, const size_t F, const char in_scale[], const char out_scale[]);
int convert_freqs_d (double *frqs, const size_t F, const char in_scale[], const char out_scale[]);

int conv_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_fft_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_fft_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_fft_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);

int cosinewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int cosinewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int cosinewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int cosinewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);

int dct_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);

int dct_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);

int dct_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);

int dct_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);

int dft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);
int dft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);
int dft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);
int dft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);

int dft_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);
int dft_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);
int dft_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);
int dft_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc);

int dst_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);

int dst_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);

int dst_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int dst_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);

int fft_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int fft_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int fft_kiss_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_kiss_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_kiss_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_kiss_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int fft_power_c (float *Y, const float *X, const size_t N);
int fft_power_z (double *Y, const double *X, const size_t N);

int fft_rad2_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_rad2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_rad2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_rad2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int fft_squared_s (float *Y, const float *X, const int iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_squared_d (double *Y, const double *X, const int iscolmajor, const int R, const int C, const int dim, const int nfft);

int filter_s (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filter_d (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filter_c (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filter_z (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);

int filtfilt_s (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filtfilt_d (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filtfilt_c (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filtfilt_z (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);

int fir_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);

int fir_fft_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_fft_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_fft_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_fft_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);

int flattop_s (float *Y, const size_t L, const size_t norm);
int flattop_d (double *Y, const size_t L, const size_t norm);
int flattop_c (float *Y, const size_t L, const size_t norm);
int flattop_z (double *Y, const size_t L, const size_t norm);

int frame_univar_s (float *Y, const float *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int frame_univar_d (double *Y, const double *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int frame_univar_c (float *Y, const float *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int frame_univar_z (double *Y, const double *X, const size_t N, const size_t L, const size_t stp, const int snip_edges);

int gauss_s (float *Y, const size_t L, const float r, const size_t norm);
int gauss_d (double *Y, const size_t L, const double r, const size_t norm);
int gauss_c (float *Y, const size_t L, const float r, const size_t norm);
int gauss_z (double *Y, const size_t L, const double r, const size_t norm);

int hamming_s (float *Y, const size_t L, const size_t norm);
int hamming_d (double *Y, const size_t L, const size_t norm);
int hamming_c (float *Y, const size_t L, const size_t norm);
int hamming_z (double *Y, const size_t L, const size_t norm);

int hann_s (float *Y, const size_t L, const size_t norm);
int hann_d (double *Y, const size_t L, const size_t norm);
int hann_c (float *Y, const size_t L, const size_t norm);
int hann_z (double *Y, const size_t L, const size_t norm);

int hilbert_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);
int hilbert_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);

int idct_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);

int idct_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);

int idct_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int idct_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);

int idft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);

int idft_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);

int idst_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);

int idst_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);

int idst_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);

int ifft_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int ifft_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int ifft_kiss_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_kiss_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_kiss_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_kiss_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int ifft_rad2_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_rad2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_rad2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_rad2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);

int iir_s (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_d (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_c (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_z (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_inplace_s (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_inplace_d (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_inplace_c (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_inplace_z (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);

int inst_freq_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);
int inst_freq_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);

int inst_phase_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);
int inst_phase_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft);

int interp1q_s (float *Yi, const float *X, const float *Y, const float *Xi, const size_t N, const size_t Ni, const int decreasing);
int interp1q_d (double *Yi, const double *X, const double *Y, const double *Xi, const size_t N, const size_t Ni, const int decreasing);
int interp1q_c (float *Yi, const float *X, const float *Y, const float *Xi, const size_t N, const size_t Ni, const int decreasing);
int interp1q_z (double *Yi, const double *X, const double *Y, const double *Xi, const size_t N, const size_t Ni, const int decreasing);

int lcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const float lvl, const int causal);
int lcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const double lvl, const int causal);

int lcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going, const float lvl);
int lcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going, const double lvl);

int lcs_s (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going, const float lvl);
int lcs_d (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going, const double lvl);

int mcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int mcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);

int mcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going);
int mcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going);

int mcs_s (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);
int mcs_d (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);

int planck_s (float *Y, const size_t L, const float epsilon, const size_t norm);
int planck_d (double *Y, const size_t L, const double epsilon, const size_t norm);
int planck_c (float *Y, const size_t L, const float epsilon, const size_t norm);
int planck_z (double *Y, const size_t L, const double epsilon, const size_t norm);

int poly2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int poly2psd_s (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2psd_d (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2psd_c (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2psd_z (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int poly2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);

int poly2roots_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2roots_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2roots_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int poly2roots_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int povey_s (float *Y, const size_t L, const size_t norm);
int povey_d (double *Y, const size_t L, const size_t norm);
int povey_c (float *Y, const size_t L, const size_t norm);
int povey_z (double *Y, const size_t L, const size_t norm);

int pow_compress_s (float *Y, const float *X, const size_t N, const float p, const float preg);
int pow_compress_d (double *Y, const double *X, const size_t N, const double p, const double preg);
int pow_compress_inplace_s (float *X, const size_t N, const float p, const float preg);
int pow_compress_inplace_d (double*X, const size_t N, const double p, const double preg);

int pulsewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs, const float dty);
int pulsewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty);
int pulsewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs, const float dty);
int pulsewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs, const double dty);

int rc2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rc2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rc2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rc2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int rc2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rc2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rc2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int rc2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int rectangular_s (float *Y, const size_t L, const size_t norm);
int rectangular_d (double *Y, const size_t L, const size_t norm);
int rectangular_c (float *Y, const size_t L, const size_t norm);
int rectangular_z (double *Y, const size_t L, const size_t norm);

int red_s (float *Y, const size_t N, const float std, const int zmn);
int red_d (double *Y, const size_t N, const double std, const int zmn);
int red_c (float *Y, const size_t N, const float std, const int zmn);
int red_z (double *Y, const size_t N, const double std, const int zmn);

int roots2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int roots2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int roots2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int roots2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int sawwave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sawwave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int sawwave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sawwave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);

int sig2ac_s (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_d (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_c (float *Y, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_z (double *Y, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);

int sig2ac_fft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_fft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);
int sig2ac_fft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Ly, const int mnz, const int unbiased, const int corr);

int sig2ar_burg_s (float *Y, float *V, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mean0);
int sig2ar_burg_d (double *Y, double *V, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mean0);

int sig2ar_s (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2ar_d (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2ar_c (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2ar_z (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);

int sig2poly_burg_s (float *Y, float *V, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mean0);
int sig2poly_burg_d (double *Y, double *V, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mean0);

int sig2poly_s (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased);
int sig2poly_d (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased);
int sig2poly_c (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased);
int sig2poly_z (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased);

int sig2psd_s (float *Y, float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2psd_d (double *Y, double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2psd_c (float *Y, float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2psd_z (double *Y, double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);

int sinewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sinewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int sinewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int sinewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);

int smooth_diff_s (float *Y, const size_t N, const size_t n);
int smooth_diff_d (double *Y, const size_t N, const size_t n);

int smooth_diffdiff_s (float *Y, const size_t N, const size_t n);
int smooth_diffdiff_d (double *Y, const size_t N, const size_t n);

int spencer_s (float *Y);
int spencer_d (double *Y);

int squarewave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int squarewave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int squarewave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int squarewave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);

int stft_flt_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t W, const size_t nfft, const float c0, const float stp, const int mn0, const int amp, const int lg);
int stft_flt_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t W, const size_t nfft, const double c0, const double stp, const int mn0, const int amp, const int lg);

int tkeo_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int tkeo_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int tkeo_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int tkeo_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);

int tkeo_smooth_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly);
int tkeo_smooth_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly);
int tkeo_smooth_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly);
int tkeo_smooth_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int cs0, const size_t str, const size_t dil, const size_t Ly);

int triangular_s (float *Y, const size_t L, const size_t norm);
int triangular_d (double *Y, const size_t L, const size_t norm);
int triangular_c (float *Y, const size_t L, const size_t norm);
int triangular_z (double *Y, const size_t L, const size_t norm);

int triwave_s (float *Y, const size_t N, const float amp, const float frq, const float phs);
int triwave_d (double *Y, const size_t N, const double amp, const double frq, const double phs);
int triwave_c (float *Y, const size_t N, const float amp, const float frq, const float phs);
int triwave_z (double *Y, const size_t N, const double amp, const double frq, const double phs);

int tukey_s (float *Y, const size_t L, const float r, const size_t norm);
int tukey_d (double *Y, const size_t L, const double r, const size_t norm);
int tukey_c (float *Y, const size_t L, const float r, const size_t norm);
int tukey_z (double *Y, const size_t L, const double r, const size_t norm);

int white_s (float *Y, const size_t N, const float std, const int uni, const int zmn);
int white_d (double *Y, const size_t N, const double std, const int uni, const int zmn);
int white_c (float *Y, const size_t N, const float std, const int uni, const int zmn);
int white_z (double *Y, const size_t N, const double std, const int uni, const int zmn);

int window_univar_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int window_univar_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int window_univar_c (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int window_univar_z (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges);

int window_univar_flt_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t W, const float c0, const float stp);
int window_univar_flt_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t W, const double c0, const double stp);
int window_univar_flt_c (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t W, const float c0, const float stp);
int window_univar_flt_z (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t W, const double c0, const double stp);

int xcorr1_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);
int xcorr1_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);
int xcorr1_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);
int xcorr1_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int es0, const size_t str, const size_t dil, const size_t Ly, const size_t dim);

int xcorr1d_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const int pad, const size_t str, const size_t dil, const size_t dim);

int xcorr1d_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_fft_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_fft_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);
int xcorr1d_fft_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim);

int xcorr_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int xcorr_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int xcorr_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int xcorr_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);

int xcorr_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int xcorr_fft_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int xcorr_fft_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int xcorr_fft_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);

int zcr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int zcr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int zcr_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);
int zcr_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t Lw, const int going, const int causal);

int zcr_windowed_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going);
int zcr_windowed_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going);
int zcr_windowed_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const float c0, const float stp, const int going);
int zcr_windowed_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Lw, const size_t W, const size_t dim, const double c0, const double stp, const int going);

int zcs_s (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);
int zcs_d (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);
int zcs_c (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);
int zcs_z (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);

#ifdef __cplusplus
}
}
#endif
