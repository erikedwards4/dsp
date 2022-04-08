//1D cross-correlation of each vector in X1f by X2 along dim.
//Each vector in X1f has length L1. X2 has length L2.

//Note that some "convolution" functions actually do cross-correlation.

//This version uses the FFT for computation.
//It is faster than xcorr only if L1+L2 is a power-of-2, or just less than a power-of-2.
//Even then, it is slower including the FFT set-up;
//so would only make sense when one could do the FFT the set-up once and xcorr repeatedly.

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fftw3.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int xcorr_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr_fft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L1<1u) { fprintf(stderr,"error in xcorr_fft_s: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr_fft_s: L2 (length X2) must be positive\n"); return 1; }

    //Set Ly, ss, es, zss according to shape
    size_t Ly, zss;
    int ss, es;
    if (strncmp(shape,"full",4u)==0)
    {
        Ly = L1 + L2 - 1u;
        es = 0; ss = 1 - (int)L2;
        zss = 0u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        Ly = L1;
        es = (int)L2/2; ss = es - (int)L2 + 1;
        zss = L2/2u;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        Ly = (L1<L2) ? 0u : L1-L2+1u;
        ss = 0; es = (int)L2 - 1;
        zss = L2 - 1u;
    }
    else
    {
        fprintf(stderr,"error in xcorr_fft_s: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (N==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while (nfft<L1+L2) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;
        //if (nfft%2u) { fprintf(stderr,"error in xcorr_fft_s: nfft must be even"); return 1; }

        //Set scale
        const float sc = 2.0f / (float)nfft;

        //Initialize FFTs
        float *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (float *)fftwf_malloc(nfft*sizeof(float));
        X2f = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan fplan1 = fftwf_plan_dft_r2c_1d((int)nfft,X1f,(fftwf_complex *)Y1f,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in xcorr_fft_s: problem creating fftw plan"); return 1; }
        fftwf_plan fplan2 = fftwf_plan_dft_r2c_1d((int)nfft,X2f,(fftwf_complex *)Y2f,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in xcorr_fft_s: problem creating fftw plan"); return 1; }
        for (size_t n=L1; n<nfft; ++n) { X1f[n] = 0.0f; }

        //Initialize IFFT
        float *Z1f;
        Z1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1f,(fftwf_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in xcorr_fft_s: problem creating fftw plan"); return 1; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1f[n] = 0.0f; }

        //Get Y2f (scaled FFT of reversed X2f)
        X2 += L2;
        for (size_t n=L2; n>0u; --n, ++X2f) { *X2f = *--X2; }
        for (size_t n=L2; n<nfft; ++n, ++X2f) { *X2f = 0.0f; }
        X2f -= nfft;
        fftwf_execute(fplan2);
        *Y2f++ /= (float)nfft; ++Y2f;
        for (size_t l=nfft-2u; l>0u; --l, ++Y2f) { *Y2f *= sc; }
        *Y2f /= (float)nfft;
        Y2f -= nfft;

        if (L1==N)
        {
            //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

            //Copy input
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
            X1f -= L1;

            //FFT
            fftwf_execute(fplan1);

            //Multiply FFTs
            *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
            for (size_t f=F-2u; f>0u; --f)
            {
                y1r = *Y1f++; y1i = *Y1f--;
                y2r = *Y2f++; y2i = *Y2f++;
                *Y1f++ = y2r*y1r - y2i*y1i;
                *Y1f++ = y2r*y1i + y2i*y1r;
            }
            *Y1f *= *Y2f;
            Y1f -= nfft; Y2f -= nfft;

            //IFFT
            fftwf_execute(iplan);

            //Copy output
            Z1f += 2u*zss;
            for (size_t l=Ly; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
            Z1f -= 2u*(zss+Ly);

            //clock_gettime(CLOCK_REALTIME,&toc);
            //fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    //Copy input
                    for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
                    X1f -= L1;

                    //FFT
                    fftwf_execute(fplan1);

                    //Multiply FFTs
                    *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                    for (size_t f=F-2u; f>0u; --f)
                    {
                        y1r = *Y1f++; y1i = *Y1f--;
                        y2r = *Y2f++; y2i = *Y2f++;
                        *Y1f++ = y2r*y1r - y2i*y1i;
                        *Y1f++ = y2r*y1i + y2i*y1r;
                    }
                    *Y1f *= *Y2f;
                    Y1f -= nfft; Y2f -= nfft;

                    //IFFT
                    fftwf_execute(iplan);

                    //Copy output
                    Z1f += 2u*zss;
                    for (size_t l=Ly; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
                    Z1f -= 2u*(zss+Ly);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X1-=K*L1-1u, Y-=K*Ly-1u)
                    {
                        //Copy input
                        for (size_t l=L1; l>0u; --l, X1+=K, ++X1f) { *X1f = *X1; }
                        X1f -= L1;

                        //FFT
                        fftwf_execute(fplan1);

                        //Multiply FFTs
                        *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                        for (size_t f=F-2u; f>0u; --f)
                        {
                            y1r = *Y1f++; y1i = *Y1f--;
                            y2r = *Y2f++; y2i = *Y2f++;
                            *Y1f++ = y2r*y1r - y2i*y1i;
                            *Y1f++ = y2r*y1i + y2i*y1r;
                        }
                        *Y1f *= *Y2f;
                        Y1f -= nfft; Y2f -= nfft;

                        //IFFT
                        fftwf_execute(iplan);

                        //Copy output
                        Z1f += 2u*zss;
                        for (size_t l=Ly; l>0u; --l, Z1f+=2, Y+=K) { *Y = *Z1f; }
                        Z1f -= 2u*(zss+Ly);
                    }
                }
            }
        }

        //Free
        fftwf_free(X1f); fftwf_free(Y2f); fftwf_free(Y1f); fftwf_free(Z1f);
        fftwf_destroy_plan(fplan1); fftwf_destroy_plan(fplan2); fftwf_destroy_plan(iplan);
    }

    return 0;
}


int xcorr_fft_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr_fft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L1<1u) { fprintf(stderr,"error in xcorr_fft_d: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr_fft_d: L2 (length X2) must be positive\n"); return 1; }

    //Set Ly, ss, es, zss according to shape
    size_t Ly, zss;
    int ss, es;
    if (strncmp(shape,"full",4u)==0)
    {
        Ly = L1 + L2 - 1u;
        es = 0; ss = 1 - (int)L2;
        zss = 0u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        Ly = L1;
        es = (int)L2/2; ss = es - (int)L2 + 1;
        zss = L2/2u;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        Ly = (L1<L2) ? 0u : L1-L2+1u;
        ss = 0; es = (int)L2 - 1;
        zss = L2 - 1u;
    }
    else
    {
        fprintf(stderr,"error in xcorr_fft_d: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (N==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while (nfft<L1+L2) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;
        //if (nfft%2u) { fprintf(stderr,"error in xcorr_fft_d: nfft must be even"); return 1; }

        //Set scale
        const double sc = 2.0 / (double)nfft;

        //Initialize FFTs
        double *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (double *)fftw_malloc(nfft*sizeof(double));
        X2f = (double *)fftw_malloc(nfft*sizeof(double));
        Y1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan fplan1 = fftw_plan_dft_r2c_1d((int)nfft,X1f,(fftw_complex *)Y1f,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in xcorr_fft_d: problem creating fftw plan"); return 1; }
        fftw_plan fplan2 = fftw_plan_dft_r2c_1d((int)nfft,X2f,(fftw_complex *)Y2f,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in xcorr_fft_d: problem creating fftw plan"); return 1; }
        for (size_t n=L1; n<nfft; ++n) { X1f[n] = 0.0; }

        //Initialize IFFT
        double *Z1f;
        Z1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1f,(fftw_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in xcorr_fft_d: problem creating fftw plan"); return 1; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1f[n] = 0.0; }

        //Get Y2f (scaled FFT of reversed X2f)
        X2 += L2;
        for (size_t n=L2; n>0u; --n, ++X2f) { *X2f = *--X2; }
        for (size_t n=L2; n<nfft; ++n, ++X2f) { *X2f = 0.0; }
        X2f -= nfft;
        fftw_execute(fplan2);
        *Y2f++ /= (double)nfft; ++Y2f;
        for (size_t l=nfft-2u; l>0u; --l, ++Y2f) { *Y2f *= sc; }
        *Y2f /= (double)nfft;
        Y2f -= nfft;

        if (L1==N)
        {
            //Copy input
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
            X1f -= L1;

            //FFT
            fftw_execute(fplan1);

            //Multiply FFTs
            *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
            for (size_t f=F-2u; f>0u; --f)
            {
                y1r = *Y1f++; y1i = *Y1f--;
                y2r = *Y2f++; y2i = *Y2f++;
                *Y1f++ = y2r*y1r - y2i*y1i;
                *Y1f++ = y2r*y1i + y2i*y1r;
            }
            *Y1f *= *Y2f;
            Y1f -= nfft; Y2f -= nfft;

            //IFFT
            fftw_execute(iplan);

            //Copy output
            Z1f += 2u*zss;
            for (size_t l=Ly; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
            Z1f -= 2u*(zss+Ly);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    //Copy input
                    for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
                    X1f -= L1;

                    //FFT
                    fftw_execute(fplan1);

                    //Multiply FFTs
                    *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                    for (size_t f=F-2u; f>0u; --f)
                    {
                        y1r = *Y1f++; y1i = *Y1f--;
                        y2r = *Y2f++; y2i = *Y2f++;
                        *Y1f++ = y2r*y1r - y2i*y1i;
                        *Y1f++ = y2r*y1i + y2i*y1r;
                    }
                    *Y1f *= *Y2f;
                    Y1f -= nfft; Y2f -= nfft;

                    //IFFT
                    fftw_execute(iplan);

                    //Copy output
                    Z1f += 2u*zss;
                    for (size_t l=Ly; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
                    Z1f -= 2u*(zss+Ly);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X1-=K*L1-1u, Y-=K*Ly-1u)
                    {
                        //Copy input
                        for (size_t l=L1; l>0u; --l, X1+=K, ++X1f) { *X1f = *X1; }
                        X1f -= L1;

                        //FFT
                        fftw_execute(fplan1);

                        //Multiply FFTs
                        *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                        for (size_t f=F-2u; f>0u; --f)
                        {
                            y1r = *Y1f++; y1i = *Y1f--;
                            y2r = *Y2f++; y2i = *Y2f++;
                            *Y1f++ = y2r*y1r - y2i*y1i;
                            *Y1f++ = y2r*y1i + y2i*y1r;
                        }
                        *Y1f *= *Y2f;
                        Y1f -= nfft; Y2f -= nfft;

                        //IFFT
                        fftw_execute(iplan);

                        //Copy output
                        Z1f += 2u*zss;
                        for (size_t l=Ly; l>0u; --l, Z1f+=2, Y+=K) { *Y = *Z1f; }
                        Z1f -= 2u*(zss+Ly);
                    }
                }
            }
        }

        //Free
        fftw_free(X1f); fftw_free(Y2f); fftw_free(Y1f); fftw_free(Z1f);
        fftw_destroy_plan(fplan1); fftw_destroy_plan(fplan2); fftw_destroy_plan(iplan);
    }

    return 0;
}


int xcorr_fft_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr_fft_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L1<1u) { fprintf(stderr,"error in xcorr_fft_c: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr_fft_c: L2 (length X2) must be positive\n"); return 1; }

    //Set Ly, ss, es, zss according to shape
    size_t Ly, zss;
    int ss, es;
    if (strncmp(shape,"full",4u)==0)
    {
        Ly = L1 + L2 - 1u;
        es = 0; ss = 1 - (int)L2;
        zss = 0u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        Ly = L1;
        es = (int)L2/2; ss = es - (int)L2 + 1;
        zss = L2/2u;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        Ly = (L1<L2) ? 0u : L1-L2+1u;
        ss = 0; es = (int)L2 - 1;
        zss = L2 - 1u;
    }
    else
    {
        fprintf(stderr,"error in xcorr_fft_c: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (N==0u) {}
    else
    {
        //Set nfft
        size_t nfft = 1u;
        while (nfft<L1+L2) { nfft *= 2u; }
        //if (nfft%2u) { fprintf(stderr,"error in xcorr_fft_c: nfft must be even"); return 1; }

        //Initialize FFTs
        float *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        X2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan fplan1 = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1f,(fftwf_complex *)Y1f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in xcorr_fft_c: problem creating fftw plan"); return 1; }
        fftwf_plan fplan2 = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X2f,(fftwf_complex *)Y2f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in xcorr_fft_c: problem creating fftw plan"); return 1; }
        for (size_t n=2u*L1; n<2u*nfft; ++n) { X1f[n] = 0.0f; }

        //Initialize IFFT
        float *Z1f;
        Z1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1f,(fftwf_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in xcorr_fft_c: problem creating fftw plan"); return 1; }

        //Get Y2f (scaled FFT of reversed X2f)
        X2 += 2u*L2;
        for (size_t n=L2; n>0u; --n, ++X2f) { X2-=2; *X2f = *X2; *++X2f = *(X2+1); }
        for (size_t n=L2; n<nfft; ++n, ++X2f) { *X2f = 0.0f; *++X2f = 0.0f; }
        X2f -= 2u*nfft;
        fftwf_execute(fplan2);
        for (size_t l=2u*nfft; l>0u; --l, ++Y2f) { *Y2f /= (float)nfft; }
        Y2f -= 2u*nfft;

        if (L1==N)
        {
            //Copy input
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; *++X1f = *++X1; }
            X1f -= 2u*L1;

            //FFT
            fftwf_execute(fplan1);

            //Multiply FFTs
            for (size_t f=nfft; f>0u; --f)
            {
                y1r = *Y1f++; y1i = *Y1f--;
                y2r = *Y2f++; y2i = *Y2f++;
                *Y1f++ = y2r*y1r - y2i*y1i;
                *Y1f++ = y2r*y1i + y2i*y1r;
            }
            Y1f -= 2u*nfft; Y2f -= 2u*nfft;

            //IFFT
            fftwf_execute(iplan);

            //Copy output
            Z1f += 2u*zss;
            for (size_t l=Ly; l>0u; --l, ++Z1f, ++Y) { *Y = *Z1f; *++Y = *++Z1f; }
            Z1f -= 2u*(zss+Ly);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    //Copy input
                    for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; *++X1f = *++X1; }
                    X1f -= 2u*L1;

                    //FFT
                    fftwf_execute(fplan1);

                    //Multiply FFTs
                    for (size_t f=nfft; f>0u; --f)
                    {
                        y1r = *Y1f++; y1i = *Y1f--;
                        y2r = *Y2f++; y2i = *Y2f++;
                        *Y1f++ = y2r*y1r - y2i*y1i;
                        *Y1f++ = y2r*y1i + y2i*y1r;
                    }
                    Y1f -= 2u*nfft; Y2f -= 2u*nfft;

                    //IFFT
                    fftwf_execute(iplan);

                    //Copy output
                    Z1f += 2u*zss;
                    for (size_t l=Ly; l>0u; --l, ++Z1f, ++Y) { *Y = *Z1f; *++Y = *++Z1f; }
                    Z1f -= 2u*(zss+Ly);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X1-=2u*K*L1-2u, Y-=2u*K*Ly-2u)
                    {
                        //Copy input
                        for (size_t l=L1; l>0u; --l, X1+=2u*K, ++X1f) { *X1f = *X1; *++X1f = *(X1+1); }
                        X1f -= 2u*L1;

                        //FFT
                        fftwf_execute(fplan1);

                        //Multiply FFTs
                        for (size_t f=nfft; f>0u; --f)
                        {
                            y1r = *Y1f++; y1i = *Y1f--;
                            y2r = *Y2f++; y2i = *Y2f++;
                            *Y1f++ = y2r*y1r - y2i*y1i;
                            *Y1f++ = y2r*y1i + y2i*y1r;
                        }
                        Y1f -= 2u*nfft; Y2f -= 2u*nfft;

                        //IFFT
                        fftwf_execute(iplan);

                        //Copy output
                        Z1f += 2u*zss;
                        for (size_t l=Ly; l>0u; --l, ++Z1f, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *++Z1f; }
                        Z1f -= 2u*(zss+Ly);
                    }
                }
            }
        }

        //Free
        fftwf_free(X1f); fftwf_free(Y2f); fftwf_free(Y1f); fftwf_free(Z1f);
        fftwf_destroy_plan(fplan1); fftwf_destroy_plan(fplan2); fftwf_destroy_plan(iplan);
    }

    return 0;
}


int xcorr_fft_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in xcorr_fft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L1<1u) { fprintf(stderr,"error in xcorr_fft_z: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in xcorr_fft_z: L2 (length X2) must be positive\n"); return 1; }

    //Set Ly, ss, es, zss according to shape
    size_t Ly, zss;
    int ss, es;
    if (strncmp(shape,"full",4u)==0)
    {
        Ly = L1 + L2 - 1u;
        es = 0; ss = 1 - (int)L2;
        zss = 0u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        Ly = L1;
        es = (int)L2/2; ss = es - (int)L2 + 1;
        zss = L2/2u;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        Ly = (L1<L2) ? 0u : L1-L2+1u;
        ss = 0; es = (int)L2 - 1;
        zss = L2 - 1u;
    }
    else
    {
        fprintf(stderr,"error in xcorr_fft_z: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (N==0u) {}
    else
    {
        //Set nfft
        size_t nfft = 1u;
        while (nfft<L1+L2) { nfft *= 2u; }
        //if (nfft%2u) { fprintf(stderr,"error in xcorr_fft_z: nfft must be even"); return 1; }

        //Initialize FFTs
        double *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        X2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan fplan1 = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1f,(fftw_complex *)Y1f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in xcorr_fft_z: problem creating fftw plan"); return 1; }
        fftw_plan fplan2 = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X2f,(fftw_complex *)Y2f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in xcorr_fft_z: problem creating fftw plan"); return 1; }
        for (size_t n=2u*L1; n<2u*nfft; ++n) { X1f[n] = 0.0; }

        //Initialize IFFT
        double *Z1f;
        Z1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1f,(fftw_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in xcorr_fft_z: problem creating fftw plan"); return 1; }

        //Get Y2f (scaled FFT of reversed X2f)
        X2 += 2u*L2;
        for (size_t n=L2; n>0u; --n, ++X2f) { X2-=2; *X2f = *X2; *++X2f = *(X2+1); }
        for (size_t n=L2; n<nfft; ++n, ++X2f) { *X2f = 0.0; *++X2f = 0.0; }
        X2f -= 2u*nfft;
        fftw_execute(fplan2);
        for (size_t l=2u*nfft; l>0u; --l, ++Y2f) { *Y2f /= (double)nfft; }
        Y2f -= 2u*nfft;

        if (L1==N)
        {
            //Copy input
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; *++X1f = *++X1; }
            X1f -= 2u*L1;

            //FFT
            fftw_execute(fplan1);

            //Multiply FFTs
            for (size_t f=nfft; f>0u; --f)
            {
                y1r = *Y1f++; y1i = *Y1f--;
                y2r = *Y2f++; y2i = *Y2f++;
                *Y1f++ = y2r*y1r - y2i*y1i;
                *Y1f++ = y2r*y1i + y2i*y1r;
            }
            Y1f -= 2u*nfft; Y2f -= 2u*nfft;

            //IFFT
            fftw_execute(iplan);

            //Copy output
            Z1f += 2u*zss;
            for (size_t l=Ly; l>0u; --l, ++Z1f, ++Y) { *Y = *Z1f; *++Y = *++Z1f; }
            Z1f -= 2u*(zss+Ly);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    //Copy input
                    for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; *++X1f = *++X1; }
                    X1f -= 2u*L1;

                    //FFT
                    fftw_execute(fplan1);

                    //Multiply FFTs
                    for (size_t f=nfft; f>0u; --f)
                    {
                        y1r = *Y1f++; y1i = *Y1f--;
                        y2r = *Y2f++; y2i = *Y2f++;
                        *Y1f++ = y2r*y1r - y2i*y1i;
                        *Y1f++ = y2r*y1i + y2i*y1r;
                    }
                    Y1f -= 2u*nfft; Y2f -= 2u*nfft;

                    //IFFT
                    fftw_execute(iplan);

                    //Copy output
                    Z1f += 2u*zss;
                    for (size_t l=Ly; l>0u; --l, ++Z1f, ++Y) { *Y = *Z1f; *++Y = *++Z1f; }
                    Z1f -= 2u*(zss+Ly);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X1-=2u*K*L1-2u, Y-=2u*K*Ly-2u)
                    {
                        //Copy input
                        for (size_t l=L1; l>0u; --l, X1+=2u*K, ++X1f) { *X1f = *X1; *++X1f = *(X1+1); }
                        X1f -= 2u*L1;

                        //FFT
                        fftw_execute(fplan1);

                        //Multiply FFTs
                        for (size_t f=nfft; f>0u; --f)
                        {
                            y1r = *Y1f++; y1i = *Y1f--;
                            y2r = *Y2f++; y2i = *Y2f++;
                            *Y1f++ = y2r*y1r - y2i*y1i;
                            *Y1f++ = y2r*y1i + y2i*y1r;
                        }
                        Y1f -= 2u*nfft; Y2f -= 2u*nfft;

                        //IFFT
                        fftw_execute(iplan);

                        //Copy output
                        Z1f += 2u*zss;
                        for (size_t l=Ly; l>0u; --l, ++Z1f, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *++Z1f; }
                        Z1f -= 2u*(zss+Ly);
                    }
                }
            }
        }

        //Free
        fftw_free(X1f); fftw_free(Y2f); fftw_free(Y1f); fftw_free(Z1f);
        fftw_destroy_plan(fplan1); fftw_destroy_plan(fplan2); fftw_destroy_plan(iplan);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
