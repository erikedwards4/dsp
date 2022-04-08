//1D convolution of each vector in X1 by X2.
//Each vector in X1 has length L1. X2 has length L2 (kernel_size).
//This gives identical output to conv1d, but uses FFT internally.
//See conv1d for more info.

//I could not find any condition where this is faster than conv1d.
//conv1d is also more numerically accurate.

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int conv1d_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_fft_s: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_fft_s: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_fft_s: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_fft_s: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_fft_s: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_fft_s: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L1+2u*pad<=dil*(L2-1u)) { fprintf(stderr,"error in conv1d_fft_s: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    //Set Ly (W)
    const int N1 = (int)(L1+2u*pad);            //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const size_t Ly = 1u + (size_t)(N1-N2)/str;

    if (Ly==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while ((int)nfft<N1+N2) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;
        //if (nfft%2u) { fprintf(stderr,"error in conv1d_fft_s: nfft must be even"); return 1; }

        //Set scale
        const float sc = 2.0f / (float)nfft;

        //Initialize FFTs
        float *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (float *)fftwf_malloc(nfft*sizeof(float));
        X2f = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan fplan1 = fftwf_plan_dft_r2c_1d((int)nfft,X1f,(fftwf_complex *)Y1f,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in conv1d_fft_s: problem creating fftw plan"); return 1; }
        fftwf_plan fplan2 = fftwf_plan_dft_r2c_1d((int)nfft,X2f,(fftwf_complex *)Y2f,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in conv1d_fft_s: problem creating fftw plan"); return 1; }
        for (size_t n=nfft; n>0u; --n, ++X1f) { *X1f = 0.0f; }
        X1f -= nfft;

        //Initialize IFFT
        float *Z1f;
        Z1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1f,(fftwf_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv1d_fft_s: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++Y1f) { *Y1f = 0.0f; }
        Y1f -= 2u*nfft;

        //Get Y2f (scaled FFT of X2f with dilation holes)
        for (size_t n=L2; n>0u; --n, ++X2, ++X2f)
        {
            *X2f = *X2;
            for (size_t i=dil; i>1u; --i) { *++X2f = 0.0f; }
        }
        for (size_t n=L2*dil; n<nfft; ++n, ++X2f) { *X2f = 0.0f; }
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
            X1f += pad;
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
            X1f -= pad + L1;

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
            Z1f += 2u*(size_t)(N2-1);
            if ((int)pad>=N2)
            {
                const size_t nz = 1u + (pad-(size_t)N2)/str;
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0f; }
                Z1f += 2u*(1u + pad - (size_t)N2);
                for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; }
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0f; }
                Z1f -= 2u*(str*(Ly-2u*nz) + pad);
            }
            else
            {
                for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; }
                Z1f -= 2u*((size_t)N2-1u+Ly*str);
            }

            //clock_gettime(CLOCK_REALTIME,&toc);
            //fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=K*L1-1u, Y-=K*Ly-1u)
                {
                    //Copy input
                    X1f += pad;
                    for (size_t l=L1; l>0u; --l, X1+=K, ++X1f) { *X1f = *X1; }
                    X1f -= pad + L1;

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
                    Z1f += 2u*(size_t)(N2-1);
                    if ((int)pad>=N2)
                    {
                        const size_t nz = 1u + (pad-(size_t)N2)/str;
                        for (size_t l=nz; l>0u; --l, Y+=K) { *Y = 0.0f; }
                        Z1f += 2u*(1u + pad - (size_t)N2);
                        for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, Y+=K) { *Y = *Z1f; }
                        for (size_t l=nz; l>0u; --l, Y+=K) { *Y = 0.0f; }
                        Z1f -= 2u*(str*(Ly-2u*nz) + pad);
                    }
                    else
                    {
                        for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, Y+=K) { *Y = *Z1f; }
                        Z1f -= 2u*((size_t)N2-1u+Ly*str);
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


int conv1d_fft_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_fft_d: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_fft_d: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_fft_d: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_fft_d: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_fft_d: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_fft_d: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L1+2u*pad<=dil*(L2-1u)) { fprintf(stderr,"error in conv1d_fft_d: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    //Set Ly (W)
    const int N1 = (int)(L1+2u*pad);            //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const size_t Ly = 1u + (size_t)(N1-N2)/str;

    if (Ly==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while ((int)nfft<N1+N2) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;
        //if (nfft%2u) { fprintf(stderr,"error in conv1d_fft_s: nfft must be even"); return 1; }

        //Set scale
        const double sc = 2.0 / (double)nfft;

        //Initialize FFTs
        double *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (double *)fftw_malloc(nfft*sizeof(double));
        X2f = (double *)fftw_malloc(nfft*sizeof(double));
        Y1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan fplan1 = fftw_plan_dft_r2c_1d((int)nfft,X1f,(fftw_complex *)Y1f,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in conv1d_fft_d: problem creating fftw plan"); return 1; }
        fftw_plan fplan2 = fftw_plan_dft_r2c_1d((int)nfft,X2f,(fftw_complex *)Y2f,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in conv1d_fft_d: problem creating fftw plan"); return 1; }
        for (size_t n=nfft; n>0u; --n, ++X1f) { *X1f = 0.0; }
        X1f -= nfft;

        //Initialize IFFT
        double *Z1f;
        Z1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1f,(fftw_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv1d_fft_d: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++Y1f) { *Y1f = 0.0; }
        Y1f -= 2u*nfft;

        //Get Y2f (scaled FFT of X2f with dilation holes)
        for (size_t n=L2; n>0u; --n, ++X2, ++X2f)
        {
            *X2f = *X2;
            for (size_t i=dil; i>1u; --i) { *++X2f = 0.0; }
        }
        for (size_t n=L2*dil; n<nfft; ++n, ++X2f) { *X2f = 0.0; }
        X2f -= nfft;
        fftw_execute(fplan2);
        *Y2f++ /= (double)nfft; ++Y2f;
        for (size_t l=nfft-2u; l>0u; --l, ++Y2f) { *Y2f *= sc; }
        *Y2f /= (double)nfft;
        Y2f -= nfft;

        if (L1==N)
        {
            //Copy input
            X1f += pad;
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
            X1f -= pad + L1;

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
            Z1f += 2u*(size_t)(N2-1);
            if ((int)pad>=N2)
            {
                const size_t nz = 1u + (pad-(size_t)N2)/str;
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0; }
                Z1f += 2u*(1u + pad - (size_t)N2);
                for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; }
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0; }
                Z1f -= 2u*(str*(Ly-2u*nz) + pad);
            }
            else
            {
                for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; }
                Z1f -= 2u*((size_t)N2-1u+Ly*str);
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=K*L1-1u, Y-=K*Ly-1u)
                {
                    //Copy input
                    X1f += pad;
                    for (size_t l=L1; l>0u; --l, X1+=K, ++X1f) { *X1f = *X1; }
                    X1f -= pad + L1;

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
                    Z1f += 2u*(size_t)(N2-1);
                    if ((int)pad>=N2)
                    {
                        const size_t nz = 1u + (pad-(size_t)N2)/str;
                        for (size_t l=nz; l>0u; --l, Y+=K) { *Y = 0.0; }
                        Z1f += 2u*(1u + pad - (size_t)N2);
                        for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, Y+=K) { *Y = *Z1f; }
                        for (size_t l=nz; l>0u; --l, Y+=K) { *Y = 0.0; }
                        Z1f -= 2u*(str*(Ly-2u*nz) + pad);
                    }
                    else
                    {
                        for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, Y+=K) { *Y = *Z1f; }
                        Z1f -= 2u*((size_t)N2-1u+Ly*str);
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


int conv1d_fft_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_fft_c: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_fft_c: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_fft_c: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_fft_c: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_fft_c: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_fft_c: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L1+2u*pad<=dil*(L2-1u)) { fprintf(stderr,"error in conv1d_fft_c: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    //Set Ly (W)
    const int N1 = (int)(L1+2u*pad);            //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const size_t Ly = 1u + (size_t)(N1-N2)/str;

    if (Ly==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while ((int)nfft<N1+N2) { nfft *= 2u; }
        //if (nfft%2u) { fprintf(stderr,"error in conv1d_fft_c: nfft must be even"); return 1; }

        //Initialize FFTs
        float *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        X2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan fplan1 = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1f,(fftwf_complex *)Y1f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in conv1d_fft_c: problem creating fftw plan"); return 1; }
        fftwf_plan fplan2 = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X2f,(fftwf_complex *)Y2f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in conv1d_fft_c: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1f) { *X1f = 0.0f; }
        X1f -= 2u*nfft;

        //Initialize IFFT
        float *Z1f;
        Z1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1f,(fftwf_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv1d_fft_c: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++Y1f) { *Y1f = 0.0f; }
        Y1f -= 2u*nfft;

        //Get Y2f (scaled FFT of X2f with dilation holes)
        for (size_t n=L2; n>0u; --n)
        {
            *X2f++ = *X2++; *X2f++ = *X2++;
            for (size_t i=dil; i>1u; --i) { *X2f++ = 0.0f; *X2f++ = 0.0f; }
        }
        for (size_t n=2u*L2*dil; n<2u*nfft; ++n, ++X2f) { *X2f = 0.0f; }
        X2f -= 2u*nfft;
        fftwf_execute(fplan2);
        for (size_t l=2u*nfft; l>0u; --l, ++Y2f) { *Y2f /= (float)nfft; }
        Y2f -= 2u*nfft;

        if (L1==N)
        {
            //Copy input
            X1f += 2u*pad;
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; *++X1f = *++X1; }
            X1f -= 2u*(pad+L1);

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
            Z1f += 2u*(size_t)(N2-1);
            if ((int)pad>=N2)
            {
                const size_t nz = 1u + (pad-(size_t)N2)/str;
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                Z1f += 2u*(1u + pad - (size_t)N2);
                for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; *++Y = *(Z1f+1); }
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0f; *++Y = 0.0f; }
                Z1f -= 2u*(str*(Ly-2u*nz) + pad);
            }
            else
            {
                for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; *++Y = *(Z1f+1); }
                Z1f -= 2u*((size_t)N2-1u+Ly*str);
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=2u*K*L1-2u, Y-=2u*K*Ly-2u)
                {
                    //Copy input
                    X1f += 2u*pad;
                    for (size_t l=L1; l>0u; --l, X1+=2u*K, ++X1f) { *X1f = *X1; *++X1f = *(X1+1); }
                    X1f -= 2u*(pad+L1);

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
                    Z1f += 2u*(size_t)(N2-1);
                    if ((int)pad>=N2)
                    {
                        const size_t nz = 1u + (pad-(size_t)N2)/str;
                        for (size_t l=nz; l>0u; --l, Y+=2u*K) { *Y = 0.0f; *(Y+1) = 0.0f; }
                        Z1f += 2u*(1u + pad - (size_t)N2);
                        for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *(Z1f+1); }
                        for (size_t l=nz; l>0u; --l, Y+=2u*K) { *Y = 0.0f; *(Y+1) = 0.0f; }
                        Z1f -= 2u*(str*(Ly-2u*nz) + pad);
                    }
                    else
                    {
                        for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *(Z1f+1); }
                        Z1f -= 2u*((size_t)N2-1u+Ly*str);
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


int conv1d_fft_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t pad, const size_t str, const size_t dil, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv1d_fft_z: dim must be in [0 3]\n"); return 1; }
    if (str<1u) { fprintf(stderr,"error in conv1d_fft_z: str (stride) must be positive\n"); return 1; }
    if (dil<1u) { fprintf(stderr,"error in conv1d_fft_z: dil (dilation) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv1d_fft_z: L2 (length of X2) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N<1u) { fprintf(stderr,"error in conv1d_fft_z: N (total length of X1) must be positive\n"); return 1; }
    if (L1<1u) { fprintf(stderr,"error in conv1d_fft_z: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L1+2u*pad<=dil*(L2-1u)) { fprintf(stderr,"error in conv1d_fft_z: L1+2*pad must be > dil*(L2-1)\n"); return 1; }

    //Set Ly (W)
    const int N1 = (int)(L1+2u*pad);            //full length of vecs in X1 including padding
    const int N2 = (int)(dil*(L2-1u)) + 1;      //full length of vec in X2 including dilation
    const size_t Ly = 1u + (size_t)(N1-N2)/str;

    if (Ly==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while ((int)nfft<N1+N2) { nfft *= 2u; }
        //if (nfft%2u) { fprintf(stderr,"error in conv1d_fft_c: nfft must be even"); return 1; }

        //Initialize FFTs
        double *X1f, *X2f, *Y1f, *Y2f, y1r, y1i, y2r, y2i;
        X1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        X2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan fplan1 = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1f,(fftw_complex *)Y1f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan1) { fprintf(stderr,"error in conv1d_fft_z: problem creating fftw plan"); return 1; }
        fftw_plan fplan2 = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X2f,(fftw_complex *)Y2f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!fplan2) { fprintf(stderr,"error in conv1d_fft_z: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1f) { *X1f = 0.0; }
        X1f -= 2u*nfft;

        //Initialize IFFT
        double *Z1f;
        Z1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1f,(fftw_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv1d_fft_z: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++Y1f) { *Y1f = 0.0; }
        Y1f -= 2u*nfft;

        //Get Y2f (scaled FFT of X2f with dilation holes)
        for (size_t n=L2; n>0u; --n)
        {
            *X2f++ = *X2++; *X2f++ = *X2++;
            for (size_t i=dil; i>1u; --i) { *X2f++ = 0.0; *X2f++ = 0.0; }
        }
        for (size_t n=2u*L2*dil; n<2u*nfft; ++n, ++X2f) { *X2f = 0.0; }
        X2f -= 2u*nfft;
        fftw_execute(fplan2);
        for (size_t l=2u*nfft; l>0u; --l, ++Y2f) { *Y2f /= (double)nfft; }
        Y2f -= 2u*nfft;

        if (L1==N)
        {
            //Copy input
            X1f += 2u*pad;
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; *++X1f = *++X1; }
            X1f -= 2u*(pad+L1);

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
            Z1f += 2u*(size_t)(N2-1);
            if ((int)pad>=N2)
            {
                const size_t nz = 1u + (pad-(size_t)N2)/str;
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0; *++Y = 0.0; }
                Z1f += 2u*(1u + pad - (size_t)N2);
                for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; *++Y = *(Z1f+1); }
                for (size_t l=nz; l>0u; --l, ++Y) { *Y = 0.0; *++Y = 0.0; }
                Z1f -= 2u*(str*(Ly-2u*nz) + pad);
            }
            else
            {
                for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, ++Y) { *Y = *Z1f; *++Y = *(Z1f+1); }
                Z1f -= 2u*((size_t)N2-1u+Ly*str);
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            for (size_t g=G; g>0u; --g, X1+=2u*B*(L1-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=B; b>0u; --b, X1-=2u*K*L1-2u, Y-=2u*K*Ly-2u)
                {
                    //Copy input
                    X1f += 2u*pad;
                    for (size_t l=L1; l>0u; --l, X1+=2u*K, ++X1f) { *X1f = *X1; *++X1f = *(X1+1); }
                    X1f -= 2u*(pad+L1);

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
                    Z1f += 2u*(size_t)(N2-1);
                    if ((int)pad>=N2)
                    {
                        const size_t nz = 1u + (pad-(size_t)N2)/str;
                        for (size_t l=nz; l>0u; --l, Y+=2u*K) { *Y = 0.0; *(Y+1) = 0.0; }
                        Z1f += 2u*(1u + pad - (size_t)N2);
                        for (size_t l=Ly-2u*nz; l>0u; --l, Z1f+=2u*str, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *(Z1f+1); }
                        for (size_t l=nz; l>0u; --l, Y+=2u*K) { *Y = 0.0; *(Y+1) = 0.0; }
                        Z1f -= 2u*(str*(Ly-2u*nz) + pad);
                    }
                    else
                    {
                        for (size_t l=Ly; l>0u; --l, Z1f+=2u*str, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *(Z1f+1); }
                        Z1f -= 2u*((size_t)N2-1u+Ly*str);
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
