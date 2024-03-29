//Does 1-D IFFT (inverse fast Fourier transform) of each vector in X along dim.
//The input X is complex-valued.
//The output Y has the same size as X, except along dim, where Y has length nfft.
//Y is real-valued for ifft_fftw_s and ifft_fftw_d,
//and complex-valued for ifft_fftw_c and ifft_fftw_z.
//In the former case, X has only nonnegative freqs, so Lx = nfft/2 + 1.

//If sc, then scales Y by 2.0*sqrt(2.0)/nfft, so that invertible with ifft.

//This uses the FFTW library (see also ifft.rad2.c, ifft.ffts.c).

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int ifft_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_fftw_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float s = sc ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_fftw_s: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize fftwf
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in ifft_fftw_s: problem creating fftw plan"); return 1; }
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
            X -= 2u + 2u*(1u-nfft%2u);
            for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
            X1 -= 2u*nfft;
            fftwf_execute(iplan);
            for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y1-=2u*nfft)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X -= 2u + 2u*(1u-nfft%2u);
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                    X1 -= 2u*nfft;
                    fftwf_execute(iplan);
                    for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y1-=2u*nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1u); }
                        X -= K*(2u + 2u*(1u-nfft%2u));
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                        X1 -= 2u*nfft;
                        fftwf_execute(iplan);
                        for (size_t l=nfft; l>0u; --l, Y1+=2u, Y+=K) { *Y = *Y1 * s; }
                    }
                }
            }
        }
        fftwf_destroy_plan(iplan); fftwf_free(X1); fftwf_free(Y1);
    }

    return 0;
}


int ifft_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_fftw_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double s = sc ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_fftw_d: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize fftw
        double *X1, *Y1;
        X1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in ifft_fftw_d: problem creating fftw plan"); return 1; }
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
            X -= 2u + 2u*(1u-nfft%2u);
            for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
            X1 -= 2u*nfft;
            fftw_execute(iplan);
            for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y1-=2u*nfft)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
                    X -= 2u + 2u*(1u-nfft%2u);
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                    X1 -= 2u*nfft;
                    fftw_execute(iplan);
                    for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y1-=2u*nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1u); }
                        X -= K*(2u + 2u*(1u-nfft%2u));
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                        X1 -= 2u*nfft;
                        fftw_execute(iplan);
                        for (size_t l=nfft; l>0u; --l, Y1+=2u, Y+=K) { *Y = *Y1 * s; }
                    }
                }
            }
        }
        fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
    }

    return 0;
}


int ifft_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_fftw_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float s = sc ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_fftw_c: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize fftwf
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in ifft_fftw_c: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*nfft;
            fftwf_execute(iplan);
            for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= 2u*nfft;
                fftwf_execute(iplan);
                for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2u*nfft;
                for (size_t v=1u; v<V; ++v, Y1-=2u*nfft)
                {
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*nfft;
                    fftwf_execute(iplan);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1u); }
                        X1 -= 2u*nfft;
                        fftwf_execute(iplan);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K-1u) { *Y = *Y1 * s; *++Y = *++Y1 * s; }
                    }
                }
            }
        }
        fftwf_destroy_plan(iplan); fftwf_free(X1); fftwf_free(Y1);
    }

    return 0;
}


int ifft_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_fftw_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double s = sc ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_fftw_z: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize fftwf
        double *X1, *Y1;
        X1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in ifft_fftw_z: problem creating fftw plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*nfft;
            fftw_execute(iplan);
            for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                X1 -= 2u*nfft;
                fftw_execute(iplan);
                for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2u*nfft;
                for (size_t v=1u; v<V; ++v, Y1-=2u*nfft)
                {
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*nfft;
                    fftw_execute(iplan);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1u); }
                        X1 -= 2u*nfft;
                        fftw_execute(iplan);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K-1u) { *Y = *Y1 * s; *++Y = *++Y1 * s; }
                    }
                }
            }
        }
        fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
