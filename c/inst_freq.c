//Gets instantaneous frequency (IF) for each vector in X along dim.
//The output Y is real-valued and has the same size as X.
//This is the phase angle in (-pi pi], equal to atan2(xi,xr),
//and the IF is the diff of this.
//Could just as easily use a centered diff.

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int inst_freq_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft)
{
    if (dim>3u) { fprintf(stderr,"error in inst_freq_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (nfft<L) { fprintf(stderr,"error in inst_freq_s: nfft must be >= L (vec length)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = 0.0f; }
    }
    else
    {
        const float sc = 2.0f / (float)nfft;
        const size_t isodd = nfft%2u;

        //Initialize fftwf
        float *X1, *Y1, *Z1;
        X1 = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Z1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan fplan = fftwf_plan_dft_r2c_1d((int)nfft,X1,(fftwf_complex *)Y1,FFTW_ESTIMATE);
        if (!fplan) { fprintf(stderr,"error in inst_freq_s: problem creating fftw plan"); return 1; }
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1,(fftwf_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in inst_freq_s: problem creating fftw plan"); return 1; }
        for (size_t n=L; n<nfft; ++n) { X1[n] = 0.0f; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1[n] = 0.0f; }
    
        if (L==N)
        {
            for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            //for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= L; //X1 -= nfft;
            fftwf_execute(fplan);
            *Y1++ /= (float)nfft; ++Y1;
            for (size_t n=nfft+isodd-2u; n>0u; --n, ++Y1) { *Y1 *= sc; }
            if (!isodd) { *Y1++ /= (float)nfft; }
            for (size_t l=nfft+1u; l<2u*nfft; ++l, ++Y1) { *Y1 = 0.0f; }
            Y1 -= 2u*nfft;
            fftwf_execute(iplan);
            for (size_t l=L; l>0u; --l, Z1+=2u, ++Y) { *Y = atan2f(*(Z1+1),*Z1); }
            for (size_t l=1u; l<L; ++l) { --Y; *Y -= *(Y-1); }
            --Y; *Y = *(Y+1);
            Z1 -= 2u*L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Z1-=2u*L)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    //for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
                    X1 -= L; //X1 -= nfft;
                    fftwf_execute(fplan);
                    *Y1++ /= (float)nfft; ++Y1;
                    for (size_t n=nfft+isodd-2u; n>0u; --n, ++Y1) { *Y1 *= sc; }
                    if (!isodd) { *Y1++ /= (float)nfft; }
                    Y1 -= nfft + 1u;
                    fftwf_execute(iplan);
                    for (size_t l=L; l>0u; --l, Z1+=2u, ++Y) { *Y = atan2f(*(Z1+1),*Z1); }
                    for (size_t l=1u; l<L; ++l) { --Y; *Y -= *(Y-1); }
                    --Y; *Y = *(Y+1); Y += L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Z1-=2u*L, ++Y)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                        //for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
                        X1 -= L; //X1 -= nfft;
                        fftwf_execute(fplan);
                        *Y1++ /= (float)nfft; ++Y1;
                        for (size_t n=nfft+isodd-2u; n>0u; --n, ++Y1) { *Y1 *= sc; }
                        if (!isodd) { *Y1++ /= (float)nfft; }
                        Y1 -= nfft + 1u;
                        fftwf_execute(iplan);
                        for (size_t l=L; l>0u; --l, Z1+=2u, Y+=K) { *Y = atan2f(*(Z1+1),*Z1); }
                        for (size_t l=1u; l<L; ++l) { Y-=K; *Y -= *(Y-K); }
                        Y-=K; *Y = *(Y+K);
                    }
                }
            }
        }
        fftwf_free(X1); fftwf_free(Y1); fftwf_free(Z1);
        fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan);
    }
    return 0;
}


int inst_freq_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft)
{
    if (dim>3u) { fprintf(stderr,"error in inst_freq_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (nfft<L) { fprintf(stderr,"error in inst_freq_d: nfft must be >= L (vec length)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; *++Y = 0.0; }
    }
    else
    {
        const double sc = 2.0 / (double)nfft;
        const size_t isodd = nfft%2u;

        //Initialize fftw
        double *X1, *Y1, *Z1;
        X1 = (double *)fftw_malloc(nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Z1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan fplan = fftw_plan_dft_r2c_1d((int)nfft,X1,(fftw_complex *)Y1,FFTW_ESTIMATE);
        if (!fplan) { fprintf(stderr,"error in inst_freq_d: problem creating fftw plan"); return 1; }
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1,(fftw_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in inst_freq_d: problem creating fftw plan"); return 1; }
        for (size_t n=L; n<nfft; ++n) { X1[n] = 0.0; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1[n] = 0.0; }
    
        if (L==N)
        {
            for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= L;
            fftw_execute(fplan);
            *Y1++ /= (double)nfft; ++Y1;
            for (size_t n=nfft+isodd-2u; n>0u; --n, ++Y1) { *Y1 *= sc; }
            if (!isodd) { *Y1++ /= (double)nfft; }
            Y1 -= nfft + 1u;
            fftw_execute(iplan);
            for (size_t l=L; l>0u; --l, Z1+=2u, ++Y) { *Y = atan2(*(Z1+1),*Z1); }
            for (size_t l=1u; l<L; ++l) { --Y; *Y -= *(Y-1); }
            --Y; *Y = *(Y+1);
            Z1 -= 2u*L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Z1-=2u*L)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(fplan);
                    *Y1++ /= (double)nfft; ++Y1;
                    for (size_t n=nfft+isodd-2u; n>0u; --n, ++Y1) { *Y1 *= sc; }
                    if (!isodd) { *Y1++ /= (double)nfft; }
                    Y1 -= nfft + 1u;
                    fftw_execute(iplan);
                    for (size_t l=L; l>0u; --l, Z1+=2u, ++Y) { *Y = atan2(*(Z1+1),*Z1); }
                    for (size_t l=1u; l<L; ++l) { --Y; *Y -= *(Y-1); }
                    --Y; *Y = *(Y+1);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Z1-=2u*L, ++Y)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        fftw_execute(fplan);
                        *Y1++ /= (double)nfft; ++Y1;
                        for (size_t n=nfft+isodd-2u; n>0u; --n, ++Y1) { *Y1 *= sc; }
                        if (!isodd) { *Y1++ /= (double)nfft; }
                        Y1 -= nfft + 1u;
                        fftw_execute(iplan);
                        for (size_t l=L; l>0u; --l, Z1+=2u, Y+=K) { *Y = atan2(*(Z1+1),*Z1); }
                        for (size_t l=1u; l<L; ++l) { Y-=K; *Y -= *(Y-K); }
                        Y-=K; *Y = *(Y+K);
                    }
                }
            }
        }
        fftw_free(X1); fftw_free(Y1); fftw_free(Z1);
        fftw_destroy_plan(fplan); fftw_destroy_plan(iplan);
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
