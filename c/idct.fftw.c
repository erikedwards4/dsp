//This computes "the" IDCT (inverse discrete cosine transformation) along dim of matrix X.
//This uses fftw3 to compute the IDCT-II.

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include "codee_dsp.h"

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int idct_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_fftw_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_fftw_s: ndct must be >= L (vec length)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u && ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const float xsc = (sc) ? (float)(2.0*M_SQRT1_2) : 1.0f;
        const float ysc = (sc) ? 1.0f/sqrtf((float)(2u*ndct)) : 1.0f/(float)(2u*ndct);

        //Initialize fftwf
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(ndct*sizeof(float));
        Y1 = (float *)fftwf_malloc(ndct*sizeof(float));
        fftwf_plan plan = fftwf_plan_r2r_1d((int)ndct,X1,Y1,FFTW_REDFT01,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in idct_fftw_s: problem creating fftw plan"); return 1; }
    
        if (L==N)
        {
            *X1++ = *X++ * xsc;
            for (size_t l=1u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= ndct;
            fftwf_execute(plan);
            for (size_t l=ndct; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
            Y1 -= ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                *X1++ = *X++ * xsc;
                for (size_t l=1u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0f; }
                X1 -= ndct;
                fftwf_execute(plan);
                for (size_t l=ndct; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
                Y1 -= ndct;
                for (size_t v=1u; v<V; ++v, Y1-=ndct)
                {
                    *X1++ = *X++ * xsc;
                    for (size_t l=1u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftwf_execute(plan);
                    for (size_t l=ndct; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
                }
            }
            else
            {
                X1 += L;
                for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0f; }
                X1 -= ndct;
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y1-=ndct, Y-=K*ndct-1u)
                    {
                        *X1++ = *X * xsc; X += K;
                        for (size_t l=1u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        fftwf_execute(plan);
                        for (size_t l=ndct; l>0u; --l, ++Y1, Y+=K) { *Y = *Y1 * ysc; }
                    }
                }
            }
        }
        fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    }

    return 0;
}


int idct_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_fftw_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_fftw_d: ndct must be >= L (vec length)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u && ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const double xsc = (sc) ? 2.0*M_SQRT1_2 : 1.0;
        const double ysc = (sc) ? 1.0/sqrt((double)(2u*ndct)) : 1.0/(double)(2u*ndct);

        //Initialize fftw
        double *X1, *Y1;
        X1 = (double *)fftw_malloc(ndct*sizeof(double));
        Y1 = (double *)fftw_malloc(ndct*sizeof(double));
        fftw_plan plan = fftw_plan_r2r_1d((int)ndct,X1,Y1,FFTW_REDFT01,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in idct_fftw_d: problem creating fftw plan"); return 1; }
    
        if (L==N)
        {
            *X1++ = *X++ * xsc;
            for (size_t l=1u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0; }
            X1 -= ndct;
            fftw_execute(plan);
            for (size_t l=ndct; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
            Y1 -= ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                *X1++ = *X++ * xsc;
                for (size_t l=1u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0; }
                X1 -= ndct;
                fftw_execute(plan);
                for (size_t l=ndct; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
                Y1 -= ndct;
                for (size_t v=1u; v<V; ++v, Y1-=ndct)
                {
                    *X1++ = *X++ * xsc;
                    for (size_t l=1u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(plan);
                    for (size_t l=ndct; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
                }
            }
            else
            {
                X1 += L;
                for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0; }
                X1 -= ndct;
                for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*L-1u, Y1-=ndct, Y-=K*ndct-1u)
                    {
                        *X1++ = *X * xsc; X += K;
                        for (size_t l=1u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        fftw_execute(plan);
                        for (size_t l=ndct; l>0u; --l, ++Y1, Y+=K) { *Y = *Y1 * ysc; }
                    }
                }
            }
        }
        fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    }

    return 0;
}


int idct_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_fftw_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_fftw_c: ndct must be >= L (vec length)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u && ndct==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const float xsc = (sc) ? (float)(2.0*M_SQRT1_2) : 1.0f;
        const float ysc = (sc) ? 1.0f/sqrtf((float)(2u*ndct)) : 1.0f/(float)(2u*ndct);

        //Initialize fftwf
        float *X1r, *Y1r, *X1i, *Y1i;
        X1r = (float *)fftwf_malloc(ndct*sizeof(float));
        X1i = (float *)fftwf_malloc(ndct*sizeof(float));
        Y1r = (float *)fftwf_malloc(ndct*sizeof(float));
        Y1i = (float *)fftwf_malloc(ndct*sizeof(float));
        fftwf_plan rplan = fftwf_plan_r2r_1d((int)ndct,X1r,Y1r,FFTW_REDFT01,FFTW_ESTIMATE);
        if (!rplan) { fprintf(stderr,"error in idct_fftw_c: problem creating fftw plan"); return 1; }
        fftwf_plan iplan = fftwf_plan_r2r_1d((int)ndct,X1i,Y1i,FFTW_REDFT01,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in idct_fftw_c: problem creating fftw plan"); return 1; }
    
        if (L==N)
        {
            *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
            for (size_t l=1u; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
            for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
            X1r -= ndct; X1i -= ndct;
            fftwf_execute(rplan); fftwf_execute(iplan);
            for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
            Y1r -= ndct; Y1i -= ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
                for (size_t l=1u; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
                X1r -= ndct; X1i -= ndct;
                fftwf_execute(rplan); fftwf_execute(iplan);
                for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                Y1r -= ndct; Y1i -= ndct;
                for (size_t v=1u; v<V; ++v, Y1r-=ndct, Y1i-=ndct)
                {
                    *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
                    for (size_t l=1u; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftwf_execute(rplan); fftwf_execute(iplan);
                    for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                }
            }
            else
            {
                X1r += L; X1i += L;
                for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
                X1r -= ndct; X1i -= ndct;
                for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y1r-=ndct, Y1i-=ndct, Y-=2u*K*ndct-2u)
                    {
                        *X1r++ = *X * xsc; *X1i++ = *++X * xsc; X += 2u*K-1u;
                        for (size_t l=1u; l<L; ++l, X+=2u*K-1u, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                        X1r -= L; X1i -= L;
                        fftwf_execute(rplan); fftwf_execute(iplan);
                        for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, Y+=2u*K-1u) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                    }
                }
            }
        }
        fftwf_destroy_plan(rplan); fftwf_free(X1r); fftwf_free(Y1r);
        fftwf_destroy_plan(iplan); fftwf_free(X1i); fftwf_free(Y1i);
    }

    return 0;
}


int idct_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_fftw_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_fftw_z: ndct must be >= L (vec length)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u && ndct==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const double xsc = (sc) ? 2.0*M_SQRT1_2 : 1.0;
        const double ysc = (sc) ? 1.0/sqrt((double)(2u*ndct)) : 1.0/(double)(2u*ndct);

        //Initialize fftw
        double *X1r, *Y1r, *X1i, *Y1i;
        X1r = (double *)fftw_malloc(ndct*sizeof(double));
        X1i = (double *)fftw_malloc(ndct*sizeof(double));
        Y1r = (double *)fftw_malloc(ndct*sizeof(double));
        Y1i = (double *)fftw_malloc(ndct*sizeof(double));
        fftw_plan rplan = fftw_plan_r2r_1d((int)ndct,X1r,Y1r,FFTW_REDFT01,FFTW_ESTIMATE);
        if (!rplan) { fprintf(stderr,"error in idct_fftw_z: problem creating fftw plan"); return 1; }
        fftw_plan iplan = fftw_plan_r2r_1d((int)ndct,X1i,Y1i,FFTW_REDFT01,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in idct_fftw_z: problem creating fftw plan"); return 1; }
    
        if (L==N)
        {
            *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
            for (size_t l=1u; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
            for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
            X1r -= ndct; X1i -= ndct;
            fftw_execute(rplan); fftw_execute(iplan);
            for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
            Y1r -= ndct; Y1i -= ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
                for (size_t l=1u; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
                X1r -= ndct; X1i -= ndct;
                fftw_execute(rplan); fftw_execute(iplan);
                for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                Y1r -= ndct; Y1i -= ndct;
                for (size_t v=1u; v<V; ++v, Y1r-=ndct, Y1i-=ndct)
                {
                    *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
                    for (size_t l=1u; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftw_execute(rplan); fftw_execute(iplan);
                    for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                }
            }
            else
            {
                X1r += L; X1i += L;
                for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
                X1r -= ndct; X1i -= ndct;
                for (size_t g=G; g>0u; --g, X+=2u*B*(L-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*L-2u, Y1r-=ndct, Y1i-=ndct, Y-=2u*K*ndct-2u)
                    {
                        *X1r++ = *X * xsc; *X1i++ = *++X * xsc; X += 2u*K-1u;
                        for (size_t l=1u; l<L; ++l, X+=2u*K-1u, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                        X1r -= L; X1i -= L;
                        fftw_execute(rplan); fftw_execute(iplan);
                        for (size_t l=ndct; l>0u; --l, ++Y1r, ++Y1i, Y+=2u*K-1u) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                    }
                }
            }
        }
        fftw_destroy_plan(rplan); fftw_free(X1r); fftw_free(Y1r);
        fftw_destroy_plan(iplan); fftw_free(X1i); fftw_free(Y1i);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
