//This computes "the" IDCT (inverse discrete cosine transformation) along dim of matrix X.
//This uses fftw3 to compute the IDCT-II.

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int idct_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc);
int idct_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc);
int idct_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc);
int idct_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc);


int idct_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_s: ndct must be >= L (vec length)\n"); return 1; }
    const float xsc = (sc) ? (float)(2.0*M_SQRT1_2) : 1.0f;
    const float ysc = (sc) ? (float)(1.0/sqrt(2*ndct)) : 1.0f/(2*ndct);

    //Initialize fftwf
    float *X1, *Y1;
    X1 = (float *)fftwf_malloc(ndct*sizeof(float));
    Y1 = (float *)fftwf_malloc(ndct*sizeof(float));
    fftwf_plan plan = fftwf_plan_r2r_1d((int)ndct,X1,Y1,FFTW_REDFT01,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in idct_s: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndct==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *X1++ = *X++ * xsc;
        for (size_t l=1; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0f; }
        X1 -= ndct;
        fftwf_execute(plan);
        for (size_t l=0; l<ndct; ++l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
        Y1 -= ndct;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            *X1++ = *X++ * xsc;
            for (size_t l=1; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= ndct;
            fftwf_execute(plan);
            for (size_t l=0; l<ndct; ++l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
            Y1 -= ndct;
            for (size_t v=1; v<V; ++v, Y1-=ndct)
            {
                *X1++ = *X++ * xsc;
                for (size_t l=1; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                fftwf_execute(plan);
                for (size_t l=0; l<ndct; ++l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
            }
        }
        else
        {
            X1 += L;
            for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= ndct;
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(ndct-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y1-=ndct, Y-=K*ndct-1)
                {
                    *X1++ = *X * xsc; X += K;
                    for (size_t l=1; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftwf_execute(plan);
                    for (size_t l=0; l<ndct; ++l, ++Y1, Y+=K) { *Y = *Y1 * ysc; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int idct_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_d: ndct must be >= L (vec length)\n"); return 1; }
    const double xsc = (sc) ? 2.0*M_SQRT1_2 : 1.0;
    const double ysc = (sc) ? 1.0/sqrt(2*ndct) : 1.0/(2*ndct);

    //Initialize fftw
    double *X1, *Y1;
    X1 = (double *)fftw_malloc(ndct*sizeof(double));
    Y1 = (double *)fftw_malloc(ndct*sizeof(double));
    fftw_plan plan = fftw_plan_r2r_1d((int)ndct,X1,Y1,FFTW_REDFT01,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in idct_d: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndct==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *X1++ = *X++ * xsc;
        for (size_t l=1; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0; }
        X1 -= ndct;
        fftw_execute(plan);
        for (size_t l=0; l<ndct; ++l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
        Y1 -= ndct;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            *X1++ = *X++ * xsc;
            for (size_t l=1; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0; }
            X1 -= ndct;
            fftw_execute(plan);
            for (size_t l=0; l<ndct; ++l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
            Y1 -= ndct;
            for (size_t v=1; v<V; ++v, Y1-=ndct)
            {
                *X1++ = *X++ * xsc;
                for (size_t l=1; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                fftw_execute(plan);
                for (size_t l=0; l<ndct; ++l, ++Y1, ++Y) { *Y = *Y1 * ysc; }
            }
        }
        else
        {
            X1 += L;
            for (size_t l=L; l<ndct; ++l, ++X1) { *X1 = 0.0; }
            X1 -= ndct;
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(ndct-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y1-=ndct, Y-=K*ndct-1)
                {
                    *X1++ = *X * xsc; X += K;
                    for (size_t l=1; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(plan);
                    for (size_t l=0; l<ndct; ++l, ++Y1, Y+=K) { *Y = *Y1 * ysc; }
                }
            }
        }
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


int idct_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_c: ndct must be >= L (vec length)\n"); return 1; }
    const float xsc = (sc) ? (float)(2.0*M_SQRT1_2) : 1.0f;
    const float ysc = (sc) ? (float)(1.0/sqrt(2*ndct)) : 1.0f/(2*ndct);

    //Initialize fftwf
    float *X1r, *Y1r, *X1i, *Y1i;
    X1r = (float *)fftwf_malloc(ndct*sizeof(float));
    X1i = (float *)fftwf_malloc(ndct*sizeof(float));
    Y1r = (float *)fftwf_malloc(ndct*sizeof(float));
    Y1i = (float *)fftwf_malloc(ndct*sizeof(float));
    fftwf_plan rplan = fftwf_plan_r2r_1d((int)ndct,X1r,Y1r,FFTW_REDFT01,FFTW_ESTIMATE);
    if (!rplan) { fprintf(stderr,"error in idct_c: problem creating fftw plan"); return 1; }
    fftwf_plan iplan = fftwf_plan_r2r_1d((int)ndct,X1i,Y1i,FFTW_REDFT01,FFTW_ESTIMATE);
    if (!iplan) { fprintf(stderr,"error in idct_c: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndct==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
        for (size_t l=1; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
        for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
        X1r -= ndct; X1i -= ndct;
        fftwf_execute(rplan); fftwf_execute(iplan);
        for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
        Y1r -= ndct; Y1i -= ndct;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
            for (size_t l=1; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
            for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
            X1r -= ndct; X1i -= ndct;
            fftwf_execute(rplan); fftwf_execute(iplan);
            for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
            Y1r -= ndct; Y1i -= ndct;
            for (size_t v=1; v<V; ++v, Y1r-=ndct, Y1i-=ndct)
            {
                *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
                for (size_t l=1; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                X1r -= L; X1i -= L;
                fftwf_execute(rplan); fftwf_execute(iplan);
                for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
            }
        }
        else
        {
            X1r += L; X1i += L;
            for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
            X1r -= ndct; X1i -= ndct;
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(ndct-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y1r-=ndct, Y1i-=ndct, Y-=2*K*ndct-2)
                {
                    *X1r++ = *X * xsc; *X1i++ = *++X * xsc; X += 2*K-1;
                    for (size_t l=1; l<L; ++l, X+=2*K-1, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftwf_execute(rplan); fftwf_execute(iplan);
                    for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, Y+=2*K-1) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(rplan); fftwf_free(X1r); fftwf_free(Y1r);
    fftwf_destroy_plan(iplan); fftwf_free(X1i); fftwf_free(Y1i);
    return 0;
}


int idct_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndct, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndct<L) { fprintf(stderr,"error in idct_z: ndct must be >= L (vec length)\n"); return 1; }
    const double xsc = (sc) ? 2.0*M_SQRT1_2 : 1.0;
    const double ysc = (sc) ? 1.0/sqrt(2*ndct) : 1.0/(2*ndct);

    //Initialize fftw
    double *X1r, *Y1r, *X1i, *Y1i;
    X1r = (double *)fftw_malloc(ndct*sizeof(double));
    X1i = (double *)fftw_malloc(ndct*sizeof(double));
    Y1r = (double *)fftw_malloc(ndct*sizeof(double));
    Y1i = (double *)fftw_malloc(ndct*sizeof(double));
    fftw_plan rplan = fftw_plan_r2r_1d((int)ndct,X1r,Y1r,FFTW_REDFT01,FFTW_ESTIMATE);
    if (!rplan) { fprintf(stderr,"error in idct_z: problem creating fftw plan"); return 1; }
    fftw_plan iplan = fftw_plan_r2r_1d((int)ndct,X1i,Y1i,FFTW_REDFT01,FFTW_ESTIMATE);
    if (!iplan) { fprintf(stderr,"error in idct_z: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndct==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
        for (size_t l=1; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
        for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
        X1r -= ndct; X1i -= ndct;
        fftw_execute(rplan); fftw_execute(iplan);
        for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
        Y1r -= ndct; Y1i -= ndct;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
            for (size_t l=1; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
            for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
            X1r -= ndct; X1i -= ndct;
            fftw_execute(rplan); fftw_execute(iplan);
            for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
            Y1r -= ndct; Y1i -= ndct;
            for (size_t v=1; v<V; ++v, Y1r-=ndct, Y1i-=ndct)
            {
                *X1r++ = *X++ * xsc; *X1i++ = *X++ * xsc;
                for (size_t l=1; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                X1r -= L; X1i -= L;
                fftw_execute(rplan); fftw_execute(iplan);
                for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
            }
        }
        else
        {
            X1r += L; X1i += L;
            for (size_t l=L; l<ndct; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
            X1r -= ndct; X1i -= ndct;
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(ndct-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y1r-=ndct, Y1i-=ndct, Y-=2*K*ndct-2)
                {
                    *X1r++ = *X * xsc; *X1i++ = *++X * xsc; X += 2*K-1;
                    for (size_t l=1; l<L; ++l, X+=2*K-1, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftw_execute(rplan); fftw_execute(iplan);
                    for (size_t l=0; l<ndct; ++l, ++Y1r, ++Y1i, Y+=2*K-1) { *Y = *Y1r*ysc; *++Y = *Y1i*ysc; }
                }
            }
        }
    }
    
    fftw_destroy_plan(rplan); fftw_free(X1r); fftw_free(Y1r);
    fftw_destroy_plan(iplan); fftw_free(X1i); fftw_free(Y1i);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
