//This computes "the" DST (discrete sine transformation) along dim of matrix X.
//This uses fftw3 to compute the DST-I.

//For complex input X, the output Y is complex, and the DST is just the DST of the
//real and imaginary parts separately (following Octave convention).

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int dst_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc);
int dst_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc);
int dst_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc);
int dst_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc);


int dst_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndst<L) { fprintf(stderr,"error in dst_s: ndst must be >= L (vec length)\n"); return 1; }

    //Initialize fftwf
    float *X1, *Y1;
    X1 = (float *)fftwf_malloc(ndst*sizeof(float));
    Y1 = (float *)fftwf_malloc(ndst*sizeof(float));
    fftwf_plan plan = fftwf_plan_r2r_1d((int)ndst,X1,Y1,FFTW_RODFT00,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in dst_s: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndst==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=L; l<ndst; ++l, ++X1) { *X1 = 0.0f; }
        X1 -= ndst;
        fftwf_execute(plan);
        if (sc)
        {
            const float s = 0.5f; //(float)(1.0/sqrt(ndst))
            *Y++ = *Y1++ * 0.5f;
            for (size_t l=1; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        }
        else
        {
            for (size_t l=0; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1; }
        }
        Y1 -= ndst;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<ndst; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= ndst;
            fftwf_execute(plan);
            if (sc)
            {
                const float s = 0.5f, dcsc = 0.5f;
                *Y++ = *Y1++ * dcsc;
                for (size_t l=1; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= ndst;
                for (size_t v=1; v<V; ++v, Y1-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftwf_execute(plan);
                    *Y++ = *Y1++ * dcsc;
                    for (size_t l=1; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t l=0; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1; }
                Y1 -= ndst;
                for (size_t v=1; v<V; ++v, Y1-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftwf_execute(plan);
                    for (size_t l=0; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
            }
        }
        else
        {
            const float s = sc ? 0.5f : 1.0f;
            const float dcsc = sc ? 0.5f : 1.0f;
            X1 += L;
            for (size_t l=L; l<ndst; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= ndst;
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(ndst-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y1-=ndst, Y-=K*ndst-1)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftwf_execute(plan);
                    *Y = *Y1++ * dcsc; Y += K;
                    for (size_t l=1; l<ndst; ++l, ++Y1, Y+=K) { *Y = *Y1*s; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int dst_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndst<L) { fprintf(stderr,"error in dst_d: ndst must be >= L (vec length)\n"); return 1; }

    //Initialize fftw
    double *X1, *Y1;
    X1 = (double *)fftw_malloc(ndst*sizeof(double));
    Y1 = (double *)fftw_malloc(ndst*sizeof(double));
    fftw_plan plan = fftw_plan_r2r_1d((int)ndst,X1,Y1,FFTW_RODFT00,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in dst_d: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndst==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=L; l<ndst; ++l, ++X1) { *X1 = 0.0; }
        X1 -= ndst;
        fftw_execute(plan);
        if (sc)
        {
            const double s = 0.5;
            *Y++ = *Y1++ * 0.5;
            for (size_t l=1; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        }
        else
        {
            for (size_t l=0; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1; }
        }
        Y1 -= ndst;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<ndst; ++l, ++X1) { *X1 = 0.0; }
            X1 -= ndst;
            fftw_execute(plan);
            if (sc)
            {
                const double s = 0.5, dcsc = 0.5;
                *Y++ = *Y1++ * dcsc;
                for (size_t l=1; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= ndst;
                for (size_t v=1; v<V; ++v, Y1-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(plan);
                    *Y++ = *Y1++ * dcsc;
                    for (size_t l=1; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t l=0; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1; }
                Y1 -= ndst;
                for (size_t v=1; v<V; ++v, Y1-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(plan);
                    for (size_t l=0; l<ndst; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
            }
        }
        else
        {
            const double s = sc ? 0.5 : 1.0;
            const double dcsc = sc ? 0.5 : 1.0;
            X1 += L;
            for (size_t l=L; l<ndst; ++l, ++X1) { *X1 = 0.0; }
            X1 -= ndst;
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=B*(ndst-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Y1-=ndst, Y-=K*ndst-1)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(plan);
                    *Y = *Y1++ * dcsc; Y += K;
                    for (size_t l=1; l<ndst; ++l, ++Y1, Y+=K) { *Y = *Y1*s; }
                }
            }
        }
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


int dst_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndst<L) { fprintf(stderr,"error in dst_c: ndst must be >= L (vec length)\n"); return 1; }

    //Initialize fftwf
    float *X1r, *X1i, *Y1r, *Y1i;
    X1r = (float *)fftwf_malloc(ndst*sizeof(float));
    X1i = (float *)fftwf_malloc(ndst*sizeof(float));
    Y1r = (float *)fftwf_malloc(ndst*sizeof(float));
    Y1i = (float *)fftwf_malloc(ndst*sizeof(float));
    fftwf_plan rplan = fftwf_plan_r2r_1d((int)ndst,X1r,Y1r,FFTW_RODFT00,FFTW_ESTIMATE);
    if (!rplan) { fprintf(stderr,"error in dst_c: problem creating fftw plan"); return 1; }
    fftwf_plan iplan = fftwf_plan_r2r_1d((int)ndst,X1i,Y1i,FFTW_RODFT00,FFTW_ESTIMATE);
    if (!iplan) { fprintf(stderr,"error in dst_c: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndst==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
        for (size_t l=L; l<ndst; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
        X1r -= ndst; X1i -= ndst;
        fftwf_execute(rplan); fftwf_execute(iplan);
        if (sc)
        {
            const float s = 0.5f;
            *Y++ = *Y1r++ * 0.5f;
            *Y++ = *Y1i++ * 0.5f;
            for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r * s; *++Y = *Y1i * s; }
        }
        else
        {
            for (size_t l=0; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r; *++Y = *Y1i; }
        }
        Y1r -= ndst; Y1i -= ndst;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
            for (size_t l=L; l<ndst; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
            X1r -= ndst; X1i -= ndst;
            fftwf_execute(rplan); fftwf_execute(iplan);
            if (sc)
            {
                const float s = 0.5f, dcsc = 0.5f;
                *Y++ = *Y1r++ * dcsc; *Y++ = *Y1i++ * dcsc;
                for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r * s; *++Y = *Y1i * s; }
                Y1r -= ndst; Y1i -= ndst;
                for (size_t v=1; v<V; ++v, Y1r-=ndst, Y1i-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftwf_execute(rplan); fftwf_execute(iplan);
                    *Y++ = *Y1r++ * dcsc; *Y++ = *Y1i++ * dcsc;
                    for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r * s; *++Y = *Y1i * s; }
                }
            }
            else
            {
                for (size_t l=0; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r; *++Y = *Y1i; }
                Y1r -= ndst; Y1i -= ndst;
                for (size_t v=1; v<V; ++v, Y1r-=ndst, Y1i-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftwf_execute(rplan); fftwf_execute(iplan);
                    for (size_t l=0; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r; *++Y = *Y1i; }
                }
            }
        }
        else
        {
            const float s = sc ? 0.5f : 1.0f;
            const float dcsc = sc ? 0.5f : 1.0f;
            X1r += L; X1i += L;
            for (size_t l=L; l<ndst; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0f; }
            X1r -= ndst; X1i -= ndst;
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(ndst-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y1r-=ndst, Y1i-=ndst, Y-=2*K*ndst-2)
                {
                    for (size_t l=0; l<L; ++l, X+=2*K-1, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftwf_execute(rplan); fftwf_execute(iplan);
                    *Y = *Y1r++ * dcsc; *++Y = *Y1i++ * dcsc; Y += 2*K-1;
                    for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, Y+=2*K-1) { *Y = *Y1r*s; *++Y = *Y1i*s; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(rplan); fftwf_free(X1r); fftwf_free(Y1r);
    fftwf_destroy_plan(iplan); fftwf_free(X1i); fftwf_free(Y1i);
    return 0;
}


int dst_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t ndst, const char sc)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (ndst<L) { fprintf(stderr,"error in dst_z: ndst must be >= L (vec length)\n"); return 1; }

    //Initialize fftw
    double *X1r, *X1i, *Y1r, *Y1i;
    X1r = (double *)fftw_malloc(ndst*sizeof(double));
    X1i = (double *)fftw_malloc(ndst*sizeof(double));
    Y1r = (double *)fftw_malloc(ndst*sizeof(double));
    Y1i = (double *)fftw_malloc(ndst*sizeof(double));
    fftw_plan rplan = fftw_plan_r2r_1d((int)ndst,X1r,Y1r,FFTW_RODFT00,FFTW_ESTIMATE);
    if (!rplan) { fprintf(stderr,"error in dst_z: problem creating fftw plan"); return 1; }
    fftw_plan iplan = fftw_plan_r2r_1d((int)ndst,X1i,Y1i,FFTW_RODFT00,FFTW_ESTIMATE);
    if (!iplan) { fprintf(stderr,"error in dst_z: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1 && ndst==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
        for (size_t l=L; l<ndst; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
        X1r -= ndst; X1i -= ndst;
        fftw_execute(rplan); fftw_execute(iplan);
        if (sc)
        {
            const double s = 0.5;
            *Y++ = *Y1r++ * 0.5;
            *Y++ = *Y1i++ * 0.5;
            for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r * s; *++Y = *Y1i * s; }
        }
        else
        {
            for (size_t l=0; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r; *++Y = *Y1i; }
        }
        Y1r -= ndst; Y1i -= ndst;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
            for (size_t l=L; l<ndst; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
            X1r -= ndst; X1i -= ndst;
            fftw_execute(rplan); fftw_execute(iplan);
            if (sc)
            {
                const double s = 0.5, dcsc = 0.5;
                *Y++ = *Y1r++ * dcsc; *Y++ = *Y1i++ * dcsc;
                for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r * s; *++Y = *Y1i * s; }
                Y1r -= ndst; Y1i -= ndst;
                for (size_t v=1; v<V; ++v, Y1r-=ndst, Y1i-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftw_execute(rplan); fftw_execute(iplan);
                    *Y++ = *Y1r++ * dcsc; *Y++ = *Y1i++ * dcsc;
                    for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r * s; *++Y = *Y1i * s; }
                }
            }
            else
            {
                for (size_t l=0; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r; *++Y = *Y1i; }
                Y1r -= ndst; Y1i -= ndst;
                for (size_t v=1; v<V; ++v, Y1r-=ndst, Y1i-=ndst)
                {
                    for (size_t l=0; l<L; ++l, ++X, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftw_execute(rplan); fftw_execute(iplan);
                    for (size_t l=0; l<ndst; ++l, ++Y1r, ++Y1i, ++Y) { *Y = *Y1r; *++Y = *Y1i; }
                }
            }
        }
        else
        {
            const double s = sc ? 0.5 : 1.0;
            const double dcsc = sc ? 0.5 : 1.0;
            X1r += L; X1i += L;
            for (size_t l=L; l<ndst; ++l, ++X1r, ++X1i) { *X1r = *X1i = 0.0; }
            X1r -= ndst; X1i -= ndst;
            for (size_t g=0; g<G; ++g, X+=2*B*(L-1), Y+=2*B*(ndst-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*L-2, Y1r-=ndst, Y1i-=ndst, Y-=2*K*ndst-2)
                {
                    for (size_t l=0; l<L; ++l, X+=2*K-1, ++X1r, ++X1i) { *X1r = *X; *X1i = *++X; }
                    X1r -= L; X1i -= L;
                    fftw_execute(rplan); fftw_execute(iplan);
                    *Y = *Y1r++ * dcsc; *++Y = *Y1i++ * dcsc; Y += 2*K-1;
                    for (size_t l=1; l<ndst; ++l, ++Y1r, ++Y1i, Y+=2*K-1) { *Y = *Y1r*s; *++Y = *Y1i*s; }
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
