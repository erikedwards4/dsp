//IIR filtering of each vector in X along dim.
//Filter coefficients are given in vector A with length P+1,
//where P is the IIR filter order (P=0 means only a0; P=1 means a0 and a1; etc.).

//For each vector: a0*Y[t] = X[t] - a1*Y[t-1] - ... - aP*Y[t-P]

//The calling program must ensure that the sizes are correct, the filter is stable, etc.

#include <stdio.h>
#include <math.h>
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int iir_s (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_d (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_c (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_z (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);

int iir_inplace_s (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_inplace_d (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_inplace_c (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);
int iir_inplace_z (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim);


int iir_s (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually 1), and negate a1 to aP
    //Also, initialize Y with X/a0
    const float a0 = *A++;
    if (a0==1.0f)
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X / a0; }
    }
    A -= P; Y -= N;

    if (N==0u || P==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; A-=l, ++l, ++Y)
        {
            for (size_t p=1u; p<=l; ++p, ++A) { *Y = fmaf(*A,*(Y-p),*Y); }
            //for (size_t p=1u; p<=l; ++p, ++A) { *Y += *A * *(Y-p); }
        }
        for (size_t l=P; l<L; ++l, A-=P, ++Y)
        {
            for (size_t p=1u; p<=P; ++p, ++A) { *Y = fmaf(*A,*(Y-p),*Y); }
            //for (size_t p=1u; p<=P; ++p, ++A) { *Y += *A * *(Y-p); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; A-=l, ++l, ++Y)
                {
                    for (size_t p=1u; p<=l; ++p, ++A) { *Y = fmaf(*A,*(Y-p),*Y); }
                }
                for (size_t l=P; l<L; ++l, A-=P, ++Y)
                {
                    for (size_t p=1u; p<=P; ++p, ++A) { *Y = fmaf(*A,*(Y-p),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<P; A-=l, ++l, Y+=K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A) { *Y = fmaf(*A,*(Y-K*p),*Y); }
                    }
                    for (size_t l=P; l<L; ++l, A-=P, Y+=K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A) { *Y = fmaf(*A,*(Y-K*p),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_d (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually 1) and negate a1 to aP
    //Also, initialize Y
    const double a0 = *A++;
    if (a0==1.0)
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X / a0; }
    }
    A -= P; Y -= N;

    if (N==0u || P==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; A-=l, ++l, ++Y)
        {
            for (size_t p=1u; p<=l; ++p, ++A) { *Y = fma(*A,*(Y-p),*Y); }
        }
        for (size_t l=P; l<L; ++l, A-=P, ++Y)
        {
            for (size_t p=1u; p<=P; ++p, ++A) { *Y = fma(*A,*(Y-p),*Y); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; A-=l, ++l, ++Y)
                {
                    for (size_t p=1u; p<=l; ++p, ++A) { *Y = fma(*A,*(Y-p),*Y); }
                }
                for (size_t l=P; l<L; ++l, A-=P, ++Y)
                {
                    for (size_t p=1u; p<=P; ++p, ++A) { *Y = fma(*A,*(Y-p),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<P; A-=l, ++l, Y+=K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A) { *Y = fma(*A,*(Y-K*p),*Y); }
                    }
                    for (size_t l=P; l<L; ++l, A-=P, Y+=K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A) { *Y = fma(*A,*(Y-K*p),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_c (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float ar, ai, xr, xi, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aP
    //Also, initialize Y
    const float a0r = *A++, a0i = *A++;
    const float a02 = a0r*a0r + a0i*a0i;
    if (a0r==1.0f && a0i==0.0f)
    {
        for (size_t p=1u; p<2u*P; ++p, ++A) { *A = -*A; }
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        for (size_t p=1u; p<=P; ++p, ++A)
        {
            ar = -*A; ai = -*(A+1u);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0u; n<N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = (xr*a0r + xi*a0i) / a02;
            *++Y = (xi*a0r - xr*a0i) / a02;
        }
    }
    A -= 2u*P; Y -= 2u*N;

    if (N==0u || L==1u || P==0u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u)
        {
            for (size_t p=1u; p<=l; ++p, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
        for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u)
        {
            for (size_t p=1u; p<=P; ++p, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u)
                {
                    for (size_t p=1u; p<=l; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
                for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u)
                {
                    for (size_t p=1u; p<=P; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=2u*BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, Y-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*p); yi = *(Y-2u*K*p+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*p); yi = *(Y-2u*K*p+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_z (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double ar, ai, xr, xi, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aP
    //Also, initialize Y
    const double a0r = *A++, a0i = *A++;
    if (a0r==1.0 && a0i==0.0)
    {
        for (size_t p=0u; p<2u*P; ++p, ++A) { *A = -*A; }
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const double a02 = a0r*a0r + a0i*a0i; 
        for (size_t p=1u; p<=P; ++p, ++A)
        {
            ar = -*A; ai = -*(A+1u);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0u; n<N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = (xr*a0r + xi*a0i) / a02;
            *++Y = (xi*a0r - xr*a0i) / a02;
        }
    }
    A -= 2u*P; Y -= 2u*N;

    if (N==0u || P==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u)
        {
            for (size_t p=1u; p<=l; ++p, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
        for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u)
        {
            for (size_t p=1u; p<=P; ++p, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u)
                {
                    for (size_t p=1u; p<=l; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
                for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u)
                {
                    for (size_t p=1u; p<=P; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=2u*BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, Y-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*p); yi = *(Y-2u*K*p+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*p); yi = *(Y-2u*K*p+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_inplace_s (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    //Deal with a0 (usually a0==1)
    const float a0 = *A++;
    if (a0==1.0f)
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A; }
    }
    else
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X) { *X /= a0; }
        X -= N;
    }
    A -= P;

    if (N==0u || P==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; A-=l, ++l, ++X)
        {
            for (size_t p=1u; p<=l; ++p, ++A) { *X = fmaf(*A,*(X-p),*X); }
        }
        for (size_t l=P; l<L; ++l, ++X, A-=P)
        {
            for (size_t p=1u; p<=P; ++p, ++A) { *X = fmaf(*A,*(X-p),*X); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; A-=l, ++l, ++X)
                {
                    for (size_t p=1u; p<=l; ++p, ++A) { *X = fmaf(*A,*(X-p),*X); }
                }
                for (size_t l=P; l<L; ++l, A-=P, ++X)
                {
                    for (size_t p=1u; p<=P; ++p, ++A) { *X = fmaf(*A,*(X-p),*X); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, X-=K*L-1u)
                {
                    for (size_t l=0u; l<P; A-=l, ++l, X+=K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A) { *X = fmaf(*A,*(X-K*p),*X); }
                    }
                    for (size_t l=P; l<L; ++l, A-=P, X+=K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A) { *X = fmaf(*A,*(X-K*p),*X); }
                    }
                }
            }
        }
    }
    
    clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);

    return 0;
}


int iir_inplace_d (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually a0==1)
    const double a0 = *A++;
    if (a0==1.0)
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A; }
    }
    else
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X) { *X /= a0; }
        X -= N;
    }
    A -= P;

    if (N==0u || P==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; ++l, ++X, A-=l-1u)
        {
            for (size_t p=1u; p<=l; ++p, ++A) { *X = fma(*A,*(X-p),*X); }
        }
        for (size_t l=P; l<L; ++l, ++X, A-=P)
        {
            for (size_t p=1u; p<=P; ++p, ++A) { *X = fma(*A,*(X-p),*X); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; A-=l, ++l, ++X)
                {
                    for (size_t p=1u; p<=l; ++p, ++A) { *X = fma(*A,*(X-p),*X); }
                }
                for (size_t l=P; l<L; ++l, A-=P, ++X)
                {
                    for (size_t p=1u; p<=P; ++p, ++A) { *X = fma(*A,*(X-p),*X); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, X-=K*L-1u)
                {
                    for (size_t l=0u; l<P; A-=l, ++l, X+=K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A) { *X = fma(*A,*(X-K*p),*X); }
                    }
                    for (size_t l=P; l<L; ++l, A-=P, X+=K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A) { *X = fma(*A,*(X-K*p),*X); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_inplace_c (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float ar, ai, xr, xi;

    //Deal with a0 (usually 1) and negate a1 to aP
    const float a0r = *A++, a0i = *A++;
    const float a02 = a0r*a0r + a0i*a0i;
    if (a0r==1.0f && a0i==0.0f)
    {
        for (size_t p=0u; p<2u*P; ++p, ++A) { *A = -*A; }
    }
    else
    {
        for (size_t p=1u; p<=P; ++p, ++A)
        {
            ar = -*A; ai = -*(A+1u);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0u; n<N; ++n, ++X)
        {
            xr = *X; xi = *(X+1u);
            *X = (xr*a0r + xi*a0i) / a02;
            *++X = (xi*a0r - xr*a0i) / a02;
        }
        X -= 2u*N;
    }
    A -= 2u*P;

    if (N==0u || P==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; ++l, A-=2u*l-2u, X+=2u)
        {
            for (size_t p=1u; p<=l; ++p, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2u*p); xi = *(X-2u*p+1u);
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
        for (size_t l=P; l<L; ++l, A-=2u*P, X+=2u)
        {
            for (size_t p=1u; p<=P; ++p, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2u*p); xi = *(X-2u*p+1u);
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; ++l, A-=2u*l-2u, X+=2u)
                {
                    for (size_t p=1u; p<=l; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2u*p); xi = *(X-2u*p+1u);
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
                for (size_t l=P; l<L; ++l, A-=2u*P, X+=2u)
                {
                    for (size_t p=1u; p<=P; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2u*p); xi = *(X-2u*p+1u);
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, X-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<P; ++l, A-=2u*l-2u, X+=2u*K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*p); xi = *(X-2u*K*p+1u);
                            *X += ar*xr - ai*xi;
                            *(X+1u) += ar*xi + ai*xr;
                        }
                    }
                    for (size_t l=P; l<L; ++l, A-=2u*P, X+=2u*K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*p); xi = *(X-2u*K*p+1u);
                            *X += ar*xr - ai*xi;
                            *(X+1u) += ar*xi + ai*xr;
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_inplace_z (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double ar, ai, xr, xi;

    //Deal with a0 (usually 1) and negate a1 to aP
    const double a0r = *A++, a0i = *A++;
    const double a02 = a0r*a0r + a0i*a0i;
    if (a0r==1.0 && a0i==0.0)
    {
        for (size_t p=0u; p<2u*P; ++p, ++A) { *A = -*A; }
    }
    else
    {
        for (size_t p=1u; p<=P; ++p, ++A)
        {
            ar = -*A; ai = -*(A+1u);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0u; n<N; ++n, ++X)
        {
            xr = *X; xi = *(X+1u);
            *X = (xr*a0r + xi*a0i) / a02;
            *++X = (xi*a0r - xr*a0i) / a02;
        }
        X -= 2u*N;
    }
    A -= 2u*P;

    if (N==0u || P==0u || L==1u) {}
    else if (L==N)
    {
        for (size_t l=0u; l<P; ++l, A-=2u*l-2u, X+=2u)
        {
            for (size_t p=1u; p<=l; ++p, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2u*p); xi = *(X-2u*p+1u);
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
        for (size_t l=P; l<L; ++l, A-=2u*P, X+=2u)
        {
            for (size_t p=1u; p<=P; ++p, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2u*p); xi = *(X-2u*p+1u);
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t l=0u; l<P; ++l, A-=2u*l-2u, X+=2u)
                {
                    for (size_t p=1u; p<=l; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2u*p); xi = *(X-2u*p+1u);
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
                for (size_t l=P; l<L; ++l, A-=2u*P, X+=2u)
                {
                    for (size_t p=1u; p<=P; ++p, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2u*p); xi = *(X-2u*p+1u);
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*BS*(L-1u))
            {
                for (size_t b=0u; b<BS; ++b, X-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<P; ++l, A-=2u*l-2u, X+=2u*K)
                    {
                        for (size_t p=1u; p<=l; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*p); xi = *(X-2u*K*p+1u);
                            *X += ar*xr - ai*xi;
                            *(X+1u) += ar*xi + ai*xr;
                        }
                    }
                    for (size_t l=P; l<L; ++l, A-=2u*P, X+=2u*K)
                    {
                        for (size_t p=1u; p<=P; ++p, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*p); xi = *(X-2u*K*p+1u);
                            *X += ar*xr - ai*xi;
                            *(X+1u) += ar*xi + ai*xr;
                        }
                    }
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
