//IIR filtering of each vector in X along dim.
//Filter coefficients are given in vector A with length Q+1,
//where Q is the IIR filter order (Q=0 means only a0; Q=1 means a0 and a1; etc.).

//For each vector: a0*Y[t] = X[t] - a1*Y[t-1] - ... - aQ*Y[t-Q]

//The calling program must ensure that the sizes are correct, the filter is stable, etc.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int iir_s (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);
int iir_d (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);
int iir_c (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);
int iir_z (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);

int iir_inplace_s (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);
int iir_inplace_d (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);
int iir_inplace_c (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);
int iir_inplace_z (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim);


int iir_s (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0f)
    {
        ++A;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
        A -= Q; Y -= N;
    }
    else
    {
        const float a0 = *A++;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X / a0; }
        A -= Q; Y -= N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, A-=t-1u, ++Y)
        {
            for (size_t q=0u; q<t; ++q, ++A) { *Y = fmaf(*A,*(Y-(q+1u)),*Y); }
        }
        for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
        {
            for (size_t q=0u; q<Q; ++q, ++A) { *Y = fmaf(*A,*(Y-(q+1u)),*Y); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=t-1u, ++Y)
                {
                    for (size_t q=0u; q<t; ++q, ++A) { *Y = fmaf(*A,*(Y-(q+1u)),*Y); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
                {
                    for (size_t q=0u; q<Q; ++q, ++A) { *Y = fmaf(*A,*(Y-(q+1u)),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, Y-=K*T-1u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=t-1u, Y+=K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A) { *Y = fmaf(*A,*(Y-K*(q+1u)),*Y); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, Y+=K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A) { *Y = fmaf(*A,*(Y-K*(q+1u)),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_d (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0)
    {
        ++A;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; }
        A -= Q; Y -= N;
    }
    else
    {
        const double a0 = *A++;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X / a0; }
        A -= Q; Y -= N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, A-=t-1u, ++Y)
        {
            for (size_t q=0u; q<t; ++q, ++A) { *Y = fma(*A,*(Y-(q+1u)),*Y); }
        }
        for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
        {
            for (size_t q=0u; q<Q; ++q, ++A) { *Y = fma(*A,*(Y-(q+1u)),*Y); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=t-1u, ++Y)
                {
                    for (size_t q=0u; q<t; ++q, ++A) { *Y = fma(*A,*(Y-(q+1u)),*Y); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
                {
                    for (size_t q=0u; q<Q; ++q, ++A) { *Y = fma(*A,*(Y-(q+1u)),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, Y-=K*T-1u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=t-1u, Y+=K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A) { *Y = fma(*A,*(Y-K*(q+1u)),*Y); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, Y+=K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A) { *Y = fma(*A,*(Y-K*(q+1u)),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_c (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float ar, ai, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0f && *(A+1u)==0.0f)
    {
        A += 2u;
        for (size_t q=0u; q<2u*Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
        A -= 2u*Q; Y -= 2u*N;
    }
    else
    {
        const float a0r = *A++, a0i = *A++;
        const float a02 = a0r*a0r + a0i*a0i;
        float xr, xi;
        for (size_t q=0u; q<Q; ++q, ++A)
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
        A -= 2u*Q; Y -= 2u*N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, Y+=2u)
        {
            for (size_t q=0u; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2u*Q, Y+=2u)
        {
            for (size_t q=0u; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, Y+=2u)
                {
                    for (size_t q=0u; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2u*Q, Y+=2u)
                {
                    for (size_t q=0u; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=2u*B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, Y-=2u*K*T-2u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, Y+=2u*K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*(q+1u)); yi = *(Y-2u*K*(q+1u)+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2u*Q, Y+=2u*K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*(q+1u)); yi = *(Y-2u*K*(q+1u)+1u);
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


int iir_z (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in iir_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double ar, ai, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0 && *(A+1u)==0.0)
    {
        A += 2u;
        for (size_t q=0u; q<2u*Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
        A -= 2u*Q; Y -= 2u*N;
    }
    else
    {
        const double a0r = *A++, a0i = *A++;
        const double a02 = a0r*a0r + a0i*a0i;
        double xr, xi;
        for (size_t q=0u; q<Q; ++q, ++A)
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
        A -= 2u*Q; Y -= 2u*N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, Y+=2u)
        {
            for (size_t q=0u; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2u*Q, Y+=2u)
        {
            for (size_t q=0u; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, Y+=2u)
                {
                    for (size_t q=0u; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2u*Q, Y+=2u)
                {
                    for (size_t q=0u; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-(2u*q+2u)); yi = *(Y-(2u*q+1u));
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, Y+=2u*B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, Y-=2u*K*T-2u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, Y+=2u*K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*(q+1u)); yi = *(Y-2u*K*(q+1u)+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2u*Q, Y+=2u*K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2u*K*(q+1u)); yi = *(Y-2u*K*(q+1u)+1u);
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


int iir_inplace_s (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually a0==1)
    if (*A==1.0f)
    {
        ++A;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A; }
        A -= Q;
    }
    else
    {
        const float a0 = *A++;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X) { *X /= a0; }
        A -= Q; X -= N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, ++X, A-=t-1)
        {
            for (size_t q=0u; q<t; ++q, ++A) { *X = fmaf(*A,*(X-q-1u),*X); }
        }
        for (size_t t=Q; t<T; ++t, ++X, A-=Q)
        {
            for (size_t q=0u; q<Q; ++q, ++A) { *X = fmaf(*A,*(X-q-1u),*X); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=t-1u, ++X)
                {
                    for (size_t q=0u; q<t; ++q, ++A) { *X = fmaf(*A,*(X-q-1u),*X); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++X)
                {
                    for (size_t q=0u; q<Q; ++q, ++A) { *X = fmaf(*A,*(X-q-1u),*X); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*T-1u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=t-1u, X+=K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A) { *X = fmaf(*A,*(X-K*(q+1u)),*X); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, X+=K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A) { *X = fmaf(*A,*(X-K*(q+1u)),*X); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_inplace_d (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually a0==1)
    if (*A==1.0)
    {
        ++A;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A; }
        A -= Q;
    }
    else
    {
        const double a0 = *A++;
        for (size_t q=0u; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0u; n<N; ++n, ++X) { *X /= a0; }
        A -= Q; X -= N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, ++X, A-=t-1)
        {
            for (size_t q=0u; q<t; ++q, ++A) { *X = fma(*A,*(X-q-1u),*X); }
        }
        for (size_t t=Q; t<T; ++t, ++X, A-=Q)
        {
            for (size_t q=0u; q<Q; ++q, ++A) { *X = fma(*A,*(X-q-1u),*X); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=t-1u, ++X)
                {
                    for (size_t q=0u; q<t; ++q, ++A) { *X = fma(*A,*(X-q-1u),*X); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++X)
                {
                    for (size_t q=0u; q<Q; ++q, ++A) { *X = fma(*A,*(X-q-1u),*X); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*T-1u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=t-1u, X+=K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A) { *X = fma(*A,*(X-K*(q+1u)),*X); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, X+=K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A) { *X = fma(*A,*(X-K*(q+1u)),*X); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_inplace_c (float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float ar, ai, xr, xi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    if (*A==1.0f && *(A+1u)==0.0f)
    {
        A += 2u;
        for (size_t q=0u; q<2u*Q; ++q, ++A) { *A = -*A; }
        A -= 2u*Q;
    }
    else
    {
        const float a0r = *A++, a0i = *A++;
        const float a02 = a0r*a0r + a0i*a0i;
        float x1r, x1i;
        for (size_t q=0u; q<Q; ++q, ++A)
        {
            ar = -*A; ai = -*(A+1u);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0u; n<N; ++n, ++X)
        {
            x1r = *X; x1i = *(X+1u);
            *X = (x1r*a0r + x1i*a0i) / a02;
            *++X = (x1i*a0r - x1r*a0i) / a02;
        }
        A -= 2u*Q; X -= 2u*N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, X+=2u)
        {
            for (size_t q=0u; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2u*Q, X+=2u)
        {
            for (size_t q=0u; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, X+=2u)
                {
                    for (size_t q=0u; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2u*Q, X+=2u)
                {
                    for (size_t q=0u; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*T-2u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, X+=2u*K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*(q+1u)); xi = *(X-2u*K*(q+1u)+1u);
                            *X += ar*xr - ai*xi;
                            *(X+1u) += ar*xi + ai*xr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2u*Q, X+=2u*K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*(q+1u)); xi = *(X-2u*K*(q+1u)+1u);
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


int iir_inplace_z (double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double ar, ai, xr, xi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    if (*A==1.0 && *(A+1u)==0.0)
    {
        A += 2u;
        for (size_t q=0u; q<2u*Q; ++q, ++A) { *A = -*A; }
        A -= 2u*Q;
    }
    else
    {
        const double a0r = *A++, a0i = *A++;
        const double a02 = a0r*a0r + a0i*a0i;
        double x1r, x1i;
        for (size_t q=0u; q<Q; ++q, ++A)
        {
            ar = -*A; ai = -*(A+1u);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0u; n<N; ++n, ++X)
        {
            x1r = *X; x1i = *(X+1u);
            *X = (x1r*a0r + x1i*a0i) / a02;
            *++X = (x1i*a0r - x1r*a0i) / a02;
        }
        A -= 2u*Q; X -= 2u*N;
    }

    if (N==0u || Q==0u || T==1u) {}
    else if (T==N)
    {
        for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, X+=2u)
        {
            for (size_t q=0u; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2u*Q, X+=2u)
        {
            for (size_t q=0u; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                *X += ar*xr - ai*xi;
                *(X+1u) += ar*xi + ai*xr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, X+=2u)
                {
                    for (size_t q=0u; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2u*Q, X+=2u)
                {
                    for (size_t q=0u; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-(2u*q+2u)); xi = *(X-(2u*q+1u));
                        *X += ar*xr - ai*xi;
                        *(X+1u) += ar*xi + ai*xr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(T-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*T-2u)
                {
                    for (size_t t=0u; t<Q; ++t, A-=2u*t-2u, X+=2u*K)
                    {
                        for (size_t q=0u; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*(q+1u)); xi = *(X-2u*K*(q+1u)+1u);
                            *X += ar*xr - ai*xi;
                            *(X+1u) += ar*xi + ai*xr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2u*Q, X+=2u*K)
                    {
                        for (size_t q=0u; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2u*K*(q+1u)); xi = *(X-2u*K*(q+1u)+1u);
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
