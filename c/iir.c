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
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0f)
    {
        ++A;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
        A -= Q; Y -= N;
    }
    else
    {
        const float a0 = *A++;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X / a0; }
        A -= Q; Y -= N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, A-=t-1, ++Y)
        {
            for (size_t q=0; q<t; ++q, ++A) { *Y = fmaf(*A,*(Y-q-1),*Y); }
        }
        for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
        {
            for (size_t q=0; q<Q; ++q, ++A) { *Y = fmaf(*A,*(Y-q-1),*Y); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=t-1, ++Y)
                {
                    for (size_t q=0; q<t; ++q, ++A) { *Y = fmaf(*A,*(Y-q-1),*Y); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
                {
                    for (size_t q=0; q<Q; ++q, ++A) { *Y = fmaf(*A,*(Y-q-1),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, Y+=B*(T-1))
            {
                for (size_t b=0; b<B; ++b, Y-=K*T-1)
                {
                    for (size_t t=0; t<Q; ++t, A-=t-1, Y+=K)
                    {
                        for (size_t q=0; q<t; ++q, ++A) { *Y = fmaf(*A,*(Y-K*(q+1)),*Y); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, Y+=K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A) { *Y = fmaf(*A,*(Y-K*(q+1)),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_d (double *Y, const double *X, double *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0)
    {
        ++A;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; }
        A -= Q; Y -= N;
    }
    else
    {
        const double a0 = *A++;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X / a0; }
        A -= Q; Y -= N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, A-=t-1, ++Y)
        {
            for (size_t q=0; q<t; ++q, ++A) { *Y = fma(*A,*(Y-q-1),*Y); }
        }
        for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
        {
            for (size_t q=0; q<Q; ++q, ++A) { *Y = fma(*A,*(Y-q-1),*Y); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=t-1, ++Y)
                {
                    for (size_t q=0; q<t; ++q, ++A) { *Y = fma(*A,*(Y-q-1),*Y); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++Y)
                {
                    for (size_t q=0; q<Q; ++q, ++A) { *Y = fma(*A,*(Y-q-1),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, Y+=B*(T-1))
            {
                for (size_t b=0; b<B; ++b, Y-=K*T-1)
                {
                    for (size_t t=0; t<Q; ++t, A-=t-1, Y+=K)
                    {
                        for (size_t q=0; q<t; ++q, ++A) { *Y = fma(*A,*(Y-K*(q+1)),*Y); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, Y+=K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A) { *Y = fma(*A,*(Y-K*(q+1)),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int iir_c (float *Y, const float *X, float *A, const size_t R, const size_t C, const size_t S, const size_t H, const size_t Q, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float ar, ai, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0f && *(A+1)==0.0f)
    {
        A += 2;
        for (size_t q=0; q<2*Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
        A -= 2*Q; Y -= 2*N;
    }
    else
    {
        const float a0r = *A++, a0i = *A++;
        const float a02 = a0r*a0r + a0i*a0i;
        float xr, xi;
        for (size_t q=0; q<Q; ++q, ++A)
        {
            ar = -*A; ai = -*(A+1);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0; n<N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = (xr*a0r + xi*a0i) / a02;
            *++Y = (xi*a0r - xr*a0i) / a02;
        }
        A -= 2*Q; Y -= 2*N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, A-=2*t-2, Y+=2)
        {
            for (size_t q=0; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                *Y += ar*yr - ai*yi;
                *(Y+1) += ar*yi + ai*yr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2*Q, Y+=2)
        {
            for (size_t q=0; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                *Y += ar*yr - ai*yi;
                *(Y+1) += ar*yi + ai*yr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=2*t-2, Y+=2)
                {
                    for (size_t q=0; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                        *Y += ar*yr - ai*yi;
                        *(Y+1) += ar*yi + ai*yr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2*Q, Y+=2)
                {
                    for (size_t q=0; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                        *Y += ar*yr - ai*yi;
                        *(Y+1) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, Y+=2*B*(T-1))
            {
                for (size_t b=0; b<B; ++b, Y-=2*K*T-2)
                {
                    for (size_t t=0; t<Q; ++t, A-=2*t-2, Y+=2*K)
                    {
                        for (size_t q=0; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2*K*(q+1)); yi = *(Y-2*K*(q+1)+1);
                            *Y += ar*yr - ai*yi;
                            *(Y+1) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2*Q, Y+=2*K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2*K*(q+1)); yi = *(Y-2*K*(q+1)+1);
                            *Y += ar*yr - ai*yi;
                            *(Y+1) += ar*yi + ai*yr;
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
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double ar, ai, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    //Also, initialize Y (copy X to Y)
    if (*A==1.0 && *(A+1)==0.0)
    {
        A += 2;
        for (size_t q=0; q<2*Q; ++q, ++A) { *A = -*A; }
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
        A -= 2*Q; Y -= 2*N;
    }
    else
    {
        const double a0r = *A++, a0i = *A++;
        const double a02 = a0r*a0r + a0i*a0i;
        double xr, xi;
        for (size_t q=0; q<Q; ++q, ++A)
        {
            ar = -*A; ai = -*(A+1);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0; n<N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = (xr*a0r + xi*a0i) / a02;
            *++Y = (xi*a0r - xr*a0i) / a02;
        }
        A -= 2*Q; Y -= 2*N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, A-=2*t-2, Y+=2)
        {
            for (size_t q=0; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                *Y += ar*yr - ai*yi;
                *(Y+1) += ar*yi + ai*yr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2*Q, Y+=2)
        {
            for (size_t q=0; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                *Y += ar*yr - ai*yi;
                *(Y+1) += ar*yi + ai*yr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=2*t-2, Y+=2)
                {
                    for (size_t q=0; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                        *Y += ar*yr - ai*yi;
                        *(Y+1) += ar*yi + ai*yr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2*Q, Y+=2)
                {
                    for (size_t q=0; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        yr = *(Y-2*q-2); yi = *(Y-2*q-1);
                        *Y += ar*yr - ai*yi;
                        *(Y+1) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, Y+=2*B*(T-1))
            {
                for (size_t b=0; b<B; ++b, Y-=2*K*T-2)
                {
                    for (size_t t=0; t<Q; ++t, A-=2*t-2, Y+=2*K)
                    {
                        for (size_t q=0; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2*K*(q+1)); yi = *(Y-2*K*(q+1)+1);
                            *Y += ar*yr - ai*yi;
                            *(Y+1) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2*Q, Y+=2*K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            yr = *(Y-2*K*(q+1)); yi = *(Y-2*K*(q+1)+1);
                            *Y += ar*yr - ai*yi;
                            *(Y+1) += ar*yi + ai*yr;
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
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Deal with a0 (usually a0==1)
    if (*A==1.0f)
    {
        ++A;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A; }
        A -= Q;
    }
    else
    {
        const float a0 = *A++;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0; n<N; ++n, ++X) { *X /= a0; }
        A -= Q; X -= N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, ++X, A-=t-1)
        {
            for (size_t q=0; q<t; ++q, ++A) { *X = fmaf(*A,*(X-q-1),*X); }
        }
        for (size_t t=Q; t<T; ++t, ++X, A-=Q)
        {
            for (size_t q=0; q<Q; ++q, ++A) { *X = fmaf(*A,*(X-q-1),*X); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=t-1, ++X)
                {
                    for (size_t q=0; q<t; ++q, ++A) { *X = fmaf(*A,*(X-q-1),*X); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++X)
                {
                    for (size_t q=0; q<Q; ++q, ++A) { *X = fmaf(*A,*(X-q-1),*X); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(T-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*T-1)
                {
                    for (size_t t=0; t<Q; ++t, A-=t-1, X+=K)
                    {
                        for (size_t q=0; q<t; ++q, ++A) { *X = fmaf(*A,*(X-K*(q+1)),*X); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, X+=K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A) { *X = fmaf(*A,*(X-K*(q+1)),*X); }
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
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Deal with a0 (usually a0==1)
    if (*A==1.0)
    {
        ++A;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A; }
        A -= Q;
    }
    else
    {
        const double a0 = *A++;
        for (size_t q=0; q<Q; ++q, ++A) { *A = -*A / a0; }
        for (size_t n=0; n<N; ++n, ++X) { *X /= a0; }
        A -= Q; X -= N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, ++X, A-=t-1)
        {
            for (size_t q=0; q<t; ++q, ++A) { *X = fma(*A,*(X-q-1),*X); }
        }
        for (size_t t=Q; t<T; ++t, ++X, A-=Q)
        {
            for (size_t q=0; q<Q; ++q, ++A) { *X = fma(*A,*(X-q-1),*X); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=t-1, ++X)
                {
                    for (size_t q=0; q<t; ++q, ++A) { *X = fma(*A,*(X-q-1),*X); }
                }
                for (size_t t=Q; t<T; ++t, A-=Q, ++X)
                {
                    for (size_t q=0; q<Q; ++q, ++A) { *X = fma(*A,*(X-q-1),*X); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=B*(T-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*T-1)
                {
                    for (size_t t=0; t<Q; ++t, A-=t-1, X+=K)
                    {
                        for (size_t q=0; q<t; ++q, ++A) { *X = fma(*A,*(X-K*(q+1)),*X); }
                    }
                    for (size_t t=Q; t<T; ++t, A-=Q, X+=K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A) { *X = fma(*A,*(X-K*(q+1)),*X); }
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
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float ar, ai, xr, xi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    if (*A==1.0f && *(A+1)==0.0f)
    {
        A += 2;
        for (size_t q=0; q<2*Q; ++q, ++A) { *A = -*A; }
        A -= 2*Q;
    }
    else
    {
        const float a0r = *A++, a0i = *A++;
        const float a02 = a0r*a0r + a0i*a0i;
        float x1r, x1i;
        for (size_t q=0; q<Q; ++q, ++A)
        {
            ar = -*A; ai = -*(A+1);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0; n<N; ++n, ++X)
        {
            x1r = *X; x1i = *(X+1);
            *X = (x1r*a0r + x1i*a0i) / a02;
            *++X = (x1i*a0r - x1r*a0i) / a02;
        }
        A -= 2*Q; X -= 2*N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, A-=2*t-2, X+=2)
        {
            for (size_t q=0; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2*q-2); xi = *(X-2*q-1);
                *X += ar*xr - ai*xi;
                *(X+1) += ar*xi + ai*xr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2*Q, X+=2)
        {
            for (size_t q=0; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2*q-2); xi = *(X-2*q-1);
                *X += ar*xr - ai*xi;
                *(X+1) += ar*xi + ai*xr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=2*t-2, X+=2)
                {
                    for (size_t q=0; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2*q-2); xi = *(X-2*q-1);
                        *X += ar*xr - ai*xi;
                        *(X+1) += ar*xi + ai*xr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2*Q, X+=2)
                {
                    for (size_t q=0; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2*q-2); xi = *(X-2*q-1);
                        *X += ar*xr - ai*xi;
                        *(X+1) += ar*xi + ai*xr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(T-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*T-2)
                {
                    for (size_t t=0; t<Q; ++t, A-=2*t-2, X+=2*K)
                    {
                        for (size_t q=0; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2*K*(q+1)); xi = *(X-2*K*(q+1)+1);
                            *X += ar*xr - ai*xi;
                            *(X+1) += ar*xi + ai*xr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2*Q, X+=2*K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2*K*(q+1)); xi = *(X-2*K*(q+1)+1);
                            *X += ar*xr - ai*xi;
                            *(X+1) += ar*xi + ai*xr;
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
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double ar, ai, xr, xi;

    //Deal with a0 (usually 1) and negate a1 to aQ
    if (*A==1.0 && *(A+1)==0.0)
    {
        A += 2;
        for (size_t q=0; q<2*Q; ++q, ++A) { *A = -*A; }
        A -= 2*Q;
    }
    else
    {
        const double a0r = *A++, a0i = *A++;
        const double a02 = a0r*a0r + a0i*a0i;
        double x1r, x1i;
        for (size_t q=0; q<Q; ++q, ++A)
        {
            ar = -*A; ai = -*(A+1);
            *A = (ar*a0r + ai*a0i) / a02;
            *++A = (ai*a0r - ar*a0i) / a02;
        }
        for (size_t n=0; n<N; ++n, ++X)
        {
            x1r = *X; x1i = *(X+1);
            *X = (x1r*a0r + x1i*a0i) / a02;
            *++X = (x1i*a0r - x1r*a0i) / a02;
        }
        A -= 2*Q; X -= 2*N;
    }

    if (N==0 || Q==0 || T==1) {}
    else if (T==N)
    {
        for (size_t t=0; t<Q; ++t, A-=2*t-2, X+=2)
        {
            for (size_t q=0; q<t; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2*q-2); xi = *(X-2*q-1);
                *X += ar*xr - ai*xi;
                *(X+1) += ar*xi + ai*xr;
            }
        }
        for (size_t t=Q; t<T; ++t, A-=2*Q, X+=2)
        {
            for (size_t q=0; q<Q; ++q, ++A)
            {
                ar = *A; ai = *++A;
                xr = *(X-2*q-2); xi = *(X-2*q-1);
                *X += ar*xr - ai*xi;
                *(X+1) += ar*xi + ai*xr;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t v=0; v<V; ++v)
            {
                for (size_t t=0; t<Q; ++t, A-=2*t-2, X+=2)
                {
                    for (size_t q=0; q<t; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2*q-2); xi = *(X-2*q-1);
                        *X += ar*xr - ai*xi;
                        *(X+1) += ar*xi + ai*xr;
                    }
                }
                for (size_t t=Q; t<T; ++t, A-=2*Q, X+=2)
                {
                    for (size_t q=0; q<Q; ++q, ++A)
                    {
                        ar = *A; ai = *++A;
                        xr = *(X-2*q-2); xi = *(X-2*q-1);
                        *X += ar*xr - ai*xi;
                        *(X+1) += ar*xi + ai*xr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(T-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*T-2)
                {
                    for (size_t t=0; t<Q; ++t, A-=2*t-2, X+=2*K)
                    {
                        for (size_t q=0; q<t; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2*K*(q+1)); xi = *(X-2*K*(q+1)+1);
                            *X += ar*xr - ai*xi;
                            *(X+1) += ar*xi + ai*xr;
                        }
                    }
                    for (size_t t=Q; t<T; ++t, A-=2*Q, X+=2*K)
                    {
                        for (size_t q=0; q<Q; ++q, ++A)
                        {
                            ar = *A; ai = *++A;
                            xr = *(X-2*K*(q+1)); xi = *(X-2*K*(q+1)+1);
                            *X += ar*xr - ai*xi;
                            *(X+1) += ar*xi + ai*xr;
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
