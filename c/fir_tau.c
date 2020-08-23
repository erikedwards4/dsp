//Causal FIR filtering of each vector in X along dim.
//FIR impulse response is given in vector B with length L.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fir_tau_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim);
int fir_tau_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim);
int fir_tau_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim);
int fir_tau_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim);


int fir_tau_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Initialize Y to 0
    for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0f; }
    Y -= N;

    if (N==0 || L==0 || D>=(int)T || D<=-(int)(T+L)) {}
    else if (T==1)
    {
        const float b = *(B+D);
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    }
    else if (T==N)
    {
        if (T<30000)
        {
            float b;
            if (D<0)
            {
                for (size_t l=0; l<L; ++l, X-=T-l+1, ++B, Y-=T-l)
                {
                    b = *B;
                    for (size_t t=l; t<T; ++t, ++X, ++Y) { *Y = fmaf(b,*X,*Y); }
                }
            }
        }
        else
        {
            if (D<0)
            {
                for (size_t l=0; l<L; ++l, ++B, ++Y) { cblas_saxpy((int)(T-l),*B,X,1,Y,1); }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (T<30000)
            {
                float b;
                for (size_t v=0; v<V; ++v, X+=L-1, B-=L)
                {
                    for (size_t l=0; l<L; ++l, ++B)
                    {
                        b = *B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y) { *Y = fmaf(b,*X,*Y); }
                        if (l<L-1) { X -= T-l; Y -= T-l-1; }
                    }
                }
            }
            else
            {
                for (size_t v=0; v<V; ++v, X+=T, B-=L, Y+=T-L)
                {
                    for (size_t l=0; l<L; ++l, ++B, ++Y) { cblas_saxpy((int)(T-l),*B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=BS*(T-1), Y+=BS*(T-1))
            {
                for (size_t b=0; b<BS; ++b, ++X, B-=L, Y-=K*L-1)
                {
                    for (size_t l=0; l<L; ++l, ++B, Y+=K) { cblas_saxpy((int)(T-l),*B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_tau_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;

    //Initialize Y to 0
    for (size_t n=0; n<N; ++n, ++Y) { *Y = 0.0; }
    Y -= N;

    if (N==0 || L==0) {}
    else if (T==1)
    {
        const double b = *B;
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    }
    else if (T==N)
    {
        if (T<30000)
        {
            double b;
            for (size_t l=0; l<L; ++l, X-=T-l+1, ++B, Y-=T-l)
            {
                b = *B;
                for (size_t t=l; t<T; ++t, ++X, ++Y) { *Y = fma(b,*X,*Y); }
            }
        }
        else
        {
            for (size_t l=0; l<L; ++l, ++B, ++Y) { cblas_daxpy((int)(T-l),*B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (T<30000)
            {
                double b;
                for (size_t v=0; v<V; ++v, X+=L-1, B-=L)
                {
                    for (size_t l=0; l<L; ++l, ++B)
                    {
                        b = *B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y) { *Y = fma(b,*X,*Y); }
                        if (l<L-1) { X -= T-l; Y -= T-l-1; }
                    }
                }
            }
            else
            {
                for (size_t v=0; v<V; ++v, X+=T, B-=L, Y+=T-L)
                {
                    for (size_t l=0; l<L; ++l, ++B, ++Y) { cblas_daxpy((int)(T-l),*B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=BS*(T-1), Y+=BS*(T-1))
            {
                for (size_t b=0; b<BS; ++b, ++X, B-=L, Y-=K*L-1)
                {
                    for (size_t l=0; l<L; ++l, ++B, Y+=K) { cblas_daxpy((int)(T-l),*B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_tau_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    float xr, xi, br, bi;

    //Initialize Y to 0
    for (size_t n=0; n<2*N; ++n, ++Y) { *Y = 0.0f; }
    Y -= 2*N;

    if (N==0 || L==0) {}
    else if (T==1)
    {
        br = *B; bi = *++B;
        for (size_t n=0; n<2*N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = br*xr - bi*xi;
            *++Y = br*xi + bi*xr;
        }
    }
    else if (T==N)
    {
        if (T<30000)
        {
            for (size_t l=0; l<L; ++l, X-=2*(T-l+1), ++B, Y-=2*(T-l))
            {
                br = *B; bi = *++B;
                for (size_t t=l; t<T; ++t, ++X, ++Y)
                {
                    xr = *X; xi = *++X;
                    *Y += br*xr - bi*xi;
                    *++Y += br*xi + bi*xr;
                }
            }
        }
        else
        {
            for (size_t l=0; l<L; ++l, B+=2, Y+=2) { cblas_caxpy((int)(T-l),B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (T<30000)
            {
                for (size_t v=0; v<V; ++v, X+=2*L-2, B-=2*L)
                {
                    for (size_t l=0; l<L; ++l, ++B)
                    {
                        br = *B; bi = *++B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y)
                        {
                            xr = *X; xi = *++X;
                            *Y += br*xr - bi*xi;
                            *++Y += br*xi + bi*xr;
                        }
                        if (l<L-1) { X -= 2*(T-l); Y -= 2*(T-l-1); }
                    }
                }
            }
            else
            {
                for (size_t v=0; v<V; ++v, X+=2*T, B-=2*L, Y+=2*(T-L))
                {
                    for (size_t l=0; l<L; ++l, B+=2, Y+=2) { cblas_caxpy((int)(T-l),B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*BS*(T-1), Y+=2*BS*(T-1))
            {
                for (size_t b=0; b<BS; ++b, X+=2, B-=2*L, Y-=2*K*L-2)
                {
                    for (size_t l=0; l<L; ++l, B+=2, Y+=2*K) { cblas_caxpy((int)(T-l),B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_tau_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim)
{
    const size_t N = R*C*S*H;
    const size_t T = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    double xr, xi, br, bi;

    //Initialize Y to 0
    for (size_t n=0; n<2*N; ++n, ++Y) { *Y = 0.0; }
    Y -= 2*N;

    if (N==0 || L==0) {}
    else if (T==1)
    {
        br = *B; bi = *++B;
        for (size_t n=0; n<2*N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = br*xr - bi*xi;
            *++Y = br*xi + bi*xr;
        }
    }
    else if (T==N)
    {
        if (T<30000)
        {
            for (size_t l=0; l<L; ++l, X-=2*(T-l+1), ++B, Y-=2*(T-l))
            {
                br = *B; bi = *++B;
                for (size_t t=l; t<T; ++t, ++X, ++Y)
                {
                    xr = *X; xi = *++X;
                    *Y += br*xr - bi*xi;
                    *++Y += br*xi + bi*xr;
                }
            }
        }
        else
        {
            for (size_t l=0; l<L; ++l, B+=2, Y+=2) { cblas_zaxpy((int)(T-l),B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t BS = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1 && (G==1 || BS==1))
        {
            if (T<30000)
            {
                for (size_t v=0; v<V; ++v, X+=2*L-2, B-=2*L)
                {
                    for (size_t l=0; l<L; ++l, ++B)
                    {
                        br = *B; bi = *++B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y)
                        {
                            xr = *X; xi = *++X;
                            *Y += br*xr - bi*xi;
                            *++Y += br*xi + bi*xr;
                        }
                        if (l<L-1) { X -= 2*(T-l); Y -= 2*(T-l-1); }
                    }
                }
            }
            else
            {
                for (size_t v=0; v<V; ++v, X+=2*T, B-=2*L, Y+=2*(T-L))
                {
                    for (size_t l=0; l<L; ++l, B+=2, Y+=2) { cblas_zaxpy((int)(T-l),B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*BS*(T-1), Y+=2*BS*(T-1))
            {
                for (size_t b=0; b<BS; ++b, X+=2, B-=2*L, Y-=2*K*L-2)
                {
                    for (size_t l=0; l<L; ++l, B+=2, Y+=2*K) { cblas_zaxpy((int)(T-l),B,X,(int)K,Y,(int)K); }
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
