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
    if (dim>3u) { fprintf(stderr,"error in fir_tau_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Initialize Y to 0
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0f; }
    Y -= N;

    if (N==0u || L==0u || D>=(int)T || D<=-(int)(T+L)) {}
    else if (T==1u)
    {
        const float b = *(B+D);
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    }
    else if (T==N)
    {
        if (T<30000u)
        {
            float b;
            if (D<0)
            {
                for (size_t l=0u; l<L; ++l, X-=T-l+1u, ++B, Y-=T-l)
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
                for (size_t l=0u; l<L; ++l, ++B, ++Y) { cblas_saxpy((int)(T-l),*B,X,1,Y,1); }
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (T<30000u)
            {
                float b;
                for (size_t v=0u; v<V; ++v, X+=L-1u, B-=L)
                {
                    for (size_t l=0u; l<L; ++l, ++B)
                    {
                        b = *B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y) { *Y = fmaf(b,*X,*Y); }
                        if (l<L-1u) { X -= T-l; Y -= T-l-1u; }
                    }
                }
            }
            else
            {
                for (size_t v=0u; v<V; ++v, X+=T, B-=L, Y+=T-L)
                {
                    for (size_t l=0u; l<L; ++l, ++B, ++Y) { cblas_saxpy((int)(T-l),*B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=BS*(T-1u), Y+=BS*(T-1u))
            {
                for (size_t b=0u; b<BS; ++b, ++X, B-=L, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, ++B, Y+=K) { cblas_saxpy((int)(T-l),*B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_tau_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_tau_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Initialize Y to 0
    for (size_t n=0u; n<N; ++n, ++Y) { *Y = 0.0; }
    Y -= N;

    if (N==0u || L==0u) {}
    else if (T==1u)
    {
        const double b = *B;
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    }
    else if (T==N)
    {
        if (T<30000u)
        {
            double b;
            for (size_t l=0u; l<L; ++l, X-=T-l+1u, ++B, Y-=T-l)
            {
                b = *B;
                for (size_t t=l; t<T; ++t, ++X, ++Y) { *Y = fma(b,*X,*Y); }
            }
        }
        else
        {
            for (size_t l=0u; l<L; ++l, ++B, ++Y) { cblas_daxpy((int)(T-l),*B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (T<30000u)
            {
                double b;
                for (size_t v=0u; v<V; ++v, X+=L-1u, B-=L)
                {
                    for (size_t l=0u; l<L; ++l, ++B)
                    {
                        b = *B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y) { *Y = fma(b,*X,*Y); }
                        if (l<L-1u) { X -= T-l; Y -= T-l-1u; }
                    }
                }
            }
            else
            {
                for (size_t v=0u; v<V; ++v, X+=T, B-=L, Y+=T-L)
                {
                    for (size_t l=0u; l<L; ++l, ++B, ++Y) { cblas_daxpy((int)(T-l),*B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=BS*(T-1u), Y+=BS*(T-1u))
            {
                for (size_t b=0u; b<BS; ++b, ++X, B-=L, Y-=K*L-1u)
                {
                    for (size_t l=0u; l<L; ++l, ++B, Y+=K) { cblas_daxpy((int)(T-l),*B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_tau_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_tau_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float xr, xi, br, bi;

    //Initialize Y to 0
    for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0f; }
    Y -= 2u*N;

    if (N==0u || L==0u) {}
    else if (T==1u)
    {
        br = *B; bi = *++B;
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = br*xr - bi*xi;
            *++Y = br*xi + bi*xr;
        }
    }
    else if (T==N)
    {
        if (T<30000u)
        {
            for (size_t l=0u; l<L; ++l, X-=2u*(T-l+1u), ++B, Y-=2u*(T-l))
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
            for (size_t l=0u; l<L; ++l, B+=2u, Y+=2u) { cblas_caxpy((int)(T-l),B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (T<30000u)
            {
                for (size_t v=0u; v<V; ++v, X+=2u*L-2u, B-=2u*L)
                {
                    for (size_t l=0u; l<L; ++l, ++B)
                    {
                        br = *B; bi = *++B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y)
                        {
                            xr = *X; xi = *++X;
                            *Y += br*xr - bi*xi;
                            *++Y += br*xi + bi*xr;
                        }
                        if (l<L-1u) { X -= 2u*(T-l); Y -= 2u*(T-l-1u); }
                    }
                }
            }
            else
            {
                for (size_t v=0u; v<V; ++v, X+=2u*T, B-=2u*L, Y+=2u*(T-L))
                {
                    for (size_t l=0u; l<L; ++l, B+=2u, Y+=2u) { cblas_caxpy((int)(T-l),B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*BS*(T-1u), Y+=2u*BS*(T-1u))
            {
                for (size_t b=0u; b<BS; ++b, X+=2u, B-=2u*L, Y-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<L; ++l, B+=2u, Y+=2u*K) { cblas_caxpy((int)(T-l),B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_tau_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const size_t L, const int D, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_tau_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t T = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double xr, xi, br, bi;

    //Initialize Y to 0
    for (size_t n=0u; n<2u*N; ++n, ++Y) { *Y = 0.0; }
    Y -= 2u*N;

    if (N==0u || L==0u) {}
    else if (T==1u)
    {
        br = *B; bi = *++B;
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y)
        {
            xr = *X; xi = *++X;
            *Y = br*xr - bi*xi;
            *++Y = br*xi + bi*xr;
        }
    }
    else if (T==N)
    {
        if (T<30000u)
        {
            for (size_t l=0u; l<L; ++l, X-=2u*(T-l+1u), ++B, Y-=2u*(T-l))
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
            for (size_t l=0u; l<L; ++l, B+=2u, Y+=2u) { cblas_zaxpy((int)(T-l),B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/T, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (T<30000u)
            {
                for (size_t v=0u; v<V; ++v, X+=2u*L-2u, B-=2u*L)
                {
                    for (size_t l=0u; l<L; ++l, ++B)
                    {
                        br = *B; bi = *++B;
                        for (size_t t=l; t<T; ++t, ++X, ++Y)
                        {
                            xr = *X; xi = *++X;
                            *Y += br*xr - bi*xi;
                            *++Y += br*xi + bi*xr;
                        }
                        if (l<L-1u) { X -= 2u*(T-l); Y -= 2u*(T-l-1u); }
                    }
                }
            }
            else
            {
                for (size_t v=0u; v<V; ++v, X+=2u*T, B-=2u*L, Y+=2u*(T-L))
                {
                    for (size_t l=0u; l<L; ++l, B+=2u, Y+=2u) { cblas_zaxpy((int)(T-l),B,X,1,Y,1); }
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*BS*(T-1u), Y+=2u*BS*(T-1u))
            {
                for (size_t b=0u; b<BS; ++b, X+=2u, B-=2u*L, Y-=2u*K*L-2u)
                {
                    for (size_t l=0u; l<L; ++l, B+=2u, Y+=2u*K) { cblas_zaxpy((int)(T-l),B,X,(int)K,Y,(int)K); }
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
