//Causal FIR filtering of each vector in X along dim.
//FIR impulse response is given in vector B with length Q+1.
//(I use P for IIR filter order, since same as polynomial order.)

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fir_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);


int fir_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Initialize Y
    float b = *B++;
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    X -= N; Y -= N - 1u;

    if (N==0u || L==1u || Q==0u) {}
    else if (L==N)
    {
        if (L<30000u)
        {
            for (size_t q=1u; q<=Q; ++q, ++B, X-=L-q+1u, Y-=L-q)
            {
                b = *B;
                for (size_t l=q; l<L; ++l, ++X, ++Y) { *Y = fmaf(b,*X,*Y); }
            }
        }
        else
        {
            for (size_t q=1u; q<=Q; ++q, ++B, ++Y) { cblas_saxpy((int)(L-q),*B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (L<30000u)
            {
                for (size_t v=V; v>0u; --v, B-=Q, X+=Q, ++Y)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B)
                    {
                        b = *B;
                        for (size_t l=q; l<L; ++l, ++X, ++Y) { *Y = fmaf(b,*X,*Y); }
                        if (q<Q) { X -= L-q; Y -= L-q-1u; }
                    }
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, B-=Q, X+=L, Y+=L-Q)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B, ++Y) { cblas_saxpy((int)(L-q),*B,X,1,Y,1); }
                }
            }
        }
        else
        {
            Y += K - 1u;
            for (size_t g=G; g>0u; --g, X+=BS*(L-1u), Y+=BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, ++X, B-=Q, Y-=K*Q-1u)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B, Y+=K) { cblas_saxpy((int)(L-q),*B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Initialize Y
    double b = *B++;
    for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    X -= N; Y -= N - 1u;

    if (N==0u || L==1u || Q==0u) {}
    else if (L==N)
    {
        if (L<30000u)
        {
            for (size_t q=1u; q<=Q; ++q, ++B, X-=L-q+1u, Y-=L-q)
            {
                b = *B;
                for (size_t l=q; l<L; ++l, ++X, ++Y) { *Y = fma(b,*X,*Y); }
            }
        }
        else
        {
            for (size_t q=1u; q<=Q; ++q, ++B, ++Y) { cblas_daxpy((int)(L-q),*B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (L<30000u)
            {
                for (size_t v=V; v>0u; --v, B-=Q, X+=Q, ++Y)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B)
                    {
                        b = *B;
                        for (size_t l=q; l<L; ++l, ++X, ++Y) { *Y = fma(b,*X,*Y); }
                        if (q<Q) { X -= L-q; Y -= L-q-1u; }
                    }
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, B-=Q, X+=L, Y+=L-Q)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B, ++Y) { cblas_daxpy((int)(L-q),*B,X,1,Y,1); }
                }
            }
        }
        else
        {
            Y += K - 1u;
            for (size_t g=G; g>0u; --g, X+=BS*(L-1u), Y+=BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, ++X, B-=Q, Y-=K*Q-1u)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B, Y+=K) { cblas_daxpy((int)(L-q),*B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float xr, xi, br, bi;

    //Initialize Y
    br = *B++; bi = *B++;
    for (size_t n=0u; n<N; ++n)
    {
        xr = *X++; xi = *X++;
        *Y++ = br*xr - bi*xi;
        *Y++ = br*xi + bi*xr;
    }
    X -= 2u*N; Y -= 2u*N - 2u;

    if (N==0u || L==1u || Q==0u) {}
    else if (L==N)
    {
        if (L<30000u)
        {
            for (size_t q=1u; q<=Q; ++q, X-=2u*(L-q+1u), Y-=2u*(L-q))
            {
                br = *B++; bi = *B++;
                for (size_t l=q; l<L; ++l)
                {
                    xr = *X++; xi = *X++;
                    *Y++ += br*xr - bi*xi;
                    *Y++ += br*xi + bi*xr;
                }
            }
        }
        else
        {
            for (size_t q=1u; q<=Q; ++q, B+=2u, Y+=2u) { cblas_caxpy((int)(L-q),B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (L<30000u)
            {
                for (size_t v=V; v>0u; --v, B-=2u*Q, X+=2u*Q, Y+=2u)
                {
                    for (size_t q=1u; q<=Q; ++q)
                    {
                        br = *B++; bi = *B++;
                        for (size_t l=q; l<L; ++l)
                        {
                            xr = *X++; xi = *X++;
                            *Y++ += br*xr - bi*xi;
                            *Y++ += br*xi + bi*xr;
                        }
                        if (q<Q) { X -= 2u*(L-q); Y -= 2u*(L-q-1u); }
                    }
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, X+=2u*L, B-=2u*Q, Y+=2u*(L-Q-1u))
                {
                    for (size_t q=0u; q<=Q; ++q, B+=2u, Y+=2u) { cblas_caxpy((int)(L-q),B,X,1,Y,1); }
                }
            }
        }
        else
        {
            Y += 2u*K - 2u;
            for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, X+=2u, B-=2u*Q, Y-=2u*K*Q-2u)
                {
                    for (size_t q=1u; q<=Q; ++q, B+=2u, Y+=2u*K) { cblas_caxpy((int)(L-q),B,X,(int)K,Y,(int)K); }
                }
            }
        }
    }

    return 0;
}


int fir_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double xr, xi, br, bi;

    //Initialize Y
    br = *B++; bi = *B++;
    for (size_t n=0u; n<N; ++n)
    {
        xr = *X++; xi = *X++;
        *Y++ = br*xr - bi*xi;
        *Y++ = br*xi + bi*xr;
    }
    X -= 2u*N; Y -= 2u*N - 2u;

    if (N==0u || L==1u || Q==0u) {}
    else if (L==N)
    {
        if (L<30000u)
        {
            for (size_t q=1u; q<=Q; ++q, X-=2u*(L-q+1u), Y-=2u*(L-q))
            {
                br = *B++; bi = *B++;
                for (size_t l=q; l<L; ++l)
                {
                    xr = *X++; xi = *X++;
                    *Y++ += br*xr - bi*xi;
                    *Y++ += br*xi + bi*xr;
                }
            }
        }
        else
        {
            for (size_t q=1u; q<=Q; ++q, B+=2u, Y+=2u) { cblas_zaxpy((int)(L-q),B,X,1,Y,1); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            if (L<30000u)
            {
                for (size_t v=V; v>0u; --v, B-=2u*Q, X+=2u*Q, Y+=2u)
                {
                    for (size_t q=1u; q<=Q; ++q)
                    {
                        br = *B++; bi = *B++;
                        for (size_t l=q; l<L; ++l)
                        {
                            xr = *X++; xi = *X++;
                            *Y++ += br*xr - bi*xi;
                            *Y++ += br*xi + bi*xr;
                        }
                        if (q<Q) { X -= 2u*(L-q); Y -= 2u*(L-q-1u); }
                    }
                }
            }
            else
            {
                for (size_t v=V; v>0u; --v, X+=2u*L, B-=2u*Q, Y+=2u*(L-Q-1u))
                {
                    for (size_t q=0u; q<=Q; ++q, B+=2u, Y+=2u) { cblas_zaxpy((int)(L-q),B,X,1,Y,1); }
                }
            }
        }
        else
        {
            Y += 2u*K - 2u;
            for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, X+=2u, B-=2u*Q, Y-=2u*K*Q-2u)
                {
                    for (size_t q=1u; q<=Q; ++q, B+=2u, Y+=2u*K) { cblas_zaxpy((int)(L-q),B,X,(int)K,Y,(int)K); }
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
