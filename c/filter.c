//Filters each vector in X along dim.
//Filter IIR coefficients are given in vector A with length P+1,
//where P is the IIR filter order (P=0 means only a0).
//Filter FIR coefficients are given in vector B with length Q+1,
//where Q is the FIR filter order (Q=0 means only b0).

//For each vector: a0*Y[t] = (b0*X[t] + b1*X[t-1] + ... + bQ*X[t-Q]) - (a1*Y[t-1] + ... + aP*Y[t-P])

//The calling program must ensure that the sizes are correct, the filter is stable, etc.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int filter_s (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filter_d (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filter_c (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);
int filter_z (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim);


int filter_s (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in filter_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually 1) and b0; also initialize Y
    float a = *A++, b = *B++;
    if (a==1.0f)
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    }
    else
    {
        const float b_a = b / a;
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A / a; }
        for (size_t q=1u; q<=Q; ++q, ++B) { *B /= a; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b_a * *X; }
        B -= Q;
    }
    A -= P; X -= N; Y -= N;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        //FIR
        if (L<3000u)
        {
            for (size_t q=1u; q<=Q; ++q, ++B, X-=L-q+1u, Y-=L-q)
            {
                b = *B;
                for (size_t l=q; l<L; ++l, ++X) { ++Y; *Y = fmaf(b,*X,*Y); }
            }
        }
        else
        {
            for (size_t q=1u; q<=Q; ++q, ++B) { ++Y; cblas_saxpy((int)(L-q),*B,X,1,Y,1); }
        }
        Y -= Q;

        //IIR
        for (size_t l=0u; l<P; A-=l, ++l, ++Y)
        {
            for (size_t p=1u; p<=l; ++p, ++A) { *Y = fmaf(*A,*(Y-p),*Y); }
        }
        for (size_t l=P; l<L; ++l, A-=P, ++Y)
        {
            for (size_t p=1u; p<=P; ++p, ++A) { *Y = fmaf(*A,*(Y-p),*Y); }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/BS;

        if (K==1u && (G==1u || BS==1u))
        {
            for (size_t v=V; v>0u; --v)
            {
                //FIR
                if (L<30000u)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B)
                    {
                        b = *B;
                        for (size_t l=q; l<L; ++l, ++X) { ++Y; *Y = fmaf(b,*X,*Y); }
                        if (q<Q) { X -= L-q; Y -= L-q-1u; }
                        else { X += Q; Y -= L-1u; }
                    }
                }
                else
                {
                    for (size_t q=1u; q<=Q; ++q, ++B) { ++Y; cblas_saxpy((int)(L-q),*B,X,1,Y,1); }
                    X += L; Y -= Q;
                }
                B -= Q;

                //IIR
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
            for (size_t g=G; g>0u; --g, X+=BS*(L-1u), Y+=BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, ++X, Y-=K*L-1u)
                {
                    //FIR
                    for (size_t q=1u; q<=Q; ++q, ++B) { Y+=K; cblas_saxpy((int)(L-q),*B,X,(int)K,Y,(int)K); }
                    Y -= Q*K; B -= Q;

                    //IIR
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


int filter_d (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in filter_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    //Deal with a0 (usually 1) and b0; also initialize Y
    double a = *A++, b = *B++;
    if (a==1.0)
    {
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b * *X; }
    }
    else
    {
        const double b_a = b / a;
        for (size_t p=1u; p<=P; ++p, ++A) { *A = -*A / a; }
        for (size_t q=1u; q<=Q; ++q, ++B) { *B /= a; }
        for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = b_a * *X; }
        B -= Q;
    }
    A -= P; X -= N; Y -= N;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        //FIR
        if (L<3000u)
        {
            for (size_t q=1u; q<=Q; ++q, ++B, X-=L-q+1u, Y-=L-q)
            {
                b = *B;
                for (size_t l=q; l<L; ++l, ++X) { ++Y; *Y = fma(b,*X,*Y); }
            }
        }
        else
        {
            for (size_t q=1u; q<=Q; ++q, ++B) { ++Y; cblas_daxpy((int)(L-q),*B,X,1,Y,1); }
        }
        Y -= Q;

        //IIR
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
            for (size_t v=V; v>0u; --v)
            {
                //FIR
                if (L<30000u)
                {
                    for (size_t q=1u; q<=Q; ++q, ++B)
                    {
                        b = *B;
                        for (size_t l=q; l<L; ++l, ++X) { ++Y; *Y = fma(b,*X,*Y); }
                        if (q<Q) { X -= L-q; Y -= L-q-1u; }
                        else { X += Q; Y -= L-1u; }
                    }
                }
                else
                {
                    for (size_t q=1u; q<=Q; ++q, ++B) { ++Y; cblas_daxpy((int)(L-q),*B,X,1,Y,1); }
                    X += L; Y -= Q;
                }
                B -= Q;

                //IIR
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
            for (size_t g=G; g>0u; --g, X+=BS*(L-1u), Y+=BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, ++X, Y-=K*L-1u)
                {
                    //FIR
                    for (size_t q=1u; q<=Q; ++q, ++B) { Y+=K; cblas_daxpy((int)(L-q),*B,X,(int)K,Y,(int)K); }
                    Y -= Q*K; B -= Q;

                    //IIR
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


int filter_c (float *Y, const float *X, float *A, float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in filter_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    float ar, ai, br, bi, xr, xi, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aP
    //Also, initialize Y
    const float b0r = *B++, b0i = *B++;
    const float a0r = *A++, a0i = *A++, a02 = a0r*a0r + a0i*a0i;
    const float b_ar = (b0r*a0r + b0i*a0i) / a02;
    const float b_ai = (b0i*a0r - b0r*a0i) / a02;
    for (size_t p=1u; p<=P; ++p)
    {
        ar = -*A; ai = -*(A+1u);
        *A++ = (ar*a0r + ai*a0i) / a02;
        *A++ = (ai*a0r - ar*a0i) / a02;
    }
    for (size_t q=1u; q<=Q; ++q)
    {
        br = *B++; bi = *B--;
        *B++ = (br*a0r + bi*a0i) / a02;
        *B++ = (bi*a0r - br*a0i) / a02;
    }
    for (size_t n=0u; n<N; ++n)
    {
        xr = *X++; xi = *X++;
        *Y++ = b_ar*xr - b_ai*xi;
        *Y++ = b_ar*xi + b_ai*xr;
    }
    A -= 2u*P; B -= 2u*Q; X -= 2u*N; Y -= 2u*N;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        //FIR
        for (size_t q=1u; q<=Q; ++q, X-=2u*(L-q+1u), Y-=2u*(L-q))
        {
            br = *B++; bi = *B++;
            for (size_t l=q; l<L; ++l)
            {
                Y += 2u;
                xr = *X++; xi = *X++;
                *Y++ += br*xr - bi*xi;
                *Y-- += br*xi + bi*xr;
            }
        }
        Y -= 2u*Q;

        //IIR
        for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2)
        {
            for (size_t p=1u; p<=l; ++p)
            {
                ar = *A++; ai = *A++;
                yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
        for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2)
        {
            for (size_t p=1u; p<=P; ++p)
            {
                ar = *A++; ai = *A++;
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
            for (size_t v=V; v>0u; --v)
            {
                //FIR
                for (size_t q=1u; q<=Q; ++q)
                {
                    br = *B++; bi = *B++;
                    for (size_t l=q; l<L; ++l)
                    {
                        Y += 2u;
                        xr = *X++; xi = *X++;
                        *Y++ += br*xr - bi*xi;
                        *Y-- += br*xi + bi*xr;
                    }
                    if (q<Q) { X -= 2u*(L-q); Y -= 2u*(L-q-1u); }
                    else { X += 2u*Q; Y -= 2u*L - 2u; }
                }
                B -= 2u*Q;

                //IIR
                for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2)
                {
                    for (size_t p=1u; p<=l; ++p)
                    {
                        ar = *A++; ai = *A++;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
                for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2)
                {
                    for (size_t p=1u; p<=P; ++p)
                    {
                        ar = *A++; ai = *A++;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                {
                    //FIR
                    //for (size_t q=1u; q<=Q; X-=2u*K*(L-q), ++q)
                    for (size_t q=1u; q<=Q; ++q)
                    {
                        br = *B++; bi = *B++;
                        for (size_t l=q; l<L; ++l, X+=2u*K)
                        {
                            Y += 2u*K;
                            xr = *X; xi = *(X+1);
                            *Y++ += br*xr - bi*xi;
                            *Y-- += br*xi + bi*xr;
                        }
                        if (q<Q) { X -= 2u*K*(L-q); Y -= 2u*K*(L-q-1u); }
                        else { X += 2u*K*Q; Y -= 2u*K*(L-1u); }
                    }
                    B -= 2u*Q;

                    //IIR
                    for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=l; ++p)
                        {
                            ar = *A++; ai = *A++;
                            yr = *(Y-2u*K*p); yi = *(Y-2u*K*p+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=P; ++p)
                        {
                            ar = *A++; ai = *A++;
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


int filter_z (double *Y, const double *X, double *A, double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t P, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in filter_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    double ar, ai, br, bi, xr, xi, yr, yi;

    //Deal with a0 (usually 1) and negate a1 to aP
    //Also, initialize Y
    const double b0r = *B++, b0i = *B++;
    const double a0r = *A++, a0i = *A++, a02 = a0r*a0r + a0i*a0i;
    const double b_ar = (b0r*a0r + b0i*a0i) / a02;
    const double b_ai = (b0i*a0r - b0r*a0i) / a02;
    for (size_t p=1u; p<=P; ++p)
    {
        ar = -*A; ai = -*(A+1u);
        *A++ = (ar*a0r + ai*a0i) / a02;
        *A++ = (ai*a0r - ar*a0i) / a02;
    }
    for (size_t q=1u; q<=Q; ++q)
    {
        br = *B++; bi = *B--;
        *B++ = (br*a0r + bi*a0i) / a02;
        *B++ = (bi*a0r - br*a0i) / a02;
    }
    for (size_t n=0u; n<N; ++n)
    {
        xr = *X++; xi = *X++;
        *Y++ = b_ar*xr - b_ai*xi;
        *Y++ = b_ar*xi + b_ai*xr;
    }
    A -= 2u*P; B -= 2u*Q; X -= 2u*N; Y -= 2u*N;

    if (N==0u || L==1u) {}
    else if (L==N)
    {
        //FIR
        for (size_t q=1u; q<=Q; ++q, X-=2u*(L-q+1u), Y-=2u*(L-q))
        {
            br = *B++; bi = *B++;
            for (size_t l=q; l<L; ++l)
            {
                Y += 2u;
                xr = *X++; xi = *X++;
                *Y++ += br*xr - bi*xi;
                *Y-- += br*xi + bi*xr;
            }
        }
        Y -= 2u*Q;

        //IIR
        for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2)
        {
            for (size_t p=1u; p<=l; ++p)
            {
                ar = *A++; ai = *A++;
                yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                *Y += ar*yr - ai*yi;
                *(Y+1u) += ar*yi + ai*yr;
            }
        }
        for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2)
        {
            for (size_t p=1u; p<=P; ++p)
            {
                ar = *A++; ai = *A++;
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
            for (size_t v=V; v>0u; --v)
            {
                //FIR
                for (size_t q=1u; q<=Q; ++q)
                {
                    br = *B++; bi = *B++;
                    for (size_t l=q; l<L; ++l)
                    {
                        Y += 2u;
                        xr = *X++; xi = *X++;
                        *Y++ += br*xr - bi*xi;
                        *Y-- += br*xi + bi*xr;
                    }
                    if (q<Q) { X -= 2u*(L-q); Y -= 2u*(L-q-1u); }
                    else { X += 2u*Q; Y -= 2u*L - 2u; }
                }
                B -= 2u*Q;

                //IIR
                for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2)
                {
                    for (size_t p=1u; p<=l; ++p)
                    {
                        ar = *A++; ai = *A++;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
                for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2)
                {
                    for (size_t p=1u; p<=P; ++p)
                    {
                        ar = *A++; ai = *A++;
                        yr = *(Y-2u*p); yi = *(Y-2u*p+1u);
                        *Y += ar*yr - ai*yi;
                        *(Y+1u) += ar*yi + ai*yr;
                    }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
            {
                for (size_t bs=0u; bs<BS; ++bs, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                {
                    //FIR
                    for (size_t q=1u; q<=Q; ++q)
                    {
                        br = *B++; bi = *B++;
                        for (size_t l=q; l<L; ++l, X+=2u*K)
                        {
                            Y += 2u*K;
                            xr = *X; xi = *(X+1);
                            *Y++ += br*xr - bi*xi;
                            *Y-- += br*xi + bi*xr;
                        }
                        if (q<Q) { X -= 2u*K*(L-q); Y -= 2u*K*(L-q-1u); }
                        else { X += 2u*K*Q; Y -= 2u*K*(L-1u); }
                    }
                    B -= 2u*Q;

                    //IIR
                    for (size_t l=0u; l<P; ++l, A-=2u*l-2u, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=l; ++p)
                        {
                            ar = *A++; ai = *A++;
                            yr = *(Y-2u*K*p); yi = *(Y-2u*K*p+1u);
                            *Y += ar*yr - ai*yi;
                            *(Y+1u) += ar*yi + ai*yr;
                        }
                    }
                    for (size_t l=P; l<L; ++l, A-=2u*P, Y+=2u*K)
                    {
                        for (size_t p=1u; p<=P; ++p)
                        {
                            ar = *A++; ai = *A++;
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


#ifdef __cplusplus
}
}
#endif
