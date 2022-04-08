//Gets polynomial params (Y), and error variance (E), for each vector in X.
//Works along dim, such that each signal (vector) in X is converted
//to one scalar in E and one vector (of length L) in Y.

//By default, the means of X are NOT subtracted.
//The mnz option allows the mean to be zeroed first for each vec in X.
//However, this also required the removal of the const qualifier for input *X.

//This starts by getting the autocovariance (AC) for lags 0 to L-1 for each vector in X.
//The "biased" version of AC uses N in the denominator,
//whereas the "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower, has worse MSE, and doesn't match FFT estimate.

//After getting the AC, Levinson-Durbin recursion is used to get the polynomial params (Y)
//and the error variance (E).

//This matches Octave in tests. For complex case, the sign of the imaginary part depends
//on whether the positive or negative lags of the ACF are taken.

#include <stdio.h>
#include <stdlib.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int sig2poly_s (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2poly_s: dim must be in [0 3]\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in sig2poly_s: L must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t P = L - 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2poly_s: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        if (Lx==N)
        {
            if (mnz)
            {
                float mn = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                mn /= (float)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
            }
            float sm = 0.0f;
            for (size_t n=Lx; n>0u; --n, ++X) { sm += *X * *X; }
            *E = sm / (float)Lx;
            *Y = 1.0f;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    if (mnz)
                    {
                        float mn = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                        mn /= (float)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
                    }
                    float sm = 0.0f;
                    for (size_t n=Lx; n>0u; --n, ++X) { sm += *X * *X; }
                    *E++ = sm / (float)Lx;
                    *Y++ = 1.0f;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u)
                    {
                        if (mnz)
                        {
                            float mn = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=K) { mn += *X; }
                            mn /= (float)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=K; *X -= mn; }
                        }
                        float sm = 0.0f;
                        for (size_t n=Lx; n>0u; --n, X+=K) { sm += *X * *X; }
                        *E++ = sm / (float)Lx;
                        *Y++ = 1.0f;
                    }
                }
            }
        }
    }
    else //L>1u
    {
        float *AC, *A, a, sm, e;
        if (!(AC=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sig2poly_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in sig2poly_s: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            //Subtract mean
            if (mnz)
            {
                float mn = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                mn /= (float)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
            }

            //Get ACF
            for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++AC)
            {
                sm = 0.0f;
                for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                *AC = sm;
            }

            //Normalize ACF
            if (unbiased)
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (float)(Lx-l+1u); }
            }
            else
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (float)Lx; }
            }

            //Get poly params and error var (Lev-Durb)
            *Y++ = 1.0f;
            e = *AC++;
            a = -*AC / e; *Y = a;
            e += a * *AC++;
            for (size_t p=1u; p<P; ++p, AC+=p)
            {
                a = *AC;
                for (size_t q=p; q>0u; --q, ++Y) { --AC; a += *AC * *Y; }
                a /= -e; *Y = a;
                for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                Y += p;
                for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                e *= 1.0f - a*a;
            }
            AC -= L;
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Lx, Y+=P, ++E)
                {
                    //Subtract mean
                    if (mnz)
                    {
                        float mn = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                        mn /= (float)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
                    }

                    //Get ACF
                    for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++AC)
                    {
                        sm = 0.0f;
                        for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                        *AC = sm;
                    }

                    //Normalize ACF
                    if (unbiased)
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (float)(Lx-l+1u); }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (float)Lx; }
                    }

                    //Get poly params and error var (Lev-Durb)
                    *Y++ = 1.0f;
                    e = *AC++;
                    a = -*AC / e; *Y = a;
                    e += a * *AC++;
                    for (size_t p=1u; p<P; ++p, AC+=p)
                    {
                        a = *AC;
                        for (size_t q=p; q>0u; --q, ++Y) { --AC; a += *AC * *Y; }
                        a /= -e; *Y = a;
                        for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                        Y += p;
                        for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                        e *= 1.0f - a*a;
                    }
                    AC -= L;
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K-1u, ++E)
                    {
                        //Subtract mean
                        if (mnz)
                        {
                            float mn = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=K) { mn += *X; }
                            mn /= (float)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=K; *X -= mn; }
                        }

                        //Get ACF
                        for (size_t l=0u; l<L; ++l, X-=K*(Lx-l+1u), ++AC)
                        {
                            sm = 0.0f;
                            for (size_t n=Lx-l; n>0u; --n, X+=K) { sm += *X * *(X+l*K); }
                            *AC = sm;
                        }

                        //Normalize ACF
                        if (unbiased)
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (float)(Lx-l+1u); }
                        }
                        else
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (float)Lx; }
                        }

                        //Get poly params and error var (Lev-Durb)
                        *Y = 1.0f; Y += K;
                        e = *AC++;
                        a = -*AC / e; *Y = a;
                        e += a * *AC++;
                        for (size_t p=1u; p<P; ++p, AC+=p)
                        {
                            a = *AC;
                            for (size_t q=p; q>0u; --q, Y+=K) { --AC; a += *AC * *Y; }
                            a /= -e; *Y = a;
                            for (size_t q=p; q>0u; --q, ++A) { Y-=K; *A = *Y; }
                            Y += p*K;
                            for (size_t q=p; q>0u; --q) { --A; Y-=K; *Y += a * *A; }
                            e *= 1.0f - a*a;
                        }
                        AC -= L;
                        *E = e;
                    }
                }
            }
        }
        free(AC); free(A);
    }

    return 0;
}


int sig2poly_d (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2poly_d: dim must be in [0 3]\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in sig2ar_d: L must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t P = L - 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2poly_d: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        if (Lx==N)
        {
            if (mnz)
            {
                double mn = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                mn /= (double)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
            }
            double sm = 0.0;
            for (size_t n=Lx; n>0u; --n, ++X) { sm += *X * *X; }
            *E = sm / (double)Lx;
            *Y = 1.0;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    if (mnz)
                    {
                        double mn = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                        mn /= (double)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
                    }
                    double sm = 0.0;
                    for (size_t n=Lx; n>0u; --n, ++X) { sm += *X * *X; }
                    *E++ = sm / (double)Lx;
                    *Y++ = 1.0;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u)
                    {
                        if (mnz)
                        {
                            double mn = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=K) { mn += *X; }
                            mn /= (double)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=K; *X -= mn; }
                        }
                        double sm = 0.0;
                        for (size_t n=Lx; n>0u; --n, X+=K) { sm += *X * *X; }
                        *E++ = sm / (double)Lx;
                        *Y++ = 1.0;
                    }
                }
            }
        }
    }
    else //L>1u
    {
        double *AC, *A, a, sm, e;
        if (!(AC=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in sig2poly_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in sig2poly_d: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            //Subtract mean
            if (mnz)
            {
                double mn = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                mn /= (double)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
            }

            //Get ACF
            for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++AC)
            {
                sm = 0.0;
                for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                *AC = sm;
            }

            //Normalize ACF
            if (unbiased)
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (double)(Lx-l+1u); }
            }
            else
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (double)Lx; }
            }

            //Get poly params and error var (Lev-Durb)
            *Y++ = 1.0;
            e = *AC++;
            a = -*AC / e; *Y = a;
            e += a * *AC++;
            for (size_t p=1u; p<P; ++p, AC+=p)
            {
                a = *AC;
                for (size_t q=p; q>0u; --q, ++Y) { --AC; a += *AC * *Y; }
                a /= -e; *Y = a;
                for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                Y += p;
                for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                e *= 1.0 - a*a;
            }
            AC -= L;
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Lx, Y+=P, ++E)
                {
                    //Subtract mean
                    if (mnz)
                    {
                        double mn = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X) { mn += *X; }
                        mn /= (double)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mn; }
                    }

                    //Get ACF
                    for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++AC)
                    {
                        sm = 0.0;
                        for (size_t n=Lx-l; n>0u; --n, ++X) { sm += *X * *(X+l); }
                        *AC = sm;
                    }

                    //Normalize ACF
                    if (unbiased)
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (double)(Lx-l+1u); }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (double)Lx; }
                    }

                    //Get poly params and error var (Lev-Durb)
                    *Y++ = 1.0;
                    e = *AC++;
                    a = -*AC / e; *Y = a;
                    e += a * *AC++;
                    for (size_t p=1u; p<P; ++p, AC+=p)
                    {
                        a = *AC;
                        for (size_t q=p; q>0u; --q, ++Y) { --AC; a += *AC * *Y; }
                        a /= -e; *Y = a;
                        for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                        Y += p;
                        for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                        e *= 1.0 - a*a;
                    }
                    AC -= L;
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(L-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K-1u, ++E)
                    {
                        //Subtract mean
                        if (mnz)
                        {
                            double mn = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=K) { mn += *X; }
                            mn /= (double)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=K; *X -= mn; }
                        }

                        //Get ACF
                        for (size_t l=0u; l<L; ++l, X-=K*(Lx-l+1u), ++AC)
                        {
                            sm = 0.0;
                            for (size_t n=Lx-l; n>0u; --n, X+=K) { sm += *X * *(X+l*K); }
                            *AC = sm;
                        }

                        //Normalize ACF
                        if (unbiased)
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (double)(Lx-l+1u); }
                        }
                        else
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (double)Lx; }
                        }

                        //Get poly params and error var (Lev-Durb)
                        *Y = 1.0; Y += K;
                        e = *AC++;
                        a = -*AC / e; *Y = a;
                        e += a * *AC++;
                        for (size_t p=1u; p<P; ++p, AC+=p)
                        {
                            a = *AC;
                            for (size_t q=p; q>0u; --q, Y+=K) { --AC; a += *AC * *Y; }
                            a /= -e; *Y = a;
                            for (size_t q=p; q>0u; --q, ++A) { Y-=K; *A = *Y; }
                            Y += p*K;
                            for (size_t q=p; q>0u; --q) { --A; Y-=K; *Y += a * *A; }
                            e *= 1.0 - a*a;
                        }
                        AC -= L;
                        *E = e;
                    }
                }
            }
        }
        free(AC); free(A);
    }

    return 0;
}


int sig2poly_c (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2poly_c: dim must be in [0 3]\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in sig2ar_c: L must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t P = L - 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2poly_c: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        if (Lx==N)
        {
            if (mnz)
            {
                float mnr = 0.0f, mni = 0.0f;
                for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                mnr /= (float)Lx; mni /= (float)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
            }
            float sm = 0.0f;
            for (size_t n=2u*Lx; n>0u; --n, ++X) { sm += *X * *X; }
            *E = sm / (float)Lx;
            *Y++ = 1.0f; *Y = 0.0f;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    if (mnz)
                    {
                        float mnr = 0.0f, mni = 0.0f;
                        for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                        mnr /= (float)Lx; mni /= (float)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
                    }
                    float sm = 0.0f;
                    for (size_t n=2u*Lx; n>0u; --n, ++X) { sm += *X * *X; }
                    *E++ = sm / (float)Lx;
                    *Y++ = 1.0f; *Y++ = 0.0f;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u)
                    {
                        if (mnz)
                        {
                            float mnr = 0.0f, mni = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K) { mnr += *X; mni += *(X+1); }
                            mnr /= (float)Lx; mni /= (float)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=2u*K; *X -= mnr; *(X+1) -= mni; }
                        }
                        float sm = 0.0f;
                        for (size_t n=Lx; n>0u; --n, X+=2u*K) { sm += *X**X + *(X+1)**(X+1); }
                        *E++ = sm / (float)Lx;
                        *Y++ = 1.0f; *Y++ = 0.0f;
                    }
                }
            }
        }
    }
    else //L>1u
    {
        float *AC, *A, ar, ai, smr, smi, den, e;
        if (!(AC=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in sig2poly_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in sig2poly_c: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            //Subtract mean
            if (mnz)
            {
                float mnr = 0.0f, mni = 0.0f;
                for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                mnr /= (float)Lx; mni /= (float)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
            }

            //Get ACF
            for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
            {
                smr = smi = 0.0f;
                for (size_t n=Lx-l; n>0u; --n, X+=2)
                {
                    smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                    smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                }
                *AC++ = smr; *AC++ = smi;
            }

            //Normalize ACF
            if (unbiased)
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (float)(Lx-l+1u); *--AC /= (float)(Lx-l+1u); }
            }
            else
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (float)Lx; *--AC /= (float)Lx; }
            }

            //Get poly params and error var (Lev-Durb)
            *Y++ = 1.0f; *Y++ = 0.0f;
            den = *AC**AC + *(AC+1)**(AC+1);
            ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
            ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
            *Y = ar; *(Y+1) = ai;
            e = *AC * (1.0f - (ar*ar+ai*ai));
            AC += 4;
            for (size_t p=1u; p<P; ++p, AC+=2u*p)
            {
                ar = *AC; ai = *(AC+1);
                for (size_t q=p; q>0u; --q, Y+=2)
                {
                    AC -= 2;
                    ar += *AC**Y - *(AC+1)**(Y+1);
                    ai += *(AC+1)**Y + *AC**(Y+1);
                }
                ar /= -e; ai /= -e;
                *Y = ar; *(Y+1) = ai;
                for (size_t q=p; q>0u; --q, A+=2) { Y-=2; *A = *Y; *(A+1) = -*(Y+1); }
                Y += 2u*p;
                for (size_t q=p; q>0u; --q)
                {
                    A -= 2; Y -= 2;
                    *Y += ar**A - ai**(A+1);
                    *(Y+1) += ar**(A+1) + ai**A;
                }
                e *= 1.0f - (ar*ar + ai*ai);
            }
            AC -= 2u*L;
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y+=2u*P, ++E)
                {
                    //Subtract mean
                    if (mnz)
                    {
                        float mnr = 0.0f, mni = 0.0f;
                        for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                        mnr /= (float)Lx; mni /= (float)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
                    }

                    //Get ACF
                    for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
                    {
                        smr = smi = 0.0f;
                        for (size_t n=Lx-l; n>0u; --n, X+=2)
                        {
                            smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                            smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                        }
                        *AC++ = smr; *AC++ = smi;
                    }

                    //Normalize ACF
                    if (unbiased)
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (float)(Lx-l+1u); *--AC /= (float)(Lx-l+1u); }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (float)Lx; *--AC /= (float)Lx; }
                    }

                    //Get poly params and error var (Lev-Durb)
                    *Y++ = 1.0f; *Y++ = 0.0f;
                    den = *AC**AC + *(AC+1)**(AC+1);
                    ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                    ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                    *Y = ar; *(Y+1) = ai;
                    e = *AC * (1.0f - (ar*ar+ai*ai));
                    AC += 4;
                    for (size_t p=1u; p<P; ++p, AC+=2u*p)
                    {
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=p; q>0u; --q, Y+=2)
                        {
                            AC -= 2;
                            ar += *AC**Y - *(AC+1)**(Y+1);
                            ai += *(AC+1)**Y + *AC**(Y+1);
                        }
                        ar /= -e; ai /= -e;
                        *Y = ar; *(Y+1) = ai;
                        for (size_t q=p; q>0u; --q, A+=2) { Y-=2; *A = *Y; *(A+1) = -*(Y+1); }
                        Y += 2u*p;
                        for (size_t q=p; q>0u; --q)
                        {
                            A -= 2; Y -= 2;
                            *Y += ar**A - ai**(A+1);
                            *(Y+1) += ar**(A+1) + ai**A;
                        }
                        e *= 1.0f - (ar*ar + ai*ai);
                    }
                    AC -= 2u*L;
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y+=2, ++E)
                    {
                        //Subtract mean
                        if (mnz)
                        {
                            float mnr = 0.0f, mni = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K) { mnr += *X; mni += *(X+1); }
                            mnr /= (float)Lx; mni /= (float)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=2u*K; *X -= mni; *(X+1) -= mnr; }
                        }

                        //Get ACF
                        for (size_t l=0u; l<L; ++l, X-=2u*K*(Lx-l+1u))
                        {
                            smr = smi = 0.0f;
                            for (size_t n=Lx-l; n>0u; --n, X+=2u*K)
                            {
                                smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
                                smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
                            }
                            *AC++ = smr; *AC++ = smi;
                        }

                        //Normalize ACF
                        if (unbiased)
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (float)(Lx-l+1u); *--AC /= (float)(Lx-l+1u); }
                        }
                        else
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (float)Lx; *--AC /= (float)Lx; }
                        }

                        //Get poly params and error var (Lev-Durb)
                        *Y = 1.0f; *(Y+1) = 0.0f; Y += 2u*K;
                        den = *AC**AC + *(AC+1)**(AC+1);
                        ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                        ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                        *Y = ar; *(Y+1) = ai;
                        e = *AC * (1.0f - (ar*ar+ai*ai));
                        AC += 4;
                        for (size_t p=1u; p<P; ++p, AC+=2u*p)
                        {
                            ar = *AC; ai = *(AC+1);
                            for (size_t q=p; q>0u; --q, Y+=2u*K)
                            {
                                AC -= 2;
                                ar += *AC**Y - *(AC+1)**(Y+1);
                                ai += *(AC+1)**Y + *AC**(Y+1);
                            }
                            ar /= -e; ai /= -e;
                            *Y = ar; *(Y+1) = ai;
                            for (size_t q=p; q>0u; --q, A+=2) { Y-=2u*K; *A = *Y; *(A+1) = -*(Y+1); }
                            Y += 2u*p*K;
                            for (size_t q=p; q>0u; --q)
                            {
                                A -= 2; Y -= 2u*K;
                                *Y += ar**A - ai**(A+1);
                                *(Y+1) += ar**(A+1) + ai**A;
                            }
                            e *= 1.0f - (ar*ar + ai*ai);
                        }
                        AC -= 2u*L;
                        *E = e;
                    }
                }
            }
        }
        free(AC); free(A);
    }

    return 0;
}


int sig2poly_z (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t L, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2poly_z: dim must be in [0 3]\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in sig2ar_z: L must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t P = L - 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2poly_z: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else if (L==1u)
    {
        if (Lx==N)
        {
            if (mnz)
            {
                double mnr = 0.0, mni = 0.0;
                for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                mnr /= (double)Lx; mni /= (double)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
            }
            double sm = 0.0;
            for (size_t n=2u*Lx; n>0u; --n, ++X) { sm += *X * *X; }
            *E = sm / (double)Lx;
            *Y++ = 1.0; *Y = 0.0;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    if (mnz)
                    {
                        double mnr = 0.0, mni = 0.0;
                        for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                        mnr /= (double)Lx; mni /= (double)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
                    }
                    double sm = 0.0;
                    for (size_t n=2u*Lx; n>0u; --n, ++X) { sm += *X * *X; }
                    *E++ = sm / (double)Lx;
                    *Y++ = 1.0; *Y++ = 0.0;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u)
                    {
                        if (mnz)
                        {
                            double mnr = 0.0, mni = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K) { mnr += *X; mni += *(X+1); }
                            mnr /= (double)Lx; mni /= (double)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=2u*K; *X -= mnr; *(X+1) -= mni; }
                        }
                        double sm = 0.0;
                        for (size_t n=Lx; n>0u; --n, X+=2u*K) { sm += *X**X + *(X+1)**(X+1); }
                        *E++ = sm / (double)Lx;
                        *Y++ = 1.0; *Y++ = 0.0;
                    }
                }
            }
        }
    }
    else //L>1u
    {
        double *AC, *A, ar, ai, smr, smi, den, e;
        if (!(AC=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in sig2poly_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in sig2poly_z: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            //Subtract mean
            if (mnz)
            {
                double mnr = 0.0, mni = 0.0;
                for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                mnr /= (double)Lx; mni /= (double)Lx;
                for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
            }

            //Get ACF
            for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
            {
                smr = smi = 0.0;
                for (size_t n=Lx-l; n>0u; --n, X+=2)
                {
                    smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                    smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                }
                *AC++ = smr; *AC++ = smi;
            }

            //Normalize ACF
            if (unbiased)
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (double)(Lx-l+1u); *--AC /= (double)(Lx-l+1u); }
            }
            else
            {
                for (size_t l=L; l>0u; --l) { *--AC /= (double)Lx; *--AC /= (double)Lx; }
            }

            //Get poly params and error var (Lev-Durb)
            *Y++ = 1.0; *Y++ = 0.0;
            den = *AC**AC + *(AC+1)**(AC+1);
            ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
            ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
            *Y = ar; *(Y+1) = ai;
            e = *AC * (1.0 - (ar*ar+ai*ai));
            AC += 4;
            for (size_t p=1u; p<P; ++p, AC+=2u*p)
            {
                ar = *AC; ai = *(AC+1);
                for (size_t q=p; q>0u; --q, Y+=2)
                {
                    AC -= 2;
                    ar += *AC**Y - *(AC+1)**(Y+1);
                    ai += *(AC+1)**Y + *AC**(Y+1);
                }
                ar /= -e; ai /= -e;
                *Y = ar; *(Y+1) = ai;
                for (size_t q=p; q>0u; --q, A+=2) { Y-=2; *A = *Y; *(A+1) = -*(Y+1); }
                Y += 2u*p;
                for (size_t q=p; q>0u; --q)
                {
                    A -= 2; Y -= 2;
                    *Y += ar**A - ai**(A+1);
                    *(Y+1) += ar**(A+1) + ai**A;
                }
                e *= 1.0 - (ar*ar + ai*ai);
            }
            AC -= 2u*L;
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y+=2u*P, ++E)
                {
                    //Subtract mean
                    if (mnz)
                    {
                        double mnr = 0.0, mni = 0.0;
                        for (size_t l=Lx; l>0u; --l) { mnr += *X++; mni += *X++; }
                        mnr /= (double)Lx; mni /= (double)Lx;
                        for (size_t l=Lx; l>0u; --l) { *--X -= mni; *--X -= mnr; }
                    }

                    //Get ACF
                    for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
                    {
                        smr = smi = 0.0;
                        for (size_t n=Lx-l; n>0u; --n, X+=2)
                        {
                            smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
                            smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
                        }
                        *AC++ = smr; *AC++ = smi;
                    }

                    //Normalize ACF
                    if (unbiased)
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (double)(Lx-l+1u); *--AC /= (double)(Lx-l+1u); }
                    }
                    else
                    {
                        for (size_t l=L; l>0u; --l) { *--AC /= (double)Lx; *--AC /= (double)Lx; }
                    }

                    //Get poly params and error var (Lev-Durb)
                    *Y++ = 1.0; *Y++ = 0.0;
                    den = *AC**AC + *(AC+1)**(AC+1);
                    ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                    ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                    *Y = ar; *(Y+1) = ai;
                    e = *AC * (1.0 - (ar*ar+ai*ai));
                    AC += 4;
                    for (size_t p=1u; p<P; ++p, AC+=2u*p)
                    {
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=p; q>0u; --q, Y+=2)
                        {
                            AC -= 2;
                            ar += *AC**Y - *(AC+1)**(Y+1);
                            ai += *(AC+1)**Y + *AC**(Y+1);
                        }
                        ar /= -e; ai /= -e;
                        *Y = ar; *(Y+1) = ai;
                        for (size_t q=p; q>0u; --q, A+=2) { Y-=2; *A = *Y; *(A+1) = -*(Y+1); }
                        Y += 2u*p;
                        for (size_t q=p; q>0u; --q)
                        {
                            A -= 2; Y -= 2;
                            *Y += ar**A - ai**(A+1);
                            *(Y+1) += ar**(A+1) + ai**A;
                        }
                        e *= 1.0 - (ar*ar + ai*ai);
                    }
                    AC -= 2u*L;
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y+=2, ++E)
                    {
                        //Subtract mean
                        if (mnz)
                        {
                            double mnr = 0.0, mni = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K) { mnr += *X; mni += *(X+1); }
                            mnr /= (double)Lx; mni /= (double)Lx;
                            for (size_t l=Lx; l>0u; --l) { X-=2u*K; *X -= mni; *(X+1) -= mnr; }
                        }

                        //Get ACF
                        for (size_t l=0u; l<L; ++l, X-=2u*K*(Lx-l+1u))
                        {
                            smr = smi = 0.0;
                            for (size_t n=Lx-l; n>0u; --n, X+=2u*K)
                            {
                                smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
                                smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
                            }
                            *AC++ = smr; *AC++ = smi;
                        }

                        //Normalize ACF
                        if (unbiased)
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (double)(Lx-l+1u); *--AC /= (double)(Lx-l+1u); }
                        }
                        else
                        {
                            for (size_t l=L; l>0u; --l) { *--AC /= (double)Lx; *--AC /= (double)Lx; }
                        }

                        //Get poly params and error var (Lev-Durb)
                        *Y = 1.0; *(Y+1) = 0.0; Y += 2u*K;
                        den = *AC**AC + *(AC+1)**(AC+1);
                        ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                        ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                        *Y = ar; *(Y+1) = ai;
                        e = *AC * (1.0 - (ar*ar+ai*ai));
                        AC += 4;
                        for (size_t p=1u; p<P; ++p, AC+=2u*p)
                        {
                            ar = *AC; ai = *(AC+1);
                            for (size_t q=p; q>0u; --q, Y+=2u*K)
                            {
                                AC -= 2;
                                ar += *AC**Y - *(AC+1)**(Y+1);
                                ai += *(AC+1)**Y + *AC**(Y+1);
                            }
                            ar /= -e; ai /= -e;
                            *Y = ar; *(Y+1) = ai;
                            for (size_t q=p; q>0u; --q, A+=2) { Y-=2u*K; *A = *Y; *(A+1) = -*(Y+1); }
                            Y += 2u*p*K;
                            for (size_t q=p; q>0u; --q)
                            {
                                A -= 2; Y -= 2u*K;
                                *Y += ar**A - ai**(A+1);
                                *(Y+1) += ar**(A+1) + ai**A;
                            }
                            e *= 1.0 - (ar*ar + ai*ai);
                        }
                        AC -= 2u*L;
                        *E = e;
                    }
                }
            }
        }
        free(AC); free(A);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
