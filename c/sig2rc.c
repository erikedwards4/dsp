//Gets P reflection coeffs (Y), and error variance (E), for each vector in X.
//Works along dim, such that each signal (vector) in X is converted
//to one scalar in E and one vector (of length P) in Y.

//By default, the means of X are NOT subtracted.
//The mnz option allows the mean to be zeroed first for each vec in X.
//However, this also required the removal of the const qualifier for input *X.

//This starts by getting the autocovariance (AC) for lags 0 to P for each vector in X.
//The "biased" version of AC uses N in the denominator,
//whereas the "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower, has worse MSE, and doesn't match FFT estimate.

//After getting the AC, Levinson-Durbin recursion is used to get the RCs (Y)
//and the error variance (E).

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2rc_s (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2rc_d (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2rc_c (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2rc_z (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);


int sig2rc_s (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2rc_s: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in sig2rc_s: P must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2rc_s: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        float *AC, *A1, *A2, a, sm, e;
        if (!(AC=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sig2rc_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(float *)malloc(P*sizeof(float)))) { fprintf(stderr,"error in sig2rc_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in sig2rc_s: problem with malloc. "); perror("malloc"); return 1; }

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

            //Get RCs and error var (Lev-Durb)
            a = -*(AC+1) / *AC;
            *A1 = a; *Y++ = a;
            e = *AC++; e += a * *AC++;
            for (size_t p=1u; p<P-1u; ++p, AC+=p)
            {
                a = *AC;
                for (size_t q=p; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                a /= -e;
                *A1 = a; *Y++ = a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0f - a*a;
            }
            a = *AC;
            for (size_t q=P; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
            A1 -= P;
            a /= -e; *Y = a;
            *E = e * (1.0f-a*a);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Lx, ++E)
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

                    //Get RCs and error var (Lev-Durb)
                    a = -*(AC+1) / *AC;
                    *A1 = a; *Y++ = a;
                    e = *AC++; e += a * *AC++;
                    for (size_t p=1u; p<P-1u; ++p, AC+=p)
                    {
                        a = *AC;
                        for (size_t q=p; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                        a /= -e;
                        *A1 = a; *Y++ = a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0f - a*a;
                    }
                    a = *AC;
                    for (size_t q=P; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                    A1 -= P;
                    a /= -e; *Y++ = a;
                    *E = e * (1.0f-a*a);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*P-K-1u, ++E)
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

                        //Get RCs and error var (Lev-Durb)
                        a = -*(AC+1) / *AC;
                        *A1 = a;
                        *Y = a; Y += K;
                        e = *AC++; e += a * *AC++;
                        for (size_t p=1u; p<P-1u; ++p, AC+=p)
                        {
                            a = *AC;
                            for (size_t q=p; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                            a /= -e; *A1 = a;
                            *Y = a; Y += K;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0f - a*a;
                        }
                        a = *AC;
                        for (size_t q=P; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                        A1 -= P;
                        a /= -e; *Y = a;
                        *E = e * (1.0f-a*a);
                    }
                }
            }
        }
        free(AC); free(A1); free(A2);
    }

    return 0;
}


int sig2rc_d (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2rc_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2rc_d: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        double *AC, *A1, *A2, a, sm, e;
        if (!(AC=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in sig2rc_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(double *)malloc(P*sizeof(double)))) { fprintf(stderr,"error in sig2rc_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in sig2rc_d: problem with malloc. "); perror("malloc"); return 1; }

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

            //Get RCs and error var (Lev-Durb)
            a = -*(AC+1) / *AC;
            *A1 = a; *Y++ = a;
            e = *AC++; e += a * *AC++;
            for (size_t p=1u; p<P-1u; ++p, AC+=p)
            {
                a = *AC;
                for (size_t q=p; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                a /= -e;
                *A1 = a; *Y++ = a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0 - a*a;
            }
            a = *AC;
            for (size_t q=P; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
            A1 -= P;
            a /= -e; *Y = a;
            *E = e * (1.0-a*a);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Lx, ++E)
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

                    //Get RCs and error var (Lev-Durb)
                    a = -*(AC+1) / *AC;
                    *A1 = a; *Y++ = a;
                    e = *AC++; e += a * *AC++;
                    for (size_t p=1u; p<P-1u; ++p, AC+=p)
                    {
                        a = *AC;
                        for (size_t q=p; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                        a /= -e;
                        *A1 = a; *Y++ = a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0 - a*a;
                    }
                    a = *AC;
                    for (size_t q=P; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                    A1 -= P;
                    a /= -e; *Y++ = a;
                    *E = e * (1.0-a*a);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*P-K-1u, ++E)
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

                        //Get RCs and error var (Lev-Durb)
                        a = -*(AC+1) / *AC;
                        *A1 = a;
                        *Y = a; Y += K;
                        e = *AC++; e += a * *AC++;
                        for (size_t p=1u; p<P-1u; ++p, AC+=p)
                        {
                            a = *AC;
                            for (size_t q=p; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                            a /= -e; *A1 = a;
                            *Y = a; Y += K;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0 - a*a;
                        }
                        a = *AC;
                        for (size_t q=P; q>0u; --q, ++A1) { --AC; a += *AC * *A1; }
                        A1 -= P;
                        a /= -e; *Y = a;
                        *E = e * (1.0-a*a);
                    }
                }
            }
        }
        free(AC); free(A1); free(A2);
    }

    return 0;
}


int sig2rc_c (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2rc_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2rc_c: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        float *AC, *A1, *A2, ar, ai, smr, smi, den, e;
        if (!(AC=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in sig2rc_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(float *)malloc(2u*P*sizeof(float)))) { fprintf(stderr,"error in sig2rc_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in sig2rc_c: problem with malloc. "); perror("malloc"); return 1; }

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

            //Get RCs and error var (Lev-Durb)
            den = *AC**AC + *(AC+1)**(AC+1);
            ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
            ai = (*(AC+1)**(AC+2) - *(AC+3)**AC) / den;
            *A1 = ar; *(A1+1) = ai;
            *Y = ar; *(Y+1) = ai;
            e = *AC * (1.0f - (ar*ar+ai*ai));
            AC += 4; Y += 2;
            for (size_t p=1u; p<P-1u; ++p, AC+=2u*p, Y+=2)
            {
                ar = *AC; ai = *(AC+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    AC -= 2;
                    ar += *AC**A1 - *(AC+1)**(A1+1);
                    ai += *(AC+1)**A1 + *AC**(A1+1);
                }
                ar /= -e; ai /= -e;
                *A1 = ar; *(A1+1) = ai;
                *Y = ar; *(Y+1) = ai;
                for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                A1 += 2u*p;
                for (size_t q=p; q>0u; --q)
                {
                    A2 -= 2; A1 -= 2;
                    *A1 += ar**A2 - ai**(A2+1);
                    *(A1+1) += ar**(A2+1) + ai**A2;
                }
                e *= 1.0f - (ar*ar + ai*ai);
            }
            ar = *AC; ai = *(AC+1);
            for (size_t q=P; q>0u; --q, A1+=2)
            {
                AC -= 2;
                ar += *AC**A1 - *(AC+1)**(A1+1);
                ai += *(AC+1)**A1 + *AC**(A1+1);
            }
            A1 -= 2u*P;
            ar /= -e; ai /= -e;
            *Y++ = ar; *Y = ai;
            *E = e * (1.0f-ar*ar-ai*ai);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, ++E)
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

                    //Get RCs and error var (Lev-Durb)
                    den = *AC**AC + *(AC+1)**(AC+1);
                    ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                    ai = (*(AC+1)**(AC+2) - *(AC+3)**AC) / den;
                    *A1 = ar; *(A1+1) = ai;
                    *Y = ar; *(Y+1) = ai;
                    e = *AC * (1.0f - (ar*ar+ai*ai));
                    AC += 4; Y += 2;
                    for (size_t p=1u; p<P-1u; ++p, AC+=2u*p, Y+=2)
                    {
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            AC -= 2;
                            ar += *AC**A1 - *(AC+1)**(A1+1);
                            ai += *(AC+1)**A1 + *AC**(A1+1);
                        }
                        ar /= -e; ai /= -e;
                        *A1 = ar; *(A1+1) = ai;
                        *Y = ar; *(Y+1) = ai;
                        for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                        A1 += 2u*p;
                        for (size_t q=p; q>0u; --q)
                        {
                            A2 -= 2; A1 -= 2;
                            *A1 += ar**A2 - ai**(A2+1);
                            *(A1+1) += ar**(A2+1) + ai**A2;
                        }
                        e *= 1.0f - (ar*ar + ai*ai);
                    }
                    ar = *AC; ai = *(AC+1);
                    for (size_t q=P; q>0u; --q, A1+=2)
                    {
                        AC -= 2;
                        ar += *AC**A1 - *(AC+1)**(A1+1);
                        ai += *(AC+1)**A1 + *AC**(A1+1);
                    }
                    A1 -= 2u*P;
                    ar /= -e; ai /= -e;
                    *Y++ = ar; *Y++ = ai;
                    *E = e * (1.0f-ar*ar-ai*ai);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*P-2u, ++E)
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

                        //Get RCs and error var (Lev-Durb)
                        den = *AC**AC + *(AC+1)**(AC+1);
                        ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                        ai = (*(AC+1)**(AC+2) - *(AC+3)**AC) / den;
                        *A1 = ar; *(A1+1) = ai;
                        *Y = ar; *(Y+1) = ai;
                        e = *AC * (1.0f - (ar*ar+ai*ai));
                        AC += 4; Y += 2u*K;
                        for (size_t p=1u; p<P-1u; ++p, AC+=2u*p, Y+=2u*K)
                        {
                            ar = *AC; ai = *(AC+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                AC -= 2;
                                ar += *AC**A1 - *(AC+1)**(A1+1);
                                ai += *(AC+1)**A1 + *AC**(A1+1);
                            }
                            ar /= -e; ai /= -e;
                            *A1 = ar; *(A1+1) = ai;
                            *Y = ar; *(Y+1) = ai;
                            for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                            A1 += 2u*p;
                            for (size_t q=p; q>0u; --q)
                            {
                                A2 -= 2; A1 -= 2;
                                *A1 += ar**A2 - ai**(A2+1);
                                *(A1+1) += ar**(A2+1) + ai**A2;
                            }
                            e *= 1.0f - (ar*ar + ai*ai);
                        }
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=P; q>0u; --q, A1+=2)
                        {
                            AC -= 2;
                            ar += *AC**A1 - *(AC+1)**(A1+1);
                            ai += *(AC+1)**A1 + *AC**(A1+1);
                        }
                        A1 -= 2u*P;
                        ar /= -e; ai /= -e;
                        *Y = ar; *(Y+1) = ai; Y += 2u*K;
                        *E = e * (1.0f-ar*ar-ai*ai);
                    }
                }
            }
        }
        free(AC); free(A1); free(A2);
    }

    return 0;
}


int sig2rc_z (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2rc_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2rc_z: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        double *AC, *A1, *A2, ar, ai, smr, smi, den, e;
        if (!(AC=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in sig2rc_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(double *)malloc(2u*P*sizeof(double)))) { fprintf(stderr,"error in sig2rc_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in sig2rc_z: problem with malloc. "); perror("malloc"); return 1; }

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

            //Get RCs and error var (Lev-Durb)
            den = *AC**AC + *(AC+1)**(AC+1);
            ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
            ai = (*(AC+1)**(AC+2) - *(AC+3)**AC) / den;
            *A1 = ar; *(A1+1) = ai;
            *Y = ar; *(Y+1) = ai;
            e = *AC * (1.0 - (ar*ar+ai*ai));
            AC += 4; Y += 2;
            for (size_t p=1u; p<P-1u; ++p, AC+=2u*p, Y+=2)
            {
                ar = *AC; ai = *(AC+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    AC -= 2;
                    ar += *AC**A1 - *(AC+1)**(A1+1);
                    ai += *(AC+1)**A1 + *AC**(A1+1);
                }
                ar /= -e; ai /= -e;
                *A1 = ar; *(A1+1) = ai;
                *Y = ar; *(Y+1) = ai;
                for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                A1 += 2u*p;
                for (size_t q=p; q>0u; --q)
                {
                    A2 -= 2; A1 -= 2;
                    *A1 += ar**A2 - ai**(A2+1);
                    *(A1+1) += ar**(A2+1) + ai**A2;
                }
                e *= 1.0 - (ar*ar + ai*ai);
            }
            ar = *AC; ai = *(AC+1);
            for (size_t q=P; q>0u; --q, A1+=2)
            {
                AC -= 2;
                ar += *AC**A1 - *(AC+1)**(A1+1);
                ai += *(AC+1)**A1 + *AC**(A1+1);
            }
            A1 -= 2u*P;
            ar /= -e; ai /= -e;
            *Y++ = ar; *Y = ai;
            *E = e * (1.0-ar*ar-ai*ai);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, ++E)
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

                    //Get RCs and error var (Lev-Durb)
                    den = *AC**AC + *(AC+1)**(AC+1);
                    ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                    ai = (*(AC+1)**(AC+2) - *(AC+3)**AC) / den;
                    *A1 = ar; *(A1+1) = ai;
                    *Y = ar; *(Y+1) = ai;
                    e = *AC * (1.0 - (ar*ar+ai*ai));
                    AC += 4; Y += 2;
                    for (size_t p=1u; p<P-1u; ++p, AC+=2u*p, Y+=2)
                    {
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            AC -= 2;
                            ar += *AC**A1 - *(AC+1)**(A1+1);
                            ai += *(AC+1)**A1 + *AC**(A1+1);
                        }
                        ar /= -e; ai /= -e;
                        *A1 = ar; *(A1+1) = ai;
                        *Y = ar; *(Y+1) = ai;
                        for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                        A1 += 2u*p;
                        for (size_t q=p; q>0u; --q)
                        {
                            A2 -= 2; A1 -= 2;
                            *A1 += ar**A2 - ai**(A2+1);
                            *(A1+1) += ar**(A2+1) + ai**A2;
                        }
                        e *= 1.0 - (ar*ar + ai*ai);
                    }
                    ar = *AC; ai = *(AC+1);
                    for (size_t q=P; q>0u; --q, A1+=2)
                    {
                        AC -= 2;
                        ar += *AC**A1 - *(AC+1)**(A1+1);
                        ai += *(AC+1)**A1 + *AC**(A1+1);
                    }
                    A1 -= 2u*P;
                    ar /= -e; ai /= -e;
                    *Y++ = ar; *Y++ = ai;
                    *E = e * (1.0-ar*ar-ai*ai);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*P-2u, ++E)
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

                        //Get RCs and error var (Lev-Durb)
                        den = *AC**AC + *(AC+1)**(AC+1);
                        ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                        ai = (*(AC+1)**(AC+2) - *(AC+3)**AC) / den;
                        *A1 = ar; *(A1+1) = ai;
                        *Y = ar; *(Y+1) = ai;
                        e = *AC * (1.0 - (ar*ar+ai*ai));
                        AC += 4; Y += 2u*K;
                        for (size_t p=1u; p<P-1u; ++p, AC+=2u*p, Y+=2u*K)
                        {
                            ar = *AC; ai = *(AC+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                AC -= 2;
                                ar += *AC**A1 - *(AC+1)**(A1+1);
                                ai += *(AC+1)**A1 + *AC**(A1+1);
                            }
                            ar /= -e; ai /= -e;
                            *A1 = ar; *(A1+1) = ai;
                            *Y = ar; *(Y+1) = ai;
                            for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                            A1 += 2u*p;
                            for (size_t q=p; q>0u; --q)
                            {
                                A2 -= 2; A1 -= 2;
                                *A1 += ar**A2 - ai**(A2+1);
                                *(A1+1) += ar**(A2+1) + ai**A2;
                            }
                            e *= 1.0 - (ar*ar + ai*ai);
                        }
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=P; q>0u; --q, A1+=2)
                        {
                            AC -= 2;
                            ar += *AC**A1 - *(AC+1)**(A1+1);
                            ai += *(AC+1)**A1 + *AC**(A1+1);
                        }
                        A1 -= 2u*P;
                        ar /= -e; ai /= -e;
                        *Y = ar; *(Y+1) = ai; Y += 2u*K;
                        *E = e * (1.0-ar*ar-ai*ai);
                    }
                }
            }
        }
        free(AC); free(A1); free(A2);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
