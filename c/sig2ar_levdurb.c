//Gets P autoregressive (AR) params (Y), and error variance (V), for each vector in X.
//Works along dim, such that each signal (vector) in X is converted
//to one scalar in V and one vector (of length P) in Y.

//By default, the means of X are NOT subtracted.
//The mnz option allows the mean to be zeroed first for each vec in X.
//However, this also required the removal of the const qualifier for input *X.

//This starts by getting the autocovariance (AC) for lags 0 to P for each vector in X.
//The "biased" version of AC uses N in the denominator,
//whereas the "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower, has worse MSE, and doesn't match FFT estimate.

//After getting the AC, Levinson-Durbin (levdurb) recursion is used to get the AR params (Y)
//and the error variance (V).

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2ar_levdurb_s (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2ar_levdurb_d (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2ar_levdurb_c (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2ar_levdurb_z (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);


int sig2ar_levdurb_s (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2ar_levdurb_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (P>=Lx) { fprintf(stderr,"error in sig2ar_levdurb_s: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else
    {
        float *AC, *A2, a, e, sm;
        if (!(AC=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sig2ar_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in sig2ar_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            if (mnz)
            {
                float mn = 0.0f;
                for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
                mn /= (float)Lx;
                for (size_t l=0u; l<Lx; ++l) { *--X -= mn; }
            }

            for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++AC)
            {
                sm = 0.0f;
                for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
                *AC = sm;
            }
            AC -= L;

            if (unbiased)
            {
                for (size_t l=0u; l<L; ++l, ++AC) { *AC /= (float)(Lx-l); }
                AC -= L;
            }
            for (size_t l=0u; l<L; ++l, ++AC) { fprintf(stderr,"AC[%lu] = %g \n",l,(double)AC[l]); }
            AC -= L;

            a = -*(AC+1) / *AC;
            *Y = -a;
            e = *AC++; e += a * *AC++;
            for (size_t p=1u; p<P; ++p, AC+=p)
            {
                a = *AC;
                for (size_t q=0u; q<p; ++q, ++Y) { --AC; a -= *AC * *Y; }
                a /= -e;
                *Y = -a;
                for (size_t q=0u; q<p; ++q, ++A2) { --Y; *A2 = *Y; }
                Y += p;
                for (size_t q=0u; q<p; ++q) { --A2; --Y; *Y += a * *A2; }
                e *= 1.0f - a*a;
            }
            AC -= P+1u;
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=P, Y+=P, ++E)
                {
                    if (mnz)
                    {
                        float mn = 0.0f;
                        for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
                        mn /= (float)Lx;
                        for (size_t l=0u; l<Lx; ++l) { *--X -= mn; }
                    }

                    for (size_t l=0u; l<P; ++l, X-=Lx-l+1u, ++AC)
                    {
                        sm = 0.0f;
                        for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
                        *AC = sm;
                    }
                    sm = 0.0f;
                    for (size_t n=0u; n<Lx-P; ++n, ++X) { sm += *X * *(X+P); }
                    *AC = sm; AC -= P;

                    if (unbiased)
                    {
                        for (size_t l=0u; l<L; ++l, ++AC) { *AC /= (float)(Lx-l); }
                        AC -= L;
                    }

                    a = -*(AC+1) / *AC;
                    *Y = -a;
                    e = *AC++; e += a * *AC++;
                    for (size_t p=1u; p<P; ++p, AC+=p)
                    {
                        a = *AC;
                        for (size_t q=0u; q<p; ++q, ++Y) { --AC; a -= *AC * *Y; }
                        a /= -e;
                        *Y = -a;
                        for (size_t q=0u; q<p; ++q, ++A2) { --Y; *A2 = *Y; }
                        Y += p;
                        for (size_t q=0u; q<p; ++q) { --A2; --Y; *Y += a * *A2; }
                        e *= 1.0f - a*a;
                    }
                    *E = e;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(P-1u))
                {
                    for (size_t b=0u; b<B; ++b, ++X, ++Y, ++E)
                    {
                        if (mnz)
                        {
                            float mn = 0.0f;
                            for (size_t l=0u; l<Lx; ++l, X+=K) { mn += *X; }
                            mn /= (float)Lx;
                            for (size_t l=0u; l<Lx; ++l) { X-=K; *X -= mn; }
                        }

                        for (size_t l=0u; l<L; ++l, X-=K*(Lx-l+1u), ++AC)
                        {
                            sm = 0.0f;
                            for (size_t n=0u; n<Lx-l; ++n, X+=K) { sm += *X * *(X+l*K); }
                            *AC = sm;
                        }
                        AC -= L;

                        if (unbiased)
                        {
                            for (size_t l=0u; l<L; ++l, ++AC) { *AC /= (float)(Lx-l); }
                            AC -= L;
                        }

                        a = -*(AC+1) / *AC;
                        *Y = -a;
                        e = *AC; ++AC;
                        e += a * *AC; ++AC;
                        for (size_t p=1u; p<P; ++p, AC+=p)
                        {
                            a = *AC;
                            for (size_t q=0u; q<p; ++q, Y+=K) { --AC; a -= *AC * *Y; }
                            a /= -e;
                            *Y = -a;
                            for (size_t q=0u; q<p; ++q, ++A2) { Y-=K; *A2 = *Y; }
                            Y += p*K;
                            for (size_t q=0u; q<p; ++q) { --A2; Y-=K; *Y += a * *A2; }
                            e *= 1.0f - a*a;
                        }
                        *E = e;
                    }
                }
            }
        }
        free(AC); free(A2);
    }

    return 0;
}


// int sig2ar_levdurb_d (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
// {
//     if (dim>3u) { fprintf(stderr,"error in sig2ar_levdurb_d: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (P>=Lx) { fprintf(stderr,"error in sig2ar_levdurb_d: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

//     if (N==0u) {}
//     else
//     {
//         if (Lx==N)
//         {
//             if (mnz)
//             {
//                 double mn = 0.0;
//                 for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
//                 mn /= (double)Lx;
//                 for (size_t l=0u; l<Lx; ++l) { *--X -= mn; }
//             }

//             if (L>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
//             {
//                 double x;
//                 for (size_t l=0u; l<L; ++l, ++Y) { *Y = 0.0; }
//                 Y -= L;
//                 for (size_t n=0u; n<N; ++n, ++X)
//                 {
//                     x = *X;
//                     for (size_t l=0u; l<L && l<N-n; ++l) { Y[l] += x * X[l]; }
//                 }
//             }
//             else
//             {
//                 double sm;
//                 for (size_t l=0u; l<L; ++l, X-=Lx-l+1u, ++Y)
//                 {
//                     sm = 0.0;
//                     for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
//                     *Y = sm;
//                 }
//                 Y -= L;
//             }

//             if (corr)
//             {
//                 const double y0 = *Y;
//                 *Y++ = 1.0;
//                 for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
//             }
//             else if (unbiased)
//             {
//                 for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)(Lx-l); }
//             }
//             else //xcorr leaves this blank
//             {
//                 for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)Lx; }
//             }
//         }
//         else
//         {
//             const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;
//             double sm;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=L-1u)
//                 {
//                     if (mnz)
//                     {
//                         double mn = 0.0;
//                         for (size_t l=0u; l<Lx; ++l, ++X) { mn += *X; }
//                         mn /= (double)Lx;
//                         for (size_t l=0u; l<Lx; ++l) { *--X -= mn; }
//                     }

//                     for (size_t l=0u; l<L-1u; ++l, X-=Lx-l+1u, ++Y)
//                     {
//                         sm = 0.0;
//                         for (size_t n=0u; n<Lx-l; ++n, ++X) { sm += *X * *(X+l); }
//                         *Y = sm;
//                     }
//                     sm = 0.0;
//                     for (size_t n=0u; n<Lx-L+1u; ++n, ++X) { sm += *X * *(X+L-1u); }
//                     *Y = sm; Y -= L-1u;

//                     if (corr)
//                     {
//                         const double y0 = *Y;
//                         *Y++ = 1.0;
//                         for (size_t l=1u; l<L; ++l, ++Y) { *Y /= y0; }
//                     }
//                     else if (unbiased)
//                     {
//                         for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)(Lx-l); }
//                     }
//                     else
//                     {
//                         for (size_t l=0u; l<L; ++l, ++Y) { *Y /= (double)Lx; }
//                     }
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(L-1u))
//                 {
//                     for (size_t b=0u; b<B; ++b, ++X, Y-=K*L-1u)
//                     {
//                         if (mnz)
//                         {
//                             double mn = 0.0;
//                             for (size_t l=0u; l<Lx; ++l, X+=K) { mn += *X; }
//                             mn /= (double)Lx;
//                             for (size_t l=0u; l<Lx; ++l) { X-=K; *X -= mn; }
//                         }

//                         for (size_t l=0u; l<L; ++l, X-=K*(Lx-l+1u), Y+=K)
//                         {
//                             sm = 0.0;
//                             for (size_t n=0u; n<Lx-l; ++n, X+=K) { sm += *X * *(X+l*K); }
//                             *Y = sm;
//                         }
//                         Y -= K*L;

//                         if (corr)
//                         {
//                             const double y0 = *Y;
//                             *Y = 1.0; Y += K;
//                             for (size_t l=1u; l<L; ++l, Y+=K) { *Y /= y0; }
//                         }
//                         else if (unbiased)
//                         {
//                             for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (double)(Lx-l); }
//                         }
//                         else
//                         {
//                             for (size_t l=0u; l<L; ++l, Y+=K) { *Y /= (double)Lx; }
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return 0;
// }


// int sig2ar_levdurb_c (float *Y, float *E, float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
// {
//     if (dim>3u) { fprintf(stderr,"error in sig2ar_levdurb_c: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (P>=Lx) { fprintf(stderr,"error in sig2ar_levdurb_c: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

//     if (N==0u) {}
//     else
//     {
//         if (Lx==N)
//         {
//             if (mnz)
//             {
//                 float mnr = 0.0f, mni = 0.0f;
//                 for (size_t l=0u; l<Lx; ++l) { mnr += *X++; mni += *X++; }
//                 mnr /= (float)Lx; mni /= (float)Lx;
//                 for (size_t l=0u; l<Lx; ++l) { *--X -= mni; *--X -= mnr; }
//             }

//             if (L>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
//             {
//                 float xr, xi;
//                 for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y = 0.0f; }
//                 Y -= 2u*L;
//                 for (size_t n=0u; n<N; ++n)
//                 {
//                     xr = *X++; xi = *X++;
//                     for (size_t l=0u; l<L && l<N-n; ++l)
//                     {
//                         Y[2u*l] += xr*X[2u*l] + xi*X[2u*l+1u];
//                         Y[2u*l+1u] += xi*X[2u*l] - xr*X[2u*l+1u];
//                     }
//                 }
//             }
//             else
//             {
//                 float smr, smi;
//                 for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
//                 {
//                     smr = smi = 0.0f;
//                     for (size_t n=0u; n<Lx-l; ++n, X+=2)
//                     {
//                         smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
//                         smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
//                     }
//                     *Y++ = smr; *Y++ = smi;
//                 }
//                 Y -= 2u*L;
//             }
            
//             if (corr)
//             {
//                 const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
//                 float yr, yi;
//                 *Y++ = 1.0f; *Y++ = 0.0f;
//                 for (size_t l=1u; l<L; ++l)
//                 {
//                     yr = *Y; yi = *(Y+1);
//                     *Y++ = (yr*y0r+yi*y0i) / y0a;
//                     *Y++ = (yi*y0r-yr*y0i) / y0a;
//                 }
//             }
//             else if (unbiased)
//             {
//                 for (size_t l=0u; l<L; ++l) { *Y++ /= (float)(Lx-l); *Y++ /= (float)(Lx-l); }
//             }
//             else
//             {
//                 for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (float)Lx; }
//             }
//         }
//         else
//         {
//             const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;
//             float smr, smi;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=2u*L-2u)
//                 {
//                     if (mnz)
//                     {
//                         float mnr = 0.0f, mni = 0.0f;
//                         for (size_t l=0u; l<Lx; ++l) { mnr += *X++; mni += *X++; }
//                         mnr /= (float)Lx; mni /= (float)Lx;
//                         for (size_t l=0u; l<Lx; ++l) { *--X -= mni; *--X -= mnr; }
//                     }

//                     for (size_t l=0u; l<L-1u; ++l, X-=2u*(Lx-l+1u))
//                     {
//                         smr = smi = 0.0f;
//                         for (size_t n=0u; n<Lx-l; ++n, X+=2)
//                         {
//                             smr += *X**(X+2u*l) + *(X+1u)**(X+2u*l+1u);
//                             smi += *(X+1u)**(X+2u*l) - *X**(X+2u*l+1u);
//                         }
//                         *Y++ = smr; *Y++ = smi;
//                     }
//                     smr = smi = 0.0f;
//                     for (size_t n=0u; n<Lx-L+1u; ++n, X+=2)
//                     {
//                         smr += *X**(X+2u*L-2u) + *(X+1u)**(X+2u*L-1u);
//                         smi += *(X+1u)**(X+2u*L-2u) - *X**(X+2u*L-1u);
//                     }
//                     *Y = smr; *(Y+1) = smi; Y -= 2u*L-2u;

//                     if (corr)
//                     {
//                         const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
//                         float yr, yi;
//                         *Y++ = 1.0f; *Y++ = 0.0f;
//                         for (size_t l=1u; l<L; ++l)
//                         {
//                             yr = *Y; yi = *(Y+1);
//                             *Y++ = (yr*y0r+yi*y0i) / y0a;
//                             *Y++ = (yi*y0r-yr*y0i) / y0a;
//                         }
//                     }
//                     else if (unbiased)
//                     {
//                         for (size_t l=0u; l<L; ++l) { *Y++ /= (float)(Lx-l); *Y++ /= (float)(Lx-l); }
//                     }
//                     else
//                     {
//                         for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (float)Lx; }
//                     }
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(L-1u))
//                 {
//                     for (size_t b=0u; b<B; ++b, X+=2, Y-=2u*K*L-2u)
//                     {
//                         if (mnz)
//                         {
//                             float mnr = 0.0f, mni = 0.0f;
//                             for (size_t l=0u; l<Lx; ++l, X+=2u*K) { mnr += *X; mni += *(X+1); }
//                             mnr /= (float)Lx; mni /= (float)Lx;
//                             for (size_t l=0u; l<Lx; ++l) { X-=2u*K; *X -= mnr; *(X+1) -= mni; }
//                         }

//                         for (size_t l=0u; l<L; ++l, X-=2u*K*(Lx-l+1u), Y+=2u*K)
//                         {
//                             smr = smi = 0.0f;
//                             for (size_t n=0u; n<Lx-l; ++n, X+=2u*K)
//                             {
//                                 smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
//                                 smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
//                             }
//                             *Y = smr; *(Y+1) = smi;
//                         }
//                         Y -= 2u*K*L;

//                         if (corr)
//                         {
//                             const float y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
//                             float yr, yi;
//                             *Y = 1.0f; *(Y+1) = 0.0f; Y += 2u*K;
//                             for (size_t l=1u; l<L; ++l, Y+=2u*K)
//                             {
//                                 yr = *Y; yi = *(Y+1);
//                                 *Y = (yr*y0r+yi*y0i) / y0a;
//                                 *(Y+1) = (yi*y0r-yr*y0i) / y0a;
//                             }
//                         }
//                         else if (unbiased)
//                         {
//                             for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (float)(Lx-l); *(Y+1) /= (float)(Lx-l); }
//                         }
//                         else
//                         {
//                             for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (float)Lx; *(Y+1) /= (float)Lx; }
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return 0;
// }


// int sig2ar_levdurb_z (double *Y, double *E, double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
// {
//     if (dim>3u) { fprintf(stderr,"error in sig2ar_levdurb_z: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (P>=Lx) { fprintf(stderr,"error in sig2ar_levdurb_z: P (polynomial order) must be < Lx (length of vecs in X)\n"); return 1; }

//     if (N==0u) {}
//     else
//     {
//         if (Lx==N)
//         {
//             if (mnz)
//             {
//                 double mnr = 0.0, mni = 0.0;
//                 for (size_t l=0u; l<Lx; ++l) { mnr += *X++; mni += *X++; }
//                 mnr /= (double)Lx; mni /= (double)Lx;
//                 for (size_t l=0u; l<Lx; ++l) { *--X -= mni; *--X -= mnr; }
//             }

//             if (L>850u && Lx>80000u) //I also leave this solution since it is closer to real-time use
//             {
//                 double xr, xi;
//                 for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y = 0.0; }
//                 Y -= 2u*L;
//                 for (size_t n=0u; n<N; ++n)
//                 {
//                     xr = *X++; xi = *X++;
//                     for (size_t l=0u; l<L && l<N-n; ++l)
//                     {
//                         Y[2u*l] += xr*X[2u*l] + xi*X[2u*l+1u];
//                         Y[2u*l+1u] += xi*X[2u*l] - xr*X[2u*l+1u];
//                     }
//                 }
//             }
//             else
//             {
//                 double smr, smi;
//                 for (size_t l=0u; l<L; ++l, X-=2u*(Lx-l+1u))
//                 {
//                     smr = smi = 0.0;
//                     for (size_t n=0u; n<Lx-l; ++n, X+=2)
//                     {
//                         smr += *X**(X+2u*l) + *(X+1)**(X+2u*l+1u);
//                         smi += *(X+1)**(X+2u*l) - *X**(X+2u*l+1u);
//                     }
//                     *Y++ = smr; *Y++ = smi;
//                 }
//                 Y -= 2u*L;
//             }
            
//             if (corr)
//             {
//                 const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
//                 double yr, yi;
//                 *Y++ = 1.0; *Y++ = 0.0;
//                 for (size_t l=1u; l<L; ++l)
//                 {
//                     yr = *Y; yi = *(Y+1);
//                     *Y++ = (yr*y0r+yi*y0i) / y0a;
//                     *Y++ = (yi*y0r-yr*y0i) / y0a;
//                 }
//             }
//             else if (unbiased)
//             {
//                 for (size_t l=0u; l<L; ++l) { *Y++ /= (double)(Lx-l); *Y++ /= (double)(Lx-l); }
//             }
//             else
//             {
//                 for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (double)Lx; }
//             }
//         }
//         else
//         {
//             const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;
//             double smr, smi;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=2u*L-2u)
//                 {
//                     if (mnz)
//                     {
//                         double mnr = 0.0, mni = 0.0;
//                         for (size_t l=0u; l<Lx; ++l) { mnr += *X++; mni += *X++; }
//                         mnr /= (double)Lx; mni /= (double)Lx;
//                         for (size_t l=0u; l<Lx; ++l) { *--X -= mni; *--X -= mnr; }
//                     }

//                     for (size_t l=0u; l<L-1u; ++l, X-=2u*(Lx-l+1u))
//                     {
//                         smr = smi = 0.0;
//                         for (size_t n=0u; n<Lx-l; ++n, X+=2)
//                         {
//                             smr += *X**(X+2u*l) + *(X+1u)**(X+2u*l+1u);
//                             smi += *(X+1u)**(X+2u*l) - *X**(X+2u*l+1u);
//                         }
//                         *Y++ = smr; *Y++ = smi;
//                     }
//                     smr = smi = 0.0;
//                     for (size_t n=0u; n<Lx-L+1u; ++n, X+=2)
//                     {
//                         smr += *X**(X+2u*L-2u) + *(X+1u)**(X+2u*L-1u);
//                         smi += *(X+1u)**(X+2u*L-2u) - *X**(X+2u*L-1u);
//                     }
//                     *Y = smr; *(Y+1) = smi; Y -= 2u*L-2u;

//                     if (corr)
//                     {
//                         const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
//                         double yr, yi;
//                         *Y++ = 1.0; *Y++ = 0.0;
//                         for (size_t l=1u; l<L; ++l)
//                         {
//                             yr = *Y; yi = *(Y+1);
//                             *Y++ = (yr*y0r+yi*y0i) / y0a;
//                             *Y++ = (yi*y0r-yr*y0i) / y0a;
//                         }
//                     }
//                     else if (unbiased)
//                     {
//                         for (size_t l=0u; l<L; ++l) { *Y++ /= (double)(Lx-l); *Y++ /= (double)(Lx-l); }
//                     }
//                     else
//                     {
//                         for (size_t l=0u; l<2u*L; ++l, ++Y) { *Y /= (double)Lx; }
//                     }
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(L-1u))
//                 {
//                     for (size_t b=0u; b<B; ++b, X+=2, Y-=2u*K*L-2u)
//                     {
//                         if (mnz)
//                         {
//                             double mnr = 0.0, mni = 0.0;
//                             for (size_t l=0u; l<Lx; ++l, X+=2u*K) { mnr += *X; mni += *(X+1); }
//                             mnr /= (double)Lx; mni /= (double)Lx;
//                             for (size_t l=0u; l<Lx; ++l) { X-=2u*K; *X -= mnr; *(X+1) -= mni; }
//                         }

//                         for (size_t l=0u; l<L; ++l, X-=2u*K*(Lx-l+1u), Y+=2u*K)
//                         {
//                             smr = smi = 0.0;
//                             for (size_t n=0u; n<Lx-l; ++n, X+=2u*K)
//                             {
//                                 smr += *X**(X+2u*l*K) + *(X+1)**(X+2u*l*K+1u);
//                                 smi += *(X+1)**(X+2u*l*K) - *X**(X+2u*l*K+1u);
//                             }
//                             *Y = smr; *(Y+1) = smi;
//                         }
//                         Y -= 2u*K*L;

//                         if (corr)
//                         {
//                             const double y0r = *Y, y0i = *(Y+1), y0a = y0r*y0r + y0i*y0i;
//                             double yr, yi;
//                             *Y = 1.0; *(Y+1) = 0.0; Y += 2u*K;
//                             for (size_t l=1u; l<L; ++l, Y+=2u*K)
//                             {
//                                 yr = *Y; yi = *(Y+1);
//                                 *Y = (yr*y0r+yi*y0i) / y0a;
//                                 *(Y+1) = (yi*y0r-yr*y0i) / y0a;
//                             }
//                         }
//                         else if (unbiased)
//                         {
//                             for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (double)(Lx-l); *(Y+1) /= (double)(Lx-l); }
//                         }
//                         else
//                         {
//                             for (size_t l=0u; l<L; ++l, Y+=2u*K) { *Y /= (double)Lx; *(Y+1) /= (double)Lx; }
//                         }
//                     }
//                 }
//             }
//         }
//     }

//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
