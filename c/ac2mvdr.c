//Gets minimum-variance distortionless response (MVDR) starting from the autocorrelation (AC) function.
//Starts with a Levinson-Durbin recursion from the AC values.
//Then does FFT according to Marple (use only real part of output, so I use FFTW_R2HC).
//Also note that Marple used the Burg method (here I use Levinson-Durbin recursion from the AC).

//Input F is the usual nfft/2 + 1, and is the size of output MVDR along dimension dim.
//Marple [p. 359] uses a large nfft of fixed size 4096, and requires nfft is a power of 2.
//I believe that nfft must be > 2*P, where P = nlags-1 is the AR model order.
//Thus, set nfft to nextpow2(2*Lx).
//Then take usual F = nfft/2 + 1, and set MVDR size to be FxC or RxF.
//The output is then on a linear frequency scale starting at 0 Hz.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ac2mvdr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const float preg);
int ac2mvdr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const double preg);
int ac2mvdr_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const float preg);
int ac2mvdr_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const double preg);


int ac2mvdr_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const float preg)
{
    if (dim>3u) { fprintf(stderr,"error in ac2mvdr_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (F<Lx) { fprintf(stderr,"error in ac2mvdr_s: F must be >= Lx (length of vecs in X)\n"); return 1; }
    if (F%2u!=1u) { fprintf(stderr,"error in ac2mvdr_s: F must be odd (since F=nfft/2+1)\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        float *A1, *A2, a, e;
        if (!(A1=(float *)malloc(P*sizeof(float)))) { fprintf(stderr,"error in ac2mvdr_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2mvdr_s: problem with malloc. "); perror("malloc"); return 1; }
        
        //Initialize FFT
        float *X1, *Y1;
        fftwf_plan fplan;
        size_t nfft = 1u;
        while (nfft<F) { nfft *= 2u; }  //assumes F=nfft/2+1, where nfft is a power of 2
        X1 = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(nfft*sizeof(float));
        fplan = fftwf_plan_r2r_1d((int)nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
        if (!fplan) { fprintf(stderr,"error in ac2mvdr_s: problem creating fftw plan"); return 1; }
        for (size_t nf=0u; nf<nfft; ++nf) { X1[nf] = 0.0f; }

        if (Lx==N)
        {
            //Levinson-Durbin from AC
            a = -*(X+1) / *X;
            *A1 = a;
            e = *X++; e += a * *X++;
            for (size_t p=1u; p<P; ++p, X+=p)
            {
                a = *X;
                for (size_t q=0u; q<p; ++q, ++A1) { --X; a += *X * *A1; }
                a /= -e;
                *A1 = a;
                for (size_t q=0u; q<p; ++q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=0u; q<p; ++q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0f - a*a;
            }
            for (size_t p=0u; p<P; ++p) { fprintf(stderr,"A1[%lu]=%g\n",p,(double)A1[p]); }

            //Enter mus directly into X1
            //(but this can't be right, because A1[q+p] can exceed P, so check again against Marple)
            for (size_t p=0u; p<=P; ++p)
            {
                for (size_t q=0u; q<=P-p; ++q)
                {
                    X1[p] += (float)(P+1u-p-2u*q) * A1[q] * A1[q+p];
                }
                if (p>0u) { X1[nfft-p] = X1[p]; }
            }
            for (size_t nf=0u; nf<nfft; ++nf) { fprintf(stderr,"X1[%lu]=%g\n",nf,(double)X1[nf]); }

            //Marple FFT part (but may need to scale by fs and/or nfft)
            fftwf_execute(fplan);
            
            //Power (from fftw half-complex)
            //(This was not part of my original code, so check again vs. Marple)
            for (size_t f=0u; f<F; ++f, ++Y1, ++Y) { *Y = *Y1 * *Y1; }
            Y -= 2u;
            for (size_t f=1u; f<F-1u; ++f, ++Y1, --Y) { *Y += *Y1 * *Y1; }
            Y1 -= nfft;

            //MVDR
            for (size_t f=0u; f<F; ++f, ++Y) { *Y = e / (*Y+preg); }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=Lx-F, ++Y)
                {
                    a = -*(X+1) / *X;
                    *A1 = a;
                    e = *X++; e += a * *X++;
                    for (size_t p=1u; p<P; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=0u; q<p; ++q, ++A1) { --X; a += *X * *A1; }
                        a /= -e;
                        *A1 = a;
                        for (size_t q=0u; q<p; ++q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=0u; q<p; ++q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0f - a*a;
                    }
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*F-1u, Y-=K*P-1u)
                    {
                        a = -*(X+K) / *X;
                        *A1 = a;
                        e = *X; X += K;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=0u; q<p; ++q, ++A1) { X-=K; a += *X * *A1; }
                            a /= -e;
                            *A1 = a;
                            for (size_t q=0u; q<p; ++q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=0u; q<p; ++q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0f - a*a;
                        }
                    }
                }
            }
        }
        free(A1); free(A2);
        fftwf_destroy_plan(fplan); fftwf_free(X1); fftwf_free(Y1);
    }
	
	return 0;
}

// int ac2mvdr_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const double preg)
// {
// 	if (dim>3u) { fprintf(stderr,"error in ac2mvdr_d: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

//     if (N==0u || Lx<2u) {}
//     else if (Lx==2u)
//     {
//         if (Lx==N) { *Y++ = -*(X+1) / *X; }
//         else
//         {
//             const size_t F = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=2) { *Y++ = -*(X+1) / *X; }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=B)
//                 {
//                     for (size_t b=0u; b<B; ++b, ++X, ++Y) { *Y = -*(X+K) / *X; }
//                 }
//             }
//         }
//     }
//     else
//     {
//         const size_t P = Lx - 1u;
//         double *A1, *A2, a, e;
//         if (!(A1=(double *)malloc(P*sizeof(double)))) { fprintf(stderr,"error in ac2mvdr_d: problem with malloc. "); perror("malloc"); return 1; }
//         if (!(A2=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2mvdr_d: problem with malloc. "); perror("malloc"); return 1; }
        
//         if (Lx==N)
//         {
//             a = -*(X+1) / *X;
//             *A1 = a; *Y++ = a;
//             e = *X++; e += a * *X++;
//             for (size_t p=1u; p<P-1u; ++p, X+=p)
//             {
//                 a = *X;
//                 for (size_t q=0u; q<p; ++q, ++A1) { --X; a += *X * *A1; }
//                 a /= -e;
//                 *A1 = a; *Y++ = a;
//                 for (size_t q=0u; q<p; ++q, ++A2) { --A1; *A2 = *A1; }
//                 A1 += p;
//                 for (size_t q=0u; q<p; ++q) { --A2; --A1; *A1 += a * *A2; }
//                 e *= 1.0 - a*a;
//             }
//             a = *X;
//             for (size_t q=0u; q<P; ++q, ++A1) { --X; a += *X * *A1; }
//             A1 -= P;
//             *Y++ = a/-e;
//         }
//         else
//         {
//             const size_t F = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=Lx, ++Y)
//                 {
//                     a = -*(X+1) / *X;
//                     *A1 = a; *Y++ = a;
//                     e = *X++; e += a * *X++;
//                     for (size_t p=1u; p<P-1u; ++p, X+=p)
//                     {
//                         a = *X;
//                         for (size_t q=0u; q<p; ++q, ++A1) { --X; a += *X * *A1; }
//                         a /= -e;
//                         *A1 = a; *Y++ = a;
//                         for (size_t q=0u; q<p; ++q, ++A2) { --A1; *A2 = *A1; }
//                         A1 += p;
//                         for (size_t q=0u; q<p; ++q) { --A2; --A1; *A1 += a * *A2; }
//                         e *= 1.0 - a*a;
//                     }
//                     a = *X;
//                     for (size_t q=0u; q<P; ++q, ++A1) { --X; a += *X * *A1; }
//                     A1 -= P;
//                     *Y = a/-e;
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(P-1u))
//                 {
//                     for (size_t b=0u; b<B; ++b, X-=K*Lx-K-1u, Y-=K*P-K-1u)
//                     {
//                         a = -*(X+1) / *X;
//                         *A1 = a; *Y = a; Y += K;
//                         e = *X; X += K;
//                         e += a * *X; X += K;
//                         for (size_t p=1u; p<P-1u; ++p, X+=p*K)
//                         {
//                             a = *X;
//                             for (size_t q=0u; q<p; ++q, ++A1) { X-=K; a += *X * *A1; }
//                             a /= -e;
//                             *A1 = a; *Y = a; Y += K;
//                             for (size_t q=0u; q<p; ++q, ++A2) { --A1; *A2 = *A1; }
//                             A1 += p;
//                             for (size_t q=0u; q<p; ++q) { --A2; --A1; *A1 += a * *A2; }
//                             e *= 1.0 - a*a;
//                         }
//                         a = *X;
//                         for (size_t q=0u; q<P; ++q, ++A1) { X-=K; a += *X * *A1; }
//                         A1 -= P;
//                         *Y = a/-e;
//                     }
//                 }
//             }
//         }
//         free(A1); free(A2);
//     }
	
// 	return 0;
// }


// int ac2mvdr_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const float preg)
// {
//     if (dim>3u) { fprintf(stderr,"error in ac2mvdr_c: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

//     if (N==0u || Lx<2u) {}
//     else if (Lx==2u)
//     {
//         if (Lx==N)
//         {
//             const float den = *X**X + *(X+1)**(X+1);
//             *Y++ = -(*X**(X+2) + *(X+1)**(X+3)) / den;
//             *Y++ = (*(X+1)**(X+2) - *X**(X+3)) / den;
//         }
//         else
//         {
//             const size_t F = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;
//             float den;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=4)
//                 {
//                     den = *X**X + *(X+1)**(X+1);
//                     *Y++ = -(*X**(X+2) + *(X+1)**(X+3)) / den;
//                     *Y++ = (*(X+1)**(X+2) - *X**(X+3)) / den;
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=2u*B)
//                 {
//                     for (size_t b=0u; b<B; ++b, X+=2u, Y+=2u)
//                     {
//                         den = *X**X + *(X+1)**(X+1);
//                         *Y = -(*X**(X+2u*K) + *(X+1)**(X+2u*K+1u)) / den;
//                         *(Y+1) = (*(X+1)**(X+2u*K) - *X**(X+2u*K+1u)) / den;
//                     }
//                 }
//             }
//         }
//     }
//     else
//     {
//         const size_t P = Lx - 1u;
//         float *A1, *A2, ar, ai, e, den;
//         if (!(A1=(float *)malloc(2u*P*sizeof(float)))) { fprintf(stderr,"error in ac2mvdr_c: problem with malloc. "); perror("malloc"); return 1; }
//         if (!(A2=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2mvdr_c: problem with malloc. "); perror("malloc"); return 1; }
        
//         if (Lx==N)
//         {
//             den = *X**X + *(X+1)**(X+1);
//             ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
//             ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
//             *A1 = ar; *(A1+1) = ai;
//             *Y = ar; *(Y+1) = ai;
//             e = *X * (1.0f - (ar*ar+ai*ai));
//             X += 4; Y += 2;
//             for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
//             {
//                 ar = *X; ai = *(X+1);
//                 for (size_t q=0u; q<p; ++q, A1+=2)
//                 {
//                     X -= 2;
//                     ar += *X**A1 - *(X+1)**(A1+1);
//                     ai += *(X+1)**A1 + *X**(A1+1);
//                 }
//                 ar /= -e; ai /= -e;
//                 *A1 = ar; *(A1+1) = ai;
//                 *Y = ar; *(Y+1) = ai;
//                 for (size_t q=0u; q<p; ++q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
//                 A1 += 2u*p;
//                 for (size_t q=0u; q<p; ++q)
//                 {
//                     A2 -= 2; A1 -= 2;
//                     *A1 += ar**A2 - ai**(A2+1);
//                     *(A1+1) += ar**(A2+1) + ai**A2;
//                 }
//                 e *= 1.0f - (ar*ar + ai*ai);
//             }
//             ar = *X; ai = *(X+1);
//             for (size_t q=0u; q<P; ++q, A1+=2)
//             {
//                 X -= 2;
//                 ar += *X**A1 - *(X+1)**(A1+1);
//                 ai += *(X+1)**A1 + *X**(A1+1);
//             }
//             A1 -= 2u*P;
//             *Y++ = ar/-e; *Y = ai/-e;
//         }
//         else
//         {
//             const size_t F = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=2u*Lx)
//                 {
//                     den = *X**X + *(X+1)**(X+1);
//                     ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
//                     ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
//                     *A1 = ar; *(A1+1) = ai;
//                     *Y = ar; *(Y+1) = ai;
//                     e = *X * (1.0f - (ar*ar+ai*ai));
//                     X += 4; Y += 2;
//                     for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
//                     {
//                         ar = *X; ai = *(X+1);
//                         for (size_t q=0u; q<p; ++q, A1+=2)
//                         {
//                             X -= 2;
//                             ar += *X**A1 - *(X+1)**(A1+1);
//                             ai += *(X+1)**A1 + *X**(A1+1);
//                         }
//                         ar /= -e; ai /= -e;
//                         *A1 = ar; *(A1+1) = ai;
//                         *Y = ar; *(Y+1) = ai;
//                         for (size_t q=0u; q<p; ++q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
//                         A1 += 2u*p;
//                         for (size_t q=0u; q<p; ++q)
//                         {
//                             A2 -= 2; A1 -= 2;
//                             *A1 += ar**A2 - ai**(A2+1);
//                             *(A1+1) += ar**(A2+1) + ai**A2;
//                         }
//                         e *= 1.0f - (ar*ar + ai*ai);
//                     }
//                     ar = *X; ai = *(X+1);
//                     for (size_t q=0u; q<P; ++q, A1+=2)
//                     {
//                         X -= 2;
//                         ar += *X**A1 - *(X+1)**(A1+1);
//                         ai += *(X+1)**A1 + *X**(A1+1);
//                     }
//                     A1 -= 2u*P;
//                     *Y++ = ar/-e; *Y++ = ai/-e;
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
//                 {
//                     for (size_t b=0u; b<B; ++b, X+=2, Y-=2u*(K*P-K-1u))
//                     {
//                         den = *X**X + *(X+1)**(X+1);
//                         ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
//                         ai = (*(X+1)**(X+2u*K) - *(X+2u*K+1u)**X) / den;
//                         *A1 = ar; *(A1+1) = ai;
//                         *Y = ar; *(Y+1) = ai;
//                         e = *X * (1.0f - (ar*ar+ai*ai));
//                         X += 4u*K; Y += 2u*K;
//                         for (size_t p=1u; p<P-1u; ++p, X+=2u*p*K, Y+=2u*K)
//                         {
//                             ar = *X; ai = *(X+1);
//                             for (size_t q=0u; q<p; ++q, A1+=2)
//                             {
//                                 X -= 2u*K;
//                                 ar += *X**A1 - *(X+1)**(A1+1);
//                                 ai += *(X+1)**A1 + *X**(A1+1);
//                             }
//                             ar /= -e; ai /= -e;
//                             *A1 = ar; *(A1+1) = ai;
//                             *Y = ar; *(Y+1) = ai;
//                             for (size_t q=0u; q<p; ++q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
//                             A1 += 2u*p;
//                             for (size_t q=0u; q<p; ++q)
//                             {
//                                 A2 -= 2; A1 -= 2;
//                                 *A1 += ar**A2 - ai**(A2+1);
//                                 *(A1+1) += ar**(A2+1) + ai**A2;
//                             }
//                             e *= 1.0f - (ar*ar + ai*ai);
//                         }
//                         ar = *X; ai = *(X+1);
//                         for (size_t q=0u; q<P; ++q, A1+=2)
//                         {
//                             X -= 2u*K;
//                             ar += *X**A1 - *(X+1)**(A1+1);
//                             ai += *(X+1)**A1 + *X**(A1+1);
//                         }
//                         A1 -= 2u*P;
//                         *Y = ar/-e; *(Y+1) = ai/-e;
//                     }
//                 }
//             }
//         }
//         free(A1); free(A2);
//     }
	
// 	return 0;
// }


// int ac2mvdr_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t F, const double preg)
// {
//     if (dim>3u) { fprintf(stderr,"error in ac2mvdr_z: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

//     if (N==0u || Lx<2u) {}
//     else if (Lx==2u)
//     {
//         if (Lx==N)
//         {
//             const double den = *X**X + *(X+1)**(X+1);
//             *Y++ = -(*X**(X+2) + *(X+1)**(X+3)) / den;
//             *Y++ = (*(X+1)**(X+2) - *X**(X+3)) / den;
//         }
//         else
//         {
//             const size_t F = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;
//             double den;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=4)
//                 {
//                     den = *X**X + *(X+1)**(X+1);
//                     *Y++ = -(*X**(X+2) + *(X+1)**(X+3)) / den;
//                     *Y++ = (*(X+1)**(X+2) - *X**(X+3)) / den;
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=2u*B)
//                 {
//                     for (size_t b=0u; b<B; ++b, X+=2u, Y+=2u)
//                     {
//                         den = *X**X + *(X+1)**(X+1);
//                         *Y = -(*X**(X+2u*K) + *(X+1)**(X+2u*K+1u)) / den;
//                         *(Y+1) = (*(X+1)**(X+2u*K) - *X**(X+2u*K+1u)) / den;
//                     }
//                 }
//             }
//         }
//     }
//     else
//     {
//         const size_t P = Lx - 1u;
//         double *A1, *A2, ar, ai, e, den;
//         if (!(A1=(double *)malloc(2u*P*sizeof(double)))) { fprintf(stderr,"error in ac2mvdr_z: problem with malloc. "); perror("malloc"); return 1; }
//         if (!(A2=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2mvdr_z: problem with malloc. "); perror("malloc"); return 1; }
        
//         if (Lx==N)
//         {
//             den = *X**X + *(X+1)**(X+1);
//             ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
//             ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
//             *A1 = ar; *(A1+1) = ai;
//             *Y = ar; *(Y+1) = ai;
//             e = *X * (1.0 - (ar*ar+ai*ai));
//             X += 4; Y += 2;
//             for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
//             {
//                 ar = *X; ai = *(X+1);
//                 for (size_t q=0u; q<p; ++q, A1+=2)
//                 {
//                     X -= 2;
//                     ar += *X**A1 - *(X+1)**(A1+1);
//                     ai += *(X+1)**A1 + *X**(A1+1);
//                 }
//                 ar /= -e; ai /= -e;
//                 *A1 = ar; *(A1+1) = ai;
//                 *Y = ar; *(Y+1) = ai;
//                 for (size_t q=0u; q<p; ++q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
//                 A1 += 2u*p;
//                 for (size_t q=0u; q<p; ++q)
//                 {
//                     A2 -= 2; A1 -= 2;
//                     *A1 += ar**A2 - ai**(A2+1);
//                     *(A1+1) += ar**(A2+1) + ai**A2;
//                 }
//                 e *= 1.0 - (ar*ar + ai*ai);
//             }
//             ar = *X; ai = *(X+1);
//             for (size_t q=0u; q<P; ++q, A1+=2)
//             {
//                 X -= 2;
//                 ar += *X**A1 - *(X+1)**(A1+1);
//                 ai += *(X+1)**A1 + *X**(A1+1);
//             }
//             A1 -= 2u*P;
//             *Y++ = ar/-e; *Y = ai/-e;
//         }
//         else
//         {
//             const size_t F = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=2u*Lx)
//                 {
//                     den = *X**X + *(X+1)**(X+1);
//                     ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
//                     ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
//                     *A1 = ar; *(A1+1) = ai;
//                     *Y = ar; *(Y+1) = ai;
//                     e = *X * (1.0 - (ar*ar+ai*ai));
//                     X += 4; Y += 2;
//                     for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
//                     {
//                         ar = *X; ai = *(X+1);
//                         for (size_t q=0u; q<p; ++q, A1+=2)
//                         {
//                             X -= 2;
//                             ar += *X**A1 - *(X+1)**(A1+1);
//                             ai += *(X+1)**A1 + *X**(A1+1);
//                         }
//                         ar /= -e; ai /= -e;
//                         *A1 = ar; *(A1+1) = ai;
//                         *Y = ar; *(Y+1) = ai;
//                         for (size_t q=0u; q<p; ++q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
//                         A1 += 2u*p;
//                         for (size_t q=0u; q<p; ++q)
//                         {
//                             A2 -= 2; A1 -= 2;
//                             *A1 += ar**A2 - ai**(A2+1);
//                             *(A1+1) += ar**(A2+1) + ai**A2;
//                         }
//                         e *= 1.0 - (ar*ar + ai*ai);
//                     }
//                     ar = *X; ai = *(X+1);
//                     for (size_t q=0u; q<P; ++q, A1+=2)
//                     {
//                         X -= 2;
//                         ar += *X**A1 - *(X+1)**(A1+1);
//                         ai += *(X+1)**A1 + *X**(A1+1);
//                     }
//                     A1 -= 2u*P;
//                     *Y++ = ar/-e; *Y++ = ai/-e;
//                 }
//             }
//             else
//             {
//                 for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
//                 {
//                     for (size_t b=0u; b<B; ++b, X+=2, Y-=2u*(K*P-K-1u))
//                     {
//                         den = *X**X + *(X+1)**(X+1);
//                         ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
//                         ai = (*(X+1)**(X+2u*K) - *(X+2u*K+1u)**X) / den;
//                         *A1 = ar; *(A1+1) = ai;
//                         *Y = ar; *(Y+1) = ai;
//                         e = *X * (1.0 - (ar*ar+ai*ai));
//                         X += 4u*K; Y += 2u*K;
//                         for (size_t p=1u; p<P-1u; ++p, X+=2u*p*K, Y+=2u*K)
//                         {
//                             ar = *X; ai = *(X+1);
//                             for (size_t q=0u; q<p; ++q, A1+=2)
//                             {
//                                 X -= 2u*K;
//                                 ar += *X**A1 - *(X+1)**(A1+1);
//                                 ai += *(X+1)**A1 + *X**(A1+1);
//                             }
//                             ar /= -e; ai /= -e;
//                             *A1 = ar; *(A1+1) = ai;
//                             *Y = ar; *(Y+1) = ai;
//                             for (size_t q=0u; q<p; ++q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
//                             A1 += 2u*p;
//                             for (size_t q=0u; q<p; ++q)
//                             {
//                                 A2 -= 2; A1 -= 2;
//                                 *A1 += ar**A2 - ai**(A2+1);
//                                 *(A1+1) += ar**(A2+1) + ai**A2;
//                             }
//                             e *= 1.0 - (ar*ar + ai*ai);
//                         }
//                         ar = *X; ai = *(X+1);
//                         for (size_t q=0u; q<P; ++q, A1+=2)
//                         {
//                             X -= 2u*K;
//                             ar += *X**A1 - *(X+1)**(A1+1);
//                             ai += *(X+1)**A1 + *X**(A1+1);
//                         }
//                         A1 -= 2u*P;
//                         *Y = ar/-e; *(Y+1) = ai/-e;
//                     }
//                 }
//             }
//         }
//         free(A1); free(A2);
//     }
	
// 	return 0;
// }


#ifdef __cplusplus
}
}
#endif
