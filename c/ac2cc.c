//Gets cepstral coefficients (CCs) from the autocorrelation (AC) function for each vector in X.

//Starts with linear prediction (LP) by Levinson-Durbin recursion from the AC values.
//The order P of the LP is determined as Lx-1, where Lx is the length of each vector in X.

//The first CC returned is log(V+preg), where V is the prediction error variance from the LP,
//and preg is a small power regulation constant provided as input (e.g., FLT_EPSILON).

//A total of K CCs are returned, such that each vector in Y has length K.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ac2cc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t K, const double preg);
int ac2cc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t K, const double preg);
int ac2cc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t K, const double preg);
int ac2cc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t K, const double preg);


int ac2cc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t K, const double preg)
{
    if (dim>3u) { fprintf(stderr,"error in ac2cc_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2cc_s: Lx (length of vecs in X) must be > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        float *AR, *A2, a, e, sm;
        if (!(AR=(float *)malloc(P*sizeof(float)))) { fprintf(stderr,"error in ac2cc_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2cc_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            //AC-to-AR
            e = *X++;
            *Y = *X / e; a = -*Y;
            e += a * *X++;
            for (size_t p=1u; p<P; ++p, X+=p)
            {
                a = *X;
                for (size_t q=p; q>0u; --q, ++AR) { --X; a -= *X * *AR; }
                a /= -e; *AR = -a;
                for (size_t q=p; q>0u; --q, ++A2) { --AR; *A2 = *AR; }
                AR += p;
                for (size_t q=p; q>0u; --q) { --A2; --AR; *AR += a * *A2; }
                e *= 1.0f - a*a;
            }

            //AR-to-CC
            *Y++ = logf(e+preg);
            for (size_t k=1u; k<K && k<=P; ++k, ++Y, ++AR)
            {
                sm = 0.0f;
                for (size_t j=1u; j<k; ++j) { sm += (float)j * *(Y-k+j) * *(AR-j); }
                *Y = *AR - sm/(float)k;
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, AR+=P)
                {
                    //AC-to-AR
                    *AR++ = 1.0f;
                    a = -*(X+1) / *X;
                    *AR = a;
                    e = *X++; e += a * *X++;
                    for (size_t p=1u; p<P; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=p; q>0u; --q, ++AR) { --X; a += *X * *AR; }
                        a /= -e; *AR = a;
                        for (size_t q=p; q>0u; --q, ++A2) { --AR; *A2 = *AR; }
                        AR += p;
                        for (size_t q=p; q>0u; --q) { --A2; --AR; *AR += a * *A2; }
                        e *= 1.0f - a*a;
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), AR+=B*(Lx-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, AR-=K-1u)
                    {
                        //AC-to-AR
                        *AR = 1.0f; AR += K;
                        a = -*(X+K) / *X;
                        *AR = a;
                        e = *X; X += K;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=p; q>0u; --q, AR+=K) { X-=K; a += *X * *AR; }
                            a /= -e; *AR = a;
                            for (size_t q=p; q>0u; --q, ++A2) { AR-=K; *A2 = *AR; }
                            AR += p*K;
                            for (size_t q=p; q>0u; --q) { --A2; AR-=K; *AR += a * *A2; }
                            e *= 1.0f - a*a;
                        }
                    }
                }
            }
        }
        free(AR); free(A2);
    }
	
	return 0;
}


// int ac2cc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t K, const double preg)
// {
// 	if (dim>3u) { fprintf(stderr,"error in ac2cc_d: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

//     if (N==0u) {}
//     else if (Lx==1u)
//     {
//         for (size_t n=N; n>0u; --n, ++X, ++Y, ++E) { *Y = 1.0; *E = *X; }
//     }
//     else
//     {
//         const size_t P = Lx - 1u;
//         double *A, a, e;
//         if (!(A=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2cc_d: problem with malloc. "); perror("malloc"); return 1; }
        
//         if (Lx==N)
//         {
//             *Y++ = 1.0;
//             a = -*(X+1) / *X;
//             *Y = a;
//             e = *X++; e += a * *X++;
//             for (size_t p=1u; p<P; ++p, X+=p)
//             {
//                 a = *X;
//                 for (size_t q=p; q>0u; --q, ++Y) { --X; a += *X * *Y; }
//                 a /= -e;
//                 *Y = a;
//                 for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
//                 Y += p;
//                 for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
//                 e *= 1.0 - a*a;
//             }
//             *E = e;
//         }
//         else
//         {
//             const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=V; v>0u; --v, Y+=P, ++E)
//                 {
//                     *Y++ = 1.0;
//                     a = -*(X+1) / *X;
//                     *Y = a;
//                     e = *X++; e += a * *X++;
//                     for (size_t p=1u; p<P; ++p, X+=p)
//                     {
//                         a = *X;
//                         for (size_t q=p; q>0u; --q, ++Y) { --X; a += *X * *Y; }
//                         a /= -e; *Y = a;
//                         for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
//                         Y += p;
//                         for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
//                         e *= 1.0 - a*a;
//                     }
//                     *E = e;
//                 }
//             }
//             else
//             {
//                 for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Lx-1u))
//                 {
//                     for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K-1u, ++E)
//                     {
//                         *Y = 1.0; Y += K;
//                         a = -*(X+K) / *X;
//                         *Y = a;
//                         e = *X; X += K;
//                         e += a * *X; X += K;
//                         for (size_t p=1u; p<P; ++p, X+=p*K)
//                         {
//                             a = *X;
//                             for (size_t q=p; q>0u; --q, Y+=K) { X-=K; a += *X * *Y; }
//                             a /= -e; *Y = a;
//                             for (size_t q=p; q>0u; --q, ++A) { Y-=K; *A = *Y; }
//                             Y += p*K;
//                             for (size_t q=p; q>0u; --q) { --A; Y-=K; *Y += a * *A; }
//                             e *= 1.0 - a*a;
//                         }
//                         *E = e;
//                     }
//                 }
//             }
//         }
//         free(A);
//     }
	
// 	return 0;
// }


#ifdef __cplusplus
}
}
#endif
