//Gets autoregressive (AR) params from the autocorrelation (AC) function for each vector in X.
//Uses a Levinson-Durbin recursion from the AC values.

//The complex-valued case here uses multiply and divide with no conjugation.
//This matches the real-valued case (i.e., make a complex-valued input with 0 imaginary parts),
//and matches Octave's levinson.m by all current tests.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int ac2ar_s (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2ar_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2ar_s: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        float *A, a, e;
        if (!(A=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2ar_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            e = *X++;
            *Y = *X / e; a = -*Y;
            e += a * *X++;
            for (size_t p=1u; p<P; ++p, X+=p)
            {
                a = *X;
                for (size_t q=p; q>0u; --q, ++Y) { --X; a -= *X * *Y; }
                a /= -e; *Y = -a;
                for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                Y += p;
                for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                e *= 1.0f - a*a;
            }
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=P, ++E)
                {
                    e = *X++;
                    *Y = *X / e; a = -*Y;
                    e += a * *X++;
                    for (size_t p=1u; p<P; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=p; q>0u; --q, ++Y) { --X; a -= *X * *Y; }
                        a /= -e;
                        *Y = -a;
                        for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                        Y += p;
                        for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                        e *= 1.0f - a*a;
                    }
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*P, Y+=B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, ++Y, ++E)
                    {
                        e = *X; X += K;
                        *Y = *X / e; a = -*Y;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=p; q>0u; --q, Y+=K) { X-=K; a -= *X * *Y; }
                            a /= -e;
                            *Y = -a;
                            for (size_t q=p; q>0u; --q, ++A) { Y-=K; *A = *Y; }
                            Y += p*K;
                            for (size_t q=p; q>0u; --q) { --A; Y-=K; *Y += a * *A; }
                            e *= 1.0f - a*a;
                        }
                        *E = e;
                    }
                }
            }
        }
        free(A);
    }
	
	return 0;
}


int ac2ar_d (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
	if (dim>3u) { fprintf(stderr,"error in ac2ar_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2ar_d: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        double *A, a, e;
        if (!(A=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2ar_d: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            e = *X++;
            *Y = *X / e; a = -*Y;
            e += a * *X++;
            for (size_t p=1u; p<P; ++p, X+=p)
            {
                a = *X;
                for (size_t q=p; q>0u; --q, ++Y) { --X; a -= *X * *Y; }
                a /= -e;
                *Y = -a;
                for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                Y += p;
                for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                e *= 1.0 - a*a;
            }
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=P, ++E)
                {
                    e = *X++;
                    *Y = *X / e; a = -*Y;
                    e += a * *X++;
                    for (size_t p=1u; p<P; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=p; q>0u; --q, ++Y) { --X; a -= *X * *Y; }
                        a /= -e;
                        *Y = -a;
                        for (size_t q=p; q>0u; --q, ++A) { --Y; *A = *Y; }
                        Y += p;
                        for (size_t q=p; q>0u; --q) { --A; --Y; *Y += a * *A; }
                        e *= 1.0 - a*a;
                    }
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*P, Y+=B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, ++Y, ++E)
                    {
                        e = *X; X += K;
                        *Y = *X / e; a = -*Y;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=p; q>0u; --q, Y+=K) { X-=K; a -= *X * *Y; }
                            a /= -e;
                            *Y = -a;
                            for (size_t q=p; q>0u; --q, ++A) { Y-=K; *A = *Y; }
                            Y += p*K;
                            for (size_t q=p; q>0u; --q) { --A; Y-=K; *Y += a * *A; }
                            e *= 1.0 - a*a;
                        }
                        *E = e;
                    }
                }
            }
        }
        free(A);
    }
	
	return 0;
}


int ac2ar_c (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2ar_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2ar_c: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        float *A, ar, ai, e, den;
        if (!(A=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2ar_c: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            den = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
            ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
            *Y = -ar; *(Y+1) = -ai;
            e = *X * (1.0f - (ar*ar+ai*ai));
            X += 4;
            for (size_t p=1u; p<P; ++p, X+=2u*p)
            {
                ar = *X; ai = *(X+1);
                for (size_t q=p; q>0u; --q, Y+=2)
                {
                    X -= 2;
                    ar -= *X**Y - *(X+1)**(Y+1);
                    ai -= *(X+1)**Y + *X**(Y+1);
                }
                ar /= -e; ai /= -e;
                *Y = -ar; *(Y+1) = -ai;
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
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=2u*P, ++E)
                {
                    den = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
                    ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
                    *Y = -ar; *(Y+1) = -ai;
                    e = *X * (1.0f - (ar*ar+ai*ai));
                    X += 4;
                    for (size_t p=1u; p<P; ++p, X+=2u*p)
                    {
                        ar = *X; ai = *(X+1);
                        for (size_t q=p; q>0u; --q, Y+=2)
                        {
                            X -= 2;
                            ar -= *X**Y - *(X+1)**(Y+1);
                            ai -= *(X+1)**Y + *X**(Y+1);
                        }
                        ar /= -e; ai /= -e;
                        *Y = -ar; *(Y+1) = -ai;
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
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*P, Y+=2u*B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y+=2, ++E)
                    {
                        den = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
                        ai = -(*(X+2u*K+1u)**X - *(X+1)**(X+2u*K)) / den;
                        *Y = -ar; *(Y+1) = -ai;
                        e = *X * (1.0f - (ar*ar+ai*ai));
                        X += 4u*K;
                        for (size_t p=1u; p<P; ++p, X+=2u*p*K)
                        {
                            ar = *X; ai = *(X+1);
                            for (size_t q=p; q>0u; --q, Y+=2u*K)
                            {
                                X -= 2u*K;
                                ar -= *X**Y - *(X+1)**(Y+1);
                                ai -= *(X+1)**Y + *X**(Y+1);
                            }
                            ar /= -e; ai /= -e;
                            *Y = -ar; *(Y+1) = -ai;
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
                        *E = e;
                    }
                }
            }
        }
        free(A);
    }
	
	return 0;
}


int ac2ar_z (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2ar_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2ar_z: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        double *A, ar, ai, e, den;
        if (!(A=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2ar_z: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            den = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
            ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
            *Y = -ar; *(Y+1) = -ai;
            e = *X * (1.0 - (ar*ar+ai*ai));
            X += 4;
            for (size_t p=1u; p<P; ++p, X+=2u*p)
            {
                ar = *X; ai = *(X+1);
                for (size_t q=p; q>0u; --q, Y+=2)
                {
                    X -= 2;
                    ar -= *X**Y - *(X+1)**(Y+1);
                    ai -= *(X+1)**Y + *X**(Y+1);
                }
                ar /= -e; ai /= -e;
                *Y = -ar; *(Y+1) = -ai;
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
            *E = e;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=2u*P, ++E)
                {
                    den = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
                    ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
                    *Y = -ar; *(Y+1) = -ai;
                    e = *X * (1.0 - (ar*ar+ai*ai));
                    X += 4;
                    for (size_t p=1u; p<P; ++p, X+=2u*p)
                    {
                        ar = *X; ai = *(X+1);
                        for (size_t q=p; q>0u; --q, Y+=2)
                        {
                            X -= 2;
                            ar -= *X**Y - *(X+1)**(Y+1);
                            ai -= *(X+1)**Y + *X**(Y+1);
                        }
                        ar /= -e; ai /= -e;
                        *Y = -ar; *(Y+1) = -ai;
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
                    *E = e;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*P, Y+=2u*B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y+=2, ++E)
                    {
                        den = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
                        ai = -(*(X+2u*K+1u)**X - *(X+1)**(X+2u*K)) / den;
                        *Y = -ar; *(Y+1) = -ai;
                        e = *X * (1.0 - (ar*ar+ai*ai));
                        X += 4u*K;
                        for (size_t p=1u; p<P; ++p, X+=2u*p*K)
                        {
                            ar = *X; ai = *(X+1);
                            for (size_t q=p; q>0u; --q, Y+=2u*K)
                            {
                                X -= 2u*K;
                                ar -= *X**Y - *(X+1)**(Y+1);
                                ai -= *(X+1)**Y + *X**(Y+1);
                            }
                            ar /= -e; ai /= -e;
                            *Y = -ar; *(Y+1) = -ai;
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
                        *E = e;
                    }
                }
            }
        }
        free(A);
    }
	
	return 0;
}


#ifdef __cplusplus
}
}
#endif
