//Gets reflection coefficients (RCs) from the autocorrelation (AC) function for each vector in X.
//Uses a Levinson-Durbin recursion from the AC values.

//This adopts levinson.m from Octave's signal package as the correct answer, including the sign convention.
//This even matches Octave for some pathological cases giving Inf and/or NaN.
//The function from the tsa package may be wrong for the complex case!

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ac2rc_s (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_d (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_c (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ac2rc_z (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int ac2rc_s (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2rc_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2rc_s: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        float *A1, *A2, a, e;
        if (!(A1=(float *)malloc(P*sizeof(float)))) { fprintf(stderr,"error in ac2rc_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2rc_s: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            a = -*(X+1) / *X;
            *A1 = a; *Y++ = a;
            e = *X++; e += a * *X++;
            for (size_t p=1u; p<P-1u; ++p, X+=p)
            {
                a = *X;
                for (size_t q=p; q>0u; --q, ++A1) { --X; a += *X * *A1; }
                a /= -e;
                *A1 = a; *Y++ = a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0f - a*a;
            }
            a = *X;
            for (size_t q=P; q>0u; --q, ++A1) { --X; a += *X * *A1; }
            A1 -= P;
            a /= -e; *Y = a;
            *E = e * (1.0f-a*a);

            //This has approx. same speed, but assembly code definitely longer
            // A1[0] = 1.0f;
            // A1[1] = A2[0] = a = -X[1]/X[0];
            // Y[0] = -a;
            // e = fmaf(X[1],a,X[0]);
            // for (size_t p=2u; p<P; ++p)
            // {
            //     a = X[p];
            //     for (size_t q=1u; q<p; ++q) { a = fmaf(X[q],A1[p-q],a); }
            //     A1[p] = a = -a/e;
            //     Y[p-1u] = -a;
            //     for (size_t q=1u; q<p; ++q) { A1[q] = fmaf(a,A2[p-q-1u],A1[q]); }
            //     for (size_t q=p; q>0u; --q) { A2[q] = A1[q+1u]; }
            //     e *= fmaf(a,-a,1.0f);
            // }
            // a = X[P];
            // for (size_t q=1u; q<P; ++q) { a = fmaf(X[q],A1[P-q],a); }
            // Y[P-1u] = a/e; E[0] = e*fmaf(a,-a,1.0f);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Lx, ++Y, ++E)
                {
                    a = -*(X+1) / *X;
                    *A1 = a; *Y++ = a;
                    e = *X++; e += a * *X++;
                    for (size_t p=1u; p<P-1u; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=p; q>0u; --q, ++A1) { --X; a += *X * *A1; }
                        a /= -e;
                        *A1 = a; *Y++ = a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0f - a*a;
                    }
                    a = *X;
                    for (size_t q=P; q>0u; --q, ++A1) { --X; a += *X * *A1; }
                    A1 -= P;
                    a /= -e; *Y = a;
                    *E = e * (1.0f-a*a);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*P-K-1u, ++E)
                    {
                        a = -*(X+K) / *X;
                        *A1 = a; *Y = a; Y += K;
                        e = *X; X += K;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P-1u; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=p; q>0u; --q, ++A1) { X-=K; a += *X * *A1; }
                            a /= -e;
                            *A1 = a; *Y = a; Y += K;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0f - a*a;
                        }
                        a = *X;
                        for (size_t q=P; q>0u; --q, ++A1) { X-=K; a += *X * *A1; }
                        A1 -= P;
                        a /= -e; *Y = a;
                        *E = e * (1.0f-a*a);
                    }
                }
            }
        }
        free(A1); free(A2);
    }
	
	return 0;
}


int ac2rc_d (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
	if (dim>3u) { fprintf(stderr,"error in ac2rc_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2rc_d: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        double *A1, *A2, a, e;
        if (!(A1=(double *)malloc(P*sizeof(double)))) { fprintf(stderr,"error in ac2rc_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2rc_d: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            a = -*(X+1) / *X;
            *A1 = a; *Y++ = a;
            e = *X++; e += a * *X++;
            for (size_t p=1u; p<P-1u; ++p, X+=p)
            {
                a = *X;
                for (size_t q=p; q>0u; --q, ++A1) { --X; a += *X * *A1; }
                a /= -e;
                *A1 = a; *Y++ = a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0 - a*a;
            }
            a = *X;
            for (size_t q=P; q>0u; --q, ++A1) { --X; a += *X * *A1; }
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
                for (size_t v=V; v>0u; --v, X+=Lx, ++Y, ++E)
                {
                    a = -*(X+1) / *X;
                    *A1 = a; *Y++ = a;
                    e = *X++; e += a * *X++;
                    for (size_t p=1u; p<P-1u; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=p; q>0u; --q, ++A1) { --X; a += *X * *A1; }
                        a /= -e;
                        *A1 = a; *Y++ = a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0 - a*a;
                    }
                    a = *X;
                    for (size_t q=P; q>0u; --q, ++A1) { --X; a += *X * *A1; }
                    A1 -= P;
                    a /= -e; *Y = a;
                    *E = e * (1.0-a*a);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(P-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*P-K-1u, ++E)
                    {
                        a = -*(X+K) / *X;
                        *A1 = a; *Y = a; Y += K;
                        e = *X; X += K;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P-1u; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=p; q>0u; --q, ++A1) { X-=K; a += *X * *A1; }
                            a /= -e;
                            *A1 = a; *Y = a; Y += K;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0 - a*a;
                        }
                        a = *X;
                        for (size_t q=P; q>0u; --q, ++A1) { X-=K; a += *X * *A1; }
                        A1 -= P;
                        a /= -e; *Y = a;
                        *E = e * (1.0-a*a);
                    }
                }
            }
        }
        free(A1); free(A2);
    }
	
	return 0;
}


int ac2rc_c (float *Y, float *E, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2rc_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2rc_c: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        float *A1, *A2, ar, ai, e, den;
        if (!(A1=(float *)malloc(2u*P*sizeof(float)))) { fprintf(stderr,"error in ac2rc_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2rc_c: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            den = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
            ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
            *A1 = ar; *(A1+1) = ai;
            *Y = ar; *(Y+1) = ai;
            e = *X * (1.0f - (ar*ar+ai*ai));
            X += 4; Y += 2;
            for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
            {
                ar = *X; ai = *(X+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    X -= 2;
                    ar += *X**A1 - *(X+1)**(A1+1);
                    ai += *(X+1)**A1 + *X**(A1+1);
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
            ar = *X; ai = *(X+1);
            for (size_t q=P; q>0u; --q, A1+=2)
            {
                X -= 2;
                ar += *X**A1 - *(X+1)**(A1+1);
                ai += *(X+1)**A1 + *X**(A1+1);
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
                    den = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
                    ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
                    *A1 = ar; *(A1+1) = ai;
                    *Y = ar; *(Y+1) = ai;
                    e = *X * (1.0f - (ar*ar+ai*ai));
                    X += 4; Y += 2;
                    for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
                    {
                        ar = *X; ai = *(X+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            X -= 2;
                            ar += *X**A1 - *(X+1)**(A1+1);
                            ai += *(X+1)**A1 + *X**(A1+1);
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
                    ar = *X; ai = *(X+1);
                    for (size_t q=P; q>0u; --q, A1+=2)
                    {
                        X -= 2;
                        ar += *X**A1 - *(X+1)**(A1+1);
                        ai += *(X+1)**A1 + *X**(A1+1);
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
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*(K*P-K-1u), ++E)
                    {
                        den = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
                        ai = (*(X+1)**(X+2u*K) - *(X+2u*K+1u)**X) / den;
                        *A1 = ar; *(A1+1) = ai;
                        *Y = ar; *(Y+1) = ai;
                        e = *X * (1.0f - (ar*ar+ai*ai));
                        X += 4u*K; Y += 2u*K;
                        for (size_t p=1u; p<P-1u; ++p, X+=2u*p*K, Y+=2u*K)
                        {
                            ar = *X; ai = *(X+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                X -= 2u*K;
                                ar += *X**A1 - *(X+1)**(A1+1);
                                ai += *(X+1)**A1 + *X**(A1+1);
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
                        ar = *X; ai = *(X+1);
                        for (size_t q=P; q>0u; --q, A1+=2)
                        {
                            X -= 2u*K;
                            ar += *X**A1 - *(X+1)**(A1+1);
                            ai += *(X+1)**A1 + *X**(A1+1);
                        }
                        A1 -= 2u*P;
                        ar /= -e; ai /= -e;
                        *Y++ = ar; *(Y+1) = ai;
                        *E = e * (1.0f-ar*ar-ai*ai);
                    }
                }
            }
        }
        free(A1); free(A2);
    }
	
	return 0;
}


int ac2rc_z (double *Y, double *E, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2rc_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (Lx<2u) { fprintf(stderr,"error in ac2rc_z: ACF must have length > 1\n"); return 1; }

    if (N==0u) {}
    else
    {
        const size_t P = Lx - 1u;
        double *A1, *A2, ar, ai, e, den;
        if (!(A1=(double *)malloc(2u*P*sizeof(double)))) { fprintf(stderr,"error in ac2rc_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2rc_z: problem with malloc. "); perror("malloc"); return 1; }
        
        if (Lx==N)
        {
            den = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
            ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
            *A1 = ar; *(A1+1) = ai;
            *Y = ar; *(Y+1) = ai;
            e = *X * (1.0 - (ar*ar+ai*ai));
            X += 4; Y += 2;
            for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
            {
                ar = *X; ai = *(X+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    X -= 2;
                    ar += *X**A1 - *(X+1)**(A1+1);
                    ai += *(X+1)**A1 + *X**(A1+1);
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
            ar = *X; ai = *(X+1);
            for (size_t q=P; q>0u; --q, A1+=2)
            {
                X -= 2;
                ar += *X**A1 - *(X+1)**(A1+1);
                ai += *(X+1)**A1 + *X**(A1+1);
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
                    den = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
                    ai = (*(X+1)**(X+2) - *(X+3)**X) / den;
                    *A1 = ar; *(A1+1) = ai;
                    *Y = ar; *(Y+1) = ai;
                    e = *X * (1.0 - (ar*ar+ai*ai));
                    X += 4; Y += 2;
                    for (size_t p=1u; p<P-1u; ++p, X+=2u*p, Y+=2)
                    {
                        ar = *X; ai = *(X+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            X -= 2;
                            ar += *X**A1 - *(X+1)**(A1+1);
                            ai += *(X+1)**A1 + *X**(A1+1);
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
                    ar = *X; ai = *(X+1);
                    for (size_t q=P; q>0u; --q, A1+=2)
                    {
                        X -= 2;
                        ar += *X**A1 - *(X+1)**(A1+1);
                        ai += *(X+1)**A1 + *X**(A1+1);
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
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*(K*P-K-1u), ++E)
                    {
                        den = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
                        ai = (*(X+1)**(X+2u*K) - *(X+2u*K+1u)**X) / den;
                        *A1 = ar; *(A1+1) = ai;
                        *Y = ar; *(Y+1) = ai;
                        e = *X * (1.0 - (ar*ar+ai*ai));
                        X += 4u*K; Y += 2u*K;
                        for (size_t p=1u; p<P-1u; ++p, X+=2u*p*K, Y+=2u*K)
                        {
                            ar = *X; ai = *(X+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                X -= 2u*K;
                                ar += *X**A1 - *(X+1)**(A1+1);
                                ai += *(X+1)**A1 + *X**(A1+1);
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
                        ar = *X; ai = *(X+1);
                        for (size_t q=P; q>0u; --q, A1+=2)
                        {
                            X -= 2u*K;
                            ar += *X**A1 - *(X+1)**(A1+1);
                            ai += *(X+1)**A1 + *X**(A1+1);
                        }
                        A1 -= 2u*P;
                        ar /= -e; ai /= -e;
                        *Y++ = ar; *(Y+1) = ai;
                        *E = e * (1.0-ar*ar-ai*ai);
                    }
                }
            }
        }
        free(A1); free(A2);
    }
	
	return 0;
}


#ifdef __cplusplus
}
}
#endif
