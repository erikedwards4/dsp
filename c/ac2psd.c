//Gets power spectral densities (PSDs) from autocovariance (AC) for each vec in X.
//The 2nd input is a vector W of F freqs (in radians) at which to get the PSD.
//Note that Octave uses freqs in Hz (with default of Fs=1 Hz).

//Following convention of Octave signal package ar_psd.m, I double the power for real-valued X.
//See Eq. (2.38) of Kay and Marple [1981].
//I have confirmed that this matches Octave output for real and complex.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int ac2psd_s (float *Y, const float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2psd_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ac2psd_s: autocovariance input (X) empty\n"); return 1; }

    if (F==0u) {}
    else if (Lx==1u)
    {
        float v2;
        for (size_t n=N; n>0u; --n, ++X)
        {
            v2 = 2.0f * *X;
            for (size_t f=F; f>0u; --f, ++Y) { *Y = v2; }
        }
    }
    else
    {
        //Initialize AC-to-AR
        const size_t P = Lx - 1u;
        float *A1, *A2, a, e;
        if (!(A1=(float *)malloc(P*sizeof(float)))) { fprintf(stderr,"error in ac2psd_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2psd_s: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        float *Er, *Ei, yr, yi, wp;
        if (!(Er=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ac2psd_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ac2psd_s: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t f=F; f>0u; --f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (float)(p+1u);
                *Er = -cosf(wp);
                *Ei = sinf(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            //AC-to-AR
            e = *X++;
            *A1 = *X / e; a = -*A1;
            e += a * *X++;
            for (size_t p=1u; p<P; ++p, X+=p)
            {
                a = *X;
                for (size_t q=p; q>0u; --q, ++A1) { --X; a -= *X * *A1; }
                a /= -e; *A1 = -a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0f - a*a;
            }
            
            //AR-to-PSD
            e *= 2.0f;
            for (size_t f=F; f>0u; --f, A1-=P, ++Y)
            {
                yr = 1.0f; yi = 0.0f;
                for (size_t p=P; p>0u; --p, ++A1, ++Er, ++Ei)
                {
                    yr += *Er * *A1;
                    yi += *Ei * *A1;
                }
                *Y = e / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
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
                    //AC-to-AR
                    e = *X++;
                    *A1 = *X / e; a = -*A1;
                    e += a * *X++;
                    for (size_t p=1u; p<P; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=p; q>0u; --q, ++A1) { --X; a -= *X * *A1; }
                        a /= -e; *A1 = -a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0f - a*a;
                    }
                    
                    //AR-to-PSD
                    e *= 2.0f;
                    for (size_t f=F; f>0u; --f, A1-=P, ++Y)
                    {
                        yr = 1.0f; yi = 0.0f;
                        for (size_t p=P; p>0u; --p, ++A1, ++Er, ++Ei)
                        {
                            yr += *A1 * *Er;
                            yi += *A1 * *Ei;
                        }
                        *Y = e / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K*F-1u)
                    {
                        //AC-to-AR
                        e = *X; X += K;
                        *A1 = *X / e; a = -*A1;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=p; q>0u; --q, ++A1) { X-=K; a -= *X * *A1; }
                            a /= -e; *A1 = -a;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0f - a*a;
                        }

                        //AR-to-PSD
                        e *= 2.0f;
                        for (size_t f=F; f>0u; --f, A1-=P, Y+=K)
                        {
                            yr = 1.0f; yi = 0.0f;
                            for (size_t p=P; p>0u; --p, ++A1, ++Er, ++Ei)
                            {
                                yr += *A1 * *Er;
                                yi += *A1 * *Ei;
                            }
                            *Y = e / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(A1); free(A2); free(Er); free(Ei);
    }

    return 0;
}


int ac2psd_d (double *Y, const double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2psd_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ac2psd_d: autocovariance input (X) empty\n"); return 1; }

    if (F==0u) {}
    else if (Lx==1u)
    {
        double v2;
        for (size_t n=N; n>0u; --n, ++X)
        {
            v2 = 2.0 * *X;
            for (size_t f=F; f>0u; --f, ++Y) { *Y = v2; }
        }
    }
    else
    {
        //Initialize AC-to-AR
        const size_t P = Lx - 1u;
        double *A1, *A2, a, e;
        if (!(A1=(double *)malloc(P*sizeof(double)))) { fprintf(stderr,"error in ac2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2psd_d: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        double *Er, *Ei, yr, yi, wp;
        if (!(Er=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ac2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ac2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t f=F; f>0u; --f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (double)(p+1u);
                *Er = -cos(wp);
                *Ei = sin(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            //AC-to-AR
            e = *X++;
            *A1 = *X / e; a = -*A1;
            e += a * *X++;
            for (size_t p=1u; p<P; ++p, X+=p)
            {
                a = *X;
                for (size_t q=p; q>0u; --q, ++A1) { --X; a -= *X * *A1; }
                a /= -e; *A1 = -a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0 - a*a;
            }
            
            //AR-to-PSD
            e *= 2.0;
            for (size_t f=F; f>0u; --f, A1-=P, ++Y)
            {
                yr = 1.0; yi = 0.0;
                for (size_t p=P; p>0u; --p, ++A1, ++Er, ++Ei)
                {
                    yr += *Er * *A1;
                    yi += *Ei * *A1;
                }
                *Y = e / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
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
                    //AC-to-AR
                    e = *X++;
                    *A1 = *X / e; a = -*A1;
                    e += a * *X++;
                    for (size_t p=1u; p<P; ++p, X+=p)
                    {
                        a = *X;
                        for (size_t q=p; q>0u; --q, ++A1) { --X; a -= *X * *A1; }
                        a /= -e; *A1 = -a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0 - a*a;
                    }
                    
                    //AR-to-PSD
                    e *= 2.0;
                    for (size_t f=F; f>0u; --f, A1-=P, ++Y)
                    {
                        yr = 1.0; yi = 0.0;
                        for (size_t p=P; p>0u; --p, ++A1, ++Er, ++Ei)
                        {
                            yr += *A1 * *Er;
                            yi += *A1 * *Ei;
                        }
                        *Y = e / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K*F-1u)
                    {
                        //AC-to-AR
                        e = *X; X += K;
                        *A1 = *X / e; a = -*A1;
                        e += a * *X; X += K;
                        for (size_t p=1u; p<P; ++p, X+=p*K)
                        {
                            a = *X;
                            for (size_t q=p; q>0u; --q, ++A1) { X-=K; a -= *X * *A1; }
                            a /= -e; *A1 = -a;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0 - a*a;
                        }

                        //AR-to-PSD
                        e *= 2.0;
                        for (size_t f=F; f>0u; --f, A1-=P, Y+=K)
                        {
                            yr = 1.0; yi = 0.0;
                            for (size_t p=P; p>0u; --p, ++A1, ++Er, ++Ei)
                            {
                                yr += *A1 * *Er;
                                yi += *A1 * *Ei;
                            }
                            *Y = e / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(A1); free(A2); free(Er); free(Ei);
    }

    return 0;
}


int ac2psd_c (float *Y, const float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2psd_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ac2psd_c: autocovariance input (X) empty\n"); return 1; }

    if (F==0u) {}
    else if (Lx==1u)
    {
        float v2;
        for (size_t n=N; n>0u; --n, X+=2)
        {
            v2 = 2.0f * *X;
            for (size_t f=F; f>0u; --f, ++Y) { *Y = v2; }
        }
    }
    else
    {
        //Initialize AC-to-AR
        const size_t P = Lx - 1u;
        float *A1, *A2, ar, ai, e, den;
        if (!(A1=(float *)malloc(2u*P*sizeof(float)))) { fprintf(stderr,"error in ac2psd_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in ac2psd_c: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        float *Er, *Ei, yr, yi, wp;
        if (!(Er=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ac2psd_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ac2psd_c: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t f=F; f>0u; --f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (float)(p+1u);
                *Er = -cosf(wp);
                *Ei = sinf(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            //AC-to-AR
            den = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
            ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
            *A1 = -ar; *(A1+1) = -ai;
            e = *X * (1.0f - (ar*ar+ai*ai));
            X += 4;
            for (size_t p=1u; p<P; ++p, X+=2u*p)
            {
                ar = *X; ai = *(X+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    X -= 2;
                    ar -= *X**A1 - *(X+1)**(A1+1);
                    ai -= *(X+1)**A1 + *X**(A1+1);
                }
                ar /= -e; ai /= -e;
                *A1 = -ar; *(A1+1) = -ai;
                for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                A1 += 2u*p;
                for (size_t q=p; q>0u; --q)
                {
                    A1 -= 2; A2 -= 2;
                    *A1 += ar**A2 - ai**(A2+1);
                    *(A1+1) += ar**(A2+1) + ai**A2;
                }
                e *= 1.0f - (ar*ar + ai*ai);
            }

            //AR-to-PSD
            for (size_t f=F; f>0u; --f, A1-=2u*P, ++Y)
            {
                yr = 1.0f; yi = 0.0f;
                for (size_t p=P; p>0u; --p, A1+=2, ++Er, ++Ei)
                {
                    yr += *Er**A1 - *Ei**(A1+1);
                    yi += *Ei**A1 + *Er**(A1+1);
                }
                *Y = e / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
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
                    //AC-to-AR
                    den = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
                    ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
                    *A1 = -ar; *(A1+1) = -ai;
                    e = *X * (1.0f - (ar*ar+ai*ai));
                    X += 4;
                    for (size_t p=1u; p<P; ++p, X+=2u*p)
                    {
                        ar = *X; ai = *(X+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            X -= 2;
                            ar -= *X**A1 - *(X+1)**(A1+1);
                            ai -= *(X+1)**A1 + *X**(A1+1);
                        }
                        ar /= -e; ai /= -e;
                        *A1 = -ar; *(A1+1) = -ai;
                        for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                        A1 += 2u*p;
                        for (size_t q=p; q>0u; --q)
                        {
                            A1 -= 2; A2 -= 2;
                            *A1 += ar**A2 - ai**(A2+1);
                            *(A1+1) += ar**(A2+1) + ai**A2;
                        }
                        e *= 1.0f - (ar*ar + ai*ai);
                    }

                    //AR-to-PSD
                    for (size_t f=F; f>0u; --f, A1-=2u*P, ++Y)
                    {
                        yr = 1.0f; yi = 0.0f;
                        for (size_t p=P; p>0u; --p, A1+=2, ++Er, ++Ei)
                        {
                            yr += *Er**A1 - *Ei**(A1+1);
                            yi += *Ei**A1 + *Er**(A1+1);
                        }
                        *Y = e / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=K*F-1u)
                    {
                        //AC-to-AR
                        den = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
                        ai = -(*(X+2u*K+1u)**X - *(X+1)**(X+2u*K)) / den;
                        *A1 = -ar; *(A1+1) = -ai;
                        e = *X * (1.0f - (ar*ar+ai*ai));
                        X += 4u*K;
                        for (size_t p=1u; p<P; ++p, X+=2u*p*K)
                        {
                            ar = *X; ai = *(X+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                X -= 2u*K;
                                ar -= *X**A1 - *(X+1)**(A1+1);
                                ai -= *(X+1)**A1 + *X**(A1+1);
                            }
                            ar /= -e; ai /= -e;
                            *A1 = -ar; *(A1+1) = -ai;
                            for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                            A1 += 2u*p;
                            for (size_t q=p; q>0u; --q)
                            {
                                A1 -= 2; A2 -= 2;
                                *A1 += ar**A2 - ai**(A2+1);
                                *(A1+1) += ar**(A2+1) + ai**A2;
                            }
                            e *= 1.0f - (ar*ar + ai*ai);
                        }
                        
                        //AR-to-PSD
                        for (size_t f=F; f>0u; --f, A1-=2u*P, Y+=K)
                        {
                            yr = 1.0f; yi = 0.0f;
                            for (size_t p=P; p>0u; --p, A1+=2, ++Er, ++Ei)
                            {
                                yr += *Er**A1 - *Ei**(A1+1);
                                yi += *Ei**A1 + *Er**(A1+1);
                            }
                            *Y = e / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(A1); free(A2); free(Er); free(Ei);
    }

    return 0;
}


int ac2psd_z (double *Y, const double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ac2psd_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ac2psd_z: autocovariance input (X) empty\n"); return 1; }

    if (F==0u) {}
    else if (Lx==1u)
    {
        double v2;
        for (size_t n=N; n>0u; --n, X+=2)
        {
            v2 = 2.0 * *X;
            for (size_t f=F; f>0u; --f, ++Y) { *Y = v2; }
        }
    }
    else
    {
        //Initialize AC-to-AR
        const size_t P = Lx - 1u;
        double *A1, *A2, ar, ai, e, den;
        if (!(A1=(double *)malloc(2u*P*sizeof(double)))) { fprintf(stderr,"error in ac2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in ac2psd_z: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        double *Er, *Ei, yr, yi, wp;
        if (!(Er=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ac2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ac2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t f=F; f>0u; --f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (double)(p+1u);
                *Er = -cos(wp);
                *Ei = sin(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            //AC-to-AR
            den = *X**X + *(X+1)**(X+1);
            ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
            ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
            *A1 = -ar; *(A1+1) = -ai;
            e = *X * (1.0 - (ar*ar+ai*ai));
            X += 4;
            for (size_t p=1u; p<P; ++p, X+=2u*p)
            {
                ar = *X; ai = *(X+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    X -= 2;
                    ar -= *X**A1 - *(X+1)**(A1+1);
                    ai -= *(X+1)**A1 + *X**(A1+1);
                }
                ar /= -e; ai /= -e;
                *A1 = -ar; *(A1+1) = -ai;
                for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                A1 += 2u*p;
                for (size_t q=p; q>0u; --q)
                {
                    A1 -= 2; A2 -= 2;
                    *A1 += ar**A2 - ai**(A2+1);
                    *(A1+1) += ar**(A2+1) + ai**A2;
                }
                e *= 1.0 - (ar*ar + ai*ai);
            }

            //AR-to-PSD
            for (size_t f=F; f>0u; --f, A1-=2u*P, ++Y)
            {
                yr = 1.0; yi = 0.0;
                for (size_t p=P; p>0u; --p, A1+=2, ++Er, ++Ei)
                {
                    yr += *Er**A1 - *Ei**(A1+1);
                    yi += *Ei**A1 + *Er**(A1+1);
                }
                *Y = e / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
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
                    //AC-to-AR
                    den = *X**X + *(X+1)**(X+1);
                    ar = -(*(X+2)**X + *(X+3)**(X+1)) / den;
                    ai = -(*(X+3)**X - *(X+1)**(X+2)) / den;
                    *A1 = -ar; *(A1+1) = -ai;
                    e = *X * (1.0 - (ar*ar+ai*ai));
                    X += 4;
                    for (size_t p=1u; p<P; ++p, X+=2u*p)
                    {
                        ar = *X; ai = *(X+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            X -= 2;
                            ar -= *X**A1 - *(X+1)**(A1+1);
                            ai -= *(X+1)**A1 + *X**(A1+1);
                        }
                        ar /= -e; ai /= -e;
                        *A1 = -ar; *(A1+1) = -ai;
                        for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                        A1 += 2u*p;
                        for (size_t q=p; q>0u; --q)
                        {
                            A1 -= 2; A2 -= 2;
                            *A1 += ar**A2 - ai**(A2+1);
                            *(A1+1) += ar**(A2+1) + ai**A2;
                        }
                        e *= 1.0 - (ar*ar + ai*ai);
                    }

                    //AR-to-PSD
                    for (size_t f=F; f>0u; --f, A1-=2u*P, ++Y)
                    {
                        yr = 1.0; yi = 0.0;
                        for (size_t p=P; p>0u; --p, A1+=2, ++Er, ++Ei)
                        {
                            yr += *Er**A1 - *Ei**(A1+1);
                            yi += *Ei**A1 + *Er**(A1+1);
                        }
                        *Y = e / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=K*F-1u)
                    {
                        //AC-to-AR
                        den = *X**X + *(X+1)**(X+1);
                        ar = -(*(X+2u*K)**X + *(X+2u*K+1u)**(X+1)) / den;
                        ai = -(*(X+2u*K+1u)**X - *(X+1)**(X+2u*K)) / den;
                        *A1 = -ar; *(A1+1) = -ai;
                        e = *X * (1.0 - (ar*ar+ai*ai));
                        X += 4u*K;
                        for (size_t p=1u; p<P; ++p, X+=2u*p*K)
                        {
                            ar = *X; ai = *(X+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                X -= 2u*K;
                                ar -= *X**A1 - *(X+1)**(A1+1);
                                ai -= *(X+1)**A1 + *X**(A1+1);
                            }
                            ar /= -e; ai /= -e;
                            *A1 = -ar; *(A1+1) = -ai;
                            for (size_t q=p; q>0u; --q, A2+=2) { A1-=2; *A2 = *A1; *(A2+1) = -*(A1+1); }
                            A1 += 2u*p;
                            for (size_t q=p; q>0u; --q)
                            {
                                A1 -= 2; A2 -= 2;
                                *A1 += ar**A2 - ai**(A2+1);
                                *(A1+1) += ar**(A2+1) + ai**A2;
                            }
                            e *= 1.0 - (ar*ar + ai*ai);
                        }
                        
                        //AR-to-PSD
                        for (size_t f=F; f>0u; --f, A1-=2u*P, Y+=K)
                        {
                            yr = 1.0; yi = 0.0;
                            for (size_t p=P; p>0u; --p, A1+=2, ++Er, ++Ei)
                            {
                                yr += *Er**A1 - *Ei**(A1+1);
                                yi += *Ei**A1 + *Er**(A1+1);
                            }
                            *Y = e / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(A1); free(A2); free(Er); free(Ei);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
