//Gets power spectral densities (PSDs) from autoregressive (AR) params for each vec in X.
//The 2nd input is the vector E of variances (prediction errors) for each vec in X.
//The 3rd input is a vector W of F freqs (in radians) at which to get the PSD.

//Following convention of Octave signal package ar_psd.m, I double the power for real-valued X.
//See Eq. (2.38) of Kay and Marple [1981].

//According to Octave:
//This function is intended for use with [a,v,k] = arburg(x,poles,criterion);
//which uses the Burg (1968) method to calculate a maximum-entropy AR model of X.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ar2psd_s (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2psd_d (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2psd_c (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2psd_z (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int ar2psd_s (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2psd_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2psd_s: polynomial input (X) empty\n"); return 1; }

    if (F==0u) {}
    else
    {
        const size_t P = Lx, FP = F*P;
        float *Er, *Ei, yr, yi, wp;
        if (!(Er=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ar2psd_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ar2psd_s: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        for (size_t f=0u; f<F; ++f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (float)(p+1u);
                *Er = cosf(wp);
                *Ei = -sinf(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            const float v2 = 2.0f**E;
            for (size_t f=0u; f<F; ++f, X-=P, ++Y)
            {
                yr = 1.0f; yi = 0.0f;
                for (size_t p=0u; p<P; ++p, ++X, ++Er, ++Ei)
                {
                    yr += *Er * *X;
                    yi += *Ei * *X;
                }
                *Y = v2 / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            float v2;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=P)
                {
                    v2 = 2.0f**E++;
                    for (size_t f=0u; f<F; ++f, X-=P, ++Y)
                    {
                        yr = 1.0f; yi = 0.0f;
                        for (size_t p=0u; p<P; ++p, ++X, ++Er, ++Ei)
                        {
                            yr += *X * *Er;
                            yi += *X * *Ei;
                        }
                        *Y = v2 / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=0u; b<B; ++b, ++X, Y-=K*F-1u)
                    {
                        v2 = 2.0f**E++;
                        for (size_t f=0u; f<F; ++f, X-=K*P, Y+=K)
                        {
                            yr = 1.0f; yi = 0.0f;
                            for (size_t p=0u; p<P; ++p, X+=K, ++Er, ++Ei)
                            {
                                yr += *X * *Er;
                                yi += *X * *Ei;
                            }
                            *Y = v2 / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(Er); free(Ei);
    }

    return 0;
}


int ar2psd_d (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2psd_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2psd_d: polynomial input (X) empty\n"); return 1; }

    if (F==0u) {}
    else
    {
        const size_t P = Lx, FP = F*P;
        double *Er, *Ei, yr, yi, wp;
        if (!(Er=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ar2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ar2psd_d: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        for (size_t f=0u; f<F; ++f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (double)(p+1u);
                *Er = cos(wp);
                *Ei = -sin(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            const double v2 = 2.0**E;
            for (size_t f=0u; f<F; ++f, X-=P, ++Y)
            {
                yr = 1.0; yi = 0.0;
                for (size_t p=0u; p<P; ++p, ++X, ++Er, ++Ei)
                {
                    yr += *Er * *X;
                    yi += *Ei * *X;
                }
                *Y = v2 / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            double v2;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=P)
                {
                    v2 = 2.0**E++;
                    for (size_t f=0u; f<F; ++f, X-=P, ++Y)
                    {
                        yr = 1.0; yi = 0.0;
                        for (size_t p=0u; p<P; ++p, ++X, ++Er, ++Ei)
                        {
                            yr += *X * *Er;
                            yi += *X * *Ei;
                        }
                        *Y = v2 / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=0u; b<B; ++b, ++X, Y-=K*F-1u)
                    {
                        v2 = 2.0**E++;
                        for (size_t f=0u; f<F; ++f, X-=K*P, Y+=K)
                        {
                            yr = 1.0; yi = 0.0;
                            for (size_t p=0u; p<P; ++p, X+=K, ++Er, ++Ei)
                            {
                                yr += *X * *Er;
                                yi += *X * *Ei;
                            }
                            *Y = v2 / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(Er); free(Ei);
    }

    return 0;
}


int ar2psd_c (float *Y, const float *X, const float *E, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2psd_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2psd_c: polynomial input (X) empty\n"); return 1; }

    if (F==0u) {}
    else
    {
        const size_t P = Lx, FP = F*P;
        float *Er, *Ei, yr, yi, wp;
        if (!(Er=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ar2psd_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in ar2psd_c: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        for (size_t f=0u; f<F; ++f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (float)(p+1u);
                *Er = cosf(wp);
                *Ei = -sinf(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            const float v2 = *E;
            for (size_t f=0u; f<F; ++f, X-=2u*P, ++Y)
            {
                yr = 1.0f; yi = 0.0f;
                for (size_t p=0u; p<P; ++p, X+=2, ++Er, ++Ei)
                {
                    yr += *Er**X - *Ei**(X+1);
                    yi += *Ei**X + *Er**(X+1);
                }
                *Y = v2 / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            float v2;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*P)
                {
                    v2 = *E++;
                    for (size_t f=0u; f<F; ++f, X-=2u*P, ++Y)
                    {
                        yr = 1.0f; yi = 0.0f;
                        for (size_t p=0u; p<P; ++p, X+=2, ++Er, ++Ei)
                        {
                            yr += *Er**X - *Ei**(X+1);
                            yi += *Ei**X + *Er**(X+1);
                        }
                        *Y = v2 / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=0u; b<B; ++b, X+=2, Y-=K*F-1u)
                    {
                        v2 = *E++;
                        for (size_t f=0u; f<F; ++f, X-=2u*K*P, Y+=K)
                        {
                            yr = 1.0f; yi = 0.0f;
                            for (size_t p=0u; p<P; ++p, X+=2u*K, ++Er, ++Ei)
                            {
                                yr += *Er**X - *Ei**(X+1);
                                yi += *Ei**X + *Er**(X+1);
                            }
                            *Y = v2 / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(Er); free(Ei);
    }

    return 0;
}


int ar2psd_z (double *Y, const double *X, const double *E, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2psd_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2psd_z: polynomial input (X) empty\n"); return 1; }

    if (F==0u) {}
    else
    {
        const size_t P = Lx, FP = F*P;
        double *Er, *Ei, yr, yi, wp;
        if (!(Er=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ar2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in ar2psd_z: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        for (size_t f=0u; f<F; ++f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (double)(p+1u);
                *Er = cos(wp);
                *Ei = -sin(wp);
            }
        }
        Er -= FP; Ei -= FP;

        if (Lx==N)
        {
            const double v2 = *E;
            for (size_t f=0u; f<F; ++f, X-=2u*P, ++Y)
            {
                yr = 1.0; yi = 0.0;
                for (size_t p=0u; p<P; ++p, X+=2, ++Er, ++Ei)
                {
                    yr += *Er**X - *Ei**(X+1);
                    yi += *Ei**X + *Er**(X+1);
                }
                *Y = v2 / (yr*yr+yi*yi);
            }
            Er -= FP; Ei -= FP;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;
            double v2;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*P)
                {
                    v2 = *E++;
                    for (size_t f=0u; f<F; ++f, X-=2u*P, ++Y)
                    {
                        yr = 1.0; yi = 0.0;
                        for (size_t p=0u; p<P; ++p, X+=2, ++Er, ++Ei)
                        {
                            yr += *Er**X - *Ei**(X+1);
                            yi += *Ei**X + *Er**(X+1);
                        }
                        *Y = v2 / (yr*yr+yi*yi);
                    }
                    Er -= FP; Ei -= FP;
                }
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=B*(F-1u))
                {
                    for (size_t b=0u; b<B; ++b, X+=2, Y-=K*F-1u)
                    {
                        v2 = *E++;
                        for (size_t f=0u; f<F; ++f, X-=2u*K*P, Y+=K)
                        {
                            yr = 1.0; yi = 0.0;
                            for (size_t p=0u; p<P; ++p, X+=2u*K, ++Er, ++Ei)
                            {
                                yr += *Er**X - *Ei**(X+1);
                                yi += *Ei**X + *Er**(X+1);
                            }
                            *Y = v2 / (yr*yr+yi*yi);
                        }
                        Er -= FP; Ei -= FP;
                    }
                }
            }
        }
        free(Er); free(Ei);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
