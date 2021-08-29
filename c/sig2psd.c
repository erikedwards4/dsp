//Gets power spectral densities (PSDs) for each signal vec in X.
//The 2nd input is a vector W of F freqs (in radians) at which to get the PSD.
//Note that Octave uses freqs in Hz (with default of Fs=1 Hz).

//Following convention of Octave signal package ar_psd.m, I double the power for real-valued X.
//See Eq. (2.38) of Kay and Marple [1981].
//I have confirmed that this matches Octave output for real and complex.

//This first gets the AC (autocovariance) of each sig,
//and then obtains AR (autoregressive) coeffs by Levinson-Durbin recursion,
//and finally uses the complex E matrix to get the parametric PSD.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2psd_s (float *Y, float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2psd_d (double *Y, double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2psd_c (float *Y, float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);
int sig2psd_z (double *Y, double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased);


int sig2psd_s (float *Y, float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2psd_s: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in sig2psd_s: P (polynomial order) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (N==0u) { fprintf(stderr,"error in sig2psd_s: input (X) empty\n"); return 1; }
    if (L>Lx) { fprintf(stderr,"error in sig2psd_s: L (num lags in AC) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (F==0u) {}
    else
    {
        //Initialize AC and ac2ar
        float *AC, *A1, *A2, a, sm, e;
        if (!(AC=(float *)malloc(L*sizeof(float)))) { fprintf(stderr,"error in sig2psd_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(float *)malloc(P*sizeof(float)))) { fprintf(stderr,"error in sig2psd_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc((P-1u)*sizeof(float)))) { fprintf(stderr,"error in sig2psd_s: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        float *Er, *Ei, yr, yi, wp;
        if (!(Er=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in sig2psd_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in sig2psd_s: problem with malloc. "); perror("malloc"); return 1; }
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

            //AC-to-AR
            a = -*(AC+1) / *AC;
            *A1 = -a;
            e = *AC++; e += a * *AC++;
            for (size_t p=1u; p<P; ++p, AC+=p)
            {
                a = *AC;
                for (size_t q=p; q>0u; --q, ++A1) { --AC; a -= *AC * *A1; }
                a /= -e; *A1 = -a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0f - a*a;
            }
            AC -= L;
            
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
                for (size_t v=V; v>0u; --v, X+=Lx)
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

                    //AC-to-AR
                    a = -*(AC+1) / *AC;
                    *A1 = -a;
                    e = *AC++; e += a * *AC++;
                    for (size_t p=1u; p<P; ++p, AC+=p)
                    {
                        a = *AC;
                        for (size_t q=p; q>0u; --q, ++A1) { --AC; a -= *AC * *A1; }
                        a /= -e; *A1 = -a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0f - a*a;
                    }
                    AC -= L;
                    
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
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*F-1u)
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

                        //AC-to-AR
                        a = -*(AC+1) / *AC;
                        *A1 = -a;
                        e = *AC++; e += a * *AC++;
                        for (size_t p=1u; p<P; ++p, AC+=p)
                        {
                            a = *AC;
                            for (size_t q=p; q>0u; --q, ++A1) { --AC; a -= *AC * *A1; }
                            a /= -e; *A1 = -a;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0f - a*a;
                        }
                        AC -= L;

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


int sig2psd_d (double *Y, double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2psd_d: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in sig2psd_d: P (polynomial order) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (N==0u) { fprintf(stderr,"error in sig2psd_d: input (X) empty\n"); return 1; }
    if (L>Lx) { fprintf(stderr,"error in sig2psd_d: L (num lags in AC) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (F==0u) {}
    else
    {
        //Initialize AC and ac2ar
        double *AC, *A1, *A2, a, sm, e;
        if (!(AC=(double *)malloc(L*sizeof(double)))) { fprintf(stderr,"error in sig2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(double *)malloc(P*sizeof(double)))) { fprintf(stderr,"error in sig2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc((P-1u)*sizeof(double)))) { fprintf(stderr,"error in sig2psd_d: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        double *Er, *Ei, yr, yi, wp;
        if (!(Er=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in sig2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in sig2psd_d: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t f=F; f>0u; --f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (double)(p+1u);
                *Er = -cosf(wp);
                *Ei = sinf(wp);
            }
        }
        Er -= FP; Ei -= FP;

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

            //AC-to-AR
            a = -*(AC+1) / *AC;
            *A1 = -a;
            e = *AC++; e += a * *AC++;
            for (size_t p=1u; p<P; ++p, AC+=p)
            {
                a = *AC;
                for (size_t q=p; q>0u; --q, ++A1) { --AC; a -= *AC * *A1; }
                a /= -e; *A1 = -a;
                for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                A1 += p;
                for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                e *= 1.0 - a*a;
            }
            AC -= L;
            
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
                for (size_t v=V; v>0u; --v, X+=Lx)
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

                    //AC-to-AR
                    a = -*(AC+1) / *AC;
                    *A1 = -a;
                    e = *AC++; e += a * *AC++;
                    for (size_t p=1u; p<P; ++p, AC+=p)
                    {
                        a = *AC;
                        for (size_t q=p; q>0u; --q, ++A1) { --AC; a -= *AC * *A1; }
                        a /= -e; *A1 = -a;
                        for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                        A1 += p;
                        for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                        e *= 1.0 - a*a;
                    }
                    AC -= L;
                    
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
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*F-1u)
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

                        //AC-to-AR
                        a = -*(AC+1) / *AC;
                        *A1 = -a;
                        e = *AC++; e += a * *AC++;
                        for (size_t p=1u; p<P; ++p, AC+=p)
                        {
                            a = *AC;
                            for (size_t q=p; q>0u; --q, ++A1) { --AC; a -= *AC * *A1; }
                            a /= -e; *A1 = -a;
                            for (size_t q=p; q>0u; --q, ++A2) { --A1; *A2 = *A1; }
                            A1 += p;
                            for (size_t q=p; q>0u; --q) { --A2; --A1; *A1 += a * *A2; }
                            e *= 1.0 - a*a;
                        }
                        AC -= L;

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


int sig2psd_c (float *Y, float *X, const float *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2psd_c: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in sig2psd_c: P (polynomial order) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (N==0u) { fprintf(stderr,"error in sig2psd_c: input (X) empty\n"); return 1; }
    if (L>Lx) { fprintf(stderr,"error in sig2psd_c: L (num lags in AC) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (F==0u) {}
    else
    {
        //Initialize AC and ac2ar
        float *AC, *A1, *A2, ar, ai, smr, smi, den, e;
        if (!(AC=(float *)malloc(2u*L*sizeof(float)))) { fprintf(stderr,"error in sig2psd_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(float *)malloc(2u*P*sizeof(float)))) { fprintf(stderr,"error in sig2psd_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(float *)malloc(2u*(P-1u)*sizeof(float)))) { fprintf(stderr,"error in sig2psd_c: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        float *Er, *Ei, yr, yi, wp;
        if (!(Er=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in sig2psd_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(float *)malloc(FP*sizeof(float)))) { fprintf(stderr,"error in sig2psd_c: problem with malloc. "); perror("malloc"); return 1; }
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

            //AC-to-AR
            den = *AC**AC + *(AC+1)**(AC+1);
            ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
            ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
            *A1 = -ar; *(A1+1) = -ai;
            e = *AC * (1.0f - (ar*ar+ai*ai));
            AC += 4;
            for (size_t p=1u; p<P; ++p, AC+=2u*p)
            {
                ar = *AC; ai = *(AC+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    AC -= 2;
                    ar -= *AC**A1 - *(AC+1)**(A1+1);
                    ai -= *(AC+1)**A1 + *AC**(A1+1);
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
            AC -= 2u*L;

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
                for (size_t v=V; v>0u; --v, X+=2u*Lx)
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

                    //AC-to-AR
                    den = *AC**AC + *(AC+1)**(AC+1);
                    ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                    ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                    *A1 = -ar; *(A1+1) = -ai;
                    e = *AC * (1.0f - (ar*ar+ai*ai));
                    AC += 4;
                    for (size_t p=1u; p<P; ++p, AC+=2u*p)
                    {
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            AC -= 2;
                            ar -= *AC**A1 - *(AC+1)**(A1+1);
                            ai -= *(AC+1)**A1 + *AC**(A1+1);
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
                    AC -= 2u*L;

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
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=K*F-1u)
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

                        //AC-to-AR
                        den = *AC**AC + *(AC+1)**(AC+1);
                        ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                        ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                        *A1 = -ar; *(A1+1) = -ai;
                        e = *AC * (1.0f - (ar*ar+ai*ai));
                        AC += 4;
                        for (size_t p=1u; p<P; ++p, AC+=2u*p)
                        {
                            ar = *AC; ai = *(AC+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                AC -= 2;
                                ar -= *AC**A1 - *(AC+1)**(A1+1);
                                ai -= *(AC+1)**A1 + *AC**(A1+1);
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
                        AC -= 2u*L;
                        
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


int sig2psd_z (double *Y, double *X, const double *W, const size_t F, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t P, const int mnz, const int unbiased)
{
    if (dim>3u) { fprintf(stderr,"error in sig2psd_z: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in sig2psd_z: P (polynomial order) must be positive\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t L = P + 1u;
    if (N==0u) { fprintf(stderr,"error in sig2psd_z: input (X) empty\n"); return 1; }
    if (L>Lx) { fprintf(stderr,"error in sig2psd_z: L (num lags in AC) must be <= Lx (length of vecs in X)\n"); return 1; }

    if (F==0u) {}
    else
    {
        //Initialize AC and ac2ar
        double *AC, *A1, *A2, ar, ai, smr, smi, den, e;
        if (!(AC=(double *)malloc(2u*L*sizeof(double)))) { fprintf(stderr,"error in sig2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A1=(double *)malloc(2u*P*sizeof(double)))) { fprintf(stderr,"error in sig2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(A2=(double *)malloc(2u*(P-1u)*sizeof(double)))) { fprintf(stderr,"error in sig2psd_z: problem with malloc. "); perror("malloc"); return 1; }

        //Make complex-valued E matrix
        const size_t FP = F*P;
        double *Er, *Ei, yr, yi, wp;
        if (!(Er=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in sig2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Ei=(double *)malloc(FP*sizeof(double)))) { fprintf(stderr,"error in sig2psd_z: problem with malloc. "); perror("malloc"); return 1; }
        for (size_t f=F; f>0u; --f, ++W)
        {
            for (size_t p=0u; p<P; ++p, ++Er, ++Ei)
            {
                wp = *W * (double)(p+1u);
                *Er = -cosf(wp);
                *Ei = sinf(wp);
            }
        }
        Er -= FP; Ei -= FP;

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

            //AC-to-AR
            den = *AC**AC + *(AC+1)**(AC+1);
            ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
            ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
            *A1 = -ar; *(A1+1) = -ai;
            e = *AC * (1.0 - (ar*ar+ai*ai));
            AC += 4;
            for (size_t p=1u; p<P; ++p, AC+=2u*p)
            {
                ar = *AC; ai = *(AC+1);
                for (size_t q=p; q>0u; --q, A1+=2)
                {
                    AC -= 2;
                    ar -= *AC**A1 - *(AC+1)**(A1+1);
                    ai -= *(AC+1)**A1 + *AC**(A1+1);
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
            AC -= 2u*L;

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
                for (size_t v=V; v>0u; --v, X+=2u*Lx)
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

                    //AC-to-AR
                    den = *AC**AC + *(AC+1)**(AC+1);
                    ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                    ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                    *A1 = -ar; *(A1+1) = -ai;
                    e = *AC * (1.0 - (ar*ar+ai*ai));
                    AC += 4;
                    for (size_t p=1u; p<P; ++p, AC+=2u*p)
                    {
                        ar = *AC; ai = *(AC+1);
                        for (size_t q=p; q>0u; --q, A1+=2)
                        {
                            AC -= 2;
                            ar -= *AC**A1 - *(AC+1)**(A1+1);
                            ai -= *(AC+1)**A1 + *AC**(A1+1);
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
                    AC -= 2u*L;

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
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=K*F-1u)
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

                        //AC-to-AR
                        den = *AC**AC + *(AC+1)**(AC+1);
                        ar = -(*(AC+2)**AC + *(AC+3)**(AC+1)) / den;
                        ai = -(*(AC+3)**AC - *(AC+1)**(AC+2)) / den;
                        *A1 = -ar; *(A1+1) = -ai;
                        e = *AC * (1.0 - (ar*ar+ai*ai));
                        AC += 4;
                        for (size_t p=1u; p<P; ++p, AC+=2u*p)
                        {
                            ar = *AC; ai = *(AC+1);
                            for (size_t q=p; q>0u; --q, A1+=2)
                            {
                                AC -= 2;
                                ar -= *AC**A1 - *(AC+1)**(A1+1);
                                ai -= *(AC+1)**A1 + *AC**(A1+1);
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
                        AC -= 2u*L;
                        
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
