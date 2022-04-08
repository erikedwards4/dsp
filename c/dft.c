//This computes the DFT (discrete Fourier transformation) along dim of matrix X.
//This uses a direct matrix multiplication by the DFT matrix.

//This can be useful for smaller transform sizes, especially odd-length ones.
//It is nice that the scaling can be included in the matrix.

//This can also be useful if using only a few of the output frequencies, as given by F.
//Especially useful (with modification) if can be combined with another matrix multiply.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "codee_dsp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int dft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in dft_s: ndft must be >= Lx (length of vecs in X)\n"); return 1; }
    if (ndft<F) { fprintf(stderr,"error in dft_s: ndft must be >= F (num output freqs)\n"); return 1; }

    //Scaling
    const float s = sc ? 1.0f/sqrtf((float)(2u*ndft)) : 1.0f;

    if (N==0u || F==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        Y -= 2u*N;
    }
    else if (F==1u)
    {
        //DC term
        float sm = 0.0f;
        for (size_t l=Lx; l>0u; --l) { sm += *X++; }
        *Y++ = s * sm; *Y++ = 0.0f;
    }
    else
    {
        //Init DFT matrix multiply (real, imag parts separate)
        const size_t LF1 = Lx * (F-1u);
        const float P2_N = (float)(2.0*M_PI/(double)ndft);
        float *DFTr, *DFTi, smr, smi;
        DFTr = (float *)aligned_alloc(sizeof(float),LF1*sizeof(float));
        DFTi = (float *)aligned_alloc(sizeof(float),LF1*sizeof(float));
        if (!DFTr) { fprintf(stderr,"error in dft_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!DFTi) { fprintf(stderr,"error in dft_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        //for (size_t l=0u; l<Lx; ++l, ++DFTr, ++DFTi) { *DFTr = s; *DFTi = 0.0f; }
        for (size_t f=1u; f<F; ++f)
        {
            for (size_t l=0u; l<Lx; ++l)
            {
                *DFTr++ = s * cosf(P2_N*(float)f*(float)l);
                *DFTi++ = -s * sinf(P2_N*(float)f*(float)l);
            }
        }
        DFTr -= LF1; DFTi -= LF1;
    
        if (Lx==N)
        {
            //DC term
            smr = 0.0f;
            for (size_t l=Lx; l>0u; --l) { smr += *X++; }
            *Y++ = s*smr; *Y++ = 0.0f;
            X -= Lx;
            
            //Matrix multiply
            for (size_t f=F; f>1u; --f, X-=Lx)
            {
                smr = smi = 0.0f;
                for (size_t l=Lx; l>0u; --l, ++X)
                {
                    smr += *X * *DFTr++;
                    smi += *X * *DFTi++;
                }
                *Y++ = smr; *Y++ = smi;
            }
            DFTr -= LF1; DFTi -= LF1;

            //Enforce real Nyquist
            if (ndft%2u==0u && F>ndft/2u) { *(Y-1) = 0.0f; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    //DC term
                    smr = 0.0f;
                    for (size_t l=Lx; l>0u; --l) { smr += *X++; }
                    *Y++ = s*smr; *Y++ = 0.0f;
                    X -= Lx;
                    
                    //Matrix multiply
                    for (size_t f=F; f>1u; --f, X-=Lx)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=Lx; l>0u; --l, ++X)
                        {
                            smr += *X * *DFTr++;
                            smi += *X * *DFTi++;
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    DFTr -= LF1; DFTi -= LF1;

                    //Enforce real Nyquist
                    if (ndft%2u==0u && F>ndft/2u) { *(Y-1) = 0.0f; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=2u*K*F-2u)
                    {
                        //DC term
                        smr = 0.0f;
                        for (size_t l=Lx; l>0u; --l, X+=K) { smr += *X; }
                        *Y = s*smr; *(Y+1) = 0.0f;
                        X -= K*Lx; Y += 2u*K;
                        
                        //Matrix multiply
                        for (size_t f=F; f>1u; --f, X-=K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=K)
                            {
                                smr += *X * *DFTr++;
                                smi += *X * *DFTi++;
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DFTr -= LF1; DFTi -= LF1;

                        //Enforce real Nyquist
                        if (ndft%2u==0u && F>ndft/2u) { *(Y-2u*K+1u) = 0.0f; }
                    }
                }
            }
        }
        free(DFTr); free(DFTi);
    }

    return 0;
}


int dft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in dft_d: ndft must be >= Lx (length of vecs in X)\n"); return 1; }
    if (ndft<F) { fprintf(stderr,"error in dft_d: ndft must be >= F (num output freqs)\n"); return 1; }

    //Scaling
    const double s = sc ? 1.0/sqrt((double)(2u*ndft)) : 1.0;

    if (N==0u || F==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *Y++ = *X; *Y++ = 0.0; }
        Y -= 2u*N;
    }
    else if (F==1u)
    {
        //DC term
        double sm = 0.0;
        for (size_t l=Lx; l>0u; --l) { sm += *X++; }
        *Y++ = s * sm; *Y++ = 0.0;
    }
    else
    {
        //Init DFT matrix multiply (real, imag parts separate)
        const size_t LF1 = Lx * (F-1u);
        const double P2_N = 2.0*M_PI/(double)ndft;
        double *DFTr, *DFTi, smr, smi;
        DFTr = (double *)aligned_alloc(sizeof(double),LF1*sizeof(double));
        DFTi = (double *)aligned_alloc(sizeof(double),LF1*sizeof(double));
        if (!DFTr) { fprintf(stderr,"error in dft_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!DFTi) { fprintf(stderr,"error in dft_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        //for (size_t l=0u; l<Lx; ++l, ++DFTr, ++DFTi) { *DFTr = s; *DFTi = 0.0; }
        for (size_t f=1u; f<F; ++f)
        {
            for (size_t l=0u; l<Lx; ++l)
            {
                *DFTr++ = s * cos(P2_N*(double)f*(double)l);
                *DFTi++ = -s * sin(P2_N*(double)f*(double)l);
            }
        }
        DFTr -= LF1; DFTi -= LF1;
    
        if (Lx==N)
        {
            //DC term
            smr = 0.0;
            for (size_t l=Lx; l>0u; --l) { smr += *X++; }
            *Y++ = s*smr; *Y++ = 0.0;
            X -= Lx;
            
            //Matrix multiply
            for (size_t f=F; f>1u; --f, X-=Lx)
            {
                smr = smi = 0.0;
                for (size_t l=Lx; l>0u; --l, ++X)
                {
                    smr += *X * *DFTr++;
                    smi += *X * *DFTi++;
                }
                *Y++ = smr; *Y++ = smi;
            }
            DFTr -= LF1; DFTi -= LF1;

            //Enforce real Nyquist
            if (ndft%2u==0u && F>ndft/2u) { *(Y-1) = 0.0; }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    //DC term
                    smr = 0.0;
                    for (size_t l=Lx; l>0u; --l) { smr += *X++; }
                    *Y++ = s*smr; *Y++ = 0.0;
                    X -= Lx;
                    
                    //Matrix multiply
                    for (size_t f=F; f>1u; --f, X-=Lx)
                    {
                        smr = smi = 0.0;
                        for (size_t l=Lx; l>0u; --l, ++X)
                        {
                            smr += *X * *DFTr++;
                            smi += *X * *DFTi++;
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    DFTr -= LF1; DFTi -= LF1;

                    //Enforce real Nyquist
                    if (ndft%2u==0u && F>ndft/2u) { *(Y-1) = 0.0; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=2u*K*F-2u)
                    {
                        //DC term
                        smr = 0.0;
                        for (size_t l=Lx; l>0u; --l, X+=K) { smr += *X; }
                        *Y = s*smr; *(Y+1) = 0.0;
                        X -= K*Lx; Y += 2u*K;
                        
                        //Matrix multiply
                        for (size_t f=F; f>1u; --f, X-=K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=K)
                            {
                                smr += *X * *DFTr++;
                                smi += *X * *DFTi++;
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DFTr -= LF1; DFTi -= LF1;

                        //Enforce real Nyquist
                        if (ndft%2u==0u && F>ndft/2u) { *(Y-2u*K+1u) = 0.0; }
                    }
                }
            }
        }
        free(DFTr); free(DFTi);
    }

    return 0;
}


int dft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dft_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in dft_c: ndft must be >= Lx (length of vecs in X)\n"); return 1; }
    if (ndft<F) { fprintf(stderr,"error in dft_c: ndft must be >= F (num output freqs)\n"); return 1; }

    //Scaling
    const float s = sc ? 1.0f/sqrtf((float)(2u*ndft)) : 1.0f;

    if (N==0u || F==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
    }
    else
    {
        //Init DFT matrix multiply
        const size_t LF = Lx * F;
        const float P2_N = (float)(2.0*M_PI/(double)ndft);
        float *DFT, smr, smi, xr, xi, dr, di;
        DFT = (float *)aligned_alloc(sizeof(float),2u*LF*sizeof(float));
        if (!DFT) { fprintf(stderr,"error in dft_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t f=0u; f<F; ++f)
        {
            for (size_t l=0u; l<Lx; ++l)
            {
                *DFT++ = s * cosf(P2_N*(float)f*(float)l);
                *DFT++ = -s * sinf(P2_N*(float)f*(float)l);
            }
        }
        DFT -= 2u*LF;
    
        if (Lx==N)
        {
            //Matrix multiply
            for (size_t f=F; f>0u; --f, X-=2u*Lx)
            {
                smr = smi = 0.0f;
                for (size_t l=Lx; l>0u; --l)
                {
                    xr = *X++; xi = *X++;
                    dr = *DFT++; di = *DFT++;
                    smr += xr*dr - xi*di;
                    smi += xr*di + xi*dr;
                }
                *Y++ = smr; *Y++ = smi;
            }
            DFT -= 2u*LF;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    //Matrix multiply
                    for (size_t f=F; f>0u; --f, X-=2u*Lx)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=Lx; l>0u; --l)
                        {
                            xr = *X++; xi = *X++;
                            dr = *DFT++; di = *DFT++;
                            smr += xr*dr - xi*di;
                            smi += xr*di + xi*dr;
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    DFT -= 2u*LF;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*F-2u)
                    {
                        //Matrix multiply
                        for (size_t f=F; f>0u; --f, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0f;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K)
                            {
                                xr = *X; xi = *(X+1);
                                dr = *DFT++; di = *DFT++;
                                smr += xr*dr - xi*di;
                                smi += xr*di + xi*dr;
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DFT -= 2u*LF;
                    }
                }
            }
        }
        free(DFT);
    }

    return 0;
}


int dft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const size_t F, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in dft_z: ndft must be >= Lx (length of vecs in X)\n"); return 1; }
    if (ndft<F) { fprintf(stderr,"error in dft_z: ndft must be >= F (num output freqs)\n"); return 1; }

    //Scaling
    const double s = sc ? 1.0/sqrt((double)(2u*ndft)) : 1.0;

    if (N==0u || F==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
    }
    else
    {
        //Init DFT matrix multiply
        const size_t LF = Lx * F;
        const double P2_N = 2.0*M_PI/(double)ndft;
        double *DFT, smr, smi, xr, xi, dr, di;
        DFT = (double *)aligned_alloc(sizeof(double),2u*LF*sizeof(double));
        if (!DFT) { fprintf(stderr,"error in dft_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t f=0u; f<F; ++f)
        {
            for (size_t l=0u; l<Lx; ++l)
            {
                *DFT++ = s * cos(P2_N*(double)f*(double)l);
                *DFT++ = -s * sin(P2_N*(double)f*(double)l);
            }
        }
        DFT -= 2u*LF;
    
        if (Lx==N)
        {
            //Matrix multiply
            for (size_t f=F; f>0u; --f, X-=2u*Lx)
            {
                smr = smi = 0.0;
                for (size_t l=Lx; l>0u; --l)
                {
                    xr = *X++; xi = *X++;
                    dr = *DFT++; di = *DFT++;
                    smr += xr*dr - xi*di;
                    smi += xr*di + xi*dr;
                }
                *Y++ = smr; *Y++ = smi;
            }
            DFT -= LF; DFT -= LF;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v)
                {
                    //Matrix multiply
                    for (size_t f=F; f>0u; --f, X-=2u*Lx)
                    {
                        smr = smi = 0.0;
                        for (size_t l=Lx; l>0u; --l)
                        {
                            xr = *X++; xi = *X++;
                            dr = *DFT++; di = *DFT++;
                            smr += xr*dr - xi*di;
                            smi += xr*di + xi*dr;
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    DFT -= LF; DFT -= LF;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(F-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*F-2u)
                    {
                        //Matrix multiply
                        for (size_t f=F; f>0u; --f, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K)
                            {
                                xr = *X; xi = *(X+1);
                                dr = *DFT++; di = *DFT++;
                                smr += xr*dr - xi*di;
                                smi += xr*di + xi*dr;
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        DFT -= LF; DFT -= LF;
                    }
                }
            }
        }
        free(DFT);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
