//This computes the IDFT (inverse discrete Fourier transformation) along dim of matrix X.
//This uses a direct matrix multiplication by the IDFT matrix.

//The real-valued case is complicated (i.e., where only real part is output).
//This will only exactly invert with DFT if the DFT output all F=ndft/2+1 positive freqs.
//That is, use of this function to invert a specific DFT requires knowledge of how that DFT was run.

//Profile notes: this is faster than fft.rad2 for power-of-2 IFFT up to ndft=32,
//and faster than fft.fftw for power-of-2 IFFT up to ndft=128,
//but quickly becomes MUCH slower for larger ndft.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int idft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int idft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);


int idft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_s: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Scaling
    const float s = sc ? 2.0f*sqrtf(0.5f*(float)ndft)/(float)ndft : 1.0f/(float)ndft;

    if (N==0u || ndft==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Init IDFT matrix multiply
        const size_t NN = ndft * ndft;
        const float P2_N = (float)(2.0*M_PI/(double)ndft);
        float *IDFT, smr;
        IDFT = (float *)aligned_alloc(sizeof(float),2u*NN*sizeof(float));
        if (!IDFT) { fprintf(stderr,"error in idft_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t l=0u; l<ndft; ++l) { *IDFT++ = s; *IDFT++ = 0.0f; }
        for (size_t n=1u; n<ndft; ++n)
        {
            for (size_t l=0u; l<ndft; ++l)
            {
                *IDFT++ = s * cosf(P2_N*(float)n*(float)l);
                *IDFT++ = s * sinf(P2_N*(float)n*(float)l);
            }
        }
        IDFT -= 2u*NN;

        if (Lx==N)
        {
            //Matrix multiply
            for (size_t nf=ndft; nf>0u; --nf, X-=2u)
            {
                smr = 0.0f;
                for (size_t l=Lx; l>0u; --l) //positive freqs
                {
                    smr += *X++ * *IDFT++;
                    smr -= *X++ * *IDFT++;
                }
                if (ndft%2u==0u) { X -= 2; }
                for (size_t nf=ndft-ndft/2u; nf>1u; --nf) //negative freqs
                {
                    X -= 2;
                    smr += *X * *IDFT++;
                    smr += *(X+1) * *IDFT++;
                }
                *Y++ = smr;
            }
            IDFT -= 2u*NN;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=ndft)
                {
                    //Matrix multiply
                    for (size_t nf=ndft; nf>0u; --nf, X-=Lx)
                    {
                        smr = 0.0f;
                        for (size_t l=Lx; l>0u; --l) //positive freqs
                        {
                            smr += *X++ * *IDFT++;
                            smr -= *X++ * *IDFT++;
                        }
                        if (ndft%2u==0u) { X -= 2; }
                        for (size_t nf=ndft-ndft/2u; nf>1u; --nf) //negative freqs
                        {
                            X -= 2;
                            smr += *X * *IDFT++;
                            smr += *(X+1) * *IDFT++;
                        }
                        *Y++ = smr;
                    }
                    IDFT -= 2u*NN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u*K*(ndft-ndft/2u-ndft%2u-Lx)+2u, ++Y)
                    {
                        //Matrix multiply
                        for (size_t nf=ndft; nf>0u; --nf, X-=Lx)
                        {
                            smr = 0.0f;
                            for (size_t l=Lx; l>0u; --l) //positive freqs
                            {
                                smr += *X++ * *IDFT++;
                                smr -= *X++ * *IDFT++;
                            }
                            if (ndft%2u==0u) { X -= 2; }
                            for (size_t nf=ndft-ndft/2u; nf>1u; --nf) //negative freqs
                            {
                                X -= 2;
                                smr += *X * *IDFT++;
                                smr += *(X+1) * *IDFT++;
                            }
                            *Y++ = smr;
                        }
                        IDFT -= 2u*NN;
                    }
                }
            }
        }
        free(IDFT);
    }

    return 0;
}


int idft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_d: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Scaling
    const double s = sc ? 2.0*sqrt(0.5*(double)ndft)/(double)ndft : 1.0/(double)ndft;

    if (N==0u || ndft==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Init IDFT matrix
        const size_t NN = ndft * ndft;
        const double P2_N = 2.0*M_PI/(double)ndft;
        double *IDFT, smr;
        IDFT = (double *)aligned_alloc(sizeof(double),2u*NN*sizeof(double));
        if (!IDFT) { fprintf(stderr,"error in idft_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t l=0u; l<ndft; ++l) { *IDFT++ = s; *IDFT++ = 0.0; }
        for (size_t n=1u; n<ndft; ++n)
        {
            for (size_t l=0u; l<ndft; ++l)
            {
                *IDFT++ = s * cos(P2_N*(double)n*(double)l);
                *IDFT++ = s * sin(P2_N*(double)n*(double)l);
            }
        }
        IDFT -= 2u*NN;

        if (Lx==N)
        {
            //Matrix multiply
            for (size_t nf=ndft; nf>0u; --nf, X-=Lx)
            {
                smr = 0.0;
                for (size_t l=Lx; l>0u; --l) //positive freqs
                {
                    smr += *X++ * *IDFT++;
                    smr -= *X++ * *IDFT++;
                }
                if (ndft%2u==0u) { X -= 2; }
                for (size_t nf=ndft-ndft/2u; nf>1u; --nf) //negative freqs
                {
                    X -= 2;
                    smr += *X * *IDFT++;
                    smr += *(X+1) * *IDFT++;
                }
                *Y++ = smr;
            }
            IDFT -= 2u*NN;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=ndft)
                {
                    //Matrix multiply
                    for (size_t nf=ndft; nf>0u; --nf, X-=Lx)
                    {
                        smr = 0.0;
                        for (size_t l=Lx; l>0u; --l) //positive freqs
                        {
                            smr += *X++ * *IDFT++;
                            smr -= *X++ * *IDFT++;
                        }
                        if (ndft%2u==0u) { X -= 2; }
                        for (size_t nf=ndft-ndft/2u; nf>1u; --nf) //negative freqs
                        {
                            X -= 2;
                            smr += *X * *IDFT++;
                            smr += *(X+1) * *IDFT++;
                        }
                        *Y++ = smr;
                    }
                    IDFT -= 2u*NN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u*K*(ndft-ndft/2u-ndft%2u-Lx)+2u, ++Y)
                    {
                        //Matrix multiply
                        for (size_t nf=ndft; nf>0u; --nf, X-=Lx)
                        {
                            smr = 0.0;
                            for (size_t l=Lx; l>0u; --l) //positive freqs
                            {
                                smr += *X++ * *IDFT++;
                                smr -= *X++ * *IDFT++;
                            }
                            if (ndft%2u==0u) { X -= 2; }
                            for (size_t nf=ndft-ndft/2u; nf>1u; --nf) //negative freqs
                            {
                                X -= 2;
                                smr += *X * *IDFT++;
                                smr += *(X+1) * *IDFT++;
                            }
                            *Y++ = smr;
                        }
                        IDFT -= 2u*NN;
                    }
                }
            }
        }
        free(IDFT);
    }

    return 0;
}


int idft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_c: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Scaling
    const float s = sc ? 2.0f*sqrtf(0.5f*(float)ndft)/(float)ndft : 1.0f/(float)ndft;

    if (N==0u || ndft==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Init IDFT matrix multiply
        const size_t LN = Lx * ndft;
        const float P2_N = (float)(2.0*M_PI/(double)ndft);
        float *IDFT, smr, smi, xr, xi, dr, di;
        IDFT = (float *)aligned_alloc(sizeof(float),2u*LN*sizeof(float));
        if (!IDFT) { fprintf(stderr,"error in idft_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t l=0u; l<Lx; ++l) { *IDFT++ = s; *IDFT++ = 0.0f; }
        for (size_t n=1u; n<ndft; ++n)
        {
            for (size_t l=0u; l<Lx; ++l)
            {
                *IDFT++ = s * cosf(P2_N*(float)n*(float)l);
                *IDFT++ = s * sinf(P2_N*(float)n*(float)l);
            }
        }
        IDFT -= 2u*LN;

        if (Lx==N)
        {
            //Matrix multiply
            for (size_t nf=ndft; nf>0u; --nf, X-=2u*Lx)
            {
                smr = smi = 0.0f;
                for (size_t l=Lx; l>0u; --l)
                {
                    xr = *X++; xi = *X++;
                    dr = *IDFT++; di = *IDFT++;
                    smr += xr*dr - xi*di;
                    smi += xr*di + xi*dr;
                }
                *Y++ = smr; *Y++ = smi;
            }
            IDFT -= 2u*LN;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*Lx)
                {
                    //Matrix multiply
                    for (size_t nf=ndft; nf>0u; --nf, X-=2u*Lx)
                    {
                        smr = smi = 0.0f;
                        for (size_t l=Lx; l>0u; --l)
                        {
                            xr = *X++; xi = *X++;
                            dr = *IDFT++; di = *IDFT++;
                            smr += xr*dr - xi*di;
                            smi += xr*di + xi*dr;
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    IDFT -= 2u*LN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*ndft-2u)
                    {
                        //Matrix multiply
                        for (size_t nf=ndft; nf>0u; --nf, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K)
                            {
                                xr = *X; xi = *(X+1);
                                dr = *IDFT++; di = *IDFT++;
                                smr += xr*dr - xi*di;
                                smi += xr*di + xi*dr;
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        IDFT -= 2u*LN;
                    }
                }
            }
        }
        free(IDFT);
    }

    return 0;
}


int idft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_z: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Scaling
    const double s = sc ? 2.0*sqrt(0.5*(double)ndft)/(double)ndft : 1.0/(double)ndft;

    if (N==0u || ndft==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Init IDFT matrix
        const size_t LN = Lx * ndft;
        const double P2_N = 2.0*M_PI/(double)ndft;
        double *IDFT, smr, smi, xr, xi, dr, di;
        IDFT = (double *)aligned_alloc(sizeof(double),2u*LN*sizeof(double));
        if (!IDFT) { fprintf(stderr,"error in idft_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t l=0u; l<Lx; ++l) { *IDFT++ = s; *IDFT++ = 0.0; }
        for (size_t n=1u; n<ndft; ++n)
        {
            for (size_t l=0u; l<Lx; ++l)
            {
                *IDFT++ = s * cos(P2_N*(double)n*(double)l);
                *IDFT++ = s * sin(P2_N*(double)n*(double)l);
            }
        }
        IDFT -= 2u*LN;

        if (Lx==N)
        {
            //Matrix multiply
            for (size_t nf=ndft; nf>0u; --nf, X-=2u*Lx)
            {
                smr = smi = 0.0;
                for (size_t l=Lx; l>0u; --l)
                {
                    xr = *X++; xi = *X++;
                    dr = *IDFT++; di = *IDFT++;
                    smr += xr*dr - xi*di;
                    smi += xr*di + xi*dr;
                }
                *Y++ = smr; *Y++ = smi;
            }
            IDFT -= 2u*LN;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*Lx)
                {
                    //Matrix multiply
                    for (size_t nf=ndft; nf>0u; --nf, X-=2u*Lx)
                    {
                        smr = smi = 0.0;
                        for (size_t l=Lx; l>0u; --l)
                        {
                            xr = *X++; xi = *X++;
                            dr = *IDFT++; di = *IDFT++;
                            smr += xr*dr - xi*di;
                            smi += xr*di + xi*dr;
                        }
                        *Y++ = smr; *Y++ = smi;
                    }
                    IDFT -= 2u*LN;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y-=2u*K*ndft-2u)
                    {
                        //Matrix multiply
                        for (size_t nf=ndft; nf>0u; --nf, X-=2u*K*Lx, Y+=2u*K)
                        {
                            smr = smi = 0.0;
                            for (size_t l=Lx; l>0u; --l, X+=2u*K)
                            {
                                xr = *X; xi = *(X+1);
                                dr = *IDFT++; di = *IDFT++;
                                smr += xr*dr - xi*di;
                                smi += xr*di + xi*dr;
                            }
                            *Y = smr; *(Y+1) = smi;
                        }
                        IDFT -= 2u*LN;
                    }
                }
            }
        }
        free(IDFT);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
