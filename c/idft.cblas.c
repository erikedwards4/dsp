//This computes the IDFT (inverse discrete Fourier transformation) along dim of matrix X.
//This uses a CBLAS matrix multiplication by the IDFT matrix.

//The real-valued case is complicated (i.e., where only real part is output).
//This will only exactly invert with DFT if the DFT output all F=ndft/2+1 positive freqs.
//That is, use of this function to invert a specific DFT requires knowledge of how that DFT was run.

//Profile note: I tried splitting the matrix multiply into two halves,
//for positive and negative freqs, so that X would not need to be copied into Xc,
//but it was definitely slower.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "codee_dsp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int idft_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_cblas_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_cblas_s: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Scaling
    const float s = sc ? 2.0f*sqrtf(0.5f*(float)ndft)/(float)ndft : 1.0f/(float)ndft;

    if (N==0u || ndft==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Init IDFT matrix
        const size_t NN = ndft * ndft;
        const float P2_N = (float)(2.0*M_PI/(double)ndft);
        float *IDFT;
        IDFT = (float *)aligned_alloc(sizeof(float),2u*NN*sizeof(float));
        if (!IDFT) { fprintf(stderr,"error in idft_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
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

        //Init complex matrix multiply
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        float *Xc, *Yc;
        Xc = (float *)aligned_alloc(sizeof(float),2u*ndft*sizeof(float));
        Yc = (float *)aligned_alloc(sizeof(float),2u*ndft*sizeof(float));
        if (!Xc) { fprintf(stderr,"error in idft_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Yc) { fprintf(stderr,"error in idft_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<2u*ndft; ++n, ++Xc) { *Xc = 0.0f; }
        Xc -= 2u*ndft;

        if (Lx==N)
        {
            //Copy into Xc
            for (size_t nf=Lx; nf>0u; --nf) { *Xc++ = *X++; *Xc++ = *X++; }
            if (ndft%2u==0u) { X -= 2; }
            for (size_t nf=ndft-ndft/2u; nf>1u; --nf) { X-=2; *Xc++ = *X; *Xc++ = -*(X+1); }
            Xc -= 2u*ndft;

            //Matrix multiply
            cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)ndft,o,IDFT,(int)ndft,Xc,1,z,Yc,1);
            
            //Output real part only
            cblas_scopy((int)ndft,Yc,2,Y,1);
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
                    //Copy into Xc
                    for (size_t nf=Lx; nf>0u; --nf) { *Xc++ = *X++; *Xc++ = *X++; }
                    if (ndft%2u==0u) { X -= 2; }
                    for (size_t nf=ndft-ndft/2u; nf>1u; --nf) { X-=2; *Xc++ = *X; *Xc++ = -*(X+1); }
                    Xc -= 2u*ndft;
                    X += 2u*(ndft-ndft/2u-ndft%2u);

                    //Matrix multiply
                    cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)ndft,o,IDFT,(int)ndft,Xc,1,z,Yc,1);

                    //Output real part only
                    cblas_scopy((int)ndft,Yc,2,Y,1);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u*K*(ndft-ndft/2u-ndft%2u-Lx)+2u, ++Y)
                    {
                        //Copy into Xc
                        for (size_t nf=Lx; nf>0u; --nf, X+=2u*K) { *Xc++ = *X; *Xc++ = *(X+1); }
                        if (ndft%2u==0u) { X -= 2u*K; }
                        for (size_t nf=ndft-ndft/2u; nf>1u; --nf) { X-=2u*K; *Xc++ = *X; *Xc++ = -*(X+1); }
                        Xc -= 2u*ndft;

                        //Matrix multiply
                        cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)ndft,o,IDFT,(int)ndft,Xc,1,z,Yc,1);

                        //Output real part only
                        cblas_scopy((int)ndft,Yc,2,Y,(int)K);
                    }
                }
            }
        }
        free(IDFT); free(Xc); free(Yc);
    }

    return 0;
}


int idft_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_cblas_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_cblas_d: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

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
        double *IDFT;
        IDFT = (double *)aligned_alloc(sizeof(double),2u*NN*sizeof(double));
        if (!IDFT) { fprintf(stderr,"error in idft_cblas_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
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

        //Init complex matrix multiply
        const double z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
        double *Xc, *Yc;
        Xc = (double *)aligned_alloc(sizeof(double),2u*ndft*sizeof(double));
        Yc = (double *)aligned_alloc(sizeof(double),2u*ndft*sizeof(double));
        if (!Xc) { fprintf(stderr,"error in idft_cblas_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Yc) { fprintf(stderr,"error in idft_cblas_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<2u*ndft; ++n, ++Xc) { *Xc = 0.0; }
        Xc -= 2u*ndft;

        if (Lx==N)
        {
            //Copy into Xc
            for (size_t nf=Lx; nf>0u; --nf) { *Xc++ = *X++; *Xc++ = *X++; }
            if (ndft%2u==0u) { X -= 2; }
            for (size_t nf=ndft-ndft/2u; nf>1u; --nf) { X-=2; *Xc++ = *X; *Xc++ = -*(X+1); }
            Xc -= 2u*ndft;

            //Matrix multiply
            cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)ndft,o,IDFT,(int)ndft,Xc,1,z,Yc,1);
            
            //Output real part only
            cblas_dcopy((int)ndft,Yc,2,Y,1);
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
                    //Copy into Xc
                    for (size_t nf=Lx; nf>0u; --nf) { *Xc++ = *X++; *Xc++ = *X++; }
                    if (ndft%2u==0u) { X -= 2; }
                    for (size_t nf=ndft-ndft/2u; nf>1u; --nf) { X-=2; *Xc++ = *X; *Xc++ = -*(X+1); }
                    Xc -= 2u*ndft;
                    X += 2u*(ndft-ndft/2u-ndft%2u);

                    //Matrix multiply
                    cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)ndft,o,IDFT,(int)ndft,Xc,1,z,Yc,1);

                    //Output real part only
                    cblas_dcopy((int)ndft,Yc,2,Y,1);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u*K*(ndft-ndft/2u-ndft%2u-Lx)+2u, ++Y)
                    {
                        //Copy into Xc
                        for (size_t nf=Lx; nf>0u; --nf, X+=2u*K) { *Xc++ = *X; *Xc++ = *(X+1); }
                        if (ndft%2u==0u) { X -= 2u*K; }
                        for (size_t nf=ndft-ndft/2u; nf>1u; --nf) { X-=2u*K; *Xc++ = *X; *Xc++ = -*(X+1); }
                        Xc -= 2u*ndft;

                        //Matrix multiply
                        cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)ndft,o,IDFT,(int)ndft,Xc,1,z,Yc,1);

                        //Output real part only
                        cblas_dcopy((int)ndft,Yc,2,Y,(int)K);
                    }
                }
            }
        }
        free(IDFT); free(Xc); free(Yc);
    }

    return 0;
}


int idft_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_cblas_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_cblas_c: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Scaling
    const float s = sc ? 2.0f*sqrtf(0.5f*(float)ndft)/(float)ndft : 1.0f/(float)ndft;

    if (N==0u || ndft==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Init IDFT matrix
        const size_t LN = Lx * ndft;
        const float P2_N = (float)(2.0*M_PI/(double)ndft);
        float *IDFT;
        IDFT = (float *)aligned_alloc(sizeof(float),2u*LN*sizeof(float));
        if (!IDFT) { fprintf(stderr,"error in idft_cblas_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
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

        //Init complex matrix multiply
        const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
    
        if (Lx==N)
        {
            //Matrix multiply
            cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,o,IDFT,(int)Lx,X,1,z,Y,1);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=2u*ndft)
                {
                    //Matrix multiply
                    cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,o,IDFT,(int)Lx,X,1,z,Y,1);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y+=2)
                    {
                        //Matrix multiply
                        cblas_cgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,o,IDFT,(int)Lx,X,(int)K,z,Y,(int)K);
                    }
                }
            }
        }
        free(IDFT);
    }

    return 0;
}


int idft_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idft_cblas_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in idft_cblas_z: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

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
        double *IDFT;
        IDFT = (double *)aligned_alloc(sizeof(double),2u*LN*sizeof(double));
        if (!IDFT) { fprintf(stderr,"error in idft_cblas_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
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

        //Init complex matrix multiply
        const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
    
        if (Lx==N)
        {
            //Matrix multiply
            cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,o,IDFT,(int)Lx,X,1,z,Y,1);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y+=2u*ndft)
                {
                    //Matrix multiply
                    cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,o,IDFT,(int)Lx,X,1,z,Y,1);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y+=2)
                    {
                        //Matrix multiply
                        cblas_zgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,o,IDFT,(int)Lx,X,(int)K,z,Y,(int)K);
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
