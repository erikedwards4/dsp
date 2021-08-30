//This computes the DFT (discrete Fourier transformation) along dim of matrix X.
//This uses a CBLAS matrix multiplication by the DFT matrix.
//This can be useful for smaller transform sizes, especially odd-length ones.

//This can also be useful if using only a few of the output frequencies,
//in which case a simple modification of this script would work.
//Especially useful (with modification) if can be combined with another matrix multiply.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int dft_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int dft_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int dft_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);
int dft_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc);


int dft_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dft_cblas_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndft<Lx) { fprintf(stderr,"error in dft_cblas_s: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

    if (N==0u) {}
    else if (ndft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        Y -= 2u*N;
    }
    else
    {
        //Scaling
        const float s = sc ? 2.0f/sqrtf((float)(2u*ndft)) : 2.0f;
        const float dcsc = sc ? 1.0f/sqrtf((float)ndft) : 2.0f;

        //Initialize DFT matrices (real, imag parts separate)
        const size_t LN = Lx * ndft;
        const float P_N = (float)(M_PI/(double)ndft);
        float *DFTr, *DFTi;
        DFTr = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DFTr) { fprintf(stderr,"error in dft_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        DFTi = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DFTi) { fprintf(stderr,"error in dft_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t l=0u; l<Lx; ++l, ++DFTr) { *DFTr = dcsc; }
        for (size_t n=1u; n<ndft; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DFTr)
            {
                *DFTr = s * cosf(P_N*(0.5f+(float)l)*(float)n);
            }
        }
        DFTr -= LN;
    
        if (Lx==N)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X,1,0.0f,Y,1);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)ndft,(int)V,(int)Lx,1.0f,DFT,(int)Lx,X,(int)Lx,0.0f,Y,(int)ndft);
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndft-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, ++Y)
                    {
                        cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X,(int)K,0.0f,Y,(int)K);
                    }
                }
            }
        }
        free(DFT);
    }

    return 0;
}


// int dft_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
// {
//     if (dim>3u) { fprintf(stderr,"error in dft_cblas_d: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (ndft<Lx) { fprintf(stderr,"error in dft_cblas_d: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

//     if (ndft==0u || N==0u) {}
//     else if (ndft==1u)
//     {
//         for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
//     }
//     else
//     {
//         //Scaling
//         const double s = sc ? 2.0/sqrt((double)(2u*ndft)) : 2.0;
//         const double dcsc = sc ? 1.0/sqrt((double)ndft) : 2.0;

//         //Initialize DFT-II matrix
//         const size_t LN = Lx * ndft;
//         const double P_N = M_PI/(double)ndft;
//         double *DFT;
//         DFT = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
//         if (!DFT) { fprintf(stderr,"error in dft_cblas_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
//         for (size_t l=0u; l<Lx; ++l, ++DFT) { *DFT = dcsc; }
//         for (size_t n=1u; n<ndft; ++n)
//         {
//             for (size_t l=0u; l<Lx; ++l, ++DFT)
//             {
//                 *DFT = s * cos(P_N*(0.5+(double)l)*(double)n);
//             }
//         }
//         DFT -= LN;
    
//         if (Lx==N)
//         {
//             cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X,1,0.0,Y,1);
//         }
//         else
//         {
//             const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)ndft,(int)V,(int)Lx,1.0,DFT,(int)Lx,X,(int)Lx,0.0,Y,(int)ndft);
//             }
//             else
//             {
//                 for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndft-1u))
//                 {
//                     for (size_t b=B; b>0u; --b, ++X, ++Y)
//                     {
//                         cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X,(int)K,0.0,Y,(int)K);
//                     }
//                 }
//             }
//         }
//         free(DFT);
//     }

//     return 0;
// }


// int dft_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
// {
//     if (dim>3u) { fprintf(stderr,"error in dft_cblas_c: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (ndft<Lx) { fprintf(stderr,"error in dft_cblas_c: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

//     if (ndft==0u || N==0u) {}
//     else if (ndft==1u)
//     {
//         for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
//     }
//     else
//     {
//         //Scaling
//         const float s = sc ? 2.0f/sqrtf((float)(2u*ndft)) : 2.0f;
//         const float dcsc = sc ? 1.0f/sqrtf((float)ndft) : 2.0f;

//         //Initialize DFT-II matrix
//         const size_t LN = Lx * ndft;
//         const float P_N = (float)(M_PI/(double)ndft);
//         float *DFT;
//         DFT = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
//         if (!DFT) { fprintf(stderr,"error in dft_cblas_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
//         for (size_t l=0u; l<Lx; ++l, ++DFT) { *DFT = dcsc; }
//         for (size_t n=1u; n<ndft; ++n)
//         {
//             for (size_t l=0u; l<Lx; ++l, ++DFT)
//             {
//                 *DFT = s * cosf(P_N*(0.5f+(float)l)*(float)n);
//             }
//         }
//         DFT -= LN;
    
//         if (Lx==N)
//         {
//             cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X,2,0.0f,Y,2);
//             cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X+1,2,0.0f,Y+1,2);
//         }
//         else
//         {
//             const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=2u*Lx, Y+=2u*ndft)
//                 {
//                     cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X,2,0.0f,Y,2);
//                     cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X+1,2,0.0f,Y+1,2);
//                 }
//             }
//             else
//             {
//                 for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndft-1u))
//                 {
//                     for (size_t b=B; b>0u; --b, ++X, ++Y)
//                     {
//                         cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X,2*(int)K,0.0f,Y,2*(int)K);
//                         ++X; ++Y;
//                         cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0f,DFT,(int)Lx,X,2*(int)K,0.0f,Y,2*(int)K);
//                     }
//                 }
//             }
//         }
//         free(DFT);
//     }

//     return 0;
// }


// int dft_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndft, const int sc)
// {
//     if (dim>3u) { fprintf(stderr,"error in dft_cblas_z: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     if (ndft<Lx) { fprintf(stderr,"error in dft_cblas_z: ndft must be >= Lx (length of vecs in X)\n"); return 1; }

//     if (ndft==0u || N==0u) {}
//     else if (ndft==1u)
//     {
//         for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
//     }
//     else
//     {
//         //Scaling
//         const double s = sc ? 2.0/sqrt((double)(2u*ndft)) : 2.0;
//         const double dcsc = sc ? 1.0/sqrt((double)ndft) : 2.0;

//         //Initialize DFT-II matrix
//         const size_t LN = Lx * ndft;
//         const double P_N = M_PI/(double)ndft;
//         double *DFT;
//         DFT = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
//         if (!DFT) { fprintf(stderr,"error in dft_cblas_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
//         for (size_t l=0u; l<Lx; ++l, ++DFT) { *DFT = dcsc; }
//         for (size_t n=1u; n<ndft; ++n)
//         {
//             for (size_t l=0u; l<Lx; ++l, ++DFT)
//             {
//                 *DFT = s * cos(P_N*(0.5+(double)l)*(double)n);
//             }
//         }
//         DFT -= LN;
    
//         if (Lx==N)
//         {
//             cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X,2,0.0,Y,2);
//             cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X+1,2,0.0,Y+1,2);
//         }
//         else
//         {
//             const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
//             const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
//             const size_t V = N/Lx, G = V/B;

//             if (K==1u && (G==1u || B==1u))
//             {
//                 for (size_t v=0u; v<V; ++v, X+=2u*Lx, Y+=2u*ndft)
//                 {
//                     cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X,2,0.0,Y,2);
//                     cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X+1,2,0.0,Y+1,2);
//                 }
//             }
//             else
//             {
//                 for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndft-1u))
//                 {
//                     for (size_t b=B; b>0u; --b, ++X, ++Y)
//                     {
//                         cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X,2*(int)K,0.0,Y,2*(int)K);
//                         ++X; ++Y;
//                         cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndft,(int)Lx,1.0,DFT,(int)Lx,X,2*(int)K,0.0,Y,2*(int)K);
//                     }
//                 }
//             }
//         }
//         free(DFT);
//     }

//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
