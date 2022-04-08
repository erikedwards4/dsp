//This computes "the" IDCT (inverse discrete cosine transformation) along dim of matrix X.
//This uses a CBLAS matrix multiplication by the DCT-III matrix.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include "codee_dsp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int idct_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_cblas_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in idct_cblas_s: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float xsc = (sc) ? (float)(2.0*M_SQRT1_2) : 1.0f;
        const float ysc = (sc) ? 2.0f/sqrtf((float)(2u*ndct)) : 2.0f/(float)(2u*ndct);
        const float dcsc = 0.5f * xsc * ysc;

        //Initialize DCT-III matrix
        const size_t LN = Lx * ndct;
        const float P_N = (float)(M_PI/(double)ndct);
        float *DCT, x0;
        DCT = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DCT) { fprintf(stderr,"error in idct_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndct; ++n)
        {
            *DCT++ = xsc * ysc;
            for (size_t l=1u; l<Lx; ++l, ++DCT)
            {
                *DCT = ysc * cosf(P_N*(0.5f+(float)n)*(float)l);
            }
        }
        DCT -= LN;
    
        if (Lx==N)
        {
            x0 = dcsc * *X;
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X,1,0.0f,Y,1);
            for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0; }
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
                    x0 = dcsc * *X;
                    cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X,1,0.0f,Y,1);
                    for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndct-1u)
                    {
                        x0 = dcsc * *X;
                        cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X,(int)K,0.0f,Y,(int)K);
                        for (size_t n=ndct; n>0u; --n, Y+=K) { *Y -= x0; }
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


int idct_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_cblas_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in idct_cblas_d: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
     else
    {
        //Scaling
        const double xsc = (sc) ? 2.0*M_SQRT1_2 : 1.0;
        const double ysc = (sc) ? 2.0/sqrt((double)(2u*ndct)) : 2.0/(double)(2u*ndct);
        const double dcsc = 0.5 * xsc * ysc;

        //Initialize DCT-III matrix
        const size_t LN = Lx * ndct;
        const double P_N = M_PI/(double)ndct;
        double *DCT, x0;
        DCT = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
        if (!DCT) { fprintf(stderr,"error in idct_cblas_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndct; ++n)
        {
            *DCT++ = xsc * ysc;
            for (size_t l=1u; l<Lx; ++l, ++DCT)
            {
                *DCT = ysc * cos(P_N*(0.5+(double)n)*(double)l);
            }
        }
        DCT -= LN;
    
        if (Lx==N)
        {
            x0 = dcsc * *X;
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X,1,0.0,Y,1);
            for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0; }
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
                    x0 = dcsc * *X;
                    cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X,1,0.0,Y,1);
                    for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndct-1u)
                    {
                        x0 = dcsc * *X;
                        cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X,(int)K,0.0,Y,(int)K);
                        for (size_t n=ndct; n>0u; --n, Y+=K) { *Y -= x0; }
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


int idct_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_cblas_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in idct_cblas_c: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float xsc = (sc) ? (float)(2.0*M_SQRT1_2) : 1.0f;
        const float ysc = (sc) ? 2.0f/sqrtf((float)(2u*ndct)) : 2.0f/(float)(2u*ndct);
        const float dcsc = 0.5f * xsc * ysc;

        //Initialize DCT-III matrix
        const size_t LN = Lx * ndct;
        const float P_N = (float)(M_PI/(double)ndct);
        float *DCT, x0r, x0i;
        DCT = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DCT) { fprintf(stderr,"error in idct_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndct; ++n)
        {
            *DCT++ = xsc * ysc;
            for (size_t l=1u; l<Lx; ++l, ++DCT)
            {
                *DCT = ysc * cosf(P_N*(0.5f+(float)n)*(float)l);
            }
        }
        DCT -= LN;
    
        if (Lx==N)
        {
            x0r = dcsc * *X; x0i = dcsc * *(X+1);
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X,2,0.0f,Y,2);
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X+1,2,0.0f,Y+1,2);
            for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0r; *++Y -= x0i; }
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
                    x0r = dcsc * *X; x0i = dcsc * *(X+1);
                    cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X,2,0.0f,Y,2);
                    cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X+1,2,0.0f,Y+1,2);
                    for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0r; *++Y -= x0i; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*ndct-2u)
                    {
                        x0r = dcsc * *X; x0i = dcsc * *(X+1);
                        cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X,2*(int)K,0.0f,Y,2*(int)K);
                        cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0f,DCT,(int)Lx,X+1,2*(int)K,0.0f,Y+1,2*(int)K);
                        for (size_t n=ndct; n>0u; --n, Y+=2u*K) { *Y -= x0r; *(Y+1) -= x0i; }
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


int idct_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idct_cblas_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in idct_cblas_z: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const double xsc = (sc) ? 2.0*M_SQRT1_2 : 1.0;
        const double ysc = (sc) ? 2.0/sqrt((double)(2u*ndct)) : 2.0/(double)(2u*ndct);
        const double dcsc = 0.5 * xsc * ysc;

        //Initialize DCT-III matrix
        const size_t LN = Lx * ndct;
        const double P_N = M_PI/(double)ndct;
        double *DCT, x0r, x0i;
        DCT = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
        if (!DCT) { fprintf(stderr,"error in idct_cblas_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndct; ++n)
        {
            *DCT++ = xsc * ysc;
            for (size_t l=1u; l<Lx; ++l, ++DCT)
            {
                *DCT = ysc * cos(P_N*(0.5+(double)n)*(double)l);
            }
        }
        DCT -= LN;
    
        if (Lx==N)
        {
            x0r = dcsc * *X; x0i = dcsc * *(X+1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X,2,0.0,Y,2);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X+1,2,0.0,Y+1,2);
            for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0r; *++Y -= x0i; }
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
                    x0r = dcsc * *X; x0i = dcsc * *(X+1);
                    cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X,2,0.0,Y,2);
                    cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X+1,2,0.0,Y+1,2);
                    for (size_t n=ndct; n>0u; --n, ++Y) { *Y -= x0r; *++Y -= x0i; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*ndct-2u)
                    {
                        x0r = dcsc * *X; x0i = dcsc * *(X+1);
                        cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X,2*(int)K,0.0,Y,2*(int)K);
                        cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndct,(int)Lx,1.0,DCT,(int)Lx,X+1,2*(int)K,0.0,Y+1,2*(int)K);
                        for (size_t n=ndct; n>0u; --n, Y+=2u*K) { *Y -= x0r; *(Y+1) -= x0i; }
                    }
                }
            }
        }
        free(DCT);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
