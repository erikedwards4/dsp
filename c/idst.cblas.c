//This computes "the" DST (discrete sine transformation) along dim of matrix X.
//This uses CBLAS matrix multiplication by the DST-I matrix.
//This is the same as DST-I except scaling.

//For complex input X, the output Y is complex, and the DST is just the DST of the
//real and imaginary parts separately (following Octave convention).

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int idst_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);
int idst_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc);


int idst_cblas_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_cblas_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_cblas_s: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float s = sc ? 2.0f/(float)(ndst+1u) : 2.0f/(float)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const float P_N = (float)(M_PI/(double)(ndst+1u));
        float *DST;
        DST = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DST) { fprintf(stderr,"error in idst_cblas_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sinf(P_N*(float)(l+1u)*(float)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X,1,0.0f,Y,1);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                cblas_sgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)ndst,(int)V,(int)Lx,1.0f,DST,(int)Lx,X,(int)Lx,0.0f,Y,(int)ndst);
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, ++Y)
                    {
                        cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X,(int)K,0.0f,Y,(int)K);
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


int idst_cblas_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_cblas_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_cblas_d: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const double s = sc ? 2.0/(double)(ndst+1u) : 2.0/(double)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const double P_N = M_PI/(double)(ndst+1u);
        double *DST;
        DST = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
        if (!DST) { fprintf(stderr,"error in idst_cblas_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sin(P_N*(double)(l+1u)*(double)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X,1,0.0,Y,1);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                cblas_dgemm(CblasColMajor,CblasTrans,CblasNoTrans,(int)ndst,(int)V,(int)Lx,1.0,DST,(int)Lx,X,(int)Lx,0.0,Y,(int)ndst);
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, ++Y)
                    {
                        cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X,(int)K,0.0,Y,(int)K);
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


int idst_cblas_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_cblas_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_cblas_c: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const float s = sc ? 2.0f/(float)(ndst+1u) : 2.0f/(float)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const float P_N = (float)(M_PI/(double)(ndst+1u));
        float *DST;
        DST = (float *)aligned_alloc(sizeof(float),LN*sizeof(float));
        if (!DST) { fprintf(stderr,"error in idst_cblas_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sinf(P_N*(float)(l+1u)*(float)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X,2,0.0f,Y,2);
            cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X+1,2,0.0f,Y+1,2);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*Lx, Y+=2u*ndst)
                {
                    cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X,2,0.0f,Y,2);
                    cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X+1,2,0.0f,Y+1,2);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, ++Y)
                    {
                        cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X,2*(int)K,0.0f,Y,2*(int)K);
                        ++X; ++Y;
                        cblas_sgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0f,DST,(int)Lx,X,2*(int)K,0.0f,Y,2*(int)K);
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


int idst_cblas_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndst, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in idst_cblas_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndst<Lx) { fprintf(stderr,"error in idst_cblas_z: ndst must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndst==0u || N==0u) {}
    else if (ndst==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Scaling
        const double s = sc ? 2.0/(double)(ndst+1u) : 2.0/(double)(2u*ndst+2u);

        //Initialize DST-I matrix
        const size_t LN = Lx * ndst;
        const double P_N = M_PI/(double)(ndst+1u);
        double *DST;
        DST = (double *)aligned_alloc(sizeof(double),LN*sizeof(double));
        if (!DST) { fprintf(stderr,"error in idst_cblas_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        for (size_t n=0u; n<ndst; ++n)
        {
            for (size_t l=0u; l<Lx; ++l, ++DST)
            {
                *DST = s * sin(P_N*(double)(l+1u)*(double)(n+1u));
            }
        }
        DST -= LN;
    
        if (Lx==N)
        {
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X,2,0.0,Y,2);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X+1,2,0.0,Y+1,2);
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, X+=2u*Lx, Y+=2u*ndst)
                {
                    cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X,2,0.0,Y,2);
                    cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X+1,2,0.0,Y+1,2);
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndst-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, ++Y)
                    {
                        cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X,2*(int)K,0.0,Y,2*(int)K);
                        ++X; ++Y;
                        cblas_dgemv(CblasRowMajor,CblasNoTrans,(int)ndst,(int)Lx,1.0,DST,(int)Lx,X,2*(int)K,0.0,Y,2*(int)K);
                    }
                }
            }
        }
        free(DST);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
