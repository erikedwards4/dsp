//Causal FIR filtering of each row or col of X according to dim.

//FIR impulse responses are given in matrix B.
//For dim=0, X is NxT and B is NxL
//For dim=1, X is TxN and B is LxN
//where N is the number of neurons, T the number of time points, and L the FIR filter order.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#ifdef __cplusplus
namespace openn {
extern "C" {
#endif

int fir_s (float *Y, const float *X, const float *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor);
int fir_d (double *Y, const double *X, const double *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor);
int fir_c (float *Y, const float *X, const float *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor);
int fir_z (double *Y, const double *X, const double *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor);


int fir_s (float *Y, const float *X, const float *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor)
{
    const float z = 0.0f;
    int n, l;

    //Checks
    if (N<1) { fprintf(stderr,"error in fir_s: N (num neurons) must be positive\n"); return 1; }
    if (T<1) { fprintf(stderr,"error in fir_s: T (num time points) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_s: L (filter IR length) must be positive\n"); return 1; }
    
    //Initialize Y to 0
    cblas_scopy(N*T,&z,0,Y,1);

    if (N==1)
    {
        for (l=0; l<L; l++) { cblas_saxpy(T-l,B[l],&X[0],1,&Y[l],1); }
    }
    else if (dim==0)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_saxpy(T-l,B[n+l*N],&X[n],N,&Y[n+l*N],N); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_saxpy(T-l,B[n*L+l],&X[n*T],1,&Y[n*T+l],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_saxpy(T-l,B[n*L+l],&X[n*T],1,&Y[n*T+l],1); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_saxpy(T-l,B[n+l*N],&X[n],N,&Y[n+l*N],N); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int fir_d (double *Y, const double *X, const double *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor)
{
    const double z = 0.0;
    int n, l;

    //Checks
    if (N<1) { fprintf(stderr,"error in fir_d: N (num neurons) must be positive\n"); return 1; }
    if (T<1) { fprintf(stderr,"error in fir_d: T (num time points) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_d: L (filter IR length) must be positive\n"); return 1; }
    
    //Initialize Y to 0
    cblas_dcopy(N*T,&z,0,Y,1);

    if (N==1)
    {
        for (l=0; l<L; l++) { cblas_daxpy(T-l,B[l],&X[0],1,&Y[l],1); }
    }
    else if (dim==0)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_daxpy(T-l,B[n+l*N],&X[n],N,&Y[n+l*N],N); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_daxpy(T-l,B[n*L+l],&X[n*T],1,&Y[n*T+l],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_daxpy(T-l,B[n*L+l],&X[n*T],1,&Y[n*T+l],1); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_daxpy(T-l,B[n+l*N],&X[n],N,&Y[n+l*N],N); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int fir_c (float *Y, const float *X, const float *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor)
{
    const float z[2] = {0.0f,0.0f};
    int n, l;

    //Checks
    if (N<1) { fprintf(stderr,"error in fir_c: N (num neurons) must be positive\n"); return 1; }
    if (T<1) { fprintf(stderr,"error in fir_c: T (num time points) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_c: L (filter IR length) must be positive\n"); return 1; }
    
    //Initialize Y to 0
    cblas_ccopy(N*T,z,0,Y,1);

    if (N==1)
    {
        for (l=0; l<L; l++) { cblas_caxpy(T-l,&B[2*l],&X[0],1,&Y[2*l],1); }
    }
    else if (dim==0)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_caxpy(T-l,&B[2*(n+l*N)],&X[2*n],N,&Y[2*(n+l*N)],N); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_caxpy(T-l,&B[2*(n*L+l)],&X[2*n*T],1,&Y[2*(n*T+l)],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_caxpy(T-l,&B[2*(n*L+l)],&X[2*n*T],1,&Y[2*(n*T+l)],1); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_caxpy(T-l,&B[2*(n+l*N)],&X[2*n],N,&Y[2*(n+l*N)],N); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int fir_z (double *Y, const double *X, const double *B, const int N, const int T, const int L, const size_t dim, const char iscolmajor)
{
    const double z[2] = {0.0,0.0};
    int n, l;

    //Checks
    if (N<1) { fprintf(stderr,"error in fir_z: N (num neurons) must be positive\n"); return 1; }
    if (T<1) { fprintf(stderr,"error in fir_z: T (num time points) must be positive\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in fir_z: L (filter IR length) must be positive\n"); return 1; }
    
    //Initialize Y to 0
    cblas_zcopy(N*T,z,0,Y,1);

    if (N==1)
    {
        for (l=0; l<L; l++) { cblas_zaxpy(T-l,&B[2*l],&X[0],1,&Y[2*l],1); }
    }
    else if (dim==0)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(T-l,&B[2*(n+l*N)],&X[2*n],N,&Y[2*(n+l*N)],N); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(T-l,&B[2*(n*L+l)],&X[2*n*T],1,&Y[2*(n*T+l)],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(T-l,&B[2*(n*L+l)],&X[2*n*T],1,&Y[2*(n*T+l)],1); }
            }
        }
        else
        {
            for (n=0; n<N; n++)
            {
                for (l=0; l<L; l++) { cblas_zaxpy(T-l,&B[2*(n+l*N)],&X[2*n],N,&Y[2*(n+l*N)],N); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fir_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
