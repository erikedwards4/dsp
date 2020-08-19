//CELL (~soma) stage: IIR filtering of each row or col of X according to dim.
//The IIR filters are specified by an Nx(Q+1) or (Q+1)xN matrix A,
//where Q is the IIR filter order (Q=0 means only a0; Q=1 means a0 and a1; etc.).

//The calling program must ensure that the sizes are correct, the filter is stable, etc.

//I just started this... finish later!! (Or skip.)

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int iir_s (float *Y, const float *X, const float *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);
int iir_d (double *Y, const double *X, const double *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);
int iir_c (float *Y, const float *X, const float *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);
int iir_z (double *Y, const double *X, const double *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);

int iir_inplace_s (float *X, const float *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);
int iir_inplace_d (double *X, const double *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);
int iir_inplace_c (float *X, const float *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);
int iir_inplace_z (double *X, const double *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor);


int iir_inplace_s (float *X, const float *A, const int N, const int T, const int Q, const size_t dim, const char iscolmajor)
{
    const int M = Q - 1;
    int n, t;

    //Checks
    if (N<1) { fprintf(stderr,"error in iir_s: N (num neurons) must be positive\n"); return 1; }
    if (T<1) { fprintf(stderr,"error in iir_s: T (num time points) must be positive\n"); return 1; }
    if (Q<0) { fprintf(stderr,"error in iir_s: Q (filter order) must be nonnegative\n"); return 1; }

    if (N==1)
    {
        if (A[0]!=1.0f) { cblas_sscal(T,1.0f/A[0],A,1); cblas_sscal(N*T,1.0f/A[0],X,1); }
        for (t=1; t<M; t++) { X[t] -= cblas_sdot(t,&A[M-t],1,&X[0],1); }
        for (t=M; t<T; t++) { X[t] -= cblas_sdot(M,&A[0],1,&X[t-M],1); }
    }
    else if (dim==0)
    {
        if (N<3)
        {
            if (iscolmajor)
            {
                for (n=0; n<N; n++)
                {
                    if (A[n*Q1]!=1.0f) { cblas_sscal(T,1.0f/A[n*Q1],&A[n*Q1],1); }
                }
                for (n=1; n<M; n++)
                {
                    X[n] -= cblas_sdot(n,&A[M-n],1,&X[0],1);
                    X[n+R] -= cblas_sdot(n,&A[M-n],1,&X[R],1);
                }
                for (n=M; n<R; n++)
                {
                    X[n] -= cblas_sdot(M,&A[0],1,&X[n-M],1);
                    X[n+R] -= cblas_sdot(M,&A[0],1,&X[n-M+R],1);
                }
            }
            else
            {
                for (n=1; n<M; n++)
                {
                    X[n*C] -= cblas_sdot(n,&A[M-n],1,&X[0],C);
                    X[1+n*C] -= cblas_sdot(n,&A[M-n],1,&X[1],C);
                }
                for (n=M; n<R; n++)
                {
                    X[n*C] -= cblas_sdot(M,&A[0],1,&X[(n-M)*C],C);
                    X[1+n*C] -= cblas_sdot(M,&A[0],1,&X[1+(n-M)*C],C);
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (n=1; n<M; n++) { cblas_sgemv(CblasColMajor,CblasTrans,n,C,-1.0f,&X[0],R,&A[M-n],1,1.0f,&X[n],R); }
                for (n=M; n<R; n++) { cblas_sgemv(CblasColMajor,CblasTrans,M,C,-1.0f,&X[n-M],R,&A[0],1,1.0f,&X[n],R); }
            }
            else
            {
                for (n=1; n<M; n++) { cblas_sgemv(CblasRowMajor,CblasTrans,n,C,-1.0f,&X[0],C,&A[M-n],1,1.0f,&X[n*C],1); }
                for (n=M; n<R; n++) { cblas_sgemv(CblasRowMajor,CblasTrans,M,C,-1.0f,&X[(n-M)*C],C,&A[0],1,1.0f,&X[n*C],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (R==1)
        {
            for (n=1; n<M; n++) { X[n] -= cblas_sdot(n,&A[M-n],1,&X[0],1); }
            for (n=M; n<C; n++) { X[n] -= cblas_sdot(M,&A[0],1,&X[n-M],1); }
        }
        else if (R==2)
        {
            if (iscolmajor)
            {
                for (n=1; n<M; n++)
                {
                    X[n*R] -= cblas_sdot(n,&A[M-n],1,&X[0],R);
                    X[1+n*R] -= cblas_sdot(n,&A[M-n],1,&X[1],R);
                }
                for (n=M; n<C; n++)
                {
                    X[n*R] -= cblas_sdot(M,&A[0],1,&X[(n-M)*R],R);
                    X[1+n*R] -= cblas_sdot(M,&A[0],1,&X[1+(n-M)*R],R);
                }
            }
            else
            {
                for (n=1; n<M; n++)
                {
                    X[n] -= cblas_sdot(n,&A[M-n],1,&X[0],1);
                    X[n+C] -= cblas_sdot(n,&A[M-n],1,&X[C],1);
                }
                for (n=M; n<C; n++)
                {
                    X[n] -= cblas_sdot(M,&A[0],1,&X[n-M],1);
                    X[n+C] -= cblas_sdot(M,&A[0],1,&X[n-M+C],1);
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (n=1; n<M; n++) { cblas_sgemv(CblasColMajor,CblasNoTrans,R,n,-1.0f,&X[0],R,&A[M-n],1,1.0f,&X[n*R],1); }
                for (n=M; n<C; n++) { cblas_sgemv(CblasColMajor,CblasNoTrans,R,M,-1.0f,&X[(n-M)*R],R,&A[0],1,1.0f,&X[n*R],1); }
            }
            else
            {
                for (n=1; n<M; n++) { cblas_sgemv(CblasRowMajor,CblasNoTrans,R,n,-1.0f,&X[0],C,&A[M-n],1,1.0f,&X[n],C); }
                for (n=M; n<C; n++) { cblas_sgemv(CblasRowMajor,CblasNoTrans,R,M,-1.0f,&X[n-M],C,&A[0],1,1.0f,&X[n],C); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in iir_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int iir_d (double *X, const char iscolmajor, const int R, const int C, const double *A, const int N, const size_t dim)
{
    const int M = N - 1;
    int n;

    //Checks
    if (R<1) { fprintf(stderr,"error in iir_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in iir_d: C (ncols X) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in iir_d: N (filter order) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (C==1)
        {
            for (n=1; n<M; n++) { X[n] -= cblas_ddot(n,&A[M-n],1,&X[0],1); }
            for (n=M; n<R; n++) { X[n] -= cblas_ddot(M,&A[0],1,&X[n-M],1); }
        }
        else if (C==2)
        {
            if (iscolmajor)
            {
                for (n=1; n<M; n++)
                {
                    X[n] -= cblas_ddot(n,&A[M-n],1,&X[0],1);
                    X[n+R] -= cblas_ddot(n,&A[M-n],1,&X[R],1);
                }
                for (n=M; n<R; n++)
                {
                    X[n] -= cblas_ddot(M,&A[0],1,&X[n-M],1);
                    X[n+R] -= cblas_ddot(M,&A[0],1,&X[n-M+R],1);
                }
            }
            else
            {
                for (n=1; n<M; n++)
                {
                    X[n*C] -= cblas_ddot(n,&A[M-n],1,&X[0],C);
                    X[1+n*C] -= cblas_ddot(n,&A[M-n],1,&X[1],C);
                }
                for (n=M; n<R; n++)
                {
                    X[n*C] -= cblas_ddot(M,&A[0],1,&X[(n-M)*C],C);
                    X[1+n*C] -= cblas_ddot(M,&A[0],1,&X[1+(n-M)*C],C);
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (n=1; n<M; n++) { cblas_dgemv(CblasColMajor,CblasTrans,n,C,-1.0,&X[0],R,&A[M-n],1,1.0,&X[n],R); }
                for (n=M; n<R; n++) { cblas_dgemv(CblasColMajor,CblasTrans,M,C,-1.0,&X[n-M],R,&A[0],1,1.0,&X[n],R); }
            }
            else
            {
                for (n=1; n<M; n++) { cblas_dgemv(CblasRowMajor,CblasTrans,n,C,-1.0,&X[0],C,&A[M-n],1,1.0,&X[n*C],1); }
                for (n=M; n<R; n++) { cblas_dgemv(CblasRowMajor,CblasTrans,M,C,-1.0,&X[(n-M)*C],C,&A[0],1,1.0,&X[n*C],1); }
            }
        }
    }
    else if (dim==1)
    {
        if (R==1)
        {
            for (n=1; n<M; n++) { X[n] -= cblas_ddot(n,&A[M-n],1,&X[0],1); }
            for (n=M; n<C; n++) { X[n] -= cblas_ddot(M,&A[0],1,&X[n-M],1); }
        }
        else if (R==2)
        {
            if (iscolmajor)
            {
                for (n=1; n<M; n++)
                {
                    X[n*R] -= cblas_ddot(n,&A[M-n],1,&X[0],R);
                    X[1+n*R] -= cblas_ddot(n,&A[M-n],1,&X[1],R);
                }
                for (n=M; n<C; n++)
                {
                    X[n*R] -= cblas_ddot(M,&A[0],1,&X[(n-M)*R],R);
                    X[1+n*R] -= cblas_ddot(M,&A[0],1,&X[1+(n-M)*R],R);
                }
            }
            else
            {
                for (n=1; n<M; n++)
                {
                    X[n] -= cblas_ddot(n,&A[M-n],1,&X[0],1);
                    X[n+C] -= cblas_ddot(n,&A[M-n],1,&X[C],1);
                }
                for (n=M; n<C; n++)
                {
                    X[n] -= cblas_ddot(M,&A[0],1,&X[n-M],1);
                    X[n+C] -= cblas_ddot(M,&A[0],1,&X[n-M+C],1);
                }
            }
        }
        else
        {
            if (iscolmajor)
            {
                for (n=1; n<M; n++) { cblas_dgemv(CblasColMajor,CblasNoTrans,R,n,-1.0,&X[0],R,&A[M-n],1,1.0,&X[n*R],1); }
                for (n=M; n<C; n++) { cblas_dgemv(CblasColMajor,CblasNoTrans,R,M,-1.0,&X[(n-M)*R],R,&A[0],1,1.0,&X[n*R],1); }
            }
            else
            {
                for (n=1; n<M; n++) { cblas_dgemv(CblasRowMajor,CblasNoTrans,R,n,-1.0,&X[0],C,&A[M-n],1,1.0,&X[n],C); }
                for (n=M; n<C; n++) { cblas_dgemv(CblasRowMajor,CblasNoTrans,R,M,-1.0,&X[n-M],C,&A[0],1,1.0,&X[n],C); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in iir_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int iir_c (float *X, const char iscolmajor, const int R, const int C, const float *A, const int N, const size_t dim)
{
    const float a[2] = {-1.0f,0.0f}, b[2] = {1.0f,0.0f};
    const int M = N - 1;
    int n;

    //Checks
    if (R<1) { fprintf(stderr,"error in iir_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in iir_c: C (ncols X) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in iir_c: N (filter order) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (n=1; n<M; n++) { cblas_cgemv(CblasColMajor,CblasTrans,n,C,&a[0],&X[0],R,&A[2*(M-n)],1,&b[0],&X[2*n],R); }
            for (n=M; n<R; n++) { cblas_cgemv(CblasColMajor,CblasTrans,M,C,&a[0],&X[2*(n-M)],R,&A[0],1,&b[0],&X[2*n],R); }
        }
        else
        {
            for (n=1; n<M; n++) { cblas_cgemv(CblasRowMajor,CblasTrans,n,C,&a[0],&X[0],C,&A[2*(M-n)],1,&b[0],&X[2*n*C],1); }
            for (n=M; n<R; n++) { cblas_cgemv(CblasRowMajor,CblasTrans,M,C,&a[0],&X[2*(n-M)*C],C,&A[0],1,&b[0],&X[2*n*C],1); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=1; n<M; n++) { cblas_cgemv(CblasColMajor,CblasNoTrans,R,n,&a[0],&X[0],R,&A[2*(M-n)],1,&b[0],&X[2*n*R],1); }
            for (n=M; n<C; n++) { cblas_cgemv(CblasColMajor,CblasNoTrans,R,M,&a[0],&X[2*(n-M)*R],R,&A[0],1,&b[0],&X[2*n*R],1); }
        }
        else
        {
            for (n=1; n<M; n++) { cblas_cgemv(CblasRowMajor,CblasNoTrans,R,n,&a[0],&X[0],C,&A[2*(M-n)],1,&b[0],&X[2*n],C); }
            for (n=M; n<C; n++) { cblas_cgemv(CblasRowMajor,CblasNoTrans,R,M,&a[0],&X[2*(n-M)],C,&A[0],1,&b[0],&X[2*n],C); }
        }
    }
    else
    {
        fprintf(stderr,"error in iir_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int iir_z (double *X, const char iscolmajor, const int R, const int C, const double *A, const int N, const size_t dim)
{
    const double a[2] = {-1.0,0.0}, b[2] = {1.0,0.0};
    const int M = N - 1;
    int n;

    //Checks
    if (R<1) { fprintf(stderr,"error in iir_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in iir_z: C (ncols X) must be positive\n"); return 1; }
    if (N<1) { fprintf(stderr,"error in iir_z: N (filter order) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (n=1; n<M; n++) { cblas_zgemv(CblasColMajor,CblasTrans,n,C,&a[0],&X[0],R,&A[2*(M-n)],1,&b[0],&X[2*n],R); }
            for (n=M; n<R; n++) { cblas_zgemv(CblasColMajor,CblasTrans,M,C,&a[0],&X[2*(n-M)],R,&A[0],1,&b[0],&X[2*n],R); }
        }
        else
        {
            for (n=1; n<M; n++) { cblas_zgemv(CblasRowMajor,CblasTrans,n,C,&a[0],&X[0],C,&A[2*(M-n)],1,&b[0],&X[2*n*C],1); }
            for (n=M; n<R; n++) { cblas_zgemv(CblasRowMajor,CblasTrans,M,C,&a[0],&X[2*(n-M)*C],C,&A[0],1,&b[0],&X[2*n*C],1); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (n=1; n<M; n++) { cblas_zgemv(CblasColMajor,CblasNoTrans,R,n,&a[0],&X[0],R,&A[2*(M-n)],1,&b[0],&X[2*n*R],1); }
            for (n=M; n<C; n++) { cblas_zgemv(CblasColMajor,CblasNoTrans,R,M,&a[0],&X[2*(n-M)*R],R,&A[0],1,&b[0],&X[2*n*R],1); }
        }
        else
        {
            for (n=1; n<M; n++) { cblas_zgemv(CblasRowMajor,CblasNoTrans,R,n,&a[0],&X[0],C,&A[2*(M-n)],1,&b[0],&X[2*n],C); }
            for (n=M; n<C; n++) { cblas_zgemv(CblasRowMajor,CblasNoTrans,R,M,&a[0],&X[2*(n-M)],C,&A[0],1,&b[0],&X[2*n],C); }
        }
    }
    else
    {
        fprintf(stderr,"error in iir_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
