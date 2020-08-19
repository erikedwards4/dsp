//The "unbiased" version uses N-l in the denominator instead of N.
//It is actually just "less biased", but is slower,
//has larger mean-squared error, and doesn't match FFT estimate.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2ac_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);
int sig2ac_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);
int sig2ac_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);
int sig2ac_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased);


int sig2ac_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
{
    if (dim>3) { fprintf(stderr,"error in sig2ac_s: dim must be in [0 3]\n"); return 1; }

    float *X1;  //1 row or col of X

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in sig2ac_s: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)R*sizeof(float)))) { fprintf(stderr,"error in sig2ac_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy((int)R,&X[c*R],1,&X1[0],1);
                //m = cblas_sdot((int)R,&X1[0],1,&o,0) / R;
                //cblas_saxpy((int)R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+c*L] = cblas_sdot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal((int)C,(float)R/(R-l),&Y[l],L); } }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy((int)R,&X[c],(int)C,&X1[0],1);
                //m = cblas_sdot((int)R,&X1[0],1,&o,0) / R;
                //cblas_saxpy((int)R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[c+l*C] = cblas_sdot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal((int)C,(float)R/(R-l),&Y[l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in sig2ac_s: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)C*sizeof(float)))) { fprintf(stderr,"error in sig2ac_s: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy((int)C,&X[r],(int)R,&X1[0],1);
                //m = cblas_sdot((int)C,&X1[0],1,&o,0) / C;
                //cblas_saxpy((int)C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[r+l*R] = cblas_sdot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal((int)R,(float)C/(C-l),&Y[l*R],1); } }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy((int)C,&X[r*C],1,&X1[0],1);
                //m = cblas_sdot((int)C,&X1[0],1,&o,0) / C;
                //cblas_saxpy((int)C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+r*L] = cblas_sdot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_sscal((int)R,(float)C/(C-l),&Y[l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int sig2ac_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
{
    if (dim>3) { fprintf(stderr,"error in sig2ac_d: dim must be in [0 3]\n"); return 1; }

    double *X1;  //1 row or col of X

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in sig2ac_d: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)R*sizeof(double)))) { fprintf(stderr,"error in sig2ac_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy((int)R,&X[c*R],1,&X1[0],1);
                //m = cblas_ddot((int)R,&X1[0],1,&o,0) / R;
                //cblas_daxpy((int)R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+c*L] = cblas_ddot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal((int)C,(double)R/(R-l),&Y[l],L); } }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy((int)R,&X[c],(int)C,&X1[0],1);
                //m = cblas_ddot((int)R,&X1[0],1,&o,0) / R;
                //cblas_daxpy((int)R,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[c+l*C] = cblas_ddot(R-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal((int)C,(double)R/(R-l),&Y[l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in sig2ac_d: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)C*sizeof(double)))) { fprintf(stderr,"error in sig2ac_d: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy((int)C,&X[r],(int)R,&X1[0],1);
                //m = cblas_ddot((int)C,&X1[0],1,&o,0) / C;
                //cblas_daxpy((int)C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[r+l*R] = cblas_ddot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal((int)R,(double)C/(C-l),&Y[l*R],1); } }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy((int)C,&X[r*C],1,&X1[0],1);
                //m = cblas_ddot((int)C,&X1[0],1,&o,0) / C;
                //cblas_daxpy((int)C,-m,&o,0,&X1[0],1);
                for (l=0; l<L; l++) { Y[l+r*L] = cblas_ddot(C-l,&X1[0],1,&X1[l],1); }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_dscal((int)R,(double)C/(C-l),&Y[l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int sig2ac_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
{
    if (dim>3) { fprintf(stderr,"error in sig2ac_c: dim must be in [0 3]\n"); return 1; }

    _Complex float dot;
    float *X1;  //1 row or col of X

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in sig2ac_c: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)(2*R)*sizeof(float)))) { fprintf(stderr,"error in sig2ac_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_ccopy((int)R,&X[2*c*R],1,&X1[0],1);
                //m = -cblas_cdotu((int)R,&X1[0],1,&o[0],0) / R;
                //cblas_caxpy((int)R,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(l+c*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal((int)C,(float)R/(R-l),&Y[2*l],L); } }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_ccopy((int)R,&X[2*c],2*C,&X1[0],1);
                //m = -cblas_cdotu((int)R,&X1[0],1,&o[0],0) / R;
                //cblas_caxpy((int)R,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(c+l*C)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal((int)C,(float)R/(R-l),&Y[2*l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in sig2ac_c: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(float *)malloc((size_t)(2*C)*sizeof(float)))) { fprintf(stderr,"error in sig2ac_c: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy((int)C,&X[2*r],2*R,&X1[0],1);
                //m = -cblas_cdotu((int)C,&X1[0],1,&o[0],0) / C;
                //cblas_caxpy((int)C,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(r+l*R)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal((int)R,(float)C/(C-l),&Y[2*l*R],1); } }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_ccopy((int)C,&X[2*r*C],1,&X1[0],1);
                //m = -cblas_cdotu((int)C,&X1[0],1,&o[0],0) / C;
                //cblas_caxpy((int)C,(float *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_cdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_ccopy(1,(float *)&dot,1,&Y[2*(l+r*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_csscal((int)R,(float)C/(C-l),&Y[2*l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int sig2ac_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const int L, const char unbiased)
{
    if (dim>3) { fprintf(stderr,"error in sig2ac_z: dim must be in [0 3]\n"); return 1; }

    _Complex double dot;
    double *X1;  //1 row or col of X

    if (dim==0)
    {
        if (L>R) { fprintf(stderr,"error in sig2ac_z: nlags must be < nrows X for dim==0\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)(2*R)*sizeof(double)))) { fprintf(stderr,"error in sig2ac_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_zcopy((int)R,&X[2*c*R],1,&X1[0],1);
                //m = -cblas_zdotu((int)R,&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy((int)R,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(l+c*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal((int)C,(double)R/(R-l),&Y[2*l],L); } }
        }
        else
        {
            for (size_t c=0; c<C; c++)
            {
                cblas_zcopy((int)R,&X[2*c],2*C,&X1[0],1);
                //m = -cblas_zdotu((int)R,&X1[0],1,&o[0],0) / R;
                //cblas_zaxpy((int)R,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(R-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(c+l*C)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal((int)C,(double)R/(R-l),&Y[2*l*C],1); } }
        }
    }
    else if (dim==1)
    {
        if (L>C) { fprintf(stderr,"error in sig2ac_z: nlags must be < ncols X for dim==1\n"); return 1; }
        if (!(X1=(double *)malloc((size_t)(2*C)*sizeof(double)))) { fprintf(stderr,"error in sig2ac_z: problem with malloc. "); perror("malloc"); return 1; }
        if (iscolmajor)
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy((int)C,&X[2*r],2*R,&X1[0],1);
                //m = -cblas_zdotu((int)C,&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy((int)C,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(r+l*R)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal((int)R,(double)C/(C-l),&Y[2*l*R],1); } }
        }
        else
        {
            for (size_t r=0; r<R; r++)
            {
                cblas_zcopy((int)C,&X[2*r*C],1,&X1[0],1);
                //m = -cblas_zdotu((int)C,&X1[0],1,&o[0],0) / C;
                //cblas_zaxpy((int)C,(double *)&m,&o[0],0,&X1[0],1);
                for (l=0; l<L; l++)
                {
                    dot = cblas_zdotu(C-l,&X1[0],1,&X1[2*l],1);
                    cblas_zcopy(1,(double *)&dot,1,&Y[2*(l+r*L)],1);
                }
            }
            if (unbiased) { for (l=0; l<L; l++) { cblas_zdscal((int)R,(double)C/(C-l),&Y[2*l],L); } }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2ac_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
