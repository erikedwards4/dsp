//Applies a window W to each row or col of X according to dim
//This is just element-wise multiplication for each row or col.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int apply_win_s (float *X, const char iscolmajor, const int R, const int C, const float *W, const int dim);
int apply_win_d (double *X, const char iscolmajor, const int R, const int C, const double *W, const int dim);
int apply_win_c (float *X, const char iscolmajor, const int R, const int C, const float *W, const int dim);
int apply_win_z (double *X, const char iscolmajor, const int R, const int C, const double *W, const int dim);


int apply_win_s (float *X, const char iscolmajor, const int R, const int C, const float *W, const int dim)
{
    int l;

    //Checks
    if (R<1) { fprintf(stderr,"error in apply_win_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in apply_win_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (l=0; l<C; l++) { cblas_sscal(R,W[l],&X[l*C],1); }
        }
        else
        {
            for (l=0; l<C; l++) { cblas_sscal(R,W[l],&X[l],C); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (l=0; l<R; l++) { cblas_sscal(C,W[l],&X[l],R); }
        }
        else
        {
            for (l=0; l<R; l++) { cblas_sscal(C,W[l],&X[l*R],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in apply_win_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int apply_win_d (double *X, const char iscolmajor, const int R, const int C, const double *W, const int dim)
{
    int l;

    //Checks
    if (R<1) { fprintf(stderr,"error in apply_win_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in apply_win_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (l=0; l<C; l++) { cblas_dscal(R,W[l],&X[l*C],1); }
        }
        else
        {
            for (l=0; l<C; l++) { cblas_dscal(R,W[l],&X[l],C); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (l=0; l<R; l++) { cblas_dscal(C,W[l],&X[l],R); }
        }
        else
        {
            for (l=0; l<R; l++) { cblas_dscal(C,W[l],&X[l*R],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in apply_win_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int apply_win_c (float *X, const char iscolmajor, const int R, const int C, const float *W, const int dim)
{
    int l;

    //Checks
    if (R<1) { fprintf(stderr,"error in apply_win_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in apply_win_c: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (l=0; l<2*C; l+=2) { cblas_cscal(R,&W[l],&X[l*C],1); }
        }
        else
        {
            for (l=0; l<2*C; l+=2) { cblas_cscal(R,&W[l],&X[l],C); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (l=0; l<2*R; l+=2) { cblas_cscal(C,&W[l],&X[l],R); }
        }
        else
        {
            for (l=0; l<2*R; l+=2) { cblas_cscal(C,&W[l],&X[l*R],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in apply_win_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int apply_win_z (double *X, const char iscolmajor, const int R, const int C, const double *W, const int dim)
{
    int l;

    //Checks
    if (R<1) { fprintf(stderr,"error in apply_win_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in apply_win_z: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (l=0; l<2*C; l+=2) { cblas_zscal(R,&W[l],&X[l*C],1); }
        }
        else
        {
            for (l=0; l<2*C; l+=2) { cblas_zscal(R,&W[l],&X[l],C); }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (l=0; l<2*R; l+=2) { cblas_zscal(C,&W[l],&X[l],R); }
        }
        else
        {
            for (l=0; l<2*R; l+=2) { cblas_zscal(C,&W[l],&X[l*R],1); }
        }
    }
    else
    {
        fprintf(stderr,"error in apply_win_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

