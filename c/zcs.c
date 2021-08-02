//The signbit function was definitely slower.
//Also slower was my original while loop idea, and a few other variations.

//For complex numbers, I use zero-crossings of the imaginary part.

#include <stdio.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int zcs_s (char *Y, const float *X, const int iscolmajor, const int R, const int C, const int dim, const int going);
int zcs_d (char *Y, const double *X, const int iscolmajor, const int R, const int C, const int dim, const int going);
int zcs_c (char *Y, const float *X, const int iscolmajor, const int R, const int C, const int dim, const int going);
int zcs_z (char *Y, const double *X, const int iscolmajor, const int R, const int C, const int dim, const int going);


int zcs_s (char *Y, const float *X, const int iscolmajor, const int R, const int C, const int dim, const int going)
{
    const int N = R*C;
    int n = -1;

    //Checks
    if (R<1) { fprintf(stderr,"error in zcs_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zcs_s: C (ncols X) must be positive\n"); return 1; }

    //Raw material
    if (going>0) { while (++n<N) { Y[n] = (X[n]>=0.0f); } }
    else { while (++n<N) { Y[n] = (X[n]<0.0f); } }

    if (dim==0)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += R; }
        }
        else
        {
            if (going) { while (--n>=C) { Y[n] *= (Y[n]!=Y[n-C]); } }
            else { while (--n>=C) { Y[n] = (Y[n]!=Y[n-C]); } }
            while (n>=0) { Y[n--] = 0; }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>=R) { Y[n] *= (Y[n]!=Y[n-R]); } }
            else { while (--n>=R) { Y[n] = (Y[n]!=Y[n-R]); } }
            while (n>=0) { Y[n--] = 0; }
        }
        else
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += C; }
            //cblas_scopy(R,&z,0,&Y[0],C);
        }
    }
    else
    {
        fprintf(stderr,"error in zcs_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int zcs_d (char *Y, const double *X, const int iscolmajor, const int R, const int C, const int dim, const int going)
{
    const int N = R*C;
    int n = -1;

    //Checks
    if (R<1) { fprintf(stderr,"error in zcs_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zcs_d: C (ncols X) must be positive\n"); return 1; }

    //Raw material
    if (going>0) { while (++n<N) { Y[n] = (X[n]>=0.0); } }
    else { while (++n<N) { Y[n] = (X[n]<0.0); } }

    if (dim==0)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += R; }
        }
        else
        {
            if (going) { while (--n>=C) { Y[n] *= (Y[n]!=Y[n-C]); } }
            else { while (--n>=C) { Y[n] = (Y[n]!=Y[n-C]); } }
            while (n>=0) { Y[n--] = 0; }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>=R) { Y[n] *= (Y[n]!=Y[n-R]); } }
            else { while (--n>=R) { Y[n] = (Y[n]!=Y[n-R]); } }
            while (n>=0) { Y[n--] = 0; }
        }
        else
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += C; }
        }
    }
    else
    {
        fprintf(stderr,"error in zcs_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int zcs_c (char *Y, const float *X, const int iscolmajor, const int R, const int C, const int dim, const int going)
{
    const int N = R*C;
    int n = -1;

    //Checks
    if (R<1) { fprintf(stderr,"error in zcs_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zcs_c: C (ncols X) must be positive\n"); return 1; }

    //Raw material
    if (going>0) { while (++n<N) { Y[n] = (X[2*n+1]>=0.0f); } }
    else { while (++n<N) { Y[n] = (X[2*n+1]<0.0f); } }

    if (dim==0)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += R; }
        }
        else
        {
            if (going) { while (--n>=C) { Y[n] *= (Y[n]!=Y[n-C]); } }
            else { while (--n>=C) { Y[n] = (Y[n]!=Y[n-C]); } }
            while (n>=0) { Y[n--] = 0; }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>=R) { Y[n] *= (Y[n]!=Y[n-R]); } }
            else { while (--n>=R) { Y[n] = (Y[n]!=Y[n-R]); } }
            while (n>=0) { Y[n--] = 0; }
        }
        else
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += C; }
        }
    }
    else
    {
        fprintf(stderr,"error in zcs_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int zcs_z (char *Y, const double *X, const int iscolmajor, const int R, const int C, const int dim, const int going)
{
    const int N = R*C;
    int n = -1;

    //Checks
    if (R<1) { fprintf(stderr,"error in zcs_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zcs_z: C (ncols X) must be positive\n"); return 1; }

    //Raw material
    if (going>0) { while (++n<N) { Y[n] = (X[2*n+1]>=0.0); } }
    else { while (++n<N) { Y[n] = (X[2*n+1]<0.0); } }

    if (dim==0)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += R; }
        }
        else
        {
            if (going) { while (--n>=C) { Y[n] *= (Y[n]!=Y[n-C]); } }
            else { while (--n>=C) { Y[n] = (Y[n]!=Y[n-C]); } }
            while (n>=0) { Y[n--] = 0; }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            if (going) { while (--n>=R) { Y[n] *= (Y[n]!=Y[n-R]); } }
            else { while (--n>=R) { Y[n] = (Y[n]!=Y[n-R]); } }
            while (n>=0) { Y[n--] = 0; }
        }
        else
        {
            if (going) { while (--n>0) { Y[n] *= (Y[n]!=Y[n-1]); } }
            else { while (--n>0) { Y[n] = (Y[n]!=Y[n-1]); } }
            while (n<N) { Y[n] = 0; n += C; }
        }
    }
    else
    {
        fprintf(stderr,"error in zcs_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

