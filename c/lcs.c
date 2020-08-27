//Level-crossings

#include <stdio.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int lcs_s (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going, const float level);
int lcs_d (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going, const double level);


int lcs_s (char *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int going, const float level)
{
    const int N = R*C;
    int n = -1;

    //Checks
    if (R<1) { fprintf(stderr,"error in lcs_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lcs_s: C (ncols X) must be positive\n"); return 1; }

    //Raw material
    if (going>0) { while (++n<N) { Y[n] = (X[n]>=level); } }
    else { while (++n<N) { Y[n] = (X[n]<level); } }

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
        fprintf(stderr,"error in lcs_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int lcs_d (char *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int going, const double level)
{
    const int N = R*C;
    int n = -1;

    //Checks
    if (R<1) { fprintf(stderr,"error in lcs_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in lcs_d: C (ncols X) must be positive\n"); return 1; }

    //Raw material
    if (going>0) { while (++n<N) { Y[n] = (X[n]>=level); } }
    else { while (++n<N) { Y[n] = (X[n]<level); } }

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
        fprintf(stderr,"error in lcs_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif

