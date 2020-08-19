//Normalizes standard deviation to 1 for each row or col according to dim.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int stdev1_s (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);
int stdev1_d (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);
int stdev1_c (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);
int stdev1_z (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);


int stdev1_s (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    float s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in stdev1_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stdev1_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrtf(R); }
        else { den = sqrtf(R-1); }

        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                s = den / cblas_snrm2(R,&X[c*R],1);
                cblas_sscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                s = den / cblas_snrm2(R,&X[c],C);
                cblas_sscal(R,s,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = sqrtf(C); }
        else { den = sqrtf(C-1); }
        
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                s = den / cblas_snrm2(C,&X[r],R);
                cblas_sscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                s = den / cblas_snrm2(C,&X[r*C],1);
                cblas_sscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stdev1_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int stdev1_d (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    double s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in stdev1_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stdev1_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrt(R); }
        else { den = sqrt(R-1); }

        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                s = den / cblas_dnrm2(R,&X[c*R],1);
                cblas_dscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                s = den / cblas_dnrm2(R,&X[c],C);
                cblas_dscal(R,s,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = sqrt(C); }
        else { den = sqrt(C-1); }
        
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                s = den / cblas_dnrm2(C,&X[r],R);
                cblas_dscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                s = den / cblas_dnrm2(C,&X[r*C],1);
                cblas_dscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stdev1_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int stdev1_c (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    float s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in stdev1_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stdev1_c: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrtf(R); }
        else { den = sqrtf(R-1); }

        if (iscolmajor)
        {
            for (c=0; c<2*C; c+=2)
            {
                s = den / cblas_scnrm2(R,&X[c*R],1);
                cblas_csscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<2*C; c+=2)
            {
                s = den / cblas_scnrm2(R,&X[c],C);
                cblas_csscal(R,s,&X[c],C); 
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = sqrtf(C); }
        else { den = sqrtf(C-1); }
        
        if (iscolmajor)
        {
            for (r=0; r<2*R; r+=2)
            {
                s = den / cblas_scnrm2(C,&X[r],R);
                cblas_csscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<2*R; r+=2)
            {
                s = den / cblas_scnrm2(C,&X[r*C],1);
                cblas_csscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stdev1_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int stdev1_z (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    double s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in stdev1_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stdev1_z: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrt(R); }
        else { den = sqrt(R-1); }

        if (iscolmajor)
        {
            for (c=0; c<2*C; c+=2)
            {
                s = den / cblas_dznrm2(R,&X[c*R],1);
                cblas_zdscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<2*C; c+=2)
            {
                s = den / cblas_dznrm2(R,&X[c],C);
                cblas_zdscal(R,s,&X[c],C); 
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = sqrt(C); }
        else { den = sqrt(C-1); }
        
        if (iscolmajor)
        {
            for (r=0; r<2*R; r+=2)
            {
                s = den / cblas_dznrm2(C,&X[r],R);
                cblas_zdscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<2*R; r+=2)
            {
                s = den / cblas_dznrm2(C,&X[r*C],1);
                cblas_zdscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stdev1_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
