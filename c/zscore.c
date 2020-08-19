//Z-scores each row or col according to dim.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zscore_s (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);
int zscore_d (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);
int zscore_c (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);
int zscore_z (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased);


int zscore_s (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    const float o = 1.0f;
    float m, s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in zscore_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zscore_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrtf(R-1); }
        else { den = sqrtf(R); }

        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = cblas_sdot(R,&X[c*R],1,&o,0) / R;
                cblas_saxpy(R,-m,&o,0,&X[c*R],1);
                s = den / cblas_snrm2(R,&X[c*R],1);
                cblas_sscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = cblas_sdot(R,&X[c],C,&o,0) / R;
                cblas_saxpy(R,-m,&o,0,&X[c],C);
                s = den / cblas_snrm2(R,&X[c],C);
                cblas_sscal(R,s,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = sqrtf(C-1); }
        else { den = sqrtf(C); }
        
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = cblas_sdot(C,&X[r],R,&o,0) / C;
                cblas_saxpy(C,-m,&o,0,&X[r],R);
                s = den / cblas_snrm2(C,&X[r],R);
                cblas_sscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = cblas_sdot(C,&X[r*C],1,&o,0) / C;
                cblas_saxpy(C,-m,&o,0,&X[r*C],1);
                s = den / cblas_snrm2(C,&X[r*C],1);
                cblas_sscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in zscore_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int zscore_d (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    const double o = 1.0;
    double m, s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in zscore_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zscore_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrt(R); }
        else { den = sqrt(R-1); }

        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = cblas_ddot(R,&X[c*R],1,&o,0) / R;
                cblas_daxpy(R,-m,&o,0,&X[c*R],1);
                s = den / cblas_dnrm2(R,&X[c*R],1);
                cblas_dscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = cblas_ddot(R,&X[c],C,&o,0) / R;
                cblas_daxpy(R,-m,&o,0,&X[c],C);
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
                m = cblas_ddot(C,&X[r],R,&o,0) / C;
                cblas_daxpy(C,-m,&o,0,&X[r],R);
                s = den / cblas_dnrm2(C,&X[r],R);
                cblas_dscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = cblas_ddot(C,&X[r*C],1,&o,0) / C;
                cblas_daxpy(C,-m,&o,0,&X[r*C],1);
                s = den / cblas_dnrm2(C,&X[r*C],1);
                cblas_dscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in zscore_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int zscore_c (float *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    const float o[2] =  {1.0,0.0};
    _Complex float m;
    float s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in zscore_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zscore_c: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrtf(R-1); }
        else { den = sqrtf(R); }

        if (iscolmajor)
        {
            for (c=0; c<2*C; c+=2)
            {
                m = -cblas_cdotu(R,&X[c*R],1,&o[0],0) / R;
                cblas_caxpy(R,(float *)&m,&o[0],0,&X[c*R],1);
                s = den / cblas_scnrm2(R,&X[c*R],1);
                cblas_csscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<2*C; c+=2)
            {
                m = -cblas_cdotu(R,&X[c],C,&o[0],0) / R;
                cblas_caxpy(R,(float *)&m,&o[0],0,&X[c],C);
                s = den / cblas_scnrm2(R,&X[c],C);
                cblas_csscal(R,s,&X[c],C); 
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = sqrtf(C-1); }
        else { den = sqrtf(C); }
        
        if (iscolmajor)
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_cdotu(C,&X[r],R,&o[0],0) / C;
                cblas_caxpy(C,(float *)&m,&o[0],0,&X[r],R);
                s = den / cblas_scnrm2(C,&X[r],R);
                cblas_csscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_cdotu(C,&X[r*C],1,&o[0],0) / C;
                cblas_caxpy(C,(float *)&m,&o[0],0,&X[r*C],1);
                s = den / cblas_scnrm2(C,&X[r*C],1);
                cblas_csscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in zscore_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int zscore_z (double *X, const char iscolmajor, const int R, const int C, const size_t dim, const char biased)
{
    const double o[2] =  {1.0,0.0};
    _Complex double m;
    double s, den;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in zscore_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in zscore_z: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (biased) { den = sqrt(R-1); }
        else { den = sqrt(R); }

        if (iscolmajor)
        {
            for (c=0; c<2*C; c+=2)
            {
                m = -cblas_zdotu(R,&X[c*R],1,&o[0],0) / R;
                cblas_zaxpy(R,(double *)&m,&o[0],0,&X[c*R],1);
                s = den / cblas_dznrm2(R,&X[c*R],1);
                cblas_zdscal(R,s,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<2*C; c+=2)
            {
                m = -cblas_zdotu(R,&X[c],C,&o[0],0) / R;
                cblas_zaxpy(R,(double *)&m,&o[0],0,&X[c],C);
                s = den / cblas_dznrm2(R,&X[c],C);
                cblas_zdscal(R,s,&X[c],C); 
            }
        }
    }
    else if (dim==1)
    {
        if (biased) { den = sqrt(C-1); }
        else { den = sqrt(C); }
        
        if (iscolmajor)
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_zdotu(C,&X[r],R,&o[0],0) / C;
                cblas_zaxpy(C,(double *)&m,&o[0],0,&X[r],R);
                s = den / cblas_dznrm2(C,&X[r],R);
                cblas_zdscal(C,s,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_zdotu(C,&X[r*C],1,&o[0],0) / C;
                cblas_zaxpy(C,(double *)&m,&o[0],0,&X[r*C],1);
                s = den / cblas_dznrm2(C,&X[r*C],1);
                cblas_zdscal(C,s,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in zscore_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
