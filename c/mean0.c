//Zeroes mean of each row or col according to dim.

#include <stdio.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int mean0_s (float *X, const char iscolmajor, const int R, const int C, const size_t dim);
int mean0_d (double *X, const char iscolmajor, const int R, const int C, const size_t dim);
int mean0_c (float *X, const char iscolmajor, const int R, const int C, const size_t dim);
int mean0_z (double *X, const char iscolmajor, const int R, const int C, const size_t dim);


int mean0_s (float *X, const char iscolmajor, const int R, const int C, const size_t dim)
{
    const float o = 1.0f;
    float m;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean0_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean0_s: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = cblas_sdot(R,&X[c*R],1,&o,0) / R;
                cblas_saxpy(R,-m,&o,0,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = cblas_sdot(R,&X[c],C,&o,0) / R;
                cblas_saxpy(R,-m,&o,0,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = cblas_sdot(C,&X[r],R,&o,0) / C;
                cblas_saxpy(C,-m,&o,0,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = cblas_sdot(C,&X[r*C],1,&o,0) / C;
                cblas_saxpy(C,-m,&o,0,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int mean0_d (double *X, const char iscolmajor, const int R, const int C, const size_t dim)
{
    const double o = 1.0;
    double m;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean0_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean0_d: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                m = cblas_ddot(R,&X[c*R],1,&o,0) / R;
                cblas_daxpy(R,-m,&o,0,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                m = cblas_ddot(R,&X[c],C,&o,0) / R;
                cblas_daxpy(R,-m,&o,0,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                m = cblas_ddot(C,&X[r],R,&o,0) / C;
                cblas_daxpy(C,-m,&o,0,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                m = cblas_ddot(C,&X[r*C],1,&o,0) / C;
                cblas_daxpy(C,-m,&o,0,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int mean0_c (float *X, const char iscolmajor, const int R, const int C, const size_t dim)
{
    //const float o = 1.0f; float m;
    const float o[2] =  {1.0f,0.0f};
    _Complex float m;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean0_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean0_c: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<2*C; c+=2)
            {
                //m = cblas_sdot(R,&X[c*R],2,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X[c*R],2);
                //m = cblas_sdot(R,&X[c*R+1],2,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X[c*R+1],2);

                //cblas_saxpy(R,-cblas_sdotu(R,&X[c*R],2,&o,0)/R,&o,0,&X[c*R],2);
                //cblas_saxpy(R,-cblas_sdotu(R,&X[c*R+1],2,&o,0)/R,&o,0,&X[c*R+1],2);

                m = -cblas_cdotu(R,&X[c*R],1,&o[0],0) / R;
                cblas_caxpy(R,(float *)&m,&o[0],0,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<2*C; c+=2)
            {
                //m = cblas_sdot(R,&X[2*c],2*C,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X[2*c],2*C);
                //m = cblas_sdot(R,&X[2*c+1],2*C,&o,0) / R;
                //cblas_saxpy(R,-m,&o,0,&X[2*c+1],2*C);
                
                m = -cblas_cdotu(R,&X[c],C,&o[0],0) / R;
                cblas_caxpy(R,(float *)&m,&o[0],0,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_cdotu(C,&X[r],R,&o[0],0) / C;
                cblas_caxpy(C,(float *)&m,&o[0],0,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_cdotu(C,&X[r*C],1,&o[0],0) / C;
                cblas_caxpy(C,(float *)&m,&o[0],0,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int mean0_z (double *X, const char iscolmajor, const int R, const int C, const size_t dim)
{
    const double o[2] =  {1.0,0.0};
    _Complex double m;
    int r, c;

    //Checks
    if (R<1) { fprintf(stderr,"error in mean0_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in mean0_z: C (ncols X) must be positive\n"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<2*C; c+=2)
            {
                m = -cblas_zdotu(R,&X[c*R],1,&o[0],0) / R;
                cblas_zaxpy(R,(double *)&m,&o[0],0,&X[c*R],1);
            }
        }
        else
        {
            for (c=0; c<2*C; c+=2)
            {
                m = -cblas_zdotu(R,&X[c],C,&o[0],0) / R;
                cblas_zaxpy(R,(double *)&m,&o[0],0,&X[c],C);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_zdotu(C,&X[r],R,&o[0],0) / C;
                cblas_zaxpy(C,(double *)&m,&o[0],0,&X[r],R);
            }
        }
        else
        {
            for (r=0; r<2*R; r+=2)
            {
                m = -cblas_zdotu(C,&X[r*C],1,&o[0],0) / C;
                cblas_zaxpy(C,(double *)&m,&o[0],0,&X[r*C],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in mean0_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
