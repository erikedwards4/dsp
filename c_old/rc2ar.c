//Gets autoregressive (AR) parameters from reflection coefficients (RCs) along rows or cols of X.
//Input (X) and output (Y) have the same size.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rc2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rc2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rc2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rc2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int rc2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2ar_s: dim must be in [0 3]\n"); return 1; }

    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    float sc;
    float *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2ar_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2ar_s: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in rc2ar_s: problem with malloc. "); perror("malloc"); return 1; }
    cblas_scopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[c*R],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*R+q] -= sc * y[p-1-q]; }
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[c],C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c+q*C] -= sc * y[p-1-q]; }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[r],R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r+q*R] -= sc * y[p-1-q]; }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[r*C],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*C+q] -= sc * y[p-1-q]; }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2ar_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


int rc2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2ar_d: dim must be in [0 3]\n"); return 1; }

    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    double sc;
    double *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2ar_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2ar_d: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in rc2ar_d: problem with malloc. "); perror("malloc"); return 1; }
    cblas_dcopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[c*R],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*R+q] -= sc * y[p-1-q]; }
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[c],C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c+q*C] -= sc * y[p-1-q]; }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[r],R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r+q*R] -= sc * y[p-1-q]; }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[r*C],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*C+q] -= sc * y[p-1-q]; }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2ar_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


int rc2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2ar_c: dim must be in [0 3]\n"); return 1; }

    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    float sc[2] = {0.0f,0.0f};
    float *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2ar_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2ar_c: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(float *)malloc((size_t)(2*P)*sizeof(float)))) { fprintf(stderr,"error in rc2ar_c: problem with malloc. "); perror("malloc"); return 1; }
    cblas_ccopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*c*R],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*c],C,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c+q*C)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c+q*C)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*r],R,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r+q*R)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r+q*R)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*r*C],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2ar_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


int rc2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2ar_z: dim must be in [0 3]\n"); return 1; }

    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    double sc[2] = {0.0,0.0};
    double *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2ar_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2ar_z: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(double *)malloc((size_t)(2*P)*sizeof(double)))) { fprintf(stderr,"error in rc2ar_z: problem with malloc. "); perror("malloc"); return 1; }
    cblas_zcopy(R*C,X,1,Y,1);

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*c*R],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*c],C,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c+q*C)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c+q*C)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*r],R,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r+q*R)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r+q*R)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*r*C],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2ar_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
