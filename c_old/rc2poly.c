//Gets polynomial coefficients from reflection coefficients (RCs) along rows or cols of X.
//Input (X) and output (Y) have the same size.

#include <stdio.h>
#include <stdlib.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int rc2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rc2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rc2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int rc2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int rc2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2poly_s: dim must be in [0 3]\n"); return 1; }

    const float o = 1.0f;
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    float sc;
    float *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2poly_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2poly_s: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in rc2poly_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(C,&o,0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&Y[c*(R+1)+1],1);
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[c*R+1],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*R+q+1] -= sc * y[p-1-q]; }
                }
                cblas_sscal(R,-1.0f,&Y[c*(R+1)+1],1);
            }
        }
        else
        {
            cblas_scopy(C,&o,0,Y,1);
            cblas_scopy(R*C,X,1,&Y[C],1);
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[C+c],C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[C+c+q*C] -= sc * y[p-1-q]; }
                }
            }
            cblas_sscal(R*C,-1.0f,&Y[C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R,&o,0,Y,1);
            cblas_scopy(R*C,X,1,&Y[R],1);
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[R+r],R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[R+r+q*R] -= sc * y[p-1-q]; }
                }
            }
            cblas_sscal(R*C,-1.0f,&Y[R],1);
        }
        else
        {
            cblas_scopy(R,&o,0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&Y[r*(C+1)+1],1);
                for (p=1; p<P; p++)
                {
                    cblas_scopy(p+1,&Y[r*C+1],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*C+q+1] -= sc * y[p-1-q]; }
                }
                cblas_sscal(C,-1.0f,&Y[r*(C+1)+1],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2poly_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


int rc2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2poly_d: dim must be in [0 3]\n"); return 1; }

    const double o = 1.0;
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    double sc;
    double *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2poly_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2poly_d: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in rc2poly_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(C,&o,0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&Y[c*(R+1)+1],1);
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[c*R+1],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*R+q+1] -= sc * y[p-1-q]; }
                }
                cblas_dscal(R,-1.0,&Y[c*(R+1)+1],1);
            }
        }
        else
        {
            cblas_dcopy(C,&o,0,Y,1);
            cblas_dcopy(R*C,X,1,&Y[C],1);
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[C+c],C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[C+c+q*C] -= sc * y[p-1-q]; }
                }
            }
            cblas_dscal(R*C,-1.0,&Y[C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R,&o,0,Y,1);
            cblas_dcopy(R*C,X,1,&Y[R],1);
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[R+r],R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[R+r+q*R] -= sc * y[p-1-q]; }
                }
            }
            cblas_dscal(R*C,-1.0,&Y[R],1);
        }
        else
        {
            cblas_dcopy(R,&o,0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&Y[r*(C+1)+1],1);
                for (p=1; p<P; p++)
                {
                    cblas_dcopy(p+1,&Y[r*C+1],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*C+q+1] -= sc * y[p-1-q]; }
                }
                cblas_dscal(C,-1.0,&Y[r*(C+1)+1],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2poly_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


int rc2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2poly_c: dim must be in [0 3]\n"); return 1; }

    const float o[2] =  {1.0f,0.0f};
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    float sc[2] = {0.0f,0.0f};
    float *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2poly_c: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2poly_c: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(float *)malloc((size_t)(2*P)*sizeof(float)))) { fprintf(stderr,"error in rc2poly_c: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_ccopy(C,&o[0],0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c*R],1,&Y[2*c*(R+1)+2],1);
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*c*R+2],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q+1)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q+1)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
                cblas_csscal(R,-1.0f,&Y[2*c*(R+1)+2],1);
            }
        }
        else
        {
            cblas_ccopy(C,&o[0],0,Y,1);
            cblas_ccopy(R*C,X,1,&Y[2*C],1);
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c*R],1,&Y[2*c*(R+1)+2],1);
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*(C+c)],C,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(C+c+q*C)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(C+c+q*C)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
            cblas_csscal(R*C,-1.0f,&Y[2*C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_ccopy(R,&o[0],0,Y,1);
            cblas_ccopy(R*C,X,1,&Y[2*R],1);
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*(r+R)],R,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(R+r+q*R)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(R+r+q*R)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
            cblas_csscal(R*C,-1.0f,&Y[2*R],1);
        }
        else
        {
            cblas_ccopy(R,&o[0],0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r*C],1,&Y[2*r*(C+1)+2],1);
                for (p=1; p<P; p++)
                {
                    cblas_ccopy(p+1,&Y[2*r*C+2],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q+1)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q+1)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
                cblas_csscal(C,-1.0f,&Y[2*r*(C+1)+2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2poly_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


int rc2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in rc2poly_z: dim must be in [0 3]\n"); return 1; }

    const double o[2] =  {1.0,0.0};
    const int P = (dim==0) ? R : C;
    int r, c, p, q;
    double sc[2] = {0.0,0.0};
    double *y;

    //Checks
    if (R<1) { fprintf(stderr,"error in rc2poly_z: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in rc2poly_z: ncols X must be positive\n"); return 1; }

    //Initialize
    if (!(y=(double *)malloc((size_t)(2*P)*sizeof(double)))) { fprintf(stderr,"error in rc2poly_z: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_zcopy(C,&o[0],0,Y,R+1);
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c*R],1,&Y[2*c*(R+1)+2],1);
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*c*R+2],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q+1)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q+1)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
                cblas_zdscal(R,-1.0,&Y[2*c*(R+1)+2],1);
            }
        }
        else
        {
            cblas_zcopy(C,&o[0],0,Y,1);
            cblas_zcopy(R*C,X,1,&Y[2*C],1);
            for (c=0; c<C; c++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*(C+c)],C,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(C+c+q*C)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(C+c+q*C)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
            cblas_zdscal(R*C,-1.0,&Y[2*C],1);
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_zcopy(R,&o[0],0,Y,1);
            cblas_zcopy(R*C,X,1,&Y[2*R],1);
            for (r=0; r<R; r++)
            {
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*(R+r)],R,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(R+r+q*R)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(R+r+q*R)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
            }
            cblas_zdscal(R*C,-1.0,&Y[2*R],1);
        }
        else
        {
            cblas_zcopy(R,&o[0],0,Y,C+1);
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r*C],1,&Y[2*r*(C+1)+2],1);
                for (p=1; p<P; p++)
                {
                    cblas_zcopy(p+1,&Y[2*r*C+2],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q+1)] -= sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q+1)+1] -= sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                }
                cblas_zdscal(C,-1.0,&Y[2*r*(C+1)+2],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in rc2poly_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(y);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
