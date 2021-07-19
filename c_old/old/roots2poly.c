//Polynomial from roots along cols or rows of X.

#include <stdio.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int roots2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int roots2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int roots2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int roots2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int roots2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in roots2poly_s: dim must be in [0 3]\n"); return 1; }

    const float z = 0.0f, o = 1.0f;
    const size_t P = (dim==0) ? R+1 : C+1;

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy(P*C,&z,0,Y,1); cblas_scopy(C,&o,0,Y,P);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--) { Y[p+c*P] = fmaf(-X[r+c*R],Y[p-1+c*P],Y[p+c*P]); }
                }
            }
        }
        else
        {
            cblas_scopy(P*C,&z,0,Y,1); cblas_scopy(C,&o,0,Y,1);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--) { Y[c+p*C] = fmaf(-X[c+r*C],Y[c+(p-1)*C],Y[c+p*C]); }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy(R*P,&z,0,Y,1); cblas_scopy((int)R,&o,0,Y,1);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--) { Y[r+p*R] = fmaf(-X[r+c*R],Y[r+(p-1)*R],Y[r+p*R]); }
                }
            }
        }
        else
        {
            cblas_scopy(R*P,&z,0,Y,1); cblas_scopy((int)R,&o,0,Y,P);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--) { Y[p+r*P] = fmaf(-X[c+r*C],Y[p-1+r*P],Y[p+r*P]); }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in roots2poly_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int roots2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in roots2poly_d: dim must be in [0 3]\n"); return 1; }

    const double z = 0.0, o = 1.0;
    const size_t P = (dim==0) ? R+1 : C+1;

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy(P*C,&z,0,Y,1); cblas_dcopy(C,&o,0,Y,P);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--) { Y[p+c*P] = fma(-X[r+c*R],Y[p-1+c*P],Y[p+c*P]); }
                }
            }
        }
        else
        {
            cblas_dcopy(P*C,&z,0,Y,1); cblas_dcopy(C,&o,0,Y,1);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--) { Y[c+p*C] = fma(-X[c+r*C],Y[c+(p-1)*C],Y[c+p*C]); }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R*P,&z,0,Y,1); cblas_dcopy((int)R,&o,0,Y,1);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--) { Y[r+p*R] = fma(-X[r+c*R],Y[r+(p-1)*R],Y[r+p*R]); }
                }
            }
        }
        else
        {
            cblas_dcopy(R*P,&z,0,Y,1); cblas_dcopy((int)R,&o,0,Y,P);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--) { Y[p+r*P] = fma(-X[c+r*C],Y[p-1+r*P],Y[p+r*P]); }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in roots2poly_d: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int roots2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in roots2poly_c: dim must be in [0 3]\n"); return 1; }

    const float z[2] =  {0.0,0.0}, o[2] =  {1.0,0.0};
    const size_t P = (dim==0) ? R+1 : C+1;

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_ccopy(P*C,z,0,Y,1); cblas_ccopy(C,o,0,Y,P);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        Y[2*(p+c*P)] = fmaf(-X[2*(r+c*R)],Y[2*(p-1+c*P)],fmaf(X[2*(r+c*R)+1],Y[2*(p-1+c*P)+1],Y[2*(p+c*P)]));
                        Y[2*(p+c*P)+1] = fmaf(-X[2*(r+c*R)],Y[2*(p-1+c*P)+1],fmaf(-X[2*(r+c*R)+1],Y[2*(p-1+c*P)],Y[2*(p+c*P)+1]));
                    }
                }
            }
        }
        else
        {
            cblas_ccopy(P*C,z,0,Y,1); cblas_ccopy(C,o,0,Y,1);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        Y[2*(c+p*C)] = fmaf(-X[2*(c+r*C)],Y[2*(c+(p-1)*C)],fmaf(X[2*(c+r*C)+1],Y[2*(c+(p-1)*C)+1],Y[2*(c+p*C)]));
                        Y[2*(c+p*C)+1] = fmaf(-X[2*(c+r*C)],Y[2*(c+(p-1)*C)+1],fmaf(-X[2*(c+r*C)+1],Y[2*(c+(p-1)*C)],Y[2*(c+p*C)+1]));
                    }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_ccopy(R*P,z,0,Y,1); cblas_ccopy((int)R,o,0,Y,1);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        Y[2*(r+p*R)] = fmaf(-X[2*(r+c*R)],Y[2*(r+(p-1)*R)],fmaf(X[2*(R+c*R)+1],Y[2*(r+(p-1)*R)+1],Y[2*(r+p*R)]));
                        Y[2*(r+p*R)+1] = fmaf(-X[2*(r+c*R)],Y[2*(r+(p-1)*R)+1],fmaf(-X[2*(r+c*R)+1],Y[2*(r+(p-1)*R)],Y[2*(r+p*R)+1]));
                    }
                }
            }
        }
        else
        {
            cblas_ccopy(R*P,z,0,Y,1); cblas_ccopy((int)R,o,0,Y,P);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        Y[2*(p+r*P)] = fmaf(-X[2*(c+r*C)],Y[2*(p-1+r*P)],fmaf(X[2*(c+r*C)+1],Y[2*(p-1+r*P)+1],Y[2*(p+r*P)]));
                        Y[2*(p+r*P)+1] = fmaf(-X[2*(c+r*C)],Y[2*(p-1+r*P)+1],fmaf(-X[2*(c+r*C)+1],Y[2*(p-1+r*P)],Y[2*(p+r*P)+1]));
                    }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in roots2poly_c: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int roots2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3) { fprintf(stderr,"error in roots2poly_z: dim must be in [0 3]\n"); return 1; }

    const double z[2] =  {0.0,0.0}, o[2] =  {1.0,0.0};
    const size_t P = (dim==0) ? R+1 : C+1;

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_zcopy(P*C,z,0,Y,1); cblas_zcopy(C,o,0,Y,P);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        Y[2*(p+c*P)] = fma(-X[2*(r+c*R)],Y[2*(p-1+c*P)],fma(X[2*(r+c*R)+1],Y[2*(p-1+c*P)+1],Y[2*(p+c*P)]));
                        Y[2*(p+c*P)+1] = fma(-X[2*(r+c*R)],Y[2*(p-1+c*P)+1],fma(-X[2*(r+c*R)+1],Y[2*(p-1+c*P)],Y[2*(p+c*P)+1]));
                    }
                }
            }
        }
        else
        {
            cblas_zcopy(P*C,z,0,Y,1); cblas_zcopy(C,o,0,Y,1);
            for (size_t c=0; c<C; c++)
            {
                for (size_t r=0; r<R; r++)
                {
                    for (p=r+1; p>0; p--)
                    {
                        Y[2*(c+p*C)] = fma(-X[2*(c+r*C)],Y[2*(c+(p-1)*C)],fma(X[2*(c+r*C)+1],Y[2*(c+(p-1)*C)+1],Y[2*(c+p*C)]));
                        Y[2*(c+p*C)+1] = fma(-X[2*(c+r*C)],Y[2*(c+(p-1)*C)+1],fma(-X[2*(c+r*C)+1],Y[2*(c+(p-1)*C)],Y[2*(c+p*C)+1]));
                    }
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_zcopy(R*P,z,0,Y,1); cblas_zcopy((int)R,o,0,Y,1);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        Y[2*(r+p*R)] = fma(-X[2*(r+c*R)],Y[2*(r+(p-1)*R)],fma(X[2*(R+c*R)+1],Y[2*(r+(p-1)*R)+1],Y[2*(r+p*R)]));
                        Y[2*(r+p*R)+1] = fma(-X[2*(r+c*R)],Y[2*(r+(p-1)*R)+1],fma(-X[2*(r+c*R)+1],Y[2*(r+(p-1)*R)],Y[2*(r+p*R)+1]));
                    }
                }
            }
        }
        else
        {
            cblas_zcopy(R*P,z,0,Y,1); cblas_zcopy((int)R,o,0,Y,P);
            for (size_t r=0; r<R; r++)
            {
                for (size_t c=0; c<C; c++)
                {
                    for (p=c+1; p>0; p--)
                    {
                        Y[2*(p+r*P)] = fma(-X[2*(c+r*C)],Y[2*(p-1+r*P)],fma(X[2*(c+r*C)+1],Y[2*(p-1+r*P)+1],Y[2*(p+r*P)]));
                        Y[2*(p+r*P)+1] = fma(-X[2*(c+r*C)],Y[2*(p-1+r*P)+1],fma(-X[2*(c+r*C)+1],Y[2*(p-1+r*P)],Y[2*(p+r*P)+1]));
                    }
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in roots2poly_z: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
