//Gets reflection coefficients (RCs) from polynomials along rows or cols of X.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int poly2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int poly2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2rc_s: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in poly2rc_s: P (num rcs) must be positive\n"); return 1; }

    const size_t P = (dim==0u) ? R-1u : C-1u;
    int r, c, p, q;
    float sc;
    float *y;
    if (!(y=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in poly2rc_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_scopy(P,&X[c*R+1],1,&Y[c*P],1);
                cblas_sscal(P,-1.0f/X[c*R],&Y[c*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[c*P],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*P+q] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[c*P],1);
                }
            }
        }
        else
        {
            cblas_scopy(P*C,&X[C],1,Y,1);
            for (size_t c=0u; c<C; ++c)
            {
                cblas_sscal(P,-1.0f/X[c],&Y[c],(int)C);
                for (p=P-1u; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[c],(int)C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c+q*C] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[c],(int)C);
                }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            cblas_scopy(R*P,&X[R],1,Y,1);
            for (size_t r=0u; r<R; ++r)
            {
                cblas_sscal(P,-1.0f/X[r],&Y[r],(int)R);
                for (p=P-1u; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[r],(int)R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r+q*R] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[r],(int)R);
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_scopy(P,&X[r*C+1],1,&Y[r*P],1);
                cblas_sscal(P,-1.0f/X[r*C],&Y[r*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_scopy(p+1,&Y[r*P],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*P+q] += sc * y[p-1-q]; }
                    cblas_sscal(p,1.0f/fmaf(sc,-sc,1.0f),&Y[r*P],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2rc_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int poly2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2rc_d: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in poly2rc_d: P (num rcs) must be positive\n"); return 1; }

    const size_t P = (dim==0u) ? R : C;
    int r, c, p, q;
    double sc;
    double *y;
    if (!(y=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in poly2rc_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_dcopy(P,&X[c*R+1],1,&Y[c*P],1);
                cblas_dscal(P,-1.0/X[c*R],&Y[c*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[c*P],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c*P+q] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[c*P],1);
                }
            }
        }
        else
        {
            cblas_dcopy(P*C,&X[C],1,Y,1);
            for (size_t c=0u; c<C; ++c)
            {
                cblas_dscal(P,-1.0/X[c],&Y[c],(int)C);
                for (p=P-1u; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[c],(int)C,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[c+q*C] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[c],(int)C);
                }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            cblas_dcopy(R*P,&X[R],1,Y,1);
            for (size_t r=0u; r<R; ++r)
            {
                cblas_dscal(P,-1.0/X[r],&Y[r],(int)R);
                for (p=P-1u; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[r],(int)R,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r+q*R] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[r],(int)R);
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_dcopy(P,&X[r*C+1],1,&Y[r*P],1);
                cblas_dscal(P,-1.0/X[r*C],&Y[r*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_dcopy(p+1,&Y[r*P],1,y,1);
                    sc = y[p];
                    for (q=0; q<p; q++) { Y[r*P+q] += sc * y[p-1-q]; }
                    cblas_dscal(p,1.0/fma(sc,-sc,1.0),&Y[r*P],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2rc_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int poly2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2rc_c: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in poly2rc_c: P (num rcs) must be positive\n"); return 1; }

    const size_t P = (dim==0u) ? R-1u : C-1u;
    int r, c, p, q;
    float sc[2] = {0.0f,0.0f};
    float *y;
    if (!(y=(float *)malloc((size_t)(2*P)*sizeof(float)))) { fprintf(stderr,"error in poly2rc_c: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_ccopy(P,&X[2*c*R+2],1,&Y[2*c*P],1);
                sc[0] = 1.0f / (X[2*c*R]*X[2*c*R]+X[2*c*R+1]*X[2*c*R+1]);
                sc[1] = sc[0] * X[2*c*R+1]; sc[0] *= -X[2*c*R];
                cblas_cscal(P,&sc[0],&Y[2*c*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*c*R],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c*R],1);
                }
            }
        }
        else
        {
            cblas_ccopy(P*C,&X[2*C],1,Y,1);
            for (size_t c=0u; c<C; ++c)
            {
                sc[0] = 1.0f / (X[2*c]*X[2*c]+X[2*c+1]*X[2*c+1]);
                sc[1] = sc[0] * X[2*c+1]; sc[0] *= -X[2*c];
                cblas_cscal(P,&sc[0],&Y[2*c],2*(int)C);
                for (p=P-1u; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*c],(int)C,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c+q*C)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c+q*C)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c],(int)C);
                }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            cblas_ccopy(R*P,&X[2*R],1,Y,1);
            for (size_t r=0u; r<R; ++r)
            {
                sc[0] = 1.0f / (X[2*r]*X[2*r]+X[2*r+1]*X[2*r+1]);
                sc[1] = sc[0] * X[2*r+1]; sc[0] *= -X[2*r];
                cblas_cscal(P,&sc[0],&Y[2*r],(int)R);
                for (p=P-1u; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*r],(int)R,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r+q*R)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r+q*R)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r],(int)R);
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_ccopy(P,&X[2*r*C+2],1,&Y[2*r*P],1);
                sc[0] = 1.0f / (X[2*r*C]*X[2*r*C]+X[2*r*C+1]*X[2*r*C+1]);
                sc[1] = sc[0] * X[2*r*C+1]; sc[0] *= -X[2*r*C];
                cblas_cscal(P,&sc[0],&Y[2*r*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_ccopy(p+1,&Y[2*r*C],1,y,1);
                    cblas_ccopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_csscal(p,1.0f/(1.0f-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r*C],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2rc_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


int poly2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2rc_z: dim must be in [0 3]\n"); return 1; }
    if (P<1u) { fprintf(stderr,"error in poly2rc_z: P (num rcs) must be positive\n"); return 1; }
    
    const size_t P = (dim==0u) ? R-1u : C-1u;
    int r, c, p, q;
    double sc[2] = {0.0,0.0};
    double *y;
    if (!(y=(double *)malloc((size_t)(2*P)*sizeof(double)))) { fprintf(stderr,"error in poly2rc_z: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_zcopy(P,&X[2*c*R+2],1,&Y[2*c*P],1);
                sc[0] = 1.0 / (X[2*c*R]*X[2*c*R]+X[2*c*R+1]*X[2*c*R+1]);
                sc[1] = sc[0] * X[2*c*R+1]; sc[0] *= -X[2*c*R];
                cblas_zscal(P,&sc[0],&Y[2*c*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*c*R],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c*R+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c*R+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c*R],1);
                }
            }
        }
        else
        {
            cblas_zcopy(P*C,&X[2*C],1,Y,1);
            for (size_t c=0u; c<C; ++c)
            {
                sc[0] = 1.0 / (X[2*c]*X[2*c]+X[2*c+1]*X[2*c+1]);
                sc[1] = sc[0] * X[2*c+1]; sc[0] *= -X[2*c];
                cblas_zscal(P,&sc[0],&Y[2*c],2*(int)C);
                for (p=P-1u; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*c],(int)C,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(c+q*C)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(c+q*C)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*c],(int)C);
                }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            cblas_zcopy(R*P,&X[2*R],1,Y,1);
            for (size_t r=0u; r<R; ++r)
            {
                sc[0] = 1.0 / (X[2*r]*X[2*r]+X[2*r+1]*X[2*r+1]);
                sc[1] = sc[0] * X[2*r+1]; sc[0] *= -X[2*r];
                cblas_zscal(P,&sc[0],&Y[2*r],(int)R);
                for (p=P-1u; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*r],(int)R,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r+q*R)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r+q*R)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r],(int)R);
                }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_zcopy(P,&X[2*r*C+2],1,&Y[2*r*P],1);
                sc[0] = 1.0 / (X[2*r*C]*X[2*r*C]+X[2*r*C+1]*X[2*r*C+1]);
                sc[1] = sc[0] * X[2*r*C+1]; sc[0] *= -X[2*r*C];
                cblas_zscal(P,&sc[0],&Y[2*r*P],1);
                for (p=P-1u; p>0; p--)
                {
                    cblas_zcopy(p+1,&Y[2*r*C],1,y,1);
                    cblas_zcopy(1,&y[2*p],1,&sc[0],1);
                    for (q=0; q<p; q++)
                    {
                        Y[2*(r*C+q)] += sc[0]*y[2*(p-1-q)] - sc[1]*y[2*(p-1-q)+1];
                        Y[2*(r*C+q)+1] += sc[0]*y[2*(p-1-q)+1] + sc[1]*y[2*(p-1-q)];
                    }
                    cblas_zdscal(p,1.0/(1.0-(sc[0]*sc[0]+sc[1]*sc[1])),&Y[2*r*C],1);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2rc_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    return 0;
}


#ifdef __cplusplus
}
}
#endif
