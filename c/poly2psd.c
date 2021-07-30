//Gets power spectral densities (PSDs) from polynomial coeffs along rows or cols of X.
//The 2nd input is the vector V of variances (prediction errors) for each row or col of X.
//The 3rd input is a vector W of F freqs (in radians) at which to get the PSD.

//Following convention of Octave signal package ar_psd.m, I double the power for real-valued X.
//See Eq. (2.38) of Kay and Marple [1981].

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int poly2psd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const float *V, const float *W, const size_t F, const size_t dim);
int poly2psd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const double *V, const double *W, const size_t F, const size_t dim);
int poly2psd_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const float *V, const float *W, const size_t F, const size_t dim);
int poly2psd_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const double *V, const double *W, const size_t F, const size_t dim);


int poly2psd_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const float *V, const float *W, const size_t F, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2psd_s: dim must be in [0 3]\n"); return 1; }
    if (F<1u) { fprintf(stderr,"error in poly2psd_d: F (length W) must be positive\n"); return 1; }

    const size_t P = (dim==0u) ? R-1u : C-1u;
    int r, c, p, f, n;
    float *Er, *Ei, *yr, *yi;
    if (!(yr=(float *)malloc((size_t)(F)*sizeof(float)))) { fprintf(stderr,"error in poly2psd_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(yi=(float *)malloc((size_t)(F)*sizeof(float)))) { fprintf(stderr,"error in poly2psd_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Er=(float *)malloc((size_t)(F*P)*sizeof(float)))) { fprintf(stderr,"error in poly2psd_s: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Ei=(float *)malloc((size_t)(F*P)*sizeof(float)))) { fprintf(stderr,"error in poly2psd_s: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        Er[n] = cosf(W[f]*(p+1)); Ei[n] = -sinf(W[f]*(p+1));
    }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Er,F,&X[c*R+1],1,0.0f,yr,1);
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Ei,F,&X[c*R+1],1,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[c*F+f] = 2.0f*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Er,P,&X[c+C],(int)C,0.0f,yr,1);
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Ei,P,&X[c+C],(int)C,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[c+f*C] = 2.0f*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Er,F,&X[r+R],(int)R,0.0f,yr,1);
                cblas_sgemv(CblasColMajor,CblasNoTrans,F,P,1.0f,Ei,F,&X[r+R],(int)R,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[r+f*R] = 2.0f*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Er,P,&X[r*C+1],1,0.0f,yr,1);
                cblas_sgemv(CblasRowMajor,CblasNoTrans,F,P,1.0f,Ei,P,&X[r*C+1],1,0.0f,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0f; Y[r*F+f] = 2.0f*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2psd_s: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(Er); free(Ei); free(yr); free(yi);
    return 0;
}


int poly2psd_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const double *V, const double *W, const size_t F, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2psd_d: dim must be in [0 3]\n"); return 1; }
    if (F<1u) { fprintf(stderr,"error in poly2psd_d: F (length W) must be positive\n"); return 1; }

    const size_t P = (dim==0u) ? R-1u : C-1u;
    int r, c, p, f, n;
    double *Er, *Ei, *yr, *yi;

    if (!(yr=(double *)malloc((size_t)(F)*sizeof(double)))) { fprintf(stderr,"error in poly2psd_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(yi=(double *)malloc((size_t)(F)*sizeof(double)))) { fprintf(stderr,"error in poly2psd_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Er=(double *)malloc((size_t)(F*P)*sizeof(double)))) { fprintf(stderr,"error in poly2psd_d: problem with malloc. "); perror("malloc"); return 1; }
    if (!(Ei=(double *)malloc((size_t)(F*P)*sizeof(double)))) { fprintf(stderr,"error in poly2psd_d: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        Er[n] = cos(W[f]*(p+1)); Ei[n] = -sin(W[f]*(p+1));
    }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Er,F,&X[c*R+1],1,0.0,yr,1);
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Ei,F,&X[c*R+1],1,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[c*F+f] = 2.0*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Er,P,&X[c+C],(int)C,0.0,yr,1);
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Ei,P,&X[c+C],(int)C,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[c+f*C] = 2.0*V[c]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Er,F,&X[r+R],(int)R,0.0,yr,1);
                cblas_dgemv(CblasColMajor,CblasNoTrans,F,P,1.0,Ei,F,&X[r+R],(int)R,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[r+f*R] = 2.0*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Er,P,&X[r*C+1],1,0.0,yr,1);
                cblas_dgemv(CblasRowMajor,CblasNoTrans,F,P,1.0,Ei,P,&X[r*C+1],1,0.0,yi,1);
                for (f=0; f<F; f++) { yr[f] += 1.0; Y[r*F+f] = 2.0*V[r]/(yr[f]*yr[f]+yi[f]*yi[f]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2psd_d: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(Er); free(Ei); free(yr); free(yi);
    return 0;
}


int poly2psd_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const float *V, const float *W, const size_t F, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2psd_c: dim must be in [0 3]\n"); return 1; }
    if (F<1u) { fprintf(stderr,"error in poly2psd_c: F (length W) must be positive\n"); return 1; }

    const float z[2] = {0.0f,0.0f}, o[2] = {1.0f,0.0f};
    const size_t P = (dim==0u) ? R-1u : C-1u;
    int r, c, p, f, n;
    float *E, *y;
    if (!(y=(float *)malloc((size_t)(2*F)*sizeof(float)))) { fprintf(stderr,"error in poly2psd_c: problem with malloc. "); perror("malloc"); return 1; }
    if (!(E=(float *)malloc((size_t)(2*F*P)*sizeof(float)))) { fprintf(stderr,"error in poly2psd_c: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        E[2*n] = cosf(W[f]*(p+1)); E[2*n+1] = -sinf(W[f]*(p+1));
    }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_cgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*c*R+2],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[c*F+f] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_cgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*(c+C)],(int)C,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[c+f*C] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_cgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*(r+R)],(int)R,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[r+f*R] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_cgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*r*C+2],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0f; Y[r*F+f] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2psd_c: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(E); free(y);
    return 0;
}


int poly2psd_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const double *V, const double *W, const size_t F, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2psd_z: dim must be in [0 3]\n"); return 1; }
    if (F<1u) { fprintf(stderr,"error in poly2psd_z: F (length W) must be positive\n"); return 1; }

    const double z[2] = {0.0,0.0}, o[2] = {1.0,0.0};
    const size_t P = (dim==0u) ? R-1u : C-1u;
    int r, c, p, f, n;
    double *E, *y;
    if (!(y=(double *)malloc((size_t)(2*F)*sizeof(double)))) { fprintf(stderr,"error in poly2psd_z: problem with malloc. "); perror("malloc"); return 1; }
    if (!(E=(double *)malloc((size_t)(2*F*P)*sizeof(double)))) { fprintf(stderr,"error in poly2psd_z: problem with malloc. "); perror("malloc"); return 1; }

    //Make complex-valued E matrix
    for (n=0; n<F*P; n++)
    {
        if (iscolmajor) { f = n%F; p = n/F; } else { f = n/P; p = n%P; }
        E[2*n] = cos(W[f]*(p+1)); E[2*n+1] = -sin(W[f]*(p+1));
    }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_zgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*c*R+2],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[c*F+f] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (size_t c=0u; c<C; ++c)
            {
                cblas_zgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*(c+C)],(int)C,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[c+f*C] = V[c] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_zgemv(CblasColMajor,CblasNoTrans,F,P,o,E,F,&X[2*(r+R)],(int)R,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[r+f*R] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
        else
        {
            for (size_t r=0u; r<R; ++r)
            {
                cblas_zgemv(CblasRowMajor,CblasNoTrans,F,P,o,E,P,&X[2*r*C+2],1,z,y,1);
                for (f=0; f<F; f++) { y[2*f] += 1.0; Y[r*F+f] = V[r] / (y[2*f]*y[2*f]+y[2*f+1]*y[2*f+1]); }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in poly2psd_z: dim must be 0 or 1.\n"); return 1;
    }

    //Exit
    free(E); free(y);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
