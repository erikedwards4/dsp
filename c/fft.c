//This takes the output of window_univar, and applies the 1D FFT to each row or col.
//The input dim is which dim to take the 1D FFT along (0->along cols, 1->along rows).
//The output Y is complex-valued.
//For real-valued X, the size of Y must be FxC if dim==0, and RxF if dim==1, where F = nfft/2 + 1.
//For complex-valued X, the size of Y must be nfftxC if dim==0, and Rxnfft if dim==1.

//I tried parallel versions with OpenMP, but much slower (have to make fftw_plan P times!).

#include <stdio.h>
#include <cblas.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fft_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft);
int fft_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft);
int fft_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft);
int fft_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft);


int fft_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft)
{
    const float z = 0.0f;
    const int F = nfft/2 + 1;
    int r, c;
    float *X1;
    fftwf_complex *Y1;
    fftwf_plan plan;

    //Checks
    if (R<1) { fprintf(stderr,"error in fft_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fft_s: C (ncols X) must be positive\n"); return 1; }
    if (nfft<R && dim==0) { fprintf(stderr,"error in fft_s: nfft must be >= R (winlength)\n"); return 1; }
    if (nfft<C && dim==1) { fprintf(stderr,"error in fft_s: nfft must be >= C (winlength)\n"); return 1; }

    //Initialize fftwf
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_complex((size_t)F);
    plan = fftwf_plan_dft_r2c_1d(nfft,X1,Y1,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_s: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                fftwf_execute(plan);
                //cblas_scopy(2*F,&Y[2*c*F],1,(float *)&Y1[0],1);
                cblas_ccopy(F,(float *)&Y1[0],1,&Y[2*c*F],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,&X1[0],1);
                fftwf_execute(plan);
                //cblas_scopy(F,(float *)&Y1[0],2,&Y[2*c],2*C); cblas_scopy(F,(float *)&Y1[1],2,&Y[2*c+1],2*C);
                cblas_ccopy(F,(float *)&Y1[0],1,&Y[2*c],C);
            }
        }
    }
    else if (dim==1)
    {
        cblas_scopy(nfft-C,&z,0,&X1[C],1); //zero-pad
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r],R,&X1[0],1);
                fftwf_execute(plan);
                cblas_ccopy(F,(float *)&Y1[0],1,&Y[2*r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                fftwf_execute(plan);
                cblas_ccopy(F,(float *)&Y1[0],1,&Y[2*r*F],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fft_s: dim must be 0 or 1.\n"); return 1;
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int fft_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft)
{
    const double z = 0.0;
    const int F = nfft/2 + 1;
    int r, c;
    double *X1;
    fftw_complex *Y1;
    fftw_plan plan;

    //Checks
    if (R<1) { fprintf(stderr,"error in fft_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fft_d: C (ncols X) must be positive\n"); return 1; }
    if (nfft<R && dim==0) { fprintf(stderr,"error in fft_d: nfft must be >= R (winlength)\n"); return 1; }
    if (nfft<C && dim==1) { fprintf(stderr,"error in fft_d: nfft must be >= C (winlength)\n"); return 1; }

    //Initialize fftw
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_complex((size_t)F);
    plan = fftw_plan_dft_r2c_1d(nfft,X1,Y1,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_d: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(F,(double *)&Y1[0],1,&Y[2*c*F],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(F,(double *)&Y1[0],1,&Y[2*c],C);
            }
        }
    }
    else if (dim==1)
    {
        cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r],R,&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(F,(double *)&Y1[0],1,&Y[2*r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(F,(double *)&Y1[0],1,&Y[2*r*F],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fft_d: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


int fft_c (float *Y, const float *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft)
{
    const float z[2] = {0.0f,0.0f};
    int r, c;
    fftwf_complex *X1, *Y1;
    fftwf_plan plan;

    //Checks
    if (R<1) { fprintf(stderr,"error in fft_c: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fft_c: C (ncols X) must be positive\n"); return 1; }
    if (nfft<R && dim==0) { fprintf(stderr,"error in fft_c: nfft must be >= R (winlength)\n"); return 1; }
    if (nfft<C && dim==1) { fprintf(stderr,"error in fft_c: nfft must be >= C (winlength)\n"); return 1; }

    //Initialize fftwf
    X1 = fftwf_alloc_complex((size_t)nfft);
    Y1 = fftwf_alloc_complex((size_t)nfft);
    plan = fftwf_plan_dft_1d(nfft,X1,Y1,FFTW_FORWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_c: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        cblas_ccopy(nfft-R,&z[0],0,(float *)&X1[R],1); //zero-pad
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c*R],1,(float *)&X1[0],1);
                fftwf_execute(plan);
                cblas_ccopy(nfft,(float *)&Y1[0],1,&Y[2*c*nfft],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_ccopy(R,&X[2*c],C,(float *)&X1[0],1);
                fftwf_execute(plan);
                cblas_ccopy(nfft,(float *)&Y1[0],1,&Y[2*c],C);
            }
        }
    }
    else if (dim==1)
    {
        cblas_ccopy(nfft-C,&z[0],0,(float *)&X1[C],1);
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r],R,(float *)&X1[0],1);
                fftwf_execute(plan);
                cblas_ccopy(nfft,(float *)&Y1[0],1,&Y[2*r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_ccopy(C,&X[2*r*C],1,(float *)&X1[0],1);
                fftwf_execute(plan);
                cblas_ccopy(nfft,(float *)&Y1[0],1,&Y[2*r*nfft],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fft_c: dim must be 0 or 1.\n"); return 1;
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int fft_z (double *Y, const double *X, const char iscolmajor, const int R, const int C, const size_t dim, const int nfft)
{
    const double z[2] = {0.0,0.0};
    int r, c;
    fftw_complex *X1, *Y1;
    fftw_plan plan;

    //Checks
    if (R<1) { fprintf(stderr,"error in fft_z: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fft_z: C (ncols X) must be positive\n"); return 1; }
    if (nfft<R && dim==0) { fprintf(stderr,"error in fft_z: nfft must be >= R (winlength)\n"); return 1; }
    if (nfft<C && dim==1) { fprintf(stderr,"error in fft_z: nfft must be >= C (winlength)\n"); return 1; }

    //Initialize fftw
    X1 = fftw_alloc_complex((size_t)nfft);
    Y1 = fftw_alloc_complex((size_t)nfft);
    plan = fftw_plan_dft_1d(nfft,X1,Y1,FFTW_FORWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_z: problem creating fftw plan"); return 1; }

    if (dim==0)
    {
        cblas_zcopy(nfft-R,&z[0],0,(double *)&X1[R],1); //zero-pad
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c*R],1,(double *)&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(nfft,(double *)&Y1[0],1,&Y[2*c*nfft],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_zcopy(R,&X[2*c],C,(double *)&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(nfft,(double *)&Y1[0],1,&Y[2*c],C);
            }
        }
    }
    else if (dim==1)
    {
        cblas_zcopy(nfft-C,&z[0],0,(double *)&X1[C],1);
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r],R,(double *)&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(nfft,(double *)&Y1[0],1,&Y[2*r],R);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_zcopy(C,&X[2*r*C],1,(double *)&X1[0],1);
                fftw_execute(plan);
                cblas_zcopy(nfft,(double *)&Y1[0],1,&Y[2*r*nfft],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fft_z: dim must be 0 or 1.\n"); return 1;
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
