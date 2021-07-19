//This takes the output of window_univar, and applies the 1D FFT to each row or col.
//The input dim is which dim to take the 1D FFT along (0->along cols, 1->along rows).
//The output Y is power (real-valued), i.e. the FFT squared, Y = X*conj(X) element-wise.
//The size of Y must be FxC if dim==0, and RxF if dim==1, where F = nfft/2 + 1.

//I tried parallel versions with OpenMP, but much slower (have to make fftw_plan P times!).

#include <stdio.h>
#include <cblas.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int fft_squared_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);
int fft_squared_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft);


int fft_squared_s (float *Y, const float *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft)
{
    const float z = 0.0f;
    const int F = nfft/2 + 1;
    int r, c, f;
    float *X1, *Y1;
    fftwf_plan plan;
    //struct timespec tic, toc;

    //Checks
    if (R<1) { fprintf(stderr,"error in fft_squared_s: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fft_squared_s: C (ncols X) must be positive\n"); return 1; }
    if (nfft<R && dim==0) { fprintf(stderr,"error in fft_squared_s: nfft must be >= R (winlength)\n"); return 1; }
    if (nfft<C && dim==1) { fprintf(stderr,"error in fft_squared_s: nfft must be >= C (winlength)\n"); return 1; }

    //Initialize fftwf
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    plan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_squared_s: problem creating fftw plan\n"); return 1; }
    cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad

    if (dim==0u)
    {
        if (iscolmajor)
        {
            //clock_gettime(CLOCK_REALTIME,&tic);
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c*R],1,&X1[0],1);
                fftwf_execute(plan);
                
                //n = c*F; Y[n++] = fmaf(Y1[0],Y1[0],0.0f); //this was ~25% slower
                //for (f=1; f<F-1; f++, n++) { Y[n] = fmaf(Y1[f],Y1[f],fmaf(Y1[nfft-f],Y1[nfft-f],0.0f)); }
                //Y[n] = (nfft%2) ? fmaf(Y1[f],Y1[f],fmaf(Y1[nfft-f],Y1[nfft-f],0.0f)) : fmaf(Y1[f],Y1[f],0.0f);
                
                //n = c*F; Y[n++] = Y1[0]*Y1[0];       //this was also ~25% slower
                //for (f=1; f<F-1; f++, n++) { Y[n] = Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f]; }
                //Y[n] = (nfft%2) ? Y1[f]*Y1[f] + Y1[nfft-f]*Y1[nfft-f] : Y1[f]*Y1[f];

                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_scopy(F,&Y[c*F],1,&Y1[0],1);
            }
            //clock_gettime(CLOCK_REALTIME,&toc);
            //fprintf(stderr,"elapsed time = %f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_scopy(R,&X[c],C,&X1[0],1);
                fftwf_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_scopy(F,&Y[c],C,&Y1[0],1);
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r],R,&X1[0],1);
                fftwf_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_scopy(F,&Y[r],R,&Y1[0],1);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_scopy(C,&X[r*C],1,&X1[0],1);
                fftwf_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_scopy(F,&Y[r*F],1,&Y1[0],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fft_squared_s: dim must be 0 or 1.\n"); return 1;
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int fft_squared_d (double *Y, const double *X, const char iscolmajor, const int R, const int C, const int dim, const int nfft)
{
    const double z = 0.0;
    const int F = nfft/2 + 1;
    int r, c, f;
    double *X1, *Y1;
    fftw_plan plan;
    
    //Checks
    if (R<1) { fprintf(stderr,"error in fft_squared_d: R (nrows X) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in fft_squared_d: C (ncols X) must be positive\n"); return 1; }
    if (nfft<R && dim==0) { fprintf(stderr,"error in fft_squared_s: nfft must be >= R (winlength)\n"); return 1; }
    if (nfft<C && dim==1) { fprintf(stderr,"error in fft_squared_s: nfft must be >= C (winlength)\n"); return 1; }


    //Initialize fftw
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    plan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_squared_d: problem creating fftw plan\n"); return 1; }
    cblas_dcopy(nfft,&z,0,&X1[0],1); //zero-pad

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c*R],1,&X1[0],1);
                fftw_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_dcopy(F,&Y[c*F],1,&Y1[0],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                cblas_dcopy(R,&X[c],C,&X1[0],1);
                fftw_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_dcopy(F,&Y[c],C,&Y1[0],1);
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r],R,&X1[0],1);
                fftw_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_dcopy(F,&Y[r],R,&Y1[0],1);
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                cblas_dcopy(C,&X[r*C],1,&X1[0],1);
                fftw_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_dcopy(F,&Y[r*F],1,&Y1[0],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in fft_squared_d: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

