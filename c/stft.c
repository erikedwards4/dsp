//This takes univariate X, windows X, and then applies the 1D FFT to each frame.
//The output Y is power (real-valued), i.e. the FFT squared.
//The size of Y is FxC if dim==0, and RxF if dim==1, where F = nfft/2 + 1.

#include <stdio.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int stft_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const int dim, const int c0, const float stp, const char mn0, const int nfft);
int stft_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const int dim, const int c0, const double stp, const char mn0, const int nfft);


int stft_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int N, const float *W, const int L, const int dim, const int c0, const float stp, const char mn0, const int nfft)
{
    const float z = 0.0f, o = 1.0f;
    const int F = nfft/2 + 1;
    const int Lpre = L/2;      //nsamps before center samp
    int ss = c0 - Lpre;        //start samp of current frame
    int r, c, f;
    float *X1, *Y1;
    fftwf_plan plan;

    //Checks
    if (N<1) { fprintf(stderr,"error in stft_s: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in stft_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stft_s: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in stft_s: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in stft_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in stft_s: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in stft_s: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0f) { fprintf(stderr,"error in stft_s: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in stft_s: nfft must be >= L (winlength)\n"); return 1; }
    if (dim==0 && F!=R) { fprintf(stderr,"error in stft_s: F must equal nrows in X for dim==0\n"); return 1; }
    if (dim==1 && F!=C) { fprintf(stderr,"error in stft_s: F must equal ncols in X for dim==1\n"); return 1; }

    //Initialize fftwf
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    plan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in stft_s: problem creating fftw plan\n"); return 1; }
    cblas_scopy(nfft,&z,0,&X1[0],1); //zero-pad

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                ss = (int)(roundf(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_s: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
                fftwf_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_scopy(F,&Y[c*F],1,&Y1[0],1);
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                ss = (int)(roundf(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_s: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
                fftwf_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_scopy(F,&Y[c],C,&Y1[0],1);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                ss = (int)(roundf(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasColMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasColMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_s: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
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
                ss = (int)(roundf(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_scopy(-ss,&z,0,&X1[0],1);
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0f,&W[-ss],1,&X[0],1,0.0f,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,L,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_ssbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0f,&W[0],1,&X[ss],1,0.0f,&X1[0],1);
                    cblas_scopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_scopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_s: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_saxpy(L,-cblas_sdot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
                fftwf_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_scopy(F,&Y[r*F],1,&Y1[0],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stft_s: dim must be 0 or 1.\n"); return 1;
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int stft_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int N, const double *W, const int L, const int dim, const int c0, const double stp, const char mn0, const int nfft)
{
    const double z = 0.0, o = 1.0;
    const int F = nfft/2 + 1;
    const int Lpre = L/2;      //nsamps before center samp
    int ss = c0 - Lpre;        //start samp of current frame
    int r, c, f;
    double *X1, *Y1;
    fftw_plan plan;

    //Checks
    if (N<1) { fprintf(stderr,"error in stft_d: N (length X) must be positive\n"); return 1; }
    if (R<1) { fprintf(stderr,"error in stft_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in stft_d: C (ncols Y) must be positive\n"); return 1; }
    if (c0<0) { fprintf(stderr,"error in stft_d: c0 (center samp of 1st frame) must be nonnegative\n"); return 1; }
    if (c0>N-1) { fprintf(stderr,"error in stft_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (L<1) { fprintf(stderr,"error in stft_d: L (winlength) must be positive\n"); return 1; }
    if (L>=N) { fprintf(stderr,"error in stft_d: L (winlength) must be < N (length X)\n"); return 1; }
    if (stp<=0.0) { fprintf(stderr,"error in stft_d: stp (step size) must be positive\n"); return 1; }
    if (nfft<L) { fprintf(stderr,"error in stft_d: nfft must be >= L (winlength)\n"); return 1; }
    if (dim==0 && F!=R) { fprintf(stderr,"error in stft_d: F must equal nrows in X for dim==0\n"); return 1; }
    if (dim==1 && F!=C) { fprintf(stderr,"error in stft_d: F must equal ncols in X for dim==1\n"); return 1; }

    //Initialize fftw
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    plan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in stft_d: problem creating fftw plan\n"); return 1; }
    cblas_dcopy(nfft,&z,0,&X1[0],1); //zero-pad

    if (dim==0)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                ss = (int)(round(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_d: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
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
                ss = (int)(round(c*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_d: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
                fftw_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_dcopy(F,&Y[c],C,&Y1[0],1);
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                ss = (int)(round(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasColMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasColMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_d: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
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
                ss = (int)(round(r*stp)) + c0 - Lpre;
                if (ss<0)
                {
                    cblas_dcopy(-ss,&z,0,&X1[0],1);
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L+ss,0,1.0,&W[-ss],1,&X[0],1,0.0,&X1[-ss],1);
                }
                else if (ss+L<=N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,L,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                }
                else if (ss<N)
                {
                    cblas_dsbmv(CblasRowMajor,CblasUpper,N-ss,0,1.0,&W[0],1,&X[ss],1,0.0,&X1[0],1);
                    cblas_dcopy(L-N+ss,&z,0,&X1[N-ss],1);
                }
                else
                {
                    cblas_dcopy(L,&z,0,&X1[0],1);
                    //fprintf(stderr,"error in stft_d: frame start samp out of range"); return 1;
                }
                if (mn0) { cblas_daxpy(L,-cblas_ddot(L,&X1[0],1,&o,0)/L,&o,0,&X1[0],1); }
                fftw_execute(plan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1+nfft%2; f++) { Y1[f] += Y1[nfft-f]; }
                cblas_dcopy(F,&Y[r*F],1,&Y1[0],1);
            }
        }
    }
    else
    {
        fprintf(stderr,"error in stft_d: dim must be 0 or 1.\n"); return 1;
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif

