//This gets polynomial coeffs for each row/col of X.
//The autocov is the biased version (e.g. default of Octave, etc.).

//An opt mean0 is added to zero the mean of each row/col of X first.
//In this case, this is mean0 -> autocov_fft -> lev_durb.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int sig2poly_levdurb_s (float *Y, float *V, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const char mean0);
int sig2poly_levdurb_d (double *Y, double *V, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const char mean0);


int sig2poly_levdurb_s (float *Y, float *V, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const char mean0)
{
    if (dim>3) { fprintf(stderr,"error in sig2poly_levdurb_s: dim must be in [0 3]\n"); return 1; }

    const float z = 0.0f, o = 1.0f;
    const size_t P1 = P + 1; //size of Y along dim
    float m, g, sc;
    int r, c, f, nfft, F, p, q;
    float *X1, *Y1, *AStmp;
    fftwf_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2poly_levdurb_s: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2poly_levdurb_s: ncols X must be positive\n"); return 1; }
    if (dim==0 && P>=R) { fprintf(stderr,"error in sig2poly_levdurb_s: P must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in sig2poly_levdurb_s: P must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+P1 : C+P1;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }
    F = nfft/2 + 1;
    sc = 1.0f/nfft;

    //Initialize
    X1 = fftwf_alloc_real((size_t)nfft);
    Y1 = fftwf_alloc_real((size_t)nfft);
    fplan = fftwf_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    iplan = fftwf_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in sig2poly_levdurb_s: problem creating fftw plan"); return 1; }
    if (!(AStmp=(float *)malloc((size_t)(P)*sizeof(float)))) { fprintf(stderr,"error in sig2poly_levdurb_s: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_scopy((int)C,&o,0,Y,P1);
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_scopy((int)R,&X[c*R],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot((int)R,&X1[0],1,&o,0) / R;
                    cblas_saxpy((int)R,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                Y[c*P1+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P1,sc,&X1[0],1);
                V[c] = fmaf(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fmaf(X1[q+1],Y[p-q+c*P1],g); }
                    Y[p+1+c*P1] = g = -g/V[c];
                    for (q=0; q<p; q++) { Y[q+1+c*P1] = fmaf(g,AStmp[p-q-1],Y[q+1+c*P1]); }
                    cblas_scopy(p,&Y[1+c*P1],1,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy((int)C,&o,0,Y,1);
            for (size_t c=0; c<C; c++)
            {
                cblas_scopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_scopy((int)R,&X[c],(int)C,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot((int)R,&X1[0],1,&o,0) / R;
                    cblas_saxpy((int)R,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                Y[c+C] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P1,sc,&X1[0],1);
                V[c] = fmaf(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fmaf(X1[q+1],Y[c+(p-q)*C],g); }
                    Y[c+(p+1)*C] = g = -g/V[c];
                    for (q=0; q<p; q++) { Y[c+(q+1)*C] = fmaf(g,AStmp[p-q-1],Y[c+(q+1)*C]); }
                    cblas_scopy(p,&Y[c+C],(int)C,&AStmp[0],1);
                    V[c] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_scopy((int)R,&o,0,Y,1);
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_scopy((int)C,&X[r],(int)R,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot((int)C,&X1[0],1,&o,0) / C;
                    cblas_saxpy((int)C,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                Y[r+R] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P1,sc,&X1[0],1);
                V[r] = fmaf(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fmaf(X1[q+1],Y[r+(p-q)*R],g); }
                    Y[r+(p+1)*R] = g = -g/V[r];
                    for (q=0; q<p; q++) { Y[r+(q+1)*R] = fmaf(g,AStmp[p-q-1],Y[r+(q+1)*R]); }
                    cblas_scopy(p,&Y[r+R],(int)R,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
        else
        {
            cblas_scopy((int)R,&o,0,Y,P1);
            for (size_t r=0; r<R; r++)
            {
                cblas_scopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_scopy((int)C,&X[r*C],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_sdot((int)C,&X1[0],1,&o,0) / C;
                    cblas_saxpy((int)C,-m,&o,0,&X1[0],1);
                }
                fftwf_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftwf_execute(iplan);
                Y[r*P1+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_sscal(P1,sc,&X1[0],1);
                V[r] = fmaf(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fmaf(X1[q+1],Y[p-q+r*P1],g); }
                    Y[p+1+r*P1] = g = -g/V[r];
                    for (q=0; q<p; q++) { Y[q+1+r*P1] = fmaf(g,AStmp[p-q-1],Y[q+1+r*P1]); }
                    cblas_scopy(p,&Y[1+r*P],1,&AStmp[0],1);
                    V[r] *= fmaf(g,-g,1.0f);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2poly_levdurb_s: dim must be 0 or 1.\n"); return 1;
    }

    fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int sig2poly_levdurb_d (double *Y, double *V, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t P, const char mean0)
{
    if (dim>3) { fprintf(stderr,"error in sig2poly_levdurb_d: dim must be in [0 3]\n"); return 1; }

    const double z = 0.0, o = 1.0;
    const size_t P1 = P + 1; //size of Y along dim
    double m, g, sc;
    int r, c, f, nfft, F, p, q;
    double *X1, *Y1, *AStmp;
    fftw_plan fplan, iplan;

    //Checks
    if (R<1) { fprintf(stderr,"error in sig2poly_levdurb_d: nrows X must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in sig2poly_levdurb_d: ncols X must be positive\n"); return 1; }
    if (dim==0 && P>=R) { fprintf(stderr,"error in sig2poly_levdurb_d: P must be < nrows X for dim==0\n"); return 1; }
    if (dim==1 && P>=C) { fprintf(stderr,"error in sig2poly_levdurb_d: P must be < ncols X for dim==1\n"); return 1; }

    //Get nfft
    nfft = (dim==0) ? R+P1 : C+P1;
    if (nfft>16384) { nfft += nfft%2; }
    else { f = 1; while (f<nfft) { f *= 2; } nfft = f; }
    F = nfft/2 + 1;
    sc = 1.0/nfft;

    //Initialize
    X1 = fftw_alloc_real((size_t)nfft);
    Y1 = fftw_alloc_real((size_t)nfft);
    fplan = fftw_plan_r2r_1d(nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
    iplan = fftw_plan_r2r_1d(nfft,Y1,X1,FFTW_R2HC,FFTW_ESTIMATE);
    if (!fplan || !iplan) { fprintf(stderr,"error in sig2poly_levdurb_d: problem creating fftw plan"); return 1; }
    if (!(AStmp=(double *)malloc((size_t)(P)*sizeof(double)))) { fprintf(stderr,"error in sig2poly_levdurb_d: problem with malloc. "); perror("malloc"); return 1; }

    if (dim==0)
    {
        if (iscolmajor)
        {
            cblas_dcopy((int)C,&o,0,Y,P1);
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_dcopy((int)R,&X[c*R],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot((int)R,&X1[0],1,&o,0) / R;
                    cblas_daxpy((int)R,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                Y[c*P1+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P1,sc,&X1[0],1);
                V[c] = fma(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fma(X1[q+1],Y[p-q+c*P1],g); }
                    Y[p+1+c*P1] = g = -g/V[c];
                    for (q=0; q<p; q++) { Y[q+1+c*P1] = fma(g,AStmp[p-q-1],Y[q+1+c*P1]); }
                    cblas_dcopy(p,&Y[1+c*P1],1,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy((int)C,&o,0,Y,1);
            for (size_t c=0; c<C; c++)
            {
                cblas_dcopy(nfft-R,&z,0,&X1[R],1); //zero-pad
                cblas_dcopy((int)R,&X[c],(int)C,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot((int)R,&X1[0],1,&o,0) / R;
                    cblas_daxpy((int)R,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                Y[c+C] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P1,sc,&X1[0],1);
                V[c] = fma(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fma(X1[q+1],Y[c+(p-q)*C],g); }
                    Y[c+(p+1)*C] = g = -g/V[c];
                    for (q=0; q<p; q++) { Y[c+(q+1)*C] = fma(g,AStmp[p-q-1],Y[c+(q+1)*C]); }
                    cblas_dcopy(p,&Y[c+C],(int)C,&AStmp[0],1);
                    V[c] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else if (dim==1)
    {
        if (iscolmajor)
        {
            cblas_dcopy((int)R,&o,0,Y,1);
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_dcopy((int)C,&X[r],(int)R,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot((int)C,&X1[0],1,&o,0) / C;
                    cblas_daxpy((int)C,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                Y[r+R] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P1,sc,&X1[0],1);
                V[r] = fma(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fma(X1[q+1],Y[r+(p-q)*R],g); }
                    Y[r+(p+1)*R] = g = -g/V[r];
                    for (q=0; q<p; q++) { Y[r+(q+1)*R] = fma(g,AStmp[p-q-1],Y[r+(q+1)*R]); }
                    cblas_dcopy(p,&Y[r+R],(int)R,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
        else
        {
            cblas_dcopy((int)R,&o,0,Y,P1);
            for (size_t r=0; r<R; r++)
            {
                cblas_dcopy(nfft-C,&z,0,&X1[C],1); //zero-pad
                cblas_dcopy((int)C,&X[r*C],1,&X1[0],1);
                if (mean0)
                {
                    m = cblas_ddot((int)C,&X1[0],1,&o,0) / C;
                    cblas_daxpy((int)C,-m,&o,0,&X1[0],1);
                }
                fftw_execute(fplan);
                for (f=0; f<nfft; f++) { Y1[f] *= Y1[f]; }
                for (f=1; f<F-1; f++) { Y1[f] += Y1[nfft-f]; Y1[nfft-f] = Y1[f]; }
                fftw_execute(iplan);
                Y[r*P1+1] = AStmp[0] = g = -X1[1]/X1[0];
                cblas_dscal(P1,sc,&X1[0],1);
                V[r] = fma(X1[1],g,X1[0]);
                for (p=1; p<P; p++)
                {
                    g = X1[p+1];
                    for (q=0; q<p; q++) { g = fma(X1[q+1],Y[p-q+r*P1],g); }
                    Y[p+1+r*P1] = g = -g/V[r];
                    for (q=0; q<p; q++) { Y[q+1+r*P1] = fma(g,AStmp[p-q-1],Y[q+1+r*P1]); }
                    cblas_dcopy(p,&Y[1+r*P1],1,&AStmp[0],1);
                    V[r] *= fma(g,-g,1.0);
                }
            }
        }
    }
    else
    {
        fprintf(stderr,"error in sig2poly_levdurb_d: dim must be 0 or 1.\n"); return 1;
    }

    fftw_destroy_plan(fplan); fftw_destroy_plan(iplan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
