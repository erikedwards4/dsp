//This squares half-complex (hc) input X element-wise,
//equivalent to: |X|.^2, i.e. Xr*Xr + Xi*Xi.
//X is in the format of the output of fft_hc.

//The half-complex format is a more efficient format for the
//FFT output given by fftw3 library, using FFTW_R2HC setting.

//This does not use in-place, since size of output is different from input.
//The input (X) has frames with nfft points, whereas
//the output (Y) has frames with F points,
//where F = nfft/2 + 1 = num non-negative FFT freqs.

#include <stdio.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace ov {
extern "C" {
#endif

int hc_square_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int nfft, const int dim);
int hc_square_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int nfft, const int dim);


int hc_square_s (float *Y, const char iscolmajor, const int R, const int C, const float *X, const int nfft, const int dim)
{
    const int F = nfft/2 + 1;
    int r, c, f, nx, ny;

    //Checks
    if (R<1) { fprintf(stderr,"error in hc_square_s: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in hc_square_s: C (ncols Y) must be positive\n"); return 1; }
    if (nfft<1) { fprintf(stderr,"error in hc_square_s: nfft must be positive\n"); return 1; }
    if (dim==0 && F!=R) { fprintf(stderr,"error in hc_square_d: nrows in Y must = F = nfft/2+1\n"); return 1; }
    if (dim==1 && F!=C) { fprintf(stderr,"error in hc_square_d: ncols in Y must = F = nfft/2+1\n"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                nx = c*nfft; ny = c*F; Y[ny++] = X[nx]*X[nx];
                for (f=1; f<F-1; f++, ny++) { Y[ny] = X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f]; }
                Y[ny] = (nfft%2) ? X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f] : X[nx+f]*X[nx+f];
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                nx = ny = c; Y[ny] = X[nx]*X[nx]; ny += C;
                for (f=1; f<F-1; f++, ny+=C) { Y[ny] = X[nx+f*C]*X[nx+f*C] + X[nx+(nfft-f)*C]*X[nx+(nfft-f)*C]; }
                Y[ny] = (nfft%2) ? X[nx+f*C]*X[nx+f*C] + X[nx+(nfft-f)*C]*X[nx+(nfft-f)*C] : X[nx+f*C]*X[nx+f*C];
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                nx = ny = r; Y[ny] = X[nx]*X[nx]; ny += R;
                for (f=1; f<F-1; f++, ny+=R) { Y[ny] = X[nx+f*R]*X[nx+f*R] + X[nx+(nfft-f)*R]*X[nx+(nfft-f)*R]; }
                Y[ny] = (nfft%2) ? X[nx+f*R]*X[nx+f*R] + X[nx+(nfft-f)*R]*X[nx+(nfft-f)*R] : X[nx+f*R]*X[nx+f*R];
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                nx = r*nfft; ny = r*F; Y[ny++] = X[nx]*X[nx];
                for (f=1; f<F-1; f++, ny++) { Y[ny] = X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f]; }
                Y[ny] = (nfft%2) ? X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f] : X[nx+f]*X[nx+f];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in hc_square_s: dim must be 0 or 1.\n"); return 1;
    }

    return 0;
}


int hc_square_d (double *Y, const char iscolmajor, const int R, const int C, const double *X, const int nfft, const int dim)
{
    const int F = nfft/2 + 1;
    int r, c, f, nx, ny;

    //Checks
    if (R<1) { fprintf(stderr,"error in hc_square_d: R (nrows Y) must be positive\n"); return 1; }
    if (C<1) { fprintf(stderr,"error in hc_square_d: C (ncols Y) must be positive\n"); return 1; }
    if (nfft<1) { fprintf(stderr,"error in hc_square_s: nfft must be positive\n"); return 1; }
    if (dim==0 && F!=R) { fprintf(stderr,"error in hc_square_d: nrows in Y must equal F = nfft/2+1\n"); return 1; }
    if (dim==1 && F!=C) { fprintf(stderr,"error in hc_square_d: ncols in Y must equal F = nfft/2+1\n"); return 1; }

    if (dim==0u)
    {
        if (iscolmajor)
        {
            for (c=0; c<C; c++)
            {
                nx = c*nfft; ny = c*F; Y[ny++] = X[nx]*X[nx];
                for (f=1; f<F-1; f++, ny++) { Y[ny] = X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f]; }
                Y[ny] = (nfft%2) ? X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f] : X[nx+f]*X[nx+f];
            }
        }
        else
        {
            for (c=0; c<C; c++)
            {
                nx = ny = c; Y[ny] = X[nx]*X[nx]; ny += C;
                for (f=1; f<F-1; f++, ny+=C) { Y[ny] = X[nx+f*C]*X[nx+f*C] + X[nx+(nfft-f)*C]*X[nx+(nfft-f)*C]; }
                Y[ny] = (nfft%2) ? X[nx+f*C]*X[nx+f*C] + X[nx+(nfft-f)*C]*X[nx+(nfft-f)*C] : X[nx+f*C]*X[nx+f*C];
            }
        }
    }
    else if (dim==1u)
    {
        if (iscolmajor)
        {
            for (r=0; r<R; r++)
            {
                nx = ny = r; Y[ny] = X[nx]*X[nx]; ny += R;
                for (f=1; f<F-1; f++, ny+=R) { Y[ny] = X[nx+f*R]*X[nx+f*R] + X[nx+(nfft-f)*R]*X[nx+(nfft-f)*R]; }
                Y[ny] = (nfft%2) ? X[nx+f*R]*X[nx+f*R] + X[nx+(nfft-f)*R]*X[nx+(nfft-f)*R] : X[nx+f*R]*X[nx+f*R];
            }
        }
        else
        {
            for (r=0; r<R; r++)
            {
                nx = r*nfft; ny = r*F; Y[ny++] = X[nx]*X[nx];
                for (f=1; f<F-1; f++, ny++) { Y[ny] = X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f]; }
                Y[ny] = (nfft%2) ? X[nx+f]*X[nx+f] + X[nx+nfft-f]*X[nx+nfft-f] : X[nx+f]*X[nx+f];
            }
        }
    }
    else
    {
        fprintf(stderr,"error in hc_square_d: dim must be 0 or 1.\n"); return 1;
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif

