//STFT (short-term Fourier transform) of univariate time series X1.

//This takes a continuous univariate time series (vector X1) of length N,
//and a window (vector X2) of length L,
//and outputs the STFT (matrix Y) at each of W windows.
//If Y is row-major, then it has size W x F.
//If Y is col-major, then it has size F x W.
//where F is nfft/2+1, and nfft is the next-pow-2 of L.

//This follows Kaldi conventions for compatibility.
//For framing conventions, see window_univar.c.
//For more flexibility, see also window_univar_float.c.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int stft_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges);
int stft_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges);


int stft_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges)
{
    if (L<1u) { fprintf(stderr,"error in stft_s: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in stft_s: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in stft_s: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while (nfft<L) { nfft *= 2u; }
        size_t F = nfft/2u + 1u;

        //Initialize fftwf
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(nfft*sizeof(float));
        fftwf_plan plan = fftwf_plan_r2r_1d((int)nfft,X1,Y1,FFTW_R2HC,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in fft_fftw_r2hc_s: problem creating fftw plan"); return 1; }

        if (snip_edges)
        {
            const int xd = (int)L - (int)stp;
            for (size_t w=0u; w<W; ++w, X1-=xd, X2-=L)
            {
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y) { *Y = *X1 * *X2; }
            }
        }
        else
        {
            const int xd = (int)L - (int)stp;           //X inc after each frame
            const size_t Lpre = L/2u;                   //nsamps before center samp
            int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
            int n, prev_n = 0;                          //current/prev samps in X

            for (size_t w=0u; w<W; ++w, ss+=stp, X2-=L)
            {
                if (ss<0 || ss>(int)N-(int)L)
                {
                    for (int s=ss; s<ss+(int)L; ++s, ++X2, ++Y)
                    {
                        n = s; //This ensures extrapolation by signal reversal to any length
                        while (n<0 || n>=(int)N) { n = (n<0) ? -n : (n<(int)N) ? n : 2*(int)N-2-n; }
                        X1 += n - prev_n;
                        *Y = *X1 * *X2;
                        prev_n = n;
                    }
                }
                else
                {
                    X1 += ss - prev_n;
                    for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y) { *Y = *X1 * *X2; }
                    X1 -= xd; prev_n = ss + stp;
                }
            }
        }
    }

    return 0;
}


int stft_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t stp, const int snip_edges)
{
    if (L<1u) { fprintf(stderr,"error in stft_d: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in stft_d: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in stft_d: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else if (snip_edges)
    {
        const int xd = (int)L - (int)stp;
        for (size_t w=0u; w<W; ++w, X1-=xd, X2-=L)
        {
            for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y) { *Y = *X1 * *X2; }
        }
    }
    else
    {
        const int xd = (int)L - (int)stp;           //X inc after each frame
        const size_t Lpre = L/2u;                   //nsamps before center samp
        int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
        int n, prev_n = 0;                          //current/prev samps in X

        for (size_t w=0u; w<W; ++w, ss+=stp, X2-=L)
        {
            if (ss<0 || ss>(int)N-(int)L)
            {
                for (int s=ss; s<ss+(int)L; ++s, ++X2, ++Y)
                {
                    n = s; //This ensures extrapolation by signal reversal to any length
                    while (n<0 || n>=(int)N) { n = (n<0) ? -n : (n<(int)N) ? n : 2*(int)N-2-n; }
                    X1 += n - prev_n;
                    *Y = *X1 * *X2;
                    prev_n = n;
                }
            }
            else
            {
                X1 += ss - prev_n;
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Y) { *Y = *X1 * *X2; }
                X1 -= xd; prev_n = ss + stp;
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
