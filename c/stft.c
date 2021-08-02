//STFT (short-term Fourier transform) of univariate time series X1,
//using window X2, and outputing power in Y.

//This takes a continuous univariate time series (vector X1) of length N,
//and a window (vector X2) of length L,
//and outputs the STFT (matrix Y) at each of W windows.
//If Y is row-major, then it has size W x F.
//If Y is col-major, then it has size F x W.
//where F is nfft/2+1, and nfft is the next-pow-2 of L.

//This follows Kaldi conventions for compatibility.
//For framing conventions, see window_univar.c.
//For more flexibility, see also window_univar_float.c.

//The following boolean (int) options are added:
//mn0: subtract mean from each frame just after windowing.
//amp: take sqrt of each element of Y just after getting power.
//lg:  take log of each element of Y just before output.

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int stft_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t nfft, const size_t stp, const int snip_edges, const int mn0, const int amp, const int lg);
int stft_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t nfft, const size_t stp, const int snip_edges, const int mn0, const int amp, const int lg);


int stft_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t nfft, const size_t stp, const int snip_edges, const int mn0, const int amp, const int lg)
{
    if (L<1u) { fprintf(stderr,"error in stft_s: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in stft_s: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in stft_s: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else
    {
        const size_t F = nfft/2u + 1u;          //Num non-negative FFT freqs
        const int xd = (int)L - (int)stp;       //X1 inc after each frame
        float mn;                               //Mean of Xw (one window of X)
        
        //Initialize FFT
        float *Xw, *Yw;
        Xw = (float *)fftwf_malloc(nfft*sizeof(float));
        Yw = (float *)fftwf_malloc(nfft*sizeof(float));
        fftwf_plan plan = fftwf_plan_r2r_1d((int)nfft,Xw,Yw,FFTW_R2HC,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in stft_s: problem creating fftw plan"); return 1; }
        for (size_t l=0u; l<nfft; ++l, ++Xw) { *Xw = 0.0f; }
        Xw -= nfft;

        if (snip_edges)
        {
            for (size_t w=0u; w<W; ++w)
            {
                //Window
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                X1 -= xd; X2 -= L;
                
                //Zero mean
                if (mn0)
                {
                    mn = 0.0f;
                    for (size_t l=0u; l<L; ++l) { mn += *--Xw; }
                    mn /= (float)L;
                    for (size_t l=0u; l<L; ++l, ++Xw) { *Xw -= mn; }
                }
                
                //FFT
                Xw -= L;
                fftwf_execute(plan);
                
                //Power
                for (size_t f=0u; f<F; ++f, ++Yw, ++Y) { *Y = *Yw * *Yw; }
                Y -= 2u;
                for (size_t f=1u; f<F-1u; ++f, ++Yw, --Y) { *Y += *Yw * *Yw; }
                Yw -= nfft;

                //Amplitude and/or log
                if (amp)
                {
                    if (lg) { for (size_t f=0u; f<F; ++f, ++Y) { *Y = logf(sqrtf(*Y)); } }
                    else { for (size_t f=0u; f<F; ++f, ++Y) { *Y = sqrtf(*Y); } }
                }
                else if (lg)
                {
                    for (size_t f=0u; f<F; ++f, ++Y) { *Y = logf(*Y); }
                }
                else { Y += F; }
            }
        }
        else
        {
            const size_t Lpre = L/2u;                   //nsamps before center samp
            int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
            int n, prev_n = 0;                          //current/prev samps in X

            for (size_t w=0u; w<W; ++w)
            {
                //Window
                if (ss<0 || ss>(int)N-(int)L)
                {
                    for (int s=ss; s<ss+(int)L; ++s, ++X2, ++Xw)
                    {
                        //Window
                        n = s; //This ensures extrapolation by signal reversal to any length
                        while (n<0 || n>=(int)N) { n = (n<0) ? -n-1 : (n<(int)N) ? n : 2*(int)N-1-n; }
                        X1 += n - prev_n; prev_n = n;
                        *Xw = *X1 * *X2;
                    }
                }
                else
                {
                    X1 += ss - prev_n;
                    for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    X1 -= xd; prev_n = ss + (int)stp;
                }
                ss += stp; X2 -= L;

                //Zero mean
                if (mn0)
                {
                    mn = 0.0f;
                    for (size_t l=0u; l<L; ++l) { mn += *--Xw; }
                    mn /= (float)L;
                    for (size_t l=0u; l<L; ++l, ++Xw) { *Xw -= mn; }
                }

                //FFT
                Xw -= L;
                fftwf_execute(plan);
                
                //Power
                for (size_t f=0u; f<F; ++f, ++Yw, ++Y) { *Y = *Yw * *Yw; }
                Y -= 2u;
                for (size_t f=1u; f<F-1u; ++f, ++Yw, --Y) { *Y += *Yw * *Yw; }
                Yw -= nfft;
                
                //Amplitude and/or log
                if (amp)
                {
                    if (lg) { for (size_t f=0u; f<F; ++f, ++Y) { *Y = logf(sqrtf(*Y)); } }
                    else { for (size_t f=0u; f<F; ++f, ++Y) { *Y = sqrtf(*Y); } }
                }
                else if (lg)
                {
                    for (size_t f=0u; f<F; ++f, ++Y) { *Y = logf(*Y); }
                }
                else { Y += F; }
            }
        }
        fftwf_destroy_plan(plan); fftwf_free(Xw); fftwf_free(Yw);
    }

    return 0;
}


int stft_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t nfft, const size_t stp, const int snip_edges, const int mn0, const int amp, const int lg)
{
    if (L<1u) { fprintf(stderr,"error in stft_d: L must be positive\n"); return 1; }
    if (stp<1u) { fprintf(stderr,"error in stft_d: stp must be positive\n"); return 1; }
    if (snip_edges && L>N) { fprintf(stderr,"error in stft_d: L must be < N if snip_edges\n"); return 1; }

    //Set number of frames (W)
    const size_t W = (snip_edges) ? 1u+(N-L)/stp : (N+stp/2u)/stp;

    if (W==0u) {}
    else
    {
        const size_t F = nfft/2u + 1u;          //Num non-negative FFT freqs
        const int xd = (int)L - (int)stp;       //X1 inc after each frame
        double mn;                              //Mean of Xw (one window of X)

        //Initialize FFT
        double *Xw, *Yw;
        Xw = (double *)fftw_malloc(nfft*sizeof(double));
        Yw = (double *)fftw_malloc(nfft*sizeof(double));
        fftw_plan plan = fftw_plan_r2r_1d((int)nfft,Xw,Yw,FFTW_R2HC,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in stft_d: problem creating fftw plan"); return 1; }
        for (size_t l=0u; l<nfft; ++l, ++Xw) { *Xw = 0.0; }
        Xw -= nfft;

        if (snip_edges)
        {
            //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
            for (size_t w=0u; w<W; ++w)
            {
                //Window
                for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                X1 -= xd; X2 -= L;
                
                //Zero mean
                if (mn0)
                {
                    mn = 0.0;
                    for (size_t l=0u; l<L; ++l) { mn += *--Xw; }
                    mn /= (double)L;
                    for (size_t l=0u; l<L; ++l, ++Xw) { *Xw -= mn; }
                }
                
                //FFT
                Xw -= L;
                fftw_execute(plan);
                
                //Power
                for (size_t f=0u; f<F; ++f, ++Yw, ++Y) { *Y = *Yw * *Yw; }
                Y -= 2u;
                for (size_t f=1u; f<F-1u; ++f, ++Yw, --Y) { *Y += *Yw * *Yw; }
                Yw -= nfft;
                // for (size_t n=0u; n<nfft; ++n, ++Yw) { *Yw *= *Yw; }
                // for (size_t f=0u; f<F-1u; ++f) { *++Y = *--Yw; }
                // *Y-- = *Yw--;
                // for (size_t f=1u; f<F-1u; ++f, --Yw, --Y) { *Y += *Yw; }
                // *Y = *Yw;

                //Amplitude and/or log
                if (amp)
                {
                    if (lg) { for (size_t f=0u; f<F; ++f, ++Y) { *Y = log(sqrt(*Y)); } }
                    else { for (size_t f=0u; f<F; ++f, ++Y) { *Y = sqrt(*Y); } }
                }
                else if (lg)
                {
                    for (size_t f=0u; f<F; ++f, ++Y) { *Y = log(*Y); }
                }
                else { Y += F; }
            }
            //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
        }
        else
        {
            const size_t Lpre = L/2u;                   //nsamps before center samp
            int ss = (int)(stp/2u) - (int)Lpre;         //start-samp of current frame
            int n, prev_n = 0;                          //current/prev samps in X

            for (size_t w=0u; w<W; ++w)
            {
                //Window
                if (ss<0 || ss>(int)N-(int)L)
                {
                    for (int s=ss; s<ss+(int)L; ++s, ++X2, ++Xw)
                    {
                        //Window
                        n = s; //This ensures extrapolation by signal reversal to any length
                        while (n<0 || n>=(int)N) { n = (n<0) ? -n-1 : (n<(int)N) ? n : 2*(int)N-1-n; }
                        X1 += n - prev_n; prev_n = n;
                        *Xw = *X1 * *X2;
                    }
                }
                else
                {
                    X1 += ss - prev_n;
                    for (size_t l=0u; l<L; ++l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    X1 -= xd; prev_n = ss + (int)stp;
                }
                ss += stp; X2 -= L;

                //Zero mean
                if (mn0)
                {
                    mn = 0.0;
                    for (size_t l=0u; l<L; ++l) { mn += *--Xw; }
                    mn /= (double)L;
                    for (size_t l=0u; l<L; ++l, ++Xw) { *Xw -= mn; }
                }

                //FFT
                Xw -= L;
                fftw_execute(plan);
                
                //Power
                for (size_t f=0u; f<F; ++f, ++Yw, ++Y) { *Y = *Yw * *Yw; }
                Y -= 2u;
                for (size_t f=1u; f<F-1u; ++f, ++Yw, --Y) { *Y += *Yw * *Yw; }
                Yw -= nfft;

                //Amplitude and/or log
                if (amp)
                {
                    if (lg) { for (size_t f=0u; f<F; ++f, ++Y) { *Y = log(sqrt(*Y)); } }
                    else { for (size_t f=0u; f<F; ++f, ++Y) { *Y = sqrt(*Y); } }
                }
                else if (lg)
                {
                    for (size_t f=0u; f<F; ++f, ++Y) { *Y = log(*Y); }
                }
                else { Y += F; }
            }
        }
        fftw_destroy_plan(plan); fftw_free(Xw); fftw_free(Yw);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
