//STFT (short-term Fourier transform) of univariate time series X1,
//using window X2, and outputing power in Y.

//This takes a continuous univariate time series (vector X1) of length N,
//and a window (vector X2) of length L,
//and outputs the STFT (matrix Y) at each of W windows.
//If Y is row-major, then it has size W x F.
//If Y is col-major, then it has size F x W.
//where F is nfft/2+1, and nfft is the next-pow-2 of L.

//This uses the framing/windowing conventions of frame_univar_flt and window_univar_flt.

//The following boolean (int) options are added:
//mn0: subtract mean from each frame just after windowing.
//amp: take sqrt of each element of Y just after getting power.
//lg:  take log of each element of Y just before output.

#include <stdio.h>
#include <float.h>
#include <math.h>
#include <fftw3.h>
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int stft_flt_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t W, const size_t nfft, const float c0, const float stp, const int mn0, const int amp, const int lg);
int stft_flt_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t W, const size_t nfft, const double c0, const double stp, const int mn0, const int amp, const int lg);


int stft_flt_s (float *Y, const float *X1, const float *X2, const size_t N, const size_t L, const size_t W, const size_t nfft, const float c0, const float stp, const int mn0, const int amp, const int lg)
{
    if (W>N) { fprintf(stderr,"error in stft_flt_s: W must be <= N (length X)\n"); return 1; }
    if (L>N) { fprintf(stderr,"error in stft_flt_s: L must be <= N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in stft_flt_s: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in stft_flt_s: N (length X) must be positive\n"); return 1; }
    if (c0>(float)(N-1u)) { fprintf(stderr,"error in stft_flt_s: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<FLT_EPSILON) { fprintf(stderr,"error in stft_flt_s: stp (step size) must be positive\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t F = nfft/2u + 1u;          //Num non-negative FFT freqs
        float mn;                               //Mean of Xw (one window of X)
        
        //Initialize FFT
        float *Xw, *Yw;
        Xw = (float *)fftwf_malloc(nfft*sizeof(float));
        Yw = (float *)fftwf_malloc(nfft*sizeof(float));
        fftwf_plan plan = fftwf_plan_r2r_1d((int)nfft,Xw,Yw,FFTW_R2HC,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in stft_flt_s: problem creating fftw plan"); return 1; }
        for (size_t n=nfft; n>0u; --n, ++Xw) { *Xw = 0.0f; }
        Xw -= nfft;

        const size_t Lpre = L/2u;                   //nsamps before center samp
        const size_t Lpost = L-Lpre-1u;             //nsamps after center samp
        float cc = c0;                              //current exact center-samp
        int cs = (int)roundf(c0);                   //current rounded center-samp
        int ss = cs-(int)Lpre, es = cs+(int)Lpost;  //current rounded start-samp, end-samp

        for (size_t w=W; w>0u; --w)
        {
            //Window
            if (es<0 || ss>(int)N-1)
            {
                for (size_t l=L; l>0u; --l, ++Xw) { *Xw = 0.0f; }
            }
            else
            {
                if (ss<0)
                {
                    for (size_t l=(size_t)(-ss); l>0u; --l, ++X2, ++Xw) { *Xw = 0.0f; }
                    for (size_t l=(size_t)(-ss); l<L; ++l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    if (cs<(int)Lpre) { X1 -= (int)L + ss; } else { X1 -= (int)L; }
                }
                else if (es<(int)N)
                {
                    for (size_t l=L; l>0u; --l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    X1 += (int)roundf(cc+stp) - cs - (int)L;
                }
                else //if (ss<(int)N)
                {
                    for (size_t l=N-(size_t)ss; l>0u; --l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    for (size_t l=N-(size_t)ss; l<L; ++l, ++X2, ++Xw) { *Xw = 0.0f; }
                    X1 += (int)roundf(cc+stp) - cs - (int)N + ss;
                }
                X2 -= L;
            }
            cc += stp; cs = (int)roundf(cc);
            ss = cs - (int)Lpre; es = cs + (int)Lpost;

            //Zero mean
            if (mn0)
            {
                mn = 0.0f;
                for (size_t l=L; l>0u; --l) { mn += *--Xw; }
                mn /= (float)L;
                for (size_t l=L; l>0u; --l, ++Xw) { *Xw -= mn; }
            }

            //FFT
            Xw -= L;
            fftwf_execute(plan);
            
            //Power
            for (size_t f=F; f>0u; --f, ++Yw, ++Y) { *Y = *Yw * *Yw; }
            Y -= 2u;
            for (size_t f=F-2u; f>0u; --f, ++Yw, --Y) { *Y += *Yw * *Yw; }
            Yw -= nfft;
            
            //Amplitude and/or log
            if (amp)
            {
                if (lg) { for (size_t f=F; f>0u; --f, ++Y) { *Y = logf(sqrtf(*Y)); } }
                else { for (size_t f=F; f>0u; --f, ++Y) { *Y = sqrtf(*Y); } }
            }
            else if (lg)
            {
                for (size_t f=F; f>0u; --f, ++Y) { *Y = logf(*Y); }
            }
            else { Y += F; }
        }
        fftwf_destroy_plan(plan); fftwf_free(Xw); fftwf_free(Yw);
    }

    return 0;
}


int stft_flt_d (double *Y, const double *X1, const double *X2, const size_t N, const size_t L, const size_t W, const size_t nfft, const double c0, const double stp, const int mn0, const int amp, const int lg)
{
    if (W>N) { fprintf(stderr,"error in stft_flt_d: W must be <= N (length X)\n"); return 1; }
    if (L>N) { fprintf(stderr,"error in stft_flt_d: L must be <= N (length X)\n"); return 1; }
    if (L<1u) { fprintf(stderr,"error in stft_flt_d: L must be positive\n"); return 1; }
    if (N<1u) { fprintf(stderr,"error in stft_flt_d: N (length X) must be positive\n"); return 1; }
    if (c0>(double)(N-1u)) { fprintf(stderr,"error in stft_flt_d: c0 (center samp of 1st frame) must be < N (length X)\n"); return 1; }
    if (stp<(double)FLT_EPSILON) { fprintf(stderr,"error in stft_flt_d: stp (step size) must be positive\n"); return 1; }

    if (W==0u) {}
    else
    {
        const size_t F = nfft/2u + 1u;          //Num non-negative FFT freqs
        double mn;                              //Mean of Xw (one window of X)

        //Initialize FFT
        double *Xw, *Yw;
        Xw = (double *)fftw_malloc(nfft*sizeof(double));
        Yw = (double *)fftw_malloc(nfft*sizeof(double));
        fftw_plan plan = fftw_plan_r2r_1d((int)nfft,Xw,Yw,FFTW_R2HC,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in stft_flt_d: problem creating fftw plan"); return 1; }
        for (size_t n=nfft; n>0u; --n, ++Xw) { *Xw = 0.0; }
        Xw -= nfft;

        const size_t Lpre = L/2u;                   //nsamps before center samp
        const size_t Lpost = L-Lpre-1u;             //nsamps after center samp
        double cc = c0;                             //current exact center-samp
        int cs = (int)round(c0);                    //current rounded center-samp
        int ss = cs-(int)Lpre, es = cs+(int)Lpost;  //current rounded start-samp, end-samp

        for (size_t w=W; w>0u; --w)
        {
            //Window
            if (es<0 || ss>(int)N-1)
            {
                for (size_t l=L; l>0u; --l, ++Xw) { *Xw = 0.0; }
            }
            else
            {
                if (ss<0)
                {
                    for (size_t l=(size_t)(-ss); l>0u; --l, ++X2, ++Xw) { *Xw = 0.0; }
                    for (size_t l=(size_t)(-ss); l<L; ++l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    if (cs<(int)Lpre) { X1 -= (int)L + ss; } else { X1 -= (int)L; }
                }
                else if (es<(int)N)
                {
                    for (size_t l=L; l>0u; --l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    X1 += (int)round(cc+stp) - cs - (int)L;
                }
                else //if (ss<(int)N)
                {
                    for (size_t l=N-(size_t)ss; l>0u; --l, ++X1, ++X2, ++Xw) { *Xw = *X1 * *X2; }
                    for (size_t l=N-(size_t)ss; l<L; ++l, ++X2, ++Xw) { *Xw = 0.0; }
                    X1 += (int)round(cc+stp) - cs - (int)N + ss;
                }
                X2 -= L;
            }
            cc += stp; cs = (int)round(cc);
            ss = cs - (int)Lpre; es = cs + (int)Lpost;

            //Zero mean
            if (mn0)
            {
                mn = 0.0;
                for (size_t l=L; l>0u; --l) { mn += *--Xw; }
                mn /= (double)L;
                for (size_t l=L; l>0u; --l, ++Xw) { *Xw -= mn; }
            }

            //FFT
            Xw -= L;
            fftw_execute(plan);
            
            //Power
            for (size_t f=F; f>0u; --f, ++Yw, ++Y) { *Y = *Yw * *Yw; }
            Y -= 2u;
            for (size_t f=F-2u; f>0u; --f, ++Yw, --Y) { *Y += *Yw * *Yw; }
            Yw -= nfft;

            //Amplitude and/or log
            if (amp)
            {
                if (lg) { for (size_t f=F; f>0u; --f, ++Y) { *Y = log(sqrt(*Y)); } }
                else { for (size_t f=F; f>0u; --f, ++Y) { *Y = sqrt(*Y); } }
            }
            else if (lg)
            {
                for (size_t f=F; f>0u; --f, ++Y) { *Y = log(*Y); }
            }
            else { Y += F; }
        }
        fftw_destroy_plan(plan); fftw_free(Xw); fftw_free(Yw);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
