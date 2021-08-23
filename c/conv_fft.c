//1D convolution of each vector in X1f by X2 along dim.
//Each vector in X1f has length L1. X2 has length L2.

//This version uses the FFT for computation.

//FIR filtering is similar, except FIR is causal and conv is non-causal.
//Note that some "convolution" functions actually do cross-correlation.
//For actual cross-correlation (no flip of X2), see xcorr_fft.

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int conv_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_fft_d (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_fft_c (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);
int conv_fft_z (double *Y, const double *X1, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim);


int conv_fft_s (float *Y, const float *X1, const float *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const char shape[], const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_fft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L1 = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (L1<1u) { fprintf(stderr,"error in conv_fft_s: L1 (length of vecs in X1) must be positive\n"); return 1; }
    if (L2<1u) { fprintf(stderr,"error in conv_fft_s: L2 (length X2) must be positive\n"); return 1; }

    //Set ss, es, Ly according to shape
    size_t Ly;
    int ss, es;
    if (strncmp(shape,"full",4u)==0)
    {
        es = 0; ss = 1 - (int)L2;
        Ly = L1 + L2 - 1u;
    }
    else if (strncmp(shape,"same",4u)==0)
    {
        es = (int)L2/2; ss = es - (int)L2 + 1;
        Ly = L1;
    }
    else if (strncmp(shape,"valid",5u)==0)
    {
        ss = 0; es = (int)L2 - 1;
        Ly = (L1<L2) ? 0u : L1-L2+1u;
    }
    else
    {
        fprintf(stderr,"error in conv_fft_s: shape string must be 'full', 'same' or 'valid'\n"); return 1;
    }

    if (N==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while (nfft<L1+L2) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;

        //Set scale
        const float sc = 2.0f / (float)nfft;
        const size_t isodd = nfft%2u;

        //Initialize fftw
        float *X1f, *X2f, *Y1f, *Y2f, *Z1f, yr, yi, dr, di;
        X1f = (float *)fftwf_malloc(nfft*sizeof(float));
        X2f = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Z1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan xplan = fftwf_plan_dft_r2c_1d((int)nfft,X1f,(fftwf_complex *)Y1f,FFTW_ESTIMATE);
        fftwf_plan bplan = fftwf_plan_dft_r2c_1d((int)nfft,X2f,(fftwf_complex *)Y2f,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in conv_fft_s: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in conv_fft_s: problem creating fftw plan"); return 1; }
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1f,(fftwf_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv_fft_s: problem creating fftw plan"); return 1; }
        for (size_t n=L1; n<nfft; ++n) { X1f[n] = 0.0f; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1f[n] = 0.0f; }

        //Get Y2f (FFT of X2f, then scale)
        for (size_t n=0u; n<=L2; ++n, ++X2, ++X2f) { *X2f = *X2; }
        for (size_t n=L2+1u; n<nfft; ++n, ++X2f) { *X2f = 0.0f; }
        X2f -= nfft;
        fftwf_execute(bplan);
        *Y2f++ /= (float)nfft; ++Y2f;
        for (size_t l=nfft+isodd-2u; l>0u; --l, ++Y2f) { *Y2f *= sc; }
        if (!isodd) { *Y2f++ /= (float)nfft; }
        Y2f -= nfft + 1u;

        if (L1==N)
        {
            for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
            X1f -= L1;
            fftwf_execute(xplan);
            *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
            for (size_t f=F-2u; f>0u; --f)
            {
                yr = *Y1f++; yi = *Y1f--;
                dr = *Y2f++; di = *Y2f++;
                *Y1f++ = dr*yr - di*yi;
                *Y1f++ = dr*yi + di*yr;
            }
            *Y1f *= *Y2f;
            Y2f -= 2u*F-2u; Y1f -= 2u*F-2u;
            fftwf_execute(iplan);
            for (size_t l=L1; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
            Z1f -= 2u*L1;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L1, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    for (size_t l=L1; l>0u; --l, ++X1, ++X1f) { *X1f = *X1; }
                    X1f -= L1;
                    fftwf_execute(xplan);
                    *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                    for (size_t f=F-2u; f>0u; --f)
                    {
                        yr = *Y1f++; yi = *Y1f--;
                        dr = *Y2f++; di = *Y2f++;
                        *Y1f++ = dr*yr - di*yi;
                        *Y1f++ = dr*yi + di*yr;
                    }
                    *Y1f *= *Y2f;
                    Y2f -= 2u*F-2u; Y1f -= 2u*F-2u;
                    fftwf_execute(iplan);
                    for (size_t l=L1; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
                    Z1f -= 2u*L1;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X1+=B*(L1-1u), Y+=B*(Ly-1u))
                {
                    for (size_t bs=0u; bs<B; ++bs, X1-=K*L1-1u, Y-=K*Ly-1u)
                    {
                        for (size_t l=L1; l>0u; --l, X1+=K, ++X1f) { *X1f = *X1; }
                        X1f -= L1;
                        fftwf_execute(xplan);
                        *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                        for (size_t f=F-2u; f>0u; --f)
                        {
                            yr = *Y1f++; yi = *Y1f--;
                            dr = *Y2f++; di = *Y2f++;
                            *Y1f++ = dr*yr - di*yi;
                            *Y1f++ = dr*yi + di*yr;
                        }
                        *Y1f *= *Y2f;
                        Y2f -= 2u*F-2u; Y1f -= 2u*F-2u;
                        fftwf_execute(iplan);
                        for (size_t l=Ly; l>0u; --l, Z1f+=2, Y+=K) { *Y = *Z1f; }
                        Z1f -= 2u*Ly;
                    }
                }
            }
        }
        fftwf_free(X1f); fftwf_free(Y2f); fftwf_free(Y1f); fftwf_free(Z1f);
        fftwf_destroy_plan(xplan); fftwf_destroy_plan(bplan); fftwf_destroy_plan(iplan);
    }

    return 0;
}


int conv_fft_d (double *Y, const double *X, const double *X2, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_fft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while (nfft<L+L2) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;

        //Set scale
        const double sc = 2.0 / (double)nfft;
        const size_t isodd = nfft%2u;

        //Initialize fftw
        double *X1f, *X2f, *Y1f, *Y2f, *Z1f, yr, yi, dr, di;
        X1f = (double *)fftw_malloc(nfft*sizeof(double));
        X2f = (double *)fftw_malloc(nfft*sizeof(double));
        Y1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Z1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan xplan = fftw_plan_dft_r2c_1d((int)nfft,X1f,(fftw_complex *)Y1f,FFTW_ESTIMATE);
        fftw_plan bplan = fftw_plan_dft_r2c_1d((int)nfft,X2f,(fftw_complex *)Y2f,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in conv_fft_d: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in conv_fft_d: problem creating fftw plan"); return 1; }
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1f,(fftw_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv_fft_d: problem creating fftw plan"); return 1; }
        for (size_t n=L; n<nfft; ++n) { X1f[n] = 0.0; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1f[n] = 0.0; }

        //Get Y2f (FFT of X2f, then scale)
        for (size_t n=0u; n<=L2; ++n, ++B, ++X2f) { *X2f = *B; }
        for (size_t n=L2+1u; n<nfft; ++n, ++X2f) { *X2f = 0.0; }
        X2f -= nfft;
        fftw_execute(bplan);
        *Y2f++ /= (double)nfft; ++Y2f;
        for (size_t l=nfft+isodd-2u; l>0u; --l, ++Y2f) { *Y2f *= sc; }
        if (!isodd) { *Y2f++ /= (double)nfft; }
        Y2f -= nfft + 1u;

        if (L==N)
        {
            for (size_t l=L; l>0u; --l, ++X, ++X1f) { *X1f = *X; }
            X1f -= L;
            fftw_execute(xplan);
            *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
            for (size_t f=F-2u; f>0u; --f)
            {
                yr = *Y1f++; yi = *Y1f--;
                dr = *Y2f++; di = *Y2f++;
                *Y1f++ = dr*yr - di*yi;
                *Y1f++ = dr*yi + di*yr;
            }
            *Y1f *= *Y2f;
            Y2f -= 2u*F-2u; Y1f -= 2u*F-2u;
            fftw_execute(iplan);
            for (size_t l=L; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
            Z1f -= 2u*L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/BS;

            if (K==1u && (G==1u || BS==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    for (size_t l=L; l>0u; --l, ++X, ++X1f) { *X1f = *X; }
                    X1f -= L;
                    fftw_execute(xplan);
                    *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                    for (size_t f=F-2u; f>0u; --f)
                    {
                        yr = *Y1f++; yi = *Y1f--;
                        dr = *Y2f++; di = *Y2f++;
                        *Y1f++ = dr*yr - di*yi;
                        *Y1f++ = dr*yi + di*yr;
                    }
                    *Y1f *= *Y2f;
                    Y2f -= 2u*F-2u; Y1f -= 2u*F-2u;
                    fftw_execute(iplan);
                    for (size_t l=L; l>0u; --l, Z1f+=2, ++Y) { *Y = *Z1f; }
                    Z1f -= 2u*L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=BS*(L-1u), Y+=BS*(L-1u))
                {
                    for (size_t bs=BS; bs>0u; --bs, X-=K*L-1u, Y-=K*L-1u)
                    {
                        for (size_t l=L; l>0u; --l, X+=K, ++X1f) { *X1f = *X; }
                        X1f -= L;
                        fftw_execute(xplan);
                        *Y1f++ *= *Y2f++; ++Y1f; ++Y2f;
                        for (size_t f=F-2u; f>0u; --f)
                        {
                            yr = *Y1f++; yi = *Y1f--;
                            dr = *Y2f++; di = *Y2f++;
                            *Y1f++ = dr*yr - di*yi;
                            *Y1f++ = dr*yi + di*yr;
                        }
                        *Y1f *= *Y2f;
                        Y2f -= 2u*F-2u; Y1f -= 2u*F-2u;
                        fftw_execute(iplan);
                        for (size_t l=L; l>0u; --l, Z1f+=2, Y+=K) { *Y = *Z1f; }
                        Z1f -= 2u*L;
                    }
                }
            }
        }
        fftw_free(X1f); fftw_free(Y2f); fftw_free(Y1f); fftw_free(Z1f);
        fftw_destroy_plan(xplan); fftw_destroy_plan(bplan); fftw_destroy_plan(iplan);
    }

    return 0;
}


int conv_fft_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_fft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        //Set nfft
        size_t nfft = 1u;
        while (nfft<L+L2) { nfft *= 2u; }

        //Initialize fftw
        float *X1f, *X2f, *Y1f, *Y2f, *Z1f, yr, yi, dr, di;
        X1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        X2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y2f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Z1f = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan xplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1f,(fftwf_complex *)Y1f,FFTW_FORWARD,FFTW_ESTIMATE);
        fftwf_plan bplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X2f,(fftwf_complex *)Y2f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in conv_fft_c: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in conv_fft_c: problem creating fftw plan"); return 1; }
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1f,(fftwf_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv_fft_c: problem creating fftw plan"); return 1; }
        for (size_t n=2u*L; n<2u*nfft; ++n) { X1f[n] = 0.0f; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1f[n] = 0.0f; }

        //Get Y2f (FFT of X2f, then scale)
        for (size_t n=0u; n<=L2; ++n, ++B, ++X2f) { *X2f = *B; *++X2f = *++B; }
        for (size_t n=L2+1u; n<nfft; ++n, ++X2f) { *X2f = 0.0f; *++X2f = 0.0f; }
        X2f -= 2u*nfft;
        fftwf_execute(bplan);
        for (size_t l=0u; l<2u*nfft; ++l, ++Y2f) { *Y2f /= (float)nfft; }
        Y2f -= 2u*nfft;

        if (L==N)
        {
            for (size_t l=0u; l<2u*L; ++l, ++X, ++X1f) { *X1f = *X; }
            X1f -= 2u*L;
            fftwf_execute(xplan);
            for (size_t n=0u; n<nfft; ++n)
            {
                yr = *Y1f++; yi = *Y1f--;
                dr = *Y2f++; di = *Y2f++;
                *Y1f++ = dr*yr - di*yi;
                *Y1f++ = dr*yi + di*yr;
            }
            Y2f -= 2u*nfft; Y1f -= 2u*nfft;
            fftwf_execute(iplan);
            for (size_t l=0u; l<2u*L; ++l, ++Z1f, ++Y) { *Y = *Z1f; }
            Z1f -= 2u*L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/BS;

            if (K==1u && (G==1u || BS==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    for (size_t l=0u; l<2u*L; ++l, ++X, ++X1f) { *X1f = *X; }
                    X1f -= 2u*L;
                    fftwf_execute(xplan);
                    for (size_t n=0u; n<nfft; ++n)
                    {
                        yr = *Y1f++; yi = *Y1f--;
                        dr = *Y2f++; di = *Y2f++;
                        *Y1f++ = dr*yr - di*yi;
                        *Y1f++ = dr*yi + di*yr;
                    }
                    Y2f -= 2u*nfft; Y1f -= 2u*nfft;
                    fftwf_execute(iplan);
                    for (size_t l=0u; l<2u*L; ++l, ++Z1f, ++Y) { *Y = *Z1f; }
                    Z1f -= 2u*L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
                {
                    for (size_t bs=BS; bs>0u; --bs, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                    {
                        for (size_t l=L; l>0u; --l, X+=2u*K, ++X1f) { *X1f = *X; *++X1f = *(X+1); }
                        X1f -= 2u*L;
                        fftwf_execute(xplan);
                        for (size_t n=0u; n<nfft; ++n)
                        {
                            yr = *Y1f++; yi = *Y1f--;
                            dr = *Y2f++; di = *Y2f++;
                            *Y1f++ = dr*yr - di*yi;
                            *Y1f++ = dr*yi + di*yr;
                        }
                        Y2f -= 2u*nfft; Y1f -= 2u*nfft;
                        fftwf_execute(iplan);
                        for (size_t l=L; l>0u; --l, ++Z1f, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *++Z1f; }
                        Z1f -= 2u*L;
                    }
                }
            }
        }
        fftwf_free(X1f); fftwf_free(Y2f); fftwf_free(Y1f); fftwf_free(Z1f);
        fftwf_destroy_plan(xplan); fftwf_destroy_plan(bplan); fftwf_destroy_plan(iplan);
    }

    return 0;
}


int conv_fft_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t L2, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in conv_fft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        //Set nfft
        size_t nfft = 1u;
        while (nfft<L+L2) { nfft *= 2u; }

        //Initialize fftw
        double *X1f, *X2f, *Y1f, *Y2f, *Z1f, yr, yi, dr, di;
        X1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        X2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y2f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Z1f = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan xplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1f,(fftw_complex *)Y1f,FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_plan bplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X2f,(fftw_complex *)Y2f,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in conv_fft_z: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in conv_fft_z: problem creating fftw plan"); return 1; }
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1f,(fftw_complex *)Z1f,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in conv_fft_z: problem creating fftw plan"); return 1; }
        for (size_t n=2u*L; n<2u*nfft; ++n) { X1f[n] = 0.0; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1f[n] = 0.0; }

        //Get Y2f (FFT of X2f, then scale)
        for (size_t n=0u; n<=L2; ++n, ++B, ++X2f) { *X2f = *B; *++X2f = *++B; }
        for (size_t n=L2+1u; n<nfft; ++n, ++X2f) { *X2f = 0.0; *++X2f = 0.0; }
        X2f -= 2u*nfft;
        fftw_execute(bplan);
        for (size_t l=0u; l<2u*nfft; ++l, ++Y2f) { *Y2f /= (double)nfft; }
        Y2f -= 2u*nfft;

        if (L==N)
        {
            for (size_t l=0u; l<2u*L; ++l, ++X, ++X1f) { *X1f = *X; }
            X1f -= 2u*L;
            fftw_execute(xplan);
            for (size_t n=0u; n<nfft; ++n)
            {
                yr = *Y1f++; yi = *Y1f--;
                dr = *Y2f++; di = *Y2f++;
                *Y1f++ = dr*yr - di*yi;
                *Y1f++ = dr*yi + di*yr;
            }
            Y2f -= 2u*nfft; Y1f -= 2u*nfft;
            fftw_execute(iplan);
            for (size_t l=0u; l<2u*L; ++l, ++Z1f, ++Y) { *Y = *Z1f; }
            Z1f -= 2u*L;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t BS = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/L, G = V/BS;

            if (K==1u && (G==1u || BS==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    for (size_t l=0u; l<2u*L; ++l, ++X, ++X1f) { *X1f = *X; }
                    X1f -= 2u*L;
                    fftw_execute(xplan);
                    for (size_t n=0u; n<nfft; ++n)
                    {
                        yr = *Y1f++; yi = *Y1f--;
                        dr = *Y2f++; di = *Y2f++;
                        *Y1f++ = dr*yr - di*yi;
                        *Y1f++ = dr*yi + di*yr;
                    }
                    Y2f -= 2u*nfft; Y1f -= 2u*nfft;
                    fftw_execute(iplan);
                    for (size_t l=0u; l<2u*L; ++l, ++Z1f, ++Y) { *Y = *Z1f; }
                    Z1f -= 2u*L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
                {
                    for (size_t bs=BS; bs>0u; --bs, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                    {
                        for (size_t l=L; l>0u; --l, X+=2u*K, ++X1f) { *X1f = *X; *++X1f = *(X+1); }
                        X1f -= 2u*L;
                        fftw_execute(xplan);
                        for (size_t n=0u; n<nfft; ++n)
                        {
                            yr = *Y1f++; yi = *Y1f--;
                            dr = *Y2f++; di = *Y2f++;
                            *Y1f++ = dr*yr - di*yi;
                            *Y1f++ = dr*yi + di*yr;
                        }
                        Y2f -= 2u*nfft; Y1f -= 2u*nfft;
                        fftw_execute(iplan);
                        for (size_t l=L; l>0u; --l, ++Z1f, Y+=2u*K) { *Y = *Z1f; *(Y+1) = *++Z1f; }
                        Z1f -= 2u*L;
                    }
                }
            }
        }
        fftw_free(X1f); fftw_free(Y2f); fftw_free(Y1f); fftw_free(Z1f);
        fftw_destroy_plan(xplan); fftw_destroy_plan(bplan); fftw_destroy_plan(iplan);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
