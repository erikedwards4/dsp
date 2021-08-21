//Causal FIR filtering of each vector in X along dim.
//FIR impulse response is given in vector B with length Q+1.
//(I use P for IIR filter order, since same as polynomial order.)

//This uses essentially: Y = ifft(fft(X,nfft).*fft(B,nfft),nfft).
//Note that this is NOT a OLA (overlap and add) method, just one big FFT.

//nfft is set to nextpow2 of L+Q

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fir_fft_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_fft_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_fft_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);
int fir_fft_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim);


int fir_fft_s (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_fft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while (nfft<L+Q) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;

        //Set scale
        const float sc = 2.0f / (float)nfft;
        const size_t isodd = nfft%2u;

        //Initialize fftw
        float *X1, *B1, *Y1, *D1, *Z1, yr, yi, dr, di;
        X1 = (float *)fftwf_malloc(nfft*sizeof(float));
        B1 = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        D1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Z1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan xplan = fftwf_plan_dft_r2c_1d((int)nfft,X1,(fftwf_complex *)Y1,FFTW_ESTIMATE);
        fftwf_plan bplan = fftwf_plan_dft_r2c_1d((int)nfft,B1,(fftwf_complex *)D1,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in fir_fft_s: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in fir_fft_s: problem creating fftw plan"); return 1; }
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1,(fftwf_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in fir_fft_s: problem creating fftw plan"); return 1; }
        for (size_t n=L; n<nfft; ++n) { X1[n] = 0.0f; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1[n] = 0.0f; }

        //Get D1 (FFT of B1, then scale)
        for (size_t n=0u; n<=Q; ++n, ++B, ++B1) { *B1 = *B; }
        for (size_t n=Q+1u; n<nfft; ++n, ++B1) { *B1 = 0.0f; }
        B1 -= nfft;
        fftwf_execute(bplan);
        *D1++ /= (float)nfft; ++D1;
        for (size_t l=2u; l<nfft+isodd; ++l, ++D1) { *D1 *= sc; }
        if (!isodd) { *D1++ /= (float)nfft; }
        D1 -= nfft + 1u;

        if (L==N)
        {
            for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= L;
            fftwf_execute(xplan);
            *Y1++ *= *D1++; ++Y1; ++D1;
            for (size_t f=1u; f<F-1u; ++f)
            {
                yr = *Y1++; yi = *Y1--;
                dr = *D1++; di = *D1++;
                *Y1++ = dr*yr - di*yi;
                *Y1++ = dr*yi + di*yr;
            }
            *Y1 *= *D1;
            D1 -= 2u*F-2u; Y1 -= 2u*F-2u;
            fftwf_execute(iplan);
            for (size_t l=0u; l<L; ++l, Z1+=2, ++Y) { *Y = *Z1; }
            Z1 -= 2u*L;
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
                    for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftwf_execute(xplan);
                    *Y1++ *= *D1++; ++Y1; ++D1;
                    for (size_t f=1u; f<F-1u; ++f)
                    {
                        yr = *Y1++; yi = *Y1--;
                        dr = *D1++; di = *D1++;
                        *Y1++ = dr*yr - di*yi;
                        *Y1++ = dr*yi + di*yr;
                    }
                    *Y1 *= *D1;
                    D1 -= 2u*F-2u; Y1 -= 2u*F-2u;
                    fftwf_execute(iplan);
                    for (size_t l=0u; l<L; ++l, Z1+=2, ++Y) { *Y = *Z1; }
                    Z1 -= 2u*L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=BS*(L-1u), Y+=BS*(L-1u))
                {
                    for (size_t bs=0u; bs<BS; ++bs, X-=K*L-1u, Y-=K*L-1u)
                    {
                        for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        fftwf_execute(xplan);
                        *Y1++ *= *D1++; ++Y1; ++D1;
                        for (size_t f=1u; f<F-1u; ++f)
                        {
                            yr = *Y1++; yi = *Y1--;
                            dr = *D1++; di = *D1++;
                            *Y1++ = dr*yr - di*yi;
                            *Y1++ = dr*yi + di*yr;
                        }
                        *Y1 *= *D1;
                        D1 -= 2u*F-2u; Y1 -= 2u*F-2u;
                        fftwf_execute(iplan);
                        for (size_t l=0u; l<L; ++l, Z1+=2, Y+=K) { *Y = *Z1; }
                        Z1 -= 2u*L;
                    }
                }
            }
        }
        fftwf_free(X1); fftwf_free(D1); fftwf_free(Y1); fftwf_free(Z1);
        fftwf_destroy_plan(xplan); fftwf_destroy_plan(bplan); fftwf_destroy_plan(iplan);
    }

    return 0;
}


int fir_fft_d (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_fft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        //Set nfft and F
        size_t nfft = 1u;
        while (nfft<L+Q) { nfft *= 2u; }
        const size_t F = nfft/2u + 1u;

        //Set scale
        const double sc = 2.0 / (double)nfft;
        const size_t isodd = nfft%2u;

        //Initialize fftw
        double *X1, *B1, *Y1, *D1, *Z1, yr, yi, dr, di;
        X1 = (double *)fftw_malloc(nfft*sizeof(double));
        B1 = (double *)fftw_malloc(nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        D1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Z1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan xplan = fftw_plan_dft_r2c_1d((int)nfft,X1,(fftw_complex *)Y1,FFTW_ESTIMATE);
        fftw_plan bplan = fftw_plan_dft_r2c_1d((int)nfft,B1,(fftw_complex *)D1,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in fir_fft_d: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in fir_fft_d: problem creating fftw plan"); return 1; }
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1,(fftw_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in fir_fft_d: problem creating fftw plan"); return 1; }
        for (size_t n=L; n<nfft; ++n) { X1[n] = 0.0; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1[n] = 0.0; }

        //Get D1 (FFT of B1, then scale)
        for (size_t n=0u; n<=Q; ++n, ++B, ++B1) { *B1 = *B; }
        for (size_t n=Q+1u; n<nfft; ++n, ++B1) { *B1 = 0.0; }
        B1 -= nfft;
        fftw_execute(bplan);
        *D1++ /= (double)nfft; ++D1;
        for (size_t l=2u; l<nfft+isodd; ++l, ++D1) { *D1 *= sc; }
        if (!isodd) { *D1++ /= (double)nfft; }
        D1 -= nfft + 1u;

        if (L==N)
        {
            for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= L;
            fftw_execute(xplan);
            *Y1++ *= *D1++; ++Y1; ++D1;
            for (size_t f=1u; f<F-1u; ++f)
            {
                yr = *Y1++; yi = *Y1--;
                dr = *D1++; di = *D1++;
                *Y1++ = dr*yr - di*yi;
                *Y1++ = dr*yi + di*yr;
            }
            *Y1 *= *D1;
            D1 -= 2u*F-2u; Y1 -= 2u*F-2u;
            fftw_execute(iplan);
            for (size_t l=0u; l<L; ++l, Z1+=2, ++Y) { *Y = *Z1; }
            Z1 -= 2u*L;
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
                    for (size_t l=0u; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(xplan);
                    *Y1++ *= *D1++; ++Y1; ++D1;
                    for (size_t f=1u; f<F-1u; ++f)
                    {
                        yr = *Y1++; yi = *Y1--;
                        dr = *D1++; di = *D1++;
                        *Y1++ = dr*yr - di*yi;
                        *Y1++ = dr*yi + di*yr;
                    }
                    *Y1 *= *D1;
                    D1 -= 2u*F-2u; Y1 -= 2u*F-2u;
                    fftw_execute(iplan);
                    for (size_t l=0u; l<L; ++l, Z1+=2, ++Y) { *Y = *Z1; }
                    Z1 -= 2u*L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=BS*(L-1u), Y+=BS*(L-1u))
                {
                    for (size_t bs=0u; bs<BS; ++bs, X-=K*L-1u, Y-=K*L-1u)
                    {
                        for (size_t l=0u; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= L;
                        fftw_execute(xplan);
                        *Y1++ *= *D1++; ++Y1; ++D1;
                        for (size_t f=1u; f<F-1u; ++f)
                        {
                            yr = *Y1++; yi = *Y1--;
                            dr = *D1++; di = *D1++;
                            *Y1++ = dr*yr - di*yi;
                            *Y1++ = dr*yi + di*yr;
                        }
                        *Y1 *= *D1;
                        D1 -= 2u*F-2u; Y1 -= 2u*F-2u;
                        fftw_execute(iplan);
                        for (size_t l=0u; l<L; ++l, Z1+=2, Y+=K) { *Y = *Z1; }
                        Z1 -= 2u*L;
                    }
                }
            }
        }
        fftw_free(X1); fftw_free(D1); fftw_free(Y1); fftw_free(Z1);
        fftw_destroy_plan(xplan); fftw_destroy_plan(bplan); fftw_destroy_plan(iplan);
    }

    return 0;
}


int fir_fft_c (float *Y, const float *X, const float *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_fft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        //Set nfft
        size_t nfft = 1u;
        while (nfft<L+Q) { nfft *= 2u; }

        //Initialize fftw
        float *X1, *B1, *Y1, *D1, *Z1, yr, yi, dr, di;
        X1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        B1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        D1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Z1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan xplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        fftwf_plan bplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)B1,(fftwf_complex *)D1,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in fir_fft_c: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in fir_fft_c: problem creating fftw plan"); return 1; }
        fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1,(fftwf_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in fir_fft_c: problem creating fftw plan"); return 1; }
        for (size_t n=2u*L; n<2u*nfft; ++n) { X1[n] = 0.0f; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1[n] = 0.0f; }

        //Get D1 (FFT of B1, then scale)
        for (size_t n=0u; n<=Q; ++n, ++B, ++B1) { *B1 = *B; *++B1 = *++B; }
        for (size_t n=Q+1u; n<nfft; ++n, ++B1) { *B1 = 0.0f; *++B1 = 0.0f; }
        B1 -= 2u*nfft;
        fftwf_execute(bplan);
        for (size_t l=0u; l<2u*nfft; ++l, ++D1) { *D1 /= (float)nfft; }
        D1 -= 2u*nfft;

        if (L==N)
        {
            for (size_t l=0u; l<2u*L; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*L;
            fftwf_execute(xplan);
            for (size_t n=0u; n<nfft; ++n)
            {
                yr = *Y1++; yi = *Y1--;
                dr = *D1++; di = *D1++;
                *Y1++ = dr*yr - di*yi;
                *Y1++ = dr*yi + di*yr;
            }
            D1 -= 2u*nfft; Y1 -= 2u*nfft;
            fftwf_execute(iplan);
            for (size_t l=0u; l<2u*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
            Z1 -= 2u*L;
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
                    for (size_t l=0u; l<2u*L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*L;
                    fftwf_execute(xplan);
                    for (size_t n=0u; n<nfft; ++n)
                    {
                        yr = *Y1++; yi = *Y1--;
                        dr = *D1++; di = *D1++;
                        *Y1++ = dr*yr - di*yi;
                        *Y1++ = dr*yi + di*yr;
                    }
                    D1 -= 2u*nfft; Y1 -= 2u*nfft;
                    fftwf_execute(iplan);
                    for (size_t l=0u; l<2u*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
                    Z1 -= 2u*L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
                {
                    for (size_t bs=0u; bs<BS; ++bs, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                    {
                        for (size_t l=0u; l<L; ++l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                        X1 -= 2u*L;
                        fftwf_execute(xplan);
                        for (size_t n=0u; n<nfft; ++n)
                        {
                            yr = *Y1++; yi = *Y1--;
                            dr = *D1++; di = *D1++;
                            *Y1++ = dr*yr - di*yi;
                            *Y1++ = dr*yi + di*yr;
                        }
                        D1 -= 2u*nfft; Y1 -= 2u*nfft;
                        fftwf_execute(iplan);
                        for (size_t l=0u; l<L; ++l, ++Z1, Y+=2u*K) { *Y = *Z1; *(Y+1) = *++Z1; }
                        Z1 -= 2u*L;
                    }
                }
            }
        }
        fftwf_free(X1); fftwf_free(D1); fftwf_free(Y1); fftwf_free(Z1);
        fftwf_destroy_plan(xplan); fftwf_destroy_plan(bplan); fftwf_destroy_plan(iplan);
    }

    return 0;
}


int fir_fft_z (double *Y, const double *X, const double *B, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t Q, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in fir_fft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else
    {
        //Set nfft
        size_t nfft = 1u;
        while (nfft<L+Q) { nfft *= 2u; }

        //Initialize fftw
        double *X1, *B1, *Y1, *D1, *Z1, yr, yi, dr, di;
        X1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        B1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        D1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Z1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan xplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        fftw_plan bplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)B1,(fftw_complex *)D1,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!xplan) { fprintf(stderr,"error in fir_fft_z: problem creating fftw plan"); return 1; }
        if (!bplan) { fprintf(stderr,"error in fir_fft_z: problem creating fftw plan"); return 1; }
        fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1,(fftw_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
        if (!iplan) { fprintf(stderr,"error in fir_fft_z: problem creating fftw plan"); return 1; }
        for (size_t n=2u*L; n<2u*nfft; ++n) { X1[n] = 0.0; }
        for (size_t n=nfft; n<2u*nfft; ++n) { Y1[n] = 0.0; }

        //Get D1 (FFT of B1, then scale)
        for (size_t n=0u; n<=Q; ++n, ++B, ++B1) { *B1 = *B; *++B1 = *++B; }
        for (size_t n=Q+1u; n<nfft; ++n, ++B1) { *B1 = 0.0; *++B1 = 0.0; }
        B1 -= 2u*nfft;
        fftw_execute(bplan);
        for (size_t l=0u; l<2u*nfft; ++l, ++D1) { *D1 /= (double)nfft; }
        D1 -= 2u*nfft;

        if (L==N)
        {
            for (size_t l=0u; l<2u*L; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*L;
            fftw_execute(xplan);
            for (size_t n=0u; n<nfft; ++n)
            {
                yr = *Y1++; yi = *Y1--;
                dr = *D1++; di = *D1++;
                *Y1++ = dr*yr - di*yi;
                *Y1++ = dr*yi + di*yr;
            }
            D1 -= 2u*nfft; Y1 -= 2u*nfft;
            fftw_execute(iplan);
            for (size_t l=0u; l<2u*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
            Z1 -= 2u*L;
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
                    for (size_t l=0u; l<2u*L; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*L;
                    fftw_execute(xplan);
                    for (size_t n=0u; n<nfft; ++n)
                    {
                        yr = *Y1++; yi = *Y1--;
                        dr = *D1++; di = *D1++;
                        *Y1++ = dr*yr - di*yi;
                        *Y1++ = dr*yi + di*yr;
                    }
                    D1 -= 2u*nfft; Y1 -= 2u*nfft;
                    fftw_execute(iplan);
                    for (size_t l=0u; l<2u*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
                    Z1 -= 2u*L;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*BS*(L-1u), Y+=2u*BS*(L-1u))
                {
                    for (size_t bs=0u; bs<BS; ++bs, X-=2u*K*L-2u, Y-=2u*K*L-2u)
                    {
                        for (size_t l=0u; l<L; ++l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                        X1 -= 2u*L;
                        fftw_execute(xplan);
                        for (size_t n=0u; n<nfft; ++n)
                        {
                            yr = *Y1++; yi = *Y1--;
                            dr = *D1++; di = *D1++;
                            *Y1++ = dr*yr - di*yi;
                            *Y1++ = dr*yi + di*yr;
                        }
                        D1 -= 2u*nfft; Y1 -= 2u*nfft;
                        fftw_execute(iplan);
                        for (size_t l=0u; l<L; ++l, ++Z1, Y+=2u*K) { *Y = *Z1; *(Y+1) = *++Z1; }
                        Z1 -= 2u*L;
                    }
                }
            }
        }
        fftw_free(X1); fftw_free(D1); fftw_free(Y1); fftw_free(Z1);
        fftw_destroy_plan(xplan); fftw_destroy_plan(bplan); fftw_destroy_plan(iplan);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
