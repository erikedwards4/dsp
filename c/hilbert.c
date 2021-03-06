//Does 1-D Hilbert transform of each vector in X along dim.
//The output Y is complex-valued and has the same size as X.
//Y is the analytic signal, with X in the real part,
//and the actual Hilbert transform in the imaginary part.

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int hilbert_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft);
int hilbert_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft);


int hilbert_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft)
{
    if (dim>3) { fprintf(stderr,"error in hilbert_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (nfft<L) { fprintf(stderr,"error in hilbert_s: nfft must be >= L (vec length)\n"); return 1; }
    const float sc = 2.0f / nfft;
    const size_t isodd = nfft%2;

    //Initialize fftwf
    float *X1, *Y1, *Z1;
    X1 = (float *)fftwf_malloc(nfft*sizeof(float));
    Y1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    Z1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    fftwf_plan fplan = fftwf_plan_dft_r2c_1d((int)nfft,X1,(fftwf_complex *)Y1,FFTW_ESTIMATE);
    if (!fplan) { fprintf(stderr,"error in hilbert_s: problem creating fftw plan"); return 1; }
    fftwf_plan iplan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)Y1,(fftwf_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!iplan) { fprintf(stderr,"error in hilbert_s: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; *++Y = 0.0f; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
        X1 -= nfft;
        fftwf_execute(fplan);
        *Y1++ /= nfft; ++Y1;
        for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
        if (!isodd) { *Y1++ /= nfft; }
        for (size_t l=nfft+1; l<2*nfft; ++l, ++Y1) { *Y1 = 0.0f; }
        Y1 -= 2*nfft;
        fftwf_execute(iplan);
        for (size_t l=0; l<2*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
        Z1 -= 2*L;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= nfft;
            fftwf_execute(fplan);
            *Y1++ /= nfft; ++Y1;
            for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
            if (!isodd) { *Y1++ /= nfft; }
            for (size_t l=nfft+1; l<2*nfft; ++l, ++Y1) { *Y1 = 0.0f; }
            Y1 -= 2*nfft;
            fftwf_execute(iplan);
            for (size_t l=0; l<2*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
            Z1 -= 2*L;
            for (size_t v=1; v<V; ++v, Z1-=2*L)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                fftwf_execute(fplan);
                *Y1++ /= nfft; ++Y1;
                for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
                if (!isodd) { *Y1++ /= nfft; }
                Y1 -= nfft + 1;
                fftwf_execute(iplan);
                for (size_t l=0; l<2*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
            }
        }
        else
        {
            X1 += L; Y1 += nfft + 1;
            for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
            for (size_t l=nfft+1; l<2*nfft; ++l, ++Y1) { *Y1 = 0.0f; }
            X1 -= nfft; Y1 -= 2*nfft;
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Z1-=2*L, Y-=2*K*L-2)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftwf_execute(fplan);
                    *Y1++ /= nfft; ++Y1;
                    for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
                    if (!isodd) { *Y1++ /= nfft; }
                    Y1 -= nfft + 1;
                    fftwf_execute(iplan);
                    for (size_t l=0; l<L; ++l, ++Z1, Y+=2*K-1) { *Y = *Z1; *++Y = *++Z1; }
                }
            }
        }
    }
    
    fftwf_free(X1); fftwf_free(Y1); fftwf_free(Z1);
    fftwf_destroy_plan(fplan); fftwf_destroy_plan(iplan);
    return 0;
}


int hilbert_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft)
{
    if (dim>3) { fprintf(stderr,"error in hilbert_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (nfft<L) { fprintf(stderr,"error in hilbert_d: nfft must be >= L (vec length)\n"); return 1; }
    const double sc = 2.0 / nfft;
    const size_t isodd = nfft%2;

    //Initialize fftw
    double *X1, *Y1, *Z1;
    X1 = (double *)fftw_malloc(nfft*sizeof(double));
    Y1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    Z1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    fftw_plan fplan = fftw_plan_dft_r2c_1d((int)nfft,X1,(fftw_complex *)Y1,FFTW_ESTIMATE);
    if (!fplan) { fprintf(stderr,"error in hilbert_d: problem creating fftw plan"); return 1; }
    fftw_plan iplan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)Y1,(fftw_complex *)Z1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!iplan) { fprintf(stderr,"error in hilbert_d: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (L==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; *++Y = 0.0; }
    }
    else if (L==N)
    {
        for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0; }
        X1 -= nfft;
        fftw_execute(fplan);
        *Y1++ /= nfft; ++Y1;
        for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
        if (!isodd) { *Y1++ /= nfft; }
        for (size_t l=nfft+1; l<2*nfft; ++l, ++Y1) { *Y1 = 0.0; }
        Y1 -= 2*nfft;
        fftw_execute(iplan);
        for (size_t l=0; l<2*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
        Z1 -= 2*L;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0; }
            X1 -= nfft;
            fftw_execute(fplan);
            *Y1++ /= nfft; ++Y1;
            for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
            if (!isodd) { *Y1++ /= nfft; }
            for (size_t l=nfft+1; l<2*nfft; ++l, ++Y1) { *Y1 = 0.0; }
            Y1 -= 2*nfft;
            fftw_execute(iplan);
            for (size_t l=0; l<2*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
            Z1 -= 2*L;
            for (size_t v=1; v<V; ++v, Z1-=2*L)
            {
                for (size_t l=0; l<L; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= L;
                fftw_execute(fplan);
                *Y1++ /= nfft; ++Y1;
                for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
                if (!isodd) { *Y1++ /= nfft; }
                Y1 -= nfft + 1;
                fftw_execute(iplan);
                for (size_t l=0; l<2*L; ++l, ++Z1, ++Y) { *Y = *Z1; }
            }
        }
        else
        {
            X1 += L; Y1 += nfft + 1;
            for (size_t l=L; l<nfft; ++l, ++X1) { *X1 = 0.0; }
            for (size_t l=nfft+1; l<2*nfft; ++l, ++Y1) { *Y1 = 0.0; }
            X1 -= nfft; Y1 -= 2*nfft;
            for (size_t g=0; g<G; ++g, X+=B*(L-1), Y+=2*B*(L-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*L-1, Z1-=2*L, Y-=2*K*L-2)
                {
                    for (size_t l=0; l<L; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= L;
                    fftw_execute(fplan);
                    *Y1++ /= nfft; ++Y1;
                    for (size_t l=2; l<nfft+isodd; ++l, ++Y1) { *Y1 *= sc; }
                    if (!isodd) { *Y1++ /= nfft; }
                    Y1 -= nfft + 1;
                    fftw_execute(iplan);
                    for (size_t l=0; l<L; ++l, ++Z1, Y+=2*K-1) { *Y = *Z1; *++Y = *++Z1; }
                }
            }
        }
    }
    
    fftw_free(X1); fftw_free(Y1); fftw_free(Z1);
    fftw_destroy_plan(fplan); fftw_destroy_plan(iplan);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
