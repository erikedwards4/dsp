//Does 1-D FFT (fast Fourier transform) of each vector in X along dim.
//The output Y is complex-valued and has the same size as X,
//except along dim, where Y has length Ly = nfft/2 + 1 for real-valued X,
//and length Ly = nfft for complex-valued X.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

//I tried parallel versions with OpenMP, but much slower (have to make fftw_plan P times!).

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <time.h>

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifndef M_SQRT1_2f
    #define M_SQRT1_2f 0.707106781186547524401f
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int fft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int fft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);


int fft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in fft_s: dim must be in [0 3]\n"); return 1; }
    struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = nfft/2 + 1;
    if (nfft<Lx) { fprintf(stderr,"error in fft_s: nfft must be >= L (vec length)\n"); return 1; }

    //Initialize fftwf
    float *X1, *Y1;
    X1 = (float *)fftwf_malloc(nfft*sizeof(float));
    Y1 = (float *)fftwf_malloc(2*Ly*sizeof(float));
    fftwf_plan plan = fftwf_plan_dft_r2c_1d((int)nfft,X1,(fftwf_complex *)Y1,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_s: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && Ly==1)
    {
        if (sc)
        {
            for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X * M_SQRT1_2f; *++Y = 0.0f; }
        }
        else
        {
            for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; *++Y = 0.0f; }
        }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
        X1 -= nfft;
        fftwf_execute(plan);
        if (sc)
        {
            const float s = (float)(1.0/sqrt(2*Lx));
            for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        }
        else
        {
            for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
        }
        Y1 -= 2*Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= nfft;
            fftwf_execute(plan);
            if (sc)
            {
                const float s = (float)(1.0/sqrt(2*Lx));
                for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2*Ly;
                for (size_t v=1; v<V; ++v, Y1-=2*Ly)
                {
                    for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftwf_execute(plan);
                    for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
                Y1 -= 2*Ly;
                for (size_t v=1; v<V; ++v, Y1-=2*Ly)
                {
                    for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftwf_execute(plan);
                    for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
            }
        }
        else
        {
            const float s = sc ? (float)(1.0/sqrt(2*Lx)) : 1.0f;
            X1 += Lx;
            for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= nfft;
            for (size_t g=0; g<G; ++g, X+=B*(Lx-1), Y+=2*B*(Ly-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*Lx-1, Y1-=2*Ly, Y-=2*K*Ly-2)
                {
                    for (size_t l=0; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftwf_execute(plan);
                    for (size_t l=0; l<Ly; ++l, ++Y1, Y+=2*K-1) { *Y = *Y1*s; *++Y = *++Y1*s; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


int fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in fft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const size_t Ly = nfft/2 + 1;
    if (nfft<Lx) { fprintf(stderr,"error in fft_d: nfft must be >= L (vec length)\n"); return 1; }

    //Initialize fftw
    double *X1, *Y1;
    X1 = (double *)fftw_malloc(nfft*sizeof(double));
    Y1 = (double *)fftw_malloc(2*Ly*sizeof(double));
    fftw_plan plan = fftw_plan_dft_r2c_1d((int)nfft,X1,(fftw_complex *)Y1,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_d: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && Ly==1)
    {
        if (sc)
        {
            for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X * M_SQRT1_2; *++Y = 0.0; }
        }
        else
        {
            for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X; *++Y = 0.0; }
        }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0; }
        X1 -= nfft;
        fftw_execute(plan);
        if (sc)
        {
            const double s = 1.0/sqrt(2*Lx);
            for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        }
        else
        {
            for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
        }
        Y1 -= 2*Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0; }
            X1 -= nfft;
            fftw_execute(plan);
            if (sc)
            {
                const double s = 1.0/sqrt(2*Lx);
                for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2*Ly;
                for (size_t v=1; v<V; ++v, Y1-=2*Ly)
                {
                    for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftw_execute(plan);
                    for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
                Y1 -= 2*Ly;
                for (size_t v=1; v<V; ++v, Y1-=2*Ly)
                {
                    for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftw_execute(plan);
                    for (size_t l=0; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
            }
        }
        else
        {
            const double s = sc ? 1.0/sqrt(2*Lx) : 1.0;
            X1 += Lx;
            for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0; }
            X1 -= nfft;
            for (size_t g=0; g<G; ++g, X+=B*(Lx-1), Y+=2*B*(Ly-1))
            {
                for (size_t b=0; b<B; ++b, X-=K*Lx-1, Y1-=2*Ly, Y-=2*K*Ly-2)
                {
                    for (size_t l=0; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftw_execute(plan);
                    for (size_t l=0; l<Ly; ++l, ++Y1, Y+=2*K-1) { *Y = *Y1*s; *++Y = *++Y1*s; }
                }
            }
        }
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


int fft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in fft_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (nfft<Lx) { fprintf(stderr,"error in fft_c: nfft must be >= L (vec length)\n"); return 1; }

    //Initialize fftwf
    float *X1, *Y1;
    X1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    Y1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    fftwf_plan plan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_c: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && nfft==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
        if (sc)
        {
            for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X * M_SQRT1_2f; }
        }
        else
        {
            for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
        }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=2*Lx; l<2*nfft; ++l, ++X1) { *X1 = 0.0f; }
        X1 -= 2*nfft;
        fftwf_execute(plan);
        if (sc)
        {
            const float s = (float)(1.0/sqrt(2*Lx));
            for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        }
        else
        {
            for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
        }
        Y1 -= 2*nfft;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=2*Lx; l<2*nfft; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= 2*nfft;
            fftwf_execute(plan);
            if (sc)
            {
                const float s = (float)(1.0/sqrt(2*Lx));
                for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2*nfft;
                for (size_t v=1; v<V; ++v, Y1-=2*nfft)
                {
                    for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2*Lx;
                    fftwf_execute(plan);
                    for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
                Y1 -= 2*nfft;
                for (size_t v=1; v<V; ++v, Y1-=2*nfft)
                {
                    for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2*Lx;
                    fftwf_execute(plan);
                    for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
            }
        }
        else
        {
            const float s = sc ? (float)(1.0/sqrt(2*Lx)) : 1.0f;
            X1 += 2*Lx;
            for (size_t l=2*Lx; l<2*nfft; ++l, ++X1) { *X1 = 0.0f; }
            X1 -= 2*nfft;
            for (size_t g=0; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(nfft-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*Lx-1, Y-=2*K*nfft-1)
                {
                    for (size_t l=0; l<2*Lx; ++l, X+=2*K, ++X1) { *X1 = *X; }
                    X1 -= 2*Lx;
                    fftwf_execute(plan);
                    for (size_t l=0; l<2*nfft; ++l, ++Y1, Y+=2*K) { *Y = *Y1 * s; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int fft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in fft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    if (nfft<Lx) { fprintf(stderr,"error in fft_z: nfft must be >= L (vec length)\n"); return 1; }

    //Initialize fftw
    double *X1, *Y1;
    X1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    Y1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    fftw_plan plan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in fft_z: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && nfft==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
        if (sc)
        {
            for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X * M_SQRT1_2; }
        }
        else
        {
            for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X; }
        }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
        for (size_t l=2*Lx; l<2*nfft; ++l, ++X1) { *X1 = 0.0; }
        X1 -= 2*nfft;
        fftw_execute(plan);
        if (sc)
        {
            const double s = 1.0/sqrt(2*Lx);
            for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        }
        else
        {
            for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
        }
        Y1 -= 2*nfft;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
            for (size_t l=2*Lx; l<2*nfft; ++l, ++X1) { *X1 = 0.0; }
            X1 -= 2*nfft;
            fftw_execute(plan);
            if (sc)
            {
                const double s = 1.0/sqrt(2*Lx);
                for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2*nfft;
                for (size_t v=1; v<V; ++v, Y1-=2*nfft)
                {
                    for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2*Lx;
                    fftw_execute(plan);
                    for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
                Y1 -= 2*nfft;
                for (size_t v=1; v<V; ++v, Y1-=2*nfft)
                {
                    for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2*Lx;
                    fftw_execute(plan);
                    for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
            }
        }
        else
        {
            const double s = sc ? 1.0/sqrt(2*Lx) : 1.0;
            X1 += 2*Lx;
            for (size_t l=2*Lx; l<2*nfft; ++l, ++X1) { *X1 = 0.0; }
            X1 -= 2*nfft;
            for (size_t g=0; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(nfft-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*Lx-1, Y-=2*K*nfft-1)
                {
                    for (size_t l=0; l<2*Lx; ++l, X+=2*K, ++X1) { *X1 = *X; }
                    X1 -= 2*Lx;
                    fftw_execute(plan);
                    for (size_t l=0; l<2*nfft; ++l, ++Y1, Y+=2*K) { *Y = *Y1 * s; }
                }
            }
        }
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
