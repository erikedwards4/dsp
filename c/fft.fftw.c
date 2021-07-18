//Does 1-D FFT (fast Fourier transform) of each vector in X along dim.
//The output Y is complex-valued and has the same size as X,
//except along dim, where Y has length Ly = nfft/2 + 1 for real-valued X,
//and length Ly = nfft for complex-valued X.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

//I tried parallel versions with OpenMP, but much slower (have to make fftw_plan P times!).

#include <stdio.h>
#include <math.h>
#include <fftw3.h>
//#include <time.h>

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

int fft_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int fft_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int fft_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int fft_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);


int fft_fftw_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_fftw_s: dim must be in [0 3]\n"); return 1; }
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft/2u + 1u;
    if (nfft<Lx) { fprintf(stderr,"error in fft_fftw_s: nfft must be >= L (vec length)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        Y -= 2u*N;
    }
    else
    {
        //Initialize fftwf
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*Ly*sizeof(float));
        fftwf_plan plan = fftwf_plan_dft_r2c_1d((int)nfft,X1,(fftwf_complex *)Y1,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in fft_fftw_s: problem creating fftw plan"); return 1; }
        for (size_t l=0u; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
        X1 -= nfft;

        if (Lx==N)
        {
            for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= Lx;
            fftwf_execute(plan);
            for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
            Y1 -= 2u*Ly; Y -= 2u*Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y1-=2u*Ly)
                {
                    for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftwf_execute(plan);
                    for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
                Y -= 2u*Ly*V;
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= Lx;
                        fftwf_execute(plan);
                        for (size_t l=0u; l<Ly; ++l, ++Y1, Y+=2u*K-1u) { *Y = *Y1; *++Y = *++Y1; }
                    }
                }
                Y -= 2u*G*B*Ly;
            }
        }
        fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    }

    //Scale
    if (sc)
    {
        const float s = (float)(1.0/sqrt(2u*nfft));
        for (size_t l=0u; l<2u*Ly*N/Lx; ++l, ++Y) { *Y *= s; }
    }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


int fft_fftw_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_fftw_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft/2u + 1u;
    if (nfft<Lx) { fprintf(stderr,"error in fft_fftw_d: nfft must be >= L (vec length)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *Y++ = *X; *Y++ = 0.0; }
        Y -= 2u*N;
    }
    else
    {
        //Initialize fftw
        double *X1, *Y1;
        X1 = (double *)fftw_malloc(nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*Ly*sizeof(double));
        fftw_plan plan = fftw_plan_dft_r2c_1d((int)nfft,X1,(fftw_complex *)Y1,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in fft_fftw_d: problem creating fftw plan"); return 1; }
        for (size_t l=0u; l<nfft; ++l, ++X1) { *X1 = 0.0; }
        X1 -= nfft;

        if (Lx==N)
        {
            for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= nfft;
            fftw_execute(plan);
            for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
            Y1 -= 2u*Ly; Y -= 2u*Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y1-=2u*Ly)
                {
                    for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    fftw_execute(plan);
                    for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
                Y -= 2u*Ly*V;
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= Lx;
                        fftw_execute(plan);
                        for (size_t l=0u; l<Ly; ++l, ++Y1, Y+=2u*K-1u) { *Y = *Y1; *++Y = *++Y1; }
                    }
                }
                Y -= 2u*G*B*Ly;
            }
        }
        fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    }

    //Scale
    if (sc)
    {
        const double s = 1.0/sqrt(2u*nfft);
        for (size_t l=0u; l<2u*Ly*N/Lx; ++l, ++Y) { *Y *= s; }
    }

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


int fft_fftw_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_fftw_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (nfft<Lx) { fprintf(stderr,"error in fft_fftw_c: nfft must be >= L (vec length)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Initialize fftwf
        float *X1, *Y1;
        X1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        Y1 = (float *)fftwf_malloc(2u*nfft*sizeof(float));
        fftwf_plan plan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in fft_fftw_c: problem creating fftw plan"); return 1; }
        for (size_t l=0u; l<2u*nfft; ++l, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*Lx;
            fftwf_execute(plan);
            for (size_t l=0u; l<2u*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
            Y1 -= 2u*nfft; Y -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y1-=2u*nfft)
                {
                    for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*Lx;
                    fftwf_execute(plan);
                    for (size_t l=0u; l<2u*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
                Y -= 2u*nfft*V;
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-1u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                        X1 -= 2u*Lx;
                        fftwf_execute(plan);
                        for (size_t l=0u; l<nfft; ++l, ++Y1, Y+=2u*K-1u) { *Y = *Y1; *++Y = *++Y1; }
                    }
                }
                Y -= 2u*G*B*nfft;
            }
        }
        fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    }
    
    //Scale
    if (sc)
    {
        const float s = (float)(1.0f/sqrt(2u*nfft));
        for (size_t l=0u; l<2u*nfft*N/Lx; ++l, ++Y) { *Y *= s; }
        Y -= 2u*nfft;
    }
    
    return 0;
}


int fft_fftw_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_fftw_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (nfft<Lx) { fprintf(stderr,"error in fft_fftw_z: nfft must be >= L (vec length)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=0u; n<2u*N; ++n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        //Initialize fftw
        double *X1, *Y1;
        X1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        Y1 = (double *)fftw_malloc(2u*nfft*sizeof(double));
        fftw_plan plan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_FORWARD,FFTW_ESTIMATE);
        if (!plan) { fprintf(stderr,"error in fft_fftw_z: problem creating fftw plan"); return 1; }
        for (size_t l=0u; l<2u*nfft; ++l, ++X1) { *X1 = 0.0; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*Lx;
            fftw_execute(plan);
            for (size_t l=0u; l<2u*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
            Y1 -= 2u*nfft; Y -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=0u; v<V; ++v, Y1-=2u*nfft)
                {
                    for (size_t l=0u; l<2u*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*Lx;
                    fftw_execute(plan);
                    for (size_t l=0u; l<2u*nfft; ++l, ++Y1, ++Y) { *Y = *Y1; }
                }
                Y -= 2u*nfft*V;
            }
            else
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-1u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=0u; l<Lx; ++l, X+=2u*K-1u, ++X1) { *X1 = *X; *++X1 = *++X; }
                        X1 -= 2u*Lx;
                        fftw_execute(plan);
                        for (size_t l=0u; l<nfft; ++l, ++Y1, Y+=2u*K-1u) { *Y = *Y1; *++Y = *++Y1; }
                    }
                }
                Y -= 2u*G*B*nfft;
            }
        }
        fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    }
    
    //Scale
    if (sc)
    {
        const double s = 1.0/sqrt(2u*nfft);
        for (size_t l=0u; l<2u*nfft*N/Lx; ++l, ++Y) { *Y *= s; }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
