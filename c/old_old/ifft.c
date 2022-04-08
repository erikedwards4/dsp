//Does 1-D IFFT (inverse fast Fourier transform) of each vector in X along dim.
//The input X is complex-valued.
//The output Y has the same size as X, except along dim, where Y has length nfft.
//Y is real-valued for ifft_s and ifft_d,
//and complex-valued for ifft_c and ifft_z.
//In the former case, X has only nonnegative freqs, so Lx = nfft/2 + 1.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

#include <stdio.h>
#include <fftw3.h>

#ifndef M_SQRT2
    #define M_SQRT2 1.41421356237309504880
#endif

#ifndef M_SQRT2f
    #define M_SQRT2f 1.41421356237309504880f
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ifft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int ifft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int ifft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int ifft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);


int ifft_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in ifft_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float s = (sc) ? M_SQRT2f : 1.0f/nfft;
    if (Lx!=nfft/2+1) { fprintf(stderr,"error in ifft_s: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    //Initialize fftwf
    float *X1, *Y1;
    X1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    Y1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    fftwf_plan plan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in ifft_s: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && nfft==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
        X -= 2 + 2*(1-nfft%2);
        for (size_t l=Lx; l<nfft; ++l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
        X1 -= 2*nfft;
        fftwf_execute(plan);
        for (size_t l=0; l<nfft; ++l, Y1+=2, ++Y) { *Y = *Y1 * s; }
        Y1 -= 2*nfft;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
            X -= 2 + 2*(1-nfft%2);
            for (size_t l=Lx; l<nfft; ++l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
            X1 -= 2*nfft;
            fftwf_execute(plan);
            for (size_t l=0; l<nfft; ++l, Y1+=2, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2*nfft;
            for (size_t v=1; v<V; ++v, X+=2*Lx, Y1-=2*nfft)
            {
                for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
                X -= 2 + 2*(1-nfft%2);
                for (size_t l=Lx; l<nfft; ++l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
                X1 -= 2*nfft;
                fftwf_execute(plan);
                for (size_t l=0; l<nfft; ++l, Y1+=2, ++Y) { *Y = *Y1 * s; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(Lx-1), Y+=B*(nfft-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, Y1-=2*nfft, Y-=K*nfft-1)
                {
                    for (size_t l=0; l<Lx; ++l, X+=2*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                    X -= K*(2 + 2*(1-nfft%2));
                    for (size_t l=Lx; l<nfft; ++l, X-=2*K, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
                    X1 -= 2*nfft;
                    fftwf_execute(plan);
                    for (size_t l=0; l<nfft; ++l, Y1+=2, Y+=K) { *Y = *Y1 * s; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int ifft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in ifft_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double s = (sc) ? M_SQRT2 : 1.0/nfft;
    if (Lx!=nfft/2+1) { fprintf(stderr,"error in ifft_d: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    //Initialize fftw
    double *X1, *Y1;
    X1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    Y1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    fftw_plan plan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in ifft_d: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && nfft==1)
    {
        for (size_t n=0; n<N; ++n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
        X -= 2 + 2*(1-nfft%2);
        for (size_t l=Lx; l<nfft; ++l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
        X1 -= 2*nfft;
        fftw_execute(plan);
        for (size_t l=0; l<nfft; ++l, Y1+=2, ++Y) { *Y = *Y1 * s; }
        Y1 -= 2*nfft;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0) ? 1 : (dim==1) ? R : (dim==2) ? R*C : R*C*S) : ((dim==0) ? C*S*H : (dim==1) ? S*H : (dim==2) ? H : 1);
        const size_t B = (iscolmajor && dim==0) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1 && (G==1 || B==1))
        {
            for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
            X -= 2 + 2*(1-nfft%2);
            for (size_t l=Lx; l<nfft; ++l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
            X1 -= 2*nfft;
            fftw_execute(plan);
            for (size_t l=0; l<nfft; ++l, Y1+=2, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2*nfft;
            for (size_t v=1; v<V; ++v, X+=2*Lx, Y1-=2*nfft)
            {
                for (size_t l=0; l<Lx; ++l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
                X -= 2 + 2*(1-nfft%2);
                for (size_t l=Lx; l<nfft; ++l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
                X1 -= 2*nfft;
                fftw_execute(plan);
                for (size_t l=0; l<nfft; ++l, Y1+=2, ++Y) { *Y = *Y1 * s; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(Lx-1), Y+=B*(nfft-1))
            {
                for (size_t b=0; b<B; ++b, X+=2, Y1-=2*nfft, Y-=K*nfft-1)
                {
                    for (size_t l=0; l<Lx; ++l, X+=2*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                    X -= K*(2 + 2*(1-nfft%2));
                    for (size_t l=Lx; l<nfft; ++l, X-=2*K, ++X1) { *X1 = *X; *++X1 = -*(X+1); }
                    X1 -= 2*nfft;
                    fftw_execute(plan);
                    for (size_t l=0; l<nfft; ++l, Y1+=2, Y+=K) { *Y = *Y1 * s; }
                }
            }
        }
    }
    
    fftw_destroy_plan(plan); fftw_free(X1); fftw_free(Y1);
    return 0;
}


int ifft_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in ifft_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const float s = (sc) ? M_SQRT2f : 1.0f/nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_c: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    //Initialize fftwf
    float *X1, *Y1;
    X1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    Y1 = (float *)fftwf_malloc(2*nfft*sizeof(float));
    fftwf_plan plan = fftwf_plan_dft_1d((int)nfft,(fftwf_complex *)X1,(fftwf_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in ifft_c: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && nfft==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= 2*nfft;
        fftwf_execute(plan);
        for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
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
            X1 -= 2*nfft;
            fftwf_execute(plan);
            for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2*nfft;
            for (size_t v=1; v<V; ++v, Y1-=2*nfft)
            {
                for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= 2*nfft;
                fftwf_execute(plan);
                for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(nfft-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*Lx-2, Y1-=2*nfft, Y-=2*K*nfft-2)
                {
                    for (size_t l=0; l<Lx; ++l, X+=2*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                    X1 -= 2*nfft;
                    fftwf_execute(plan);
                    for (size_t l=0; l<nfft; ++l, ++Y1, Y+=2*K-1) { *Y = *Y1 * s; *++Y = *++Y1 * s; }
                }
            }
        }
    }
    
    fftwf_destroy_plan(plan); fftwf_free(X1); fftwf_free(Y1);
    return 0;
}


int ifft_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3) { fprintf(stderr,"error in ifft_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0) ? R : (dim==1) ? C : (dim==2) ? S : H;
    const double s = (sc) ? M_SQRT2 : 1.0/nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_z: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    //Initialize fftw
    double *X1, *Y1;
    X1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    Y1 = (double *)fftw_malloc(2*nfft*sizeof(double));
    fftw_plan plan = fftw_plan_dft_1d((int)nfft,(fftw_complex *)X1,(fftw_complex *)Y1,FFTW_BACKWARD,FFTW_ESTIMATE);
    if (!plan) { fprintf(stderr,"error in ifft_z: problem creating fftw plan"); return 1; }

    if (N==0) {}
    else if (Lx==1 && nfft==1)
    {
        for (size_t n=0; n<2*N; ++n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (Lx==N)
    {
        for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
        X1 -= 2*nfft;
        fftw_execute(plan);
        for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
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
            X1 -= 2*nfft;
            fftw_execute(plan);
            for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2*nfft;
            for (size_t v=1; v<V; ++v, Y1-=2*nfft)
            {
                for (size_t l=0; l<2*Lx; ++l, ++X, ++X1) { *X1 = *X; }
                X1 -= 2*nfft;
                fftw_execute(plan);
                for (size_t l=0; l<2*nfft; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
            }
        }
        else
        {
            for (size_t g=0; g<G; ++g, X+=2*B*(Lx-1), Y+=2*B*(nfft-1))
            {
                for (size_t b=0; b<B; ++b, X-=2*K*Lx-2, Y1-=2*nfft, Y-=2*K*nfft-2)
                {
                    for (size_t l=0; l<Lx; ++l, X+=2*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                    X1 -= 2*nfft;
                    fftw_execute(plan);
                    for (size_t l=0; l<nfft; ++l, ++Y1, Y+=2*K-1) { *Y = *Y1 * s; *++Y = *++Y1 * s; }
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
