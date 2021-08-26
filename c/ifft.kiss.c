//Does 1-D IFFT (inverse fast Fourier transform) of each vector in X along dim.
//The input X is complex-valued.
//The output Y has the same size as X, except along dim, where Y has length nfft.
//Y is real-valued for ifft_kiss_s and ifft_kiss_d,
//and complex-valued for ifft_kiss_c and ifft_kiss_z.
//In the former case, X has only nonnegative freqs, so Lx = nfft/2 + 1.

//If sc, then scales Y by 2.0*sqrt(2.0)/nfft, so that invertible with ifft.

//This uses the fully-open-source KISS FFTS (Keep-it-simple-stupid FFT library) library
//https://github.com/mborgerding/kissfft

//I installed into: /opt/kissfft.
//So, either use that, change the #include, and/or use the -I flag.

//This definitely requires allocating Y1 and copying Y1 into Y.
//However, this appears to work fine without a separate X1,
//for a single FFT only; but fails for multiple FFTs, so using X1 for now.

//Odd-length FFT supported, but must make X1 complex.

//Double-precision is not supported!
//So, I cast to float for X1, Y1 to support this, and issue a warning.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/opt/kissfft/kiss_fft.h"
#include "/opt/kissfft/kiss_fft.c"
#include "/opt/kissfft/kiss_fftr.h"
#include "/opt/kissfft/kiss_fftr.c"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ifft_kiss_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_kiss_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_kiss_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_kiss_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);


int ifft_kiss_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_kiss_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float s = (sc) ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_kiss_s: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize IFFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in ifft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,NULL,NULL);
        #endif
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = *X; X1->i = *++X; }
            X -= 2u + 2u*(1u-nfft%2u);
            for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { X1->r = *X; X1->i = -*(X+1u); }
            X1 -= nfft;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = Y1->r * s; }
            Y1 -= nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = *X; X1->i = *++X; }
                    X -= 2u + 2u*(1u-nfft%2u);
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { X1->r = *X; X1->i = -*(X+1u); }
                    X1 -= nfft;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = Y1->r * s; }
                    Y1 -= nfft;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { X1->r = *X; X1->i = *(X+1u); }
                        X -= K*(2u + 2u*(1u-nfft%2u));
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { X1->r = *X; X1->i = -*(X+1u); }
                        X1 -= nfft;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=K) { *Y = Y1->r * s; }
                        Y1 -= nfft;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }

    return 0;
}


int ifft_kiss_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_kiss_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double s = (sc) ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_kiss_d: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in ifft_kiss_d: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize IFFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in ifft_kiss_d: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_kiss_d: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,NULL,NULL);
        #endif
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*++X; }
            X -= 2u + 2u*(1u-nfft%2u);
            for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = -(kiss_fft_scalar)*(X+1u); }
            X1 -= nfft;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)Y1->r * s; }
            Y1 -= nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*++X; }
                    X -= 2u + 2u*(1u-nfft%2u);
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = -(kiss_fft_scalar)*(X+1u); }
                    X1 -= nfft;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)Y1->r * s; }
                    Y1 -= nfft;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*(X+1u); }
                        X -= K*(2u + 2u*(1u-nfft%2u));
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = -(kiss_fft_scalar)*(X+1u); }
                        X1 -= nfft;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=K) { *Y = (double)Y1->r * s; }
                        Y1 -= nfft;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }

    return 0;
}


int ifft_kiss_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_kiss_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float s = (sc) ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_kiss_c: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize IFFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in ifft_kiss_c: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_kiss_c: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,NULL,NULL);
        #endif
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = *X; X1->i = *++X; }
            X1 -= nfft;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = Y1->r * s; *++Y = Y1->i * s; }
            Y1 -= nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = *X; X1->i = *++X; }
                    X1 -= nfft;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = Y1->r * s; *++Y = Y1->i * s; }
                    Y1 -= nfft;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { X1->r = *X; X1->i = *(X+1); }
                        X1 -= nfft;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K) { *Y = Y1->r * s; *(Y+1) = Y1->i * s; }
                        Y1 -= nfft;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }

    return 0;
}


int ifft_kiss_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_kiss_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double s = (sc) ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_kiss_z: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in ifft_kiss_z: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize IFFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in ifft_kiss_z: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_kiss_z: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,1,NULL,NULL);
        #endif
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*++X; }
            X1 -= nfft;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)Y1->r * s; *++Y = (double)Y1->i * s; }
            Y1 -= nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*++X; }
                    X1 -= nfft;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)Y1->r * s; *++Y = (double)Y1->i * s; }
                    Y1 -= nfft;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*(X+1); }
                        X1 -= nfft;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K) { *Y = (double)Y1->r * s; *(Y+1) = (double)Y1->i * s; }
                        Y1 -= nfft;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
