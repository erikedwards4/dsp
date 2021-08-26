//Does 1-D FFT (fast Fourier transform) of each vector in X along dim.
//The output Y is complex-valued and has the same size as X,
//except along dim, where Y has length Ly = nfft/2 + 1 for real-valued X,
//and length Ly = nfft for complex-valued X.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

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
#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fft_kiss_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_kiss_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_kiss_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_kiss_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);


int fft_kiss_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_kiss_s: dim must be in [0 3]\n"); return 1; }
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft/2u + 1u;
    if (nfft<Lx) { fprintf(stderr,"error in fft_kiss_s: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        Y -= 2u*N;
    }
    else if (nfft%2u)
    {
        //Initialize FFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,NULL,NULL);
        #endif

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = *X; X1->i = 0.0f; }
            X1 -= Lx;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = Y1->r; *++Y = Y1->i; }
            Y1 -= Ly;
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
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = *X; X1->i = 0.0f; }
                    X1 -= Lx;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = Y1->r; *++Y = Y1->i; }
                    Y1 -= Ly;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, ++X1) { X1->r = *X; X1->i = 0.0f; }
                        X1 -= Lx;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, Y+=2u*K, ++Y1) { *Y = Y1->r; *(Y+1) = Y1->i; }
                        Y1 -= Ly;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }
    else    //nfft is a power-of-2
    {
        //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
        
        //Initialize FFT
        kiss_fft_scalar *X1;    //kiss_fft_scalar is float
        kiss_fft_cpx *Y1;
        X1 = (kiss_fft_scalar *)calloc(sizeof(kiss_fft_scalar),(nfft+2u)*sizeof(kiss_fft_scalar));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fftr_cfg kiss_fftr_state;
        #ifdef __cplusplus
            kiss_fftr_state = kiss_fftr_alloc((int)nfft,0,nullptr,nullptr);
        #else
            kiss_fftr_state = kiss_fftr_alloc((int)nfft,0,NULL,NULL);
        #endif

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= Lx;
            kiss_fftr(kiss_fftr_state,X1,Y1);
            for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = Y1->r; *++Y = Y1->i; }
            Y1 -= Ly;
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
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    kiss_fftr(kiss_fftr_state,X1,Y1);
                    for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = Y1->r; *++Y = Y1->i; }
                    Y1 -= Ly;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= Lx;
                        kiss_fftr(kiss_fftr_state,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, Y+=2u*K, ++Y1) { *Y = Y1->r; *(Y+1) = Y1->i; }
                        Y1 -= Ly;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fftr_state); free(X1); free(Y1);

        //clock_gettime(CLOCK_REALTIME,&toc);
        //fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);
    }

    //Scale
    if (sc)
    {
        const float s = 1.0f/sqrtf((float)(2u*nfft));
        for (size_t l=2u*Ly*N/Lx; l>0u; --l, ++Y) { *Y *= s; }
    }

    return 0;
}


int fft_kiss_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_kiss_d: dim must be in [0 3]\n"); return 1; }    

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft/2u + 1u;
    if (nfft<Lx) { fprintf(stderr,"error in fft_kiss_d: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in fft_kiss_d: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        Y -= 2u*N;
    }
    else if (nfft%2u)
    {
        //Initialize FFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,NULL,NULL);
        #endif

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = 0.0f; }
            X1 -= Lx;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = (double)Y1->r; *++Y = (double)Y1->i; }
            Y1 -= Ly;
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
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = 0.0f; }
                    X1 -= Lx;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = (double)Y1->r; *++Y = (double)Y1->i; }
                    Y1 -= Ly;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = 0.0f; }
                        X1 -= Lx;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, Y+=2u*K, ++Y1) { *Y = (double)Y1->r; *(Y+1) = (double)Y1->i; }
                        Y1 -= Ly;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }
    else    //nfft is a power-of-2
    {
        //Initialize FFT
        kiss_fft_scalar *X1;    //kiss_fft_scalar is float
        kiss_fft_cpx *Y1;
        X1 = (kiss_fft_scalar *)calloc(sizeof(kiss_fft_scalar),(nfft+2u)*sizeof(kiss_fft_scalar));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_kiss_s: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fftr_cfg kiss_fftr_state;
        #ifdef __cplusplus
            kiss_fftr_state = kiss_fftr_alloc((int)nfft,0,nullptr,nullptr);
        #else
            kiss_fftr_state = kiss_fftr_alloc((int)nfft,0,NULL,NULL);
        #endif

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = (kiss_fft_scalar)*X; }
            X1 -= Lx;
            kiss_fftr(kiss_fftr_state,X1,Y1);
            for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = (double)Y1->r; *++Y = (double)Y1->i; }
            Y1 -= Ly;
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
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = (kiss_fft_scalar)*X; }
                    X1 -= Lx;
                    kiss_fftr(kiss_fftr_state,X1,Y1);
                    for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = (double)Y1->r; *++Y = (double)Y1->i; }
                    Y1 -= Ly;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, ++X1) { *X1 = (kiss_fft_scalar)*X; }
                        X1 -= Lx;
                        kiss_fftr(kiss_fftr_state,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, Y+=2u*K, ++Y1) { *Y = (double)Y1->r; *(Y+1) = (double)Y1->i; }
                        Y1 -= Ly;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fftr_state); free(X1); free(Y1);
    }

    //Scale
    if (sc)
    {
        const double s = 1.0/sqrt((double)(2u*nfft));
        for (size_t l=2u*Ly*N/Lx; l>0u; --l, ++Y) { *Y *= s; }
    }

    return 0;
}


int fft_kiss_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_kiss_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft;
    if (nfft<Lx) { fprintf(stderr,"error in fft_kiss_c: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
    }
    else
    {
        //Initialize FFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in fft_kiss_c: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_kiss_c: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,NULL,NULL);
        #endif

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = *X; X1->i = *++X; }
            X1 -= Lx;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = Y1->r; *++Y = Y1->i; }
            Y1 -= Ly;
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
                    X1 -= Lx;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = Y1->r; *++Y = Y1->i; }
                    Y1 -= Ly;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { X1->r = *X; X1->i = *(X+1); }
                        X1 -= Lx;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, Y+=2u*K, ++Y1) { *Y = Y1->r; *(Y+1) = Y1->i; }
                        Y1 -= Ly;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }

    //Scale
    if (sc)
    {
        const float s = 1.0f/sqrtf((float)(2u*nfft));
        for (size_t l=2u*nfft*N/Lx; l>0u; --l, ++Y) { *Y *= s; }
    }

    return 0;
}


int fft_kiss_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_kiss_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft;
    if (nfft<Lx) { fprintf(stderr,"error in fft_kiss_z: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in fft_kiss_z: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
    }
    else
    {
        //Initialize FFT
        kiss_fft_cpx *X1, *Y1;
        X1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        Y1 = (kiss_fft_cpx *)calloc(sizeof(kiss_fft_cpx),nfft*sizeof(kiss_fft_cpx));
        if (!X1) { fprintf(stderr,"error in fft_kiss_z: problem with calloc. "); perror("calloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_kiss_z: problem with calloc. "); perror("calloc"); return 1; }
        kiss_fft_cfg kiss_fft_state;
        #ifdef __cplusplus
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,nullptr,nullptr);
        #else
            kiss_fft_state = kiss_fft_alloc((int)nfft,0,NULL,NULL);
        #endif

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*++X; }
            X1 -= Lx;
            kiss_fft(kiss_fft_state,X1,Y1);
            for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = (double)Y1->r; *++Y = (double)Y1->i; }
            Y1 -= Ly;
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
                    X1 -= Lx;
                    kiss_fft(kiss_fft_state,X1,Y1);
                    for (size_t l=Ly; l>0u; --l, ++Y, ++Y1) { *Y = (double)Y1->r; *++Y = (double)Y1->i; }
                    Y1 -= Ly;
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { X1->r = (kiss_fft_scalar)*X; X1->i = (kiss_fft_scalar)*(X+1); }
                        X1 -= Lx;
                        kiss_fft(kiss_fft_state,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, Y+=2u*K, ++Y1) { *Y = (double)Y1->r; *(Y+1) = (double)Y1->i; }
                        Y1 -= Ly;
                    }
                }
            }
        }
        
        //Free
        free(kiss_fft_state); free(X1); free(Y1);
    }

    //Scale
    if (sc)
    {
        const double s = 1.0/sqrt((double)(2u*nfft));
        for (size_t l=2u*nfft*N/Lx; l>0u; --l, ++Y) { *Y *= s; }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
