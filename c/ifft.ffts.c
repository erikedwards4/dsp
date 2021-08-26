//Does 1-D IFFT (inverse fast Fourier transform) of each vector in X along dim.
//The input X is complex-valued.
//The output Y has the same size as X, except along dim, where Y has length nfft.
//Y is real-valued for ifft_fftw_s and ifft_fftw_d,
//and complex-valued for ifft_fftw_c and ifft_fftw_z.
//In the former case, X has only nonnegative freqs, so Lx = nfft/2 + 1.

//If sc, then scales Y by 2.0*sqrt(2.0)/nfft, so that invertible with ifft.

//This uses the fully-open-source FFTS (Fastest FFT in the South) library
//https://github.com/anthonix/ffts

//I installed into their default: /usr/local/include/ffts.
//So, either use that, change the #include, and/or use the -I flag.

//Double-precision is not supported! At least not under same install as single.
//So, I cast to float for X1, Y1 to support this, and issue a warning.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/usr/local/include/ffts/ffts.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ifft_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);


int ifft_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_ffts_s: dim must be in [0 3]\n"); return 1; }
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float s = (sc) ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_fftw_s: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (nfft%2u)
    {
        fprintf(stderr,"error in ifft_ffts_s: nfft must be even\n"); return 1;
    }
    else if (nfft & (nfft-1u)) //not a power-of-2
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in ifft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *iplan = ffts_init_1d(nfft,FFTS_BACKWARD);
        if (!iplan) { fprintf(stderr,"error in ifft_ffts_s: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;
        
        if (Lx==N)
        {
            for (size_t l=Lx; l>1u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
            *X1++ = *X; *X1++ = 0.0f; X -= 2u;
            for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
            X1 -= 2u*nfft;
            ffts_execute(iplan,X1,Y1);
            for (size_t l=nfft; l>0u; --l, Y1+=2, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y1-=2u*nfft)
                {
                    for (size_t l=Lx; l>1u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
                    *X1++ = *X; *X1++ = 0.0f; X -= 2u;
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                    X1 -= 2u*nfft;
                    ffts_execute(iplan,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, Y1+=2, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y1-=2u*nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>1u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1u); }
                        *X1++ = *X; *X1++ = 0.0f; X -= 2u*K;
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                        X1 -= 2u*nfft;
                        ffts_execute(iplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, Y1+=2, Y+=K) { *Y = *Y1 * s; }
                    }
                }
                Y -= G*B*nfft;
            }
        }
        
        //Free
        ffts_free(iplan); free(X1); free(Y1);
    }
    else    //nfft is a power-of-2
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in ifft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *iplan = ffts_init_1d_real(nfft,FFTS_BACKWARD);
        if (!iplan) { fprintf(stderr,"error in ifft_ffts_s: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=Lx; l>1u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
            *X1++ = *X; *X1++ = 0.0f; X -= 2u;
            for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
            X1 -= 2u*nfft;
            ffts_execute(iplan,X1,Y1);
            for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
            Y1 -= nfft; Y -= nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=nfft)
                {
                    for (size_t l=Lx; l>1u; --l, ++X, ++X1) { *X1 = *X; *++X1 = *++X; }
                    *X1++ = *X; *X1++ = 0.0f; X -= 2u;
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                    X1 -= 2u*nfft;
                    ffts_execute(iplan,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
                Y -= V*nfft;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y1-=nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>1u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                        *X1++ = *X; *X1++ = 0.0f; X -= 2u*K;
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { *X1 = *X; *++X1 = -*(X+1u); }
                        X1 -= 2u*nfft;
                        ffts_execute(iplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=K) { *Y = *Y1 * s; }
                    }
                }
                Y -= G*B*nfft;
            }
        }
        
        //Free
        ffts_free(iplan); free(X1); free(Y1);
    }

    return 0;
}


int ifft_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_ffts_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double s = (sc) ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_ffts_d: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in ifft_ffts_d: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (nfft%2u)
    {
        fprintf(stderr,"error in ifft_ffts_d: nfft must be even\n"); return 1;
    }
    else if (nfft & (nfft-1u)) //not a power-of-2
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in ifft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *iplan = ffts_init_1d(nfft,FFTS_BACKWARD);
        if (!iplan) { fprintf(stderr,"error in ifft_ffts_d: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=Lx; l>1u; --l, X+=2, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1); }
            *X1++ = (float)*X; *X1++ = 0.0; X -= 2u;
            for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = (float)*X; *++X1 = (float)(-*(X+1u)); }
            X1 -= 2u*nfft;
            ffts_execute(iplan,X1,Y1);
            for (size_t l=nfft; l>0u; --l, Y1+=2, ++Y) { *Y = (double)*Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y1-=2u*nfft)
                {
                    for (size_t l=Lx; l>1u; --l, X+=2, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1); }
                    *X1++ = (float)*X; *X1++ = 0.0; X -= 2u;
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = (float)*X; *++X1 = (float)(-*(X+1u)); }
                    X1 -= 2u*nfft;
                    ffts_execute(iplan,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, Y1+=2, ++Y) { *Y = (double)*Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2, Y1-=2u*nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>1u; --l, X+=2u*K, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1u); }
                        *X1++ = (float)*X; *X1++ = 0.0; X -= 2u*K;
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { *X1 = (float)*X; *++X1 = (float)(-*(X+1u)); }
                        X1 -= 2u*nfft;
                        ffts_execute(iplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, Y1+=2u, Y+=K) { *Y = (double)*Y1 * s; }
                    }
                }
            }
        }
        
        //Free
        ffts_free(iplan); free(X1); free(Y1);
    }
    else    //nfft is a power-of-2
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in ifft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *iplan = ffts_init_1d(nfft,FFTS_BACKWARD);
        if (!iplan) { fprintf(stderr,"error in ifft_ffts_d: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=Lx; l>1u; --l, X+=2, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1); }
            *X1++ = (float)*X; *X1++ = 0.0; X -= 2u;
            for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = (float)*X; *++X1 = (float)(-*(X+1u)); }
            X1 -= 2u*nfft;
            ffts_execute(iplan,X1,Y1);
            for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1 * s; }
            Y1 -= nfft; Y -= nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=nfft)
                {
                    for (size_t l=Lx; l>1u; --l, X+=2, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1); }
                    *X1++ = (float)*X; *X1++ = 0.0; X -= 2u;
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2, ++X1) { *X1 = (float)*X; *++X1 = (float)(-*(X+1u)); }
                    X1 -= 2u*nfft;
                    ffts_execute(iplan,X1,Y1);
                    for (size_t l=nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1 * s; }
                }
                Y -= V*nfft;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y1-=nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>1u; --l, X+=2u*K, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1u); }
                        *X1++ = (float)*X; *X1++ = 0.0; X -= 2u*K;
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++X1) { *X1 = (float)*X; *++X1 = (float)(-*(X+1u)); }
                        X1 -= 2u*nfft;
                        ffts_execute(iplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=K) { *Y = (double)*Y1 * s; }
                    }
                }
                Y -= G*B*nfft;
            }
        }
        
        //Free
        ffts_free(iplan); free(X1); free(Y1);
    }

    return 0;
}


int ifft_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_ffts_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float s = (sc) ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (nfft<Lx) { fprintf(stderr,"error in ifft_ffts_c: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (nfft%2u) { fprintf(stderr,"error in ifft_ffts_z: nfft must be even\n"); return 1; }
    else
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in ifft_ffts_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_ffts_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *iplan = ffts_init_1d(nfft,FFTS_BACKWARD);
        if (!iplan) { fprintf(stderr,"error in ifft_ffts_c: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*nfft;
            ffts_execute(iplan,X1,Y1);
            for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=2u*nfft)
                {
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= 2u*nfft;
                    ffts_execute(iplan,X1,Y1);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1u); }
                        X1 -= 2u*nfft;
                        ffts_execute(iplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K-1u) { *Y = *Y1 * s; *++Y = *++Y1 * s; }
                    }
                }
            }
        }
        ffts_free(iplan); free(X1); free(Y1);
    }

    return 0;
}


int ifft_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_ffts_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double s = (sc) ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (nfft<Lx) { fprintf(stderr,"error in ifft_ffts_z: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in ifft_ffts_z: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else if (nfft%2u) { fprintf(stderr,"error in ifft_ffts_z: nfft must be even\n"); return 1; }
    else
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in ifft_ffts_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in ifft_ffts_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *iplan = ffts_init_1d(nfft,FFTS_BACKWARD);
        if (!iplan) { fprintf(stderr,"error in ifft_ffts_z: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = (float)*X; }
            X1 -= 2u*nfft;
            ffts_execute(iplan,X1,Y1);
            for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=2u*nfft)
                {
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = (float)*X; }
                    X1 -= 2u*nfft;
                    ffts_execute(iplan,X1,Y1);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1u); }
                        X1 -= 2u*nfft;
                        ffts_execute(iplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, Y1+=2u, Y+=2u*K) { *Y = (double)*Y1 * s; *(Y+1) = (double)*(Y1+1) * s; }
                    }
                }
            }
        }
        ffts_free(iplan); free(X1); free(Y1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
