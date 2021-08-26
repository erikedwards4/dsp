//This computes "the" DCT (type-II discrete cosine transformation) along dim of matrix X.
//This uses FFTS to compute the 4u*ndct length FFT of each vec in X.

//For complex input X, the output Y is complex, and the DCT is just the DCT of the
//real and imaginary parts separately (following Octave convention).

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/usr/local/include/ffts/ffts.h"
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int dct_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);
int dct_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc);


int dct_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_ffts_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_ffts_s: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const float s = sc ? 1.0f/sqrtf((float)(2u*ndct)) : 1.0f;
        const float dcsc = sc ? 0.5f/sqrtf((float)ndct) : 1.0f;

        //Initialize FFT
        const size_t nfft = 4u*ndct;
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in dct_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in dct_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in dct_ffts_s: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;
    
        if (Lx==N)
        {
            X1 += 2u;
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *X++; }
            X1 += 8u*(ndct-Lx);
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *--X; }
            X1 -= 2u + 8u*ndct;
            ffts_execute(fplan,X1,Y1);
            if (sc)
            {
                *Y = *Y1 * dcsc; ++Y; Y1 += 2;
                for (size_t n=ndct; n>1u; --n, ++Y, Y1+=2u) { *Y = *Y1 * s; }
            }
            else
            {
                for (size_t n=ndct; n>0u; --n, ++Y, Y1+=2u) { *Y = *Y1; }
            }
            Y1 -= 2u*ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Lx, Y1-=2u*ndct)
                {
                    X1 += 2u;
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *X++; }
                    X1 += 8u*(ndct-Lx);
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *--X; }
                    X1 -= 2u + 8u*ndct;
                    ffts_execute(fplan,X1,Y1);
                    if (sc)
                    {
                        *Y = *Y1 * dcsc; ++Y; Y1 += 2;
                        for (size_t n=ndct; n>1u; --n, ++Y, Y1+=2u) { *Y = *Y1 * s; }
                    }
                    else
                    {
                        for (size_t n=ndct; n>0u; --n, ++Y, Y1+=2u) { *Y = *Y1; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndct-1u, Y1-=2u*ndct)
                    {
                        X1 += 2u;
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *X; X+=K; }
                        X1 += 8u*(ndct-Lx);
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { X-=K; *X1 = *X; }
                        X1 -= 2u + 8u*ndct;
                        ffts_execute(fplan,X1,Y1);
                        if (sc)
                        {
                            *Y = *Y1 * dcsc; Y+=K; Y1 += 2;
                            for (size_t n=ndct; n>1u; --n, Y+=K, Y1+=2u) { *Y = *Y1 * s; }
                        }
                        else
                        {
                            for (size_t n=ndct; n>0u; --n, Y+=K, Y1+=2u) { *Y = *Y1; }
                        }
                    }
                }
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }

    return 0;
}


int dct_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_ffts_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_ffts_d: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const double s = sc ? 1.0/sqrt((double)(2u*ndct)) : 1.0;
        const double dcsc = sc ? 0.5/sqrt((double)ndct) : 1.0;

        //Initialize FFT
        const size_t nfft = 4u*ndct;
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in dct_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in dct_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in dct_ffts_d: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;
    
        if (Lx==N)
        {
            X1 += 2u;
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*X++; }
            X1 += 8u*(ndct-Lx);
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*--X; }
            X1 -= 2u + 8u*ndct;
            ffts_execute(fplan,X1,Y1);
            if (sc)
            {
                *Y = (double)*Y1 * dcsc; ++Y; Y1 += 2;
                for (size_t n=ndct; n>1u; --n, ++Y, Y1+=2u) { *Y = (double)*Y1 * s; }
            }
            else
            {
                for (size_t n=ndct; n>0u; --n, ++Y, Y1+=2u) { *Y = (double)*Y1; }
            }
            Y1 -= 2u*ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=Lx, Y1-=2u*ndct)
                {
                    X1 += 2u;
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*X++; }
                    X1 += 8u*(ndct-Lx);
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*--X; }
                    X1 -= 2u + 8u*ndct;
                    ffts_execute(fplan,X1,Y1);
                    if (sc)
                    {
                        *Y = (double)*Y1 * dcsc; ++Y; Y1 += 2;
                        for (size_t n=ndct; n>1u; --n, ++Y, Y1+=2u) { *Y = (double)*Y1 * s; }
                    }
                    else
                    {
                        for (size_t n=ndct; n>0u; --n, ++Y, Y1+=2u) { *Y = (double)*Y1; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, ++X, Y-=K*ndct-1u, Y1-=2u*ndct)
                    {
                        X1 += 2u;
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*X; X+=K; }
                        X1 += 8u*(ndct-Lx);
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { X-=K; *X1 = (float)*X; }
                        X1 -= 2u + 8u*ndct;
                        ffts_execute(fplan,X1,Y1);
                        if (sc)
                        {
                            *Y = (double)*Y1 * dcsc; Y+=K; Y1 += 2;
                            for (size_t n=ndct; n>1u; --n, Y+=K, Y1+=2u) { *Y = (double)*Y1 * s; }
                        }
                        else
                        {
                            for (size_t n=ndct; n>0u; --n, Y+=K, Y1+=2u) { *Y = (double)*Y1; }
                        }
                    }
                }
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }

    return 0;
}


int dct_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_ffts_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_ffts_c: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const float s = sc ? 1.0f/sqrtf((float)(2u*ndct)) : 1.0f;
        const float dcsc = sc ? 0.5f/sqrtf((float)ndct) : 1.0f;

        //Initialize FFT
        const size_t nfft = 4u*ndct;
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in dct_ffts_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in dct_ffts_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in dct_ffts_c: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;
    
        if (Lx==N)
        {
            X1 += 2u;
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *X++; *(X1+1) = *X++; }
            X1 += 8u*(ndct-Lx);
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *(X1+1) = *--X; *X1 = *--X; }
            X1 -= 2u + 8u*ndct;
            ffts_execute(fplan,X1,Y1);
            if (sc)
            {
                *Y = *Y1 * dcsc; ++Y; ++Y1;
                *Y = *Y1 * dcsc; ++Y; ++Y1;
                for (size_t n=2u*ndct; n>2u; --n, ++Y, ++Y1) { *Y = *Y1 * s; }
            }
            else
            {
                for (size_t n=2u*ndct; n>0u; --n, ++Y, ++Y1) { *Y = *Y1; }
            }
            Y1 -= 2u*ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y1-=2u*ndct)
                {
                    X1 += 2u;
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *X++; *(X1+1) = *X++; }
                    X1 += 8u*(ndct-Lx);
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *(X1+1) = *--X; *X1 = *--X; }
                    X1 -= 2u + 8u*ndct;
                    ffts_execute(fplan,X1,Y1);
                    if (sc)
                    {
                        *Y = *Y1 * dcsc; ++Y; ++Y1;
                        *Y = *Y1 * dcsc; ++Y; ++Y1;
                        for (size_t n=2u*ndct; n>2u; --n, ++Y, ++Y1) { *Y = *Y1 * s; }
                    }
                    else
                    {
                        for (size_t n=2u*ndct; n>0u; --n, ++Y, ++Y1) { *Y = *Y1; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*ndct-2u, Y1-=2u*ndct)
                    {
                        X1 += 2u;
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = *X; *(X1+1) = *(X+1); X+=2u*K; }
                        X1 += 8u*(ndct-Lx);
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { X-=2u*K; *X1 = *X; *(X1+1) = *(X+1); }
                        X1 -= 2u + 8u*ndct;
                        ffts_execute(fplan,X1,Y1);
                        if (sc)
                        {
                            *Y = *Y1 * dcsc; ++Y1;
                            *(Y+1) = *Y1 * dcsc; Y+=2u*K; ++Y1;
                            for (size_t n=ndct; n>1u; --n, Y+=2u*K, ++Y1) { *Y = *Y1 * s; *(Y+1) = *++Y1 * s; }
                        }
                        else
                        {
                            for (size_t n=ndct; n>0u; --n, Y+=2u*K, ++Y1) { *Y = *Y1; *(Y+1) = *++Y1; }
                        }
                    }
                }
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }

    return 0;
}


int dct_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t ndct, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in dct_ffts_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (ndct<Lx) { fprintf(stderr,"error in dct_ffts_z: ndct must be >= Lx (length of vecs in X)\n"); return 1; }

    if (ndct==0u || N==0u) {}
    else if (ndct==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const double s = sc ? 1.0/sqrt((double)(2u*ndct)) : 1.0;
        const double dcsc = sc ? 0.5/sqrt((double)ndct) : 1.0;

        //Initialize FFT
        const size_t nfft = 4u*ndct;
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in dct_ffts_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in dct_ffts_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in dct_ffts_z: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;
    
        if (Lx==N)
        {
            X1 += 2u;
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*X++; *(X1+1) = (float)*X++; }
            X1 += 8u*(ndct-Lx);
            for (size_t l=Lx; l>0u; --l, X1+=4u) { *(X1+1) = (float)*--X; *X1 = (float)*--X; }
            X1 -= 2u + 8u*ndct;
            ffts_execute(fplan,X1,Y1);
            if (sc)
            {
                *Y = (double)*Y1 * dcsc; ++Y; ++Y1;
                *Y = (double)*Y1 * dcsc; ++Y; ++Y1;
                for (size_t n=2u*ndct; n>2u; --n, ++Y, ++Y1) { *Y = (double)*Y1 * s; }
            }
            else
            {
                for (size_t n=2u*ndct; n>0u; --n, ++Y, ++Y1) { *Y = (double)*Y1; }
            }
            Y1 -= 2u*ndct;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, X+=2u*Lx, Y1-=2u*ndct)
                {
                    X1 += 2u;
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*X++; *(X1+1) = (float)*X++; }
                    X1 += 8u*(ndct-Lx);
                    for (size_t l=Lx; l>0u; --l, X1+=4u) { *(X1+1) = (float)*--X; *X1 = (float)*--X; }
                    X1 -= 2u + 8u*ndct;
                    ffts_execute(fplan,X1,Y1);
                    if (sc)
                    {
                        *Y = (double)*Y1 * dcsc; ++Y; ++Y1;
                        *Y = (double)*Y1 * dcsc; ++Y; ++Y1;
                        for (size_t n=2u*ndct; n>2u; --n, ++Y, ++Y1) { *Y = (double)*Y1 * s; }
                    }
                    else
                    {
                        for (size_t n=2u*ndct; n>0u; --n, ++Y, ++Y1) { *Y = (double)*Y1; }
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(ndct-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y-=2u*K*ndct-2u, Y1-=2u*ndct)
                    {
                        X1 += 2u;
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { *X1 = (float)*X; *(X1+1) = (float)*(X+1); X+=2u*K; }
                        X1 += 8u*(ndct-Lx);
                        for (size_t l=Lx; l>0u; --l, X1+=4u) { X-=2u*K; *X1 = (float)*X; *(X1+1) = (float)*(X+1); }
                        X1 -= 2u + 8u*ndct;
                        ffts_execute(fplan,X1,Y1);
                        if (sc)
                        {
                            *Y = (double)*Y1 * dcsc; ++Y1;
                            *(Y+1) = (double)*Y1 * dcsc; Y+=2u*K; ++Y1;
                            for (size_t n=ndct; n>1u; --n, Y+=2u*K, ++Y1) { *Y = (double)*Y1 * s; *(Y+1) = (double)*++Y1 * s; }
                        }
                        else
                        {
                            for (size_t n=ndct; n>0u; --n, Y+=2u*K, ++Y1) { *Y = (double)*Y1; *(Y+1) = (double)*++Y1; }
                        }
                    }
                }
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
