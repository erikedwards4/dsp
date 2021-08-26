//Does 1-D FFT (fast Fourier transform) of each vector in X along dim.
//The output Y is complex-valued and has the same size as X,
//except along dim, where Y has length Ly = nfft/2 + 1 for real-valued X,
//and length Ly = nfft for complex-valued X.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

//This uses the fully-open-source FFTS (Fastest FFT in the South) library
//https://github.com/anthonix/ffts

//I installed into their default: /usr/local/include/ffts.
//So, either use that, change the #include, and/or use the -I flag.

//This definitely requires allocating X1 and copying X into X1.
//However, this appears to work fine without a separate Y1.
//But the proper way is to use Y1, and I haven't done extensive tests, so using Y1 for now.

//Odd-length FFT not supported. Even-length, non-power-of-2 supported, but must make X1 complex.

//Double-precision is not supported! At least not under same install as single.
//So, I cast to float for X1, Y1 to support this, and issue a warning.

//Their memory allocation is very instructional!
//Since aligned_alloc is the most preferred, and is in <stdlib.h>, and tests to approx. same speed as malloc, I hard-code it for now.
//To use _mm_malloc, need to include <xmmintrin.h>, use signed ints as input, use _mm_free, and compile with flags -msse -mmmx -msse2.
//posix_memalign is also in <stdlib.h>, but is harder to use; and memalign requires include <malloc.h>.
//The next 2 are obscure, valloc is "obsolete", and finally malloc is lowest-ranked.
//This is their list of methods, in order of preference:
//#if defined(HAVE_ALIGNED_ALLOC)
//     p = aligned_alloc(32, size);
// #elif defined(__ICC) || defined(__INTEL_COMPILER) || defined(HAVE__MM_MALLOC)
//     p = (void*) _mm_malloc(size, 32);
// #elif defined(HAVE_POSIX_MEMALIGN)
//     if (posix_memalign(&p, 32, size)) { p = NULL; }
// #elif defined(HAVE_MEMALIGN)
//     p = memalign(32, size);
// #elif defined(__ALTIVEC__)
//     p = vec_malloc(size);
// #elif defined(_MSC_VER) || defined(WIN32)
//     p = _aligned_malloc(size, 32);
// #elif defined(HAVE_VALLOC)
//     p = valloc(size);
// #else
//     p = malloc(size);
// #endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "/usr/local/include/ffts/ffts.h"
//#include <time.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int fft_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int fft_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);


int fft_ffts_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_ffts_s: dim must be in [0 3]\n"); return 1; }
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft/2u + 1u;
    if (nfft<Lx) { fprintf(stderr,"error in fft_ffts_s: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        Y -= 2u*N;
    }
    else if (nfft%2u)
    {
        fprintf(stderr,"error in fft_ffts_s: nfft must be even\n"); return 1;
    }
    else if (nfft & (nfft-1u)) //not a power-of-2
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in fft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in fft_ffts_s: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;
        
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, X1+=2) { *X1 = *X; }
            X1 -= 2u*Lx;
            ffts_execute(fplan,X1,Y1);
            for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = *Y1; }
            Y1 -= 2u*Ly; Y -= 2u*Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=2u*Ly)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, X1+=2) { *X1 = *X; }
                    X1 -= 2u*Lx;
                    ffts_execute(fplan,X1,Y1);
                    for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = *Y1; }
                }
                Y -= 2u*Ly*V;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, X1+=2) { *X1 = *X; }
                        X1 -= 2u*Lx;
                        ffts_execute(fplan,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, ++Y1, Y+=2u*K) { *Y = *Y1; *(Y+1) = *++Y1; }
                    }
                }
                Y -= 2u*G*B*Ly;
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }
    else    //nfft is a power-of-2
    {
        //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
        
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in fft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_ffts_s: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d_real(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in fft_ffts_s: problem creating ffts plan"); return 1; }
        for (size_t n=nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= nfft;

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= Lx;
            ffts_execute(fplan,X1,Y1);
            for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = *Y1; }
            Y1 -= 2u*Ly; Y -= 2u*Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=2u*Ly)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
                    X1 -= Lx;
                    ffts_execute(fplan,X1,Y1);
                    for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = *Y1; }
                }
                Y -= 2u*Ly*V;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, ++X1) { *X1 = *X; }
                        X1 -= Lx;
                        ffts_execute(fplan,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, ++Y1, Y+=2u*K) { *Y = *Y1; *(Y+1) = *++Y1; }
                    }
                }
                Y -= 2u*G*B*Ly;
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);

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


int fft_ffts_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_ffts_d: dim must be in [0 3]\n"); return 1; }    

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft/2u + 1u;
    if (nfft<Lx) { fprintf(stderr,"error in fft_ffts_d: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in fft_ffts_d: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X) { *Y++ = *X; *Y++ = 0.0; }
        Y -= 2u*N;
    }
    else if (nfft%2u)
    {
        fprintf(stderr,"error in fft_ffts_d: nfft must be even\n"); return 1;
    }
    else if (nfft & (nfft-1u)) //not a power-of-2
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in fft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in fft_ffts_d: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, X1+=2) { *X1 = (float)*X; }
            X1 -= 2u*Lx;
            ffts_execute(fplan,X1,Y1);
            for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1; }
            Y1 -= 2u*Ly; Y -= 2u*Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=2u*Ly)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, X1+=2) { *X1 = (float)*X; }
                    X1 -= 2u*Lx;
                    ffts_execute(fplan,X1,Y1);
                    for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1; }
                }
                Y -= 2u*Ly*V;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, X1+=2) { *X1 = (float)*X; }
                        X1 -= 2u*Lx;
                        ffts_execute(fplan,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, ++Y1, Y+=2u*K) { *Y = (double)*Y1; *(Y+1) = (double)*++Y1; }
                    }
                }
                Y -= 2u*G*B*Ly;
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }
    else    //nfft is a power-of-2
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in fft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_ffts_d: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d_real(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in fft_ffts_d: problem creating ffts plan"); return 1; }
        for (size_t n=nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= nfft;

        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = (float)*X; }
            X1 -= Lx;
            ffts_execute(fplan,X1,Y1);
            for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1; }
            Y1 -= 2u*Ly; Y -= 2u*Ly;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y1-=2u*Ly)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++X1) { *X1 = (float)*X; }
                    X1 -= Lx;
                    ffts_execute(fplan,X1,Y1);
                    for (size_t l=2u*Ly; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1; }
                }
                Y -= 2u*Ly*V;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=K, ++X1) { *X1 = (float)*X; }
                        X1 -= Lx;
                        ffts_execute(fplan,X1,Y1);
                        for (size_t l=Ly; l>0u; --l, ++Y1, Y+=2u*K) { *Y = (double)*Y1; *(Y+1) = (double)*++Y1; }
                    }
                }
                Y -= 2u*G*B*Ly;
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }

    //Scale
    if (sc)
    {
        const double s = 1.0/sqrt((double)(2u*nfft));
        for (size_t l=2u*Ly*N/Lx; l>0u; --l, ++Y) { *Y *= s; }
    }

    return 0;
}


int fft_ffts_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_ffts_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (nfft<Lx) { fprintf(stderr,"error in fft_ffts_c: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
    }
    else if (nfft%2u) { fprintf(stderr,"error in fft_ffts_c: nfft must be even\n"); return 1; }
    else
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in fft_ffts_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_ffts_c: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in fft_ffts_c: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = *X; }
            X1 -= 2u*Lx;
            ffts_execute(fplan,X1,Y1);
            for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1; }
            Y1 -= 2u*nfft; Y -= 2u*nfft;
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
                    X1 -= 2u*Lx;
                    ffts_execute(fplan,X1,Y1);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = *Y1; }
                }
                Y -= 2u*nfft*V;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = *X; *++X1 = *(X+1); }
                        X1 -= 2u*Lx;
                        ffts_execute(fplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K) { *Y = *Y1; *(Y+1) = *++Y1; }
                    }
                }
                Y -= 2u*G*B*nfft;
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
    }

    //Scale
    if (sc)
    {
        const float s = 1.0f/sqrtf((float)(2u*nfft));
        for (size_t l=2u*nfft*N/Lx; l>0u; --l, ++Y) { *Y *= s; }
    }

    return 0;
}


int fft_ffts_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_ffts_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (nfft<Lx) { fprintf(stderr,"error in fft_ffts_z: nfft must be >= Lx (length of vecs in X)\n"); return 1; }

    //Double not directly supported!
    fprintf(stderr,"warning in fft_ffts_z: double precision not directly supported, casting to float for FFT part\n");

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X; }
        Y -= 2u*N;
    }
    else if (nfft%2u) { fprintf(stderr,"error in fft_ffts_z: nfft must be even\n"); return 1; }
    else
    {
        //Initialize FFT
        float *X1, *Y1;
        X1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        Y1 = (float *)aligned_alloc(sizeof(float),2u*nfft*sizeof(float));
        if (!X1) { fprintf(stderr,"error in fft_ffts_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        if (!Y1) { fprintf(stderr,"error in fft_ffts_z: problem with aligned_alloc. "); perror("aligned_alloc"); return 1; }
        ffts_plan_t *fplan = ffts_init_1d(nfft,FFTS_FORWARD);
        if (!fplan) { fprintf(stderr,"error in fft_ffts_z: problem creating ffts plan"); return 1; }
        for (size_t n=2u*nfft; n>0u; --n, ++X1) { *X1 = 0.0f; }
        X1 -= 2u*nfft;

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++X1) { *X1 = (float)*X; }
            X1 -= 2u*Lx;
            ffts_execute(fplan,X1,Y1);
            for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1; }
            Y1 -= 2u*nfft; Y -= 2u*nfft;
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
                    X1 -= 2u*Lx;
                    ffts_execute(fplan,X1,Y1);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y1, ++Y) { *Y = (double)*Y1; }
                }
                Y -= 2u*nfft*V;
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++X1) { *X1 = (float)*X; *++X1 = (float)*(X+1); }
                        X1 -= 2u*Lx;
                        ffts_execute(fplan,X1,Y1);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K) { *Y = (double)*Y1; *(Y+1) = (double)*++Y1; }
                    }
                }
                Y -= 2u*G*B*nfft;
            }
        }
        
        //Free
        ffts_free(fplan); free(X1); free(Y1);
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
