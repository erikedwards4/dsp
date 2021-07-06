//Does 1-D FFT (fast Fourier transform) of each vector in X along dim.
//The output Y is complex-valued and has the same size as X,
//except along dim, where Y has length Ly = nfft/2 + 1 for real-valued X,
//and length Ly = nfft for complex-valued X.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

//This uses the algorithm from Ch. 20 of Introduction to Algorithms [Cormen et al.].

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <time.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_PIf
    #define M_PIf 3.14159265358979323846f
#endif

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

void get_bittbl(size_t* bitrev, const size_t nfft);
void get_sintbl_s (float* sintbl, const size_t nfft);
void get_sintbl_d (double* sintbl, const size_t nfft);
int fft_algo_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
//int fft_algo_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);


void get_bittbl(size_t* bittbl, const size_t nfft)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);

    //Direct from Wikipedia article (https://en.wikipedia.org/wiki/Bit-reversal_permutation)
    // size_t K=0u, P2=1u;
    // while (P2<nfft) { P2 *= 2u; ++K; }
    // bittbl[0u] = 0u;
    // size_t L = 1u;
    // for (size_t k=0u; k<K; ++k, L*=2u)
    // {
    //     for (size_t l=1u; l<L; ++l) { bittbl[l] *= 2u; }
    //     for (size_t l=0u; l<L; ++l) { bittbl[l+L] = bittbl[l] + 1u; }
    // }

    //Other solution
    const size_t nfft2 = nfft/2u;
    size_t j=0u, k;
    *bittbl++ = 0u;
    for (size_t i=1u; i<nfft; ++i, ++bittbl)
    {
        k = nfft2;
        while (k<=j) { j -= k; k /= 2u; }
        j += k;
        *bittbl = j;
    }
    bittbl -= nfft;

    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
}


void get_sintbl_s (float* sintbl, const size_t nfft)
{
    const size_t nfft2=nfft/2u, nfft4=nfft/4u, nfft8=nfft/8u;
    float c=1.0f, s=0.0f, dc, ds, t;

    t = (float)(sin(M_PI/(double)nfft));
    dc = 2.0f * t * t;
    t = 2.0f * dc;
    ds = sqrtf(t-dc*dc);

    for (size_t i=0u; i<nfft8; ++i, ++sintbl, s+=ds, ds-=t*s) { *sintbl = s; }
    if (nfft8>0u) { *sintbl = M_SQRT1_2f; }
    sintbl += nfft4 - nfft8;
    for (size_t i=0u; i<nfft8; ++i, --sintbl, c-=dc, dc+=t*c) { *sintbl = c; }
    sintbl -= nfft4 - nfft8;
    for (size_t i=0u; i<nfft4; ++i) { sintbl[nfft2-i] = sintbl[i]; }
    for (size_t i=0u; i<nfft2+nfft4; ++i) { sintbl[i+nfft2] = -sintbl[i]; }
}


void get_sintbl_d (double* sintbl, const size_t nfft)
{
    const size_t nfft2=nfft/2u, nfft4=nfft/4u, nfft8=nfft/8u;
    double c=1.0, s=0.0, dc, ds, t;

    t = sin(M_PI/(double)nfft);
    dc = 2.0 * t * t;
    t = 2.0 * dc;
    ds = sqrt(t-dc*dc);

    for (size_t i=0u; i<nfft8; ++i, ++sintbl, s+=ds, ds-=t*s) { *sintbl = s; }
    if (nfft8>0u) { *sintbl = M_SQRT1_2; }
    sintbl += nfft4 - nfft8;
    for (size_t i=0u; i<nfft8; ++i, --sintbl, c-=dc, dc+=t*c) { *sintbl = c; }
    sintbl -= nfft4 - nfft8;
    for (size_t i=0u; i<nfft4; ++i) { sintbl[nfft2-i] = sintbl[i]; }
    for (size_t i=0u; i<nfft2+nfft4; ++i) { sintbl[i+nfft2] = -sintbl[i]; }
}


void fft_1d_s (float *X, const size_t N, const size_t *bittbl, const float *sintbl)
{
    const size_t N4 = N/4u;
    size_t m, d, ii, kk, ik;
    float x, s, c, dx, dy;

    //Bit reversal
    for (size_t n=0u; n<N; ++n, X+=2u)
    {
        m = bittbl[n];
        if (n<m) { x = *X; *X = *(X+2u*(m-n)); *(X+2u*(m-n)) = x; }
    }
    X -= 2u*N;

    //Transformation
    for (size_t k=1u; k<N; sintbl-=k*d, k=kk)
    {
        kk=k+k; d=N/kk;
        for (size_t j=0u; j<k; ++j, sintbl+=d)
        {
            s = *sintbl; c = *(sintbl+N4);
            for (size_t i=j; i<N; i+=kk)
            {
                ii = i+i; ik = ii + kk;
                dx = s*X[ik+1u] + c*X[ik];
                dy = c*X[ik+1u] - s*X[ik];
                X[ik] = X[ii] - dx;
                X[ii] += dx;
                X[ik+1u] = X[ii+1u] - dy;
                X[ii+1u] += dy;
                fprintf(stderr,"k=%lu, j=%lu, i=%lu, ik=%lu, ii=%lu \n",k,j,i,ik,ii);
            }
        }
    }
}


int fft_algo_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_algo_s: dim must be in [0 3]\n"); return 1; }
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft;
    if (nfft<Lx) { fprintf(stderr,"error in fft_algo_s: nfft must be >= L (vec length)\n"); return 1; }

    //Initialize fft
    float *X1;
    size_t bittbl[nfft];
    float sintbl[nfft+nfft/4u];
    get_bittbl(bittbl,nfft);
    get_sintbl_s(sintbl,nfft);
    X1 = (float *)malloc(2u*nfft*sizeof(float));

    if (N==0u) {}
    else if (Lx==1u && Ly==1u)
    {
        if (sc)
        {
            for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X * M_SQRT1_2f; *++Y = 0.0f; }
        }
        else
        {
            for (size_t n=0u; n<N; ++n, ++X, ++Y) { *Y = *X; *++Y = 0.0f; }
        }
    }
    else if (Lx==N)
    {
        for (size_t l=0u; l<Lx; ++l, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        for (size_t l=Lx; l<Ly; ++l) { *Y++ = 0.0f; *Y++ = 0.0f; }
        Y -= 2u*Ly;
        fft_1d_s(Y,nfft,bittbl,sintbl);
    }
    else
    {
        // const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        // const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        // const size_t V = N/Lx, G = V/B;

        // if (K==1u && (G==1u || B==1u))
        // {
        //     for (size_t l=0u; l<Lx; ++l, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        //     for (size_t l=Lx; l<Ly; ++l) { *Y++ = 0.0f; *Y++ = 0.0f; }
        //     Y -= 2u*Ly;
        //     fft_1d_s(Y,nfft,bittbl,sintbl);

        //     if (sc)
        //     {
        //         const float s = (float)(1.0/sqrt(2*Lx));
        //         for (size_t l=0u; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        //         Y1 -= 2u*Ly;
        //         for (size_t v=1; v<V; ++v, Y1-=2u*Ly)
        //         {
        //             for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
        //             X1 -= Lx;
        //             //fftwf_execute(plan);
        //             for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
        //         }
        //     }
        //     else
        //     {
        //         for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
        //         Y1 -= 2u*Ly;
        //         for (size_t v=1u; v<V; ++v, Y1-=2u*Ly)
        //         {
        //             for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
        //             X1 -= Lx;
        //             //fftwf_execute(plan);
        //             for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
        //         }
        //     }
        // }
        // else
        // {
        //     const float s = sc ? (float)(1.0/sqrt(2*Lx)) : 1.0f;
        //     X1 += Lx;
        //     for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
        //     X1 -= nfft;
        //     for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
        //     {
        //         for (size_t b=0; b<B; ++b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
        //         {
        //             for (size_t l=0u; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
        //             X1 -= Lx;
        //             //fftwf_execute(plan);
        //             for (size_t l=0u; l<Ly; ++l, ++Y1, Y+=2u*K-1u) { *Y = *Y1*s; *++Y = *++Y1*s; }
        //         }
        //     }
        // }
    }
    
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


// int fft_algo_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
// {
//     if (dim>3u) { fprintf(stderr,"error in fft_algo_d: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     const size_t Ly = nfft/2u + 1u;
//     if (nfft<Lx) { fprintf(stderr,"error in fft_algo_d: nfft must be >= L (vec length)\n"); return 1; }

//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
