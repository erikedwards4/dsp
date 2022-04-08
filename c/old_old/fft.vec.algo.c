//Does 1-D FFT (fast Fourier transform) of a single vector in X.
//The output Y is complex-valued and has length Ly = nfft/2 + 1
//for real-valued X, and length Ly = nfft for complex-valued X.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

//This uses the algorithm from Ch. 20 of Introduction to Algorithms [Cormen et al.].

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

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
void get_cstbl_s (float* cstbl, const size_t nfft);
void get_cstbl_d (double* cstbl, const size_t nfft);
int fft_vec_algo_s (float *Y, const float *X, const size_t N, const size_t nfft, const char sc);
int fft_vec_algo_d (double *Y, const double *X, const size_t N, const size_t nfft, const char sc);
int fft_vec_algo_c (float *Y, const float *X, const size_t N, const size_t nfft, const char sc);
int fft_vec_algo_z (double *Y, const double *X, const size_t N, const size_t nfft, const char sc);

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

    //Other solution (approx. same speed, but shorter assembly code)
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


void get_cstbl_s (float* cstbl, const size_t nfft)
{
    const size_t nfft2=nfft/2u, nfft4=nfft/4u, nfft8=nfft/8u;
    float c=1.0f, s=0.0f, dc, ds, t;

    t = (float)(sin(M_PI/(double)nfft));
    dc = 2.0f * t * t;
    t = 2.0f * dc;
    ds = sqrtf(t-dc*dc);

    for (size_t i=0u; i<nfft8; ++i, ++cstbl, s+=ds, ds-=t*s) { *cstbl = s; }
    if (nfft8>0u) { *cstbl = M_SQRT1_2f; }
    cstbl += nfft4 - nfft8;
    for (size_t i=0u; i<nfft8; ++i, --cstbl, c-=dc, dc+=t*c) { *cstbl = c; }
    cstbl -= nfft4 - nfft8;
    for (size_t i=0u; i<nfft4; ++i) { cstbl[nfft2-i] = cstbl[i]; }
    for (size_t i=0u; i<nfft2+nfft4; ++i) { cstbl[i+nfft2] = -cstbl[i]; }
}


void get_cstbl_d (double* cstbl, const size_t nfft)
{
    const size_t nfft2=nfft/2u, nfft4=nfft/4u, nfft8=nfft/8u;
    double c=1.0, s=0.0, dc, ds, t;

    t = sin(M_PI/(double)nfft);
    dc = 2.0 * t * t;
    t = 2.0 * dc;
    ds = sqrt(t-dc*dc);

    for (size_t i=0u; i<nfft8; ++i, ++cstbl, s+=ds, ds-=t*s) { *cstbl = s; }
    if (nfft8>0u) { *cstbl = M_SQRT1_2; }
    cstbl += nfft4 - nfft8;
    for (size_t i=0u; i<nfft8; ++i, --cstbl, c-=dc, dc+=t*c) { *cstbl = c; }
    cstbl -= nfft4 - nfft8;
    for (size_t i=0u; i<nfft4; ++i) { cstbl[nfft2-i] = cstbl[i]; }
    for (size_t i=0u; i<nfft2+nfft4; ++i) { cstbl[i+nfft2] = -cstbl[i]; }
}


int fft_vec_algo_s (float *Y, const float *X, const size_t N, const size_t nfft, const char sc)
{
    if (nfft<N) { fprintf(stderr,"error in fft_vec_algo_s: nfft must be >= N (vec length)\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in fft_vec_algo_s: nfft must be a power of 2\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        if (sc) { *Y = *X * M_SQRT1_2f; *++Y = 0.0f; }
        else { *Y = *X; *++Y = 0.0f; }
    }
    else if (nfft==2u)
    {
        *Y++ = *X; *Y++ = 0.0f; *Y++ = *X; *Y-- = 0.0f;
        if (N==2u) { *Y -= *++X; *(Y-2u) += *X; }
        if (sc)
        {
            const float s = (N==1u) ? M_SQRT1_2f : 0.5f;
            *Y *= s; *(Y-2u) *= s;
        }
    }
    else if (nfft==4u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        for (size_t n=N; n<4u; ++n) { *Y++ = 0.0f; *Y++ = 0.0f; }
        Y -= 8u;
        const float dc = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const float ny = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        Y[3u] = Y[6u] - Y[2u];
        Y[2u] = Y[0u] - Y[4u];
        Y[0u] = dc; Y[4u] = ny;
        if (sc)
        {
            const float s = (float)(1.0/sqrt(2u*N));
            for (size_t n=0u; n<8u; ++n, ++Y) { *Y *= s; }
        }
    }
    else
    {
        const size_t nfft4 = nfft/4u, nfft2 = nfft/2u;
        size_t m=0u, h=0u, d, ii, kk, ik;
        float y, sr, si, dr, di;

        //Initialize fft
        float cstbl[nfft+nfft/4u];
        get_cstbl_s(cstbl,nfft);

        //Put X into real part of Y
        for (size_t n=0u; n<N; ++n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        for (size_t n=N; n<nfft; ++n) { *Y++ = 0.0f; *Y++ = 0.0f; }
        Y -= 2u*nfft-2u;

        //Bit reverse
        for (size_t n=1u; n<nfft; ++n, Y+=2u)
        {
            size_t k = nfft2;
            while (k<=m) { m -= k; k /= 2u; }
            m += k;
            if (n<m) { y = *Y; *Y = *(Y+2u*(m-n)); *(Y+2u*(m-n)) = y; }
        }
        Y -= 2u*nfft;

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = cstbl[h]; sr = cstbl[h+nfft4];
                for (size_t i=j; i<nfft; i+=kk)
                {
                    ii = i+i; ik = ii + kk;
                    dr = si*Y[ik+1u] + sr*Y[ik];
                    di = sr*Y[ik+1u] - si*Y[ik];
                    Y[ik] = Y[ii] - dr; Y[ii] += dr;
                    Y[ik+1u] = Y[ii+1u] - di; Y[ii+1u] += di;
                }
            }
        }

        //Scale
        if (sc)
        {
            const float s = M_SQRT1_2f/sqrtf(N);
            for (size_t n=0u; n<2u*nfft; ++n, ++Y) { *Y *= s; }
        }
    }
    
    return 0;
}


int fft_vec_algo_d (double *Y, const double *X, const size_t N, const size_t nfft, const char sc)
{
    if (nfft<N) { fprintf(stderr,"error in fft_vec_algo_d: nfft must be >= N (vec length)\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in fft_vec_algo_d: nfft must be a power of 2\n"); return 1; }

    if (nfft==0u) {}
    else if (nfft==1u)
    {
        if (sc) { *Y = *X * M_SQRT1_2; *++Y = 0.0; }
        else { *Y = *X; *++Y = 0.0; }
    }
    else if (nfft==2u)
    {
        *Y++ = *X; *Y++ = 0.0; *Y++ = *X; *Y-- = 0.0;
        if (N==2u) { *Y -= *++X; *(Y-2u) += *X; }
        if (sc)
        {
            const float s = (N==1u) ? M_SQRT1_2f : 0.5;
            *Y *= s; *(Y-2u) *= s;
        }
    }
    else if (nfft==4u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *Y++ = *X; *Y++ = 0.0; }
        for (size_t n=N; n<4u; ++n) { *Y++ = 0.0; *Y++ = 0.0; }
        Y -= 8u;
        const double dc = Y[0u] + Y[2u] + Y[4u] +Y[6u];
        const double ny = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        Y[3u] = Y[6u] - Y[2u];
        Y[2u] = Y[0u] - Y[4u];
        Y[0u] = dc; Y[4u] = ny;
        if (sc)
        {
            const double s = 1.0/sqrt(2u*N);
            for (size_t n=0u; n<8u; ++n, ++Y) { *Y *= s; }
        }
    }
    else
    {
        const size_t nfft4 = nfft/4u, nfft2 = nfft/2u;
        size_t m=0u, h=0u, d, ii, kk, ik;
        double y, sr, si, dr, di;

        //Initialize fft
        double cstbl[nfft+nfft/4u];
        get_cstbl_d(cstbl,nfft);

        //Put X into real part of Y
        for (size_t n=0u; n<N; ++n, ++X) { *Y++ = *X; *Y++ = 0.0; }
        for (size_t n=N; n<nfft; ++n) { *Y++ = 0.0; *Y++ = 0.0; }
        Y -= 2u*nfft-2u;

        //Bit reverse
        for (size_t n=1u; n<nfft; ++n, Y+=2u)
        {
            size_t k = nfft2;
            while (k<=m) { m -= k; k /= 2u; }
            m += k;
            if (n<m) { y = *Y; *Y = *(Y+2u*(m-n)); *(Y+2u*(m-n)) = y; }
        }
        Y -= 2u*nfft;

        //Transform
        // for (size_t k=1u; k<nfft; cstbl-=k*d, Y-=2u*k, k=kk)
        // {
        //     kk=k+k; d=nfft/kk;
        //     for (size_t j=0u; j<k; ++j, cstbl+=d)
        //     {
        //         si = *cstbl; sr = *(cstbl+nfft4);
        //         for (size_t i=j; i<nfft; i+=kk, Y+=2u*kk-2u)
        //         {
        //             dr = s**(Y+kk+1u) + c**(Y+kk);
        //             di = c**(Y+kk+1u) - s**(Y+kk);
        //             *(Y+kk) = *Y - dr; *Y++ += dr;
        //             *(Y+kk) = *Y - di; *Y++ += di;
        //         }
        //         Y -= 2u*kk*(1u+(nfft-j-1u)/kk) - 2u;
        //     }
        // }

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = cstbl[h]; sr = cstbl[h+nfft4];
                for (size_t i=j; i<nfft; i+=kk)
                {
                    ii = i+i; ik = ii + kk;
                    dr = si*Y[ik+1u] + sr*Y[ik];
                    di = sr*Y[ik+1u] - si*Y[ik];
                    Y[ik] = Y[ii] - dr; Y[ii] += dr;
                    Y[ik+1u] = Y[ii+1u] - di; Y[ii+1u] += di;
                }
            }
        }

        //Scale
        if (sc)
        {
            const double s = M_SQRT1_2/sqrt(N);
            for (size_t n=0u; n<2u*nfft; ++n, ++Y) { *Y *= s; }
        }
    }

    return 0;
}


int fft_vec_algo_c (float *Y, const float *X, const size_t N, const size_t nfft, const char sc)
{
    if (nfft<N) { fprintf(stderr,"error in fft_vec_algo_c: nfft must be >= N (vec length)\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in fft_vec_algo_c: nfft must be a power of 2\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        if (sc) { *Y = *X * M_SQRT1_2f; *++Y = *++X * M_SQRT1_2f; }
        else { *Y = *X; *++Y = *++X; }
    }
    else if (nfft==2u)
    {
        if (N==1u) { *Y++ = *X++; *Y++ = *X--; *Y++ = *X++; *Y = *X; }
        else
        {
            *Y++ = *X + *(X+2u); ++X; *Y++ = *X + *(X+2u); ++X;
            *Y++ = *(X-2u) - *X; ++X; *Y = *(X-2u) - *X;
        }
        if (sc)
        {
            const float s = (N==1u) ? M_SQRT1_2f : 0.5f;
            *Y-- *= s; *Y-- *= s; *Y-- *= s; *Y *= s;
        }
    }
    else if (nfft==4u)
    {
        for (size_t n=0u; n<N; ++n) { *Y++ = *X++; *Y++ = *X++; }
        for (size_t n=N; n<4u; ++n) { *Y++ = 0.0f; *Y++ = 0.0f; }
        Y -= 8u;
        const float dcr = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const float dci = Y[1u] + Y[3u] + Y[5u] + Y[7u];
        const float nyr = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        const float nyi = Y[1u] - Y[3u] + Y[5u] - Y[7u];
        const float y1r = Y[0u] - Y[4u] + Y[3u] - Y[7u];
        const float y1i = Y[1u] - Y[2u] + Y[6u] - Y[5u];
        const float y3r = Y[0u] - Y[3u] + Y[7u] - Y[4u];
        const float y3i = Y[1u] - Y[5u] + Y[2u] - Y[6u];
        Y[0u] = dcr; Y[1u] = dci; Y[2u] = y1r; Y[3u] = y1i;
        Y[4u] = nyr; Y[5u] = nyi; Y[6u] = y3r; Y[7u] = y3i;
        if (sc)
        {
            const float s = (float)(1.0/sqrt(2u*N));
            for (size_t n=0u; n<8u; ++n, ++Y) { *Y *= s; }
        }
    }
    else
    {
        const size_t nfft4 = nfft/4u;
        size_t m=0u, h=0u, d, ii, kk, ik;
        float y, sr, si, dr, di;

        //Initialize fft
        float cstbl[nfft+nfft4];
        get_cstbl_s(cstbl,nfft);

        //Copy X into Y
        for (size_t n=0u; n<N; ++n) { *Y++ = *X++; *Y++ = *X++; }
        for (size_t n=N; n<nfft; ++n) { *Y++ = 0.0f; *Y++ = 0.0f; }
        Y -= 2u*nfft;

        //Bit reverse
        Y += 2u;
        for (size_t n=1u; n<nfft; ++n, Y+=2u)
        {
            size_t k = nfft/2u;
            while (k<=m) { m -= k; k /= 2u; }
            m += k;
            if (n<m)
            {
                y = *Y; *Y = *(Y+2u*(m-n)); *(Y+2u*(m-n)) = y;
                y = *(Y+1u); *(Y+1u) = *(Y+1u+2u*(m-n)); *(Y+1u+2u*(m-n)) = y;
            }
        }
        Y -= 2u*nfft;

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = cstbl[h]; sr = cstbl[h+nfft4];
                for (size_t i=j; i<nfft; i+=kk)
                {
                    ii = i+i; ik = ii + kk;
                    dr = si*Y[ik+1u] + sr*Y[ik];
                    di = sr*Y[ik+1u] - si*Y[ik];
                    Y[ik] = Y[ii] - dr; Y[ii] += dr;
                    Y[ik+1u] = Y[ii+1u] - di; Y[ii+1u] += di;
                }
            }
        }

        //Scale
        if (sc)
        {
            const float s = M_SQRT1_2f/sqrtf(N);
            for (size_t n=0u; n<2u*nfft; ++n, ++Y) { *Y *= s; }
        }
    }
    
    return 0;
}


int fft_vec_algo_z (double *Y, const double *X, const size_t N, const size_t nfft, const char sc)
{
    if (nfft<N) { fprintf(stderr,"error in fft_vec_algo_z: nfft must be >= N (vec length)\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in fft_vec_algo_z: nfft must be a power of 2\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        if (sc) { *Y = *X * M_SQRT1_2; *++Y = *++X * M_SQRT1_2; }
        else { *Y = *X; *++Y = *++X; }
    }
    else if (nfft==2u)
    {
        if (N==1u) { *Y++ = *X++; *Y++ = *X--; *Y++ = *X++; *Y = *X; }
        else
        {
            *Y++ = *X + *(X+2u); ++X; *Y++ = *X + *(X+2u); ++X;
            *Y++ = *(X-2u) - *X; ++X; *Y = *(X-2u) - *X;
        }
        if (sc)
        {
            const double s = (N==1u) ? M_SQRT1_2 : 0.5;
            *Y-- *= s; *Y-- *= s; *Y-- *= s; *Y *= s;
        }
    }
    else if (nfft==4u)
    {
        for (size_t n=0u; n<N; ++n) { *Y++ = *X++; *Y++ = *X++; }
        for (size_t n=N; n<4u; ++n) { *Y++ = 0.0; *Y++ = 0.0; }
        Y -= 8u;
        const double dcr = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const double dci = Y[1u] + Y[3u] + Y[5u] + Y[7u];
        const double nyr = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        const double nyi = Y[1u] - Y[3u] + Y[5u] - Y[7u];
        const double y1r = Y[0u] - Y[4u] + Y[3u] - Y[7u];
        const double y1i = Y[1u] - Y[2u] + Y[6u] - Y[5u];
        const double y3r = Y[0u] - Y[3u] + Y[7u] - Y[4u];
        const double y3i = Y[1u] - Y[5u] + Y[2u] - Y[6u];
        Y[0u] = dcr; Y[1u] = dci; Y[2u] = y1r; Y[3u] = y1i;
        Y[4u] = nyr; Y[5u] = nyi; Y[6u] = y3r; Y[7u] = y3i;
        if (sc)
        {
            const double s = 1.0/sqrt(2u*N);
            for (size_t n=0u; n<8u; ++n, ++Y) { *Y *= s; }
        }
    }
    else
    {
        const size_t nfft4 = nfft/4u;
        size_t m=0u, h=0u, d, ii, kk, ik;
        double y, sr, si, dr, di;

        //Initialize fft
        double cstbl[nfft+nfft4];
        get_cstbl_d(cstbl,nfft);

        //Copy X into Y
        for (size_t n=0u; n<N; ++n) { *Y++ = *X++; *Y++ = *X++; }
        for (size_t n=N; n<nfft; ++n) { *Y++ = 0.0; *Y++ = 0.0; }
        Y -= 2u*nfft;

        //Bit reverse
        Y += 2u;
        for (size_t n=1u; n<nfft; ++n, Y+=2u)
        {
            size_t k = nfft/2u;
            while (k<=m) { m -= k; k /= 2u; }
            m += k;
            if (n<m)
            {
                y = *Y; *Y = *(Y+2u*(m-n)); *(Y+2u*(m-n)) = y;
                y = *(Y+1u); *(Y+1u) = *(Y+1u+2u*(m-n)); *(Y+1u+2u*(m-n)) = y;
            }
        }
        Y -= 2u*nfft;

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = cstbl[h]; sr = cstbl[h+nfft4];
                for (size_t i=j; i<nfft; i+=kk)
                {
                    ii = i+i; ik = ii + kk;
                    dr = si*Y[ik+1u] + sr*Y[ik];
                    di = sr*Y[ik+1u] - si*Y[ik];
                    Y[ik] = Y[ii] - dr; Y[ii] += dr;
                    Y[ik+1u] = Y[ii+1u] - di; Y[ii+1u] += di;
                }
            }
        }

        //Scale
        if (sc)
        {
            const double s = M_SQRT1_2/sqrt(N);
            for (size_t n=0u; n<2u*nfft; ++n, ++Y) { *Y *= s; }
        }
    }
    
    return 0;
}


#ifdef __cplusplus
}
}
#endif
