//Does 1-D IFFT (inverse fast Fourier transform) of each vector in X along dim.
//The input X is complex-valued.
//The output Y has the same size as X, except along dim, where Y has length nfft.
//Y is real-valued for ifft_rad2_s and ifft_rad2_d,
//and complex-valued for ifft_rad2_c and ifft_rad2_z.
//In the former case, X has only nonnegative freqs, so Lx = nfft/2 + 1.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

//This uses a variant of the radix-2 algorithm, and so only valid for power-of-2 nfft.
//See also ifft.fftw.c.

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#ifndef M_SQRT1_2
    #define M_SQRT1_2 0.707106781186547524401
#endif

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

static void get_bittbl(size_t* bittbl, const size_t nfft);
static void get_cstbl_s (float* cstbl, const size_t nfft);
static void get_cstbl_d (double* cstbl, const size_t nfft);
static void ifft_1d_s (float *Y, const size_t nfft, const size_t *bittbl, const float *cstbl);
static void ifft_1d_d (double *Y, const size_t nfft, const size_t *bittbl, const double *cstbl);
static void ifft_1d_c (float *Y, const size_t nfft, const size_t *bittbl, const float *cstbl);
static void ifft_1d_z (double *Y, const size_t nfft, const size_t *bittbl, const double *cstbl);

int ifft_rad2_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_rad2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_rad2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);
int ifft_rad2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc);


static void get_bittbl(size_t* bittbl, const size_t nfft)
{
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
}


static void get_cstbl_s (float* cstbl, const size_t nfft)
{
    const size_t nfft2=nfft/2u, nfft4=nfft/4u, nfft8=nfft/8u;
    float c=1.0f, s=0.0f, dc, ds, t;

    t = (float)(sin(M_PI/(double)nfft));
    dc = 2.0f * t * t;
    t = 2.0f * dc;
    ds = sqrtf(t-dc*dc);

    for (size_t i=nfft8; i>0u; --i, ++cstbl, s+=ds, ds-=t*s) { *cstbl = s; }
    if (nfft8>0u) { *cstbl = (float)M_SQRT1_2; }
    cstbl += nfft4 - nfft8;
    for (size_t i=nfft8; i>0u; --i, --cstbl, c-=dc, dc+=t*c) { *cstbl = c; }
    cstbl -= nfft4 - nfft8;
    for (size_t i=0u; i<nfft4; ++i) { cstbl[nfft2-i] = cstbl[i]; }
    for (size_t i=0u; i<nfft2+nfft4; ++i) { cstbl[i+nfft2] = -cstbl[i]; }
}


static void get_cstbl_d (double* cstbl, const size_t nfft)
{
    const size_t nfft2=nfft/2u, nfft4=nfft/4u, nfft8=nfft/8u;
    double c=1.0, s=0.0, dc, ds, t;

    t = sin(M_PI/(double)nfft);
    dc = 2.0 * t * t;
    t = 2.0 * dc;
    ds = sqrt(t-dc*dc);

    for (size_t i=nfft8; i>0u; --i, ++cstbl, s+=ds, ds-=t*s) { *cstbl = s; }
    if (nfft8>0u) { *cstbl = M_SQRT1_2; }
    cstbl += nfft4 - nfft8;
    for (size_t i=nfft8; i>0u; --i, --cstbl, c-=dc, dc+=t*c) { *cstbl = c; }
    cstbl -= nfft4 - nfft8;
    for (size_t i=0u; i<nfft4; ++i) { cstbl[nfft2-i] = cstbl[i]; }
    for (size_t i=0u; i<nfft2+nfft4; ++i) { cstbl[i+nfft2] = -cstbl[i]; }
}


static void ifft_1d_s (float *Y, const size_t nfft, const size_t *bittbl, const float *cstbl)
{
    if (nfft<2u) {}
    else if (nfft==2u)
    {
        const float y1r = Y[0u] - Y[2u];
        Y[0u] += Y[2u]; Y[2u] = y1r;
    }
    else if (nfft==4u)
    {
        const float y0r = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const float y1r = Y[0u] - Y[3u] + Y[7u] - Y[4u];
        const float y2r = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        const float y3r = Y[0u] - Y[4u] + Y[3u] - Y[7u];
        Y[0u] = y0r; Y[2u] = y1r; Y[4u] = y2r; Y[6u] = y3r;
    }
    else
    {
        const size_t nfft4 = nfft/4u;
        size_t h=0u, d, ii, kk, ik, b2;
        float y, sr, si, dr, di;

        //Bit reverse
        for (size_t n=1u, n2=2u; n<nfft; ++n, n2+=2u)
        {
            b2 = 2u * bittbl[n];
            if (n2<b2)
            {
                y = Y[n2]; Y[n2] = Y[b2]; Y[b2] = y;
                y = Y[n2+1u]; Y[n2+1u] = Y[b2+1u]; Y[b2+1u] = y;
            }
        }

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = -cstbl[h]; sr = cstbl[h+nfft4];
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
    }
}


static void ifft_1d_d (double *Y, const size_t nfft, const size_t *bittbl, const double *cstbl)
{
    if (nfft<2u) {}
    else if (nfft==2u)
    {
        const double y1r = Y[0u] - Y[2u];
        Y[0u] += Y[2u]; Y[2u] = y1r;
    }
    else if (nfft==4u)
    {
        const double y0r = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const double y1r = Y[0u] - Y[3u] + Y[7u] - Y[4u];
        const double y2r = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        const double y3r = Y[0u] - Y[4u] + Y[3u] - Y[7u];
        Y[0u] = y0r; Y[2u] = y1r; Y[4u] = y2r; Y[6u] = y3r;
    }
    else
    {
        const size_t nfft4 = nfft/4u;
        size_t h=0u, d, ii, kk, ik, b2;
        double y, sr, si, dr, di;

        //Bit reverse
        for (size_t n=1u, n2=2u; n<nfft; ++n, n2+=2u)
        {
            b2 = 2u * bittbl[n];
            if (n2<b2)
            {
                y = Y[n2]; Y[n2] = Y[b2]; Y[b2] = y;
                y = Y[n2+1u]; Y[n2+1u] = Y[b2+1u]; Y[b2+1u] = y;
            }
        }

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = -cstbl[h]; sr = cstbl[h+nfft4];
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
    }
}


static void ifft_1d_c (float *Y, const size_t nfft, const size_t *bittbl, const float *cstbl)
{
    if (nfft<2u) {}
    else if (nfft==2u)
    {
        const float y1r = Y[0u] - Y[2u];
        const float y1i = Y[1u] - Y[3u];
        Y[0u] += Y[2u]; Y[1u] += Y[3u];
        Y[2u] = y1r; Y[3u] = y1i;
    }
    else if (nfft==4u)
    {
        const float y0r = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const float y0i = Y[1u] + Y[3u] + Y[5u] + Y[7u];
        const float y1r = Y[0u] - Y[3u] + Y[7u] - Y[4u];
        const float y1i = Y[1u] - Y[5u] + Y[2u] - Y[6u];
        const float y2r = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        const float y2i = Y[1u] - Y[3u] + Y[5u] - Y[7u];
        const float y3r = Y[0u] - Y[4u] + Y[3u] - Y[7u];
        const float y3i = Y[1u] - Y[2u] + Y[6u] - Y[5u];
        Y[0u] = y0r; Y[1u] = y0i; Y[2u] = y1r; Y[3u] = y1i;
        Y[4u] = y2r; Y[5u] = y2i; Y[6u] = y3r; Y[7u] = y3i;
    }
    else
    {
        const size_t nfft4 = nfft/4u;
        size_t h=0u, d, ii, kk, ik, b2;
        float y, sr, si, dr, di;

        //Bit reverse
        for (size_t n=1u, n2=2u; n<nfft; ++n, n2+=2u)
        {
            b2 = 2u * bittbl[n];
            if (n2<b2)
            {
                y = Y[n2]; Y[n2] = Y[b2]; Y[b2] = y;
                y = Y[n2+1u]; Y[n2+1u] = Y[b2+1u]; Y[b2+1u] = y;
            }
        }

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = -cstbl[h]; sr = cstbl[h+nfft4];
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
    }
}


static void ifft_1d_z (double *Y, const size_t nfft, const size_t *bittbl, const double *cstbl)
{
    if (nfft<2u) {}
    else if (nfft==2u)
    {
        const double y1r = Y[0u] - Y[2u];
        const double y1i = Y[1u] - Y[3u];
        Y[0u] += Y[2u]; Y[1u] += Y[3u];
        Y[2u] = y1r; Y[3u] = y1i;
    }
    else if (nfft==4u)
    {
        const double y0r = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const double y0i = Y[1u] + Y[3u] + Y[5u] + Y[7u];
        const double y1r = Y[0u] - Y[3u] + Y[7u] - Y[4u];
        const double y1i = Y[1u] - Y[5u] + Y[2u] - Y[6u];
        const double y2r = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        const double y2i = Y[1u] - Y[3u] + Y[5u] - Y[7u];
        const double y3r = Y[0u] - Y[4u] + Y[3u] - Y[7u];
        const double y3i = Y[1u] - Y[2u] + Y[6u] - Y[5u];
        Y[0u] = y0r; Y[1u] = y0i; Y[2u] = y1r; Y[3u] = y1i;
        Y[4u] = y2r; Y[5u] = y2i; Y[6u] = y3r; Y[7u] = y3i;
    }
    else
    {
        const size_t nfft4 = nfft/4u;
        size_t h=0u, d, ii, kk, ik, b2;
        double y, sr, si, dr, di;

        //Bit reverse
        for (size_t n=1u, n2=2u; n<nfft; ++n, n2+=2u)
        {
            b2 = 2u * bittbl[n];
            if (n2<b2)
            {
                y = Y[n2]; Y[n2] = Y[b2]; Y[b2] = y;
                y = Y[n2+1u]; Y[n2+1u] = Y[b2+1u]; Y[b2+1u] = y;
            }
        }

        //Transform
        for (size_t k=1u; k<nfft; k=kk, h=0u)
        {
            kk=k+k; d=nfft/kk;
            for (size_t j=0u; j<k; ++j, h+=d)
            {
                si = -cstbl[h]; sr = cstbl[h+nfft4];
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
    }
}


int ifft_rad2_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_rad2_s: dim must be in [0 3]\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in ifft_rad2_s: nfft must be a power of 2\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const float s = (sc) ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_rad2_s: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize ifft
        size_t *bittbl; float *cstbl, *Y1;
        if (!(bittbl=(size_t *)malloc(nfft*sizeof(size_t)))) { fprintf(stderr,"error in ifft_rad2_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(cstbl=(float *)malloc((nfft+nfft/4u)*sizeof(float)))) { fprintf(stderr,"error in ifft_rad2_s: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Y1=(float *)malloc(2u*nfft*sizeof(float)))) { fprintf(stderr,"error in ifft_rad2_s: problem with malloc. "); perror("malloc"); return 1; }
        get_bittbl(bittbl,nfft);
        get_cstbl_s(cstbl,nfft);
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
            X -= 2u + 2u*(1u-nfft%2u);
            for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
            Y1 -= 2u*nfft;
            ifft_1d_s(Y1,nfft,bittbl,cstbl);
            for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t l=Lx; l>0u; --l, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
                X -= 2u + 2u*(1u-nfft%2u);
                for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
                X+=2u*Lx; Y1 -= 2u*nfft;
                ifft_1d_s(Y1,nfft,bittbl,cstbl);
                for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2u*nfft;
                for (size_t v=1u; v<V; ++v, X+=2u*Lx, Y1-=2u*nfft)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
                    X -= 2u + 2u*(1u-nfft%2u);
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
                    Y1 -= 2u*nfft;
                    ifft_1d_s(Y1,nfft,bittbl,cstbl);
                    for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y1-=2u*nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++Y1) { *Y1 = *X; *++Y1 = *(X+1u); }
                        X -= K*(2u + 2u*(1u-nfft%2u));
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
                        Y1 -= 2u*nfft;
                        ifft_1d_s(Y1,nfft,bittbl,cstbl);
                        for (size_t l=nfft; l>0u; --l, Y1+=2u, Y+=K) { *Y = *Y1 * s; }
                    }
                }
            }
        }
        free(bittbl); free(cstbl); free(Y1);
    }

    return 0;
}


int ifft_rad2_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_rad2_d: dim must be in [0 3]\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in ifft_rad2_d: nfft must be a power of 2\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const double s = (sc) ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (Lx!=nfft/2u+1u) { fprintf(stderr,"error in ifft_rad2_d: nfrqs (vec length in X) must equal nfft/2+1\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize ifft
        size_t *bittbl; double *cstbl, *Y1;
        if (!(bittbl=(size_t *)malloc(nfft*sizeof(size_t)))) { fprintf(stderr,"error in ifft_rad2_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(cstbl=(double *)malloc((nfft+nfft/4u)*sizeof(double)))) { fprintf(stderr,"error in ifft_rad2_d: problem with malloc. "); perror("malloc"); return 1; }
        if (!(Y1=(double *)malloc(2u*nfft*sizeof(double)))) { fprintf(stderr,"error in ifft_rad2_d: problem with malloc. "); perror("malloc"); return 1; }
        get_bittbl(bittbl,nfft);
        get_cstbl_d(cstbl,nfft);
    
        if (Lx==N)
        {
            for (size_t l=Lx; l>0u; --l, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
            X -= 2u + 2u*(1u-nfft%2u);
            for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
            Y1 -= 2u*nfft;
            ifft_1d_d(Y1,nfft,bittbl,cstbl);
            for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
            Y1 -= 2u*nfft;
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t l=Lx; l>0u; --l, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
                X -= 2u + 2u*(1u-nfft%2u);
                for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
                X+=2u*Lx; Y1 -= 2u*nfft;
                ifft_1d_d(Y1,nfft,bittbl,cstbl);
                for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
                Y1 -= 2u*nfft;
                for (size_t v=1u; v<V; ++v, X+=2u*Lx, Y1-=2u*nfft)
                {
                    for (size_t l=Lx; l>0u; --l, ++X, ++Y1) { *Y1 = *X; *++Y1 = *++X; }
                    X -= 2u + 2u*(1u-nfft%2u);
                    for (size_t l=nfft-Lx; l>0u; --l, X-=2u, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
                    Y1 -= 2u*nfft;
                    ifft_1d_d(Y1,nfft,bittbl,cstbl);
                    for (size_t l=nfft; l>0u; --l, Y1+=2u, ++Y) { *Y = *Y1 * s; }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X+=2u, Y1-=2u*nfft, Y-=K*nfft-1u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++Y1) { *Y1 = *X; *++Y1 = *(X+1u); }
                        X -= K*(2u + 2u*(1u-nfft%2u));
                        for (size_t l=nfft-Lx; l>0u; --l, X-=2u*K, ++Y1) { *Y1 = *X; *++Y1 = -*(X+1u); }
                        Y1 -= 2u*nfft;
                        ifft_1d_d(Y1,nfft,bittbl,cstbl);
                        for (size_t l=nfft; l>0u; --l, Y1+=2u, Y+=K) { *Y = *Y1 * s; }
                    }
                }
            }
        }
        free(bittbl); free(cstbl); free(Y1);
    }

    return 0;
}


int ifft_rad2_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_rad2_c: dim must be in [0 3]\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in ifft_rad2_c: nfft must be a power of 2\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft;
    const float s = (sc) ? 2.0f*sqrtf(0.5f*(float)nfft)/(float)nfft : 1.0f/(float)nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_rad2_c: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize ifft
        size_t *bittbl; float *cstbl;
        if (!(bittbl=(size_t *)malloc(nfft*sizeof(size_t)))) { fprintf(stderr,"error in ifft_rad2_c: problem with malloc. "); perror("malloc"); return 1; }
        if (!(cstbl=(float *)malloc((nfft+nfft/4u)*sizeof(float)))) { fprintf(stderr,"error in ifft_rad2_c: problem with malloc. "); perror("malloc"); return 1; }
        get_bittbl(bittbl,nfft);
        get_cstbl_s(cstbl,nfft);

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
            Y -= 2u*nfft;
            ifft_1d_c(Y,nfft,bittbl,cstbl);
            for (size_t l=2u*nfft; l>0u; --l, ++Y) { *Y *= s; }
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
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    Y -= 2u*nfft;
                    ifft_1d_c(Y,nfft,bittbl,cstbl);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y) { *Y *= s; }
                }
            }
            else
            {
                float *Y1;
                if (!(Y1=(float *)malloc(2u*Ly*sizeof(float)))) { fprintf(stderr,"error in ifft_rad2_c: problem with malloc. "); perror("malloc"); return 1; }
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++Y1) { *Y1 = *X; *++Y1 = *(X+1u); }
                        Y1 -= 2u*nfft;
                        ifft_1d_c(Y1,nfft,bittbl,cstbl);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K-1u) { *Y = *Y1 * s; *++Y = *++Y1 * s; }
                    }
                }
                free(Y1);
            }
        }
        free(bittbl); free(cstbl);
    }

    return 0;
}


int ifft_rad2_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const size_t nfft, const int sc)
{
    if (dim>3u) { fprintf(stderr,"error in ifft_rad2_z: dim must be in [0 3]\n"); return 1; }
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in ifft_rad2_z: nfft must be a power of 2\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft;
    const double s = (sc) ? 2.0*sqrt(0.5*(double)nfft)/(double)nfft : 1.0/(double)nfft;
    if (Lx!=nfft) { fprintf(stderr,"error in ifft_rad2_z: nfrqs (vec length in X) must equal nfft\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=2u*N; n>0u; --n, ++X, ++Y) { *Y = *X * s; }
    }
    else
    {
        //Initialize ifft
        size_t *bittbl; double *cstbl;
        if (!(bittbl=(size_t *)malloc(nfft*sizeof(size_t)))) { fprintf(stderr,"error in ifft_rad2_z: problem with malloc. "); perror("malloc"); return 1; }
        if (!(cstbl=(double *)malloc((nfft+nfft/4u)*sizeof(double)))) { fprintf(stderr,"error in ifft_rad2_z: problem with malloc. "); perror("malloc"); return 1; }
        get_bittbl(bittbl,nfft);
        get_cstbl_d(cstbl,nfft);

        if (Lx==N)
        {
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
            Y -= 2u*nfft;
            ifft_1d_z(Y,nfft,bittbl,cstbl);
            for (size_t l=2u*nfft; l>0u; --l, ++Y) { *Y *= s; }
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
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    Y -= 2u*nfft;
                    ifft_1d_z(Y,nfft,bittbl,cstbl);
                    for (size_t l=2u*nfft; l>0u; --l, ++Y) { *Y *= s; }
                }
            }
            else
            {
                double *Y1;
                if (!(Y1=(double *)malloc(2u*Ly*sizeof(double)))) { fprintf(stderr,"error in ifft_rad2_z: problem with malloc. "); perror("malloc"); return 1; }
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(nfft-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y1-=2u*nfft, Y-=2u*K*nfft-2u)
                    {
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, ++Y1) { *Y1 = *X; *++Y1 = *(X+1u); }
                        Y1 -= 2u*nfft;
                        ifft_1d_z(Y1,nfft,bittbl,cstbl);
                        for (size_t l=nfft; l>0u; --l, ++Y1, Y+=2u*K-1u) { *Y = *Y1 * s; *++Y = *++Y1 * s; }
                    }
                }
                free(Y1);
            }
        }
        free(bittbl); free(cstbl);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
