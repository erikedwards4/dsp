//Does 1-D FFT (fast Fourier transform) of each vector in X along dim.
//The output Y is complex-valued and has the same size as X,
//except along dim, where Y has length Ly = nfft.

//If sc, then scales Y by sqrt(0.5/n) so that invertible with ifft.

//This uses the algorithm from Ch 12.3 of Numerical Recipes in C, 2nd Ed. [Press et al. 1992].
//...Actually, this is clearly a variant of the algorithm in fft.rad2.c,
//so sticking with that one for now. This one may have the advantage of giving
//only the positive-half of the spectrum, but it does that with a 2nd transformation,
//so I don't think it would be more efficient than copy operations.

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

void fft_1d_s (float *Y, const size_t nfft);
void fft_1d_d (double *Y, const size_t nfft);
//void fft_1d_c (float *Y, const size_t nfft, const size_t *bittbl, const float *cstbl);
//void fft_1d_z (double *Y, const size_t nfft, const size_t *bittbl, const double *cstbl);

int fft_nrc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
int fft_nrc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
//int fft_nrc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
//int fft_nrc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);


void fft_1d_s (float *Y, const size_t nfft, const size_t *bittbl, const float *cstbl)
{
    if (nfft<2u) {}
    else if (nfft==2u)
    {
        const float ny = Y[0u] - Y[2u];
        Y[0u] += Y[2u]; Y[2u] = ny;
    }
    else if (nfft==4u)
    {
        const float dc = Y[0u] + Y[2u] + Y[4u] + Y[6u];
        const float ny = Y[0u] - Y[2u] + Y[4u] - Y[6u];
        Y[3u] = Y[6u] - Y[2u];
        Y[2u] = Y[0u] - Y[4u];
        Y[0u] = dc; Y[4u] = ny;
    }
    else
    {
        const size_t np3 = nfft + 3u;
        size_t i, i1, i2, i3, i4;
        const float c1 = 0.5f, c2 = -0.5f;
        float h1r, h1i, h2r, h2i;
        const double theta = M_PI/(double)(nfft>>1u);
        double wtmp = sin(0.5*theta);
        double wpr = -2.0*wtmp*wtmp, wpi = sin(theta);
        double wr = 1.0+wpr, wi = wpi, np3 = n+3;

        for (size_t i=2u; i<=(nfft>>2u); ++i)
        {
            i1 = i + i - 1u; i2 = i + 1u;
            i3 = np3 - i2; i4 = 1u + i3;
            h1r = c1 * (Y[i1] + Y[i3]);
            h1i = c1 * (Y[i2] - Y[i4]);
            h2r = -c2 * (Y[i2] + Y[i4]);
            h2i =  c2 * (Y[i1] - Y[i3]);
            Y[i1] =  h1r + wr*h2r - wi*h2i;
            Y[i2] =  h1i + wr*h2i + wi*h2r;
            Y[i3] =  h1r - wr*h2r + wi*h2i;
            Y[i4] = -h1i + wr*h2i + wi*h2r;
            wtmp = wr;
            wr += wtmp*wpr - wi*wpi;
            wi += wi*wpr + wtmp*wpi;
        }
        h1r = Y[1u];
        Y[1u] = h1r + Y[2u];
        Y[2u] = h1r - Y[2u];
    }
}


int fft_nrc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    if (nfft>0u && (nfft & (nfft-1u))) { fprintf(stderr,"error in fft_nrc_s: nfft must be a power of 2\n"); return 1; }
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft;
    if (nfft<Lx) { fprintf(stderr,"error in fft_nrc_s: nfft must be >= L (vec length)\n"); return 1; }

    if (nfft==0u || N==0u) {}
    else if (nfft==1u)
    {
        for (size_t n=0u; n<N; ++n, ++X) { *Y++ = *X; *Y++ = 0.0f; }
    }
    else if (Lx==N)
    {
        for (size_t l=0u; l<Lx; ++l, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        for (size_t l=Lx; l<Ly; ++l) { *Y++ = 0.0f; *Y++ = 0.0f; }
        Y -= 2u*Ly;
        fft_1d_s(Y,nfft,bittbl,cstbl);
        free(bittbl); free(cstbl);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v, Y+=2u*Ly)
            {
                for (size_t l=0u; l<Lx; ++l, ++X) { *Y++ = *X; *Y++ = 0.0f; }
                for (size_t l=Lx; l<Ly; ++l) { *Y++ = 0.0f; *Y++ = 0.0f; }
                Y -= 2u*Ly;
                fft_1d_s(Y,nfft,bittbl,cstbl);
            }
        }
        else
        {
            float *Y1;
            if (!(Y1=(float *)malloc(2u*Ly*sizeof(float)))) { fprintf(stderr,"error in fft_nrc_s: problem with malloc. "); perror("malloc"); return 1; }
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=0; b<B; ++b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
                {
                    for (size_t l=0u; l<Lx; ++l, X+=K) { *Y1++ = *X; *Y1++ = 0.0f; }
                    for (size_t l=Lx; l<Ly; ++l) { *Y1++ = 0.0f; *Y1++ = 0.0f; }
                    Y1 -= 2u*Ly;
                    fft_1d_s(Y1,nfft,bittbl,cstbl);
                    for (size_t l=0u; l<Ly; ++l, ++Y1, Y+=2u*K-1u) { *Y = *Y1; *++Y = *++Y1; }
                }
            }
            free(Y1);
        }
        free(bittbl); free(cstbl);
    }

    //Scale
    if (sc)
    {
        const float s = M_SQRT1_2f/sqrtf(Lx);
        for (size_t l=0u; l<2u*Ly*N/Lx; ++l, ++Y) { *Y *= s; }
    }
    
    //clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


#ifdef __cplusplus
}
}
#endif
