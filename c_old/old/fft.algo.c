//Does 1-D FFT (fast Fourier transform) of each vector in X along dim.
//The output Y is complex-valued and has the same size as X,
//except along dim, where Y has length Ly = nfft/2 + 1 for real-valued X,
//and length Ly = nfft for complex-valued X.

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

void make_bitrev(size_t* bitrev, const size_t N);
void make_sintbl_s (float* sintbl, const size_t N);
void make_sintbl_d (double* sintbl, const size_t N);
int fft_algo_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);
//int fft_algo_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc);


void make_bitrev(size_t* bitrev, const size_t N)
{
    const size_t N2 = N/2u;
    size_t j=0u, k;
    *bitrev++ = 0u;
    for (size_t i=1u; i<N; ++i, ++bitrev)
    {
        k = N2;
        while (k<=j) { j -= k; k /= 2u; }
        j += k;
        *bitrev = j;
    }
    bitrev -= N;
}


void make_sintbl_s (float* sintbl, const size_t N)
{
    const size_t N2=N/2u, N4=N/4u, N8=N/8u;
    float c=1.0f, s=0.0f, dc, ds, t;

    t = (float)(sin(M_PI/(double)N));
    dc = 2.0f * t * t;
    t = 2.0f * dc;
    ds = sqrtf(dc*(2.0f-dc));

    for (size_t i=0u; i<N8; ++i, ++sintbl, s+=ds, ds-=t*s) { *sintbl = s; }
    if (N8>0u) { *sintbl = M_SQRT1_2f; }
    sintbl += N4 - N8;
    for (size_t i=0u; i<N8; ++i, --sintbl, c-=dc, dc+=t*c) { *sintbl = c; }
    sintbl -= N4 - N8;
    for (size_t i=0u; i<N4; ++i) { sintbl[N2-i] = sintbl[i]; }
    for (size_t i=0u; i<N2+N4; ++i) { sintbl[i+N2] = -sintbl[i]; }
}


void make_sintbl_d (double* sintbl, const size_t N)
{
    const size_t N2=N/2u, N4=N/4u, N8=N/8u;
    double c=1.0, s=0.0, dc, ds, t;

    t = sin(M_PI/(double)N);
    dc = 2.0 * t * t;
    t = 2.0 * dc;
    ds = sqrt(dc*(2.0-dc));

    for (size_t i=0u; i<N8; ++i, ++sintbl, s+=ds, ds-=t*s) { *sintbl = s; }
    if (N8>0u) { *sintbl = M_SQRT1_2; }
    sintbl += N4 - N8;
    for (size_t i=0u; i<N8; ++i, --sintbl, c-=dc, dc+=t*c) { *sintbl = c; }
    sintbl -= N4 - N8;
    for (size_t i=0u; i<N4; ++i) { sintbl[N2-i] = sintbl[i]; }
    for (size_t i=0u; i<N2+N4; ++i) { sintbl[i+N2] = -sintbl[i]; }
}


int fft_algo_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
{
    if (dim>3u) { fprintf(stderr,"error in fft_s: dim must be in [0 3]\n"); return 1; }
    struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
    
    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = nfft;
    //const size_t Ly = nfft/2u + 1u;
    if (nfft<Lx) { fprintf(stderr,"error in fft_s: nfft must be >= L (vec length)\n"); return 1; }

    //Initialize fft
    size_t N4=N/4u, j, ik, h, d, k2;
    float t, s, c, dx, dy;
    size_t bitrev[N];
    float sintbl[N+N4];
    make_sintbl_s(sintbl,N); make_bitrev(bitrev,N);

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
        //Copy X to Y
        for (size_t l=0u; l<Lx; ++l, ++X) { *Y++ = *X; *Y++ = 0.0f; }
        for (size_t l=Lx; l<Ly; ++l) { *Y++ = 0.0f; *Y++ = 0.0f; }
        Y -= 2u*Ly;

        //Bit reversal
        // for (size_t i=0u; i<N; ++i, ++bitrev, ++Xr)
        // {
        //     j = *bitrev;
        //     if (i<j) { t = *Xr; *Xr = *(Xr-i+j); *(Xr-i+j) = t; }
        // }
        // bitrev -= N; Xr -= N;

        //Transformation
        // for (size_t k=1u; k<N; k=k2)
        // {
        //     h=0u; k2=k+k; d=N/k2;
        //     for (size_t j=0u; j<k; ++j)
        //     {
        //         c = sintbl[h+N4]; s = sintbl[h];
        //         for (size_t i=j; i<N; i+=k2)
        //         {
        //             ik = i + k;
        //             dx = s*Xi[ik] + c*Xr[ik];
        //             dy = c*Xi[ik] - s*Xr[ik];
        //             Xr[ik] = Xr[i] - dx;
        //             Xr[i] += dx;
        //             Xi[ik] = Xi[i] - dy;
        //             Xi[i] += dy;
        //         }
        //         h += d;
        //     }
        // }

        //Bit reversal
        // for (size_t l=0u; l<N; ++l, ++bitrev, Y+=2u)
        // {
        //     j = *bitrev;
        //     if (l<j) { t = *Y; *Y = *(Y+2u*(j-l)); *(Y+2u*(j-l)) = t; }
        // }
        // bitrev -= N; Y -= 2u*N;
        for (size_t l=0u; l<N; ++l, Y+=2u)
        {
            j = bitrev[l];
            if (l<j) { t = *Y; *Y = *(Y+2u*(j-l)); *(Y+2u*(j-l)) = t; }
        }
        Y -= 2u*N;

        //Transformation
        for (size_t k=1u; k<N; k=k2)
        {
            h=0u; k2=k+k; d=N/k2;
            for (size_t j=0u; j<k; ++j)
            {
                c = sintbl[h+N4]; s = sintbl[h];
                for (size_t i=j; i<N; i+=k2)
                {
                    ik = 2u*(i+k);
                    dx = s*Y[ik+1u] + c*Y[ik];
                    dy = c*Y[ik+1u] - s*Y[ik];
                    Y[ik] = Y[2u*i] - dx;
                    Y[2u*i] += dx;
                    Y[ik+1u] = Y[2u*i+1u] - dy;
                    Y[2u*i+1u] += dy;
                }
                h += d;
            }
        }

        //Finish
        //for (size_t l=0u; l<Ly; ++l, ++Xr, ++Xi) { *Y++ = *Xr; *Y++ = *Xi; }
        //Xr -= Ly; Xi -= Ly; Y -= 2u*Ly;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            // for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            // for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
            // X1 -= nfft;
            // //fftwf_execute(plan);
            // if (sc)
            // {
            //     const float s = (float)(1.0/sqrt(2*Lx));
            //     for (size_t l=0u; l<2*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
            //     Y1 -= 2u*Ly;
            //     for (size_t v=1; v<V; ++v, Y1-=2u*Ly)
            //     {
            //         for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            //         X1 -= Lx;
            //         //fftwf_execute(plan);
            //         for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1 * s; }
            //     }
            // }
            // else
            // {
            //     for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
            //     Y1 -= 2u*Ly;
            //     for (size_t v=1u; v<V; ++v, Y1-=2u*Ly)
            //     {
            //         for (size_t l=0u; l<Lx; ++l, ++X, ++X1) { *X1 = *X; }
            //         X1 -= Lx;
            //         //fftwf_execute(plan);
            //         for (size_t l=0u; l<2u*Ly; ++l, ++Y1, ++Y) { *Y = *Y1; }
            //     }
            // }
        }
        else
        {
            // const float s = sc ? (float)(1.0/sqrt(2*Lx)) : 1.0f;
            // X1 += Lx;
            // for (size_t l=Lx; l<nfft; ++l, ++X1) { *X1 = 0.0f; }
            // X1 -= nfft;
            // for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=2u*B*(Ly-1u))
            // {
            //     for (size_t b=0; b<B; ++b, X-=K*Lx-1u, Y1-=2u*Ly, Y-=2u*K*Ly-2u)
            //     {
            //         for (size_t l=0u; l<Lx; ++l, X+=K, ++X1) { *X1 = *X; }
            //         X1 -= Lx;
            //         //fftwf_execute(plan);
            //         for (size_t l=0u; l<Ly; ++l, ++Y1, Y+=2u*K-1u) { *Y = *Y1*s; *++Y = *++Y1*s; }
            //     }
            // }
        }
    }
    
    //free(bitrev); free(sintbl); free(Xr); free(Xi);
    clock_gettime(CLOCK_REALTIME,&toc); fprintf(stderr,"elapsed time = %.6f ms\n",(toc.tv_sec-tic.tv_sec)*1e3+(toc.tv_nsec-tic.tv_nsec)/1e6);
    return 0;
}


// int fft_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim, const size_t nfft, const char sc)
// {
//     if (dim>3u) { fprintf(stderr,"error in fft_d: dim must be in [0 3]\n"); return 1; }

//     const size_t N = R*C*S*H;
//     const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
//     const size_t Ly = nfft/2u + 1u;
//     if (nfft<Lx) { fprintf(stderr,"error in fft_d: nfft must be >= L (vec length)\n"); return 1; }

//     return 0;
// }


#ifdef __cplusplus
}
}
#endif
