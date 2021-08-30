//Basic Teager-Kaiser Energy Operator (TKEO):
//y[n] = x[n]*x[n] - x[n-1]*x[n+1]

//One can also use smooth_diff and smooth_diffdiff, and then:
//y[n] = dx[n]*dx[n] - x[n]*ddx[n]
//See de Matos MC. 2018. Seismic attributes from the complex Teager-Kaiser energy.
//and its excellent reference:
//Holoborodko P. 2008. Smooth noise robust differentiators. www.holoborodko.com.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int tkeo_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int tkeo_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int tkeo_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==N)
    {
        //struct timespec tic, toc; clock_gettime(CLOCK_REALTIME,&tic);
        *Y = *X * *X;
        ++X; ++Y;
        for (size_t l=L-2u; l>0u; --l, ++X, ++Y) { *Y = *X**X + *(X-1)**(X+1); }
        *Y = *X * *X;
        //clock_gettime(CLOCK_REALTIME,&toc);
        //fprintf(stderr,"elapsed time = %.6f ms\n",(double)(toc.tv_sec-tic.tv_sec)*1e3+(double)(toc.tv_nsec-tic.tv_nsec)/1e6);
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y)
            {
                *Y = *X * *X;
                ++X; ++Y;
                for (size_t l=L-2u; l>0u; --l, ++X, ++Y) { *Y = *X**X + *(X-1)**(X+1); }
                *Y = *X * *X;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-K-1u, Y-=K*L-K-1u)
                {
                    *Y = *X * *X;
                    X += K; Y += K;
                    for (size_t l=L-2u; l>0u; --l, X+=K, Y+=K) { *Y = *X**X + *(X-K)**(X+K); }
                    *Y = *X * *X;
                }
            }
        }
    }

    return 0;
}


int tkeo_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in tkeo_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (L==N)
    {
        *Y = *X * *X;
        ++X; ++Y;
        for (size_t l=L-2u; l>0u; --l, ++X, ++Y) { *Y = *X**X + *(X-1)**(X+1); }
        *Y = *X * *X;
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, ++X, ++Y)
            {
                *Y = *X * *X;
                ++X; ++Y;
                for (size_t l=L-2u; l>0u; --l, ++X, ++Y) { *Y = *X**X + *(X-1)**(X+1); }
                *Y = *X * *X;
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(L-1u), Y+=B*(L-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*L-K-1u, Y-=K*L-K-1u)
                {
                    *Y = *X * *X;
                    X += K; Y += K;
                    for (size_t l=L-2u; l>0u; --l, X+=K, Y+=K) { *Y = *X**X + *(X-K)**(X+K); }
                    *Y = *X * *X;
                }
            }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
