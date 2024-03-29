//Gets polynomial coefficients from reflection coefficients (RCs) for each vector in X.
//Input (X) and output (Y) have the same size,
//except Y has an extra 1.0 at the lead of each vector.

#include <stdio.h>
#include <stdlib.h>
#include "codee_dsp.h"

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif


int rc2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rc2poly_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx + 1u;

    if (N==0u) {}
    else
    {
        float sc, *y;
        if (!(y=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in rc2poly_s: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            *Y++ = 1.0f;
            for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
            Y -= Lx;
            for (size_t l=1u; l<Lx; ++l)
            {
                for (size_t q=l+1u; q>0u; --q, ++Y, ++y) { *y = *Y; }
                Y -= l + 1u; sc = *--y;
                for (size_t q=l; q>0u; --q, ++Y) { --y; *Y += sc * *y; }
                Y -= l;
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=Lx)
                {
                    *Y++ = 1.0f;
                    for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    Y -= Lx;
                    for (size_t l=1u; l<Lx; ++l)
                    {
                        for (size_t q=l+1u; q>0u; --q, ++Y, ++y) { *y = *Y; }
                        Y -= l + 1u; sc = *--y;
                        for (size_t q=l; q>0u; --q, ++Y) { --y; *Y += sc * *y; }
                        Y -= l;
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K-1u)
                    {
                        *Y = 1.0f; Y += K;
                        for (size_t l=Lx; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                        Y -= K*Lx;
                        for (size_t l=1u; l<Lx; ++l)
                        {
                            for (size_t q=l+1u; q>0u; --q, Y+=K, ++y) { *y = *Y; }
                            Y -= K*(l+1u); sc = *--y;
                            for (size_t q=l; q>0u; --q, Y+=K) { --y; *Y += sc * *y; }
                            Y -= K*l;
                        }
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


int rc2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rc2poly_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx + 1u;

    if (N==0u) {}
    else
    {
        double sc, *y;
        if (!(y=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in rc2poly_d: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            *Y++ = 1.0;
            for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
            Y -= Lx;
            for (size_t l=1u; l<Lx; ++l)
            {
                for (size_t q=l+1u; q>0u; --q, ++Y, ++y) { *y = *Y; }
                Y -= l + 1u; sc = *--y;
                for (size_t q=l; q>0u; --q, ++Y) { --y; *Y += sc * *y; }
                Y -= l;
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=Lx)
                {
                    *Y++ = 1.0;
                    for (size_t l=Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    Y -= Lx;
                    for (size_t l=1u; l<Lx; ++l)
                    {
                        for (size_t q=l+1u; q>0u; --q, ++Y, ++y) { *y = *Y; }
                        Y -= l + 1u; sc = *--y;
                        for (size_t q=l; q>0u; --q, ++Y) { --y; *Y += sc * *y; }
                        Y -= l;
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K-1u)
                    {
                        *Y = 1.0; Y += K;
                        for (size_t l=Lx; l>0u; --l, X+=K, Y+=K) { *Y = *X; }
                        Y -= K*Lx;
                        for (size_t l=1u; l<Lx; ++l)
                        {
                            for (size_t q=l+1u; q>0u; --q, Y+=K, ++y) { *y = *Y; }
                            Y -= K*(l+1u); sc = *--y;
                            for (size_t q=l; q>0u; --q, Y+=K) { --y; *Y += sc * *y; }
                            Y -= K*l;
                        }
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


int rc2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rc2poly_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx + 1u;

    if (N==0u) {}
    else
    {
        float scr, sci, *y;
        if (!(y=(float *)malloc(Ly*sizeof(float)))) { fprintf(stderr,"error in rc2poly_c: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            *Y++ = 1.0f; *Y++ = 0.0f;
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
            Y -= 2u*Lx;
            for (size_t l=1u; l<Lx; ++l)
            {
                for (size_t q=2u*l+2u; q>0u; --q, ++y, ++Y) { *y = *Y; }
                Y -= 2u*l + 2u; sci = *--y; scr = *--y;
                for (size_t q=l; q>0u; --q)
                {
                    y -= 2;
                    *Y++ += scr**y - sci**(y+1);
                    *Y++ += scr**(y+1) + sci**y;
                }
                Y -= 2u*l;
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=2u*Lx)
                {
                    *Y++ = 1.0f; *Y++ = 0.0f;
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    Y -= 2u*Lx;
                    for (size_t l=1u; l<Lx; ++l)
                    {
                        for (size_t q=2u*l+2u; q>0u; --q, ++Y, ++y) { *y = *Y; }
                        Y -= 2u*l + 2u; sci = *--y; scr = *--y;
                        for (size_t q=l; q>0u; --q)
                        {
                            y -= 2;
                            *Y++ += scr**y - sci**(y+1);
                            *Y++ += scr**(y+1) + sci**y;
                        }
                        Y -= 2u*l;
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K-2u)
                    {
                        *Y = 1.0f; *(Y+1) = 0.0f; Y += 2u*K;
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, Y+=2u*K) { *Y = *X; *(Y+1) = *(X+1); }
                        Y -= 2u*K*Lx;
                        for (size_t l=1u; l<Lx; ++l)
                        {
                            for (size_t q=l+1u; q>0u; --q, Y+=2u*K) { *y++ = *Y; *y++ = *(Y+1); }
                            Y -= 2u*K*(l+1u); sci = *--y; scr = *--y;
                            for (size_t q=l; q>0u; --q, Y+=2u*K)
                            {
                                y -= 2;
                                *Y += scr**y - sci**(y+1);
                                *(Y+1) += scr**(y+1) + sci**y;
                            }
                            Y -= 2u*K*l;
                        }
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


int rc2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in rc2poly_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx + 1u;

    if (N==0u) {}
    else
    {
        double scr, sci, *y;
        if (!(y=(double *)malloc(Ly*sizeof(double)))) { fprintf(stderr,"error in rc2poly_z: problem with malloc. "); perror("malloc"); return 1; }

        if (Lx==N)
        {
            *Y++ = 1.0; *Y++ = 0.0;
            for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
            Y -= 2u*Lx;
            for (size_t l=1u; l<Lx; ++l)
            {
                for (size_t q=2u*l+2u; q>0u; --q, ++y, ++Y) { *y = *Y; }
                Y -= 2u*l + 2u; sci = *--y; scr = *--y;
                for (size_t q=l; q>0u; --q)
                {
                    y -= 2;
                    *Y++ += scr**y - sci**(y+1);
                    *Y++ += scr**(y+1) + sci**y;
                }
                Y -= 2u*l;
            }
        }
        else
        {
            const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
            const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
            const size_t V = N/Lx, G = V/B;

            if (K==1u && (G==1u || B==1u))
            {
                for (size_t v=V; v>0u; --v, Y+=2u*Lx)
                {
                    *Y++ = 1.0; *Y++ = 0.0;
                    for (size_t l=2u*Lx; l>0u; --l, ++X, ++Y) { *Y = *X; }
                    Y -= 2u*Lx;
                    for (size_t l=1u; l<Lx; ++l)
                    {
                        for (size_t q=2u*l+2u; q>0u; --q, ++Y, ++y) { *y = *Y; }
                        Y -= 2u*l + 2u; sci = *--y; scr = *--y;
                        for (size_t q=l; q>0u; --q)
                        {
                            y -= 2;
                            *Y++ += scr**y - sci**(y+1);
                            *Y++ += scr**(y+1) + sci**y;
                        }
                        Y -= 2u*l;
                    }
                }
            }
            else
            {
                for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
                {
                    for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K-2u)
                    {
                        *Y = 1.0; *(Y+1) = 0.0; Y += 2u*K;
                        for (size_t l=Lx; l>0u; --l, X+=2u*K, Y+=2u*K) { *Y = *X; *(Y+1) = *(X+1); }
                        Y -= 2u*K*Lx;
                        for (size_t l=1u; l<Lx; ++l)
                        {
                            for (size_t q=l+1u; q>0u; --q, Y+=2u*K) { *y++ = *Y; *y++ = *(Y+1); }
                            Y -= 2u*K*(l+1u); sci = *--y; scr = *--y;
                            for (size_t q=l; q>0u; --q, Y+=2u*K)
                            {
                                y -= 2;
                                *Y += scr**y - sci**(y+1);
                                *(Y+1) += scr**(y+1) + sci**y;
                            }
                            Y -= 2u*K*l;
                        }
                    }
                }
            }
        }
        free(y);
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
