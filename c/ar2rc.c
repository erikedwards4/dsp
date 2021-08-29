//Gets reflection coefficients (RCs) from autoregressive (AR) coeffs for each vector in X.
//This is similar to poly2rc, but no need to normalize by X[0], since AR coeffs assume this.

#include <stdio.h>
#include <stdlib.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ar2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int ar2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int ar2rc_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2rc_s: polynomial input (X) empty\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in ar2rc_s: Lx (length of vecs in X) must be positive\n"); return 1; }

    const size_t P = Lx;
    float sc, *y;
    if (!(y=(float *)malloc(P*sizeof(float)))) { fprintf(stderr,"error in ar2rc_s: problem with malloc. "); perror("malloc"); return 1; }

    if (Lx==N)
    {
        for (size_t p=P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
        for (size_t p=P; p>1u; --p)
        {
            y += p - 1u;
            for (size_t q=p; q>0u; --q, --y) { *y = *--Y; }
            y += p; sc = -*y;
            for (size_t q=p; q>1u; --q, ++Y) { --y; *Y += sc**y; *Y /= 1.0f-sc*sc; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=P-1u)
            {
                for (size_t p=P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
                for (size_t p=P; p>1u; --p)
                {
                    y += p - 1u;
                    for (size_t q=p; q>0u; --q, --y) { *y = *--Y; }
                    y += p; sc = -*y;
                    for (size_t q=p; q>1u; --q, ++Y) { --y; *Y += sc**y; *Y /= 1.0f-sc*sc; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(P-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K-1u)
                {
                    for (size_t p=P; p>0u; --p, X+=K, Y+=K) { *Y = -*X; }
                    for (size_t p=P; p>1u; --p)
                    {
                        y += p - 1u;
                        for (size_t q=p; q>0u; --q, --y) { Y-=K; *y = *Y; }
                        y += p; sc = -*y;
                        for (size_t q=p; q>1u; --q, Y+=K) { --y; *Y += sc**y; *Y /= 1.0f-sc*sc; }
                    }
                }
            }
        }
    }

    //Exit
    return 0;
}


int ar2rc_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2rc_d: polynomial input (X) empty\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in ar2rc_d: Lx (length of vecs in X) must be positive\n"); return 1; }

    const size_t P = Lx;
    double sc, *y;
    if (!(y=(double *)malloc(P*sizeof(double)))) { fprintf(stderr,"error in ar2rc_d: problem with malloc. "); perror("malloc"); return 1; }

    if (Lx==N)
    {
        for (size_t p=P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
        for (size_t p=P; p>1u; --p)
        {
            y += p - 1u;
            for (size_t q=p; q>0u; --q, --y) { *y = *--Y; }
            y += p; sc = -*y;
            for (size_t q=p; q>1u; --q, ++Y) { --y; *Y += sc**y; *Y /= 1.0-sc*sc; }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=P-1u)
            {
                for (size_t p=P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
                for (size_t p=P; p>1u; --p)
                {
                    y += p - 1u;
                    for (size_t q=p; q>0u; --q, --y) { *y = *--Y; }
                    y += p; sc = -*y;
                    for (size_t q=p; q>1u; --q, ++Y) { --y; *Y += sc**y; *Y /= 1.0-sc*sc; }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*(P-1u))
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K-1u)
                {
                    for (size_t p=P; p>0u; --p, X+=K, Y+=K) { *Y = -*X; }
                    for (size_t p=P; p>1u; --p)
                    {
                        y += p - 1u;
                        for (size_t q=p; q>0u; --q, --y) { Y-=K; *y = *Y; }
                        y += p; sc = -*y;
                        for (size_t q=p; q>1u; --q, Y+=K) { --y; *Y += sc**y; *Y /= 1.0-sc*sc; }
                    }
                }
            }
        }
    }

    //Exit
    return 0;
}


int ar2rc_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2rc_s: polynomial input (X) empty\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in ar2rc_s: Lx (length of vecs in X) must be positive\n"); return 1; }

    const size_t P = Lx;
    float scr, sci, sca, *y;
    if (!(y=(float *)malloc(2u*P*sizeof(float)))) { fprintf(stderr,"error in ar2rc_s: problem with malloc. "); perror("malloc"); return 1; }

    if (Lx==N)
    {
        for (size_t p=2u*P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
        for (size_t p=P; p>1u; --p)
        {
            y += 2u*p - 2u;
            for (size_t q=p; q>0u; --q, y-=2) { Y-=2; *y = *Y; *(y+1) = *(Y+1); }
            y += 2u*p; scr = -*y; sci = -*(y+1);
            sca = 1.0f - scr*scr - sci*sci;
            for (size_t q=p; q>1u; --q)
            {
                y -= 2;
                *Y += scr**y - sci**(y+1);
                *Y++ /= sca;
                *Y += scr**(y+1) + sci**y;
                *Y++ /= sca;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=2u*P-2u)
            {
                for (size_t p=2u*P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
                for (size_t p=P; p>1u; --p)
                {
                    y += 2u*p - 2u;
                    for (size_t q=p; q>0u; --q, y-=2) { Y-=2; *y = *Y; *(y+1) = *(Y+1); }
                    y += 2u*p; scr = -*y; sci = -*(y+1);
                    sca = 1.0f - scr*scr - sci*sci;
                    for (size_t q=p; q>1u; --q)
                    {
                        y -= 2;
                        *Y += scr**y - sci**(y+1);
                        *Y++ /= sca;
                        *Y += scr**(y+1) + sci**y;
                        *Y++ /= sca;
                    }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K-2u)
                {
                    for (size_t p=P; p>0u; --p, X+=2u*K, Y+=2u*K) { *Y = -*X; *(Y+1) = -*(X+1); }
                    for (size_t p=P; p>1u; --p)
                    {
                        y += 2u*p - 2u;
                        for (size_t q=p; q>0u; --q, y-=2) { Y-=2u*K; *y = *Y; *(y+1) = *(Y+1); }
                        y += 2u*p; scr = -*y; sci = -*(y+1);
                        sca = 1.0f - scr*scr - sci*sci;
                        for (size_t q=p; q>1u; --q, Y+=2u*K)
                        {
                            y -= 2;
                            *Y += scr**y - sci**(y+1);
                            *Y /= sca;
                            *(Y+1) += scr**(y+1) + sci**y;
                            *(Y+1) /= sca;
                        }
                    }
                }
            }
        }
    }

    //Exit
    return 0;
}


int ar2rc_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2rc_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    if (N==0u) { fprintf(stderr,"error in ar2rc_z: polynomial input (X) empty\n"); return 1; }
    if (Lx<1u) { fprintf(stderr,"error in ar2rc_z: Lx (length of vecs in X) must be positive\n"); return 1; }

    const size_t P = Lx;
    double scr, sci, sca, *y;
    if (!(y=(double *)malloc(2u*P*sizeof(double)))) { fprintf(stderr,"error in ar2rc_z: problem with malloc. "); perror("malloc"); return 1; }

    if (Lx==N)
    {
        for (size_t p=2u*P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
        for (size_t p=P; p>1u; --p)
        {
            y += 2u*p - 2u;
            for (size_t q=p; q>0u; --q, y-=2) { Y-=2; *y = *Y; *(y+1) = *(Y+1); }
            y += 2u*p; scr = -*y; sci = -*(y+1);
            sca = 1.0 - scr*scr - sci*sci;
            for (size_t q=p; q>1u; --q)
            {
                y -= 2;
                *Y += scr**y - sci**(y+1);
                *Y++ /= sca;
                *Y += scr**(y+1) + sci**y;
                *Y++ /= sca;
            }
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=V; v>0u; --v, Y+=2u*P-2u)
            {
                for (size_t p=2u*P; p>0u; --p, ++X, ++Y) { *Y = -*X; }
                for (size_t p=P; p>1u; --p)
                {
                    y += 2u*p - 2u;
                    for (size_t q=p; q>0u; --q, y-=2) { Y-=2; *y = *Y; *(y+1) = *(Y+1); }
                    y += 2u*p; scr = -*y; sci = -*(y+1);
                    sca = 1.0 - scr*scr - sci*sci;
                    for (size_t q=p; q>1u; --q)
                    {
                        y -= 2;
                        *Y += scr**y - sci**(y+1);
                        *Y++ /= sca;
                        *Y += scr**(y+1) + sci**y;
                        *Y++ /= sca;
                    }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*(P-1u))
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*K-2u)
                {
                    for (size_t p=P; p>0u; --p, X+=2u*K, Y+=2u*K) { *Y = -*X; *(Y+1) = -*(X+1); }
                    for (size_t p=P; p>1u; --p)
                    {
                        y += 2u*p - 2u;
                        for (size_t q=p; q>0u; --q, y-=2) { Y-=2u*K; *y = *Y; *(y+1) = *(Y+1); }
                        y += 2u*p; scr = -*y; sci = -*(y+1);
                        sca = 1.0 - scr*scr - sci*sci;
                        for (size_t q=p; q>1u; --q, Y+=2u*K)
                        {
                            y -= 2;
                            *Y += scr**y - sci**(y+1);
                            *Y /= sca;
                            *(Y+1) += scr**(y+1) + sci**y;
                            *(Y+1) /= sca;
                        }
                    }
                }
            }
        }
    }

    //Exit
    return 0;
}


#ifdef __cplusplus
}
}
#endif
