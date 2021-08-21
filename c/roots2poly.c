//Polynomial from roots for each vector in X.

#include <stdio.h>
#include <math.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int roots2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int roots2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int roots2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int roots2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int roots2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in roots2poly_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    //const size_t Ly = Lx + 1u;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0f; *Y++ = -*X++;
        for (size_t l=1u; l<Lx; ++l, ++X, Y+=l)
        {
            *Y = -*X * *(Y-1);
            for (size_t p=0u; p<l; ++p) { --Y; *Y = fmaf(-*X,*(Y-1),*Y); }
        }
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
                *Y++ = 1.0f; *Y++ = -*X++;
                for (size_t l=1u; l<Lx; ++l, ++X, Y+=l)
                {
                    *Y = -*X * *(Y-1);
                    for (size_t p=0u; p<l; ++p) { --Y; *Y = fmaf(-*X,*(Y-1),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*Lx)
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K*Lx+K-1u)
                {
                    *Y = 1.0f; Y += K;
                    *Y = -*X; Y += K; X += K;
                    for (size_t l=1u; l<Lx; ++l, X+=K, Y+=l*K)
                    {
                        *Y = -*X * *(Y-K);
                        for (size_t p=0u; p<l; ++p) { Y-=K; *Y = fmaf(-*X,*(Y-K),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int roots2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in roots2poly_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0; *Y++ = -*X++;
        for (size_t l=1u; l<Lx; ++l, ++X, Y+=l)
        {
            *Y = -*X * *(Y-1);
            for (size_t p=0u; p<l; ++p) { --Y; *Y = fma(-*X,*(Y-1),*Y); }
        }
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
                *Y++ = 1.0; *Y++ = -*X++;
                for (size_t l=1u; l<Lx; ++l, ++X, Y+=l)
                {
                    *Y = -*X * *(Y-1);
                    for (size_t p=0u; p<l; ++p) { --Y; *Y = fma(-*X,*(Y-1),*Y); }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=B*(Lx-1u), Y+=B*Lx)
            {
                for (size_t b=B; b>0u; --b, X-=K*Lx-1u, Y-=K*Lx+K-1u)
                {
                    *Y = 1.0; Y += K;
                    *Y = -*X; Y += K; X += K;
                    for (size_t l=1u; l<Lx; ++l, X+=K, Y+=l*K)
                    {
                        *Y = -*X * *(Y-K);
                        for (size_t p=0u; p<l; ++p) { Y-=K; *Y = fma(-*X,*(Y-K),*Y); }
                    }
                }
            }
        }
    }

    return 0;
}


int roots2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in roots2poly_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0f; *Y++ = 0.0f;
        *Y++ = -*X++; *Y++ = -*X++;
        for (size_t l=1u; l<Lx; ++l, X+=2, Y+=2u*l)
        {
            *Y = -*X**(Y-2) + *(X+1)**(Y-1);
            *(Y+1) = -*X**(Y-1) - *(X+1)**(Y-2);
            for (size_t p=0u; p<l; ++p)
            {
                Y -= 2;
                *Y = fmaf(-*X,*(Y-2),fmaf(*(X+1),*(Y-1),*Y));
                *(Y+1) = fmaf(-*X,*(Y-1),fmaf(-*(X+1),*(Y-2),*(Y+1)));
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
            for (size_t v=V; v>0u; --v)
            {
                *Y++ = 1.0f; *Y++ = 0.0f;
                *Y++ = -*X++; *Y++ = -*X++;
                for (size_t l=1u; l<Lx; ++l, X+=2, Y+=2u*l)
                {
                    *Y = -*X**(Y-2) + *(X+1)**(Y-1);
                    *(Y+1) = -*X**(Y-1) - *(X+1)**(Y-2);
                    for (size_t p=0u; p<l; ++p)
                    {
                        Y -= 2;
                        *Y = fmaf(-*X,*(Y-2),fmaf(*(X+1),*(Y-1),*Y));
                        *(Y+1) = fmaf(-*X,*(Y-1),fmaf(-*(X+1),*(Y-2),*(Y+1)));
                    }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*Lx)
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*(K*Lx+K-1u))
                {
                    *Y = 1.0f; *(Y+1) = 0.0f; Y += 2u*K;
                    *Y = -*X; *(Y+1) = -*(X+1); Y += 2u*K; X += 2u*K;
                    for (size_t l=1u; l<Lx; ++l, X+=2u*K, Y+=2u*l*K)
                    {
                        *Y = -*X**(Y-2u*K) + *(X+1)**(Y-2u*K+1u);
                        *(Y+1) = -*X**(Y-2u*K+1u) - *(X+1)**(Y-2u*K);
                        for (size_t p=0u; p<l; ++p)
                        {
                            Y -= 2u*K;
                            *Y = fmaf(-*X,*(Y-2u*K),fmaf(*(X+1),*(Y-2u*K+1u),*Y));
                            *(Y+1) = fmaf(-*X,*(Y-2u*K+1),fmaf(-*(X+1),*(Y-2u*K),*(Y+1)));
                        }
                    }
                }
            }
        }
    }

    return 0;
}


int roots2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in roots2poly_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0; *Y++ = 0.0;
        *Y++ = -*X++; *Y++ = -*X++;
        for (size_t l=1u; l<Lx; ++l, X+=2, Y+=2u*l)
        {
            *Y = -*X**(Y-2) + *(X+1)**(Y-1);
            *(Y+1) = -*X**(Y-1) - *(X+1)**(Y-2);
            for (size_t p=0u; p<l; ++p)
            {
                Y -= 2;
                *Y = fma(-*X,*(Y-2),fma(*(X+1),*(Y-1),*Y));
                *(Y+1) = fma(-*X,*(Y-1),fma(-*(X+1),*(Y-2),*(Y+1)));
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
            for (size_t v=V; v>0u; --v)
            {
                *Y++ = 1.0; *Y++ = 0.0;
                *Y++ = -*X++; *Y++ = -*X++;
                for (size_t l=1u; l<Lx; ++l, X+=2, Y+=2u*l)
                {
                    *Y = -*X**(Y-2) + *(X+1)**(Y-1);
                    *(Y+1) = -*X**(Y-1) - *(X+1)**(Y-2);
                    for (size_t p=0u; p<l; ++p)
                    {
                        Y -= 2;
                        *Y = fma(-*X,*(Y-2),fma(*(X+1),*(Y-1),*Y));
                        *(Y+1) = fma(-*X,*(Y-1),fma(-*(X+1),*(Y-2),*(Y+1)));
                    }
                }
            }
        }
        else
        {
            for (size_t g=G; g>0u; --g, X+=2u*B*(Lx-1u), Y+=2u*B*Lx)
            {
                for (size_t b=B; b>0u; --b, X-=2u*K*Lx-2u, Y-=2u*(K*Lx+K-1u))
                {
                    *Y = 1.0; *(Y+1) = 0.0; Y += 2u*K;
                    *Y = -*X; *(Y+1) = -*(X+1); Y += 2u*K; X += 2u*K;
                    for (size_t l=1u; l<Lx; ++l, X+=2u*K, Y+=2u*l*K)
                    {
                        *Y = -*X**(Y-2u*K) + *(X+1)**(Y-2u*K+1u);
                        *(Y+1) = -*X**(Y-2u*K+1u) - *(X+1)**(Y-2u*K);
                        for (size_t p=0u; p<l; ++p)
                        {
                            Y -= 2u*K;
                            *Y = fma(-*X,*(Y-2u*K),fma(*(X+1),*(Y-2u*K+1u),*Y));
                            *(Y+1) = fma(-*X,*(Y-2u*K+1),fma(-*(X+1),*(Y-2u*K),*(Y+1)));
                        }
                    }
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
