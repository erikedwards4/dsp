//Gets autoregressive (AR) parameters from polynomials for each vec in X.
//If the polynomial is a0 a1 a2..., then the AR coeffs are -a1/a0 -a2/a0...
//This is invertible with ar2poly only if a0==1.

//It is not perfectly clear from online resources how to do the complex case.
//However, since a0 must end up as 1+0i, this implies divide by a0 in any case,
//so that is what is implmented here.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int poly2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);
int poly2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim);


int poly2ar_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2ar_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx - 1u;

    if (N==0u) {}
    else if (Lx==N)
    {
        const float x0 = *X++;
        for (size_t l=0u; l<Ly; ++l, ++X, ++Y) { *Y = *X / x0; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;
        float x0;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                x0 = *X++;
                for (size_t l=0u; l<Ly; ++l, ++X, ++Y) { *Y = *X / x0; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*Ly-1u)
                {
                    x0 = *X; X += K;
                    for (size_t l=0; l<Ly; ++l, X+=K, Y+=K) { *Y = *X / x0; }
                }
            }
        }
    }

    return 0;
}


int poly2ar_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2ar_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx - 1u;

    if (N==0u) {}
    else if (Lx==N)
    {
        const double x0 = *X++;
        for (size_t l=0u; l<Ly; ++l, ++X, ++Y) { *Y = *X / x0; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;
        double x0;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                x0 = *X++;
                for (size_t l=0u; l<Ly; ++l, ++X, ++Y) { *Y = *X / x0; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*Ly-1u)
                {
                    x0 = *X; X += K;
                    for (size_t l=0; l<Ly; ++l, X+=K, Y+=K) { *Y = *X / x0; }
                }
            }
        }
    }

    return 0;
}


int poly2ar_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2ar_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx - 1u;

    if (N==0u) {}
    else if (Lx==N)
    {
        const float x0r = *X++, x0i = *X++;
        const float x0a = x0r*x0r + x0i*x0i;
        for (size_t l=0u; l<Ly; ++l, X+=2)
        {
            *Y++ = (*X*x0r+*(X+1)*x0i) / x0a;
            *Y++ = (*(X+1)*x0r-*X*x0i) / x0a;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;
        float x0r, x0i, x0a;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                x0r = *X++; x0i = *X++;
                x0a = x0r*x0r + x0i*x0i;
                for (size_t l=0u; l<Ly; ++l, X+=2)
                {
                    *Y++ = (*X*x0r+*(X+1)*x0i) / x0a;
                    *Y++ = (*(X+1)*x0r-*X*x0i) / x0a;
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                {
                    x0r = *X++; x0i = *X;
                    x0a = x0r*x0r + x0i*x0i;
                    X += 2u*K - 1u;
                    for (size_t l=0u; l<Ly; ++l, X+=2u*K, Y+=2u*K-1u)
                    {
                        *Y++ = (*X*x0r+*(X+1)*x0i) / x0a;
                        *Y = (*(X+1)*x0r-*X*x0i) / x0a;
                    }
                }
            }
        }
    }

    return 0;
}


int poly2ar_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const char iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in poly2ar_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    const size_t Ly = Lx - 1u;

    if (N==0u) {}
    else if (Lx==N)
    {
        const double x0r = *X++, x0i = *X++;
        const double x0a = x0r*x0r + x0i*x0i;
        for (size_t l=0u; l<Ly; ++l, X+=2)
        {
            *Y++ = (*X*x0r+*(X+1)*x0i) / x0a;
            *Y++ = (*(X+1)*x0r-*X*x0i) / x0a;
        }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;
        double x0r, x0i, x0a;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                x0r = *X++; x0i = *X++;
                x0a = x0r*x0r + x0i*x0i;
                for (size_t l=0u; l<Ly; ++l, X+=2)
                {
                    *Y++ = (*X*x0r+*(X+1)*x0i) / x0a;
                    *Y++ = (*(X+1)*x0r-*X*x0i) / x0a;
                }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*(Ly-1u))
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*K*Ly-2u)
                {
                    x0r = *X++; x0i = *X;
                    x0a = x0r*x0r + x0i*x0i;
                    X += 2u*K - 1u;
                    for (size_t l=0u; l<Ly; ++l, X+=2u*K, Y+=2u*K-1u)
                    {
                        *Y++ = (*X*x0r+*(X+1)*x0i) / x0a;
                        *Y = (*(X+1)*x0r-*X*x0i) / x0a;
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
