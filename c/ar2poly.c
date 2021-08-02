//Gets polynomials from autoregressive (AR) parameters along rows or cols of X.
//If the polynomial is a0 a1 a2..., then the AR coeffs are -a1/a0 -a2/a0...
//Since a0 cannot be recovered from AR coeffs, a0 is always set to 1 here.
//Thus, the output is just: Y = [1 X] (i.e., same as input with leading 1).

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int ar2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);
int ar2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim);


int ar2poly_s (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2poly_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0f;
        for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;
        //const size_t Ly = Lx + 1u;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                *Y++ = 1.0f;
                for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*Lx)
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*Lx+K-1u)
                {
                    *Y = 1.0f; Y += K;
                    for (size_t l=0; l<Lx; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}


int ar2poly_d (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2poly_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0;
        for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                *Y++ = 1.0;
                for (size_t l=0u; l<Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=B*(Lx-1u), Y+=B*Lx)
            {
                for (size_t b=0u; b<B; ++b, X-=K*Lx-1u, Y-=K*Lx+K-1u)
                {
                    *Y = 1.0; Y += K;
                    for (size_t l=0; l<Lx; ++l, X+=K, Y+=K) { *Y = *X; }
                }
            }
        }
    }

    return 0;
}


int ar2poly_c (float *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2poly_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0f; *Y++ = 0.0f;
        for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                *Y++ = 1.0f; *Y++ = 0.0f;
                for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*Lx)
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*(K*Lx+K-1u))
                {
                    *Y = 1.0f; *(Y+1) = 0.0f; Y += 2u*K;
                    for (size_t l=0; l<Lx; ++l, X+=2u*K, Y+=2u*K) { *Y = *X; *(Y+1) = *(X+1); }
                }
            }
        }
    }

    return 0;
}


int ar2poly_z (double *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim)
{
    if (dim>3u) { fprintf(stderr,"error in ar2poly_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t Lx = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;

    if (N==0u) {}
    else if (Lx==N)
    {
        *Y++ = 1.0; *Y++ = 0.0;
        for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/Lx, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            for (size_t v=0u; v<V; ++v)
            {
                *Y++ = 1.0; *Y++ = 0.0;
                for (size_t l=0u; l<2u*Lx; ++l, ++X, ++Y) { *Y = *X; }
            }
        }
        else
        {
            for (size_t g=0u; g<G; ++g, X+=2u*B*(Lx-1u), Y+=2u*B*Lx)
            {
                for (size_t b=0u; b<B; ++b, X-=2u*K*Lx-2u, Y-=2u*(K*Lx+K-1u))
                {
                    *Y = 1.0; *(Y+1) = 0.0; Y += 2u*K;
                    for (size_t l=0; l<Lx; ++l, X+=2u*K, Y+=2u*K) { *Y = *X; *(Y+1) = *(X+1); }
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
