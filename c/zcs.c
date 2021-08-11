//The signbit function was definitely slower.
//Also slower was my original while loop idea, and a few other variations.

//For complex numbers, I use zero-crossings of the imaginary part,
//for which one imagines rotation of a vector around the complex plane.

#include <stdio.h>

#ifdef __cplusplus
namespace codee {
extern "C" {
#endif

int zcs_s (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);
int zcs_d (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);
int zcs_c (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);
int zcs_z (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going);


int zcs_s (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcs_s: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int s, sp; //sign info and previous sign info

    if (N==0u) {}
    else if (L==N)
    {
        if (going==0)
        {
            sp = (*X++<0.0f); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0f); *Y = (s!=sp); sp = s; }
        }
        else if (going==1)
        {
            sp = (*X++>=0.0f); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=0.0f); *Y = s*(s!=sp); sp = s; }
        }
        else if (going==-1)
        {
            sp = (*X++<0.0f); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0f); *Y = s*(s!=sp); sp = s; }
        }
        else { fprintf(stderr,"error in zcs_s: going must be in {-1,0,1}\n"); return 1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (going==0)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*X++<0.0f); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0f); *Y = (s!=sp); sp = s; }
                }
            }
            else if (going==1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*X++>=0.0f); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=0.0f); *Y = s*(s!=sp); sp = s; }
                }
            }
            else if (going==-1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*X++<0.0f); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0f); *Y = s*(s!=sp); sp = s; }
                }
            }
            else { fprintf(stderr,"error in zcs_s: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            if (going==0)
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<0.0f); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<0.0f); *Y = (s!=sp); sp = s; }
                    }
                }
            }
            else if (going==1)
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X>=0.0f); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X>=0.0f); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else if (going==-1)
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<0.0f); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<0.0f); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else { fprintf(stderr,"error in zcs_s: going must be in {-1,0,1}\n"); return 1; }
        }
    }

    return 0;
}


int zcs_d (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcs_d: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int s, sp; //sign info and previous sign info

    if (N==0u) {}
    else if (L==N)
    {
        if (going==0)
        {
            sp = (*X++<0.0); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0); *Y = (s!=sp); sp = s; }
        }
        else if (going==1)
        {
            sp = (*X++>=0.0); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=0.0); *Y = s*(s!=sp); sp = s; }
        }
        else if (going==-1)
        {
            sp = (*X++<0.0); *Y++ = 0;
            for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0); *Y = s*(s!=sp); sp = s; }
        }
        else { fprintf(stderr,"error in zcs_d: going must be in {-1,0,1}\n"); return 1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (going==0)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*X++<0.0); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0); *Y = (s!=sp); sp = s; }
                }
            }
            else if (going==1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*X++>=0.0); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X>=0.0); *Y = s*(s!=sp); sp = s; }
                }
            }
            else if (going==-1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*X++<0.0); *Y++ = 0;
                    for (size_t l=1u; l<L; ++l, ++X, ++Y) { s = (*X<0.0); *Y = s*(s!=sp); sp = s; }
                }
            }
            else { fprintf(stderr,"error in zcs_d: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            if (going==0)
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<0.0); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<0.0); *Y = (s!=sp); sp = s; }
                    }
                }
            }
            else if (going==1)
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X>=0.0); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X>=0.0); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else if (going==-1)
            {
                for (size_t g=0u; g<G; ++g, X+=B*(L-1u), Y+=B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=K*L-1u, Y-=K*L-1u)
                    {
                        sp = (*X<0.0); *Y = 0; X+=K; Y+=K;
                        for (size_t l=1u; l<L; ++l, X+=K, Y+=K) { s = (*X<0.0); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else { fprintf(stderr,"error in zcs_d: going must be in {-1,0,1}\n"); return 1; }
        }
    }

    return 0;
}


int zcs_c (int *Y, const float *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcs_c: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int s, sp; //sign info and previous sign info

    if (N==0u) {}
    else if (L==N)
    {
        if (going==0)
        {
            sp = (*++X<0.0f); *Y++ = 0; X+=2;
            for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0f); *Y = (s!=sp); sp = s; }
        }
        else if (going==1)
        {
            sp = (*++X>=0.0f); *Y++ = 0; X+=2;
            for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X>=0.0f); *Y = s*(s!=sp); sp = s; }
        }
        else if (going==-1)
        {
            sp = (*++X<0.0f); *Y++ = 0; X+=2;
            for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0f); *Y = s*(s!=sp); sp = s; }
        }
        else { fprintf(stderr,"error in zcs_c: going must be in {-1,0,1}\n"); return 1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (going==0)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*++X<0.0f); *Y++ = 0; X+=2;
                    for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0f); *Y = (s!=sp); sp = s; }
                }
            }
            else if (going==1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*++X>=0.0f); *Y++ = 0; X+=2;
                    for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X>=0.0f); *Y = s*(s!=sp); sp = s; }
                }
            }
            else if (going==-1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*++X<0.0f); *Y++ = 0; X+=2;
                    for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0f); *Y = s*(s!=sp); sp = s; }
                }
            }
            else { fprintf(stderr,"error in zcs_c: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            if (going==0)
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-1u, Y-=2u*K*L-1u)
                    {
                        sp = (*++X<0.0f); *Y = 0; X+=2u*K; Y+=2u*K;
                        for (size_t l=1u; l<L; ++l, X+=2u*K, Y+=2u*K) { s = (*X<0.0f); *Y = (s!=sp); sp = s; }
                    }
                }
            }
            else if (going==1)
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-1u, Y-=2u*K*L-1u)
                    {
                        sp = (*++X>=0.0f); *Y = 0; X+=2u*K; Y+=2u*K;
                        for (size_t l=1u; l<L; ++l, X+=2u*K, Y+=2u*K) { s = (*X>=0.0f); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else if (going==-1)
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-1u, Y-=2u*K*L-1u)
                    {
                        sp = (*++X<0.0f); *Y = 0; X+=2u*K; Y+=2u*K;
                        for (size_t l=1u; l<L; ++l, X+=2u*K, Y+=2u*K) { s = (*X<0.0f); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else { fprintf(stderr,"error in zcs_c: going must be in {-1,0,1}\n"); return 1; }
        }
    }

    return 0;
}


int zcs_z (int *Y, const double *X, const size_t R, const size_t C, const size_t S, const size_t H, const int iscolmajor, const size_t dim, const int going)
{
    if (dim>3u) { fprintf(stderr,"error in zcs_z: dim must be in [0 3]\n"); return 1; }

    const size_t N = R*C*S*H;
    const size_t L = (dim==0u) ? R : (dim==1u) ? C : (dim==2u) ? S : H;
    int s, sp; //sign info and previous sign info

    if (N==0u) {}
    else if (L==N)
    {
        if (going==0)
        {
            sp = (*++X<0.0); *Y++ = 0; X+=2;
            for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0); *Y = (s!=sp); sp = s; }
        }
        else if (going==1)
        {
            sp = (*++X>=0.0); *Y++ = 0; X+=2;
            for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X>=0.0); *Y = s*(s!=sp); sp = s; }
        }
        else if (going==-1)
        {
            sp = (*++X<0.0); *Y++ = 0; X+=2;
            for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0); *Y = s*(s!=sp); sp = s; }
        }
        else { fprintf(stderr,"error in zcs_z: going must be in {-1,0,1}\n"); return 1; }
    }
    else
    {
        const size_t K = (iscolmajor) ? ((dim==0u) ? 1u : (dim==1u) ? R : (dim==2u) ? R*C : R*C*S) : ((dim==0u) ? C*S*H : (dim==1u) ? S*H : (dim==2u) ? H : 1u);
        const size_t B = (iscolmajor && dim==0u) ? C*S*H : K;
        const size_t V = N/L, G = V/B;

        if (K==1u && (G==1u || B==1u))
        {
            if (going==0)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*++X<0.0); *Y++ = 0; X+=2;
                    for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0); *Y = (s!=sp); sp = s; }
                }
            }
            else if (going==1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*++X>=0.0); *Y++ = 0; X+=2;
                    for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X>=0.0); *Y = s*(s!=sp); sp = s; }
                }
            }
            else if (going==-1)
            {
                for (size_t v=0u; v<V; ++v)
                {
                    sp = (*++X<0.0); *Y++ = 0; X+=2;
                    for (size_t l=1u; l<L; ++l, X+=2, ++Y) { s = (*X<0.0); *Y = s*(s!=sp); sp = s; }
                }
            }
            else { fprintf(stderr,"error in zcs_z: going must be in {-1,0,1}\n"); return 1; }
        }
        else
        {
            if (going==0)
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-1u, Y-=2u*K*L-1u)
                    {
                        sp = (*++X<0.0); *Y = 0; X+=2u*K; Y+=2u*K;
                        for (size_t l=1u; l<L; ++l, X+=2u*K, Y+=2u*K) { s = (*X<0.0); *Y = (s!=sp); sp = s; }
                    }
                }
            }
            else if (going==1)
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-1u, Y-=2u*K*L-1u)
                    {
                        sp = (*++X>=0.0); *Y = 0; X+=2u*K; Y+=2u*K;
                        for (size_t l=1u; l<L; ++l, X+=2u*K, Y+=2u*K) { s = (*X>=0.0); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else if (going==-1)
            {
                for (size_t g=0u; g<G; ++g, X+=2u*B*(L-1u), Y+=2u*B*(L-1u))
                {
                    for (size_t b=0u; b<B; ++b, X-=2u*K*L-1u, Y-=2u*K*L-1u)
                    {
                        sp = (*++X<0.0); *Y = 0; X+=2u*K; Y+=2u*K;
                        for (size_t l=1u; l<L; ++l, X+=2u*K, Y+=2u*K) { s = (*X<0.0); *Y = s*(s!=sp); sp = s; }
                    }
                }
            }
            else { fprintf(stderr,"error in zcs_z: going must be in {-1,0,1}\n"); return 1; }
        }
    }

    return 0;
}


#ifdef __cplusplus
}
}
#endif
